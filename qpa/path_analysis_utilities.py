
import ndex.client as nc
from ndex.networkn import NdexGraph
import pandas
import logs
from pandas import read_table
from scipy.stats import spearmanr
from indra.databases import uniprot_client
from experiment_scoring_helper import ScoreExperiment
# from indra.databases import hgnc_client
# from indra.literature import pubmed_client
# from causal_paths.src.causpaths import DirectedPaths
import BigMechCausalUtils as cu
from causal_paths.src.path_scoring import PathScoring

log = logs.get_logger('bigmech')

# def get_all_gene_names(data):
#     gene_names = data['antibody']['Gene Name']
#     uniprot_ids = data['antibody']['UniProt ID']
#     all_genes = set()
#     invalid_genes = set()
#     for gn, upid in zip(gene_names, uniprot_ids):
#         # Some entries are lists of genes separated by commas
#         # and we also strip off extra spaces
#         names = [x.strip() for x in gn.split(',')]
#         ids = [x.strip() for x in upid.split(',')]
#         names_from_ids = [uniprot_client.get_gene_name(x) for x in ids]
#         # Find invalid gene names
#         for name in names:
#             if not hgnc_client.get_hgnc_id(name):
#                 print('Invalid or deprecated gene symbol: %s' % name)
#                 invalid_genes.add(name)
#         # Find inconsistent gene names and UniProt IDs
#         if set(names) != set(names_from_ids):
#             print('Inconsistent entries:')
#             print('- Given gene names: %s' % ','.join(names))
#             print('- Genes from uniprot IDs: %s' % ','.join(names_from_ids))
#         # Add both the gene names and the gene names derived from UniProt IDs
#         all_genes = all_genes.union(set(names)).union(set(names_from_ids))
#     # Finally remove the invalid gene names
#     all_genes = all_genes.difference(invalid_genes)
#     all_genes = sorted(list(all_genes))
#     return all_genes

def analyze_korkut_batch(network_uuid, ndex_host, path_comparison_method, korkut_structure, experiment_ids=None, use_drug_downstream=False):

    reference_network = NdexGraph(uuid=network_uuid, server=ndex_host)

    # interpret INDRA statements into causal directed edges
    # needs to specify which edges must be doubled to provide both forward and reverse
    two_way_edgetypes = ['Complex', 'Activation', 'in-complex-with']
    cu.indra_causality(reference_network, two_way_edgetypes)

    experiments = korkut_structure["experiments"]
    if not experiment_ids:
        experiment_ids = experiments.keys()

    for experiment_id in experiment_ids:
        print str(experiment_id)
        experiment = experiments[experiment_id]
        # note that making predictions alters the experiment, adding the predictions and any other parameters
        make_predictions(experiment, reference_network, path_comparison_method, use_drug_downstream=use_drug_downstream)

def make_predictions(experiment, network, path_comparison_method, use_drug_downstream=False):
    scoreExperiment = ScoreExperiment()
    sources = experiment["perturbed_protein_gene_symbols"]
    if use_drug_downstream:
        sources.extend(experiment["downstream_protein_gene_symbols"])
    targets = experiment["measured_protein_changes"].keys()

    ps = PathScoring()

    log.info("start: " + str(len(targets)))

    for pcm in path_comparison_method:
        target_to_top_path_map = {}
        target_path_score_map = {}
        no_path_targets = []
        for target in targets:

            paths = cu.get_source_target_paths(network, sources, [target.strip()])

            if len(paths):
                # rank the paths, add top path to map
                if pcm is "cross_country":
                    results_list_sorted = sorted(paths, lambda x,y: ps.cross_country_scoring(x, y))
                    target_to_top_path_map[target] = results_list_sorted[0]
                else:  # default to shortest_path
                    paths.sort(key = lambda s: len(s))
                    target_to_top_path_map[target] = paths[0]
                    target_path_score_map [target] = len(paths[0])
            else:
                no_path_targets.append(target)

        #======================================================
        # If comparison method is "cross_country" we can now
        # process the scores  i.e. it otherwise can't happen
        # until all results for this experiment are recorded
        #=====================================================
        if pcm is "cross_country":
            target_path_score_map = scoreExperiment.rank_score(target_to_top_path_map)

        # rank the targets by top path, producing a target-to-rank dict, i.e. the prediction_dict
        experiment["target_paths"] = target_to_top_path_map

        experiment["target_path_score"] = target_path_score_map

        log.info("Getting spearman score")
        v_changes = []
        v_path_scores = []
        for key, value in experiment["measured_protein_changes"].iteritems():
            if target_path_score_map.get(key):
                v_changes.append(value)
                v_path_scores.append(target_path_score_map[key])
            else:
                print "Target " + key + " has no path score. Ignoring it "

        spearman_rank = spearmanr(v_changes, v_path_scores)
        experiment['spearman_rho'] = float(spearman_rank[0])
        experiment['spearman_pvalue'] = float(spearman_rank[1])
        print "No path on targets(" + str(len(no_path_targets))+ "):" + str(no_path_targets)
        # compute a spearman comparison of the prediction_dict to the measured protein data

# The source_target_list for an experiment specifies a matrix of sources and targets for directed path search
# For each experiment, the sources are either the target proteins of the drugs or
# the union of target proteins and cannonical downstream proteins.
# The targets are the gene names of the measured proteins / phosphoproteins
#
def build_source_target_list(experiment, use_drug_downstream=False):
    sources = experiment["perturbed_protein_gene_symbols"]
    if use_drug_downstream:
        sources.extend(experiment["downstream_protein_gene_symbols"])
    targets = []
    for gene_symbol in experiment["measured_protein_changes"].keys():
        targets.append([gene_symbol])

    source_target = [{"sources": sources, "targets": targets}]

    return source_target

def analyze_korkut(network, experiments, gene_names, drug_desc): #, path_comparator, use_drug_downstream=False):
    path_response_dict=dict()
    results_list = []

    for i in experiments:
        path_response_dict[i]={}
        print(i)
        doses=i.split(',')
        request_vector=[]
        for j in doses:
            drug_dose=j.split('|')
            drug_abbrev=drug_dose[0]
            targets_string=drug_desc.loc[drug_desc['Drug Code'] == drug_abbrev]['HGNC target'].values[0]
            targets=targets_string.split(',')
            downstream_string=drug_desc.loc[drug_desc['Drug Code'] == drug_abbrev]['HGNC downstream'].values[0]
            downstream=downstream_string.split(',')
            drug_request_vector=targets+downstream
            request_vector=request_vector+drug_request_vector

        #send request to our service
        #source = ",".join(request_vector)
        #target = ",".join(gene_names)

        for g in gene_names:
            target=[g]
            max_number_of_paths = 1
           #url = directed_path_query_url + '?source=' + source + '&target=' + target + '&uuid='+net_uuid+'&server=www.ndexbio.org&pathnum=' + str(max_number_of_paths)

            #response=requests.post(url, files={"network_cx": f})

            # TODO use cross country scoring in place of length scoring
            response= cu.k_shortest_paths_multi(network, request_vector, target, npaths=max_number_of_paths)


            pathnum=len(list(response))
            print pathnum
            if pathnum>0:
                # path_response_dict[i][g]=min([len(x) for x in response])
                path_response_dict[i][g]=min([len(x)/2 + 1 for x in results_list])

                print "D: %d CC: %d" % (min([len(x) for x in response]), min([len(x)/2 + 1 for x in results_list]))
            else:
                path_response_dict[i][g]=100
            print path_response_dict[i][g]

    fh=open('path_length_prior_response.txt','w')

    fh.write('\t'+'\t'.join([k for k in path_response_dict.keys()])+'\n')

    for g in gene_names:
        fh.write(g)
        for k in path_response_dict.keys():
            try:
                fh.write('\t'+str(path_response_dict[k][g]))
            except KeyError:
                fh.write('\tNA')
        fh.write('\n')



def spearman_compare_data_to_predictions(data,experiment,prediction_dict,ab2gene_dict):
    score_vector=[]
    phospho_vector=[]
    genewise_scores=prediction_dict[experiment]
    for ab in ab2gene_dict.keys():
        try:
            phos_value=abs(data['protein'].loc[data['protein']['Sample Description (drug abbre. | dose or time-point)'] == experiment][ab].values[0])
            for gene in ab2gene_dict[ab]:
                try:
                    score=genewise_scores[gene]
                    score_vector.append(score)
                    phospho_vector.append(phos_value)
                except KeyError:
                    next
        except KeyError:
            next
    return spearmanr(score_vector, phospho_vector)

def write_spearman_results(method_name, reference_network_name, data, ab2gene_dict, path_response_dict):

    filename = "%s_%s_spearman_results.txt" % (method_name, reference_network_name)

    fh2=open(filename,'w')

    fh2.write('Drug\tSpearman\tp-val\n')

    for experiment in data['protein']['Sample Description (drug abbre. | dose or time-point)']:
        spear_tuple=spearman_compare_data_to_predictions(data,experiment,path_response_dict,ab2gene_dict)
        fh2.write(experiment+'\t'+str(spear_tuple[0])+'\t'+str(spear_tuple[1]/2)+'\n')

    fh2.close()
