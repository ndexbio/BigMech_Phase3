
import ndex.client as nc
import ndex.networkn as networkn
import pandas
from pandas import read_table
from scipy.stats import spearmanr
#from indra.databases import uniprot_client
#from indra.databases import hgnc_client
#from indra.literature import pubmed_client
from causal_paths.src.causpaths import DirectedPaths
import demo_notebooks.causal_paths.causal_utilities as cu
from causal_paths.src.path_scoring import PathScoring


host = "http://www.ndexbio.org"
directed_path_query_url = "http://general.bigmech.ndexbio.org:5603/directedpath/query"

ndex = nc.Ndex(host=host)

data_file = '../Korkut et al. Data 12122016.xlsx'


def read_data(fname):
    """Returns the data as a dictionary."""
    data = {}
    #note: using the recentered data here
    data['protein'] = pandas.read_excel('../rescaled_data/Korkut_recentered.xlsx', sheetname='Sheet2',
                                        skiprows=[0], index_col=None)
    data['phenotype'] = pandas.read_excel(data_file,
                                          sheetname='Phenotype Data',
                                          skiprows=[0], index_col=None)
    data['antibody'] = pandas.read_excel(data_file,
                                          sheetname='Antibody Data',
                                          skiprows=range(5), index_col=None)
    return data

data = read_data(data_file)

def get_all_gene_names(data):
    gene_names = data['antibody']['Gene Name']
    uniprot_ids = data['antibody']['UniProt ID']
    all_genes = set()
    invalid_genes = set()
    for gn, upid in zip(gene_names, uniprot_ids):
        # Some entries are lists of genes separated by commas                                           
        # and we also strip off extra spaces                                                            
        names = [x.strip() for x in gn.split(',')]
        ids = [x.strip() for x in upid.split(',')]
        names_from_ids = [uniprot_client.get_gene_name(x) for x in ids]
        # Find invalid gene names                                                                       
        for name in names:
            if not hgnc_client.get_hgnc_id(name):
                print('Invalid or deprecated gene symbol: %s' % name)
                invalid_genes.add(name)
        # Find inconsistent gene names and UniProt IDs                                                  
        if set(names) != set(names_from_ids):
            print('Inconsistent entries:')
            print('- Given gene names: %s' % ','.join(names))
            print('- Genes from uniprot IDs: %s' % ','.join(names_from_ids))
        # Add both the gene names and the gene names derived from UniProt IDs                           
        all_genes = all_genes.union(set(names)).union(set(names_from_ids))
    # Finally remove the invalid gene names                                                             
    all_genes = all_genes.difference(invalid_genes)
    all_genes = sorted(list(all_genes))
    return all_genes

gene_names = get_all_gene_names(data)

drug_desc=read_table('../drugs_hugo.txt',na_values='none')

experiments=data['protein']['Sample Description (drug abbre. | dose or time-point)']

target_vectors={}

headers = {'content-type': 'application/json'}
url = 'http://general.bigmech.ndexbio.org:5602/rank_entities'
heat_response_dict=dict()

#prior:net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
#high_confidence:net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'
#net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

network = networkn.NdexGraph(server='http://public.ndexbio.org', uuid=net_uuid)

host = "http://www.ndexbio.org"
directed_path_query_url = "http://general.bigmech.ndexbio.org:5603/directedpath/query"

ndex = nc.Ndex(host=host)

reference_network_cx = ndex.get_network_as_cx_stream(net_uuid)

path_response_dict={}

directedPaths=DirectedPaths()

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

def ab2gene(data):
    ab2gene_dict={}
    for r in data['antibody'].iterrows():
        ab2gene_dict[r[1]['Ab Name Reported on Dataset']]=r[1]['Gene Name'].split(',')
    return ab2gene_dict

ab2gene_dict=ab2gene(data)


def comparable_vectors(data,experiment,heat_response_dict,ab2gene_dict):
    heat_vector=[]
    phospho_vector=[]
    genewise_heats=heat_response_dict[experiment]
    for ab in ab2gene_dict.keys():
        try:
            phos_value=abs(data['protein'].loc[data['protein']['Sample Description (drug abbre. | dose or time-point)'] == experiment][ab].values[0])
            for gene in ab2gene_dict[ab]:
                try:
                    heat=genewise_heats[gene]
                    heat_vector.append(heat)
                    phospho_vector.append(phos_value)
                except KeyError:
                    next
        except KeyError:
            next
    return spearmanr(heat_vector, phospho_vector)

fh2=open('shortest_path_prir_spearman_results.txt','w')

fh2.write('Drug\tSpearman\tp-val\n')

for experiment in data['protein']['Sample Description (drug abbre. | dose or time-point)']:
    spear_tuple=comparable_vectors(data,experiment,path_response_dict,ab2gene_dict)
    fh2.write(experiment+'\t'+str(spear_tuple[0])+'\t'+str(spear_tuple[1]/2)+'\n')

fh2.close()
