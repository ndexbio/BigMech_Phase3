# create and save a json data structure describing the experiments


import qpa.path_analysis_utilities as pau
from pandas import read_table
from jsonschema import validate
import json
import pandas
from os import path

korkut_experiment_data_filename = 'Korkut et al. Data 12122016.xlsx'
korkut_protein_data_filename = 'rescaled_data/Korkut_recentered.xlsx'
drugs_filename = 'rescaled_data/drugs_hugo.txt'

def ab2gene(data):
    ab2gene_dict={}
    for r in data['antibody'].iterrows():
        gene_names = r[1]['Gene Name'].split(',')
        antibody_name_1 = r[1]['Ab Name Reported on Dataset']
        ab2gene_dict[antibody_name_1] = gene_names
        if "Protein Data ID" in r[1]:
            antibody_name_2 = r[1]['Protein Data ID']
            ab2gene_dict[antibody_name_2]= gene_names
    return ab2gene_dict

def read_korkut_data(korkut_protein_data_file, korkut_experiment_data_file):
    """Returns the data as a dictionary."""
    data = {}
    #note: using the recentered data here
    data['protein'] = pandas.read_excel(korkut_protein_data_file,
                                        sheetname='Sheet2',
                                        skiprows=[0],
                                        index_col=None)
    data['phenotype'] = pandas.read_excel(korkut_experiment_data_file,
                                          sheetname='Phenotype Data',
                                          skiprows=[0],
                                          index_col=None)
    data['antibody'] = pandas.read_excel(korkut_experiment_data_file,
                                          sheetname='Antibody Data',
                                          skiprows=range(5),
                                         index_col=None)
    return data

def build_korkut(data, drug_desc):
    antibody_to_gene_from_data = ab2gene(data)
    korkut_dict = {}
    experiment_ids = data['protein']['Sample Description (drug abbre. | dose or time-point)']
    first_experiment_data = data['protein'].loc[data['protein']['Sample Description (drug abbre. | dose or time-point)'] == experiment_ids[0]]
    antibody_to_change_map = first_experiment_data.to_dict()
    antibody_to_gene_symbol_map = {}
    unmapped_antibodies = []
    for antibody in antibody_to_change_map.keys():
        if antibody in antibody_to_gene_from_data:
            antibody_to_gene_symbol_map[antibody] = antibody_to_gene_from_data[antibody]
        if antibody is not 'Sample Description (drug abbre. | dose or time-point)' and antibody is not "experimental sets":
            unmapped_antibodies.append(antibody)
    korkut_dict["unmapped_antibodies"] = unmapped_antibodies
    print "%s unmapped antibodies" % (len(unmapped_antibodies))
    korkut_dict["antibody_to_gene_symbol_map"] = antibody_to_gene_symbol_map
    experiments = {}
    for experiment_id in experiment_ids:
        experiments[experiment_id] = build_experiment(experiment_id, data, drug_desc, antibody_to_gene_symbol_map)

    korkut_dict["experiments"] = experiments
    return korkut_dict

def build_experiment(experiment_id, data, drug_desc, antibody_to_gene_symbol_map):
    experiment = {"id": experiment_id}
    doses=experiment_id.split(',')
    perturbed = []
    downstream = []
    drugs = []
    for dose in doses:
        drug_dose=dose.split('|')
        drug_abbrev=drug_dose[0]
        drugs.append(drug_abbrev)
        targets_string=drug_desc.loc[drug_desc['Drug Code'] == drug_abbrev]['HGNC target'].values[0]
        if not targets_string is "no_gene":
            perturbed.extend(targets_string.split(','))
        downstream_string=drug_desc.loc[drug_desc['Drug Code'] == drug_abbrev]['HGNC downstream'].values[0]
        downstream.extend(downstream_string.split(','))

    experiment["drug_ids"] = drugs
    experiment["perturbed_protein_gene_symbols"]=perturbed
    experiment["downstream_protein_gene_symbols"]=downstream
    experiment["cell_viability"] = data["phenotype"].loc[data["phenotype"]['Sample Description (drug abbre. | dose or time-point)'] == experiment_id]["cellviab"].values[0]

    try:
        experiment_data = data['protein'].loc[data['protein']['Sample Description (drug abbre. | dose or time-point)'] == experiment_id]
    except KeyError:
         raise Exception("no protein data for %s" % (experiment_id))

    antibody_to_change_map = experiment_data.to_dict(orient='list')
    protein_changes = {}

    for antibody, gene_symbols in antibody_to_gene_symbol_map.iteritems():
        if antibody in antibody_to_change_map.keys():
            antibody_change = antibody_to_change_map[antibody]
            # getting numeric change and taking absolute value
            absolute_antibody_change = abs(antibody_change[0])
            for gene_symbol in gene_symbols:
                # TODO what do we do in case of conflict?
                protein_changes[gene_symbol] = absolute_antibody_change
        else:
            print "no data for antibody %s" % (antibody)

    experiment["measured_protein_changes"] = protein_changes

    return experiment


current_directory = path.abspath(path.dirname(__file__))

with open(path.join(current_directory, 'korkut_experiments_schema.json')) as json_file:
            korkut_schema = json.load(json_file)

data = read_korkut_data(korkut_protein_data_filename, korkut_experiment_data_filename)

print "read korkut data"

drug_desc=read_table(drugs_filename,na_values='none')

print "read drug descriptions"

#gene_names = pau.get_all_gene_names(data)

korkut_experiments = build_korkut(data, drug_desc)

#validate(korkut_experiments, korkut_schema)

with open(path.join(current_directory, 'korkut_experiments.json'), 'w') as korkut_file:
    json.dump(korkut_experiments, korkut_file, indent=4)


