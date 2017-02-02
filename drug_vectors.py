from pandas import read_table
import json
import requests
import math
from scipy.stats import spearmanr

execfile('process_data.py')

drug_desc=read_table('drugs_hugo.txt',na_values='none')

experiments=data['protein']['Sample Description (drug abbre. | dose or time-point)']

target_vectors={}

#prior:net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
#high_confidence:net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'
#net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

headers = {'content-type': 'application/json'}
url = 'http://general.bigmech.ndexbio.org:5602/rank_entities'
heat_response_dict=dict()

for i in experiments:
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
    request_dict={'identifier_set':request_vector,'kernel_id':net_uuid, 'normalize':'True'}
    response=requests.post(url, data=json.dumps(request_dict), headers=headers)
    heat_response_dict[i]=dict(response.json()['ranked_entities'])

#fh=open('prior_only_heat_response.txt','w')
fh=open('prior_heat_response_normalized.txt','w')

fh.write('\t'+'\t'.join([k for k in heat_response_dict.keys()])+'\n')

for g in gene_names:
    fh.write(g)
    for k in heat_response_dict.keys():

        try:
            fh.write('\t'+str(heat_response_dict[k][g]))
        except KeyError:
            fh.write('\tNA')

    fh.write('\n')

fh.close()

def ab2gene(data):
    ab2gene_dict={}
    for r in data['antibody'].iterrows():
        ab2gene_dict[r[1]['Ab Name Reported on Dataset']]=r[1]['Gene Name'].split(',')
    return ab2gene_dict

def comparable_vectors(data,experiment,heat_response_dict,ab2gene_dict):
    heat_vector=[]
    phospho_vector=[]
    genewise_heats=heat_response_dict[experiment]
    for ab in ab2gene_dict.keys():
        try:
            phos_value=abs(math.log(data['protein'].loc[data['protein']['Sample Description (drug abbre. | dose or time-point)'] == experiment][ab].values[0]))
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

#fh2=open('prior_only_spearman_results.txt','w')
fh2=open('prior_spearman_normalized_results.txt','w')

fh2.write('Drug\tSpearman\tp-val\n')

for experiment in data['protein']['Sample Description (drug abbre. | dose or time-point)']:
    spear_tuple=comparable_vectors(data,experiment,heat_response_dict,ab2gene_dict)
    fh2.write(experiment+'\t'+str(spear_tuple[0])+'\t'+str(spear_tuple[1]/2)+'\n')

fh2.close()
