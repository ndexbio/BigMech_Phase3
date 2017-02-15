# run test analysis on korkut experment structure

import qpa.path_analysis_utilities as pau
from os import path
import json
import csv
import logs

log = logs.get_logger('bigmech')

ndex_host = "http://www.ndexbio.org"
path_rank_method = ["cross_country" ,"shortest_path"]

#prior:net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
#network_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'

#high_confidence:net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'
#net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

reference_networks = {
    "prior": "84f321c6-dade-11e6-86b1-0ac135e8bacf"
  # "high_confidence": 'b04e406b-dc88-11e6-86b1-0ac135e8bacf'
    }


current_directory = path.abspath(path.dirname(__file__))


experiment_ids = []
with open (path.join(current_directory, 'experiment_id_list.csv'), 'r') as exp_list:
    tsvin = csv.reader(exp_list, delimiter='\t')
    counter = 0
    for row in tsvin:
        if counter>0:
            experiment_ids.append(row[0])
        counter +=1

print  str(counter) + " experiments ids will be processed"


for network_name, network_uuid in reference_networks.iteritems():
    log.info("processing network " + network_name)

    for pcm in path_rank_method:
        with open(path.join(current_directory, 'korkut_experiments.json'), 'r') as korkut_file:
            korkut = json.load(korkut_file)

        pau.analyze_korkut_batch(network_uuid, ndex_host, pcm, korkut, experiment_ids=experiment_ids)  # ["901|1.5,Tm|0.3", "901|1.5,HN|6", "901|1.5"])

        with open(path.join(current_directory, network_name + '_' + pcm + '_results.json'), 'w') as results_file:
            json.dump(korkut, results_file, indent=4)