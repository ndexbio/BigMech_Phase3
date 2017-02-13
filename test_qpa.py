# run test analysis on korkut experment structure

import qpa.path_analysis_utilities as pau
from os import path
import json

ndex_host = "http://www.ndexbio.org"
path_rank_method = "shortest_path"
#path_rank_method = "cross_country"

#prior:net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
network_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'

#high_confidence:net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'
#net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

current_directory = path.abspath(path.dirname(__file__))

with open(path.join(current_directory, 'korkut_experiments.json'), 'r') as korkut_file:
    korkut = json.load(korkut_file)

pau.analyze_korkut_batch(network_uuid, ndex_host, path_rank_method, korkut, experiment_ids=["901|1.5,Tm|0.3"])

with open(path.join(current_directory, 'korkut_results.json'), 'w') as results_file:
    json.dump(korkut, results_file, indent=4)