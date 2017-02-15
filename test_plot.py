# run test analysis on korkut experment structure

import qpa.path_analysis_utilities as pau
from os import path
import json
#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

ndex_host = "http://www.ndexbio.org"
path_rank_method = "shortest_path"
#path_rank_method = "cross_country"

#prior:net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
network_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'

#high_confidence:net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'
#net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

current_directory = path.abspath(path.dirname(__file__))

with open(path.join(current_directory, 'prior_cross_country_results.json'), 'r') as korkut_file:
    resultobj = json.load(korkut_file)

v = []
for key, value in resultobj['experiments'].iteritems():
    if value.get('spearman_rho'):
        v0 = value['spearman_rho']
        v.append(v0)

n, bins, patches = plt.hist(v,20, range=(-1,1))

    #y = mlab.normpdf(bins, mu, sigma)
    #l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('spearman_ranking')
plt.ylabel('Probability')
   # plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
   # plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()