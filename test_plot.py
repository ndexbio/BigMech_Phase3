# run test analysis on korkut experment structure

import qpa.path_analysis_utilities as pau
from os import path
import json
#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math

ndex_host = "http://www.ndexbio.org"
path_rank_method = "shortest_path"
#path_rank_method = "cross_country"

#prior:net_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'
network_uuid='84f321c6-dade-11e6-86b1-0ac135e8bacf'

#high_confidence:net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'
#net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

current_directory = path.abspath(path.dirname(__file__))
# analysis_title = "Prior Shortest Path"
# analysis = 'prior_shortest_path'

analysis_title = "Prior Cross Country"
analysis = 'prior_cross_country'


#input_file = 'prior_cross_country_results.json'
input_file = analysis + '_results.json'

with open(path.join(current_directory, "qpa_results", input_file), 'r') as korkut_file:
    resultobj = json.load(korkut_file)
v = []
for experiment_id, experiment in resultobj['experiments'].iteritems():
    if 'spearman_rho' in experiment:
        rho = experiment['spearman_rho']
        if not math.isnan(rho): #isinstance(rho, Number):
            v.append(rho)
        else:
            print "%s is a unexpected rho value" % (rho)

fig = plt.figure()

print "%s experiments" % (len(v))

n, bins, patches = plt.hist(v, 20, range=(-0.5,0.5), color='#33ccff')
#
# mu = np.mean(v)
# sigma = np.std(v)
#
# print "mu = " + str(mu)
# print "sigma = " + str(sigma)
#
# y = mlab.normpdf(bins, mu, sigma)
# l = plt.plot(bins, y, 'r--', linewidth=2)

plt.xlabel('experiment rho')
plt.ylabel('experiment count')
plt.title('Correlation Distribution for %s Analysis' % (analysis_title))
   # plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
   # plt.axis([40, 160, 0, 0.03])
plt.grid(True)


fig.savefig(path.join(current_directory, "qpa_results", analysis + '_plot.png'))

#plt.show()

