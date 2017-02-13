__author__ = 'danielcarlin'

import argparse
import ndex.client as nc
import ndex.networkn as networkn
import demo_notebooks.causal_paths.causal_utilities as cu
from pandas import read_table

net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

G = networkn.NdexGraph(server='http://public.ndexbio.org', uuid=net_uuid)

probe_genes_mat=read_table('../genes.txt',na_values='none', header=None)

probe_genes=list(probe_genes_mat[0])

drug_desc=read_table('../drugs_hugo.txt',na_values='none')

drug_genes=list()

for target_string in list(drug_desc['HGNC target']):
    targets=target_string.split(',')
    drug_genes.extend(targets)

for affector_string in list(drug_desc['HGNC downstream']):
    affectors=affector_string.split(',')
    drug_genes.extend(affectors)

drug_genes=list(set(drug_genes))

npaths=1
forward_paths = cu.k_shortest_paths_multi(G, drug_genes, probe_genes, npaths=npaths)
reverse_paths = cu.k_shortest_paths_multi(G, probe_genes, drug_genes, npaths=npaths)

def pathlength_cmp(G,path1,path2):
    return(len(path1)-len(path2))


