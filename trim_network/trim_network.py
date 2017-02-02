__author__ = 'danielcarlin'

import argparse
import ndex.client as nc
import ndex.networkn as networkn
from pandas import read_table
import sys
import logging
import json
import argparse
import os
import sys
import demo_notebooks.causal_paths.causal_utilities as cu

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

npaths=5
forward_paths = cu.k_shortest_paths_multi(G, drug_genes, probe_genes, npaths=npaths)
reverse_paths = cu.k_shortest_paths_multi(G, probe_genes, drug_genes, npaths=npaths)

G_prime = cu.network_from_paths(G,forward_paths,reverse_paths,drug_genes,probe_genes)

G_prime.set_name('directed path trimmed top 5')
G_prime.upload_to('http://public.ndexbio.org','decarlin','perfect6')
