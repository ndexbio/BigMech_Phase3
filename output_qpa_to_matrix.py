__author__ = 'aarongary'
from os import path
import json

current_directory = path.abspath(path.dirname(__file__))
with open(path.join(current_directory, 'korkut_results.json'), 'r') as korkut_out_file:
    korkut_results = json.load(korkut_out_file)


fh=open('prep_spearman_predictions.txt','w')

korkut_experiments = korkut_results.get("experiments")

fh.write('\t'+'\t'.join([k for k in korkut_experiments.keys()])+'\n')

gene_names = ["BCL2L1", "PTGS2", "PTK2", "STMN1", "MAP2K1", "PLK1", "CDKN1B", "CDKN1A", "ATM", "EGFR", "PRKAA1", "SMAD3", "MAPK14", "AKT1", "AKT2", "ATR", "XRCC1", "GSK3B", "TSC2", "YES1", "CAV1", "STAT3", "PIK3CA", "BID", "CCNE1", "JUN", "RB1", "MAPK3", "CTNNB1", "MAPK1", "RAD51", "PIK3R1", "CDH1", "MYC", "CDH2", "CHEK2", "CHEK1", "PCNA", "SRC", "ACACA", "SHC1", "TP53BP1", "ACACB", "PRKCA", "FYN", "PARP1", "PRKCD", "YWHAB", "AR", "MTOR", "BAX", "EIF4E", "BCL2L11", "BAD", "LYN", "RAF1", "YAP1", "CASP7", "CCNB1", "MAPK9", "TP53", "FN1", "RPS6", "CCND1", "YBX1", "INSR", "ESR1", "GSK3A", "IRS1", "CASP8", "KIT", "MAP2K2", "STAT5A", "ELK1", "EIF4EBP1", "BRAF", "FOXO3", "PDPK1", "RPS6KB1", "LCK"]

for g in gene_names:
    fh.write(g)
    for exp_name in korkut_experiments.keys():
        try:
            d = korkut_experiments.get(exp_name).get("target_paths").get(g)
            if d is not None:
                fh.write('\t'+str(len(d)))
            else:
                fh.write('\t100')
        except KeyError:
            fh.write('\tNA')
    fh.write('\n')







fh.close()

'''

for exp_name in korkut_experiments.keys():
    target_paths = korkut_experiments.get(exp_name).get("target_paths")
    for g in target_paths.keys():
        fh.write(g)



for g in gene_names:
    fh.write(g)
    for k in path_response_dict.keys():
        try:
            fh.write('\t'+str(path_response_dict[k][g]))
        except KeyError:
            fh.write('\tNA')
    fh.write('\n')
'''
