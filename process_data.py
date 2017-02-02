import pandas
from indra.databases import uniprot_client
from indra.databases import hgnc_client
from indra.literature import pubmed_client

data_file = 'Korkut et al. Data 12122016.xlsx'

def read_data(fname):
    """Returns the data as a dictionary."""
    data = {}
    data['protein'] = pandas.read_excel(data_file, sheetname='Protein Data',
                                        skiprows=[0], index_col=None)
    data['phenotype'] = pandas.read_excel(data_file,
                                          sheetname='Phenotype Data',
                                          skiprows=[0], index_col=None)
    data['antibody'] = pandas.read_excel(data_file,
                                          sheetname='Antibody Data',
                                          skiprows=range(5), index_col=None)
    return data

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

def get_gene_pmids(genes):
    all_pmids = set()
    for gene in genes:
        print(gene)
        pmids = pubmed_client.get_ids_for_gene(gene)
        all_pmids = all_pmids.union(set(pmids))
    all_pmids = sorted(list(all_pmids))
    return all_pmids

if __name__ == '__main__':
    data = read_data(data_file)
    gene_names = get_all_gene_names(data)
    pmids = get_gene_pmids(gene_names)
