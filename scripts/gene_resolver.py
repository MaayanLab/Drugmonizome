import pandas as pd
import numpy as np
import json

def gene_resolver(geneObject, columnName = 'Gene Name', fmt = 'DataFrame'):
    '''
    Verifies approved gene symbols and matches synonyms to approved gene symbols while dropping unresolved genes
    
    Parameters:
    geneObject : The input DataFrame or list with gene information
    columnName (string) : The column within the input df that contains the gene information
    fmt (string) : The format of the input gene information
    '''

    # Importing approved symbol table
    df_lookup = pd.read_csv('/Users/maayanlab/Documents/DrugSetEnrichment/mapping_files/Homo_sapiens.gene_info',
     delimiter = '\t')
    approved_symbols = df_lookup['Symbol'].tolist()

    # Gene synonym lookup
    with open('/Users/maayanlab/Documents/DrugSetEnrichment/mapping_files/gene_symbol_lookup.json', 'r') as f:
        synonym_lookup = json.load(f)

    if fmt == 'DataFrame':
        term_list = []
        for index, row in geneObject.iterrows():
            gene = row.loc[columnName]
            if gene in approved_symbols:
                term_list.append(gene)
            elif gene in synonym_lookup:
                term_list.append(synonym_lookup[gene])
            else:
                geneObject.drop(index, inplace = True)
        geneObject.loc[:,'Approved Symbol'] = pd.Series(np.array(term_list), index = geneObject.index)

    if fmt == 'list':
        df = pd.DataFrame(data = geneObject, columns = ['Gene Name'])
        term_list = []
        for index, row in df.iterrows():
            gene = row.loc['Gene Name']
            if gene in approved_symbols:
                term_list.append(gene)
            elif gene in synonym_lookup:
                term_list.append(synonym_lookup[gene])
            else:
                df.drop(index, inplace = True)
        df.loc[:,'Approved Symbol'] = pd.Series(np.array(term_list), index = df.index)
        genes = df['Approved Symbol'].tolist()
        return genes


