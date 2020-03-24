import pandas as pd
from collections import defaultdict

def transposer(dictionary, cutoff = True):
    '''
    Transposes gene-set library so that gene-set terms become set labels and set labels become set terms

    Parameters:
    dictionary (dict): gene-set library to be transposed

    Returns:
    drugsetlibrary (dict): transposed gene-set in the form of drug-set library

    '''
    d = defaultdict(list)
    for drug,genes in dictionary.items():
        for gene in genes:
            d[gene].append(drug)
    drugsetlibrary = dict(d)
    if cutoff is True:
        # Removing all entries with less than 5 drugs
        drugsetlibrary = {k:v for k,v in drugsetlibrary.items() if len(v) >= 5}
    return drugsetlibrary

def library_counts(drugsetlibrary):
    '''
    Gives counts of the unique drugs, terms, and associations in a drug-set library

    Parameters:
    drugsetlibrary (dict): the input drug-set library to be analyzed

    Returns:
    Number of unique drugs
    Number of unique terms
    Number of associations
    '''
    unique_drugs = len(set([drug for k,v in drugsetlibrary.items() for drug in v]))
    unique_terms = len(drugsetlibrary)
    unique_associations = len([drug for k,v in drugsetlibrary.items() for drug in v])
    average_interactions = unique_associations / unique_terms
    print(str(unique_drugs)+ ' unique drugs')
    print(str(unique_terms)+ ' unique association terms')
    print(str(unique_associations)+ ' unique associations')
    print(str(average_interactions)+ ' average drugs per term')

def csv_formatter(drugsetlibrary, filename = ''):
    '''
    Converts drugsetlibrary into csv format

    Parameters:
    drugsetlibrary (dict): The drugsetlibrary dictionary to be exported
    filename (string) : User-specified filename to be exported
    
    ** Directory of interest must be specified manually

    Returns:
    None
    '''
    
    df = pd.DataFrame.from_dict(drugsetlibrary, orient='index')
    df = df[df.columns[0:]].apply(
        lambda x: ';'.join(x.dropna().astype(str).astype(str)),
        axis = 1
    )
    df.to_csv(filename)

def gmt_formatter(drugsetlibrary, filename = ''):
    '''
    Converts drugsetlibrary into gmt format

    Parameters:
    drugsetlibrary (dict): The drugsetlibrary dictionary to be exported
    filename (string) : User-specified filename to be exported

    ** Directory of interest must be specified manually

    Returns:
    None
    '''
    output = []
    for term in drugsetlibrary.keys():    
        terms = drugsetlibrary[term]
        line = '{0}\t\t{1}'.format(term, '\t'.join(terms))
        output.append(line)

    gmt_output = '\n'.join(output)


    dataFile = open(filename, 'w')
    for eachitem in gmt_output:
        dataFile.write(eachitem)
    dataFile.close()