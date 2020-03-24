import pandas as pd
import numpy as np
from collections import defaultdict

def get_interactionlist(score_cutoff, interactions):
    '''
    Returns an 'interaction list' representation of a dataframe
    This is an nx2 pandas.DataFrame where each row gives a gene and corresponding drug, e.g.
        gene    drug
        MSRA    CHEMBL406270
        TIP1    CHEMBL406270
        EPN1    CHEMBL279107
        TIP1    CHEMBL47181
        ...     ...

    Parameters:
    score_cutoff (int): Specifies cut-off value for interaction_scores between chemicals and proteins
    interactions (dataframe) : input DataFrame of chemical-protein interactions
    '''

    print('removing low-confidence interactions')
    interactions = interactions.loc[interactions['combined_score']>= score_cutoff]

    print('converting STRING protein ids to gene symbols')
    # Use protein lookup table to convert STRING protein ids to Entrez Gene Symbols
    # The .merge pandas function will only combine based on common STITCH protein identifiers between the interactions dataframe and protein lookup table, thus removing
    # STRING protein ids paired to unapproved symbols
    lookup_table_protein = pd.read_csv('input/STRING_to_entrez.csv')
    interactions = lookup_table_protein.merge(interactions)

    print('returning output')
    #Keeping the combined_score column
    interactions = interactions.loc[~interactions['protein'].isnull().values,]
    interactions = interactions.loc[~interactions['drugbank_id'].isnull().values,]
    interactions = interactions.drop_duplicates(subset=['protein','drugbank_id'])

    # Calculating average number of protein interactions for each chemical
    count = interactions[['drugbank_id','gene symbol']].describe()
    print('')
    print(count)
    print('')
    interact_number= str((count.iloc[0]['drugbank_id'])/(count.iloc[1]['drugbank_id']))
    print('Average chemical-protein interactions: ' + interact_number)

    return interactions

def drugsetlibrary_converter(interactions, cutoff = 5):
    # Tupelizing the lists of gene symbols and chemical names so that duplicate protein ids paired to each compound id remain unique 
    id_dict = tuple(zip(interactions['approved symbol'].tolist(), interactions['drugbank_id'].tolist()))

    # Creating a drug-set library where gene symbols are matched to all chemicals with which they are associated 
    drugsetlibrary = defaultdict(list)
    for k, v in id_dict:
        drugsetlibrary[k].append(v)

    # Retaining drug-sets with only specified cutoff number of drugs per set
    drugsetlibrary = {k:list(set(v)) for k,v in drugsetlibrary.items() if len(set(v))>=cutoff} # Cutoff can be user-specified (default = 5)

    return drugsetlibrary
