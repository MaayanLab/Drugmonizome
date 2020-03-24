import requests
import time
from collections import defaultdict

def enrichr_id_retrieval(dictionary):
    '''
    Queries a genelist (value of input dictionary) through the Enrichr API and creates an output dictionary of your input dictionary key matched to the retrieved EnrichrID

    Parameters:
    dictionary (dict): Key (e.g. drug) matched to value (genelist)

    Returns:
    output_dict (dict): Key (e.g. drug) matched to value (EnrichrID)
    '''

    # Associating downregulated genelists with userListIDs in Enrichr #
    enrichr_url = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    failed_list = []
    output_dict = {}

    for term,genelist in dictionary.items():
        genes_str = '\n'.join(genelist)
        
        payload = {'list': (None, genes_str)}
        response = requests.post(enrichr_url, files=payload)
        
        if not response.ok:
            failed_list.append(term)
        
        user_id = response.json()['userListId']
        output_dict[term] = user_id
        time.sleep(0.5)
    
    print(str(len(failed_list))+ " genelists failed to be matched with a EnrichrID!")
    return output_dict

def enrichr_library_generator(gene_set_library, dictionary, p_value_cutoff = 0.01):
    '''
    Queries an EnrichrID and retrieves most enriched terms from specificed geneset library

    Parameters:
    gene_set_library (str): Specify the Enrichr library to pull enriched terms from
    dictionary (dict): Dictionary of terms (e.g. drugs) associated with an EnrichrID (and in effect, a genelist/signature)

    Returns:
    drugsetlibrary (dict): dictionary of Enrichr terms matched to sets of drugs signficantly associated with the term

    '''
    
    #################################################################
    ### Querying up/downregulated drug/gene lists through Enrichr ###
    #################################################################
    
    enrichr_url = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    drug_list = []
    term_list = []

    for drug, user_id in dictionary.items():
        response = requests.get(enrichr_url + query_string % (user_id, gene_set_library))
        try:
            response.json()
        except ValueError:
            continue
        data = response.json()
        # Replacing JSON key with drug name for easier data manipulation later #
        data[drug] = data.pop(gene_set_library)
        # Accessing each branch of json tree
        for k,v in data.items():
            for lists in v:
                if (lists[6]) < p_value_cutoff: # using a p value of 0.01 as a strict threshold
                    drug_list.append(k)
                    term_list.append(lists[1])

        time.sleep(0.5)
    

    ####################################
    ### Creating Drug-set library ###
    ####################################

    # Creating two tupelized lists in a dictionary format #
    drug_dict = tuple(zip(term_list,drug_list))
    # Creating a dictionary where a list of values are matched under their corresponding key #
    drugsetlibrary = defaultdict(list)
    for k,v in drug_dict:
        drugsetlibrary[k].append(v)

    # Removing all terms paired with less than 5 drugs #
    drugsetlibrary = {k : v for k,v in drugsetlibrary.items() if len(v)>=5}
    return drugsetlibrary