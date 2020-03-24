import requests
import time
import json

def geneshot_autorif_query(drug_list, filename = ''):
    '''
    Queries a drug list through the Geneshot API to retrieve gene terms associated with each drug
    and saves the results as a json in the input folder
    * Queries against AutoRIF

    Parameters:
    drug_list (list): The list of drugs that are to be queried through the API
    filename (str) : Specify the filename to be saved

    Returns:
    json file in specified folder

    '''
    i = 0
    with open('input/'+filename+'_geneset_autorif.json','a') as outfile:
        outfile.write("[")
        for drug in drug_list:
            i += 1
            try:
                res = requests.get('https://amp.pharm.mssm.edu/geneshot/api/search/auto/' + drug)
                geneset_autorif = res.json()
            except ValueError:
                pass
            json.dump(geneset_autorif, outfile, indent = 4)
            time.sleep(0.50)
            if i < (len(drug_list)):
                outfile.write(",")
            else:
                outfile.write("]")
    outfile.close()
    print(filename + "_geneset_autorif.json file created in input folder!")


def geneshot_generif_query(drug_list, filename = ''):
    '''
    Queries a drug list through the Geneshot API to retrieve gene terms associated with each drug
    and saves the results as a json in the input folder 
    *Queries against GeneRIF

    Parameters:
    drug_list (list): The list of drugs that are to be queried through the API
    filename (str) : Specify the filename to be saved

    Returns:
    json file in specified folder

    '''
    i = 0
    with open('input/'+filename+'_geneset_autorif.json','a') as outfile:
        outfile.write("[")
        for drug in drug_list:
            i += 1
            try:
                res = requests.get('https://amp.pharm.mssm.edu/geneshot/api/search/' + drug)
                geneset_autorif = res.json()
            except ValueError:
                pass
            json.dump(geneset_autorif, outfile, indent = 4)
            time.sleep(0.50)
            if i <= (len(drug_list)):
                outfile.write(",")
            else:
                outfile.write("]")
    outfile.close()
    print(filename + "_geneshot_generif.json file created in input folder!")


