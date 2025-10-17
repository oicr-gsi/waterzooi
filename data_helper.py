# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 14:47:15 2025

@author: rjovelin
"""

import json
from commons import list_case_workflows




def load_data(provenance_data_file):
    '''
    (str) -> list
    
    Returns the list of data contained in the provenance_data_file
    
    Parameters
    ----------
    - provenance_data_file (str): Path to the file with production data extracted from Shesmu
    '''

    infile = open(provenance_data_file, encoding='utf-8')
    provenance_data = json.load(infile)
    infile.close()
    
    return provenance_data


def clean_up_data(provenance_data):
    '''
    (list) -> list
    
    Returns the list of case information removing any case for which information is missing
            
    Parameters
    ----------
    - provenance_data (list): List of dictionaries, each representing the data of a single case
    '''
    
    to_remove = [i for i in provenance_data if len(i['project_info']) == 0]
    for i in to_remove:
        provenance_data.remove(i)
    
    return provenance_data


def clean_up_workflows(case_data):
    '''
    (dict) -> dict    
    
    Remove children and parent workflows of case workflows that are not in a case
    
    Parameters
    ----------
    case_data (dict): Dictionary with production data of a given case
    '''

    # make a list of case workflows
    workflows = list_case_workflows(case_data)
    
    # evaluate all children and parents workflows
    for k in ['children', 'parents']:
        for i in range(len(case_data['workflow_runs'])):
            L = json.loads(case_data['workflow_runs'][i][k])
            to_remove = [j for j in L if j[1] not in workflows]
            if len(to_remove) >= 1:
                for j in to_remove:
                    assert j in L
                    L.remove(j)
            case_data['workflow_runs'][i][k] = json.dumps(L)            

    return case_data



