# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 15:30:32 2025

@author: rjovelin
"""

import hashlib
import json
from db_helper import connect_to_db



def compute_md5(d):
    '''
    (dict) -> str
    
    Returns the md5 checksum of a dictionary d
    
    Parameters
    ----------
    d (dict): Dictionary with information parsed from FPR
    '''
    
    return hashlib.md5(json.dumps(d, sort_keys=True).encode('utf-8')).hexdigest()




def case_to_update(recorded_md5sums, case, md5sum):
    '''
    (dict, str, str) -> bool
    
    Returns True if the case information needs to be updated in the database
       
    Parameters
    ----------
    - recorded_md5sums (dict): Dictionary of recorded cases, checksum in the database
    - case (str): Case identifier
    - md5sum (str): md5sum of the dictionary containing case information from production
    '''
    
    if case not in recorded_md5sums or recorded_md5sums[case]['md5'] != md5sum:
        # need to update the case info
        return True
    else:
        return False
    
    


def get_cases_md5sum(database, table):
    '''
    (str, str) -> dict

    Returns a dictionary of cases and checksum extracted from the table in the database
           
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - table (str): Table storing the case checksum information 
    '''
        
    # connect to database, get recorded md5sums
    conn = connect_to_db(database)
    data = conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
    tables = [i['name'] for i in data]
    records = {}
    if table in tables:
        data = conn.execute('SELECT case_id, donor_id, md5, project_id FROM {0}'.format(table)).fetchall()  
        for i in data:
            records[i['case_id']] = {'md5': i['md5'], 'project_id': i['project_id'], 'donor_id': i['donor_id']}
    conn.close()
    
    return records



def convert_to_bool(S):
    '''
    (str) -> bool
    
    Returns the boolean value of the string representation of a boolean
    
    Parameters
    ----------
    - S (str): String indicating True or False
    '''
    
    if S.lower() == 'true':
        B = True
    elif S.lower() == 'false':
        B = False
    return B


def find_sequencing_attributes(limskeys, case_data):
    '''
    (list, dict) -> dict
    
    Returns a dictionary of library and barcode for each lims id in limskeys
    
    Parameters
    ----------
    - limskeys (list): list of lims ids
    - case_data (dict): Dictionary with a single case data   
    '''
    
    D = {}
        
    for i in limskeys:
        for d in case_data['sample_info']:
            if i == d['limsId']:
                assert i not in D
                D[i] = {'library': d['library'],
                        'barcode': d['barcode'],
                        'lane': d['lane'],
                        'project': d['project'],
                        'sample': d['sampleId'],
                        'donor': d['donor'],
                        'group_id': d['groupId'],
                        'library_type': d['libraryDesign'],                   
                        'tissue_origin': d['tissueOrigin'],
                        'tissue_type': d['tissueType'],
                        'run': d['run']}
    
    return D


def list_case_workflows(case_data):
    '''
    (dict) -> list
    
    Returns a list of all workflows in a case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with production data of a given case
    '''
    
    L = [d['wfrunid'] for d in case_data['workflow_runs']]
    
    return L


def get_donor_name(case_data):
    '''
    (str) -> str
    
    Returns the name of the donor in case_data
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data 
    '''

    donor = list(set([i['donor'] for i in case_data['sample_info']]))
    assert len(donor) == 1
    donor = donor[0]

    return donor

