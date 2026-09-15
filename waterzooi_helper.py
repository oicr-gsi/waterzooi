# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 19:17:48 2026

@author: rjovelin
"""


from commons import connect_to_db 

import json


def get_project_level_deliverables(nabu_cache, project_name):
    '''
    (str, str) -> list
    
    Returns a list of project deliverables 
    
    Parameters
    ----------
    - nabu_cache (str): Path to the nabu cache 
    - project_name (str): Name of the project of interest
    '''
        
    conn = connect_to_db(nabu_cache)
    data = conn.execute("SELECT DISTINCT release FROM signoff WHERE project_id = ?;", (project_name,)).fetchall()
    conn.close()
    
    L = []
    for i in data:
        release = json.loads(i['release'])
        for j in release:
            deliverables = release[j].keys()
            L.extend(deliverables)
    L = list(set(L))
    
    return L


def get_release_signoff(nabu_cache, project_name):
    '''
    (str, str) -> dict
    
    Returns a dictionary with the release and release approval signoff
    for all cases in project
    
    Parameters
    ----------
    - nabu_cache (str): Path to the nabu cache 
    - project_name (str): Name of the project of interest
    '''
    
    conn = connect_to_db(nabu_cache)
    data = conn.execute("SELECT DISTINCT case_id, release, release_approval FROM signoff WHERE project_id = ?;", (project_name,)).fetchall()
    conn.close()
    
    D = {}
    
    for i in data:
        case_id = i['case_id']
        release = json.loads(i['release'])
        approval = json.loads(i['release_approval'])
        D[case_id] = {'release': release, 'release_approval': approval}
        
    return D



def get_case_analysis_status(analysis_database, project_name=None):
    '''
    (str, str) -> dict
    
    Returns a dictionary with the analysis status of each case in a project if
    project name is specified or all projects otherwise.
    If a case has multiple analysis templates, it will return the status of a complete
    template if one exists
        
    Parameters
    ----------
    - analysis_database (str): Path to the database storing the analysis data
    - project_name (None str): Name of a specific project
    '''
    
    conn = connect_to_db(analysis_database)
    if project_name:
        data = conn.execute("SELECT project_id, case_id, valid FROM templates WHERE project_id = ?", (project_name,)).fetchall()
    else:
        data = conn.execute("SELECT project_id, case_id, valid FROM templates").fetchall()
    conn.close()
    
    D = {}
    for i in data:
        project = i['project_id']
        case_id = i['case_id']
        valid = i['valid']
        if project not in D:
            D[project] = {}
        assert case_id not in D[project] 
        D[project][case_id] = int(valid)
       
    return D


def count_completed_cases(analysis_status):
    '''
    (dict) -> dict
    
    Returns a dictionary with counts of cases with complete
    and incomplete analysis for each project     
       
    Parameters
    ----------
    - analysis_status (dict): Dictionary with analysis status of each case
                              for each project
    '''
    
    D = {}
    
    for project in analysis_status:
        complete = [case_id for case_id in analysis_status[project] if analysis_status[project][case_id] == 1]
        incomplete = [case_id for case_id in analysis_status[project] if analysis_status[project][case_id] == 0]
        D[project] = {'complete': len(complete), 'incomplete': len(incomplete)}
    
    return D


def extract_samples_libraries_per_case(project_name, database):
    '''
    (str, str) - > dict
    
    Returns a dictionary with samples sorted by tissue type and with libraries sorted by library type
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT case_id, donor_id, sample_id, tissue_type, tissue_origin, library_type, group_id, library FROM Libraries WHERE project_id = ?;", (project_name,)).fetchall()
    conn.close()

    D = {}
    
    for i in data:
        case = i['case_id']
        donor = i['donor_id']
        tissue_type = i['tissue_type']
        library_type = i['library_type']
        library = i['library']
        sample_id = i['sample_id']
        
        sample = '_'.join([donor, i['tissue_origin'] , tissue_type, library_type, i['group_id']])
        
        ### sample_id seems to have wrong format --> correct olive?
        
        if tissue_type == 'R':
            tissue = 'normal'
        else:
            tissue = 'tumor'
        
        if case not in D:
            D[case] = {}
        if 'samples' not in D[case]:
            D[case]['samples'] = {}
        if 'libraries' not in D[case]:
            D[case]['libraries'] = {}
        
        if tissue not in D[case]['samples']:
            D[case]['samples'][tissue] = set()
        if library_type not in D[case]['libraries']:
            D[case]['libraries'][library_type] = set()
            
        D[case]['samples'][tissue].add(sample)
        D[case]['libraries'][library_type].add(library)
        
    return D            
