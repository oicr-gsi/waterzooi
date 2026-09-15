# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 15:28:58 2026

@author: rjovelin
"""


import sqlite3
import requests
import argparse
import json
import os
from commons import load_data, is_case_info_complete, connect_to_db, \
    insert_multiple_records, delete_multiple_records



def ticket_format(d):
    '''
    (dict) -> list | None
    
    Transforms the comment for a file qc into a list of tickets and/or a list of comments or None
        
    Parameters
    ----------
    - d (dict): Dictionary with file qc info pulled from the fileqc nabu endpoint
    '''
    
    comment = d['comment']
    
    if comment:
        if comment.startswith('G') and '-' in comment:
            comment = comment.split('-')
            c = ['-'.join([comment[0], comment[i]]) for i in range(1, len(comment))]
        else:
            c = [d['comment']]
    else:
        c = d['comment']
    
    return c



def extract_nabu_signoff(nabu_key_file, nabu_endpoint='https://nabu.gsi.oicr.on.ca/case/sign-off'):
    '''
    (str, str) -> dict
    
    Returns a dictionary of signoffs for each case in cases
        
    Parameters
    ----------
    - nabu_key_file (str): File storing the nabu API key
    - nabu_endpoint (str): URL to access the signoffs in Nabu
    '''
    
    infile = open(nabu_key_file)
    nabu_key = infile.read().rstrip()
    infile.close()
    
    headers = {'accept': 'application/json',
               'X-API-KEY': nabu_key,}
    
    D = {}
    
    response = requests.get(nabu_endpoint, headers=headers)
    if response.status_code == 200:
        for d in response.json():
            case_id = d['caseIdentifier']
            ticket = ticket_format(d)
            d['comment'] = ticket
            if case_id not in D:
                D[case_id] = {}
            step = d['signoffStepName']
            step = ' '.join(list(map(lambda x: x.lower().capitalize(), step.split('_'))))
            if step in D[case_id]:
                D[case_id][step].append(d)
            else:
                D[case_id][step] = [d]
    return D




def get_file_signoff(project, nabu_endpoint = 'https://nabu.gsi.oicr.on.ca/get-fileqcs'):
    '''
    (str, str) -> dict

    Returns a dictionary 

    Parameters
    ----------
    - project (str): Project of interest
    - nabu_endpoint (str): URL to access the file qc in Nabu
    '''

    headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
    json_data = {'project': project}
    response = requests.post(nabu_endpoint, headers=headers, json=json_data)
    
    D = {}
    
    if response.status_code == 200:
        for d in response.json()['fileqcs']:
            fileid = d['fileid']
            filepath = d['filepath']
            if 'username' in d:
                username = d['username']
            else:
                username = 'NA'
            qcstatus = d['qcstatus']
            if 'comment' in d:
                ticket = d['comment']
            else:
                ticket = 'NA'
            
            assert fileid not in D
            D[fileid] = {
                'fileid' : fileid,
                'filepath' : filepath,
                'username' : username,
                'qcstatus' : qcstatus,
                'ticket' : ticket}
            
    return D   



def create_nabu_cache(database, fields):
    '''
    (str, dict) -> None
    
    Creates a table in database
    
    Parameters
    ----------
    - database_name (str): Name of the database
    - fields (dict): Dctionary with table names and columns and data type
    '''
    
    tables = ['fileqc', 'signoff']
    
    # connect to database
    conn = sqlite3.connect(database)
    cur = conn.cursor()

    for table in tables:
        # get the column names and types
        column_names, column_types = fields[table]['names'], fields[table]['types']    
        # define table format including constraints    
        table_format = ', '.join(list(map(lambda x: ' '.join(x), list(zip(column_names, column_types)))))
        # create table
        cmd = 'CREATE TABLE {0} ({1})'.format(table, table_format)
        cur.execute(cmd)
        conn.commit()
    
    conn.close()



def list_projects(provenance_data_file):
    '''
    (str)- > list
    
    Returns a list of project names across all valid cases
        
    Parameters
    - provenance_data_file (str): Path to the provenance reporter json file
    '''
    
    L = []
    
    # make a list of project
    # load production data
    provenance_data = load_data(provenance_data_file)
    
    for case_data in provenance_data:
        # check that case data is complete (all sections in the case dictionary are complete)
        if is_case_info_complete(case_data):
            project_ids = [case_data['project_info'][i]['project'] for i in range(len(case_data['project_info']))]
            L.extend(project_ids)
            L = list(set(L))

    return L



def map_file_swids_to_cases(provenance_data_file):
    '''
    (str)- > dict
        
    Returns a dictionary mapping each file swid to a case id and a list of
    project ids if the case belongs to multiple projects 
    
    Parameters
    - provenance_data_file (str): Path to the provenance reporter json file
    '''
    
    D ={}
    
    # load production data
    provenance_data = load_data(provenance_data_file)
    for case_data in provenance_data:
        # check that case data is complete (all sections in the case dictionary are complete)
        if is_case_info_complete(case_data):
            case_id = case_data['case']
            # get projects
            project_ids = [case_data['project_info'][i]['project'] for i in range(len(case_data['project_info']))]
            project_ids = ';'.join(sorted(project_ids))
            for d in case_data['workflow_runs']:
                files = json.loads(d['files'])
                for i in files:
                    D[i['accession']] = {'case_id': case_id, 'project_id': project_ids}

    return D                     





def map_fileqc_to_cases(file_swids, fileqc):
    '''
    (dict, dict) -> list
        
    Returns a list of lists with with file QC info to be added to the Nabu cache
    
    Parameters
    ----------
    - file_swids (dict): Dictionary mapping file swids with case id and project ids
    - fileqc (dict): Dictionary with file qc information for each file swid
    '''

    data = []

    for fileid in fileqc:
        if fileid in file_swids:
            projects = file_swids[fileid]['project_id'].split(';')
            for project in projects:
                L = [project, fileid, file_swids[fileid]['case_id'], fileqc[fileid]['filepath'],
                     fileqc[fileid]['username'], fileqc[fileid]['qcstatus'], fileqc[fileid]['ticket']]
                data.append(L)   
        
    return data



def map_cases_to_projects(provenance_data_file):
    '''
    (str) -> dict
        
    Returns a dictionary mapping each case id to a list of projects if the case belongs
    to multiple projects
        
    Parameters
    ----------
    - provenance_data_file (str): Path to the provenance report json
    '''
    
    D = {}
    
    # load production data
    provenance_data = load_data(provenance_data_file)
    for case_data in provenance_data:
        # check that case data is complete (all sections in the case dictionary are complete)
        if is_case_info_complete(case_data):
            case_id = case_data['case']
            # get projects
            project_ids = [case_data['project_info'][i]['project'] for i in range(len(case_data['project_info']))]
            project_ids = ';'.join(sorted(project_ids))
            D[case_id] = project_ids
            
    return D                     



def get_release_signoff(case_signoffs, cases_to_projects):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with the release approval and release signoff for each case
    assing the projects to the case ids
    
    Parameters
    ----------
    - case_signoffs (dict): Dictionary with case signoff extracted from Nabu
    - cases_to_projects (dict): Dictionary mapping each case id to its project(s) 
    '''

    D = {}

    for case_id in case_signoffs:
        if case_id in cases_to_projects:
            projects = cases_to_projects[case_id]
            if case_id not in D:
                D[case_id] = {'projects': projects}
            if 'Release' in case_signoffs[case_id]:
                for d in case_signoffs[case_id]['Release']:
                    deliverableType = d['deliverableType']
                    deliverable = d['deliverable']
                    username = d['username']
                    comment = d['comment']
                    if comment:
                        ';'.join(comment)
                    qcpassed = d['qcPassed']
                    if 'release' not in D[case_id]:
                        D[case_id]['release'] = {}
                    if deliverableType not in D[case_id]['release']:
                        D[case_id]['release'][deliverableType] = {}
                    assert deliverable not in D[case_id]['release'][deliverableType]
                    D[case_id]['release'][deliverableType][deliverable] = {'username': username,
                                                                           'comment': comment,
                                                                           'qcpassed': qcpassed}
            if 'Release Approval' in case_signoffs[case_id]:
                for d in case_signoffs[case_id]['Release Approval']:
                    deliverableType = d['deliverableType']
                    username = d['username']
                    comment = d['comment']
                    if comment:
                        ';'.join(comment)
                    qcpassed = d['qcPassed']
                    if 'release_approval' not in D[case_id]:
                        D[case_id]['release_approval'] = {}
                    assert deliverableType not in D[case_id]['release_approval']
                    D[case_id]['release_approval'][deliverableType] = {'username': username,
                                                                       'comment': comment,
                                                                       'qcpassed': qcpassed}
            
    return D




def organize_release_signoff(release_signoffs):
    '''
    (dict) -> list
    
    Returns a list of lists with case signoff data to be added to the Nabu cache
    
    Parameters
    ----------
    - release_signofs (dict): Dictionary with projects, release approval and release signoffs of each case
    '''

    data = []

    for case_id in release_signoffs:
        projects = release_signoffs[case_id]['projects'].split(';')
        for project_id in projects:
            if 'release' in release_signoffs[case_id]:
                data_release = release_signoffs[case_id]['release']
            else:
                data_release ={}
            if 'release_approval' in release_signoffs[case_id]:
                release_approval = release_signoffs[case_id]['release_approval']
            else:
                release_approval ={}
            L = [case_id, json.dumps(data_release), json.dumps(release_approval), project_id]
            data.append(L)
             
    return data
                                      



def update_nabu_cache(provenance_data_file, nabu_cache, nabu_key_file, nabu = 'https://nabu.gsi.oicr.on.ca'):
    '''
    (str, str, str, str)- > None
    
    Builds and updates the nabu cache with file qc and case signoff information
        
    Parameters
    ----------
    - provenance_data_file (str): Path to the provenance report json
    - nabu_cache (str): Path to the nabu cache
    - nabu_key_file (str): Path to the file with Nabu access key
    - nabu (str): URL of the Nabu API ('https://nabu.gsi.oicr.on.ca')
    '''

    
    # create dict to store column names for each table {table: [column names]}
    fields = {'fileqc': {'names': ['project_id', 'case_id', 'fileid', 'filepath', 'username', 'qcstatus', 'ticket'],
                         'types': ['VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(572)', 'TEXT', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)']},
              'signoff': {'names': ['case_id', 'release', 'release_approval', 'project_id'],
                          'types': ['VARCHAR(572)', 'TEXT', 'TEXT', 'VARCHAR(128)']}}

    # create database if file doesn't exist
    if os.path.isfile(nabu_cache) == False:
        create_nabu_cache(nabu_cache, fields)
    
    # make a list of project
    projects = list_projects(provenance_data_file)
       
    # map file swids to case_id and project_id
    file_swids = map_file_swids_to_cases(provenance_data_file)
    nabu_fileqc_endpoint = nabu + '/get-fileqcs'
    
    for project in projects:
        fileqc = get_file_signoff(project, nabu_fileqc_endpoint)
        # organize data to add to cache
        data = map_fileqc_to_cases(file_swids, fileqc)
        # open database to delete old entries and add new ones
        conn = connect_to_db(nabu_cache)
        # remove entries for project
        delete_multiple_records([project], conn, nabu_cache, 'fileqc', 'project_id')
        # add entries for project
        insert_multiple_records(data, conn, nabu_cache, 'fileqc', fields['fileqc']['names'])
        conn.close()
        
      
    # get the case signoff
    nabu_signoff_endpoint = nabu + '/case/sign-off'
    case_signoffs = extract_nabu_signoff(nabu_key_file, nabu_signoff_endpoint)
    # map case ids to projects
    cases_to_projects = map_cases_to_projects(provenance_data_file)
    release_signoffs = get_release_signoff(case_signoffs, cases_to_projects)
    # organize data to add to cache
    signoff_data = organize_release_signoff(release_signoffs)
    # open database to delete old entries and add new ones
    conn = connect_to_db(nabu_cache)
    # remove entries for project
    delete_multiple_records(list(release_signoffs.keys()), conn, nabu_cache, 'signoff', 'case_id')
    # add entries for project
    insert_multiple_records(signoff_data, conn, nabu_cache, 'signoff', fields['signoff']['names'])
    conn.close()
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = 'build_nabu_cache.py', description='Script to build a cache with data extracted from Nabu')
    parser.add_argument('-pv', '--provenance_reporter', dest = 'provenance_data_file', help = 'Path to the provenance reporter json file', required=True)
    parser.add_argument('-db', '--database', dest = 'nabu_cache', help = 'Path to nabu cache', required=True)
    parser.add_argument('-nk', '--nabu_key_file', dest = 'nabu_key_file', help = 'Path to the Nabu key file', required=True)
    parser.add_argument('-nb', '--nabu', dest = 'nabu', default = 'https://nabu.gsi.oicr.on.ca',
                        help = 'URL of the Nabu API. Path to the Nabu key file. Default is https://nabu.gsi.oicr.on.ca')
    # get arguments from the command line
    args = parser.parse_args()
    #args.func(args)
    update_nabu_cache(args.provenance_data_file, args.nabu_cache, args.nabu_key_file, args.nabu)
     
    
