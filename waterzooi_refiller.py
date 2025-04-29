# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 11:48:39 2024

@author: rjovelin
"""



import sqlite3
import json
import argparse
import time
import os
import hashlib
# import matplotlib
# matplotlib.use('agg')
# import matplotlib.pyplot as plt
# import networkx as nx
# import numpy as np
# import base64
# import io


def connect_to_db(database):
    '''
    (str) -> sqlite3.Connection
    
    Returns a connection to SqLite database prov_report.db.
    This database contains information extracted from FPR
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    '''
    
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    return conn



def define_column_names():
    '''
    (None) -> dict

    Returns a dictionary with column names for each table in database
    '''

    # create dict to store column names for each table {table: [column names]}
    column_names = {'Workflows': ['wfrun_id', 'wf', 'wfv', 'case_id', 'project_id', 'donor_id', 'file_count', 'lane_count'],
                    'Parents': ['parents_id', 'children_id', 'project_id', 'case_id', 'donor_id'],
                    'Projects': ['project_id', 'pipeline', 'last_updated', 'cases', 'samples',
                                 'library_types', 'assays', 'deliverables', 'active'],
                    'Files': ['file_swid', 'project_id', 'md5sum', 'wfrun_id', 'file', 'attributes', 'creation_date', 'limskey', 'case_id', 'donor_id'],
                    'Libraries': ['library', 'lims_id', 'case_id', 'donor_id', 'tissue_type', 'tissue_origin',
                                  'library_type', 'group_id', 'group_id_description', 'project_id'],
                    'Workflow_Inputs': ['library', 'run', 'lane', 'wfrun_id', 'limskey', 'barcode', 'platform', 'project_id', 'case_id', 'donor_id'],
                    'Samples': ['case_id', 'donor_id', 'ext_id', 'species', 'sex', 'miso', 'project_id'],
                    'Checksums': ['project_id', 'case_id', 'donor_id', 'md5']
                    }
        
    return column_names


def define_column_types():
    '''
    (None) -> dict

    Returns a dictionary with column types for each table in database
    '''
    
    # create dict to store column names for each table {table: [column names]}
    column_types = {'Workflows': ['VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)', 'INT', 'INT'],
                    'Parents': ['VARCHAR(572)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)'],
                    'Projects': ['VARCHAR(128) PRIMARY KEY NOT NULL UNIQUE', 'VARCHAR(128)',
                                 'VARCHAR(256)', 'INT', 'INT', 'VARCHAR(256)', 'VARCHAR(572)',
                                 'VARCHAR(572)', 'INT'],
                    'Files': ['VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(256)',
                              'VARCHAR(572)', 'TEXT', 'TEXT', 'INT', 'VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(128)'],
                    'Libraries': ['VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(572)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(128)'],
                    'Workflow_Inputs': ['VARCHAR(128)', 'VARCHAR(256)', 'INTEGER', 'VARCHAR(572)', 
                                        'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)'],
                    'Samples': ['VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)'],
                    'Checksums': ['VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(572)']
                    }
    
    return column_types



def create_table(database, table):
    '''
    (str, str) -> None
    
    Creates a table in database
    
    Parameters
    ----------
    - database (str): Name of the database
    - table (str): Table name
    '''
    
    # get the column names
    column_names = define_column_names()[table]
    # get the column types
    column_types = define_column_types()[table]    
    
    # define table format including constraints    
    table_format = ', '.join(list(map(lambda x: ' '.join(x), list(zip(column_names, column_types)))))

    if table  in ['Workflows', 'Parents', 'Files', 'Libraries', 'Workflow_Inputs', 'Samples', 'Checksums']:
        constraints = '''FOREIGN KEY (project_id)
            REFERENCES Projects (project_id)'''
        table_format = table_format + ', ' + constraints 
    
    if table == 'Parents':
        constraints = '''FOREIGN KEY (parents_id)
          REFERENCES Workflows (wfrun_id),
          FOREIGN KEY (children_id)
              REFERENCES Workflows (wfrun_id)''' 
        table_format = table_format + ', ' + constraints + ', PRIMARY KEY (parents_id, children_id, project_id, case_id)'
    
    if table == 'Worklows':
        table_format = table_format + ', PRIMARY KEY (wfrun_id, project_id)'
    
    if table == 'Files':
        constraints = '''FOREIGN KEY (wfrun_id)
            REFERENCES Workflows (wfrun_id)'''
        table_format = table_format + ', ' + constraints
    
    if table == 'Workflow_Inputs':
        constraints = '''FOREIGN KEY (wfrun_id)
            REFERENCES Workflows (wfrun_id),
            FOREIGN KEY (library)
              REFERENCES Libraries (library)'''
        table_format = table_format + ', ' + constraints
    
    if table == 'Samples':
        constraints = '''FOREIGN KEY (ext_id)
            REFERENCES Libraries (ext_id)'''
        table_format = table_format + ', ' + constraints

    if table == 'Libraries':
        constraints = '''FOREIGN KEY (case_id)
            REFERENCES Samples (case_id)'''
        table_format = table_format + ', ' + constraints

    # connect to database
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    # create table
    cmd = 'CREATE TABLE {0} ({1})'.format(table, table_format)
    cur.execute(cmd)
    conn.commit()
    conn.close()


def initiate_db(database):
    '''
    (str) -> None
    
    Create tables in database
    
    Parameters
    ----------
    - database (str): Path to the database file
    '''
    
    # check if table exists
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cur.fetchall()
    tables = [i[0] for i in tables]    
    conn.close()
    for i in ['Projects', 'Workflows', 'Parents', 'Files', 'Libraries',
              'Workflow_Inputs', 'Samples', 'Calculate_Contamination', 'Checksums']:
        if i not in tables:
            create_table(database, i)


def insert_data(database, table, data, column_names):
    '''
    (str, str, list, list) -> None
    
    Inserts data into the database table with column names 
    
    Parameters
    ----------
    - database (str): Path to the database file
    - table (str): Table in database
    - data (list): List of data to be inserted
    - column_names (list): List of table column names
    '''
    
    if data:
        print('inserting data in {0}'.format(table))
        # connect to db
        conn = sqlite3.connect(database)
        # add data
        vals = '(' + ','.join(['?'] * len(data[0])) + ')'
        conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), data)
        conn.commit()
        conn.close()


def delete_records(cases, database, table):
    '''
    (dict, str, str) -> None
    
    Remove all the rows from table with case_id in cases
    
    Parameters
    ----------
    - donors (dict): Dictionary with cases to remove from table
    - database (str): Path to the sqlite database
    - table (str): Table in database
    '''
    
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    id_list = list(cases.keys())
    query = "DELETE FROM {0} WHERE case_id IN ({1})".format(table, ", ".join("?" * len(id_list)))
    cur.execute(query, id_list)
    conn.commit()
    conn.close()


def compute_md5(d):
    '''
    (dict) -> str
    
    Returns the md5 checksum of a dictionary d
    
    Parameters
    ----------
    d (dict): Dictionary with information parsed from FPR
    '''
    
    return hashlib.md5(json.dumps(d, sort_keys=True).encode('utf-8')).hexdigest()


def compute_case_md5sum(provenance_data):
    '''
    (list) -> dict
    
    Returns a dictionary with the md5sum of the each case's data listed in
    provenance_data
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries, each representing the data of a single case
    '''

    D = {}
    for d in provenance_data:
        md5 = compute_md5(d)
        case = d['case']
        donor = get_donor_name(d)
        project = d['project']
        assert case not in D
        D[case] = {'md5':md5, 'project_id':project, 'donor_id':donor}
    return D
    

def get_cases_md5sum(database, table = 'Checksums'):
    '''
    (str, str) -> dict

    Returns a dictionary of cases and checksum extracted from the table Checksums in the database
           
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - table (str): Table storing the donor checksum information. Default is Checksums 
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
    

def cases_info_to_update(md5sums, recorded_md5sums):
    '''
    (dict, dict) -> dict

    Returns a dictionary of cases, checksum for which the information in the database needs to be updated
    (ie, the checksum in the database is different from the checksum of the production data,
     or the case recorded in the database is no longer in production)
            
    Parameters
    ----------
    - md5sums (dict): Dictionary of cases, checksum for production data
    - recorded_md5sums (dict): Dictionary of recorded cases, checksum in the database
    '''
        
    cases = {}
    for case in md5sums:
        # update if not already recorded
        if case not in recorded_md5sums:
            cases[case] = {'md5': md5sums[case]['md5'], 'project_id': md5sums[case]['project_id'], 'donor_id': md5sums[case]['donor_id']}
        # update if md5sums are different
        else:
            if recorded_md5sums[case]['md5'] != md5sums[case]['md5']:
                cases[case] = md5sums[case]
            
    # delete cases that are no longer recorded
    for case in recorded_md5sums:
        if case not in md5sums:
            cases[case] = 'delete'
        
    return cases


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


def add_checksums_info_to_db(database, cases_to_update, table = 'Checksums'):
    '''
    (str, dict, str) -> None
    
    Update table Checksums with the checksum of the case info 
       
    Parameters
    ----------
    - database (str): Path to the database file
    - cases_to_update (dict): Dictionary with cases for which records needs to be updated
    - table (str): Name of Table in database. Default is Checksums
    '''
    
    if cases_to_update:
        delete_records(cases_to_update, database, table)
        
        # make a list of data to insert
        newdata = []
        
        # connect to db
        conn = sqlite3.connect(database, timeout=30)
        # get column names
        data = conn.execute("SELECT * FROM {0}".format(table))
        column_names = [column[0] for column in data.description]

        # order values according to column names
        for i in cases_to_update:
            L = [cases_to_update[i]['project_id'], i, cases_to_update[i]['donor_id'],  cases_to_update[i]['md5']]
            newdata.append(L)
        vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
        conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
        conn.commit()
        conn.close()



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



# def remove_data_from_inactive_projects(provenance_data):
#     '''
#     (list) -> list

#     Returns the list of case data removing any data from inactive project
    
#     Parameters
#     ----------
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     '''
    
#     to_remove = [i for i in provenance_data if is_project_active(i) == False]
#     if to_remove:
#         for i in to_remove:
#             provenance_data.remove(i)
#     return provenance_data


# def is_project_active(case_data):
#     '''
#     (dict) -> bool
    
#     Returns True if the project is active and False otherwise
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single donor data    
#     '''
    
#     active = [convert_to_bool(i['active']) for i in case_data['pinery_project_data']]
#     return all(active)



# def remove_cases_without_cerberus_data(provenance_data):
#     '''
#     (list) -> list
    
#     Returns the list of case data removing any cases that do not have data from cerberus
    
#     Parameters
#     ----------
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     '''

#     to_remove = [i for i in provenance_data if len(i['cerberus_data']) == 0]
#     if to_remove:
#         for i in to_remove:
#             provenance_data.remove(i)
#     return provenance_data



# def remove_cases_without_pinery_data(provenance_data):
#     '''
#     (list) -> list
    
#     Returns the list of case data removing any cases that do not have project
#     information from pinery
    
#     Parameters
#     ----------
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     '''

#     to_remove = [i for i in provenance_data if len(i['pinery_project_data']) == 0]
#     if to_remove:
#         for i in to_remove:
#             provenance_data.remove(i)
#     return provenance_data




def get_donor_name(case_data):
    '''
    (str) -> str
    
    Returns the name of the donor in case_data
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data 
    '''

    donor = list(set([i['donor'] for i in case_data['pineryData']]))
    assert len(donor) == 1
    donor = donor[0]

    return donor


def get_file_timestamp(d):
    '''
    (dict) - > int
    
    Returns the time stamp of a file from the dictionary of the file from cerberus data
    for a donor
    
    Parameters
    ----------
    - d (dict): Dictionary representing a file information from cerberus
    '''    
    
    creation_date = d['timestamp']
    creation_date = ' '.join(creation_date.split('T')).replace('Z', '')
    creation_date = creation_date.split('.')
    if len(creation_date) > 1:
        creation_date = ' '.join(creation_date[:-1])
    else:
        creation_date = creation_date[0]
    pattern = '%Y-%m-%d %H:%M:%S'
    creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
    
    return creation_date



def collect_case_file_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the file information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''    
        
    D = {}
    
    for d in case_data['workflowRuns']:
        donor = get_donor_name(case_data)
        case = case_data['case']
        project_id = case_data['project']
        wfrun_id = d['wfrunid']
        limskeys = d['limsIds'].split(',')
        file_swids = list(map(lambda x: x.strip(), json.loads(d['files'])))   
        
        ## olive must be updated to collect these fields
        md5sums = ['NA'] * len(file_swids)
        creation_date = ['NA'] * len(file_swids)
        attributes = ['NA'] * len(file_swids)
        file_paths = ['NA'] * len(file_swids)
         
        # if attributes:
        #     attributes = json.loads(d['attributes'])
        #     for k in attributes:
        #         attributes[k] = attributes[k][0]
        #     if len(attributes) == 0:
        #         attributes = ''
        #     else:
        #         attributes = json.dumps(attributes)
        # else:
        #     attributes = ''
    
    
        for i in range(len(file_swids)):
            assert file_swids[i] not in D
            D[file_swids[i]] = {'file_swid': file_swids[i],
                                'md5sum': md5sums[i],
                                'creation_date': creation_date[i],
                                'attributes': attributes[i],
                                'file': file_paths[i],
                                'project_id': project_id,
                                'case_id': case,
                                'donor_id': donor,
                                'wfrun_id': wfrun_id,
                                'limskey': limskeys}  
                                
    return D


def add_file_info_to_db(database, provenance_data, cases_to_update, table = 'Files'):
    '''
    (str, list, dict, str) -> None
    
    Inserts file information in database's Files table
       
    Parameters
    ----------
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - cases_to_update (dict): Dictionary with cases for which records needs to be updated
    - table (str): Table in database storing file information. Default is Files
    '''
        
    if cases_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(cases_to_update, database, table)
        print('deleted records in {0}'.format(table))
    
        # make a list of data to insert in the database
        newdata = []
            
        for case_data in provenance_data:
            donor = get_donor_name(case_data)
            case =  case_data['case']
            # check if donor needs to be updated
            if case in cases_to_update and cases_to_update[case] != 'delete':
                file_info = collect_case_file_info(case_data)
                for file_swid in file_info:
                    for limskey in file_info[file_swid]['limskey']:
                        L = [file_info[file_swid][i] for i in column_names[:column_names.index('limskey')]]
                        L.append(limskey)
                        L.extend([file_info[file_swid][i] for i in column_names[column_names.index('limskey')+1:]])
                        newdata.append(L)             
    
        # add data
        insert_data(database, table, newdata, column_names)


def collect_case_library_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the library information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''    

    D = {}
    
    for d in case_data['pineryData']:
        donor = d['donor']
        library = d['library']
        case = case_data['case']
        tissue_origin = d['tissueOrigin']
        tissue_type = d['tissueType']
        library_type = d['libraryDesign']
        project = d['project']
        group_id = d['groupId']
        group_id_description = d['groupDesc']
        lims_id = d['limsId']
        sample = d['sampleId']
               
        d = {'library': library, 'case_id': case, 'donor_id': donor, 
             'project_id': project, 'group_id': group_id, 
             'group_id_description': group_id_description,
             'library_type': library_type, 'tissue_type': tissue_type,
             'tissue_origin': tissue_origin, 'lims_id': lims_id, 'sample': sample}
        
        if library not in D:
            D[library] = [d] 
        else:
            if d not in D[library]:
                D[library].append(d)
                   
    return D       


def add_library_info_to_db(database, provenance_data, cases_to_update, table = 'Libraries'):
    '''
    (str, list, dict, str) -> None
    
    Inserts library information in the Libraries table of the database
       
    Parameters
    ----------
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - cases_to_update (dict): Dictionary with cases for which records needs to be updated
    - table (str): Table in database storing file information. Default is Libraries
    '''                 
 
    if cases_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for cases to update
        delete_records(cases_to_update, database, table)
        print('deleted records in {0}'.format(table))

        # make a list of data to insert in the database
        newdata = []
    
        for case_data in provenance_data:
            case =  case_data['case']
            # check if donor needs to be updated
            if case in cases_to_update and cases_to_update[case] != 'delete':
                library_info = collect_case_library_info(case_data)
                for library in library_info:
                    for d in library_info[library]:
                        L = [d[i] for i in column_names]
                        newdata.append(L)             
  
        # add data
        insert_data(database, table, newdata, column_names)




def get_pipeline(case_data):
    '''
    (dict) -> str
    
    Returns the pipeline the project extracted from the pinery project data of a case 
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data    
    '''

    L = list(set([i['pipeline'] for i in case_data['pinery_project_data']]))
    assert len(L) == 1
    pipeline = L[0]    
    return pipeline



# def collect_project_info_from_db(database, table = 'Libraries'):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with project level information for each active project in provenance_data
    
#     Parameters
#     ----------
#     - database (str): Path to the database
#     - table (str): Table storing the library information
#     '''

#     D = {}

#     # connect to db
#     conn = connect_to_db(database)
#     # get column names
#     data = conn.execute('SELECT * FROM {0}'.format(table))
#     data = data.fetchall()
#     conn.close()
    
#     for i in data:
#         project = i['project_id']
#         library_design = i['library_type']
#         sample = '_'.join([i['donor_id'], i['tissue_type'], i['tissue_origin'],
#                            library_design, i['group_id']])
        
#         if project not in D:
#             D[project] = {}
#         D[project][sample] = library_design
                               
#     return D


def collect_project_info(provenance_data):
    '''
    (list) -> dict
    
    Returns a dictionary with project level information for each project in provenance_data
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries with case information 
    '''

    D = {}
    
    for case_data in provenance_data:
        project = case_data['project']
        case = case_data['case']
        assay = case_data['assay']
        assert len(case_data['projectInfo']) == 1
        active = case_data['projectInfo'][0]['isActive']
        pipeline = case_data['projectInfo'][0]['pipeline']
        deliverables = case_data['projectInfo'][0]['deliverables']
        samples, library_types = [], []
        for d in case_data['pineryData']:
            samples.append(d['sampleId'])
            library_types.append(d['libraryDesign'])
        if project not in D:
            D[project] = {'assays': [assay],
                          'cases': [case],
                          'active': active,
                          'pipeline': pipeline,
                          'deliverables': deliverables,
                          'samples': samples,
                          'library_types': library_types}
        else:
            D[project]['assays'].append(assay)
            D[project]['cases'].append(case)
            D[project]['samples'].extend(samples)
            D[project]['library_types'].extend(library_types)
                                
        D[project]['assays'] = list(set(D[project]['assays']))
        D[project]['cases'] = list(set(D[project]['cases']))
        D[project]['samples'] = list(set(D[project]['samples']))
        D[project]['library_types'] = list(set(D[project]['library_types']))
            
    return D



# def count_project_case(provenance_data):
#     '''
#     (list) -> dict
    
#     Returns a dictionary of case count per project
    
#     Parameters
#     ----------
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     '''
    
#     D = {}
    
#     for case_data in provenance_data:
#         project = get_project_name(case_data)
#         case = get_case_name(case_data)
#         if project in D:
#             D[project].append(case)
#         else:
#             D[project] = [case]
            
#     counts = {}
#     for project in D:
#         D[project] = set(D[project])
#         counts[project] = len(D[project])

#     return counts


def add_project_info_to_db(database, provenance_data, project_table = 'Projects'):
    '''
    (str, list, str) -> None
    
    Add project information into Projects table of database
       
    Parameters
    ----------
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single case
    - project_table (str): Name of Table storing project information in database. Default is Projects
    '''
    
    # get the column names
    column_names = define_column_names()[project_table]
    
    # collect project information
    project_info = collect_project_info(provenance_data)
    
    # make a list of data to insert in the database
    newdata = []
    
    # connect to db
    conn = connect_to_db(database)
       
    for project in project_info:
        last_updated = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
        
        # get the library types
        library_types = ','.join(sorted(project_info[project]['library_types']))
        samples = len(project_info[project]['samples'])
        cases = len(project_info[project]['cases'])
        assays = ','.join(sorted(project_info[project]['assays']))
        
        L = [project, project_info[project]['pipeline'], last_updated, str(cases),
             str(samples), str(library_types), assays, project_info[project]['deliverables'],
             project_info[project]['active']]
        
        newdata.append(L)
        # delete project entry
        query = 'DELETE FROM {0} WHERE project_id = \"{1}\"'.format(project_table, project)
        conn.execute(query)
        conn.commit()

    conn.close()
        
    # add data
    insert_data(database, project_table, newdata, column_names)



def collect_workflow_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the workflow information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''

    D = {}
    
    for d in case_data['workflowRuns']:
        project_id = case_data['project']
        case = case_data['case']
        donor = get_donor_name(case_data)
        wfrun_id = d['wfrunid']
        wf = d['wf']
        wfv = d['wfv']
        limskeys = d['limsIds'].split(',')
        file_count = len(d['files'])
        sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
        lane_count = len([sequencing_attributes[i]['lane'] for i in sequencing_attributes])
         
        data = {'project_id': project_id,
                'wfrun_id': wfrun_id,
                'wf': wf,
                'wfv': wfv,
                'case_id': case,
                'donor_id': donor,
                'file_count': file_count,
                'lane_count': lane_count}
                
        assert wfrun_id not in D
        D[wfrun_id] = data
                
    return D            



def add_workflows_to_db(database, provenance_data, cases_to_update, table = 'Workflows'):
    '''
    (str, list, dict, str) -> None
    
    Inserts or updates workflow information 
 
    Parameters
    ----------    
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - cases_to_update (dict): Dictionary with cases for which records needs to be updated
    - table (str): Table in database storing file information. Default is Workflows
    '''
    
    if cases_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(cases_to_update, database, table)
        print('deleted records in {0}'.format(table))

        # make a list of data to insert in the database
        newdata = []
        
        for case_data in provenance_data:
            case = case_data['case']
            donor = get_donor_name(case_data)
                      
            # check if donor needs to be updated
            if case in cases_to_update and cases_to_update[case] != 'delete':
                workflow_info = collect_workflow_info(case_data)                
                for workflow in workflow_info:
                    L = [workflow_info[workflow][i] for i in column_names]
                    newdata.append(L)             
      
        # add data
        insert_data(database, table, newdata, column_names)



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
        for d in case_data['pineryData']:
            if i == d['limsId']:
                assert i not in D
                D[i] = {'library': d['library'],
                        'barcode': d['barcode'],
                        'lane': d['lane'],
                        'project': d['project'],
                        'run': d['run']}
    
    return D


def collect_case_workflow_inputs(case_data):
    '''
    (dict) -> list
    
    Returns a list of dictionaries with the workflow input information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''
        
    L = []
    
    for d in case_data['workflowRuns']:
        case = case_data['case']
        donor = get_donor_name(case_data)
        limskeys = d['limsIds'].split(',')
        sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
        platform = 'NA'
        wfrun_id = d['wfrunid']
                
        for limskey in sequencing_attributes:
            D = {'case_id': case,
                 'donor_id': donor,
                 'library': sequencing_attributes[limskey]['library'],
                 'project_id': sequencing_attributes[limskey]['project'],
                 'barcode': sequencing_attributes[limskey]['barcode'],
                 'platform': platform,
                 'lane': sequencing_attributes[limskey]['lane'],
                 'wfrun_id': wfrun_id,
                 'run': sequencing_attributes[limskey]['run'],
                 'limskey': limskey}
            
            if D not in L:
                L.append(D)
    
    return L



def add_workflow_inputs_to_db(database, provenance_data, cases_to_update, table = 'Workflow_Inputs'):
    '''
    (str, list, dict, str) -> None
    
    Inserts or updates workflow input information 
 
    Parameters
    ----------    
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single donor
    - cases_to_update (dict): Dictionary with donors for which records needs to be updated
    - table (str): Table in database storing file information. Default is Workflow_Inputs
    '''

    if cases_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(cases_to_update, database, table)
        print('deleted records in {0}'.format(table))

        # make a list of data to insert in the database
        newdata = []
     
        for case_data in provenance_data:
            case = case_data['case']
            # check if donor needs to be updated
            if case in cases_to_update and cases_to_update[case] != 'delete':
                workflow_input_info = collect_case_workflow_inputs(case_data)
                for d in workflow_input_info:
                    L = [d[i] for i in column_names]
                    newdata.append(L)             

        # add data
        insert_data(database, table, newdata, column_names)


def map_file_to_worklow(case_data):
    '''
    (dict) -> dict

    Returns a dictionary of file swids matched to their workflow run id for a case    
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''
    
    D = {}
       
    for i in range(len(case_data['cerberus_data'])):
        wfrun_id = case_data['cerberus_data'][i]['workflow_run_accession']
        file_swid = case_data['cerberus_data'][i]['accession'] 
        if file_swid in D:
            assert D[file_swid] == wfrun_id
        else:
            D[file_swid] = wfrun_id
    
    return D



def get_workflow_inputs(case_data):
    '''
    (dict) -> dict

    Returns a dictionary of workflows and their input files for a case    
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''
    
    D = {}
       
    for i in range(len(case_data['cerberus_data'])):
        wfrun_id = case_data['cerberus_data'][i]['workflow_run_accession']
        input_files = json.loads(case_data['cerberus_data'][i]['input_files']) 
        if wfrun_id in D:
            assert D[wfrun_id] == input_files
        else:
            D[wfrun_id] = input_files
    
    return D


def identify_parent_children_workflows(workflow_inputs, files):
    '''
    (dict, dict) -> dict     
    
    Returns a dictionary of children: parents workflows relationsips for a case
        
    Parameters
    ----------
    - workflow_inouts (dict): Dictionary of workflows and their input files for a case
    - files (dict): Dictionary of file swids matched to their workflow run id for a case    
    '''
    
    # parents record child-parent workflow relationships
    D = {}
    
    for workflow in workflow_inputs:
        if workflow_inputs[workflow]:
            parent_workflows = sorted(list(set([files[i] for i in workflow_inputs[workflow] if i in files])))
        else:
            parent_workflows = ['NA']
        if workflow not in D:
            D[workflow] = parent_workflows
        else:
            assert D[workflow] == parent_workflows
    
    return D



def collect_workflow_relationships(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary of parent-children workflows
    for all workflows for a given case

    Paramaters
    ----------
    - case_data (dict): Dictionary with a single case data         
    '''
    
    D = {}
    
    for d in case_data['workflowRuns']:
        workflow = d['wfrunid']
        parents = json.loads(d['parents'])
        children = json.loads(d['children'])
        if parents:
            parent_workflows = [i[1] for i in parents]
        else:
            parent_workflows = ['NA']
        if children:
            children_workflows = [i[1] for i in children]
        else:
            children_workflows = ['NA']
        
        # record the parent-child relationship of the current workflow
        if workflow in D:
            D[workflow].extend(children_workflows)
        else:
            D[workflow] = children_workflows
        D[workflow] = list(set(D[workflow]))
        # record the parent-child relationships of each parent and current workflow 
        for workflow_run in parent_workflows:
            if workflow_run in D:
                D[workflow_run].append(workflow)
            else:
                D[workflow_run] = [workflow]
            D[workflow_run] = list(set(D[workflow_run]))

    return D



def add_workflows_relationships_to_db(database, provenance_data, cases_to_update, table = 'Parents'):
    '''
    (str, str, str, dict, str, str, str, str) -> None
    
    Inserts or updates workflow information and parent-children workflow relationships
 
    Parameters
    ----------   
    - database (str): Path to the database file
    - provenance_data (list): List of dictionaries, each representing the data of a single case
    - cases_to_update (dict): Dictionary with cases for which records needs to be updated
    - table (str): Table in database storing file information. Default is Parents
    '''

    if cases_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(cases_to_update, database, table)
        print('deleted records in {0}'.format(table))
        
        # make a list of data to insert in the database
        newdata = []
        
        for case_data in provenance_data:
            case = case_data['case']
            donor = get_donor_name(case_data)
            project = case_data['project']             
                       
            # check if donor needs to be updated
            if case in cases_to_update and cases_to_update[case] != 'delete':
                parents =  collect_workflow_relationships(case_data)
                for workflow in parents:
                    for child in parents[workflow]:
                        L = [workflow, child, project, case, donor]
                        if L not in newdata:
                            newdata.append(L)
        
        # add data
        insert_data(database, table, newdata, column_names)


def get_donor_external_id(case_data):
    '''
    (dict) -> str
    
    Returns the external id of a donor by extracting the external id from the
    cerberus data of the corresponding case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data    
    '''

    L = []

    for i in range(len(case_data['cerberus_data'])):
        L.append(case_data['cerberus_data'][i]['external_donor_id'])
    L = list(set(L))
    assert len(L) == 1
    return L[0]


def get_donor_sex(case_data):
    '''
    (dict) -> str
    
    Returns the sex of a donor by extracting the sex from the
    pinery data of a case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data    
    '''
    
    sex = []
    for i in range(len(case_data['pinery_data'])):
        sex.append(case_data['pinery_data'][i]['sex'])       
    sex = map(lambda x: x.lower(), list(set(sex)))
    
    # sex may not be defined for all samples.
    # assign donor sex to a known sample sex
    if 'female' in sex:
        assert 'male' not in sex
        sex = 'Female'
    elif 'male' in sex:
        assert 'female' not in sex
        sex = 'Male'
    else:
        sex = 'Unknown'
        
    return sex


def get_donor_species(case_data):
    '''
    (dict) -> str
    
    Returns the species of a donor by extracting the species from the
    pinery data of a case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data    
    '''
    
    species = []
    for i in range(len(case_data['pinery_data'])):
        species.append(case_data['pinery_data'][i]['organism'])       
    species = list(set(species))
    assert len(species) == 1
        
    return species[0]


def collect_case_sample_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with the sample information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''
        
    D = {}
        
    for d in case_data['pineryData']:
        case = case_data['case']
        donor = d['donor']
        project_id = d['project']
        
        ext_id = 'NA'
        sex = 'NA'
        miso = 'NA'
        species = 'NA'
        
        d = {'case_id': case, 'donor_id': donor, 'ext_id': ext_id, 'sex': sex,
             'species': species, 'project_id': project_id, 'miso': miso}
        
        if case not in D:
            D[case] = d
        else:
            assert D[case] == d
            
    return D



def add_samples_info_to_db(database, provenance_data, cases_to_update, table = 'Samples'):
    '''
    (str, str, str, dict, dict, dict) -> None
    
    Inserts samples data into Samples table of database    
    
    Parameters
    ----------
    - database (str): Path to the databae file
    - provenance_data (list): List of dictionaries, each representing the data of a single case
    - cases_to_update (dict): Dictionary with cases for which records needs to be updated
    - table (str): Name of table in database
    '''
    
    if cases_to_update:
        # get the column names
        column_names = define_column_names()[table]
        # remove rows for donors to update
        delete_records(cases_to_update, database, table)
        print('deleted records in {0}'.format(table))

        # make a list of data to insert in the database
        newdata = []
        
        for case_data in provenance_data:
            case = case_data['case']
            # check if donor needs to be updated
            if case in cases_to_update and cases_to_update[case] != 'delete':
                sample_info = collect_case_sample_info(case_data)
                L = [sample_info[case][i] for i in column_names]
                newdata.append(L)
                
        # add data
        insert_data(database, table, newdata, column_names) 









########################################################




# def add_contamination_info(database, calcontaqc_db, provenance_data, donors_to_update, table = 'Calculate_Contamination'):
#     '''
#     (str, str, dict, str) -> None
    
#     Parse the calcontaqc_db, reformat data and add information to the Calculate_Contamination
#     table of the database
    
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - calcontaqc_db (str): Path to the calculate contamination database
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Parents
#     '''

#     if donors_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(donors_to_update, database, table)
#         print('deleted records in {0}'.format(table))
        
#         # make a list of data to insert in the database
#         newdata = []
        
#         for donor_data in provenance_data:
#             donor = donor_data['donor']
#             # check if donor needs to be updated
#             if donor in donors_to_update and donors_to_update[donor] != 'delete':
#                 contamination = collect_donor_calculate_contamination_db(calcontaqc_db, donor)
#                 if contamination:
#                     for sample in contamination[donor]:
#                         for i in contamination[donor][sample]:
#                             L = [i[j] for j in column_names]
#                             newdata.append(L)
                   
#         # add data
#         insert_data(database, table, newdata, column_names)
        

# def collect_donor_calculate_contamination_db(calcontaqc_db, donor):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with information collected from the calculate contamination
#     QC-etl sqlite cache for a donor of interest
    
#     Parameters
#     ----------
#     - calcontaqc_db (str): Path to the calculate contamination database
#     - donor (str): Donor of interest
#     '''

#     # get tables in db
#     conn = sqlite3.connect(calcontaqc_db)
#     cur = conn.cursor()
#     cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
#     tables = cur.fetchall()
#     tables = [i[0] for i in tables]    
#     conn.close()
#     # connect to db and retrieve all data
#     conn = connect_to_db(calcontaqc_db)
#     data = conn.execute('SELECT * FROM {0} where Donor == \"{1}\"'.format(tables[0], donor)).fetchall()
#     conn.close()

#     D = {}
    
#     for i in data:
#         case_id = i['Donor']
#         group_id = i['Group ID']
#         library_type = i['Library Design']
#         tissue_origin = i['Tissue Origin']
#         tissue_type = i['Tissue Type']
#         contamination = i['contamination']
#         merged_limskey = i['Merged Pinery Lims ID']
#         merged_limskey = list(map(lambda x: x.strip(), merged_limskey.replace('[', '').replace(']', '').replace('\"', '').split(',')))
#         merged_limskey = ';'.join(sorted(merged_limskey))
#         sample_id = '_'.join([i['Donor'], i['Tissue Type'], i['Tissue Origin'],
#                               i['Library Design'], i['Group ID']])
        
#         if case_id not in D:
#             D[case_id] = {}
#         if sample_id not in D[case_id]:
#             D[case_id][sample_id] = []
#         d = {'case_id': case_id,
#               'group_id': group_id,
#               'library_type': library_type, 
#               'tissue_origin': tissue_origin,
#               'tissue_type': tissue_type,
#               'contamination': contamination,
#               'merged_limskey': merged_limskey,
#               'sample_id': sample_id
#               }
#         D[case_id][sample_id].append(d)
           
#     return D                         


# def add_WGS_blocks_to_db(database, provenance_data, donors_to_update, table,
#                          expected_workflows, qc_workflows, library_type = 'WG', platform = 'novaseq'):
#     '''
#     (str, list, dict, str, list, list, str, str) -> None
    
#     Inserts WGS blocks data into WGS_blocks table of database    
    
#     Parameters
#     ----------
#     - database (str): Path to the databae file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): Name of table in database
#     - expected_workflows (list): List of expected workflow names to define a complete block
#     - qc_workflows (tuple): List of QC workflows to exclude from analysis blocks
#     - library_type (str): Library design of the data forming the blocks.
#                           Accepted values are WG and EX
#     - platform (str): Sequencing platform. Default is novaseq
#     '''
    
#     if donors_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(donors_to_update, database, table)
#         print('deleted records in {0}'.format(table))
        
#         # make a list of data to insert in the database
#         newdata = []
    
#         for donor_data in provenance_data:
#             donor = donor_data['donor']
#             # check if donor needs to be updated
#             if donor in donors_to_update and donors_to_update[donor] != 'delete':
#                 blocks = find_donor_WGS_blocks(donor_data, library_type, platform, expected_workflows, qc_workflows)
#                 for samples in blocks:
#                     for block in blocks[samples]:
#                         L = [blocks[samples][block]['project_id'],
#                              blocks[samples][block]['case_id'],
#                              samples,
#                              '.'.join(list(map(lambda x: os.path.basename(x), block.split('.')))),
#                              ';'.join(list(map(lambda x: os.path.basename(x), blocks[samples][block]['workflows']))),
#                               blocks[samples][block]['name'],
#                               blocks[samples][block]['date'],
#                               blocks[samples][block]['complete'],
#                               blocks[samples][block]['clean'],
#                               blocks[samples][block]['network']]
#                         newdata.append(L)
        
#         # add data
#         insert_data(database, table, newdata, column_names)



# def check_workflow_environment(blocks):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with the environment of each anchor workflow in blocks
    
#     Parameters
#     ----------
#     - blocks (dict): Dictionary with analysis workflows grouped by sample pairs and blocks 
#     '''
    
#     D = {}
#     for sample in blocks:
#         for anchor in blocks[sample]:
#             workflows = anchor.split('.')
#             for workflow in workflows:
#                 workflow_id = os.path.basename(workflow)
#                 environment = os.path.dirname(workflow)
#                 if workflow_id in D:
#                     D[workflow_id].append(environment)
#                 else:
#                     D[workflow_id] = [environment]
#                 D[workflow_id] = list(set(D[workflow_id]))
#     return D


# def is_anchor_stale(blocks):
#     '''
#     (dict) -> bool
    
#     Returns True if any of the anchor workflows exists in multiple environments (ie, indicating stale records)
    
#     Parameters
#     ----------
#     - blocks (dict): Dictionary with analysis workflows grouped by sample pairs and blocks 
#     '''
    
#     stale = False
    
#     # get the environment of each anchor workflow
#     environments = check_workflow_environment(blocks)
#     # check if the anchor workflows have multiple environments
#     for workflow in environments:
#         if len(environments[workflow]) > 1:
#             stale = True
    
#     return stale
    

# def add_WT_blocks_to_db(database, provenance_data, donors_to_update, table,
#                         expected_workflows, qc_workflows, library_type = 'WT', platform = 'novaseq'):
#     '''
#     (str, str, str, str) -> None
    
#     Inserts WGS blocks data into WGS_blocks table of database    
    
#     Parameters
#     ----------
#     - database (str): Path to the databae file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - donors_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): Name of table in database
#     - expected_workflows (list): List of expected workflow names to define a complete block
#     - qc_workflows (tuple): List of QC workflows to exclude from analysis blocks
#     - library_type (str): Library design of the data forming the blocks. Default is WT
#     - platform (str): Sequencing platform. Default is novaseq
#     '''
    
#     if donors_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(donors_to_update, database, table)
#         print('deleted records in {0}'.format(table))
        
#         # make a list of data to insert in the database
#         newdata = []
    
#         for donor_data in provenance_data:
#             donor = donor_data['donor']
            
#             # check if donor needs to be updated
#             if donor in donors_to_update and donors_to_update[donor] != 'delete':
#                 blocks = find_donor_WT_blocks(donor_data, library_type, platform, expected_workflows, qc_workflows)
#                 # proceed if the anchor workflows have a unique environment (ie, no stale record)
#                 if is_anchor_stale(blocks) == False:
#                     for samples in blocks:
#                         for block in blocks[samples]:
#                             L = [blocks[samples][block]['project_id'],
#                                 blocks[samples][block]['case_id'],
#                                 samples,
#                                 '.'.join(list(map(lambda x: os.path.basename(x), block.split('.')))),
#                                 ';'.join(list(map(lambda x: os.path.basename(x), blocks[samples][block]['workflows']))),
#                                 blocks[samples][block]['name'],
#                                 blocks[samples][block]['date'],
#                                 blocks[samples][block]['complete'],
#                                 blocks[samples][block]['clean'],
#                                 blocks[samples][block]['network']]
#                             newdata.append(L)
                                    
#         # add data
#         insert_data(database, table, newdata, column_names)
    

# def get_donor_samples_and_library_types(donor_data):
#     '''
#     dict) -> dict

#     Returns a dictionary mapping all samples with their corresponding library types
#     for a given donor

#     Parameters
#     ----------
#     - donor_data (dict): Dictionary with a single donor data   
#     '''
    
#     D = {}
#     for i in range(len(donor_data['cerberus_data'])):
#         library_design = donor_data['cerberus_data'][i]['library_design']
#         sample = '_'.join([donor_data['donor'],
#                           donor_data['cerberus_data'][i]['tissue_type'], 
#                           donor_data['cerberus_data'][i]['tissue_origin'],
#                           library_design, donor_data['cerberus_data'][i]['group_id']])
#         D[sample] = library_design
                           
#     return D
            


# def collect_project_info(provenance_data):
#     '''
#     (list) -> dict
    
#     Returns a dictionary with project level information for each active project in provenance_data
    
#     Parameters
#     ----------
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     '''


#     D = {}


#     for i in range(len(provenance_data)):
#         project_name = get_project_name(provenance_data[i])
#         if project_name not in D:
#             D[project_name] = {'samples': [], 'library_types': []}
#         samples = get_donor_samples_and_library_types(provenance_data[i])
#         pipeline = get_pipeline(provenance_data[i])
#         D[project_name]['pipeline'] = pipeline
#         D[project_name]['samples'].extend(list(samples.keys()))
#         D[project_name]['library_types'].extend(list(samples.values()))
#         D[project_name]['samples'] = list(set(D[project_name]['samples']))
#         D[project_name]['library_types'] = list(set(D[project_name]['library_types']))
    
#     return D    
        
        


# def get_donor_bmpp(donor_data, platform, library_type):
#     '''
#     (dict, str, str) -> list
    
#     Returns a list of bmpp workflow Ids corresponding to the specific case from project 
#     with input sequences from platform and library_type
    
#     Parameters
#     ----------
#     - donor_data (dict): Dictionary with a single donor data   
#     - platform (str): Sequencing platform. Values are novaseq, miseq, nextseq and hiseq
#     - library_type (str): 2-letters code describing the type of the library (eg, WG, WT,..)
#     '''
    
#     L = []
        
#     for i in range(len(donor_data['cerberus_data'])):
#         library_design = donor_data['cerberus_data'][i]['library_design']
#         workflow = donor_data['cerberus_data'][i]['workflow']
#         instrument = donor_data['cerberus_data'][i]['instrument_model']
#         wfrun_id = donor_data['cerberus_data'][i]['workflow_run_accession']
#         if 'bammergepreprocessing' in workflow.lower() and platform in instrument.lower() and library_type == library_design:
#             L.append(wfrun_id)
#     L = list(set(L))        
            
#     return L
    

# def map_samples_to_bmpp_runs(samples_workflows, bmpp_ids):
#     '''
#     (dict, list) -> dict
    
#     Returns a dictionary with normal, tumor samples for each bmpp run id
      
#     Parameters
#     ----------
#     - samples_workflows (dict): Dictionary with workflow information for all samples of a given donor
#     - bmpp_ids (list): List of BamMergePreprocessing workflow run identifiers for a single case
#     '''

#     D = {}
    
    
#     for bmpp_run_id in bmpp_ids:
#         samples = {'normal': [], 'tumour': []}
#         for sample in samples_workflows:
#             for d in samples_workflows[sample]['workflows']:
#                 if d['wfrun_id'] == bmpp_run_id:
#                     if samples_workflows[sample]['tissue_type'] == 'R':
#                         tissue = 'normal'
#                     else:
#                         tissue = 'tumour'
#                     if sample not in samples[tissue]:
#                         samples[tissue].append(sample)
#         D[bmpp_run_id] = samples
        
#     return D
    

# def get_case_call_ready_samples(bmpp_samples):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with all normal and tumor samples processed for all
#     bamMergePreprocessing workflow run ids for a specific donor
    
#     Parameters
#     ----------
#     - bmpp_samples (dict): Dictionary with normal, tumor samples for each bmpp run id
#     '''
        
#     D = {'normal': [], 'tumour': []} 
    
#     for i in bmpp_samples:
#         D['normal'].extend(bmpp_samples[i]['normal'])
#         D['tumour'].extend(bmpp_samples[i]['tumour'])
#     return D


# def group_normal_tumor_pairs(samples):
#     '''
#     (dict) -> list
    
#     Returns a list with all possible normal, tumor pairs from the dictionary
#     of normal and tumor samples for a given donor
    
#     Parameters
#     ----------
#     - samples (dict): Dictionary with normal and tumour samples for a given donor
#     '''
    
#     pairs = []
#     if samples['normal'] and samples['tumour']:
#         for i in samples['normal']:
#             for j in samples['tumour']:
#                 pairs.append(sorted([i, j]))
#     elif samples['normal']:
#         for i in samples['normal']:
#             pairs.append([i])
#     elif samples['tumour']:
#         for i in samples['tumour']:
#             pairs.append([i])
       
#     return pairs


# def exclude_qc_workflows(samples_workflows, qc_workflows):
#     '''
#     (dict, list) -> dict
    
#     Returns a dictionary of workflow information for each sample in samples excluding
#     QC worflows from the qc_workflows list
       
#     Parameters
#     ----------
#     - samples (dict): Dictionary with workflow information for all samples of a given donor
#     - sample (str): Specific sample of donor
#     '''
        
#     D = {}
    
#     for sample in samples_workflows:
#         assert sample not in D
#         D[sample] = {'tissue_type': samples_workflows[sample]['tissue_type'], 'workflows': []}
#         for i in samples_workflows[sample]['workflows']:
#             if i['workflow'] not in qc_workflows:
#                 D[sample]['workflows'].append(i)
        
#     return D                  
            


# def collect_sample_workflows(donor_data, platform):
#     '''
#     (dict, str) -> dict
    
#     Returns a dictionary with all workflows for all samples of a given donor
        
#     Parameters
#     ----------
#     - donor_data (dict): Dictionary with a single donor data   
#     - platform (str): Sequencing platform. Values are: novaseq, nextseq, highseq and miseq
#     '''
    
#     D = {}
        
#     for i in range(len(donor_data['cerberus_data'])):
#         donor = donor_data['donor']
#         tissue_type = donor_data['cerberus_data'][i]['tissue_type']
#         tissue_origin = donor_data['cerberus_data'][i]['tissue_origin']
#         library_type = donor_data['cerberus_data'][i]['library_design']
#         group_id = donor_data['cerberus_data'][i]['group_id']
#         sample = '_'.join([donor, tissue_type, tissue_origin, library_type, group_id])
#         wfrun_id = donor_data['cerberus_data'][i]['workflow_run_accession']
#         workflow = donor_data['cerberus_data'][i]['workflow']
#         instrument = donor_data['cerberus_data'][i]['instrument_model']
        
#         if platform in instrument.lower():
#             if sample not in D:
#                 D[sample] = {'tissue_type': tissue_type, 'workflows': [{'workflow': workflow, 'wfrun_id': wfrun_id}]}
#             else:
#                 if {'workflow': workflow, 'wfrun_id': wfrun_id} not in D[sample]['workflows']:
#                     D[sample]['workflows'].append({'workflow': workflow, 'wfrun_id': wfrun_id})            
        
#     return D


# def find_sample_pairs_with_common_workflows(samples_workflows, pairs):
#     '''
#     (dict, list) -> dict
    
#     Returns a dictionary of sample pairs and common workflows
#     for pairs that have been processed together
    
#     Parameters
#     ----------    
#     - samples_workflows (dict): Dictionary of workflow information for each sample of a given donor
#     - pairs (list): List of normal/tumor samples
#     '''

#     D = {} 

#     for i in pairs:
#         j = ' | '.join(sorted(i))
#         # get the common worflows to each normal/tumor pair
#         L = [k for k in samples_workflows[i[0]]['workflows'] if k in samples_workflows[i[1]]['workflows']]
#         if L:
#             D[j] = L

#     return D



# def identify_WGTS_blocks(D, parent_workflows, bmpp):
#     '''
#     (dict, dict, list) -> list
    
#     Returns a dictionary with analysis workflows grouped by sample pairs and blocks 
#     (ie, originating from common call-ready workflows)
    
#     Parameters
#     ----------
#     - D (dict): Dictionary with workflows of each normal, tumour sample pair
#     - parent_worfklows (dict): Dictionary with parent workflows of each workflow for the normal, tumour sample pairs
#     - bmpp (list): List of bmpp workflow run id for a given donor                           
#     '''
    
#     # track blocks
#     blocks = {}

#     # record all parent workflows
#     L = []

#     # get the bmpp sort bmpp-dowsntream workflows by block and sample     
#     for samples in D:
#         for d in D[samples]:
#             # get the workflow id
#             wfrun_id = d['wfrun_id']
#             # get the parent workflows
#             parent = parent_workflows[wfrun_id]
#             # check that parent workflows are in bmpp list
#             if parent and all([i in bmpp for i in parent]):
#                 # get the anchor workflow
#                 parent_workflow = '.'.join(sorted(parent))
#                 if samples not in blocks:
#                     blocks[samples] = {}
#                 if parent_workflow not in blocks[samples]:
#                     blocks[samples][parent_workflow] = []
#                 blocks[samples][parent_workflow].append({wfrun_id: {'parent': d, 'children': []}})
#                 L.append(wfrun_id)
    
#     # sort workflows downstream of callers by block and sample, 
#     # keeping track of workflow parent-child relationshsips    
#     for samples in blocks:
#         for parent_workflow in blocks[samples]:
#             for j in D[samples]:
#                 if j['wfrun_id'] not in L:
#                     upstream = parent_workflows[j['wfrun_id']]
#                     for k in blocks[samples][parent_workflow]:
#                         for m in upstream:
#                             if m in k:
#                                 k[m]['children'].append(j)
#     return blocks                



# def get_block_analysis_date(donor_data, block_workflows):
#     '''
#     (dict) -> dict, dict
    
#     Returns a dictionary with the most recent analysis date of any workflow downstream
#     of bmpp parent workflows for each sample pair
    
#     Parameters
#     ----------
#     - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
#     - creation_dates (dict): Dictionary with creation date of each worklow in project
#     '''

#     D = {}
    
#     for i in range(len(donor_data['cerberus_data'])):
#         wfrun_id = donor_data['cerberus_data'][i]['workflow_run_accession']
#         creation_date = get_file_timestamp(donor_data['cerberus_data'][i])
#         if wfrun_id not in D:
#             D[wfrun_id] = creation_date
#         else:
#             D[wfrun_id] = max(creation_date, D[wfrun_id])


#     block_date = {}
#     for block in block_workflows:
#         block_date[block] = {}
#         for bmpp in block_workflows[block]:
#             # get the analysis date for each workflow
#             L = [D[wf] for wf in block_workflows[block][bmpp] if wf in D]
#             if L:
#                 block_date[block][bmpp] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(max(L))))
#             else:
#                 block_date[block][bmpp] = 'NA'
#     return block_date



# def list_block_workflows(blocks):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with list of workflow ids grouped by sample pair and analysis block
    
#     Parameters
#     ----------
#     - blocks (dict): Dictionary with analysis workflows grouped by sample pairs and blocks 
#     (ie, originating from common call-ready workflows)
#     '''
        
#     # list all workflows downstream of each bmpp anchors
#     W = {}
#     for block in blocks:
#         W[block] = {}
#         for bmpp in blocks[block]:
#             L = []
#             for d in blocks[block][bmpp]:
#                 for workflow in d:
#                     L.append(d[workflow]['parent']['wfrun_id'])
#                     for k in d[workflow]['children']:
#                         L.append(k['wfrun_id'])
#             L.extend(bmpp.split('.'))
#             L = sorted(list(set(L)))            
#             W[block][bmpp] = L
    
#     return W



# def map_workflow_ids_to_names(samples_workflows):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with workflow_id and workflow name key, value pairs
    
#     Parameters
#     ----------
#     - samples_workflows (dict): Dictionary with workflow information for each sample of a donor
#     '''

#     D = {}

#     for sample in samples_workflows:
#         for d in samples_workflows[sample]['workflows']:
#             wfrun_id = d['wfrun_id']
#             wf = d['workflow']
#             D[wfrun_id] = wf

#     return D


# def get_node_labels(block_workflows, workflow_names):
#     '''
#     (dict, dict) -> dict
    
#     Returns a dictionary with node labels (workflows) for each block (sample pair)
#     and bmpp parent workflows
    
#     Parameters
#     ----------
#     - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
#     - workflow_names (dict): Dictionary with workflow name for each workflow run id in project
#     '''
    
#     labels = {}
#     for block in block_workflows:
#         labels[block] = {}
#         for bmpp in block_workflows[block]:
#             labels[block][bmpp] = []
#             for workflow in block_workflows[block][bmpp]:
#                 workflow_name = workflow_names[workflow]
#                 workflow_name = workflow_name.split('_')[0]
#                 if workflow_name.lower() == 'varianteffectpredictor':
#                     workflow_name = 'VEP'
#                 elif workflow_name.lower() == 'bammergepreprocessing':
#                     workflow_name = 'bmpp'
#                 labels[block][bmpp].append(workflow_name)
       
#     return labels


# def make_adjacency_matrix(block_workflows, parent_workflows):
#     '''
#     (dict, dict) -> dict    
    
#     Returns a dictionary of matrices showing relationships among workflows for each
#     sample and bmpp parent workflows
    
#     Parameters
#     ----------
#     - block_workflows (dict): Dictionary of workflow run ids organized by sample pair and bmpp parent workflows
#     - parent_workflows (dict): Dictionary with parent-children workflow relationships
#     '''
          
#     matrix = {}
    
#     for block in block_workflows:
#         matrix[block] = {}
#         for bmpp in block_workflows[block]:
#             M = []
#             for i in block_workflows[block][bmpp]:
#                 m = []
#                 for j in block_workflows[block][bmpp]:
#                     if i not in parent_workflows and j not in parent_workflows:
#                         m.append(0)
#                     elif i == j:
#                         m.append(0)
#                     elif i in parent_workflows:
#                         if j in parent_workflows[i]:
#                             m.append(1)
#                         else:
#                             m.append(0)
#                     elif i not in parent_workflows:
#                         if i in parent_workflows[j]:
#                             m.append(1)
#                         else:
#                             m.append(0)
#                 M.append(m)
#             matrix[block][bmpp] = M

#     return matrix



# def convert_figure_to_base64(figure):
#     '''
#     (matplotlib.figure.Figure) -> str
    
#     Returns a base64 string representation of figure
    
#     Parameters
#     ----------
#     - figure (matplotlib.figure.Figure): Figure showing workflow relationships
#     '''
    
#     my_stringIObytes = io.BytesIO()
#     figure.savefig(my_stringIObytes, format='png')
    
#     my_stringIObytes.seek(0)
#     my_base64_jpgData = base64.b64encode(my_stringIObytes.read()).decode('utf-8')
      
#     return my_base64_jpgData


# def show_graph(adjacency_matrix, mylabels):
#     '''
#     (dict, dict) -> matplotlib.figure.Figure
    
#     Returns a matplotlib figure of workflow relationships
        
#     Parameters
#     ----------
#     - matrix (dict): Dictionary of matrices of workflow relationships for each
#                      sample pair (block) and bmpp parent workflows
#     - labels (dict): Dictionary with node labels (workflows) for each block (sample pair)
#                      and bmpp parent workflows
#     '''
    
#     figure = plt.figure()
    
#     # add a plot to figure (N row, N column, plot N)
#     ax = figure.add_subplot(1, 1, 1)

#     rows, cols = np.where(adjacency_matrix == 1)
#     edges = zip(rows.tolist(), cols.tolist())
#     gr = nx.Graph()
    
#     gr.add_edges_from(edges)
    
#     nodes = list(gr)
#     N = {}
#     for i in nodes:
#         N[i] = mylabels[i]
    
#     nx.draw(gr, node_size=1200, node_color='#ffe066', font_size = 14, with_labels=True, labels=N, linewidths=2)
    
#     # save figure    
#     plt.tight_layout()
#     #figure.savefig(Outputfile, bbox_inches = 'tight')
#     plt.close()

#     return figure


# def plot_workflow_network(matrix, labels):
#     '''
#     (dict, dict) -> dict
    
#     Returns a dictionary with plots of workflow relationships for each sample pair and 
#     bmpp parent workflows
    
#     Parameters
#     ----------
#     - matrix (dict): Dictionary of matrices of workflow relationships for each
#                      sample pair (block) and bmpp parent workflows
#     - labels (dict): Dictionary with node labels (workflows) for each block (sample pair)
#                      and bmpp parent workflows
#     '''
        
#     F = {}
#     # convert to numpy 2-D array
#     for block in matrix:
#         F[block] = {}
#         for bmpp in matrix[block]:
#             matrix[block][bmpp] = np.array(matrix[block][bmpp])
#             figure = show_graph(matrix[block][bmpp], mylabels=labels[block][bmpp])   
#             # convert to base64
#             F[block][bmpp] = convert_figure_to_base64(figure)
   
#     return F


# def is_block_complete(block_workflows, expected_workflows, workflow_names):
#     '''
#     (dict, list) -> dict
    
#     Returns a dictionary indicating if the downstream workflows of each parent bmpp workflows
#     are complete (ie, contains all the expected workflows) for each block (ie,sample pair)
    
#     Parameters
#     ----------
#     - block_workflows (dict): Dictionary wirh list of workflows for each sample pair and parent bmpp
#     - expected_workflows (list): List of expected generic workflows names
#     '''
    
#     D = {}
    
#     for block in block_workflows:
#         D[block] = {}
#         for bmpp in block_workflows[block]:
#             if len(block_workflows[block][bmpp]) == 0:
#                 complete = False
#             else:
#                 # make a list of caller workflows (ie, remove bmpp)
#                 call_ready = list(map(lambda x: x.strip(), bmpp.split('.')))
#                 callers = set(block_workflows[block][bmpp]).difference(set(call_ready))
#                 workflows = [workflow_names[i] for i in callers]
#                 # homogeneize workflow names by removing the matched suffix
#                 for i in range(len(workflows)):
#                     workflows[i] = workflows[i].split('_')[0]
#                 # check that all workflows are present
#                 complete = set(expected_workflows).issubset(set(workflows))
#             # record boolean as 0 or 1
#             D[block][bmpp] = int(complete)
            
#     return D



# def is_block_clean(block_workflows, expected_workflows):
#     '''
#     (dict, list) -> dict
    
#     Returns a dictionary indicating if each parent bmpp workflow has onlly expected 
#     or extra (more than expected) downstream workflows for each block (ie,sample pair)
    
#     Parameters
#     ----------
#     - block_workflows (dict): Dictionary wirh list of workflows for each sample pair and parent bmpp
#     - expected_workflows (list): List of expected generic workflows names
#     '''
    
#     D = {}
    
#     for block in block_workflows:
#         D[block] = {}
#         for bmpp in block_workflows[block]:
#             L = block_workflows[block][bmpp]
#             if len(L) == 0:
#                 clean = False
#             else:
#                 # remove call ready workflows
#                 call_ready = list(map(lambda x: x.strip(), bmpp.split('.')))
#                 callers = set(L).difference(set(call_ready))
#                 if len(callers) == len(expected_workflows):
#                     clean = True
#                 else:
#                     clean = False
#             # record boolean as 0 or 1
#             D[block][bmpp] = int(clean)
             
#     return D




# def score_blocks(blocks, amount_data, complete, clean):
#     '''
#     (dict, dict, dict, dict) -> dict
    
#     Returns a dictionary with scores (based on amount of data, and block status)
#     for each sub-block (bmpp anchor) for each sample pair (block)
       
#     Parameters
#     ----------
#     - blocks (dict): Dictionary of workflow information organized by sample pair and parent bmpp workflows
#     - amount_data (dict): Dictionary of amount of data for each workflow
#     - complete (dict): Dictionary with block completeness (ie presence of all expected workflows) for each sub-block
#     - clean (dict): Dictionary with block cleaness (ie presence of only expected workflows) for each sub-block 
#     '''
    
#     # rank blocks according to lane count
#     ranks = {}
#     for block in blocks:
#         d = {}
#         for bmpp in blocks[block]:
#             # sum all lanes for all call-ready workflows within each block
#             total = 0
#             workflows = list(map(lambda x: x.strip(), bmpp.split('.')))
#             for workflow_id in workflows:
#                 total += amount_data[os.path.basename(workflow_id)]
#             d[bmpp] = total
#         L = sorted(list(set(d.values())))
#         # get the indices (ranks) for each bmpp anchor
#         # weighted by the number of values
#         ranks[block] = {}
#         for bmpp in d:
#             for i in range(len(L)):
#                 if L[i] == d[bmpp]:
#                     ranks[block][bmpp] = (i + 1) / len(L)
    
#     # score sub-blocks within each block
#     D = {}
#     for block in blocks:
#         for bmpp in blocks[block]:
#             score = 0
#             score += complete[block][bmpp]
#             score += clean[block][bmpp]
#             score += ranks[block][bmpp]
#             if block not in D:
#                 D[block] = {}
#             D[block][bmpp] = score
    
#     return D


# def order_blocks(blocks, amount_data, complete, clean):
#     '''
#     (dict, dict, dict, dict) -> dict
    
#     Returns a dictionary with bmpp parent workflows ordered by amount of lane data
#     and block status for each sample pair
    
#     Parameters
#     ----------
#     - blocks (dict): Dictionary of workflow information organized by sample pair and parent bmpp workflows
#     - amount_data (dict): Dictionary of amount of data for each workflow
#     - complete (dict): Dictionary with block completeness (ie presence of all expected workflows) for each sub-block
#     - clean (dict): Dictionary with block cleaness (ie presence of only expected workflows) for each sub-block 
#     '''
    
#     # score the blocks
#     scores = score_blocks(blocks, amount_data, complete, clean)
    
#     D = {}
#     for block in blocks:
#         L = []
#         for bmpp in blocks[block]:
#             L.append([scores[block][bmpp], bmpp])
#         # sort sub-blocks according to score    
#         L.sort(key=lambda x: x[0])
#         D[block] = [i[1] for i in L]
#     return D



# def name_blocks(ordered_blocks, pipeline):
#     '''
#     (dict) -> dict  
    
#     Returns a dictionary with sub-blocks names for each sample pair (ie, block)
    
#     Parameters
#     ----------
#     - ordered_blocks (dict): Dictionary with bmpp parent worflows ordered by amount of data for each sample pair
#     - pipeline (str): Name of pipeline. Values are WG or WT
#     '''
        
#     names = {}
#     for block in ordered_blocks:
#         counter = 1
#         names[block] = {}
#         # loop in reverse order, list is sorted according to ascending scores
#         for i in ordered_blocks[block][::-1]:
#             k = '{0} Analysis Block {1}'.format(pipeline, counter)
#             names[block][i] = k
#             counter += 1
#     return names



# def find_donor_WGS_blocks(donor_data, library_type, platform, expected_workflows, qc_workflows):
#     '''
#     (str, str, str, list) -> dict
    
#     Returns a dictionary with the WGS blocks for case in project
    
#     Parameters
#     ----------
#     - donor_data (dict): Dictionary with a single donor data   
#     - library_type (str): 2-letters code describing the type of the library (eg, WG, WT,..)
#     - platform (str): Sequencing platform. Values are novaseq, miseq, nextseq and hiseq
#     - expected_workflows (list): List of expected workflow names to define a complete block
#     - qc_workflows (list): List of QC workflows to exclude
#     '''
    
#     # organize the WGS blocks as a dictionary
#     WGS_blocks = {}
    
#     # get all the bmpp runs for WG library type and Novaseq platform
#     bmpp = get_donor_bmpp(donor_data, platform, library_type)

#     # proceed only in bmpp ids exist
#     if bmpp:
#         # get workflows of all samples for the donor
#         samples_workflows = collect_sample_workflows(donor_data, platform)
#         bmpp_samples = map_samples_to_bmpp_runs(samples_workflows, bmpp)
#         # identify all the samples processed
#         samples = get_case_call_ready_samples(bmpp_samples)
#         # proceed only if tumor/normal samples exist
#         if samples['normal'] and samples['tumour']:
#             # get all pairs N/T samples
#             pairs = group_normal_tumor_pairs(samples)
#             # exclude QC workflows
#             samples_workflows = exclude_qc_workflows(samples_workflows, qc_workflows)
#             # find common workflows to each normal/tumor pair    
#             D = find_sample_pairs_with_common_workflows(samples_workflows, pairs)
            
#             if D:
#                 # find the parents of each workflow
#                 files = map_file_to_worklow(donor_data)
#                 workflow_inputs = get_workflow_inputs(donor_data)
#                 parent_workflows = identify_parent_children_workflows(workflow_inputs, files)
#                 # find the blocks for that donor           
#                 blocks = identify_WGTS_blocks(D, parent_workflows, bmpp)
#                 # list all workflows for each block
#                 block_workflows = list_block_workflows(blocks)
#                 # get the date of each block
#                 block_date = get_block_analysis_date(donor_data, block_workflows) 
#                 # map each workflow run id to its workflow name
#                 workflow_names = map_workflow_ids_to_names(samples_workflows)
#                 # get the workflow names
#                 block_workflow_names = get_node_labels(block_workflows, workflow_names)
#                 # convert workflow relationships to adjacency matrix for each block
#                 matrix = make_adjacency_matrix(block_workflows, parent_workflows)
#                 # create figures
#                 figures = plot_workflow_network(matrix, block_workflow_names)
#                 # check if blocks are complete
#                 complete = is_block_complete(block_workflows, expected_workflows, workflow_names)
#                 # check if blocks have extra workflows
#                 clean = is_block_clean(block_workflows, expected_workflows)
#                 # get the amount of data for each workflow
#                 workflow_info = collect_donor_workflow_info(donor_data)
#                 amount_data = {i: workflow_info[i]['lane_count'] for i in workflow_info}
#                 # order blocks based on scores
#                 ordered_blocks = order_blocks(blocks, amount_data, complete, clean)
#                 # name each block according to the selected block order
#                 names = name_blocks(ordered_blocks, 'WG')
               
#                 for samples in blocks:
#                     WGS_blocks[samples] = {}
#                     for block in blocks[samples]:
#                         WGS_blocks[samples][block] = {}
#                         # record network image
#                         WGS_blocks[samples][block]['network'] = figures[samples][block]
#                         # record all workflow ids
#                         WGS_blocks[samples][block]['workflows'] = block_workflows[samples][block]
#                         # record block date
#                         WGS_blocks[samples][block]['date'] = block_date[samples][block]
#                         # record complete status
#                         WGS_blocks[samples][block]['complete'] = complete[samples][block]
#                         # record extra workflow status
#                         WGS_blocks[samples][block]['clean'] = clean[samples][block]
#                         # record block name
#                         WGS_blocks[samples][block]['name'] = names[samples][block]
#                         # add project and case ids
#                         WGS_blocks[samples][block]['project_id'] = get_project_name(donor_data)
#                         WGS_blocks[samples][block]['case_id'] = donor_data['donor']
    
#     return WGS_blocks




# def get_donor_star(donor_data, platform, library_type):
#     '''
#     (dict, str, str) -> list
    
#     Returns a list of star workflow Ids corresponding to the specific donor  
#     with input sequences from platform and library_type
    
#     Parameters
#     ----------
#     - donor_data (dict): Dictionary with a single donor data   
#     - platform (str): Sequencing platform. Values are novaseq, miseq, nextseq and hiseq
#     - library_type (str): 2-letters code describing the type of the library (eg, WG, WT,..)
#     '''
    
#     L = []
        
#     for i in range(len(donor_data['cerberus_data'])):
#         library_design = donor_data['cerberus_data'][i]['library_design']
#         workflow = donor_data['cerberus_data'][i]['workflow']
#         instrument = donor_data['cerberus_data'][i]['instrument_model']
#         wfrun_id = donor_data['cerberus_data'][i]['workflow_run_accession']
#         if 'star_call_ready' in workflow.lower() and platform in instrument.lower() and library_type == library_design:
#             L.append(wfrun_id)
#     L = list(set(L))        
            
#     return L



# def map_samples_to_star_runs(samples_workflows, star_ids):
#     '''
#     (dict, list) -> dict
        
#     Returns a dictionary with tumor samples for each star run id
          
#     Parameters
#     ----------
#     - samples_workflows (dict): Dictionary with workflow information for all samples of a given donor
#     - star_ids (list): List of BamMergePreprocessing workflow run identifiers for a single case
#     '''

#     D = {}
    
#     for star_run_id in star_ids:
#         samples = {'tumour': []}
#         for sample in samples_workflows:
#             for d in samples_workflows[sample]['workflows']:
#                 if d['wfrun_id'] == star_run_id:
#                     if samples_workflows[sample]['tissue_type'] == 'R':
#                         tissue = 'normal'
#                     else:
#                         tissue = 'tumour'
#                     if tissue == 'tumour':
#                         if sample not in samples[tissue]:
#                             samples[tissue].append(sample)
#         D[star_run_id] = samples                
    
#     return D

    
# def get_WT_case_call_ready_samples(star_samples):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with tumor samples processed for all star_call_ready
#     workflow run ids for a specific donor
        
#     Parameters
#     ----------
#     - star_samples (dict): Dictionary with tumor samples for each star run id
#     '''
     
#     D = {'tumour': []}
#     for i in star_samples:
#         D['tumour'].extend(star_samples[i]['tumour'])
#     return D    
    

# def map_workflows_to_samples(samples_workflows, samples):
#     '''
#     (dict, dict) -> dict
    
#     Returns a dictionary with workflows processed by the WT pipeline for tumor samples
    
#     Parameters
#     ----------
#     - samples_workflows (dict): Dictionary with workflows for all samples of a given donor
#     - samples (dict): Dictionary with tumor samples
#     '''
    
#     D = {}
#     for i in samples['tumour']:
#         D[i] = samples_workflows[i]['workflows']
        
#     return D    
        
    
# def find_donor_WT_blocks(donor_data, library_type, platform, expected_workflows, qc_workflows):
#     '''
#     (str, str, str, list) -> dict
    
#     Returns a dictionary with the WGS blocks for case in project
    
#     Parameters
#     ----------
#     - donor_data (dict): Dictionary with a single donor data   
#     - library_type (str): 2-letters code describing the type of the library (eg, WG, WT,..)
#     - platform (str): Sequencing platform. Values are novaseq, miseq, nextseq and hiseq
#     - expected_workflows (list): List of expected workflow names to define a complete block
#     - qc_workflows (list): List of QC workflows to exclude
#     '''
    
    
#     WT_blocks = {}
    
#     # build the somatic calling block
#     # identify all call ready star runs for novaseq
#     star = get_donor_star(donor_data, platform, library_type)

#     if star:
#         # get workflows of all samples for the donor
#         samples_workflows = collect_sample_workflows(donor_data, platform)
#         star_samples = map_samples_to_star_runs(samples_workflows, star)
#         # identify all the samples processed
#         samples = get_WT_case_call_ready_samples(star_samples)
#         if samples['tumour']:
#             # exclude QC workflows
#             samples_workflows = exclude_qc_workflows(samples_workflows, qc_workflows)
#             D = map_workflows_to_samples(samples_workflows, samples)
#             if D:
#                 # find the parents of each workflow
#                 files = map_file_to_worklow(donor_data)
#                 workflow_inputs = get_workflow_inputs(donor_data)
#                 parent_workflows = identify_parent_children_workflows(workflow_inputs, files)
#                 # find the blocks for that donor           
#                 blocks = identify_WGTS_blocks(D, parent_workflows, star)
#                 # list all workflows for each block
#                 block_workflows = list_block_workflows(blocks)
#                 # get the date of each block
#                 block_date = get_block_analysis_date(donor_data, block_workflows) 
#                 # map each workflow run id to its workflow name
#                 workflow_names = map_workflow_ids_to_names(samples_workflows)
#                 # get the workflow names
#                 block_workflow_names = get_node_labels(block_workflows, workflow_names)
#                 # convert workflow relationships to adjacency matrix for each block
#                 matrix = make_adjacency_matrix(block_workflows, parent_workflows)
#                 # create figures
#                 figures = plot_workflow_network(matrix, block_workflow_names)
#                 # check if blocks are complete
#                 complete = is_block_complete(block_workflows, expected_workflows, workflow_names)
#                 # check if blocks have extra workflows
#                 clean = is_block_clean(block_workflows, expected_workflows)
#                 # get the amount of data for each workflow
#                 workflow_info = collect_donor_workflow_info(donor_data)
#                 amount_data = {i: workflow_info[i]['lane_count'] for i in workflow_info}
#                 # order blocks based on scores
#                 ordered_blocks = order_blocks(blocks, amount_data, complete, clean)
#                 # name each block according to the selected block order
#                 names = name_blocks(ordered_blocks, 'WT')
            
                           
#                 for samples in blocks:
#                     WT_blocks[samples] = {}
#                     for block in blocks[samples]:
#                         WT_blocks[samples][block] = {}
#                         # record network image
#                         WT_blocks[samples][block]['network'] = figures[samples][block]
#                         # record all workflow ids
#                         WT_blocks[samples][block]['workflows'] = block_workflows[samples][block]
#                         # record block date
#                         WT_blocks[samples][block]['date'] = block_date[samples][block]
#                         # record complete status
#                         WT_blocks[samples][block]['complete'] = complete[samples][block]
#                         # record extra workflow status
#                         WT_blocks[samples][block]['clean'] = clean[samples][block]
#                         # record block name
#                         WT_blocks[samples][block]['name'] = names[samples][block]
#                         # add project and case ids
#                         WT_blocks[samples][block]['project_id'] = get_project_name(donor_data)
#                         WT_blocks[samples][block]['case_id'] = donor_data['donor']
    
#     return WT_blocks
            


def generate_database(database, provenance_data_file):
    '''
    (str, str) -> None

    Generates the waterzooi database using data in the provenance_data_file and 
    calcontaqc_db database

    Parameters
    ----------
    - database (str): Path of the waterzooi database
    - provenance_data_file (str): Path to the provenance data json
    '''
    
    # create database if file doesn't exist
    if os.path.isfile(database) == False:
        initiate_db(database)
    print('initiated database')
    
    # load data from file
    provenance_data = load_data(provenance_data_file)
    print('loaded data')
    
    # collect the md5sum of each donor's data
    md5sums = compute_case_md5sum(provenance_data)
    print('computed md5sums')
    # collect the recorded md5sums of the donor data from the database
    recorded_md5sums = get_cases_md5sum(database, table = 'Checksums')
    print('pulled md5sums from database')
    # determine the donors that need an update
    cases_to_update = cases_info_to_update(md5sums, recorded_md5sums)
    print('determined donors to update')
    
    # add project information
    add_project_info_to_db(database, provenance_data, 'Projects')
    print('added project info to database')
    # add library information
    add_library_info_to_db(database, provenance_data, cases_to_update, 'Libraries')
    print('added library info to database')
    # add sample information
    add_samples_info_to_db(database, provenance_data, cases_to_update, 'Samples')
    print('added sample information to database')
    # add workflow inputs
    add_workflow_inputs_to_db(database, provenance_data, cases_to_update, 'Workflow_Inputs')
    print('added workflow inputs to database')
    # add workflow relationships
    add_workflows_relationships_to_db(database, provenance_data, cases_to_update, 'Parents')
    print('added workflow relationships to database')
    # add workflow information
    add_workflows_to_db(database, provenance_data, cases_to_update, 'Workflows')
    print('added workflow info to database')
    # add file information
    add_file_info_to_db(database, provenance_data, cases_to_update, 'Files')
    print('added file info to database')
    
    
    
    # update the checksums for donors
    # add_checksums_info_to_db(database, cases_to_update, 'Checksums')
    # print('added md5sums to database')
    
 
   
 
generate_database('waterzooi_db_case.db', 'CaseInfo_IRIS_NEOPOC.dump.json')    
     
 
    
 
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(prog = 'waterzooi_refiller.py', description='Script to add data to the waterzooi database')
#     parser.add_argument('-p', '--provenance', dest = 'provenance', help = 'Path to the provenance reporter data json', required=True)
#     parser.add_argument('-cq', '--calcontaqc', dest = 'calcontaqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest', help = 'Path to the merged rnaseq calculateContamination database. Default is /scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest')
#     parser.add_argument('-db', '--waterzooi_db', dest='database', help='Path to the waterzooi database', required=True)    
#     parser.add_argument('-r', '--review', dest='analysis_db', help='Path to analysis_review database', required=True)    
    
    
#     # get arguments from the command line
#     args = parser.parse_args()
#     #args.func(args)
    
#     generate_database(args.database, args.provenance, args.calcontaqc_db)    

