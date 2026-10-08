# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 15:30:32 2025

@author: rjovelin
"""

import hashlib
import json
import sqlite3


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


def is_case_info_complete(case_data):
    '''
    (dict) -> bool
    
    Returns True if the case information is complete
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''
    
    complete = True
    for i in case_data:
        if len(case_data[i]) == 0:
            complete = False
            break
    
    return complete


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
    
    if len(donor) != 1:
        print('donor', donor)
    
    assert len(donor) == 1
    donor = donor[0]

    return donor



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



def define_columns(database):
    '''
    (str) -> dict

    Returns a dictionary with column names and types for each table in database
 
    Parameters
    ----------
    - database (str): Name of the cache. Accepted values: waterzooi and analysis_review
    '''

    # create dict to store column names for each table {table: [column names]}
    if database == 'waterzooi':
        columns = {'Workflows': {'names': ['wfrun_id', 'wf', 'wfv', 'case_id', 'project_id', 'donor_id', 'file_count', 'lane_count'],
                                      'types': ['VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)', 'INT', 'INT']},
                        'Parents': {'names': ['parents_id', 'children_id', 'project_id', 'case_id', 'donor_id'],
                                    'types': ['VARCHAR(572)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)']},
                        'Projects': {'names': ['project_id', 'pipeline', 'last_updated', 'cases', 'samples',
                                               'library_types', 'assays', 'deliverables', 'active'],
                                     'types': ['VARCHAR(128) PRIMARY KEY NOT NULL UNIQUE', 'VARCHAR(128)',
                                               'VARCHAR(256)', 'INT', 'INT', 'VARCHAR(256)', 'VARCHAR(572)',
                                               'VARCHAR(572)', 'INT']},
                        'Files': {'names': ['file_swid', 'project_id', 'md5sum', 'wfrun_id',
                                            'file', 'attributes', 'creation_date', 'limskey',
                                            'case_id', 'donor_id'],
                                  'types': ['VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(572)',
                                            'TEXT', 'TEXT', 'INT', 'VARCHAR(256)',
                                            'VARCHAR(256)', 'VARCHAR(128)']},
                        'Libraries': {'names': ['library', 'lims_id', 'sample_id', 'case_id',
                                                'donor_id', 'tissue_type', 'tissue_origin',
                                                'library_type', 'group_id', 'group_id_description', 'project_id'],
                                      'types': ['VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(572)',
                                                'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                                'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(128)']},
                        'Workflow_Inputs': {'names': ['library', 'run', 'lane', 'wfrun_id', 'limskey',
                                                      'barcode', 'platform', 'project_id', 'case_id', 'donor_id'],
                                            'types': ['VARCHAR(128)', 'VARCHAR(256)', 'INTEGER', 'VARCHAR(572)',
                                                      'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                                      'VARCHAR(572)', 'VARCHAR(128)']},
                        'Samples': {'names': ['case_id', 'assay', 'donor_id', 'ext_id', 'species',
                                              'miso', 'project_id', 'sequencing_status'],
                                    'types': ['VARCHAR(572)', 'VARCHAR(256)',  'VARCHAR(128)', 'VARCHAR(256)',
                                              'VARCHAR(256)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)']},
                        'Checksums': {'names': ['project_id', 'case_id', 'donor_id', 'md5'],
                                      'types': ['VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(572)']}}
                                                  
                        
    elif database == 'analysis_review':
        columns = {'templates': {'names':  ['case_id', 'donor_id', 'project_id', 'assay',
                                            'template', 'valid', 'error', 'md5'],
                                'types': ['VARCHAR(572)', 'VARCHAR(572)', 'VARCHAR(572)',
                                          'VARCHAR(572)', 'TEXT', 'TEXT', 'VARCHAR(572)',
                                          'VARCHAR(572)']}}
    else:
        columns = {}
        
    return columns


def create_table(database_name, database, table):
    '''
    (str, str, str, dict) -> None
    
    Creates a table in database
    
    Parameters
    ----------
    - database_name (str): Name of the database
    - table (str): Table name
    - database (str):  Name of the cache. Accepted values: waterzooi and analysis_review
    '''
    
    # get the column names and types
    columns = define_columns(database)
    # get the column names and types
    column_names, column_types = columns[table]['names'], columns[table]['types']    
    
    # define table format including constraints    
    table_format = ', '.join(list(map(lambda x: ' '.join(x), list(zip(column_names, column_types)))))

    if database == 'waterzooi':
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
    conn = sqlite3.connect(database_name)
    cur = conn.cursor()
    # create table
    cmd = 'CREATE TABLE {0} ({1})'.format(table, table_format)
    cur.execute(cmd)
    conn.commit()
    conn.close()

def initiate_db(database_name, database, tables):
    '''
    (str, str, list) -> None
    
    Create tables in database
    
    Parameters
    ----------
    - database (str): Path to the database file
    - tables (list): List of tables in database
    '''
    
    # check if table exists
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    current_tables = cur.fetchall()
    current_tables = [i[0] for i in current_tables]    
    conn.close()
    
    for i in tables:
        if i not in current_tables:
            create_table(database_name, database, i)


def insert_multiple_records(data, conn, database, table, column_names):
    '''
    (list, sqlite3.Connection, str, str, list) -> None
    
    Inserts data into the database table with column names 
    
    Parameters
    ----------
    - data (list): List of data to be inserted
    - database (str): Path to the database file
    - conn (sqlite3.Connection): Open connection to the database
    - table (str): Table in database
    - column_names (list): List of table column names
    '''
       
    vals = '(' + ','.join(['?'] * len(data[0])) + ')'
    conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), data)
    conn.commit()


def delete_unique_record(identifier, conn, database, table, field):
    '''
    (str, sqlite3.Connection, str, str, str) -> None
     
    Remove all the rows from table with identifier in field
        
    Parameters
    ----------
    - conn (sqlite3.Connection): Open connection to the database
    - identifer (str): Item in table 
    - database (str): Path to the sqlite database
    - table (str): Table in the database
    - field (str): Column in table
    '''
        
    query = "DELETE FROM {0} WHERE {1} = \"{2}\"".format(table, field, identifier)
    conn.execute(query)
    conn.commit()


def delete_multiple_records(L, conn, database, table, field):
    '''
    (list, str, str) -> None
     
    Remove all the rows from table with items in L
        
    Parameters
    ----------
    - L (list): List of items to remove
    - conn (sqlite3.Connection): Open connection to the database
    - database (str): Path to the sqlite database
    - table (str): Table in database
    - field (str): Column in table
    '''
        
    query = "DELETE FROM {0} WHERE {1} IN ({2})".format(table, field, ", ".join("?" * len(L)))
    conn.execute(query, L)
    conn.commit()
