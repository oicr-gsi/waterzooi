# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:37:47 2025

@author: rjovelin
"""


import sqlite3


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
                                      'types': ['VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(572)']},
                        'File_qc': {'names': ['project_id', 'case_id', 'wfrun_id', 'file_swid',
                                              'filepath', 'username', 'ticket', 'qcstatus'],
                                    'types': ['VARCHAR(256)', 'VARCHAR(572)', 'VARCHAR(572)',
                                              'VARCHAR(572)', 'VARCHAR(572)', 'VARCHAR(128)',
                                              'VARCHAR(128)', 'VARCHAR(128)']}}
    elif database == 'analysis_review':
        columns = {'templates': {'names':  ['case_id', 'donor', 'project', 'assay',
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
        if table  in ['Workflows', 'Parents', 'Files', 'Libraries', 'Workflow_Inputs', 'Samples', 'Checksums', 'File_qc']:
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

        if table in ['Libraries', 'File_qc']:
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

 
    
    

