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
                    'Libraries': ['library', 'lims_id', 'sample_id', 'case_id', 'donor_id', 'tissue_type', 'tissue_origin',
                                  'library_type', 'group_id', 'group_id_description', 'project_id'],
                    'Workflow_Inputs': ['library', 'run', 'lane', 'wfrun_id', 'limskey', 'barcode', 'platform', 'project_id', 'case_id', 'donor_id'],
                    'Samples': ['case_id', 'assay', 'donor_id', 'ext_id', 'species', 'sex', 'miso', 'project_id', 'sequencing_status'],
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
                    'Libraries': ['VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(572)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)',
                                  'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(128)'],
                    'Workflow_Inputs': ['VARCHAR(128)', 'VARCHAR(256)', 'INTEGER', 'VARCHAR(572)', 
                                        'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)'],
                    'Samples': ['VARCHAR(572)', 'VARCHAR(256)',  'VARCHAR(128)', 'VARCHAR(256)', 'VARCHAR(256)', 'VARCHAR(128)', 'VARCHAR(572)', 'VARCHAR(128)', 'VARCHAR(128)'],
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
              'Workflow_Inputs', 'Samples', 'Checksums']:
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
        files = json.loads(d['files'])
        for i in range(len(files)):
            file_swid = files[i]['accession']
            file_path = files[i]['path']
            md5sum = files[i]['md5']
            attributes = files[i]['file_attributes']
            if attributes:
                attributes = json.loads(attributes)
                file_attributes = {}
                for k in attributes:
                    file_attributes[k] = attributes[k][0]
                if len(file_attributes) == 0:
                    file_attributes = ''
                else:
                    file_attributes = json.dumps(file_attributes)
            else:
                file_attributes = ''
            creation_date = files[i]['timestamp'].replace('Z', '')
            creation_date = creation_date.split('.')[0]
            creation_date = ' '.join(creation_date.split('T'))
            pattern = '%Y-%m-%d %H:%M:%S'
            creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
            
            
            assert file_swid not in D
            D[file_swid] = {'file_swid': file_swid,
                            'md5sum': md5sum,
                            'creation_date': creation_date,
                            'attributes': file_attributes,
                            'file': file_path,
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
             'tissue_origin': tissue_origin, 'lims_id': lims_id, 'sample_id': sample}
        
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
        file_count = len(json.loads(d['files']))   
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
                        'sample': d['sampleId'],
                        'donor': d['donor'],
                        'group_id': d['groupId'],
                        'library_type': d['libraryDesign'],                   
                        'tissue_origin': d['tissueOrigin'],
                        'tissue_type': d['tissueType'],
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
        platform = d['platform']
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
                 'limskey': limskey,
                 'platform': platform}
            
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
        assay = case_data['assay']
        incomplete = case_data['incomplete']
        sequencing_status = not incomplete
               
        ext_id = d['ext_id']
        sex = d['sex']
        miso = 'NA'
        species = d['organism']
        
        d = {'case_id': case, 'assay': assay, 'donor_id': donor, 'ext_id': ext_id, 'sex': sex,
             'species': species, 'project_id': project_id, 'miso': miso, 'sequencing_status': str(int(sequencing_status))}
        
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
    add_checksums_info_to_db(database, cases_to_update, 'Checksums')
    print('added md5sums to database')
    
 
   
 
#generate_database('waterzooi_db_case.db', 'CaseInfo_IRIS_NEOPOC.dump.json')    
     
generate_database('waterzooi_db_case.db', 'case_shesmu_data.dump.json')    
 
    
 
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

