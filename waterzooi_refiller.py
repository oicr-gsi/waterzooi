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
import requests

from db_helper import connect_to_db, define_columns, initiate_db, insert_multiple_records, \
    delete_unique_record, delete_multiple_records
from data_helper import load_data, clean_up_data, clean_up_workflows, is_case_info_complete
from commons import compute_md5, get_cases_md5sum, find_sequencing_attributes, \
    convert_to_bool, list_case_workflows, get_donor_name, \
    case_to_update    
     
        

    


# def add_checksums_info_to_db(database, cases_to_update, table = 'Checksums'):
#     '''
#     (str, dict, str) -> None
    
#     Update table Checksums with the checksum of the case info 
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Name of Table in database. Default is Checksums
#     '''
    
#     if cases_to_update:
#         delete_records(cases_to_update, database, table)
        
#         # make a list of data to insert
#         newdata = []
        
#         # connect to db
#         conn = sqlite3.connect(database, timeout=30)
#         # get column names
#         data = conn.execute("SELECT * FROM {0}".format(table))
#         column_names = [column[0] for column in data.description]

#         # order values according to column names
#         for i in cases_to_update:
#             L = [cases_to_update[i]['project_id'], i, cases_to_update[i]['donor_id'],  cases_to_update[i]['md5']]
#             newdata.append(L)
#         vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#         conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#         conn.commit()
#         conn.close()





# def get_file_timestamp(d):
#     '''
#     (dict) - > int
    
#     Returns the time stamp of a file from the dictionary of the file from cerberus data
#     for a donor
    
#     Parameters
#     ----------
#     - d (dict): Dictionary representing a file information from cerberus
#     '''    
    
#     creation_date = d['timestamp']
#     creation_date = ' '.join(creation_date.split('T')).replace('Z', '')
#     creation_date = creation_date.split('.')
#     if len(creation_date) > 1:
#         creation_date = ' '.join(creation_date[:-1])
#     else:
#         creation_date = creation_date[0]
#     pattern = '%Y-%m-%d %H:%M:%S'
#     creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
    
#     return creation_date



# def collect_case_file_info(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with all the file information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''    
        
#     D = {}
    
#     for d in case_data['workflow_runs']:
#         donor = get_donor_name(case_data)
#         case = case_data['case']
#         project_id = case_data['project_info'][0]['project']
#         wfrun_id = d['wfrunid']
#         limskeys = d['limsIds'].split(',')
#         files = json.loads(d['files'])
#         for i in range(len(files)):
#             file_swid = files[i]['accession']
#             file_path = files[i]['path']
#             md5sum = files[i]['md5']
#             attributes = files[i]['file_attributes']
#             if attributes:
#                 attributes = json.loads(attributes)
#                 file_attributes = {}
#                 for k in attributes:
#                     file_attributes[k] = attributes[k][0]
#                 if len(file_attributes) == 0:
#                     file_attributes = ''
#                 else:
#                     file_attributes = json.dumps(file_attributes)
#             else:
#                 file_attributes = ''
#             creation_date = files[i]['timestamp'].replace('Z', '')
#             creation_date = creation_date.split('.')[0]
#             creation_date = ' '.join(creation_date.split('T'))
#             pattern = '%Y-%m-%d %H:%M:%S'
#             creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
            
            
#             assert file_swid not in D
#             D[file_swid] = {'file_swid': file_swid,
#                             'md5sum': md5sum,
#                             'creation_date': creation_date,
#                             'attributes': file_attributes,
#                             'file': file_path,
#                             'project_id': project_id,
#                             'case_id': case,
#                             'donor_id': donor,
#                             'wfrun_id': wfrun_id,
#                             'limskey': limskeys}  
                                
#     return D


# def add_file_info_to_db(database, provenance_data, cases_to_update, table = 'Files'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts file information in database's Files table
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Files
#     '''
        
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))
    
#         # make a list of data to insert in the database
#         newdata = []
            
#         for case_data in provenance_data:
#             donor = get_donor_name(case_data)
#             case =  case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 file_info = collect_case_file_info(case_data)
#                 for file_swid in file_info:
#                     for limskey in file_info[file_swid]['limskey']:
#                         L = [file_info[file_swid][i] for i in column_names[:column_names.index('limskey')]]
#                         L.append(limskey)
#                         L.extend([file_info[file_swid][i] for i in column_names[column_names.index('limskey')+1:]])
#                         newdata.append(L)             
    
#         # add data
#         insert_data(database, table, newdata, column_names)


# def collect_case_library_info(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with all the library information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''    

#     D = {}
    
#     for d in case_data['sample_info']:
#         donor = d['donor']
#         library = d['library']
#         case = case_data['case']
#         tissue_origin = d['tissueOrigin']
#         tissue_type = d['tissueType']
#         library_type = d['libraryDesign']
#         project = d['project']
#         group_id = d['groupId']
#         group_id_description = d['groupDesc']
#         lims_id = d['limsId']
#         sample = d['sampleId']
               
#         d = {'library': library, 'case_id': case, 'donor_id': donor, 
#              'project_id': project, 'group_id': group_id, 
#              'group_id_description': group_id_description,
#              'library_type': library_type, 'tissue_type': tissue_type,
#              'tissue_origin': tissue_origin, 'lims_id': lims_id, 'sample_id': sample}
        
#         if library not in D:
#             D[library] = [d] 
#         else:
#             if d not in D[library]:
#                 D[library].append(d)
                   
#     return D       


# def add_library_info_to_db(database, provenance_data, cases_to_update, table = 'Libraries'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts library information in the Libraries table of the database
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Libraries
#     '''                 
 
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for cases to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
    
#         for case_data in provenance_data:
#             case =  case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 library_info = collect_case_library_info(case_data)
#                 for library in library_info:
#                     for d in library_info[library]:
#                         L = [d[i] for i in column_names]
#                         newdata.append(L)             
  
#         # add data
#         insert_data(database, table, newdata, column_names)




# def get_pipeline(case_data):
#     '''
#     (dict) -> str
    
#     Returns the pipeline the project extracted from the pinery project data of a case 
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data    
#     '''

#     L = list(set([i['pipeline'] for i in case_data['pinery_project_data']]))
#     assert len(L) == 1
#     pipeline = L[0]    
#     return pipeline


# def collect_project_info(provenance_data):
#     '''
#     (list) -> dict
    
#     Returns a dictionary with project level information for each project in provenance_data
    
#     Parameters
#     ----------
#     - provenance_data (list): List of dictionaries with case information 
#     '''

#     D = {}
    
#     for case_data in provenance_data:
#         #project = case_data['project']
#         project = case_data['project_info'][0]['project']
#         case = case_data['case']
#         assay = case_data['assay']
#         #assert len(case_data['project_info']) == 1
#         active = case_data['project_info'][0]['isActive']
#         pipeline = case_data['project_info'][0]['pipeline']
#         deliverables = case_data['project_info'][0]['deliverables']
#         samples, library_types = [], []
#         for d in case_data['sample_info']:
#             samples.append(d['sampleId'])
#             library_types.append(d['libraryDesign'])
#         if project not in D:
#             D[project] = {'assays': [assay],
#                           'cases': [case],
#                           'active': active,
#                           'pipeline': pipeline,
#                           'deliverables': deliverables,
#                           'samples': samples,
#                           'library_types': library_types}
#         else:
#             D[project]['assays'].append(assay)
#             D[project]['cases'].append(case)
#             D[project]['samples'].extend(samples)
#             D[project]['library_types'].extend(library_types)
                                
#         D[project]['assays'] = list(set(D[project]['assays']))
#         D[project]['cases'] = list(set(D[project]['cases']))
#         D[project]['samples'] = list(set(D[project]['samples']))
#         D[project]['library_types'] = list(set(D[project]['library_types']))
            
#     return D


# def add_project_info_to_db(database, provenance_data, project_table = 'Projects'):
#     '''
#     (str, list, str) -> None
    
#     Add project information into Projects table of database
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     - project_table (str): Name of Table storing project information in database. Default is Projects
#     '''
    
#     # get the column names
#     column_names = define_column_names()[project_table]
    
#     # collect project information
#     project_info = collect_project_info(provenance_data)
    
#     # make a list of data to insert in the database
#     newdata = []
    
#     # connect to db
#     conn = connect_to_db(database)
       
#     for project in project_info:
#         last_updated = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
        
#         # get the library types
#         library_types = ','.join(sorted(project_info[project]['library_types']))
#         samples = len(project_info[project]['samples'])
#         cases = len(project_info[project]['cases'])
#         assays = ','.join(sorted(project_info[project]['assays']))
        
#         L = [project, project_info[project]['pipeline'], last_updated, str(cases),
#              str(samples), str(library_types), assays, project_info[project]['deliverables'],
#              project_info[project]['active']]
        
#         newdata.append(L)
#         # delete project entry
#         query = 'DELETE FROM {0} WHERE project_id = \"{1}\"'.format(project_table, project)
#         conn.execute(query)
#         conn.commit()

#     conn.close()
        
#     # add data
#     insert_data(database, project_table, newdata, column_names)



# def collect_workflow_info(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with all the workflow information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''

#     D = {}
    
#     for d in case_data['workflow_runs']:
#         project_id = case_data['project_info'][0]['project']
#         case = case_data['case']
#         donor = get_donor_name(case_data)
#         wfrun_id = d['wfrunid']
#         wf = d['wf']
#         wfv = d['wfv']
#         limskeys = d['limsIds'].split(',')
#         file_count = len(json.loads(d['files']))   
#         sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
#         lane_count = len([sequencing_attributes[i]['lane'] for i in sequencing_attributes])
         
#         data = {'project_id': project_id,
#                 'wfrun_id': wfrun_id,
#                 'wf': wf,
#                 'wfv': wfv,
#                 'case_id': case,
#                 'donor_id': donor,
#                 'file_count': file_count,
#                 'lane_count': lane_count}
                
#         assert wfrun_id not in D
#         D[wfrun_id] = data
                
#     return D            



# def add_workflows_to_db(database, provenance_data, cases_to_update, table = 'Workflows'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts or updates workflow information 
 
#     Parameters
#     ----------    
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Workflows
#     '''
    
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
        
#         for case_data in provenance_data:
#             case = case_data['case']
#             donor = get_donor_name(case_data)
                      
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 workflow_info = collect_workflow_info(case_data)                
#                 for workflow in workflow_info:
#                     L = [workflow_info[workflow][i] for i in column_names]
#                     newdata.append(L)             
      
#         # add data
#         insert_data(database, table, newdata, column_names)





# def collect_case_workflow_inputs(case_data):
#     '''
#     (dict) -> list
    
#     Returns a list of dictionaries with the workflow input information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''
        
#     L = []
    
#     for d in case_data['workflow_runs']:
#         case = case_data['case']
#         donor = get_donor_name(case_data)
#         limskeys = d['limsIds'].split(',')
#         sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
#         #platform = d['platform']
        
#         platform = 'NA'
        
        
#         wfrun_id = d['wfrunid']
                
#         for limskey in sequencing_attributes:
#             D = {'case_id': case,
#                  'donor_id': donor,
#                  'library': sequencing_attributes[limskey]['library'],
#                  'project_id': sequencing_attributes[limskey]['project'],
#                  'barcode': sequencing_attributes[limskey]['barcode'],
#                  'platform': platform,
#                  'lane': sequencing_attributes[limskey]['lane'],
#                  'wfrun_id': wfrun_id,
#                  'run': sequencing_attributes[limskey]['run'],
#                  'limskey': limskey,
#                  'platform': platform}
            
#             if D not in L:
#                 L.append(D)
    
#     return L



# def add_workflow_inputs_to_db(database, provenance_data, cases_to_update, table = 'Workflow_Inputs'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts or updates workflow input information 
 
#     Parameters
#     ----------    
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Workflow_Inputs
#     '''

#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
     
#         for case_data in provenance_data:
#             case = case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 workflow_input_info = collect_case_workflow_inputs(case_data)
#                 for d in workflow_input_info:
#                     L = [d[i] for i in column_names]
#                     newdata.append(L)             

#         # add data
#         insert_data(database, table, newdata, column_names)


# def map_file_to_worklow(case_data):
#     '''
#     (dict) -> dict

#     Returns a dictionary of file swids matched to their workflow run id for a case    
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''
    
#     D = {}
       
#     for i in range(len(case_data['cerberus_data'])):
#         wfrun_id = case_data['cerberus_data'][i]['workflow_run_accession']
#         file_swid = case_data['cerberus_data'][i]['accession'] 
#         if file_swid in D:
#             assert D[file_swid] == wfrun_id
#         else:
#             D[file_swid] = wfrun_id
    
#     return D



# def get_workflow_inputs(case_data):
#     '''
#     (dict) -> dict

#     Returns a dictionary of workflows and their input files for a case    
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''
    
#     D = {}
       
#     for i in range(len(case_data['cerberus_data'])):
#         wfrun_id = case_data['cerberus_data'][i]['workflow_run_accession']
#         input_files = json.loads(case_data['cerberus_data'][i]['input_files']) 
#         if wfrun_id in D:
#             assert D[wfrun_id] == input_files
#         else:
#             D[wfrun_id] = input_files
    
#     return D


# def collect_workflow_relationships(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary of parent-children workflows
#     for all workflows for a given case

#     Paramaters
#     ----------
#     - case_data (dict): Dictionary with a single case data         
#     '''
    
#     D = {}
    
#     for d in case_data['workflow_runs']:
#         workflow = d['wfrunid']
#         parents = json.loads(d['parents'])
#         children = json.loads(d['children'])
#         if parents:
#             parent_workflows = [i[1] for i in parents]
#         else:
#             parent_workflows = ['NA']
#         if children:
#             children_workflows = [i[1] for i in children]
#         else:
#             children_workflows = ['NA']
        
#         # record the parent-child relationship of the current workflow
#         if workflow in D:
#             D[workflow].extend(children_workflows)
#         else:
#             D[workflow] = children_workflows
#         D[workflow] = list(set(D[workflow]))
#         # record the parent-child relationships of each parent and current workflow 
#         for workflow_run in parent_workflows:
#             if workflow_run in D:
#                 D[workflow_run].append(workflow)
#             else:
#                 D[workflow_run] = [workflow]
#             D[workflow_run] = list(set(D[workflow_run]))

#     return D



# def add_workflows_relationships_to_db(database, provenance_data, cases_to_update, table = 'Parents'):
#     '''
#     (str, str, str, dict, str, str, str, str) -> None
    
#     Inserts or updates workflow information and parent-children workflow relationships
 
#     Parameters
#     ----------   
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Parents
#     '''

#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))
        
#         # make a list of data to insert in the database
#         newdata = []
        
#         for case_data in provenance_data:
#             case = case_data['case']
#             donor = get_donor_name(case_data)
#             project = case_data['project_info'][0]['project']             
                       
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 parents =  collect_workflow_relationships(case_data)
#                 for workflow in parents:
#                     for child in parents[workflow]:
#                         L = [workflow, child, project, case, donor]
#                         if L not in newdata:
#                             newdata.append(L)
        
#         # add data
#         insert_data(database, table, newdata, column_names)


# def get_donor_external_id(case_data):
#     '''
#     (dict) -> str
    
#     Returns the external id of a donor by extracting the external id from the
#     cerberus data of the corresponding case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data    
#     '''

#     L = []

#     for i in range(len(case_data['cerberus_data'])):
#         L.append(case_data['cerberus_data'][i]['external_donor_id'])
#     L = list(set(L))
#     assert len(L) == 1
#     return L[0]


# def get_donor_sex(case_data):
#     '''
#     (dict) -> str
    
#     Returns the sex of a donor by extracting the sex from the
#     pinery data of a case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data    
#     '''
    
#     sex = []
#     for i in range(len(case_data['pinery_data'])):
#         sex.append(case_data['pinery_data'][i]['sex'])       
#     sex = map(lambda x: x.lower(), list(set(sex)))
    
#     # sex may not be defined for all samples.
#     # assign donor sex to a known sample sex
#     if 'female' in sex:
#         assert 'male' not in sex
#         sex = 'Female'
#     elif 'male' in sex:
#         assert 'female' not in sex
#         sex = 'Male'
#     else:
#         sex = 'Unknown'
        
#     return sex


# def get_donor_species(case_data):
#     '''
#     (dict) -> str
    
#     Returns the species of a donor by extracting the species from the
#     pinery data of a case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data    
#     '''
    
#     species = []
#     for i in range(len(case_data['pinery_data'])):
#         species.append(case_data['pinery_data'][i]['organism'])       
#     species = list(set(species))
#     assert len(species) == 1
        
#     return species[0]



# def get_sequencing_status(case_data):
#     '''
#     (dict) -> bool
    
#     Returns True if case sequencing is complete
        
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''
    
#     sequencing = json.loads(case_data['case_info']['sequencing'])
#     L = []
#     for d in sequencing:
#         if 'sequencing' in d['type'].lower():
#             L.append(d['complete'])
        
#     complete = all(L)
#     return complete


# def collect_case_sample_info(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with the sample information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''
        
#     D = {}
        
#     for d in case_data['sample_info']:
#         case = case_data['case']
#         donor = d['donor']
#         project_id = d['project']
#         assay = case_data['assay']
#         sequencing_status = get_sequencing_status(case_data)
#         external_id = d['externalId']
#         #sex = d['sex']
#         miso = 'NA'
#         species = d['organism']
        
#         d = {'case_id': case, 'assay': assay, 'donor_id': donor, 'ext_id': external_id, 
#              'species': species, 'project_id': project_id, 'miso': miso, 'sequencing_status': str(int(sequencing_status))}
        
#         if case not in D:
#             D[case] = d
#         else:
#             assert D[case] == d
            
#     return D



# def add_samples_info_to_db(database, provenance_data, cases_to_update, table = 'Samples'):
#     '''
#     (str, str, str, dict, dict, dict) -> None
    
#     Inserts samples data into Samples table of database    
    
#     Parameters
#     ----------
#     - database (str): Path to the databae file
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Name of table in database
#     '''
    
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
        
#         for case_data in provenance_data:
#             case = case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 sample_info = collect_case_sample_info(case_data)
#                 L = [sample_info[case][i] for i in column_names]
#                 newdata.append(L)
                
#         # add data
#         insert_data(database, table, newdata, column_names) 









# def get_file_signoff(project, nabu = 'https://nabu.gsi.oicr.on.ca/get-fileqcs'):
#     '''
#     (str, str) -> dict

#     Returns a dictionary 

#     Parameters
#     ----------
#     - project (str): Project of interest
#     - nabu (str): URL to access the file qc in Nabu
#     '''

#     headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
#     json_data = {'project': project}
#     response = requests.post(nabu, headers=headers, json=json_data)
    
#     D = {}
    
#     if response.status_code == 200:
#         for d in response.json()['fileqcs']:
#             fileid = d['fileid']
#             filepath = d['filepath']
#             if 'username' in d:
#                 username = d['username']
#             else:
#                 username = 'NA'
#             qcstatus = d['qcstatus']
#             if 'comment' in d:
#                 ticket = d['comment']
#             else:
#                 ticket = 'NA'
            
#             assert fileid not in D
#             D[fileid] = {
#                 'fileid' : fileid,
#                 'filepath' : filepath,
#                 'username' : username,
#                 'qcstatus' : qcstatus,
#                 'ticket' : ticket
#            }
            
#     return D   



# def get_case_file_qc(database, project, fileqc):
#     '''
#     (str, str, dict) -> dict
    
#     Returns a dictionary with file qc status from Nabu organized by workflow run id
#     for all cases in project
        
#     Parameters
#     ----------
#     - database (str): Path to the waterzooi sqlite database
#     - project (str): Project of interest
#     - fileqc (dict): Dictionary with file qc status extracted from nabu for a given project
#     '''
    
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT case_id, file_swid, wfrun_id FROM Files WHERE \
#                         project_id = ?;", (project,)).fetchall()
#     conn.close()
 
#     D = {}
    
#     for i in data:
#         case_id = i['case_id']
#         file_id = i['file_swid']
#         wfrun_id = i['wfrun_id']
        
#         # get the qc info
#         if file_id in fileqc:
#             username = fileqc[file_id]['username']
#             qcstatus = fileqc[file_id]['qcstatus']
#             filepath = fileqc[file_id]['filepath']
#             ticket = fileqc[file_id]['ticket']
#         else:
#             username = 'NA'
#             qcstatus = 'NA'
#             filepath = 'NA'
#             ticket = 'NA'

#         if case_id not in D:
#             D[case_id] = {}
#         if wfrun_id not in D[case_id]:
#             D[case_id][wfrun_id] = {}
#         D[case_id][wfrun_id][file_id] = {
#             'username': username,
#             'qcstatus': qcstatus,
#             'filepath': filepath,
#             'ticket': ticket}
            
#     return D                
            

# def list_projects_for_file_qc(provenance_data, cases_to_update):
#     '''
#     (str) -> list
    
#     Returns a list of project from the Projects table in database corresponding
#     to cases that need to be updated in the databse
    
#     Parameters
#     ----------
#     - database (str): Path to the waterzooi sqlite database
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     '''
    
#     L = []
#     for case_data in provenance_data:
#         case = case_data['case']
#         if case in cases_to_update:
#             for d in case_data['project_info']:
#                 L.append(d['project'])
    
#     L = list(set(L))
    
#     return L



# def add_file_qc(database, provenance_data, cases_to_update, table = 'File_qc'):
#     '''
#     (str, dict, str) -> None
    
    
#     Inserts file qc from Nabu information in database's File_qc table
    
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file qc information. Default is File_qc
#     '''
    
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))
               
#         # make a list of projects
#         projects = list_projects_for_file_qc(provenance_data, cases_to_update)
#         for project in projects:
#             # get the file qc for each project:
#             fileqc = get_file_signoff(project)
#             # organize file qc by case and workflow
#             case_file_qc = get_case_file_qc(database, project, fileqc)

#             # make a list of data to insert in the database
#             newdata = []
#             for case_id in case_file_qc:
#                 if case_id in cases_to_update:
#                     for wfrun_id in case_file_qc[case_id]:
#                         for file_id in case_file_qc[case_id][wfrun_id]:
#                             L = [project, case_id, wfrun_id, file_id, case_file_qc[case_id][wfrun_id][file_id]['filepath'],
#                                  case_file_qc[case_id][wfrun_id][file_id]['username'],
#                                  case_file_qc[case_id][wfrun_id][file_id]['ticket']]
#                             if case_file_qc[case_id][wfrun_id][file_id]['qcstatus'] == 'PASS':
#                                 L.append('1')
#                             elif case_file_qc[case_id][wfrun_id][file_id]['qcstatus'] == 'FAILED':
#                                 L.append('0')
#                             else:
#                                 L.append(case_file_qc[case_id][wfrun_id][file_id]['qcstatus'])
#                             newdata.append(L)    
                               
#             # insert data
#             insert_data(database, table, newdata, column_names)


# def generate_database(database, provenance_data_file):
#     '''
#     (str, str) -> None

#     Generates the waterzooi database using data in the provenance_data_file and 
#     calcontaqc_db database

#     Parameters
#     ----------
#     - database (str): Path of the waterzooi database
#     - provenance_data_file (str): Path to the provenance data json
#     '''
    
#     # create database if file doesn't exist
#     if os.path.isfile(database) == False:
#         initiate_db(database, 'waterzooi', ['Projects', 'Workflows', 'Parents',
#                     'Files', 'Libraries', 'Workflow_Inputs', 'Samples', 'Checksums', 'File_qc'])
#     print('initiated database')
    
#     # load data from file
#     provenance_data = load_data(provenance_data_file)
#     print('loaded data')
    
#     # remove cases with incomplete data
#     provenance_data = clean_up_data(provenance_data)
#     print('removed cases with incomplete data')    
    
#     # clean up workflows - remove children and parent workflows that are not in a given case        
#     for i in range(len(provenance_data)):
#         provenance_data[i] = clean_up_workflows(provenance_data[i])    
#     print('clean up workflows')
    
#     # collect the md5sum of each donor's data
#     md5sums = compute_case_md5sum(provenance_data)
    
       
#     print('computed md5sums')
#     # collect the recorded md5sums of the donor data from the database
#     recorded_md5sums = get_cases_md5sum(database, table = 'Checksums')
#     print('pulled md5sums from database')
#     # determine the donors that need an update
#     cases_to_update = cases_info_to_update(md5sums, recorded_md5sums)
#     print('determined donors to update')
    
#     # add project information
#     add_project_info_to_db(database, provenance_data, 'Projects')
#     print('added project info to database')
#     # add library information
#     add_library_info_to_db(database, provenance_data, cases_to_update, 'Libraries')
#     print('added library info to database')
#     # add sample information
#     add_samples_info_to_db(database, provenance_data, cases_to_update, 'Samples')
#     print('added sample information to database')
#     # add workflow inputs
#     add_workflow_inputs_to_db(database, provenance_data, cases_to_update, 'Workflow_Inputs')
#     print('added workflow inputs to database')
#     # add workflow relationships
#     add_workflows_relationships_to_db(database, provenance_data, cases_to_update, 'Parents')
#     print('added workflow relationships to database')
#     # add workflow information
#     add_workflows_to_db(database, provenance_data, cases_to_update, 'Workflows')
#     print('added workflow info to database')
#     # add file information
#     add_file_info_to_db(database, provenance_data, cases_to_update, 'Files')
#     print('added file info to database')
    
#     # add file qc information from Nabu
#     add_file_qc(database, provenance_data, cases_to_update, 'File_qc')
#     print('added File qc info to database')
       
#     # update the checksums for donors
#     add_checksums_info_to_db(database, cases_to_update, 'Checksums')
#     print('added md5sums to database')
    
 
   
 
#generate_database('waterzooi_db_case.db', 'CaseInfo_IRIS_NEOPOC.dump.json')    
     
#generate_database('waterzooi_db_case_small.db', 'provenance_reporter_small.json')    
 
#generate_database('waterzooi_db_case.db', 'provenance_reporter.json')    

    
 
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





##################################
##################################



# def add_checksums_info_to_db(database, cases_to_update, table = 'Checksums'):
#     '''
#     (str, dict, str) -> None
    
#     Update table Checksums with the checksum of the case info 
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Name of Table in database. Default is Checksums
#     '''
    
#     if cases_to_update:
#         delete_records(cases_to_update, database, table)
        
#         # make a list of data to insert
#         newdata = []
        
#         # connect to db
#         conn = sqlite3.connect(database, timeout=30)
#         # get column names
#         data = conn.execute("SELECT * FROM {0}".format(table))
#         column_names = [column[0] for column in data.description]

#         # order values according to column names
#         for i in cases_to_update:
#             L = [cases_to_update[i]['project_id'], i, cases_to_update[i]['donor_id'],  cases_to_update[i]['md5']]
#             newdata.append(L)
#         vals = '(' + ','.join(['?'] * len(newdata[0])) + ')'
#         conn.executemany('INSERT INTO {0} {1} VALUES {2}'.format(table, tuple(column_names), vals), newdata)
#         conn.commit()
#         conn.close()





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



# def collect_case_file_info(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with all the file information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''    
        
#     D = {}
    
#     for d in case_data['workflow_runs']:
#         donor = get_donor_name(case_data)
#         case = case_data['case']
#         project_id = case_data['project_info'][0]['project']
#         wfrun_id = d['wfrunid']
#         limskeys = d['limsIds'].split(',')
#         files = json.loads(d['files'])
#         for i in range(len(files)):
#             file_swid = files[i]['accession']
#             file_path = files[i]['path']
#             md5sum = files[i]['md5']
#             attributes = files[i]['file_attributes']
#             if attributes:
#                 attributes = json.loads(attributes)
#                 file_attributes = {}
#                 for k in attributes:
#                     file_attributes[k] = attributes[k][0]
#                 if len(file_attributes) == 0:
#                     file_attributes = ''
#                 else:
#                     file_attributes = json.dumps(file_attributes)
#             else:
#                 file_attributes = ''
#             creation_date = files[i]['timestamp'].replace('Z', '')
#             creation_date = creation_date.split('.')[0]
#             creation_date = ' '.join(creation_date.split('T'))
#             pattern = '%Y-%m-%d %H:%M:%S'
#             creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
            
            
#             assert file_swid not in D
#             D[file_swid] = {'file_swid': file_swid,
#                             'md5sum': md5sum,
#                             'creation_date': creation_date,
#                             'attributes': file_attributes,
#                             'file': file_path,
#                             'project_id': project_id,
#                             'case_id': case,
#                             'donor_id': donor,
#                             'wfrun_id': wfrun_id,
#                             'limskey': limskeys}  
                                
#     return D


# def add_file_info_to_db(database, provenance_data, cases_to_update, table = 'Files'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts file information in database's Files table
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Files
#     '''
        
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))
    
#         # make a list of data to insert in the database
#         newdata = []
            
#         for case_data in provenance_data:
#             donor = get_donor_name(case_data)
#             case =  case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 file_info = collect_case_file_info(case_data)
#                 for file_swid in file_info:
#                     for limskey in file_info[file_swid]['limskey']:
#                         L = [file_info[file_swid][i] for i in column_names[:column_names.index('limskey')]]
#                         L.append(limskey)
#                         L.extend([file_info[file_swid][i] for i in column_names[column_names.index('limskey')+1:]])
#                         newdata.append(L)             
    
#         # add data
#         insert_data(database, table, newdata, column_names)


def collect_case_library_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the library information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''    

    D = {}
    
    for d in case_data['sample_info']:
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


# def add_library_info_to_db(database, provenance_data, cases_to_update, table = 'Libraries'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts library information in the Libraries table of the database
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Libraries
#     '''                 
 
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for cases to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
    
#         for case_data in provenance_data:
#             case =  case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 library_info = collect_case_library_info(case_data)
#                 for library in library_info:
#                     for d in library_info[library]:
#                         L = [d[i] for i in column_names]
#                         newdata.append(L)             
  
#         # add data
#         insert_data(database, table, newdata, column_names)




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
        #project = case_data['project']
        project = case_data['project_info'][0]['project']
        case = case_data['case']
        assay = case_data['assay']
        #assert len(case_data['project_info']) == 1
        active = case_data['project_info'][0]['isActive']
        pipeline = case_data['project_info'][0]['pipeline']
        deliverables = case_data['project_info'][0]['deliverables']
        samples, library_types = [], []
        for d in case_data['sample_info']:
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


# def add_project_info_to_db(database, provenance_data, project_table = 'Projects'):
#     '''
#     (str, list, str) -> None
    
#     Add project information into Projects table of database
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     - project_table (str): Name of Table storing project information in database. Default is Projects
#     '''
    
#     # get the column names
#     column_names = define_column_names()[project_table]
    
#     # collect project information
#     project_info = collect_project_info(provenance_data)
    
#     # make a list of data to insert in the database
#     newdata = []
    
#     # connect to db
#     conn = connect_to_db(database)
       
#     for project in project_info:
#         last_updated = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
        
#         # get the library types
#         library_types = ','.join(sorted(project_info[project]['library_types']))
#         samples = len(project_info[project]['samples'])
#         cases = len(project_info[project]['cases'])
#         assays = ','.join(sorted(project_info[project]['assays']))
        
#         L = [project, project_info[project]['pipeline'], last_updated, str(cases),
#              str(samples), str(library_types), assays, project_info[project]['deliverables'],
#              project_info[project]['active']]
        
#         newdata.append(L)
#         # delete project entry
#         query = 'DELETE FROM {0} WHERE project_id = \"{1}\"'.format(project_table, project)
#         conn.execute(query)
#         conn.commit()

#     conn.close()
        
#     # add data
#     insert_data(database, project_table, newdata, column_names)



# def collect_workflow_info(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with all the workflow information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''

#     D = {}
    
#     for d in case_data['workflow_runs']:
#         project_id = case_data['project_info'][0]['project']
#         case = case_data['case']
#         donor = get_donor_name(case_data)
#         wfrun_id = d['wfrunid']
#         wf = d['wf']
#         wfv = d['wfv']
#         limskeys = d['limsIds'].split(',')
#         file_count = len(json.loads(d['files']))   
#         sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
#         lane_count = len([sequencing_attributes[i]['lane'] for i in sequencing_attributes])
         
#         data = {'project_id': project_id,
#                 'wfrun_id': wfrun_id,
#                 'wf': wf,
#                 'wfv': wfv,
#                 'case_id': case,
#                 'donor_id': donor,
#                 'file_count': file_count,
#                 'lane_count': lane_count}
                
#         assert wfrun_id not in D
#         D[wfrun_id] = data
                
#     return D            



# def add_workflows_to_db(database, provenance_data, cases_to_update, table = 'Workflows'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts or updates workflow information 
 
#     Parameters
#     ----------    
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Workflows
#     '''
    
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
        
#         for case_data in provenance_data:
#             case = case_data['case']
#             donor = get_donor_name(case_data)
                      
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 workflow_info = collect_workflow_info(case_data)                
#                 for workflow in workflow_info:
#                     L = [workflow_info[workflow][i] for i in column_names]
#                     newdata.append(L)             
      
#         # add data
#         insert_data(database, table, newdata, column_names)





# def collect_case_workflow_inputs(case_data):
#     '''
#     (dict) -> list
    
#     Returns a list of dictionaries with the workflow input information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''
        
#     L = []
    
#     for d in case_data['workflow_runs']:
#         case = case_data['case']
#         donor = get_donor_name(case_data)
#         limskeys = d['limsIds'].split(',')
#         sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
#         #platform = d['platform']
        
#         platform = 'NA'
        
        
#         wfrun_id = d['wfrunid']
                
#         for limskey in sequencing_attributes:
#             D = {'case_id': case,
#                  'donor_id': donor,
#                  'library': sequencing_attributes[limskey]['library'],
#                  'project_id': sequencing_attributes[limskey]['project'],
#                  'barcode': sequencing_attributes[limskey]['barcode'],
#                  'platform': platform,
#                  'lane': sequencing_attributes[limskey]['lane'],
#                  'wfrun_id': wfrun_id,
#                  'run': sequencing_attributes[limskey]['run'],
#                  'limskey': limskey,
#                  'platform': platform}
            
#             if D not in L:
#                 L.append(D)
    
#     return L



# def add_workflow_inputs_to_db(database, provenance_data, cases_to_update, table = 'Workflow_Inputs'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts or updates workflow input information 
 
#     Parameters
#     ----------    
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with donors for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Workflow_Inputs
#     '''

#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
     
#         for case_data in provenance_data:
#             case = case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 workflow_input_info = collect_case_workflow_inputs(case_data)
#                 for d in workflow_input_info:
#                     L = [d[i] for i in column_names]
#                     newdata.append(L)             

#         # add data
#         insert_data(database, table, newdata, column_names)


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


# def collect_workflow_relationships(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary of parent-children workflows
#     for all workflows for a given case

#     Paramaters
#     ----------
#     - case_data (dict): Dictionary with a single case data         
#     '''
    
#     D = {}
    
#     for d in case_data['workflow_runs']:
#         workflow = d['wfrunid']
#         parents = json.loads(d['parents'])
#         children = json.loads(d['children'])
#         if parents:
#             parent_workflows = [i[1] for i in parents]
#         else:
#             parent_workflows = ['NA']
#         if children:
#             children_workflows = [i[1] for i in children]
#         else:
#             children_workflows = ['NA']
        
#         # record the parent-child relationship of the current workflow
#         if workflow in D:
#             D[workflow].extend(children_workflows)
#         else:
#             D[workflow] = children_workflows
#         D[workflow] = list(set(D[workflow]))
#         # record the parent-child relationships of each parent and current workflow 
#         for workflow_run in parent_workflows:
#             if workflow_run in D:
#                 D[workflow_run].append(workflow)
#             else:
#                 D[workflow_run] = [workflow]
#             D[workflow_run] = list(set(D[workflow_run]))

#     return D



# def add_workflows_relationships_to_db(database, provenance_data, cases_to_update, table = 'Parents'):
#     '''
#     (str, str, str, dict, str, str, str, str) -> None
    
#     Inserts or updates workflow information and parent-children workflow relationships
 
#     Parameters
#     ----------   
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Parents
#     '''

#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))
        
#         # make a list of data to insert in the database
#         newdata = []
        
#         for case_data in provenance_data:
#             case = case_data['case']
#             donor = get_donor_name(case_data)
#             project = case_data['project_info'][0]['project']             
                       
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 parents =  collect_workflow_relationships(case_data)
#                 for workflow in parents:
#                     for child in parents[workflow]:
#                         L = [workflow, child, project, case, donor]
#                         if L not in newdata:
#                             newdata.append(L)
        
#         # add data
#         insert_data(database, table, newdata, column_names)


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



def get_sequencing_status(case_data):
    '''
    (dict) -> bool
    
    Returns True if case sequencing is complete
        
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''
    
    sequencing = json.loads(case_data['case_info']['sequencing'])
    L = []
    for d in sequencing:
        if 'sequencing' in d['type'].lower():
            L.append(d['complete'])
        
    complete = all(L)
    return complete


# def collect_case_sample_info(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with the sample information for a given case
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with a single case data   
#     '''
        
#     D = {}
        
#     for d in case_data['sample_info']:
#         case = case_data['case']
#         donor = d['donor']
#         project_id = d['project']
#         assay = case_data['assay']
#         sequencing_status = get_sequencing_status(case_data)
#         external_id = d['externalId']
#         #sex = d['sex']
#         miso = 'NA'
#         species = d['organism']
        
#         d = {'case_id': case, 'assay': assay, 'donor_id': donor, 'ext_id': external_id, 
#              'species': species, 'project_id': project_id, 'miso': miso, 'sequencing_status': str(int(sequencing_status))}
        
#         if case not in D:
#             D[case] = d
#         else:
#             assert D[case] == d
            
#     return D



# def add_samples_info_to_db(database, provenance_data, cases_to_update, table = 'Samples'):
#     '''
#     (str, str, str, dict, dict, dict) -> None
    
#     Inserts samples data into Samples table of database    
    
#     Parameters
#     ----------
#     - database (str): Path to the databae file
#     - provenance_data (list): List of dictionaries, each representing the data of a single case
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Name of table in database
#     '''
    
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
        
#         for case_data in provenance_data:
#             case = case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 sample_info = collect_case_sample_info(case_data)
#                 L = [sample_info[case][i] for i in column_names]
#                 newdata.append(L)
                
#         # add data
#         insert_data(database, table, newdata, column_names) 












            

# def list_projects_for_file_qc(provenance_data, cases_to_update):
#     '''
#     (str) -> list
    
#     Returns a list of project from the Projects table in database corresponding
#     to cases that need to be updated in the databse
    
#     Parameters
#     ----------
#     - database (str): Path to the waterzooi sqlite database
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     '''
    
#     L = []
#     for case_data in provenance_data:
#         case = case_data['case']
#         if case in cases_to_update:
#             for d in case_data['project_info']:
#                 L.append(d['project'])
    
#     L = list(set(L))
    
#     return L



# def add_file_qc(database, provenance_data, cases_to_update, table = 'File_qc'):
#     '''
#     (str, dict, str) -> None
    
    
#     Inserts file qc from Nabu information in database's File_qc table
    
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file qc information. Default is File_qc
#     '''
    
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for donors to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))
               
#         # make a list of projects
#         projects = list_projects_for_file_qc(provenance_data, cases_to_update)
#         for project in projects:
#             # get the file qc for each project:
#             fileqc = get_file_signoff(project)
#             # organize file qc by case and workflow
#             case_file_qc = get_case_file_qc(database, project, fileqc)

#             # make a list of data to insert in the database
#             newdata = []
#             for case_id in case_file_qc:
#                 if case_id in cases_to_update:
#                     for wfrun_id in case_file_qc[case_id]:
#                         for file_id in case_file_qc[case_id][wfrun_id]:
#                             L = [project, case_id, wfrun_id, file_id, case_file_qc[case_id][wfrun_id][file_id]['filepath'],
#                                  case_file_qc[case_id][wfrun_id][file_id]['username'],
#                                  case_file_qc[case_id][wfrun_id][file_id]['ticket']]
#                             if case_file_qc[case_id][wfrun_id][file_id]['qcstatus'] == 'PASS':
#                                 L.append('1')
#                             elif case_file_qc[case_id][wfrun_id][file_id]['qcstatus'] == 'FAILED':
#                                 L.append('0')
#                             else:
#                                 L.append(case_file_qc[case_id][wfrun_id][file_id]['qcstatus'])
#                             newdata.append(L)    
                               
#             # insert data
#             insert_data(database, table, newdata, column_names)



############################ PROCESS CASE BY CASE AND RECORD CASE INFO TO DB





def collect_case_project_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with case project information
        
    Parameters
    ----------
    - case_data (dict): Dictionary with case information 
    '''

    D = {}
    
    # there may be multiple projects associated  with a case
    for i in case_data['project_info']:
        project = i['project']
        case = case_data['case']
        assay = case_data['assay']
        active = i['isActive']
        pipeline = i['pipeline']
        deliverables = i['deliverables']
        samples, library_types = [], []
        for d in case_data['sample_info']:
           # check that samples correspond to the correct project
           if project == d['project']:
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



def update_project_info(project_info, case_project):
    '''
    (dict, dict) -> dict
    
    Updates the dictionary project_info with the project information of the case
    
    Parameters
    ----------
    - project_info (dict): Dictionary used to collect project level information
    - case_project (dict): Dictionary with case level project information
    '''
    
    for project in case_project:
        if project not in project_info:
            project_info[project] = case_project[project]
        else:
            project_info[project]['assays'].extend(case_project[project]['assays'])
            project_info[project]['cases'].extend(case_project[project]['cases'])
            project_info[project]['samples'].extend(case_project[project]['samples'])
            project_info[project]['library_types'].extend(case_project[project]['library_types'])
                                        
        project_info[project]['assays'] = list(set(project_info[project]['assays']))
        project_info[project]['cases'] = list(set(project_info[project]['cases']))
        project_info[project]['samples'] = list(set(project_info[project]['samples']))
        project_info[project]['library_types'] = list(set(project_info[project]['library_types']))
    
    return project_info
    
    
    
def collect_case_library_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the library information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''    

    D = {}
    
    for d in case_data['sample_info']:
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


def collect_case_sample_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with the sample information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''
        
    D = {}
        
    for d in case_data['sample_info']:
        case = case_data['case']
        donor = d['donor']
        project_id = d['project']
        assay = case_data['assay']
        sequencing_status = get_sequencing_status(case_data)
        external_id = d['externalId']
        #sex = d['sex']
        miso = 'NA'
        species = d['organism']
        
        d = {'case_id': case, 'assay': assay, 'donor_id': donor, 'ext_id': external_id, 
             'species': species, 'project_id': [project_id], 'miso': miso, 'sequencing_status': str(int(sequencing_status))}
        
        if case not in D:
            D[case] = d
        else:
            D[case]['project_id'].append(project_id)
        D[case]['project_id'] = list(set(D[case]['project_id']))
            
    return D



def get_intrument(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with the sequencing instrument for each limsId
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''
    
    D = {}
    
    for d in case_data['sample_info']:
        limskey = d['limsId']
        instrument = d['instrument']
        if limskey in D:
            assert D[limskey] == instrument
        else:
            D[limskey] = instrument
            
    return D
    

def collect_case_workflow_inputs(case_data):
    '''
    (dict) -> list
    
    Returns a list of dictionaries with the workflow input information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''

    # map instrument to each limskey
    instruments = get_intrument(case_data)
        
    L = []

    for d in case_data['workflow_runs']:
        case = case_data['case']
        donor = get_donor_name(case_data)
        limskeys = d['limsIds'].split(',')
        sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
        wfrun_id = d['wfrunid']
                
        for limskey in sequencing_attributes:
            D = {'case_id': case,
                 'donor_id': donor,
                 'library': sequencing_attributes[limskey]['library'],
                 'project_id': sequencing_attributes[limskey]['project'],
                 'barcode': sequencing_attributes[limskey]['barcode'],
                 'platform': instruments[limskey],
                 'lane': sequencing_attributes[limskey]['lane'],
                 'wfrun_id': wfrun_id,
                 'run': sequencing_attributes[limskey]['run'],
                 'limskey': limskey
                 }
            
            if D not in L:
                L.append(D)
    
    return L




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
    
    for i in case_data['project_info']:
        project = i['project']
        case_id = case_data['case']
        donor = get_donor_name(case_data) 
        for d in case_data['workflow_runs']:
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
                D[workflow]['children'].extend(children_workflows)
                D[workflow]['project'].append(project)
                assert D[workflow]['case'] == case_id and D[workflow]['donor'] == donor
            else:
                D[workflow] = {'children': children_workflows, 'case': case_id,
                               'donor': donor, 'project': [project]}
            D[workflow]['project'] = list(set(D[workflow]['project']))    
            D[workflow]['children'] = list(set(D[workflow]['children']))
            
            # record the parent-child relationships of each parent and current workflow 
            for workflow_run in parent_workflows:
                if workflow_run in D:
                    D[workflow_run]['children'].append(workflow)
                    D[workflow_run]['project'].append(project)
                    assert D[workflow_run]['case'] == case_id and D[workflow_run]['donor'] == donor
                else:
                    D[workflow_run] = {'children': [workflow], 'case': case_id,
                                       'donor': donor, 'project': [project]}
                D[workflow_run]['project'] = list(set(D[workflow_run]['project']))
                D[workflow_run]['children'] = list(set(D[workflow_run]['children']))

    return D







# def add_library_info_to_db(database, provenance_data, cases_to_update, table = 'Libraries'):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts library information in the Libraries table of the database
       
#     Parameters
#     ----------
#     - database (str): Path to the database file
#     - provenance_data (list): List of dictionaries, each representing the data of a single donor
#     - cases_to_update (dict): Dictionary with cases for which records needs to be updated
#     - table (str): Table in database storing file information. Default is Libraries
#     '''                 
 
#     if cases_to_update:
#         # get the column names
#         column_names = define_column_names()[table]
#         # remove rows for cases to update
#         delete_records(cases_to_update, database, table)
#         print('deleted records in {0}'.format(table))

#         # make a list of data to insert in the database
#         newdata = []
    
#         for case_data in provenance_data:
#             case =  case_data['case']
#             # check if donor needs to be updated
#             if case in cases_to_update and cases_to_update[case] != 'delete':
#                 library_info = collect_case_library_info(case_data)
#                 for library in library_info:
#                     for d in library_info[library]:
#                         L = [d[i] for i in column_names]
#                         newdata.append(L)             
  
#         # add data
#         insert_data(database, table, newdata, column_names)
    
    



def collect_workflow_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the workflow information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''

    D = {}
    
    for i in case_data['project_info']:
        project_id = i['project']
        for d in case_data['workflow_runs']:
            case_id = case_data['case']
            donor = get_donor_name(case_data)
            wfrun_id = d['wfrunid']
            wf = d['wf']
            wfv = d['wfv']
            limskeys = d['limsIds'].split(',')
            file_count = len(json.loads(d['files']))   
            sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
            lane_count = len([sequencing_attributes[i]['lane'] for i in sequencing_attributes])
         
            data = {'project_id': [project_id],
                    'wfrun_id': wfrun_id,
                    'wf': wf,
                    'wfv': wfv,
                    'case_id': case_id,
                    'donor_id': donor,
                    'file_count': file_count,
                    'lane_count': lane_count}
                
            if wfrun_id not in D:
                D[wfrun_id] = data
            else:
                D[wfrun_id]['project_id'].append(project_id)
            D[wfrun_id]['project_id'] = list(set(D[wfrun_id]['project_id']))    
    
    return D          




def collect_case_file_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the file information for a given case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''    
        
    D = {}
    
    for i in case_data['project_info']:
        project_id = i['project']
        for d in case_data['workflow_runs']:
            donor = get_donor_name(case_data)
            case_id = case_data['case']
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
            
                data = {'file_swid': file_swid,
                        'md5sum': md5sum,
                        'creation_date': creation_date,
                        'attributes': file_attributes,
                        'file': file_path,
                        'project_id': [project_id],
                        'case_id': case_id,
                        'donor_id': donor,
                        'wfrun_id': wfrun_id,
                        'limskey': limskeys}  
                                
                if file_swid not in D:
                    D[file_swid] = data
                else:
                    D[file_swid]['project_id'].append(project_id)
                D[file_swid]['project_id'] = list(set(D[file_swid]['project_id']))    
    
    return D



def get_file_signoff(project, nabu = 'https://nabu.gsi.oicr.on.ca/get-fileqcs'):
    '''
    (str, str) -> dict

    Returns a dictionary 

    Parameters
    ----------
    - project (str): Project of interest
    - nabu (str): URL to access the file qc in Nabu
    '''

    headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
    json_data = {'project': project}
    response = requests.post(nabu, headers=headers, json=json_data)
    
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



def get_project_file_qc(case_data, visited_projects, file_qc, nabu = 'https://nabu.gsi.oicr.on.ca/get-fileqcs'):
    '''
    (dict, list, dict, str) -> dict, list
    
    Returns a dictionary with the file QC information and a list of projects
    for which file QC has already been extracted
    
    Parameters
    ----------
    - case_data (dict): Dictionary with production data for a case
    - visited_projects (list): List of projects already evaluated
    - file_qc (dict): Dictionary recording the file QC information from nabu
    - nabu (str): Nabu file QC endpoint
    '''
    
    for i in case_data['project_info']:
        project = i['project']
        # check if file qc for that project have already been recorded
        if project not in visited_projects:
            # get the file qc for that project
            D = get_file_signoff(project, nabu)
            # update file qc with  of that projects
            file_qc.update(D)
            # update list of visited projects
            visited_projects.append(project)

    return file_qc, visited_projects


def collect_file_qc_info(case_data, file_qc):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with file QC information from Nabu
    
    Parameters
    ----------
    - case_data (dict): Dictionary with production data for a case
    - file_qc (dict): Dictionary with file QC extracted from Nabu
    '''

    D = {}

    for i in case_data['project_info']:
        project_id = i['project']
        for d in case_data['workflow_runs']:
            wfrun_id = d['wfrunid']
            case_id = case_data['case']
            files = json.loads(d['files'])
            for j in range(len(files)):
                file_swid = files[j]['accession']
                # get the qc info
                if file_swid in file_qc:
                    username = file_qc[file_swid]['username']
                    qcstatus = file_qc[file_swid]['qcstatus']
                    filepath = file_qc[file_swid]['filepath']
                    ticket = file_qc[file_swid]['ticket']
                else:
                    username = 'NA'
                    qcstatus = 'NA'
                    filepath = 'NA'
                    ticket = 'NA'
                         
                data = {'username': username, 'qcstatus': qcstatus,
                        'filepath': filepath, 'ticket': ticket, 'project_id': [project_id]}
                
                if case_id not in D:
                    D[case_id] = {}
                if wfrun_id not in D[case_id]:
                    D[case_id][wfrun_id] = {}
                if file_swid not in D[case_id][wfrun_id]:
                    D[case_id][wfrun_id][file_swid] = data
                else:
                    D[case_id][wfrun_id][file_swid]['project_id'].append(project_id)
                D[case_id][wfrun_id][file_swid]['project_id'] = list(set(D[case_id][wfrun_id][file_swid]['project_id']))
              
    return D                
                





def record_case_info(case_data, database, recorded_md5sums, project_info, file_qc, tables, columns):
    '''
    (dict, str, dict, dict, dict, dict) -> dict, dict
    
    
    
    Parameters
    ----------
    
    
    case_data
    
    
    , database
    
    
    , recorded_md5sums
    
    
    , project_info
    
    
    - file_qc (dict): Dictionary with File QC extracted from Nabu
    - tables (dict): Dictionary with database table names
    
    
    
    
    '''
    
    # compute the md5sum of the case info
    md5sum = compute_md5(case_data)
    # determine if case needs to be updated
    case_id = case_data['case']
    print(case_id)
    if case_to_update(recorded_md5sums, case_id, md5sum):
        print('case to update')
        # case info need to be updated
        # collect project info
        case_project = collect_case_project_info(case_data)
        print('collected case project info')
        # update project information
        project_info = update_project_info(project_info, case_project)
        print('updated project info')
        # collect library information
        library_info = collect_case_library_info(case_data)
        print('collected library info')
        # collect sample information
        sample_info = collect_case_sample_info(case_data)
        print('collected sample info')
        # collect workflow inputs
        workflow_inputs_info = collect_case_workflow_inputs(case_data)
        print('collected workflow inputs info')
        # collect workflow_relationships
        workflow_relationships = collect_workflow_relationships(case_data)
        print('collected workflow relationship info')
        # collect workflow information
        workflow_info = collect_workflow_info(case_data)
        print('collected workflow info')
        # collect file information
        file_info = collect_case_file_info(case_data)
        print('collected file info')
        # collect file qc information from Nabu
        case_file_qc = collect_file_qc_info(case_data, file_qc)
        
        #print(case_file_qc)
        
        
        
        print('collected file qc info')
           
        # insert data to tables
        # open database
        #conn = connect_to_db(database)

        
        # make a parallel list with table names
        table_names = [tables['workflows'], tables['files'], tables['files_qc'],
                       tables['libraries'], tables['samples'], tables['workflow_inputs'],
                       tables['parents'], tables['checksums']]
        # make a parallel list with column names
        column_names = [columns[i]['names'] for i in table_names]
        
        # make a list with data organized for insertion to the database
        L = [organize_workflow_info(workflow_info, column_names[0]),
             organize_file_info(file_info, column_names[1]),
             organize_file_qc_info(case_file_qc),
             organize_library_info(library_info, column_names[3]),
             organize_sample_info(sample_info, column_names[4]),
             organize_workflow_input_info(workflow_inputs_info, column_names[5]),
             organize_parent_info(workflow_relationships),
             organize_checksum_info(case_data, md5sum)]
           
        print('organized case data')

        # insert data to tables
        # open database
        conn = connect_to_db(database)

        for i in range(len(L)):
            # delete case records
            delete_unique_record(case_id, conn, database, table_names[i], 'case_id')
            # insert records
            insert_multiple_records(L[i], conn, database, table_names[i], column_names[i])

        print('inserted case data')


        # # remove case from recorded_md5sums
        # if case_id in recorded_md5sums:
        #     del recorded_md5sums[case_id]
        
        # close connection to database
        conn.close()
    
    # remove case from recorded_md5sums
    if case_id in recorded_md5sums:
        del recorded_md5sums[case_id]
    
    return project_info, recorded_md5sums
    





def organize_project_info(project_info):
    '''
    (dict) -> list
    
    Returns a list with project information to be addaded to the database
           
    Parameters
    ----------
    - project_info (dict): Dictionary with project information across all cases
    '''
    
    L = []
    
    for project in project_info:
        # include time the data was last updated
        last_updated = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
        # get the library types
        library_types = ','.join(sorted(project_info[project]['library_types']))
        samples = len(project_info[project]['samples'])
        cases = len(project_info[project]['cases'])
        assays = ','.join(sorted(project_info[project]['assays']))
        
        d = [project, project_info[project]['pipeline'], last_updated, str(cases),
             str(samples), str(library_types), assays, project_info[project]['deliverables'],
             project_info[project]['active']]
        
        L.append(d)
        
    return L



def organize_workflow_info(data, columns):
    '''
    (dict, list) -> list
    
    Returns a list with workflow information to be added to the waterzooi database
     
    Parameters
    ----------    
    - data (dict): Dictionary with case workflow information from production 
    - columns (list): List of column names in the workflow table of the waterzooi database
    '''
    
    L = []
    
    for wfrun_id in data:
        for project in data[wfrun_id]['project_id']:
            L.append([data[wfrun_id][i] if i != 'project_id' else project for i in columns])
                      
    return L



# def add_workflows_to_db(workflow_info, case, conn, database, table, columns, field):
#     '''
#     (str, list, dict, str) -> None
    
#     Inserts or updates workflow information 
 
#     Parameters
#     ----------    
#     - workflow_info (dict): Dictionary with case workflow information extracted from production    
#     - case (str): Case identifier
#     - conn (sqlite3.Connection): Open connection to the database
#     - database (str): Path to the waterzooi sqlite database
#     - table (str): Table storing workflow information in the database
#     - columns (list): List of column names in table
#     - field (str): Field in table
#     '''

#     # organize workflow information
#     L = organize_workflow_info(workflow_info, columns)
#     # delete case from workflow table
#     delete_unique_record(case, conn, database, table, field)
#     # insert records
#     insert_multiple_records(L, conn, database, table, columns)


def organize_file_info(data, columns):
    '''
    (dict, list) -> list
    
    Returns a list with file information to be added to the waterzooi database
     
    Parameters
    ----------    
    - data (dict): Dictionary with case file information from production 
    - columns (list): List of column names in the workflow table of the waterzooi database
    '''
    
    L = []

    for file_swid in data:
        for project in data[file_swid]['project_id']:
            for limskey in data[file_swid]['limskey']:
                d = [file_swid, project]
                d.extend([data[file_swid][i] for i in columns[columns.index('project_id')+1:columns.index('limskey')]])
                d.append(limskey)
                d.extend([data[file_swid][i] for i in columns[columns.index('limskey')+1:]])
                L.append(d)
    
        
    return L        
                         

# def add_files_to_db(file_info, case, conn, database, table, columns, field):
#     '''
#     (dict, str, sqlite3.Connection, str, str, list, str) -> None
    
#     Inserts or updates file information 
     
#     Parameters
#     ----------
#     - file_info (dict): Dictionary with case file information from production    
#     - case (str): Case identifier
#     - conn (sqlite3.Connection): Open connection to the database
#     - database (str): Path to the waterzooi sqlite database
#     - table (str): Table storing workflow information in the database
#     - columns (list): List of column names in table
#     - field (str): Field in table
#     '''
    
#     # organize file information
#     L = organize_file_info(file_info, columns)
#     # delete case from file table
#     delete_unique_record(case, conn, database, table, field)
#     # insert records
#     insert_multiple_records(L, conn, database, table, columns)



def organize_file_qc_info(data):
    '''
    (dict) -> list
    
    Returns a list with file QC information to be added to the waterzooi database
    
    Parameters
    ----------
    - data (dict): Dictionary with case file QC information from production and Nabu 
    '''
    
    L = []

    for case_id in data:
        for wfrun_id in data[case_id]:
            for file_id in data[case_id][wfrun_id]:
                for project in data[case_id][wfrun_id][file_id]['project_id']:
                    d = [project, case_id, wfrun_id, file_id, data[case_id][wfrun_id][file_id]['filepath'],
                         data[case_id][wfrun_id][file_id]['username'],
                         data[case_id][wfrun_id][file_id]['ticket']]
                    if data[case_id][wfrun_id][file_id]['qcstatus'] == 'PASS':
                        d.append('1')
                    elif data[case_id][wfrun_id][file_id]['qcstatus'] == 'FAILED':
                        d.append('0')
                    else:
                        d.append(data[case_id][wfrun_id][file_id]['qcstatus'])
                    L.append(d) 
                               
    return L


def organize_library_info(data, columns):
    '''
    (dict, list) -> list
    
    Returns a list with library information to be added to the waterzooi database
     
    Parameters
    ----------    
    - data (dict): Dictionary with library information from production 
    - columns (list): List of column names in the workflow table of the waterzooi database
    '''
    
    L = []
    
    for library in data:
        for d in data[library]:
            L.append([d[i] for i in columns])
         
    return L


def organize_sample_info(data, columns):
    '''
    (dict, list) -> list
    
    Returns a list with sample information to be added to the waterzooi database
     
    Parameters
    ----------    
    - data (dict): Dictionary with sample information from production 
    - columns (list): List of column names in the workflow table of the waterzooi database
    '''
    
    L = []
    
    for case_id in data:
        for project_id in data[case_id]['project_id']:
           d = [data[case_id][i] for i in columns[:columns.index('project_id')]]
           d.append(project_id)    
           d.extend([data[case_id][i] for i in columns[columns.index('project_id')+1:]])
           L.append(d)     
        
    return L
                        
     


def organize_workflow_input_info(data, columns):
    '''
    (dict, list) -> list
    
    Returns a list with workflow input information to be added to the waterzooi database
     
    Parameters
    ----------    
    - data (dict): Dictionary with workflow input information from production 
    - columns (list): List of column names in the workflow table of the waterzooi database
    '''
    
    L = []
    
    for d in data:
        L.append([d[i] for i in columns])
    
    return L


    
def organize_parent_info(data):
    '''
    (dict) -> list

    Returns a list with workflow relationships to be added to the waterzooi database
     
    Parameters
    ----------    
    - data (dict): Dictionary with workflow relationships information from production 
    '''
    
    
    L = []
    
    for parent in data:
        for project in data[parent]['project']:
            for child in data[parent]['children']:
                L.append([parent, child, project, data[parent]['case'], data[parent]['donor']]) 
                            
    return L        



def organize_checksum_info(case_data, md5sum):
    '''
    (dict, str) -> list

    Returns a list with case checksum to be added to the waterzooi database
     
    Parameters
    ----------    
    - case_data (dict): Dictionary of case data from production 
    - md5sum (list): md5sum of the case_data dictionary
    '''    

    L = []
    
    for i in case_data['project_info']:
        project = i['project']
        case_id = case_data['case']
        donor = get_donor_name(case_data)  
        L.append([project, case_id, donor, md5sum])        

    return L




def add_file_qc_to_db(file_qc_info, case, conn, database, table, columns, field):
    '''
    (dict, str, sqlite3.Connection, str, str, list, str) -> None
    
    Inserts or updates file information 
     
    Parameters
    ----------
    - file_qc_info (dict): Dictionary with case file QC information from production and Nabu
    - case (str): Case identifier
    - conn (sqlite3.Connection): Open connection to the database
    - database (str): Path to the waterzooi sqlite database
    - table (str): Table storing workflow information in the database
    - columns (list): List of column names in table
    - field (str): Field in table
    '''
    
    # organize file information
    L = organize_file_qc_info(file_qc_info)
    # delete case from file table
    delete_unique_record(case, conn, database, table, field)
    # insert records
    insert_multiple_records(L, conn, database, table, columns)


def add_project_info(database, project_info, table='Projects', field = 'project_id', database_name = 'waterzooi'):
    '''
    (str, dict, str, str, str) -> None 
    
    Inserts project information into the Projects table of the waterzooi database
    
    Paramaters
    ----------
    - database (str): Path to the waterzooi sqlite database
    - project_info (dict): Project information collected from all cases in production
    - table (str): Table name in database
    - field (str): Field in database table 
    - database_name (str): Name of the database. Default is waterzooi
    '''
    
    # record project information in database
    data = organize_project_info(project_info)
    # connect to database
    conn = connect_to_db(database)
    # delete records
    delete_multiple_records(list(project_info.keys()), conn, database, table, field)
    # get the columns
    columns = define_columns(database_name)[table]['names']
    # add data to table
    insert_multiple_records(data, conn, database, table, columns)
    # close connection to database 
    conn.close()




def generate_database(database, provenance_data_file, nabu = 'https://nabu.gsi.oicr.on.ca/get-fileqcs'):
    '''
    (str, str) -> None

    Generates the waterzooi database using data in the provenance_data_file and 
    calcontaqc_db database

    Parameters
    ----------
    - database (str): Path of the waterzooi database
    - provenance_data_file (str): Path to the provenance data json
    '''
    
    tables = {'workflows': 'Workflows', 'parents': 'Parents', 
              'files': 'Files', 'files_qc': 'File_qc', 
              'libraries': 'Libraries', 'workflow_inputs': 'Workflow_Inputs',
              'checksums': 'Checksums', 'samples':'Samples'}
       
    columns = define_columns('waterzooi')
    
    # create database if file doesn't exist
    if os.path.isfile(database) == False:
        initiate_db(database, 'waterzooi', list(tables.values()) + ['Projects'])
    print('initiated database')
    
    # load data from file
    provenance_data = load_data(provenance_data_file)
    print('loaded data')
    
    # collect the recorded md5sums of the donor data from the database
    recorded_md5sums = get_cases_md5sum(database, table = 'Checksums')
    print('collected recorded checksums')
    
    # project information is updated once all the cases have been evaluated
    # initiate dictionary to collect project level info
    project_info = {}
    
    # initiate file qc 
    file_qc = {}
    # initiate collector to check for which projects file qc have been recorded
    visited_projects = []

    total, processed = len(provenance_data), 0
    
    print(total, processed)
    
    for case_data in provenance_data:
        total -= 1
        print(total)
        # check that no information is missing
        if is_case_info_complete(case_data):
            processed += 1    
            # remove workflows that do not belong to the case
            case_data = clean_up_workflows(case_data)
            # collect file qc at the project level if not already recorded
            file_qc, visited_projects = get_project_file_qc(case_data, visited_projects, file_qc, nabu)
            # record case data, update the project level information
            # remove case from recorded md5sums. any remaining cases are not in production
            # and should be removed from the database
            project_info, recorded_md5sums = record_case_info(case_data, database, recorded_md5sums, project_info, file_qc, tables, columns)
 
            print(total, processed)
            print('---------')
 
    # record project information in database
    if project_info:
        add_project_info(database, project_info, 'Projects', 'project_id', 'waterzooi')
    
    # delete records 
    if recorded_md5sums:
        # remove all recorded cases not in production    
        conn = connect_to_db(database)
        for i in tables:
            delete_multiple_records(list(recorded_md5sums.keys()), conn, database, tables[i], 'case_id')
        conn.close()
    
       
    
    # # collect the md5sum of each donor's data
    # md5sums = compute_case_md5sum(provenance_data)
    
       
    # print('computed md5sums')
    # print('pulled md5sums from database')
    # # determine the donors that need an update
    # cases_to_update = cases_info_to_update(md5sums, recorded_md5sums)
    # print('determined donors to update')
    
    # # add project information
    # add_project_info_to_db(database, provenance_data, 'Projects')
    # print('added project info to database')
    # # add library information
    # add_library_info_to_db(database, provenance_data, cases_to_update, 'Libraries')
    # print('added library info to database')
    # # add sample information
    # add_samples_info_to_db(database, provenance_data, cases_to_update, 'Samples')
    # print('added sample information to database')
    # # add workflow inputs
    # add_workflow_inputs_to_db(database, provenance_data, cases_to_update, 'Workflow_Inputs')
    # print('added workflow inputs to database')
    # # add workflow relationships
    # add_workflows_relationships_to_db(database, provenance_data, cases_to_update, 'Parents')
    # print('added workflow relationships to database')
    # # add workflow information
    # add_workflows_to_db(database, provenance_data, cases_to_update, 'Workflows')
    # print('added workflow info to database')
    # # add file information
    # add_file_info_to_db(database, provenance_data, cases_to_update, 'Files')
    # print('added file info to database')
    
    # # add file qc information from Nabu
    # add_file_qc(database, provenance_data, cases_to_update, 'File_qc')
    # print('added File qc info to database')
       
    # # update the checksums for donors
    # add_checksums_info_to_db(database, cases_to_update, 'Checksums')
    # print('added md5sums to database')



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



generate_database('waterzooi_db_case_new.db', 'provenance_reporter.json')

