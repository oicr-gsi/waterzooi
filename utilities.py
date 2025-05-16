# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:04:02 2023

@author: rjovelin
"""

import sqlite3
import time
import string
import random




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


# def get_children_workflows(project_name, database):
#     '''
#     (str) -> list
    
#     Returns a dictionary with workflow name, list of workflow_ids that are all children of 
#     workflow_id (i.e immediate downstream workflow) for a given project_name
    
#     Parameters
#     ----------
#     - project_name (str): Name of project of interest
#     - bmpp_id (str): bamMergePreprocessing workflow id 
#     - database (str): Path to the sqlite database
#     '''
    
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT DISTINCT Workflows.wf, Parents.parents_id, \
#                         Parents.children_id FROM Parents JOIN Workflows \
#                         WHERE Parents.project_id = ? \
#                         AND Workflows.project_id = ? AND \
#                         Workflows.wfrun_id = Parents.children_id;", (project_name, project_name)).fetchall()
#     data= list(set(data))
#     conn.close()
    
#     D = {}
#     for i in data:
#         if i['parents_id'] not in D:
#             D[i['parents_id']] = []
#         D[i['parents_id']].append({'wf': i['wf'], 'children_id': i['children_id']})
    
#     return D


# def get_workflow_names(project_name, database):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with workflow_id and workflow name key, value pairs
    
#     Parameters
#     ----------
#     - project_name (str): Name of project of interest
#     - database (str): Path to the sqlite database
#     '''

#     conn = connect_to_db(database)
#     data = conn.execute("SELECT DISTINCT Workflows.wfrun_id, Workflows.wf, Workflows.wfv FROM \
#                         Workflows WHERE Workflows.project_id = ?;", (project_name,)).fetchall()
#     data= list(set(data))
#     conn.close()
    
#     D = {}
#     for i in data:
#         D[i['wfrun_id']] = [i['wf'], i['wfv']]
       
#     return D


def remove_non_analysis_workflows(L):
    '''
    (list) -> list
    
    Returns a list L of dictionaries with workflows, removing any non-analysis workflows
    
    Parameters
    ----------
    - L (list): List of dictionaries with workflow name and workflow ids
    '''
    
    non_analysis_workflows = ('wgsmetrics', 'insertsizemetrics', 'bamqc', 'calculatecontamination',
                              'callability', 'fastqc', 'crosscheckfingerprintscollector',
                              'fingerprintcollector', 'bwamem', 'bammergepreprocessing',
                              'ichorcna', 'tmbanalysis', 'casava', 'bcl2fastq',
                              'fileimportforanalysis', 'fileimport', 'import_fastq',
                              'dnaseqqc', 'hotspotfingerprintcollector', 'rnaseqqc')
    
    to_remove = [i for i in L if i['wf'].split('_')[0].lower() in non_analysis_workflows or i['wf'].lower() in non_analysis_workflows]
    for i in to_remove:
        L.remove(i)
    
    return L



def convert_epoch_time(epoch):
    '''
    (str) -> str
    
    Returns epoch time in readable format
    
    Parameters
    ----------
    - epoch (str)
    '''
    
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(epoch)))


# def get_miso_sample_link(project_name, case, database):
#     '''
#     (str, str, str) -> str
    
#     Returns a link to the sample MISO page
    
#     Parameters
#     ----------
#     - project_name (str): Project of interest
#     - case (str): Sample name
#     - database (str): Path to the sqlite database
#     '''
    
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT miso FROM Samples WHERE project_id = ? AND case_id = ?;", (project_name, case)).fetchall()
#     data = list(set(data))
#     miso_link = data[0]['miso']
    
#     return miso_link



def get_library_design(library_source):
    '''
    (str) -> str
    
    Returns the description of library_source as defined in MISO
    
    Parameters
    ----------
    - library_source (str): Code of the library source as defined in MISO
    '''

    library_design = {'WT': 'Whole Transcriptome', 'WG': 'Whole Genome', 'TS': 'Targeted Sequencing',
                      'TR': 'Total RNA', 'SW': 'Shallow Whole Genome', 'SM': 'smRNA', 'SC': 'Single Cell',
                      'NN': 'Unknown', 'MR': 'mRNA', 'EX': 'Exome', 'CT': 'ctDNA', 'CM': 'cfMEDIP',
                      'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq'}

    if library_source in library_design:
        return library_design[library_source]
    else:
        return None





# def get_pipelines(project_name, database):
#     '''
#     (str, str) -> list
    
#     Returns a list of pipeline names based on the library codes extracted from database for project_name
    
#     Parameters
#     ----------
#     - project_name (str) Name of the project of interest
#     - database (str): Path to the sqlite database
#     '''    
    
#     # connect to db
#     conn = connect_to_db(database)
#     # extract library source
#     library_source = conn.execute("SELECT DISTINCT library_type FROM Files WHERE project_id = ?;", (project_name,)).fetchall()
#     library_source = list(set([i['library_type'] for i in  list(set(library_source))]))
#     # get the library definitions
#     pipelines = [get_library_design(j) for j in library_source if get_library_design(j)]
#     conn.close()
    
#     return pipelines



def get_donors(project_name, database):
    '''
    (str, str) -> list
    
    Returns a list of donors for a given project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
    
    # connect to db
    conn = connect_to_db(database)
    # extract library source
    data = conn.execute("SELECT DISTINCT case_id FROM Samples WHERE project_id = ?;", (project_name,)).fetchall()
    donors = list(set([i['case_id'] for i in  data]))
    conn.close()
    
    return donors



def secret_key_generator(size=10):
    '''
    (int)
    
    Returns a random string of length size with upper and lower case characters
    and digit
    
    Parameters
    ----------
    - size (int): Length of the random string
    '''
    
    chars=string.ascii_uppercase + string.ascii_lowercase + string.digits
    s = ''.join(random.choice(chars) for i in range(size))
    
    return s
    






# def get_assays(project_name, database):
#     '''
#     (str, str) -> list
    
#     Returns a list of assays extracted from the cases in the Samples table of the database
    
#     Parameters
#     ----------
#     - project_name (str): Name of the project of interest
#     - database (str): Path to the sqlite database
#     '''
    
#     # connect to db
#     conn = connect_to_db(database)
#     # extract library source
#     data = conn.execute("SELECT DISTINCT case_id FROM Samples WHERE project_id = ?;", (project_name,)).fetchall()
#     conn.close()
#     assays = sorted(list(set(map(lambda x: x.split(':')[1], [i['case_id'] for i in  data]))))
        
#     return assays

    
def get_case_md5sums(database, project_name):
    '''
    (str, str) -> dict
    
    Returns a dictionary with case, md5sums from the waterzooi database
    
    Parameters
    ----------
    - database (str): Path to the waterzooi sqlite database
    - project_name (str): Name of the project of interest
    '''
    
    # connect to db
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT case_id, md5 FROM Checksums WHERE project_id = ?;", (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        case = i['case_id']
        md5sum = i['md5']
        D[case] = md5sum
    
    return D
    

