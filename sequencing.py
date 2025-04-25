# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 12:10:24 2023

@author: rjovelin
"""


import json
import os
from utilities import connect_to_db



def collect_sequence_info(project_name, database):
    '''
    (str, str) -> list
    
    Returns a list with sequence file information for a project of interest
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    # get sequences    
    conn = connect_to_db(database)
    files = conn.execute("SELECT Files.file, Files.workflow, Files.version, Files.wfrun_id, Files.attributes, \
                         Workflow_Inputs.run, Workflow_Inputs.lane, Workflow_Inputs.platform, \
                         Libraries.library, Libraries.case_id, Libraries.donor_id, Libraries.ext_id, \
                         Libraries.group_id, Libraries.group_id_description, \
                         Libraries.library_type, Libraries.tissue_origin, Libraries.tissue_type \
                         from Files JOIN Workflow_Inputs JOIN Libraries WHERE Files.project_id = '{0}' \
                         AND Workflow_Inputs.project_id = '{0}' AND Libraries.project_id = '{0}' \
                         AND Files.wfrun_id = Workflow_Inputs.wfrun_id AND Workflow_Inputs.library = Libraries.library \
                         AND Libraries.case_id = Files.case_id \
                         AND LOWER(Files.workflow) in ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq');".format(project_name)).fetchall()
    conn.close()

    return files



def get_sequences(L):
    '''
    (list) -> list

    Returns a list sequence file information by grouping paired fastqs    
    Pre-condition: all fastqs are paired-fastqs. Non-paired-fastqs are discarded.
    
    Parameters
    ----------
    - L (list): List of sqlite3.Row extracted from the database and containing sequence file information
    '''
    
    # sort list according to files
    L.sort(key = lambda x: x['file'])
    
    F = []
    
    for i in range(len(L)):
        # keep only read1
        if json.loads(L[i]['attributes'])['read_number'] == '1':
            case = L[i]['case_id']
            donor = L[i]['donor_id']
            sample = L[i]['ext_id']
            library =  L[i]['library']
            library_type =  L[i]['library_type']
            tissue_origin =  L[i]['tissue_origin']
            tissue_type =  L[i]['tissue_type']
            group_id = L[i]['group_id']
            group_description = L[i]['group_id_description']
            workflow = L[i]['workflow'] + '_' + L[i]['version']
            wfrun = L[i]['wfrun_id']
            file = L[i]['file']
            run = L[i]['run'] + '_' + str(L[i]['lane'])
            platform = L[i]['platform']
            read_count = json.loads(L[i]['attributes'])['read_count'] if 'read_count' in json.loads(L[i]['attributes']) else 'NA' 
            sample_id = '_'.join([donor, tissue_origin, tissue_type, group_id]) 
            readcount = '{:,}'.format(int(read_count)) if read_count != 'NA' else 'NA'
            fileprefix = os.path.basename(file)
            fileprefix = '_'.join(fileprefix.split('_')[:-1])
            d = {'case': case, 'donor':donor, 'sample': sample, 'sample_id': sample_id, 'library': library, 'run': run,
                 'read_count': readcount, 'workflow': workflow, 'prefix':fileprefix,
                 'platform': platform, 'group_id': group_id,
                 'group_description': group_description, 'tissue_type': tissue_type,
                 'library_type': library_type, 'tissue_origin': tissue_origin}
            F.append(d)
       
    F.sort(key = lambda x: x['case'])
     
    return F


def platform_name(project_name, database):
    '''
    (str, str) -> list
    
    Returns a dictionary with sequencing platform, shortname for all platforms 
    for the project of interest
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    # get sequences    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Workflow_Inputs.platform FROM Workflow_Inputs WHERE \
                         Workflow_Inputs.project_id = '{0}';".format(project_name)).fetchall()
    conn.close()

    D = {}
    
    for i in data:
        s = ''
        for j in i['platform']:
            if not j.isnumeric():
                s += j
        s = s.split('_')
        while '' in s:
            s.remove('')
        D[i['platform']] = s[-1].lower()
    
    return D


