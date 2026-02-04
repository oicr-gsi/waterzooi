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
    
    Returns a list sequence file information by grouping paired fastqs    
    Pre-condition: all fastqs are paired-fastqs. Non-paired-fastqs are discarded.
    
    
    
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    # get sequences    
    conn = connect_to_db(database)
    cmd = "SELECT DISTINCT Files.attributes, Files.case_id, Files.donor_id, \
          Files.file, Files.wfrun_id, Files.limskey, Samples.ext_id, Libraries.library, \
          Libraries.library_type, Libraries.tissue_origin, Libraries.sample_id,\
          Libraries.group_id, Libraries.group_id_description , Libraries.tissue_type, \
          Workflows.wf, Workflow_Inputs.platform, Workflow_Inputs.run, \
          Workflow_Inputs.lane FROM Files JOIN Samples JOIN Libraries JOIN \
          Workflows JOIN Workflow_Inputs WHERE Files.project_id = ? \
          AND Workflow_Inputs.project_id = ? AND Libraries.project_id = ? \
          AND Workflows.project_id = ? AND Samples.project_id = ? \
          AND Files.wfrun_id = Workflow_Inputs.wfrun_id \
          AND Files.wfrun_id = Workflows.wfrun_id AND \
          Libraries.case_id = Samples.case_id \
          AND Samples.donor_id = Files.donor_id \
          AND Workflow_Inputs.library = Libraries.library \
          AND Libraries.case_id = Files.case_id AND Libraries.case_id = Workflows.case_id \
          AND Libraries.case_id = Workflow_Inputs.case_id AND LOWER(Workflows.wf) in ('casava', 'bcl2fastq', 'fileimportforanalysis',\
          'fileimport', 'import_fastq');"    
    
    data = conn.execute(cmd, (project_name, project_name, project_name, project_name, project_name)).fetchall()
    
    conn.close()

    F = []
    
    for i in range(len(data)):
        # keep only read1
        if json.loads(data[i]['attributes'])['read_number'] == '1':
            case_id = data[i]['case_id']
            donor = data[i]['donor_id']
            sample = data[i]['ext_id']
            library =  data[i]['library']
            library_type =  data[i]['library_type']
            tissue_origin =  data[i]['tissue_origin']
            tissue_type =  data[i]['tissue_type']
            limskey = data[i]['limskey']
            group_id = data[i]['group_id']
            group_description = data[i]['group_id_description']
            workflow = data[i]['wf']
            wfrun = data[i]['wfrun_id']
            file = data[i]['file']
            run = data[i]['run'] + '_' + str(data[i]['lane'])
            platform = data[i]['platform']
            read_count = json.loads(data[i]['attributes'])['read_count'] if 'read_count' in json.loads(data[i]['attributes']) else 'NA' 
            readcount = '{:,}'.format(int(read_count)) if read_count != 'NA' else 'NA'
            sample_id = data[i]['sample_id']
            fileprefix = os.path.basename(file)
            fileprefix = '_'.join(fileprefix.split('_')[:-1])
            d = {'case_id': case_id, 'donor':donor, 'sample': sample, 'sample_id': sample_id, 'library': library, 'run': run,
                 'read_count': readcount, 'workflow': workflow, 'prefix':fileprefix,
                 'platform': platform, 'group_id': group_id,
                 'group_description': group_description, 'tissue_type': tissue_type,
                 'library_type': library_type, 'tissue_origin': tissue_origin,
                 'limskey': limskey}
            F.append(d)
       
    F.sort(key=lambda x: (x['case_id'], x['donor'], x['limskey'], x['sample_id'], x['library'], x['platform']))
        
    return F


def get_platform_shortname(project_name, database):
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
    cmd = "SELECT DISTINCT Workflow_Inputs.platform FROM Workflow_Inputs WHERE \
          Workflow_Inputs.project_id = ?;"
    data = conn.execute(cmd, (project_name,)).fetchall()
    conn.close()

    D = {}
    
    for i in data:
        instrument = ''
        platform = i['platform']
        if '_' in platform:
            s = platform.split('_')
        else:
            s = platform.split()
        for k in s:
            if 'seq' in k.lower():
                instrument = k
                break
        
        D[platform] = instrument.lower()
    
    return D


