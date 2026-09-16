# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 19:17:48 2026

@author: rjovelin
"""


from commons import connect_to_db 

import json
import os


def get_project_info(database, project_name=None):
    '''
    (str, str | None) -> list
    
    Returns a list with project information extracted from database for all projects 
    of for a single project if project_name is defined 
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - project_name (None | str): Project of interest
    '''
    
    # connect to db
    conn = connect_to_db(database)
    if project_name:
        # extract project info
        project = conn.execute('SELECT * FROM Projects WHERE project_id=?', (project_name,)).fetchall()
    else:
        project = conn.execute('SELECT * FROM Projects').fetchall()
    conn.close()
    
    return project



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


def get_fileqc(nabu_cache, project_name):
    '''
    (str, str) -> dict
    
    Returns a dictionary with the file qc status for all files in project
    
    Parameters
    ----------
    - nabu_cache (str): Path to the nabu cache 
    - project_name (str): Name of the project of interest
    '''
    
    conn = connect_to_db(nabu_cache)
    data = conn.execute("SELECT DISTINCT case_id, fileid, username, qcstatus, ticket FROM fileqc WHERE project_id = ?;", (project_name,)).fetchall()
    conn.close()
    
    D = {}
    
    for i in data:
        case_id = i['case_id']
        username = i['username']
        ticket = i['ticket']
        status = i['qcstatus']
        if status.lower() == 'pass':
            qcstatus = 1
        else:
            qcstatus = 0
        file_swid = i['fileid']
        D[file_swid] = {'case_id': case_id, 'username': username,
                        'ticket': ticket, 'qcstatus': qcstatus}
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
    cmd = "SELECT DISTINCT Files.attributes, Files.case_id, Files.donor_id, Files.file_swid, \
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

    D = {}
      
    for i in range(len(data)):
        # group sequences based on wfrunid
        # get file file swids for each file
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
        wfrun_id = data[i]['wfrun_id']
        file = data[i]['file']
        run = data[i]['run'] + '_' + str(data[i]['lane'])
        platform = data[i]['platform']
        read_count = json.loads(data[i]['attributes'])['read_count'] if 'read_count' in json.loads(data[i]['attributes']) else 'NA' 
        readcount = '{:,}'.format(int(read_count)) if read_count != 'NA' else 'NA'
        sample_id = data[i]['sample_id']
        fileprefix = os.path.basename(file)
        fileprefix = '_'.join(fileprefix.split('_')[:-1])
        file_swid = data[i]['file_swid']    
        
        d = {'case_id': case_id, 'donor':donor, 'sample': sample, 'sample_id': sample_id, 'library': library, 'run': run,
                 'read_count': readcount, 'workflow': workflow, 'prefix':fileprefix,
                 'platform': platform, 'group_id': group_id, 'wfrun_id': wfrun_id,
                 'group_description': group_description, 'tissue_type': tissue_type,
                 'library_type': library_type, 'tissue_origin': tissue_origin,
                 'limskey': limskey}
        
        if wfrun_id not in D:
            D[wfrun_id] = d
            D[wfrun_id]['file_swids'] = [file_swid]
        else:
            D[wfrun_id]['file_swids'].append(file_swid)
        
        #F.append(d)
    
        
    F = []
    for wfrunid in D:
        F.append(D[wfrunid])
    F.sort(key=lambda x: (x['case_id'], x['donor'], x['limskey'], x['sample_id'], x['library'], x['platform']))
    
    return F



def get_workflow_level_release(L, fileqc):
    '''
    (list, dict) -> bool
    
    Returns True if any file in L has been released
    
    Parameters
    ----------
    - L (list): List of file swids of a same workflow
    - fileqc (dict): Dictionary with file qc extracted from the nabu cache
    '''
    
    qc = []
    for fileid in L:
        if fileid in fileqc:
            qc.append(fileqc[fileid]['qcstatus'])
        else:
            qc.append(0)
      
    return any(qc)


def merge_qc_status_workflow(L, fileqc):
    '''
    
    
    
    '''
    
    
    username = []
    ticket = []
    qcstatus = []
       
    for fileid in L:
        if fileid in fileqc:
            username.append(fileqc[fileid]['username'])
            ticket.extend(fileqc[fileid]['ticket'])
            qcstatus.append(fileqc[fileid]['qcstatus'])
        else:
            username.append('NA')
            ticket.append('NA')
            qcstatus.append(0)
    
    while 'NA' in username:
        username.remove('NA')
    while 'NA' in ticket:
        ticket.remove('NA')
    qcstatus = any(qcstatus)
    
    return {'username': username, 'ticket': ticket, 'qcstatus': qcstatus}
    
       
    
def files_to_cases(database, project_name):
    '''


    '''


    # connect to db
    conn = connect_to_db(database)
    data = conn.execute('SELECT file_swid, case_id FROM Files WHERE project_id=?', (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        D[i['file_swid']] = i['case_id']
    
    return D






    
    
    