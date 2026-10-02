# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 19:17:48 2026

@author: rjovelin
"""


from commons import connect_to_db 

import json
import os
import time
import string
import random
import networkx as nx
import plotly.graph_objects as go


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


def get_case_release_signoff(nabu_cache, case_id, project_name):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with the release and release approval signoff
    for case_id in project
    
    Parameters
    ----------
    - nabu_cache (str): Path to the nabu cache 
    -case_id (str): Case identifier
    - project_name (str): Name of the project of interest
    '''
    
    conn = connect_to_db(nabu_cache)
    data = conn.execute("SELECT DISTINCT case_id, release, release_approval FROM signoff WHERE case_id = ? AND project_id = ?;", (case_id, project_name,)).fetchall()
    conn.close()
    
    D = {}
    
    for i in data:
        release = json.loads(i['release'])
        approval = json.loads(i['release_approval'])
        D[case_id] = {'release': release, 'release_approval': approval}
        
    return D






def get_fileqc(nabu_cache, project_name, case_id = None):
    '''
    (str, str, str | None) -> dict
    
    Returns a dictionary with the file qc status for all files in project
    
    Parameters
    ----------
    - nabu_cache (str): Path to the nabu cache 
    - project_name (str): Name of the project of interest
    - case_id (str | None): Optional case identifier
    '''
    
    conn = connect_to_db(nabu_cache)
    
    if case_id:
        data = conn.execute("SELECT DISTINCT case_id, fileid, username, qcstatus, ticket FROM fileqc WHERE case_id = ? AND project_id = ?;", (case_id, project_name,)).fetchall()
    else:
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
            ticket.append(fileqc[fileid]['ticket'])
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
    ticket = list(set(ticket))
    username = list(set(username))
        
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



def get_assays(database, project_name):
    '''
    (str, str) -> list
    
    Returns a list of all assay names in project (without the assay version)
        
    Parameters
    ----------
    - database (str): Path to the database
    - project_name (str): Name of project of interest
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT assays FROM Projects WHERE project_id = ?;", (project_name,)).fetchall()     
    conn.close()
    
    assays = []
    for i in data:
        assays.extend(i['assays'].split(','))
    assays = sorted(list(set((map(lambda x: '_'.join(x.split('_')[:-1]), assays))))) 
        
    return assays


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

    
def get_cases(project_name, database):
    '''
    (str, str) -> list
    
    Returns a list of dictionaries with case information
    
    Paramaters
    -----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT case_id, assay, donor_id, ext_id, species, miso FROM Samples WHERE project_id = ?", (project_name,)).fetchall()
    conn.close()
    
    data = [dict(i) for i in data]
         
    return data
    
    
    
def get_analysis_data(analysis_db, project_name, assay):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary of cases with analysis data corresponding to project and assay
    
    Parameters
    ----------
    - analysis_db (str): Path to the database storing analysis data
    - project_name (str): Name of the project of interest
    - assay (str): Name of the assay
    '''
    
    conn = connect_to_db(analysis_db)
    data = conn.execute("SELECT case_id, donor_id, template, valid, error FROM templates WHERE \
                        project_id = ? AND assay = ?;", (project_name,assay)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        case_id = i['case_id']
        template = json.loads(i['template'])
        valid = int(i['valid'])
        donor = i['donor_id']
        error = i['error']
        
        d = {'analysis': template, 'valid': valid, 'error': error, 'donor': donor}        
        
        assert case_id not in D
        D[case_id] = d
           
    return D    
    

def get_case_analysis_data(analysis_db, case_id, project_name, assay):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary of cases with analysis data corresponding to project and assay
    
    Parameters
    ----------
    - analysis_db (str): Path to the database storing analysis data
    - project_name (str): Name of the project of interest
    - assay (str): Name of the assay
    '''
    
    conn = connect_to_db(analysis_db)
    data = conn.execute("SELECT case_id, donor_id, template, valid, error FROM templates WHERE \
                        case_id = ? AND project_id = ? AND assay = ?;", (case_id, project_name, assay)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        template = json.loads(i['template'])
        valid = int(i['valid'])
        donor = i['donor_id']
        error = i['error']
        
        d = {'analysis': template, 'valid': valid, 'error': error, 'donor': donor}        
        
        assert case_id not in D
        D[case_id] = d
           
    return D    















def get_analysis_samples(analysis_data):
    '''
    
    
    
    '''
    
    
    D = {}
    
    for case_id in analysis_data:
        samples = []
        if 'analysis' in analysis_data[case_id] and analysis_data[case_id]['analysis']:
            for pipeline in analysis_data[case_id]['analysis']:
                if analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                    for workflow in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                        for d in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis'][workflow]:
                            if d['samples']:
                                samples.extend(d['samples'].split(','))
        samples = list(set(samples))            
        D[case_id] = samples            
                    
    return D                    
                    

















def get_analysis_workflows(analysis_data):
    '''
    
    
    
    '''
    
    
    D = {}
    
    for case_id in analysis_data:
        workflows, workflow_runs = [], []
        if 'analysis' in analysis_data[case_id] and analysis_data[case_id]['analysis']:
            for pipeline in analysis_data[case_id]['analysis']:
                if analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                    for workflow in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                        if analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis'][workflow]:
                            workflows.append(workflow)                 
                            for d in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis'][workflow]:
                                if d['wfrunid']:
                                    workflow_runs.append(d['wfrunid'])
        workflows = list(set(workflows))            
        workflow_runs = list(set(workflow_runs)) 
        D[case_id] = {'workflows': workflows,
                      'workflow_runs': workflow_runs}            
                    
    return D                    


def map_analysis_workflows(analysis_data, case_id):
    '''
    
    
    
    '''
    
    
    D = {}
    
    if 'analysis' in analysis_data[case_id] and analysis_data[case_id]['analysis']:
        for pipeline in analysis_data[case_id]['analysis']:
            if analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                for workflow in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                    for d in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis'][workflow]:
                        if d['wfrunid']:
                            D[d['wfrunid']] = workflow
                    
    return D                    




















def error_formatting(error):
    '''
    
    
    
    
    '''
    
    if 'missing workflows' in error.lower():
        message, workflows = error.split(':')
        workflows = workflows.strip().split(',')
    elif 'incomplete data' in error.lower():
        message, workflows = error.split(':')
        workflows = workflows.rstrip().replace('Workflows are missing ', '').split(',')
    else:
        message, workflows = error.split(':')
        workflows = ''
    
    message = message.strip().replace('[', '').replace(']', '').lower()
    
    err = {'message': message, 'workflows': workflows}
    
    return err 


def get_workflows_analysis_date(case_id, project_name, database):
    '''
    (str, str, str) -> dict
    
    Returns the creation date of any file for each workflow id for the case in project
           
    Parameters
    ----------
    - case_id (str): Case identifier
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
        
    # connect to db
    conn = connect_to_db(database)
    # extract project info
    data = conn.execute("SELECT DISTINCT creation_date, wfrun_id FROM Files WHERE case_id = ? AND  project_id= ?;", (case_id, project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        D[i['wfrun_id']] = i['creation_date']
        
    return D



def most_recent_analysis_workflow(analysis_data, case_id, creation_dates):
    '''
    (dict, dict) -> str
    
    Returns the most recent workflow creation in the analysis data of the case in project
           
    Parameters
    ----------
    - case_data (list): Dictionary with template information for each case
    - creation_dates (dict): Dictionary with creation dates of each workflow
    '''
        
    # get the workflow ids of all workflows
    workflows =  get_analysis_workflows(analysis_data)
    workflow_runs = workflows[case_id]['workflow_runs']
    
    L = sorted([creation_dates[wfrunid] for wfrunid in workflow_runs])
    
    try:
        most_recent = time.strftime('%Y-%m-%d', time.localtime(int(L[-1])))
    except:
        most_recent = 'NA'
        
    return most_recent
        
    
    
def map_workflows_to_fileids(case_id, project_name, database, wfrunid = None):
    '''
    (str, str, str, str | None) -> dict
    
    Returns a dictionary mapping each workflow run id of a case and project 
    to its output file swids
    
    Parameters
    ----------
    - case_id (str): Case identifier
    - project_name (str): Name of the project of interest
    - database (str): Path to the waterzooi database
    - wfrunid (str | None): Workflow run id
    '''    
        
    
    # connect to db
    conn = connect_to_db(database)
    if wfrunid:
        data = conn.execute("SELECT DISTINCT file_swid, wfrun_id FROM Files WHERE wfrun_id = ? AND case_id = ? AND  project_id= ?;", (wfrunid, case_id, project_name,)).fetchall()
    else:
        data = conn.execute("SELECT DISTINCT file_swid, wfrun_id FROM Files WHERE case_id = ? AND  project_id= ?;", (case_id, project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        if i['wfrun_id'] in D:
            D[i['wfrun_id']].append(i['file_swid'])
        else:
            D[i['wfrun_id']] = [i['file_swid']]
        
    return D

    
def organize_data(analysis_data, case_id):
    '''
    (dict, str) -> list

    '''
    
    # store the data as a list of dictionary for easier sorting and display
    data = []
    
    # map the workflow run ids to the workflow names
    analysis_workflows = {}
    
    # get the paren-children workflow relationships
    parents = {}
    
    
    
    if 'analysis' in analysis_data[case_id] and analysis_data[case_id]['analysis']:
        for pipeline in analysis_data[case_id]['analysis']:
            if analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                for workflow in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                    for d in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis'][workflow]:
                        if d['wfrunid']:
                            data.append(d)
                            analysis_workflows[d['wfrunid']] = workflow
                            if d['parents']:
                                for k in d['parents']:
                                    if k in parents:
                                        parents[k].append(d['wfrunid'])
                                    else:
                                        parents[k] = [d['wfrunid']]

    return data, analysis_workflows, parents
    
    
    
    
    
    
# def get_workflow_release_status(database, case_id):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with the release status of each workflow id of a given case
#     The release status is derived from the file qc status in Nabu of the workflow output files
    
#     Parameters
#     ----------
#     - database (str): Path to the waterzooi sqlite database
#     - case_id (str): Case identifier
#     '''

#     # get the file qc status for each output file of every workflows
#     workflow_qc = get_workflow_file_qc(database, case_id)
    
#     D = {}    

#     for workflow_id in workflow_qc:
#         if all(map(lambda x: x.isdigit(), workflow_qc[workflow_id])):
#             if all(map(lambda x: int(x), workflow_qc[workflow_id])):
#                 D[workflow_id] = True
#             elif any(map(lambda x: int(x), workflow_qc[workflow_id])):
#                 D[workflow_id] = True
#             elif all(map(lambda x: int(x), workflow_qc[workflow_id])) == False:
#                 D[workflow_id] = False
#         elif '1' in workflow_qc[workflow_id]:
#             D[workflow_id] = True
#         elif len(list(set(workflow_qc[workflow_id]))) == 1:
#             D[workflow_id] = '?'
        
         
#     return D        
    


def get_case_assay(database, project_id, case_id):
    '''
    
    
    
    '''
    
    # connect to db
    conn = connect_to_db(database)
    # extract project info
    data = conn.execute("SELECT DISTINCT assay FROM Samples WHERE case_id = ? AND  project_id= ?;", (case_id, project_id,)).fetchall()
    conn.close()
    
    assert len(data) == 1
    
    assay = data[0]['assay']
    
          
    return assay


def get_case_parent_to_children_workflows(database, case):
    '''
    (dict, str) -> dict
    
    Returns a dictionary of parent to children workflows for a single case
    
    Parameters
    ----------
    - database (str): Path to the database
    - case (str): Name of case of interest
    '''

    conn = connect_to_db(database)
    data = conn.execute("SELECT parents_id, children_id FROM Parents WHERE case_id = ?;", (case,)).fetchall()     
    conn.close()
    
    parent_to_children = {}
    for i in data:
        parent = i['parents_id']
        child = i['children_id']
        if parent in parent_to_children:
            parent_to_children[parent].append(child)
        else:
            parent_to_children[parent] = [child]
    
    return parent_to_children

    

def get_case_children_to_parents_workflows(parents_to_children):
    '''
    (dict) -> dict
    
    Returns a dictionary of child to parent workflows for a single case
    
    Parameters
    - parents_to_children (dict): Dictionary of parents to children workflow
                                  relationships for a single case
    '''
    
    child_to_parents = {}
    
    for parent in parents_to_children:
        for child in parents_to_children[parent]:
            if child in child_to_parents:
                child_to_parents[child].append(parent)
            else:
                child_to_parents[child] = [parent]
    
    return child_to_parents
    


def get_workflow_output_files(database, wfrun_id):
    '''
    (str, str) -> dict, dict
    
    Returns a dictionary with the output files of workflow with wfrun_id grouped by sample 
    and a dictionary with file paths mapped to file swids
    
    Parameters
    ----------
    - database (str): Path to the database
    - wfrun_id (str): Workflow run identifier
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Files.file, Files.file_swid, Libraries.sample_id FROM Files JOIN \
                        Workflow_Inputs JOIN Libraries WHERE Workflow_Inputs.wfrun_id = Files.wfrun_id \
                        AND Files.limskey = Workflow_Inputs.limskey AND Files.limskey = Libraries.lims_id \
                        AND Libraries.lims_id = Workflow_Inputs.limskey AND Files.wfrun_id = ?", (wfrun_id,)).fetchall()
    conn.close()   
    
    D = {}
    F = {}
    
    for i in data:
        sample = i['sample_id']
        file = i['file']
        fileswid = i['file_swid']
        if file in D:
            D[file].append(sample)
        else:
            D[file] = [sample]
        D[file] = sorted(list(set(D[file])))
        
        F[file] = fileswid
         
        
    # group samples sharing the same files
    S = {}
    for file in D:
        sample = ';'.join(D[file])
        if sample in S:
            S[sample].append(file)
        else:
            S[sample] = [file]
       
    return S, F

    
    
    # for i in data:
    #     sample = i['sample_id']
    #     file = i['file']
    #     if file in D:
    #         D[file].append(sample)
    #     else:
    #         D[file] = [sample]
    #     D[file] = sorted(list(set(D[file])))
            
    # # group samples sharing the same files
    # S = {}
    # for file in D:
    #     sample = ';'.join(D[file])
    #     if sample in S:
    #         S[sample].append(file)
    #     else:
    #         S[sample] = [file]
       
    # return S


def get_case_workflow_info(database, case):
    '''
    (str, str) -> dict
    
    Returns a dictionary of workflow name and workflow version for all workflows of a single case
        
    Parameters
    ----------
    - database (str): Path to the database
    - case (str): Case of interest
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT wfrun_id, wf, wfv FROM Workflows WHERE case_id = ?;", (case,)).fetchall()     
    conn.close()
    
    D = {}
    for i in data:
        workflow_id = i['wfrun_id']
        workflow_name = i['wf']
        version = i['wfv']
        D[workflow_id] = [workflow_name, version] 
    
    return  D



def map_limskeys_to_workflow(database, wfrun_id):
    '''
    (str, str) -> list

    Returns a list of limskeys matching workflow with identifier wfrun_id 

    Parameters
    ----------
    - database (str): Path to the waterzooi
    - wfrun_id (str): Workflow unique identifier
    '''

    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Workflow_Inputs.limskey FROM Workflow_Inputs \
                        WHERE Workflow_Inputs.wfrun_id = ?;", (wfrun_id,)).fetchall()
    conn.close()
    
    limskeys = [i['limskey'] for i in data]
    
    return limskeys




def get_input_sequences(database, case_id, wfrun_id):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with input sequences  of worflow with identifier wfrun_id
    
    Parameters
    ----------
    - database (str): Path to the waterzooi
    - case_id (str): Case identifier
    - wfrun_id (str): Workflow unique identifier
    '''
    
    # get the limskeys matching the workflow
    limskeys = map_limskeys_to_workflow(database, wfrun_id)

    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Files.file_swid, Files.file, Files.limskey, Libraries.library, \
                        Libraries.sample_id FROM Files JOIN Libraries JOIN Workflows \
                        WHERE Files.wfrun_id = Workflows.wfrun_id AND Files.limskey = Libraries.lims_id \
                        AND LOWER(Workflows.wf) IN ('casava', 'bcl2fastq', 'fileimportforanalysis', \
                        'fileimport', 'import_fastq') AND Files.case_id = ?;", (case_id,)).fetchall()
    conn.close()
    
    D = {}
    
    for i in data:
        sample = i['sample_id']
        library = i['library']
        limskey = i['limskey']
        file_swid = i['file_swid']
        file = i['file']
        
        # check that limskey match the limskeys of workflow wfrun_id
        if limskey in limskeys:
            if sample not in D:
                D[sample] = [[sample, library, limskey, file_swid, file]]
            else:
                D[sample].append([sample, library, limskey, file_swid, file])
    
    # sort according to sample and sequences
    for sample in D:
        D[sample].sort(key=lambda x: (x[0], x[2], x[-1]))
    
    return D


def add_workflow_qc_status(sequences, fileqc):
    '''
    (dict, dict) -> list
    
    Returns a list of dictionary with information of pairs of fastqs
    adding the release status of each pair at the workflow level
        
    Parameters
    ----------
    - sequences (list): List of sequence information for pairs of fastqs
    - filwqc (dict): Dictionary with file qc (ie release status) for each file
    '''
    
    # add qc status at the workflow level for each pair of fastqs in sequences
    for i in sequences:
        i['status'] = merge_qc_status_workflow(i['file_swids'], fileqc)
    
    return sequences


def get_sequences_to_download(sequences, platform_names, platforms):
    '''
    (list, dict, list) -> list
    
    Returns a dictionary with sequencing information for specific platforms
    
    Parameters
    ----------
    - sequences (list): List of sequence information for pairs of fastqs
    - platform_names (dict): Dictionary with the generic name of the sequencing platforms
    - plarforms (list): List of user-selected sequencing platforms
    '''
        
    L = []
    for i in sequences:
        d = {'Case': i['case_id'],
             'Donor': i['donor'],
             'DonorID': i['sample'],
             'SampleID': i['group_id'],
             'Sample': i['sample_id'],
             'Description': i['group_description'],
             'Library': i['library'],
             'Library Type': i['library_type'],
             'Tissue Origin': i['tissue_origin'],
             'Tissue Type': i['tissue_type'],
             'File Prefix': i['prefix']}
       
        if i['status']['qcstatus']:
            d['Released'] = 'YES'
            if i['status']['ticket']:
                d['Tickets'] = ';'.join(sorted(i['status']['ticket']))
        else:
            d['Released'] = 'NO'
            d['ticket'] = 'NA'
             

        # download all information if platforms are not selected
        if platforms:
            # check that platform is selected
            if platform_names[i['platform']] in platforms:
                L.append(d)    
        else:
            L.append(d)
                
    return L
    



def get_data_release_approval_signoff(signoffs, case_id):
    '''
    (dict, str) -> dict
    
    Returns a dictionary with True if Data release is complete for case_id and False otherwise
    
    Parameters
    ----------
    - signoffs (dict): Dictionary with case signoff extracted from the nabu cache
    - case_id (str): Case identifier
    '''
         
    if case_id in signoffs and 'release_approval' in signoffs[case_id] and \
      'Data Release' in signoffs[case_id]['release_approval'] and \
      signoffs[case_id]['release_approval']['Data Release']['qcpassed']:
        data_release = True
    else:
        data_release = False         
            
    return {case_id: data_release}    


def get_data_release_signoff(signoffs, project_deliverables, case_id, data):
    '''
    
    
    '''
    
    # collect data release deliverables and deliverable release signoffs
    if data == 'all':
        data_release_deliverables = [i for i in project_deliverables if 'pipeline' in i.lower()
                                     or 'fastq' in i.lower() or 'cbioportal' in i.lower()]
    elif data == 'cbioportal':
        data_release_deliverables = [i for i in project_deliverables if 'cbioportal' in i.lower()]
        
    elif data == 'pipeline':
        data_release_deliverables = [i for i in project_deliverables if 'pipeline' in i.lower()
                                     or 'fastq' in i.lower()]
        
    D = {}
    if case_id in signoffs:
        if 'release' in signoffs[case_id]:
            for i in signoffs[case_id]['release']:
                for j in signoffs[case_id]['release'][i]:
                    if j in data_release_deliverables:
                        assert j not in D
                        D[j] = signoffs[case_id]['release'][i][j]['qcpassed']
    # add deliverables not in release signoffs
    for i in data_release_deliverables:
        if i not in D:
            D[i] = False
        
    return {case_id: all(D.values())}




def get_output_files(database, project_id, case_id):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with the matching file paths and file swids
    
    Parameters
    ----------
    - database (str): Path to the database
    - wfrun_id (str): Workflow run identifier
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT file, file_swid FROM Files WHERE project_id = ? and case_id = ?;", (project_id, case_id)).fetchall()
    conn.close()   
    
    D = {}
    
    for i in data:
        file = i['file']
        fileswid = i['file_swid']
        assert fileswid not in D
        D[fileswid] = file
        
    return D




def get_workflow_outputs(database, project_name, case_id = None):
    '''
    (str, str, str | None) -> dict
    
    Returns a dictionary mapping each workflow run id to its output files for
    all cases in a project or for a single case 
    
    Parameters
    ----------
    - database (str): Path to the waterzooi database
    - project_name (str): Name of the project of interest
    - case_id (str | None): Case identifier
    '''    
        
    
    # connect to db
    conn = connect_to_db(database)
    if case_id:
        data = conn.execute("SELECT DISTINCT case_id, file, wfrun_id FROM Files WHERE case_id = ? AND  project_id= ?;", (case_id, project_name,)).fetchall()
    else:
        data = conn.execute("SELECT DISTINCT case_id, file, wfrun_id FROM Files WHERE project_id= ?;", (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        if i['case_id'] not in D:
            D[i['case_id']] = {}
        if i['wfrun_id'] in D[i['case_id']]:
            D[i['case_id']][i['wfrun_id']].append(i['file'])
        else:
            D[i['case_id']][i['wfrun_id']] = [i['file']]
        
    return D



def get_files_to_release(files, files_extensions):
    '''
    (list, list | None) -> list
    
    
    
    '''
    
    L = []
    
    for file in files:
        if files_extensions:
            for file_type in files_extensions:
                if file_type in file:
                    L.append(file)
    L = list(set(L))
    
    return L    
    
    




def prepare_analysis_json(analysis_data, workflow_outputs, workflow_deliverables = None):
    '''
    (dict, dict, dict | None) -> dict
    
    Returns a dictionary matching all the files to each workflow run id of the assay workflows
    Precondition: Data has passed validation and analysis_data contains all the required
    pipeline data
    
    Parameters
    ----------
    - analysis_data (dict): Dictionary with analysis data for a given assay 
    - workflow_outputs (dict): Dictionary with all the files for each workfflow run id
    - workflow_deliverables (dict | None): Dictionary with file outputs for workflows included in the pipeline deliverables
    '''
        
    D = {}
        
    for case_id in analysis_data:
        for pipeline in analysis_data[case_id]['analysis']:
            for workflow in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                for d in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis'][workflow]:
                    wfrunid = d['wfrunid']
                    # get the file paths    
                    files = workflow_outputs[case_id][wfrunid]
                    if workflow_deliverables:
                        # check that workflow is included in pipeline deliverables
                        if workflow in workflow_deliverables:
                            # get only the files that should included in the release
                            outputs = get_files_to_release(files, workflow_deliverables[workflow])
                            if outputs:
                                if case_id not in D:
                                    D[case_id] = {}
                                if workflow not in D[case_id]:
                                    D[case_id][workflow] = {}
                                D[case_id][workflow][wfrunid] = outputs
                    else:
                        if case_id not in D:
                            D[case_id] = {}
                        if workflow not in D[case_id]:
                            D[case_id][workflow] = {}
                        D[case_id][workflow][wfrunid] = files
    
    
    return D
    
    
   
def prepare_cbioportal_json(analysis_data, workflow_outputs):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with required output files for the cbioportal importer
    Precondition: Data has passed validation and analysis_data contains all the required
    pipeline data
    
    Parameters
    ----------
    - analysis_data (dict): Dictionary with analysis data for a given assay 
    - workflow_outputs (dict): Dictionary with all the files for each workfflow run id
    '''
    
    D = {}  
    
    for case_id in analysis_data:
        donor = analysis_data[case_id]['donor']
        for pipeline in analysis_data[case_id]['analysis']:
            for workflow in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis']:
                if any(['varianteffectpredictor' in workflow.lower(), 'rsem' in workflow.lower(),
                        'mavis' in workflow.lower(), 'purple' in workflow.lower(),
                        'sequenza' in workflow.lower()]):
                    for d in analysis_data[case_id]['analysis'][pipeline]['pipeline_analysis'][workflow]:
                        tumor_sample, file = '', ''
                        wfrunid = d['wfrunid']
                        # get the file paths    
                        files = workflow_outputs[case_id][wfrunid]
                        # get the samples
                        samples = d['samples'].split(',')
                        # find the tumour sample
                        for i in samples:
                            if 'Ly' not in i:
                                tumor_sample = i
                                break
                        assert tumor_sample and 'Ly' not in tumor_sample
                                                                  
                        if donor not in D:
                            D[donor] = {}
                        if tumor_sample not in D[donor]:
                            D[donor][tumor_sample] = {}
                        
                        if 'varianteffectpredictor' in workflow.lower():
                            for file in files:
                                if 'mutect2.filtered.maf.gz' in file:
                                    D[donor][tumor_sample][workflow] = file    
                                    break
                        elif 'rsem' in workflow.lower():
                            for file in files:
                                if '.genes.results' in file:
                                    D[donor][tumor_sample][workflow] = file    
                                    break
                        elif 'mavis' in workflow.lower():
                            for file in files:
                                if '.mavis_summary.tab' in file:
                                    D[donor][tumor_sample][workflow] = file    
                                    break
                        elif 'purple' in workflow.lower():
                            D[donor][tumor_sample][workflow] = {}
                            for file in files:
                                if '.purple.cnv.somatic.tsv' in file:
                                    D[donor][tumor_sample][workflow]['cnv'] = file
                                elif '.purple.purity.tsv' in file:
                                    D[donor][tumor_sample][workflow]['purity'] = file
                        elif 'sequenza' in workflow.lower():
                            for file in files:
                                if 'results.sequenza.zip' in file:
                                    D[donor][tumor_sample][workflow] = file
                                    break

    # remove samples and donors without data                    
    for donor in D:
        to_remove = [sample for sample in D[donor] if len(D[donor][sample]) == 0]
        for i in to_remove:
            del D[donor][i]
    to_remove = [donor for donor in D if len(D[donor]) == 0]
    for i in to_remove:
        del D[i]
                
    return D
    
    
def count_cases(analysis_data, data_release_approval, data_release, pipeline_signoff, cbio_signoff):
    
    '''
    
    
    
    
    
    '''
    
    complete = len([case_id for case_id in analysis_data if analysis_data[case_id]['valid']])
    incomplete = len(analysis_data) - complete
    
    # count complete cases with release approval and data release signed off
    complete_signedoff = len([case_id for case_id in analysis_data if analysis_data[case_id]['valid'] 
                          and data_release_approval[case_id] and data_release[case_id]])
            
    complete_to_release = len([case_id for case_id in analysis_data if analysis_data[case_id]['valid']
                           and data_release_approval[case_id] and data_release[case_id] == False])
                              
    pipeline_to_release = len([case_id for case_id in analysis_data if analysis_data[case_id]['valid']
                           and data_release_approval[case_id] and pipeline_signoff[case_id] == False])
                          
    cbio_to_release = len([case_id for case_id in analysis_data if analysis_data[case_id]['valid']
                           and data_release_approval[case_id] and cbio_signoff[case_id] == False])
                          
    
    return complete, incomplete, complete_signedoff, complete_to_release, pipeline_to_release, cbio_to_release
    
    
    
def plot_graph(edges, workflow_names):
    '''
    (list, dict) -> plotly.graph_objs._figure.Figure
       
    Returns  plotly figure of a graph showing the relationships among workflows
    
    Parameters
    ----------
    - edges (list): List of connected pairs of workflow ids
    - workflow_names (dict): Dictionary mapping workflow identifiers to their name
    '''
    
    # create the graph of workflow relationships
    G = nx.Graph()
    G.add_edges_from(edges)
    
    # add a graph layout and get positions
    pos = nx.spring_layout(G)
    
    # get edge positions
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    # get node positions
    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    # plot the edges
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='#888'),
        hoverinfo='none',
        mode='lines')
    
    # plot the nodes
    node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    marker=dict(
        showscale=True,
        # colorscale options
        #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
        #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
        #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
        colorscale='Viridis',
        reversescale=True,
        color=[],
        size=10,
        colorbar=dict(
            thickness=15,
            title=dict(
              text='Node Connections',
              side='right'
            ),
            xanchor='left',
        ),
        line_width=1.5))
    
    
    # color the nodes based on the number of connection
    node_adjacencies = [len(list(G.neighbors(node))) for node in G.nodes()]
    node_trace.marker.color = node_adjacencies
    
    # to change the size of the marker based on the number of connection
    #node_trace.marker.size = node_adjacencies
    
    # label the nodes with the workflow names
    node_text = [str(node) for node in G.nodes()]
    node_text = [workflow_names[i] for i in node_text]
    node_trace.text = node_text
    
    # # generate figure
    # fig = go.Figure(data=[edge_trace, node_trace],
    #              layout=go.Layout(
    #                 title='Workflow connections',
    #                 showlegend=False,
    #                 hovermode='closest',
    #                 margin=dict(b=20,l=5,r=5,t=40),
    #                 xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    #                 yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
    #                 )
    
    
    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title='Workflow connections',
                    showlegend=False,
                    hovermode='closest',
                    #margin=dict(b=20,l=5,r=5,t=40),
                    margin=dict(b=0,l=0,r=0,t=0),
                    
                    
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    height=350,  
                    width=1200,    
                    autosize=True
                    )
                    )
    
    
    
    # fig.update_layout(
    #     title='Workflow connections',
    # autosize=True,
    # margin=dict(l=0, r=0, t=0, b=0), # Strip padding for small spaces
    # height=300,  
    # width=1200,       
    # )
    
    
    
    return fig
    
    

def plot_small_graph(edges, workflow_names, wfrunid, parents, children):
    '''
    (list, dict) -> plotly.graph_objs._figure.Figure
       
    Returns  plotly figure of a graph showing the relationships among workflows
    
    Parameters
    ----------
    - edges (list): List of connected pairs of workflow ids
    - workflow_names (dict): Dictionary mapping workflow identifiers to their name
    '''
    
    # create the graph of workflow relationships
    G = nx.Graph()
    G.add_edges_from(edges)
    
    # add a graph layout and get positions
    pos = nx.spring_layout(G)
    
    # get edge positions
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    # get node positions
    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    # plot the edges
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.8, color='#888'),
        hoverinfo='none',
        mode='lines')
    
    # plot the nodes
    node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    marker=dict(
        showscale=False,
        # colorscale options
        #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
        #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
        #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
        colorscale='Viridis',
        reversescale=True,
        color=[],
        size=8,
        colorbar=dict(
            thickness=10,
            title=dict(
              text='Node Connections',
              side='right'
            ),
            xanchor='left',
        ),
        line_width=1))
    
    
    # color the nodes based on the number of connection
    # node_adjacencies = [len(list(G.neighbors(node))) for node in G.nodes()]
    # node_trace.marker.color = node_adjacencies
    
    node_colors = []
    for node in G.nodes():
        if str(node) == wfrunid:
            node_colors.append('#ff6666')
        elif str(node) in parents:
            node_colors.append('#0073e6')
        elif str(node) in children:
            node_colors.append('#2eb82e')
    node_trace.marker.color = node_colors
    
    
    
    # to change the size of the marker based on the number of connection
    #node_trace.marker.size = node_adjacencies
    
    # label the nodes with the workflow names
    node_text = [str(node) for node in G.nodes()]
    node_text = [workflow_names[i] for i in node_text]
    node_trace.text = node_text
    
    # generate figure
    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title=None,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    
    fig.update_layout(
    autosize=True,
    margin=dict(l=0, r=0, t=0, b=0), # Strip padding for small spaces
    height=150,  
    width=450,       
    )
    
#     fig.update_layout(
#     height=150,          # Force the height to match your HTML <div> element
#     autosize=True,       # Allows it to dynamically fill the 100% width of the <td>
#     margin=dict(l=10, r=10, t=10, b=10), # Minimize padding to prevent cropping
# )
    
    
    
    
    
    return fig
    



def convert_epoch_time(epoch):
    '''
    (str) -> str
    
    Returns epoch time in readable format
    
    Parameters
    ----------
    - epoch (str)
    '''
    
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(epoch)))




def get_last_sequencing(project_name, database):
    '''
    (str, str) -> str
    
    Returns the date of the last sequencing for the project of interest
    
    Paramaters
    ----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    sequencing = conn.execute("SELECT DISTINCT Files.creation_date FROM Files JOIN Workflows \
                              WHERE Files.project_id = '{0}' AND Workflows.project_id = '{0}' \
                              AND Workflows.wfrun_id = Files.wfrun_id AND LOWER(Workflows.wf) in \
                              ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq');".format(project_name)).fetchall()
    conn.close()
    
    # get the most recent creation date of fastq generating workflows
    if sequencing:
        seq_dates = sorted([i['creation_date'] for i in sequencing])
        most_recent = seq_dates[-1]
    else:
        most_recent = 'NA'
        
    try:
        most_recent = convert_epoch_time(most_recent)    
        return most_recent
    except:
        return most_recent



def rename_case_id(case_id):
    '''
    (str) -> str
    
    Returns the case id replacing the en dash with an hyphen
    
    Parameters
    ----------
    - case_id (str): Case identifier
    '''
    
    #convert en-dash to hyphen in file name
    if "\u2013" in case_id:
        case_name = case_id.replace("\u2013", '-')
    else:
        case_name = case_id
    
    return case_name

    

def create_graph_edges(workflow_ids, parent_to_children):
    '''
    (list, dict) -> list
    
    Returns a list of tuples, each with 2 workflow identifiers when there is a connection
    (ie parent to child) between these 2 workflows

    Parameters
    ----------
    - workflow_ids (list): List of all the workflow ids of a template of a case
    - parent_to_children (dict): Dictionary with parent to children workflow relationships 
    '''

    edges = []
        
    for i in workflow_ids:
        for j in workflow_ids:
            if i != j and (i in parent_to_children or j in parent_to_children):
                if i in parent_to_children:
                    if j in parent_to_children[i]:
                        edges.append((i, j))
                else:
                    if i in parent_to_children[j]:
                        edges.append((j, i))
    return edges


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


