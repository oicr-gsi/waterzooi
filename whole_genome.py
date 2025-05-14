# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:42:36 2023

@author: rjovelin
"""

import os
import itertools
import json
import time
from utilities import connect_to_db, convert_epoch_time, remove_non_analysis_workflows,\
    get_children_workflows, get_workflow_names, get_donors



    




def get_parent_workflows(project_name, database):
    '''
    (str, str) -> dict
    
    Returns a dictionary with workflow name, list of workflow_ids that are parent
    to all each workflow (i.e immediate upstream workflow) for a given project
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Workflows.wf, Parents.parents_id, Parents.children_id \
                        FROM Parents JOIN Workflows WHERE Parents.project_id = ? \
                        AND Workflows.project_id = ? AND Workflows.wfrun_id = Parents.parents_id;", (project_name, project_name)).fetchall()
    data= list(set(data))
    conn.close()
    
    D = {}
    for i in data:
        if i['children_id'] not in D:
            D[i['children_id']] = {}
        if i['wf'] not in D[i['children_id']]:
            D[i['children_id']][i['wf']] = []
        D[i['children_id']][i['wf']].append(i['parents_id'])
    return D
      




def get_workflows_analysis_date(project_name, database):
    '''
    (str, str) -> dict
    
    Returns the creation date of any file for each workflow id for the project of interest
       
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
        
    # connect to db
    conn = connect_to_db(database)
    # extract project info
    data = conn.execute("SELECT DISTINCT creation_date, wfrun_id FROM Files WHERE project_id= ?;", (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        D[i['wfrun_id']] = i['creation_date']
        
    return D



def get_workflow_file_count(project_name, database, workflow_table='Workflows'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with the number of files for each workflow in project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - workflow_table (str): Name of the table containing the workflow information in database
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT {0}.file_count, {0}.wfrun_id FROM {0} WHERE {0}.project_id = ?;".format(workflow_table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        counts[i['wfrun_id']] = i['file_count']
    
    return counts


def get_workflow_limskeys(project_name, database, workflow_input_table='Workflow_Inputs'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with list of limskeys for each workflow id in project
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - workflow_input_table (str): Name of the table with worklow input information in the database
    '''
        
    conn = connect_to_db(database)
    query = "SELECT {0}.limskey, {0}.wfrun_id FROM {0} WHERE {0}.project_id = ?;".format(workflow_input_table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['wfrun_id'] not in D:
            D[i['wfrun_id']] = []
        D[i['wfrun_id']].append(i['limskey'])
            
    return D


   
def get_amount_data(project_name, database, workflow_table='Workflows'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with the amount of data (ie, lane count) for each workflow in project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - workflow_table (str): Name of the table containing the workflow information in database
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT {0}.lane_count, {0}.wfrun_id FROM {0} WHERE {0}.project_id = ?;".format(workflow_table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        counts[i['wfrun_id']] = i['lane_count']
    
    return counts



def create_WG_block_json(database, project_name, case, blocks, block, anchor_workflow, workflow_names, selected_workflows, selection):
    '''
    (str, str, dict, str, str, dict, dict, str)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor parent workflow (bmpp or star)
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - project_name (str): Name of project of interest
    - case (str): Donor identifier 
    - blocks (dict): Dictionary with block information
    - block (str): Sample pair in blocks
    - anchor_workflow (str): bamMergePreprocessing parent workflow(s) or star_call_ready parent workflow
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - selection (str): Include files from all selected workflows or files from the standard deliverables
                       Values: standard or all
    '''
    
    
    libraries = map_limskey_to_library(project_name, database, table='Workflow_Inputs')
    sample_names = map_library_to_sample(project_name, database, table = 'Libraries')
    donors = map_library_to_case(project_name, database, table = 'Libraries')
    workflow_outputfiles = get_workflow_output(project_name, database, libraries, sample_names, donors, 'Files')
    
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
            
    # get the deliverables
    if selection == 'standard':
        deliverables = get_WGS_standard_deliverables()
    elif selection == 'all':
        deliverables = {}
    
    # organize the workflows by block and samples
    D = {}
    # re-organize the sample pair
    sample_id = '.'.join(list(map(lambda x: x.strip(), block.split('|'))))
    # get the workflow ids for that block
    for i in blocks[block]:
        if i['anchor_wf'] == anchor_workflow:
            D[sample_id] = map(lambda x: x.strip(), i['workflows'].split(';'))
    
    block_data = {}
    for sample in D:
        for workflow_id in D[sample]:
            # check if workflow is selected
            if selected_workflows[workflow_id]:
                # get workflow name and version
                workflow_name = workflow_names[workflow_id][0]
                workflow_version = workflow_names[workflow_id][1]
                # needed to sort outputs by sample pairs or by sample for call-ready workflows
                # even if all files are recorded
                #outputfiles = get_workflow_output(project_name, case, workflow_id, database, libraries, sample_names, 'Files')
                outputfiles = workflow_outputfiles[workflow_id]
                
                # check that only workflows in standard WGS deliverables are used
                if deliverables:
                    key = workflow_name.split('_')[0].lower()
                    if key in deliverables:
                        
                        for j in outputfiles:
                            # list all deliverable files
                            L = []
                            # gather all file paths for workflow and sample(s)
                            files = [i[0] for i in outputfiles[j]]
                            # map all file endings of deliverables with files
                            groups = list(itertools.product(files, deliverables[key]))
                            # determine which files are part of the deliverables
                            F = list(map(G, groups))
                            L = [groups[k][0] for k in range(len(F)) if F[k]]
                            
                            if L:
                                sample_id = j.replace(';', '.')
                                if case not in block_data:
                                    block_data[case] = {}
                                if sample_id not in block_data[case]:
                                    block_data[case][sample_id] = {}
                                if workflow_name not in block_data[case][sample_id]:
                                    block_data[case][sample_id][workflow_name] = []
                                
                                d = {'workflow_id': workflow_id,
                                     'workflow_version': workflow_version,
                                     'files': L}
                                if d not in block_data[case][sample_id][workflow_name]:
                                    block_data[case][sample_id][workflow_name].append(d)
                                    
                else:
                    for j in outputfiles:
                        sample_id = j.replace(';', '.')
                        d =  {'workflow_id': workflow_id, 
                              'workflow_version': workflow_version,
                              'files': [i[0] for i in outputfiles[j]]}
                        if case not in block_data:
                            block_data[case] = {}
                        if sample_id not in block_data[case]:
                            block_data[case][sample_id] = {}
                        if workflow_name not in block_data[case][sample_id]:
                            block_data[case][sample_id][workflow_name] = []
                        if d not in block_data[case][sample_id][workflow_name]:
                            block_data[case][sample_id][workflow_name].append(d)
                    
    return block_data                





def create_WGS_project_block_json(project_name, database, blocks, block_status, selected_workflows, workflow_names, deliverables=None):
    '''
    (str, str, dict, dict, dict, dict, None | dict)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor bmpp parent workflow
    
    Parameters
    ----------
    - project_name (None | str): None or name of project of interest
    - database (None | str): None or path to the sqlite database
    - blocks (dict): Dictionary with block information
    - block_status (dict): Dictionary with review status of each block
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
    - deliverables (None | dict): None or dictionary with file extensions of standard WGS deliverables
    '''
    
    libraries = map_limskey_to_library(project_name, database, table='Workflow_Inputs')
    sample_names = map_library_to_sample(project_name, database, table = 'Libraries')
    donors = map_library_to_case(project_name, database, table = 'Libraries')
    workflow_outputfiles = get_workflow_output(project_name, database, libraries, sample_names, donors, 'Files')
      
    
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
    
    D = {}
    for case in blocks:
        for samples in blocks[case]:
            # check the selection status of the block
            if block_status[case][samples] not in ['ready', 'review']:
                # block already reviewed and workflows selected
                anchor_wf = block_status[case][samples]
                for workflow in blocks[case][samples][anchor_wf]['workflows']:
                    workflow = os.path.basename(workflow)
                    
                    # get workflow name and version
                    workflow_name = workflow_names[workflow][0]
                    workflow_version = workflow_names[workflow][1]
                               
                    # check workflow status
                    if selected_workflows[workflow]:
                        # get workflow output files
                        # needed to sort outputs by sample pairs or by sample for call-ready workflows
                        # even if all files are recorded
                        #outputfiles = get_workflow_output(project_name, case, workflow, database, libraries, sample_names, 'Files')
                        outputfiles = workflow_outputfiles[workflow]                        
                        
                        # check that only workflows in standard WGS deliverables are used
                        if deliverables:
                            key = workflow_names[workflow][0].split('_')[0].lower()
                                                       
                            if key in deliverables:
                                for j in outputfiles:
                                    # list all deliverable files
                                    L = []
                                    # gather all file paths for workflow and sample(s)
                                    files = [i[0] for i in outputfiles[j]]
                                    # map all file endings of deliverables with files
                                    groups = list(itertools.product(files, deliverables[key]))
                                    # determine which files are part of the deliverables
                                    F = list(map(G, groups))
                                    L = [groups[k][0] for k in range(len(F)) if F[k]]
                                    
                                    if L:
                                        sample_id = j.replace(';', '.')
                                        if case not in D:
                                            D[case] = {}
                                        if sample_id not in D[case]:
                                            D[case][sample_id] = {}
                                        if workflow_name not in D[case][sample_id]:
                                            D[case][sample_id][workflow_name] = []
                                        
                                        d = {'workflow_id': workflow,
                                             'workflow_version': workflow_version,
                                             'files': L}    
                                        if d not in D[case][sample_id][workflow_name]: 
                                            D[case][sample_id][workflow_name].append(d)
                        
                        else:
                            for j in outputfiles:
                                sample_id = j.replace(';', '.')
                                d =  {'workflow_id': workflow,
                                      'workflow_version': workflow_version,
                                      'files': [i[0] for i in outputfiles[j]]}
                                if case not in D:
                                    D[case] = {}
                                if sample_id not in D[case]:
                                    D[case][sample_id] = {}
                                if workflow_name not in D[case][sample_id]:
                                    D[case][sample_id][workflow_name] = []
                                if d not in D[case][sample_id][workflow_name]:
                                    D[case][sample_id][workflow_name].append(d)
                                        
    
    return D



def get_call_ready_cases(project_name, platform, library_type, database):
    '''
    (str, str, str, str) -> dict

    Returns a dictionary with samples and libraries and bmpp and downstream workflow ids for each case in a project,
    restricting data to specified platform and library type

    Parameters
    ----------
    - project_name (str): Name of the project
    - platform (str): Name of sequencing platform.
                      Accepted values: novaseq, nextseq, hiseq, miseq
    - library_type (str): 2 letters-code indicating the type of library                   
    - database (str): Path to the sqlite database
    '''

    # get all the samples for project name 
    conn = connect_to_db(database)
    query = "SELECT Libraries.library, Libraries.case_id, Libraries.project_id, \
             Libraries.ext_id, Libraries.group_id, Libraries.library_type, \
             Libraries.tissue_type, Libraries.tissue_origin, \
             Workflows.wf, Workflows.wfrun_id, Workflow_Inputs.platform \
             from Workflow_Inputs JOIN Libraries JOIN Workflows \
             WHERE Libraries.project_id = ? AND Workflow_Inputs.project_id = ? \
             AND Workflows.project_id = ? AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id \
             AND Workflow_Inputs.library = Libraries.library AND Libraries.library_type = ?;"
    data = conn.execute(query, (project_name, project_name, project_name, library_type)).fetchall()
    conn.close()
    
    cases = {}
    for i in data:
        # select bmpp data sequenced on novaseq
        if platform in i['platform'].lower():
            if 'bammergepreprocessing' in i['wf'].lower():
                if i['case_id'] not in cases:
                    cases[i['case_id']] = {'project': i['project_id'], 'samples': set(), 'libraries': set(), 'bmpp': set()}
                cases[i['case_id']]['bmpp'].add(i['wfrun_id'])
                sample = '_'.join([i['case_id'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
                cases[i['case_id']]['samples'].add(sample)
                cases[i['case_id']]['libraries'].add(i['library'])
            
    # get parent-children workflow relationships
    parents = get_children_workflows(project_name, database)
    
    # find the bmpp downstream workflows
    for sample in cases:
        downstream = []
        for bmpp in cases[sample]['bmpp']:
            if bmpp in parents:
                # get the bmpp downstream workflows
                children = parents[bmpp]
                # removed any non-analysis workflow
                children = remove_non_analysis_workflows(children)
                # list all downtream workflows
                downstream.extend([i['children_id'] for i in children])
                # get the downstream workflows of downstream workflows
                # remove non-analysis workflows
                for workflow in downstream:
                    if workflow in parents:
                        L = remove_non_analysis_workflows(parents[workflow])
                        downstream.extend([i['children_id'] for i in L])
        cases[sample]['downstream'] = list(set(downstream)) 
    
    
    return cases




def get_call_ready_samples(project_name, bmpp_run_id, database):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with normal and tumour samples from project_name processed through bamMergePreprcessing
    workflow with bmpp_run_id 
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - bmpp_run_id (str): BamMergePreprocessing workflow run identifier
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    query = "SELECT Libraries.case_id, Libraries.group_id, Libraries.library, Libraries.tissue_type, \
                        Libraries.tissue_origin, Libraries.library_type \
                        FROM Libraries JOIN Workflow_Inputs WHERE Workflow_Inputs.library = Libraries.library \
                        AND Workflow_Inputs.wfrun_id = ? AND Libraries.project_id = ? \
                        AND Workflow_Inputs.project_id = ?"
    data = conn.execute(query, (os.path.basename(bmpp_run_id), project_name, project_name)).fetchall()
    conn.close()

    data = list(set(data))
    
    samples = {'normal': [], 'tumour': []}
    for i in data:
        if i['tissue_type'] == 'R':
            tissue = 'normal'
        else:
            tissue = 'tumour'
        sample = '_'.join([i['case_id'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
        if sample not in samples[tissue]:
            samples[tissue].append(sample)

    return samples



def map_samples_to_bmpp_runs(project_name, bmpp_ids, database):
    '''
    (str, list, str) -> dict
    
    Returns a dictionary with normal, tumor samples for each bmpp run id
      
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - bmpp_ids (list): List of BamMergePreprocessing workflow run identifiers for a single case
    - database (str): Path to the sqlite database
    '''

    D = {}
    for i in bmpp_ids:
        # initiate dictionary
        samples = get_call_ready_samples(project_name, i, database)
        D[i] = samples
    return D



def get_WGTS_blocks_info(project_name, case, database, table):
    '''
    (str, str, str, str) -> list 
    
    Returns a list of dictionaries containing WGS or WT block information for a given project and case
                
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - case (str): Donor id 
    - database (str): Path to the sqlite database
    - table (str): Table with block information: WGS_blocks, WT_blocks or EX_blocks
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT samples, anchor_wf, workflows, name, date, \
             complete, clean, network from {0} WHERE project_id = ? AND \
             case_id = ?;".format(table)
    data = conn.execute(query, (project_name, case)).fetchall() 
    conn.close()

    L = [dict(i) for i in data]

    D = {}
    # group by samples
    for i in L:
        samples = i['samples']
        if samples not in D:
            D[samples] = []
        # add call ready workflows
        call_ready = list(map(lambda x: x.strip(), i['anchor_wf'].split('.')))
        i['call_ready'] = call_ready    
        workflows = list(map(lambda x: x.strip(), i['workflows'].split(';')))
        # add caller workflows
        callers = set(workflows).difference(set(call_ready))
        i['callers'] = callers
        # map each sample to the
        bmpp_samples = map_samples_to_bmpp_runs(project_name, call_ready, database)
        i['pairs'] = bmpp_samples
        D[samples].append(i)
        # sort according to sub-block name
        D[samples].sort(key = lambda x: x['name'])

    return D    


def get_sequencing_platform(project_name, database, table = 'Workflow_Inputs'):
    '''
    (str, str, str) -> list 
    
    Returns a dictionary with the sequencing platform of input raw sequences
    for each workflow for project
                
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table with workflow input information
    '''
                
    conn = connect_to_db(database)
    query = "SELECT DISTINCT wfrun_id, platform FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall() 
    conn.close()

    D = {}
    for i in data:
        D[i['wfrun_id']] = i['platform']
    
    return D




# def get_selected_workflows(project_name, database, table = 'Workflows'):
#     '''
#     (str, str, str) -> dict 
    
#     Returns a dictionary with the selected status of each workflow for the given project
                  
#     Parameters
#     ----------
#     - project_name (str): Name of project of interest
#     - database (str): Path to the sqlite database
#     - table (str): Table with workflow information
#     '''
                
#     conn = connect_to_db(database)
#     query = "SELECT DISTINCT wfrun_id, selected FROM {0} WHERE project_id = ?;".format(table)
#     data = conn.execute(query, (project_name,)).fetchall() 
#     conn.close()

#     D = {}
#     for i in data:
#         D[i['wfrun_id']] = int(i['selected'])
    
#     return D

    
def get_case_workflows(case, database, table = 'WGS_blocks'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary of all workflows in each block and sample pair for a given case
    
    Parameters
    ----------
    - case (str): Donor identifier
    - database (str): Path to the sqlite database
    - table (str): Name of the table storing analysis blocks
    '''
        
    conn = connect_to_db(database)    
    query = "SELECT samples, anchor_wf, workflows FROM {0} WHERE case_id = ?;".format(table)
    data = conn.execute(query, (case,)).fetchall()
    conn.close()
    
    D = {}
        
    for i in data:
        samples = i['samples']
        block = i['anchor_wf']
        workflows = i['workflows'].split(';')
        if samples not in D:
            D[samples] = {}
        if block not in D[samples]:
            D[samples][block] = []
        D[samples][block].extend(workflows)
        D[samples][block] = list(set(workflows))
    
    return D


def update_wf_selection(workflows, selected_workflows, selection_status, database, table='Workflows'):
    '''
    (list, list, dict, str, str)
    
    Update the selection status of workflows 
    
    Parameters
    ----------
    - workflows (list): List of workflows across templates from a case
    - selected_workflows (list): List of selected workflows from the application form for a given case
    - selection_status (dict): Selection status of all workflows for a given project
    - database (str): Path to the sqlite database
    - table (str): Table storing workflows information
    '''
    
    # update selected status
    conn = connect_to_db(database)
    for i in workflows:
        if i in selected_workflows:
            status = 1
        else:
            status = 0
        
        # update only if status has changed
        if selection_status[i] != status:
            query = 'UPDATE Workflows SET selected = ? WHERE wfrun_id = ?;'
            conn.execute(query, (status, i))
            conn.commit()
    conn.close()



def get_wgs_blocks(project, database, table = 'WGS_blocks'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with analysis block names for each case in project
    
    Parameters
    ----------
    - project (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table with analysis blocks
    '''

    conn = connect_to_db(database)    
    query = "SELECT DISTINCT * FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project,)).fetchall()
    conn.close()
    
    D = {}
    
    for i in data:
        case, samples, anchor = i['case_id'], i['samples'], i['anchor_wf']
        workflows = i['workflows'].split(';')
        clean, complete = int(i['clean']), int(i['complete'])
        if case not in D:
            D[case] = {}
        if samples not in D[case]:
            D[case][samples] = {}
        assert anchor not in D[case][samples]
        D[case][samples][anchor] = {'workflows': workflows, 'clean': clean,
                                    'complete': complete}
            
    return D
    
             
def get_block_counts(analysis_blocks):
    '''
    (dict) -> dict
    
    Returns a dictionary with block counts for each case and sample pairs in given project
    
    Parameters
    ----------
    - analysis_blocks (dict): Dictionary with analysis blocks for a given project
    '''
    
    D = {}
    
    for i in analysis_blocks:
        for j in analysis_blocks[i]:
            if i not in D:
                D[i] = {}
            assert j not in D[i]
            D[i][j] = len(analysis_blocks[i][j])
    return D  



def review_wgs_blocks(blocks, selected_workflows):
    '''
    (dict, dict) -> dict 
    
    Returns a dictionary with status for analysis blocks for each case in project
                  
    Parameters
    ----------
    - blocks (dict): 
    - selected_workflows (dict): 
    '''
    
    D = {}
        
    for case in blocks:
        if case not in D:
            D[case] = {}
        for samples in blocks[case]:
            for anchor in blocks[case][samples]:
                # do not include call-ready workflows to determine selection/review status
                # these may be shared across multiple blocks
                L = [selected_workflows[os.path.basename(i)] for i in blocks[case][samples][anchor]['workflows'] if i not in anchor] 
                
                if any(L):
                    D[case][samples] = anchor
                    break
                else:
                    if blocks[case][samples][anchor]['clean'] and \
                      blocks[case][samples][anchor]['complete']:
                          D[case][samples] = 'ready'
                          break
                    else:
                        D[case][samples] = 'review'
    return D



def map_limskey_to_library(project_name, database, table='Workflow_Inputs'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary mapping limskey ids to library ids for each workflow in project
    
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - database (str): Path to the sqlite database
    - table (str): Table storing the workflow input information
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT library, limskey, wfrun_id FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        workflow = i['wfrun_id']
        if workflow not in D:
            D[workflow] = {}
        assert i['limskey'] not in D[workflow]    
        D[workflow][i['limskey']] = i['library']  
    
    return D




def map_library_to_sample(project_name, database, table = 'Libraries'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary mapping sample ids to library ids
        
    Parameters
    ----------
    - project_name (str): Name of the project of interest
    - database (str): Path to the sqlite database
    - table (str): Table storing the libraries information
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT library, case_id, tissue_type, tissue_origin, \
             library_type, group_id FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    
    D = {}
    for i in data:
        donor = i['case_id']
        library = i['library']
        sample = [i['case_id'], i['tissue_type'], i['tissue_origin'],
                           i['library_type'], i['group_id']]
        if not i['group_id']:
            sample = sample[:-1]
        sample = '_'.join(sample)    
        
        if donor not in D:
            D[donor] = {}
        if library in D[donor]:
            assert D[donor][library] == sample
        else:
            D[donor][library] = sample
        
    return D



def map_library_to_case(project_name, database, table = 'Libraries'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary mapping each library to its donor identifier
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table in database storing library information.
                   Default is Libraries
    '''
    
    # get the workflow output files sorted by sample
    conn = connect_to_db(database)
    query = "SELECT DISTINCT library, case_id FROM {0} WHERE project_id = ?".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        assert i['library'] not in D
        D[i['library']] = i['case_id']
          
    return D



def get_workflow_output(project_name, database, libraries, samples, donors, table = 'Files'):
    '''
    (str, str, dict, dict, dict, str) -> dict
    
    Returns a dictionary with workflow output files sorted by sample
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - libraries (dict): Dictionary mapping libraries to limskeys
    - samples (dict): Dictionary mapping libraries to samples
    - donors (dict): Dictionary mapping libraries to donors
    - table (str): Table in database storing File information
    '''

    # get the workflow output files sorted by sample
    conn = connect_to_db(database)
    query = "SELECT DISTINCT file, limskey, file_swid, wfrun_id FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall()
    conn.close()
    
    D = {}

    for i in data:
        file = i['file']
        limskeys = i['limskey'].split(';')
        fileswid = i['file_swid']
        workflow_id = i['wfrun_id']
        libs = list(set([libraries[workflow_id][j] for j in limskeys]))
        
        #sample_names = ';'.join(sorted(list(set([samples[case][j] for j in libs]))))
        
        sample_names = ';'.join(sorted(list(set([samples[donors[j]][j] for j in libs]))))
        
        if workflow_id not in D:
            D[workflow_id] = {}
        if sample_names in D[workflow_id]:
            D[workflow_id][sample_names].append([file, fileswid])
        else:
            D[workflow_id][sample_names] = [[file, fileswid]]
    return D



def map_fileswid_to_filename(project_name, database, table='Files'):
   '''


   '''

   # get the workflow output files sorted by sample
   conn = connect_to_db(database)
   query = "SELECT DISTINCT file_swid, file FROM {0} WHERE project_id = ?;".format(table)
   data = conn.execute(query, (project_name,)).fetchall()
   conn.close()

   D = {}
   for i in data:
       D[i['file_swid']] = i['file']
   
   return D



def get_WGS_standard_deliverables():
    '''
    (None) -> dict
    
    Returns a dictionary with the file extensions or file endings for each workflow
    for which output files are released as part of the standard WGS package
    
    Parameters
    ----------
     None
    
    '''
    
    deliverables = {'bammergepreprocessing': ['.bai', '.bam'],
                    'varianteffectpredictor': ['.mutect2.filtered.vep.vcf.gz',
                                               '.mutect2.filtered.vep.vcf.gz.tbi',
                                               '.mutect2.filtered.maf.gz'],
                    'delly': ['.somatic_filtered.delly.merged.vcf.gz',
                      '.somatic_filtered.delly.merged.vcf.gz.tbi'],
                    'sequenza': ['results.zip', 'summary.pdf', 'alternative_solutions.json'],
                    'mavis': ['.tab', '.zip']}
    
    return deliverables

    
    
def get_contamination(sample_id, database, table = 'Calculate_Contamination'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with call-ready contamination and merged limskey for sample_id     
    
    Parameters
    ----------
    - sample_id (str): Sample identifier
    - database (str): Path to the sqlite database
    - table (str): Table in database storing the call-ready contamination. Default is Calculate_Contamination
    '''    
   
    conn = connect_to_db(database)
    query = "SELECT DISTINCT contamination, merged_limskey FROM {0} WHERE sample_id = ?;".format(table)
    data = conn.execute(query, (sample_id,)).fetchall()
    conn.close()

    D = {}
    for i in data:
        if i['merged_limskey'] in D:
            D[i['merged_limskey']].append(i['contamination'])
        else:
            D[i['merged_limskey']] = [i['contamination']]
    for i in D:
        D[i] = max(D[i])    
    
    return D



def group_limskeys(block_limskeys):
    '''
    (list) -> list
    
    Sort the limskeys of an analysis block by sample
    
    Parameters
    ----------
    - block_limskeys (list): List of limskeys for a given block
        
    Examples
    --------
    >>> group_limskeys(['4991_1_LDI51430', '5073_4_LDI57812', '5073_3_LDI57812', '5073_2_LDI57812'])
    ['4991_1_LDI51430', '5073_4_LDI57812;5073_3_LDI57812;5073_2_LDI57812']
    '''
    
    D = {}
    for i in block_limskeys:
        j = i.split('_')[-1]
        if j in D:
            D[j].append(i)
        else:
            D[j] = [i]
        D[j].sort()
    
    L = [';'.join(D[j]) for j in D]

    return L


def get_block_level_contamination(project_name, database, blocks, sample_pair):
    '''
    (str, str,, dict, str) -> dict
    
    Returns a dictionary with contamination for each anchor workflow in the 
    analysis block for the given sample pair.
    The contamination is the maximum contamination among samples and among 
    bmpp workflows for each anhor workflow (sub-block)
    
    Parameters
    ----------
    - project_name (str): 
    - database (str): Path to the sqlite database
    - blocks (dict): Dictionary with analysis block information
    - sample_pair (str): Normal/tumor sample pair
    '''
       
    # get the block limskeys
    limskeys = get_workflow_limskeys(project_name, database, 'Workflow_Inputs')
    
    # list the bmpp anchor workflow ids
    bmpp_workflows = [blocks[sample_pair][i]['anchor_wf'] for i in range(len(blocks[sample_pair]))]
    bmpp_workflows = list(set(bmpp_workflows))
    
    # get contamination for each sample in sample pair
    contamination = {}
    for sample in sample_pair.split('|'):
        sample = sample.strip()
        contamination.update(get_contamination(sample, database, 'Calculate_Contamination'))
    
    # map each contamination to the workflow anchor id
    D = {}    
    for workflow in bmpp_workflows:
        block_conta = []
        for workflow_id in workflow.split('.'):
            # get the limskeys for that workflow
            workflow_limskeys = limskeys[os.path.basename(workflow_id)]
            # group the limskeys by sample
            workflow_limskeys = group_limskeys(workflow_limskeys)
            # use the maximum contamination of each sample 
            conta = [contamination[i] for i in workflow_limskeys if i in contamination]
            if conta:
                conta = round(max(conta) * 100, 3)
            else:
                conta = 'NA'
            block_conta.append(conta)
        # use the maximum contamination of each bmpp workflow
        if 'NA' in block_conta:
            D[workflow] = 'NA'
        else:
            D[workflow] = max(block_conta)
    
    return D       
    


def map_workflow_to_platform(project_name, database, table = 'Workflow_Inputs'):
    '''
    (str, str, str) -> dict
    
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table storing workflow_input information
    '''
    
    # connect to db
    conn = connect_to_db(database)
    # extract library source
    data = conn.execute("SELECT DISTINCT wfrun_id, limskey, platform FROM {} WHERE project_id = ?;".format(table), (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        workflow = i['wfrun_id']
        limskey = i['limskey']
        platform = i['platform']
        
        if workflow in D:
            D[workflow]['limskey'].add(limskey)
            D[workflow]['platform'].add(platform)
        
        else:
            D[workflow] = {'limskey': {limskey}, 'platform' : {platform}}
    
    return D
    

def get_sample_sequencing_amount(project_name, case, samples, database, workflow_table = 'Workflows', workflow_input = 'Workflow_Inputs', library_table = 'Libraries'):
    '''
    (str, str, str, str, str, str, str) -> dict

    Returns a dictionary with the lane counts for sequencing workflows
    for each sequencing platform for each sample in samples 

    Parameters
    ----------
    - project_name (str): Name of project of interest 
    - case (str): Donor identifier
    - samples (str): Sample or sample pair
    - database (str): Path to the sqlite database
    - workflow_table (str): Table storing the workflow information
    - workflow_input (str): Table storing the workflow input information
    - library_table (str): Table storing the library information
    '''    
        
    libraries = map_limskey_to_library(project_name, database, table = workflow_input)
    sample_names = map_library_to_sample(project_name, database, table = library_table)
    workflow_names = get_workflow_names(project_name, database)
    amount_data = get_amount_data(project_name, database, workflow_table)
    platforms = map_workflow_to_platform(project_name, database, table = workflow_input)

    sequencing_workflows = ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq')
    workflows = [i for i in workflow_names if workflow_names[i][0].lower() in sequencing_workflows]

    # find libraries for each sample
    L = {}
    for sample in samples.split('|'):
        sample = sample.strip()
        L[sample] = []
        for i in sample_names[case]:
            if sample_names[case][i] == sample:
                L[sample].append(i)
    # find sequencing workflows for each sample
    wfs = {}
    for sample in L:
        for workflow in libraries:
            for limskey in libraries[workflow]:
                if workflow in workflows and libraries[workflow][limskey] in L[sample]:
                    if sample not in wfs:
                        wfs[sample] = []
                    wfs[sample].append(workflow)
          
    # sort lane counts, workflow and limskey by sample and platform    
    lanes = {}
    for sample in wfs:
        for workflow in wfs[sample]:
            platform = platforms[workflow]['platform']
            assert len(platform) == 1
            platform = list(platform)[0]
            if sample not in lanes:
                lanes[sample] = {}
            if platform not in lanes[sample]:
                lanes[sample][platform] = {'count': 0, 'workflows': [], 'released': [], 'limskeys': []}
            lanes[sample][platform]['count'] += amount_data[workflow]    
            lanes[sample][platform]['workflows'].append(workflow)  
    
    return lanes    
       









def get_cases_with_analysis(analysis_db, project_name, assay):
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
    data = conn.execute("SELECT case_id, donor, template, valid, error, md5sum FROM templates WHERE \
                        project = ? AND assay = ?;", (project_name,assay)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        case = i['case_id']
        md5sum = i['md5sum']
        template = json.loads(i['template'])
        valid = int(i['valid'])
        donor = i['donor']
        error = i['error']
        
        d = {'md5sum': md5sum, 'template': template, 'valid': valid, 'error': error, 'donor': donor}        
        
        if case in D:
            D[case].append(d)
        else:
            D[case] = [d]
        
    return D    



def get_case_error_message(cases):
    '''
    (dict) -> dict
    
    Returns a dictionary with the error messages across all templates for each case
    
    Parameters
    ----------
    - cases (dict): Dictionary of cases with analysis data corresponding to project and assay
    '''

    D = {}
    
    for case in cases:
        for d in cases[case]:
            error = d['error']
            error = error.split(';')
            if case in D:
                D[case].extend(error)
            else:
                D[case] = error
    for case in D:
        D[case] = ';'.join(sorted(list(set(D[case]))))

    return D



def delete_cases_with_distinct_checksums(cases, md5sums):
    '''
    (dict, dict) -> dict
    
    Removes cases when databases are in sync (ie, cases have different md5sums)
    and modify cases in places
       
    Parameters
    ----------
    - cases (dict): Dictionary with analysis data from the analysis review database
    - md5sums (dict): Dictionary with md5sums for each case in the main waterzooi database
    '''
    
    # keep only cases with up to date data between resources
    to_remove = []
    for i in cases:
        for j in cases[i]:
            if j['md5sum'] != md5sums[i]:
                to_remove.append(i)
    to_remove = list(set(to_remove))
    if to_remove:
        alert = 'removing {0} cases for which data is not in sync between waterzooi and analysis databases'
        print(alert.format(len(to_remove)))
        for i in to_remove:
            del cases[i]



def map_donors_to_cases(cases):
    '''
    (dict) -> dict
    
    Returns a dictionary with case and corresponding donor identifier
        
    Parameters
    ----------
    - cases (dict): Dictionary with analysis data from the analysis review database
    '''

    # get the donor
    donors = {}
    for i in cases:
        for d in cases[i]:
            donor = d['donor']
            if i in donors:
                assert donor == donors[i]
            else:
                donors[i] = donor

    return donors



def get_case_analysis_samples(cases):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the samples with analysis data for each case
    
    Parameters
    ----------
    - cases (dict): Dictionary of project cases with analysis data from a specific assay
    '''
    
    D = {}
    for case in cases:
        samples = []
        for d in cases[case]:
            if d['template']:
                for i in d['template']['Samples']:
                    samples.append(d['template']['Samples'][i]['sample']) 
        D[case] = list(set(samples))
        
    return D
    


# def count_case_analysis_workflows(cases):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with analysis workflow ids organized for each case
    
#     Parameters
#     ----------
#     - cases (dict): Dictionary with case analysis extracted from the analysis review database
#     '''
        
#     D = {}
    
#     for case in cases:
#         if case not in D:
#             D[case] = {}
#         callready = []
#         for d in cases[case]:
#             if 'Anchors' in d['template']:
#                 for i in d['template']['Anchors']:
#                     callready.append(d['template']['Anchors'][i]['workflows'])
#         D[case]['callready'] = len(list(set(callready)))
#         downstream = []
#         for d in cases[case]:
#             if 'Analysis' in d['template']:
#                 for i in d['template']['Analysis']:
#                     for j in d['template']['Analysis'][i]:
#                         downstream.append(j['workflow_id'])
#         D[case]['downstream'] = len(list(set(downstream)))
#         analysis = {}
#         for d in cases[case]:
#             if 'Analysis' in d['template']:
#                 for i in d['template']['Analysis']:
#                     if i not in analysis:
#                         analysis[i] = []
#                     for j in d['template']['Analysis'][i]:
#                         analysis[i].append(j['workflow_id'])
#                     analysis[i] = list(set(analysis[i]))   
#         for i in analysis:
#             analysis[i] = len(analysis[i])
#         D[case]['analysis'] = analysis

#     return D






def get_case_analysis_workflows(cases):
    '''
    (dict) -> dict
    
    Returns a dictionary with analysis workflow ids organized for each case
    
    Parameters
    ----------
    - cases (dict): Dictionary with case analysis extracted from the analysis review database
    '''
        
    D = {}
    
    for case in cases:
        if case not in D:
            D[case] = []
        for d in cases[case]:
            callready = []
            if 'Anchors' in d['template']:
                for i in d['template']['Anchors']:
                    callready.append(d['template']['Anchors'][i]['workflows'])
            callready = list(set(callready))
            downstream = []
            if 'Analysis' in d['template']:
                for i in d['template']['Analysis']:
                    for j in d['template']['Analysis'][i]:
                        downstream.append(j['workflow_id'])
            downstream = list(set(downstream))
            sequencing = {}
            if 'Data' in d['template']:
                for i in d['template']['Data']['Sequencing']['workflows']:
                    workflow_name = i['workflow_name']
                    workflow_id = i['workflow_id']
                    if workflow_name in sequencing:
                        sequencing[workflow_name].append(workflow_id)
                    else:
                        sequencing[workflow_name] = [workflow_id]
                sequencing[workflow_name] = list(set(sequencing[workflow_name]))
            analysis = {}
            if 'Analysis' in d['template']:
                for i in d['template']['Analysis']:
                    if i not in analysis:
                        analysis[i] = []
                    for j in d['template']['Analysis'][i]:
                        analysis[i].append(j['workflow_id'])
                analysis[i] = list(set(analysis[i]))   
            alignments = {}
            if 'Data' in d['template']:
                for i in d['template']['Data']:
                    if i != 'Sequencing':
                        for k in d['template']['Data'][i]['workflows']:
                            wfname = k['workflow_name']
                            workflow_id = k['workflow_id']   
                            if wfname not in alignments:
                                alignments[wfname] = []
                            alignments[wfname].append(workflow_id)
                            alignments[wfname] = list(set(alignments[wfname]))
                        
            D[case].append({'callready': callready, 'downstream': downstream,
                            'sequencing': sequencing, 'analysis': analysis,
                            'alignments': alignments})
        
    return D




def count_case_analysis_workflows(case_analysis_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with analysis workflow ids organized for each case
    
    Parameters
    ----------
    - cases (dict): Dictionary with case analysis extracted from the analysis review database
    '''
        
    D = {}
    
    for case in case_analysis_data:
        if case not in D:
            D[case] = {'callready': [], 'downstream': [], 'analysis': {}}
        for d in case_analysis_data[case]:
            D[case]['callready'].extend(d['callready'])
            D[case]['downstream'].extend(d['downstream'])
            for i in d['analysis']:
                if i not in D[case]['analysis']:
                    D[case]['analysis'][i] = []
                D[case]['analysis'][i].extend(d['analysis'][i])
            
        D[case]['callready'] = len(set(D[case]['callready']))
        D[case]['downstream'] = len(set(D[case]['downstream']))
        for i in D[case]['analysis']:
            D[case]['analysis'][i] = len(set(D[case]['analysis'][i]))
            
    return D





def get_sequencing_input(database, case):
    '''
    (str, str) -> dict
        
    Returns a dictionary with limskey and library for each sequencing workflow id of a case
    
    Paramaters
    ----------
    - database (str): Path to the main waterzooi database
    - case (str): Case identifier
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Workflow_Inputs.limskey, Workflow_Inputs.library, \
                        Workflow_Inputs.wfrun_id FROM Workflow_Inputs JOIN Workflows \
                        WHERE  Workflow_Inputs.wfrun_id = Workflows.wfrun_id \
                        AND LOWER(Workflows.wf) IN ('casava', 'bcl2fastq', 'fileimportforanalysis', \
                        'fileimport', 'import_fastq', 'bwamem') AND Workflow_Inputs.case_id = ?;", (case,)).fetchall()
    conn.close()

    D = {}
    for i in data:
        limskey = i['limskey']
        library = i['library']
        wfrun_id = i['wfrun_id']
        
        assert wfrun_id not in D
        D[wfrun_id] = {'limskey': limskey, 'library': library}
       
    return D










    


def list_assay_analysis_workflows(case_data):
    '''
    (dict) -> list
        
    Returns a list of all the expected analysis workflows for an essay
    
    Parameters
    ----------
    - (case_data): Dictionary with analysis workflow ids organized for each case with the same assay
    '''
    
    # get all the analysis workflows
    analysis_workflows = []
    for case in case_data:
        for d in case_data[case]:
            analysis_workflows.extend(list(d['analysis'].keys()))
    analysis_workflows = sorted(list(set(analysis_workflows)))    
    
    return analysis_workflows




def most_recent_analysis_workflow(case_data, creation_dates):
    '''
    (list, dict) -> list
    
    Returns a list with the most recent workflow for each analysis template in a case
       
    Parameters
    ----------
    - case_data (list): List of templates with analysis data for a single case
    - creation_dates (dict): Dictionary with creation dates of each workflow
    '''
    
    most_recent = []
       
    for template in case_data:
        L = []
        for workflow_id in template['callready']:
            #L.append(creation_dates[os.path.basename(workflow_id)])
            L.append(creation_dates[workflow_id])
        for workflow_id in template['downstream']:
            L.append(creation_dates[workflow_id])
        
        L.sort()
        try:
            date = time.strftime('%Y-%m-%d', time.localtime(int(L[-1])))
        except:
            date = 'NA'
        most_recent.append(date)
        
    return most_recent


def get_analysis_workflow_name(analysis):
    '''
    (dict) -> dict
    
    Returns a dictionary with the name of workflows for each workflow id of
    an analysis group of a single case
    
    Parameters
    ----------
    - analysis (dict): Dictionary with workflow ids of analysis workflows of a single case
    '''
    
    D = {}
    
    for workflow_name in analysis:
        for workflow_id in analysis[workflow_name]:
            assert workflow_id not in D
            D[workflow_id] = workflow_name
    
    return D
    
    
def get_case_workflow_samples(database, case):
    '''
    (str, str) -> dict
    
    Returns a dictionary of workflow ids and list of corresponding samples for a single case
    
    Parameters
    ----------
    - database (str): Path to the database
    - case (str): Case of interest
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Workflow_Inputs.wfrun_id, Libraries.donor_id, Libraries.tissue_type, \
                        Libraries.tissue_origin, Libraries.library_type, Libraries.group_id FROM \
                        Workflow_Inputs JOIN Libraries WHERE Workflow_Inputs.library=Libraries.library \
                        AND Libraries.case_id = ?;", (case,)).fetchall()     
    conn.close()
    
    D = {}
    for i in data:
        donor = i['donor_id']
        tissue_origin = i['tissue_origin']
        tissue_type = i['tissue_type']
        library_type = i['library_type']
        groupid = i['group_id'] 
        wfrunid = i['wfrun_id']
        sample = '_'.join([donor, tissue_origin, tissue_type, library_type, groupid]) 
        
            
        if wfrunid in D:
            D[wfrunid].append(sample)
        else:
            D[wfrunid] = [sample]
        D[wfrunid] = list(set(D[wfrunid]))

    return D    

    
    
    
def get_assays(database, project_name):
    '''
    (str, str) -> str
    
    Returns a comma-separated list of all assays for a given project 
    
    Parameters
    ----------
    - database (str): Path to the database
    - project_name (str): Name of project of interest
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT assays FROM Projects WHERE project_id = ?;", (project_name,)).fetchall()     
    conn.close()
    
    assays = ','.join([i['assays'] for i in data])
    
    return assays


def get_missing_workflows(case_data):
    '''
    (list) -> list
    
    Returns a list of list of mising workflows for each template of a single case
    
    Parameters
    ----------
    - case_data (list): List of dictionaries with template of analysis data for a single case
    '''
    
    missing = []
    
    for d in case_data:
        L = []
        for workflow in d['analysis']:
            if len(d['analysis'][workflow]) == 0:
                L.append(workflow)
        missing.append(L)
        
    return missing
    
    
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
    


# def map_files_to_limskeys(database, case):
#     '''
    
    
    
#     '''
    
#     D = {}
    
    
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT file_swid, wfrun_id, file, limskey FROM Files \
#                         WHERE case_id = ?", (case,)).fetchall()
#     conn.close()
                        
#     for i in data:
#         limskeys = i['limskey'].split(';')                    
#         file_swid = i['file_swid']
#         file = i['file']
#         wfrunid = i['wfrun_id']                
                        
#         for limskey in limskeys:
#             if wfrunid not in D:
#                 D[wfrunid] = {}
#             if limskey in D[wfrunid]:
#                 D[wfrunid][limskey].append({'file': file, 'file_swid': file_swid})
#             else:
#                 D[wfrunid][limskey] = [{'file': file, 'file_swid': file_swid}]    
                
#     return D                        
                        
    

# def map_limskey_to_sample(database, case):
#     '''
    
    
#     '''
    
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT Workflow_Inputs.wfrun_id, Workflow_Inputs.limskey, \
#                         Libraries.donor_id, Libraries.tissue_type, \
#                         Libraries.tissue_origin, Libraries.library_type, Libraries.group_id FROM \
#                         Workflow_Inputs JOIN Libraries WHERE Workflow_Inputs.library=Libraries.library \
#                         AND Libraries.case_id = ?;", (case,)).fetchall()     
#     conn.close()
    
#     D = {}
#     for i in data:
#         limskey = i['limskey']
#         donor = i['donor_id']
#         tissue_origin = i['tissue_origin']
#         tissue_type = i['tissue_type']
#         library_type = i['library_type']
#         groupid = i['group_id'] 
#         wfrunid = i['wfrun_id']
#         sample = '_'.join([donor, tissue_origin, tissue_type, library_type, groupid]) 
        
#         if wfrunid not in D:
#             D[wfrunid] = {}
#         if limskey in D[wfrunid]:
#             assert sample == D[wfrunid][limskey]
#         else:
#             D[wfrunid][limskey] = sample
            
#     return D    









# def get_workflow_output_files(files_to_limskeys, limskey_to_sample):
#     '''
    
    
#     '''
    
#     D = {}
    
#     for wfrunid in files_to_limskeys:
#         for limskey in files_to_limskeys[wfrunid]:
            
#             # this should be always be true - let's see the new olive
#             if wfrunid in limskey_to_sample and limskey in limskey_to_sample[wfrunid]:
            
                
#                 sample = limskey_to_sample[wfrunid][limskey]
#                 if wfrunid not in D:
#                     D[wfrunid] = {}
#                 for d in files_to_limskeys[wfrunid][limskey]:
#                     file = d['file']
#                     file_swid = d['file_swid']
            
#                     if sample not in D[wfrunid]:
#                         D[wfrunid][sample] = {}
#                     if file in D[wfrunid][sample]:
#                         assert limskey not in D[wfrunid][sample][file]['limskey']
#                         D[wfrunid][sample][file]['limskey'].append(limskey)
#                     else:
#                         D[wfrunid][sample][file] = {'file': file, 'file_swid': file_swid, 'limskey': [limskey]}
   
#     return D
    
    
# def map_samples_to_files(outputfiles):
#     '''
#     (dict) -> dict
    
    
#     '''

#     D = {}
    
#     for sample in outputfiles:
#         for file in outputfiles[sample]:
#             if file in D:
#                 D[file].append(sample)
#             else:
#                 D[file] = [sample]
#             D[file] = sorted(list(set(D[file])))
    
#     return D



# def group_files_by_samples(files_to_samples):
#     '''
    
    
#     '''
    
#     D = {}
    
#     for file in files_to_samples:
#         sample = ';'.join(files_to_samples[file])
#         if sample in D:
#             D[sample].append(file)
#         else:
#             D[sample] = [file]
 
#     return D







def get_workflow_output_files(database, wfrun_id):
    '''
    (str, str) ->
    
    Returns a dictionary with the output files of workflow with wfrun_id grouped by sample 
    
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
    
    for i in data:
        sample = i['sample_id']
        file = i['file']
        if file in D:
            D[file].append(sample)
        else:
            D[file] = [sample]
        D[file] = sorted(list(set(D[file])))
            
    # group samples sharing the same files
    S = {}
    for file in D:
        sample = ';'.join(D[file])
        if sample in S:
            S[sample].append(file)
        else:
            S[sample] = [file]
       
    return S




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


def get_input_sequences(database, case, wfrun_id):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with input sequences  of worflow with identifier wfrun_id
    
    Parameters
    ----------
    - database (str): Path to the waterzooi
    - case (str): Case identifier
    - wfrun_id (str): Workflow unique identifier
    '''
    
    # get the limskeys matching the workflow
    limskeys = map_limskeys_to_workflow(database, wfrun_id)

    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Files.file_swid, Files.file, Files.limskey, Libraries.library, \
                        Libraries.sample_id FROM Files JOIN Libraries JOIN Workflows \
                        WHERE Files.wfrun_id = Workflows.wfrun_id AND Files.limskey = Libraries.lims_id \
                        AND LOWER(Workflows.wf) IN ('casava', 'bcl2fastq', 'fileimportforanalysis', \
                        'fileimport', 'import_fastq') AND Files.case_id = ?;", (case,)).fetchall()
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




def get_pipeline_standard_deliverables():
    '''
    (None) -> dict
    
    Returns a dictionary with the file extensions or file endings for each workflow
    for which output files are released as part of the standard WGS package
    
    Parameters
    ----------
     None
    
    '''
    
    deliverables = {'bammergepreprocessing': ['.bai', '.bam'],
                    'varianteffectpredictor': ['.mutect2.filtered.vep.vcf.gz',
                                               '.mutect2.filtered.vep.vcf.gz.tbi',
                                               '.mutect2.filtered.maf.gz'],
                    'delly': ['.somatic_filtered.delly.merged.vcf.gz',
                      '.somatic_filtered.delly.merged.vcf.gz.tbi'],
                    'sequenza': ['results.zip', 'summary.pdf', 'alternative_solutions.json'],
                    'mavis': ['.tab', '.zip'],
                    'star': ['.bai', '.bam'],
                    'rsem': ['.genes.results', '.isoforms.results', '.transcript.bam'],
                    'starfusion': ['.tsv'],
                    'arriba': ['.tsv', '.fusions.pdf'],
                    'hrdetect': ['.signatures.json'],
                    'msisensor': ['.msi', '.msi_germline', '.msi_somatic', '.msi.booted'],
                    'purple': ['.purple.purity.tsv', '.purple.purity.qc', '.purple.qc',
                               '.purple.purity.range.tsv', '.purple.cnv.somatic.tsv',
                               'purple.cnv.gene.tsv', '.purple.segment.tsv',
                               '.solPrimary.purple.zip', '.purple_alternates.zip',
                               '.purple.somatic.vcf.gz'],
                    'gridds': ['.purple.sv.vcf.gz']}
    
    return deliverables



### rreview beloe for downloading the json


def create_analysis_json(case_data, selected_workflows, workflow_outputfiles, deliverables=None):
    '''
    (dict, dict, dict, dict, dict, dict, None | dict)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor bmpp parent workflow
    
    Parameters
    ----------
    - case_data (dict): Dictionary with analysis templates for all cases in a project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - workflow_outputfiles (dict): Dictionary with outputfiles for each workflow run
    - deliverables (None | dict): None or dictionary with file extensions of standard deliverables
    '''
        
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
   
    D = {}
        
    for case in case_data:
        for template in case_data[case]:
            # make a list of workflows:
            callready = template['callready']
            workflows = {}
            for i in ['sequencing', 'analysis', 'alignments']:
                for workflow in template[i]:
                    if template[i][workflow]:
                        workflows[workflow] = template[i][workflow]
        
            # check that analysis workflows are selected
            # do not include call ready workflows because they can be shared across templates
            wfs = []
            for i in workflows.values():
                wfs.extend(i)
            
            analysis = [i for i in wfs if i not in callready]   
            if any(map(lambda x: x in selected_workflows, analysis)):
                for workflow in workflows:
                    for wfrunid in workflows[workflow]:
                        # check if workflow is selected
                        if selected_workflows[wfrunid]:
                            # get workflow output files
                            outputfiles = workflow_outputfiles[wfrunid]                        
                                        
                            # check that only workflows in standard eliverables are used
                            if deliverables:
                                key = workflow.split('_')[0].lower()
                                if key in deliverables:
                                    # map all file endings of deliverables with files
                                    groups = list(itertools.product(outputfiles, deliverables[key]))
                                    # determine which files are part of the deliverables
                                    F = list(map(G, groups))
                                    L = [groups[k][0] for k in range(len(F)) if F[k]]
                                    if L:
                                        if case not in D:
                                            D[case] = {}
                                        if workflow in D[case]:
                                            D[case][workflow].extend(L)
                                        else:
                                            D[case][workflow] = L
                            else:
                                if case not in D:
                                    D[case] = {}
                                if workflow in D[case]:
                                    D[case][workflow].extend(outputfiles)
                                else:
                                    D[case][workflow] = outputfiles
                            
                            D[case][workflow] = sorted(list(set(D[case][workflow])))  
    
    
    return D
                    



def create_case_analysis_json(case, case_data, selected_workflows, workflow_outputfiles, selection):
    '''
    (str, list, dict, dict, str)
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor parent workflow (bmpp or star)
    
    Parameters
    ----------
    - case (str): Case unique identifier
    - case_data (list): List of analysis templates for a single case in a project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - workflow_outputfiles (dict): Dictionary with outputfiles for each workflow run
    - selection (str): Include files from all selected workflows or files from the standard deliverables
                       Values: standard or all
    '''
    
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
            
    # get the deliverables
    if selection == 'standard':
        deliverables = get_pipeline_standard_deliverables()
    elif selection == 'all':
        deliverables = {}
    
    D = {}
            
    for template in case_data:
        # make a list of workflows:
        callready = template['callready']
        workflows = {}
        for i in ['sequencing', 'analysis', 'alignments']:
            for workflow in template[i]:
                if template[i][workflow]:
                    workflows[workflow] = template[i][workflow]
            
        # check that analysis workflows are selected
        # do not include call ready workflows because they can be shared across templates
        wfs = []
        for i in workflows.values():
            wfs.extend(i)
                
        analysis = [i for i in wfs if i not in callready]   
        if any(map(lambda x: x in selected_workflows, analysis)):
            for workflow in workflows:
                for wfrunid in workflows[workflow]:
                    # check if workflow is selected
                    if selected_workflows[wfrunid]:
                        # get workflow output files
                        outputfiles = workflow_outputfiles[wfrunid]                        
                                            
                        # check that only workflows in standard eliverables are used
                        if deliverables:
                            key = workflow.split('_')[0].lower()
                            if key in deliverables:
                                # map all file endings of deliverables with files
                                groups = list(itertools.product(outputfiles, deliverables[key]))
                                # determine which files are part of the deliverables
                                F = list(map(G, groups))
                                L = [groups[k][0] for k in range(len(F)) if F[k]]
                                if L:
                                    if case not in D:
                                        D[case] = {}
                                    if workflow in D[case]:
                                        D[case][workflow].extend(L)
                                    else:
                                        D[case][workflow] = L
                        else:
                            if case not in D:
                                D[case] = {}
                            if workflow in D[case]:
                                D[case][workflow].extend(outputfiles)
                            else:
                                D[case][workflow] = outputfiles
                        
                        D[case][workflow] = sorted(list(set(D[case][workflow])))  
                       
        
        
        return D
    
    
    
    
    


# def get_selected_workflows(project_name, database, table = 'Workflows'):
#     '''
#     (str, str, str) -> dict 
    
#     Returns a dictionary with the selected status of each workflow for the given project
                  
#     Parameters
#     ----------
#     - project_name (str): Name of project of interest
#     - database (str): Path to the sqlite database
#     - table (str): Table with workflow information
#     '''
                
#     conn = connect_to_db(database)
#     query = "SELECT DISTINCT wfrun_id, selected FROM {0} WHERE project_id = ?;".format(table)
#     data = conn.execute(query, (project_name,)).fetchall() 
#     conn.close()

#     D = {}
#     for i in data:
#         D[i['wfrun_id']] = int(i['selected'])
    
#     return D



# def review_wgs_blocks(blocks, selected_workflows):
#     '''
#     (dict, dict) -> dict 
    
#     Returns a dictionary with status for analysis blocks for each case in project
                  
#     Parameters
#     ----------
#     - blocks (dict): 
#     - selected_workflows (dict): 
#     '''
    
#     D = {}
        
#     for case in blocks:
#         if case not in D:
#             D[case] = {}
#         for samples in blocks[case]:
#             for anchor in blocks[case][samples]:
#                 # do not include call-ready workflows to determine selection/review status
#                 # these may be shared across multiple blocks
#                 L = [selected_workflows[os.path.basename(i)] for i in blocks[case][samples][anchor]['workflows'] if i not in anchor] 
                
#                 if any(L):
#                     D[case][samples] = anchor
#                     break
#                 else:
#                     if blocks[case][samples][anchor]['clean'] and \
#                       blocks[case][samples][anchor]['complete']:
#                           D[case][samples] = 'ready'
#                           break
#                     else:
#                         D[case][samples] = 'review'
#     return D




# def get_wgs_blocks(project, database, table = 'WGS_blocks'):
#     '''
#     (str, str, str) -> dict
    
#     Returns a dictionary with analysis block names for each case in project
    
#     Parameters
#     ----------
#     - project (str): Name of project of interest
#     - database (str): Path to the sqlite database
#     - table (str): Table with analysis blocks
#     '''

#     conn = connect_to_db(database)    
#     query = "SELECT DISTINCT * FROM {0} WHERE project_id = ?;".format(table)
#     data = conn.execute(query, (project,)).fetchall()
#     conn.close()
    
#     D = {}
    
#     for i in data:
#         case, samples, anchor = i['case_id'], i['samples'], i['anchor_wf']
#         workflows = i['workflows'].split(';')
#         clean, complete = int(i['clean']), int(i['complete'])
#         if case not in D:
#             D[case] = {}
#         if samples not in D[case]:
#             D[case][samples] = {}
#         assert anchor not in D[case][samples]
#         D[case][samples][anchor] = {'workflows': workflows, 'clean': clean,
#                                     'complete': complete}
            
#     return D







# def create_WGS_project_block_json(project_name, database, blocks, block_status, selected_workflows, workflow_names, deliverables=None):
#     '''
#     (str, str, dict, dict, dict, dict, None | dict)
    
#     Returns a dictionary with workflow information for a given block (ie, sample pair)
#     and anchor bmpp parent workflow
    
#     Parameters
#     ----------
#     - project_name (None | str): None or name of project of interest
#     - database (None | str): None or path to the sqlite database
#     - blocks (dict): Dictionary with block information
#     - block_status (dict): Dictionary with review status of each block
#     - selected_workflows (dict): Dictionary with selected status of each workflow in project
#     - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
#     - deliverables (None | dict): None or dictionary with file extensions of standard WGS deliverables
#     '''
    
#     libraries = map_limskey_to_library(project_name, database, table='Workflow_Inputs')
#     sample_names = map_library_to_sample(project_name, database, table = 'Libraries')
#     donors = map_library_to_case(project_name, database, table = 'Libraries')
#     workflow_outputfiles = get_workflow_output(project_name, database, libraries, sample_names, donors, 'Files')
      
    
#     # create a lambda to evaluate the deliverable files
#     # x is a pair of (file, file_ending)
#     G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
    
#     D = {}
#     for case in blocks:
#         for samples in blocks[case]:
#             # check the selection status of the block
#             if block_status[case][samples] not in ['ready', 'review']:
#                 # block already reviewed and workflows selected
#                 anchor_wf = block_status[case][samples]
#                 for workflow in blocks[case][samples][anchor_wf]['workflows']:
#                     workflow = os.path.basename(workflow)
                    
#                     # get workflow name and version
#                     workflow_name = workflow_names[workflow][0]
#                     workflow_version = workflow_names[workflow][1]
                               
#                     # check workflow status
#                     if selected_workflows[workflow]:
#                         # get workflow output files
#                         # needed to sort outputs by sample pairs or by sample for call-ready workflows
#                         # even if all files are recorded
#                         #outputfiles = get_workflow_output(project_name, case, workflow, database, libraries, sample_names, 'Files')
#                         outputfiles = workflow_outputfiles[workflow]                        
                        
#                         # check that only workflows in standard WGS deliverables are used
#                         if deliverables:
#                             key = workflow_names[workflow][0].split('_')[0].lower()
                                                       
#                             if key in deliverables:
#                                 for j in outputfiles:
#                                     # list all deliverable files
#                                     L = []
#                                     # gather all file paths for workflow and sample(s)
#                                     files = [i[0] for i in outputfiles[j]]
#                                     # map all file endings of deliverables with files
#                                     groups = list(itertools.product(files, deliverables[key]))
#                                     # determine which files are part of the deliverables
#                                     F = list(map(G, groups))
#                                     L = [groups[k][0] for k in range(len(F)) if F[k]]
                                    
#                                     if L:
#                                         sample_id = j.replace(';', '.')
#                                         if case not in D:
#                                             D[case] = {}
#                                         if sample_id not in D[case]:
#                                             D[case][sample_id] = {}
#                                         if workflow_name not in D[case][sample_id]:
#                                             D[case][sample_id][workflow_name] = []
                                        
#                                         d = {'workflow_id': workflow,
#                                              'workflow_version': workflow_version,
#                                              'files': L}    
#                                         if d not in D[case][sample_id][workflow_name]: 
#                                             D[case][sample_id][workflow_name].append(d)
                        
#                         else:
#                             for j in outputfiles:
#                                 sample_id = j.replace(';', '.')
#                                 d =  {'workflow_id': workflow,
#                                       'workflow_version': workflow_version,
#                                       'files': [i[0] for i in outputfiles[j]]}
#                                 if case not in D:
#                                     D[case] = {}
#                                 if sample_id not in D[case]:
#                                     D[case][sample_id] = {}
#                                 if workflow_name not in D[case][sample_id]:
#                                     D[case][sample_id][workflow_name] = []
#                                 if d not in D[case][sample_id][workflow_name]:
#                                     D[case][sample_id][workflow_name].append(d)
                                        
    
#     return D


def get_selected_workflows(project_name, database, table = 'Workflows'):
    '''
    (str, str, str) -> dict 
    
    Returns a dictionary with the selected status of each workflow for the given project
                  
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    - table (str): Table with workflow information
    '''
                
    conn = connect_to_db(database)
    query = "SELECT DISTINCT wfrun_id, selected FROM {0} WHERE project_id = ?;".format(table)
    data = conn.execute(query, (project_name,)).fetchall() 
    conn.close()

    D = {}
    for i in data:
        D[i['wfrun_id']] = int(i['selected'])
    
    return D


def get_workflow_outputfiles(database, project_name):
    '''
    (str, str) ->
    
    Returns a dictionary with the workflow output files organuzed by cases of project
        
    Parameters
    ----------
    - database (str): Path to the database
    - project_name (str): Project nameWorkflow run identifier
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT Files.file, Files.wfrun_id, Files.case_id, \
                        Workflows.wf FROM Files JOIN Workflows WHERE Workflows.wfrun_id = Files.wfrun_id \
                        AND Workflows.case_id = Files.case_id AND Files.project_id = Workflows.project_id \
                        AND Files.project_id = ?;", (project_name, )).fetchall()
    conn.close()   
    
    D = {}
    
    for i in data:
        wfrun_id = i['wfrun_id']
        case_id = i['case_id']
        file = i['file']
        workflow = i['wf']
        
        if wfrun_id in D:
            D[wfrun_id].append(file)
        else:
            D[wfrun_id] = [file]
        
    return D




def get_review_status(case_data, selected_workflows):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with review status for each case in project
    A case is considered to be reviewed if any worklow has been selected
    
    Parameters
    ----------
    - case_data (dict): Dictionary with analysis template for each case
    - selected_workflows (dict): Dictionary with selection status of each workflow in a project
    '''
    
    D = {}
    
    for case in case_data:
        status = 0
        for template in case_data[case]:
            for i in ['callready', 'downstream']:
                for wfrun_id in template[i]:
                    if selected_workflows[wfrun_id]:
                        status = 1
                        break
            for i in ['sequencing', 'alignments', 'analysis']:
                for workflow in template[i]:
                    for wfrun_id in template[i][workflow]:
                        if selected_workflows[wfrun_id]:
                            status = 1
                            break
        D[case] = status
    
    return D
    
