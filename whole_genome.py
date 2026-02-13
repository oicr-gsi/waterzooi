# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:42:36 2023

@author: rjovelin
"""

import os
import itertools
import json
import time
import plotly.graph_objects as go
import networkx as nx
from utilities import connect_to_db 
from project import *



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



def get_workflow_counts(case_id, database, workflow_table='Workflows'):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary with the number of files and the amount of data 
    (ie, lane count) for each workflow in a case
        
    Parameters
    ----------
    - case_id (str): Case identifier
    - database (str): Path to the sqlite database
    - workflow_table (str): Name of the table containing the workflow information in database
    '''
    
    conn = connect_to_db(database)
    query = "SELECT DISTINCT {0}.file_count, {0}.lane_count, {0}.wfrun_id FROM {0} WHERE {0}.case_id = ?;".format(workflow_table)
    data = conn.execute(query, (case_id,)).fetchall()
    conn.close()

    counts = {}
    for i in data:
        workflow_id = i['wfrun_id']
        lane_count = i['lane_count']
        file_count = i['file_count']
        counts[workflow_id] = {'file_count': file_count, 'lane_count': lane_count}
            
    return counts

   

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
    data = conn.execute("SELECT case_id, donor_id, template, valid, error, md5 FROM templates WHERE \
                        project_id = ? AND assay = ?;", (project_name,assay)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        case = i['case_id']
        md5sum = i['md5']
        template = json.loads(i['template'])
        valid = int(i['valid'])
        donor = i['donor_id']
        error = i['error']
        
        d = {'md5sum': md5sum, 'template': template, 'valid': valid, 'error': error, 'donor': donor}        
        
        if case in D:
            D[case].append(d)
        else:
            D[case] = [d]
        
    return D    



def get_case_error_message(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with the error messages across all templates for each case
    
    Parameters
    ----------
    - case_data (dict): Dictionary of cases with analysis data corresponding to project and assay
    '''

    D = {}
    
    for case_id in case_data:
        for d in case_data[case_id]:
            error = d['error']
            error = error.split(';')
            if case_id in D:
                D[case_id].extend(error)
            else:
                D[case_id] = error
    for case_id in D:
        D[case_id] = ';'.join(sorted(list(set(D[case_id]))))

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
    for case_id in cases:
        samples = []
        for d in cases[case_id]:
            if d['template']:
                for i in d['template']['Samples']:
                    samples.append(d['template']['Samples'][i]['sample']) 
        D[case_id] = list(set(samples))
        
    return D
    

def organize_analysis_workflows(case_data, parent_to_children):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with workflow ids organized in sections for each template and each case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with analyses for each case
    '''
        
    D = {}
    
    for case_id in case_data:
        for template in case_data[case_id]:
            if template['template']:
                # track sequencing workflows
                seq = []            
                # get sequencing workflows
                sequencing = []
                for d in template['template']['Data']['Sequencing']:
                    sequencing.append([d['workflow_id'], d['workflow_name']])
                    seq.append(d['workflow_id'])
                # track the anchor workflows
                anchors = []
                for i in template['template']['Anchors']:
                    anchors.extend(template['template']['Anchors'][i])
                # get the analysis and callready workflows
                callready, analyses = [], []
                for i in template['template']['Analysis']:
                    for d in template['template']['Analysis'][i]:
                        workflow_id = d['workflow_id']
                        workflow_name = d['workflow_name']
                        if workflow_id in anchors:
                            callready.append([workflow_id, workflow_name])
                        else:
                            analyses.append([workflow_id, workflow_name])
            
                # get alignment workflows
                alignments = []
                for d in template['template']['Data']['Sequencing']:
                    # get the children of the sequencing workflow
                    children = parent_to_children[d['workflow_id']]
                    for workflow_id in children:
                        # keep workflow that are in analyses
                        # but not in anchors
                        # (ie star_call_read <-- fastqs; bmpp <-- lane level bams)
                        if workflow_id not in anchors:
                            for i in analyses:
                                if workflow_id == i[0]:
                                    alignments.append([i[0], i[1]])
                                    break
                            
                # sort lists according to workflow names
                callready.sort(key=lambda x: x[1])
                alignments.sort(key=lambda x: x[1])
                sequencing.sort(key=lambda x: x[1])
                analyses.sort(key=lambda x: x[1])
            
                data = {'callready': callready, 'alignments': alignments,
                        'sequencing': sequencing, 'analyses': analyses}
            
                if case_id in D:
                    D[case_id].append(data)
                else:
                    D[case_id] = [data]

    return D



def count_case_analysis_workflows(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with analysis workflow ids organized for each case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case analysis extracted from the analysis review database
    '''
        
    D = {}
    
    for case_id in case_data:
        if case_id not in D:
            D[case_id] = {'callready': [], 'downstream': [], 'analysis': {}}
        for template in case_data[case_id]:
            if template['template']:
                for sample in template['template']['Anchors']:
                    D[case_id]['callready'].extend(template['template']['Anchors'][sample])
                for workflow in template['template']['Analysis']:
                    for d in template['template']['Analysis'][workflow]:
                        if d['workflow_id'] not in D[case_id]['callready']:
                            D[case_id]['downstream'].append(d['workflow_id'])
                        if workflow not in D[case_id]['analysis']:
                            D[case_id]['analysis'][workflow] = []
                        D[case_id]['analysis'][workflow].append(d['workflow_id'])
            
        D[case_id]['callready'] = len(set(D[case_id]['callready']))
        D[case_id]['downstream'] = len(set(D[case_id]['downstream']))
        for workflow in D[case_id]['analysis']:
            D[case_id]['analysis'][workflow] = len(set(D[case_id]['analysis'][workflow]))
          
    return D


def list_assay_analysis_workflows(workflow_counts):
    '''
    (dict) -> list
        
    Returns a list of all the expected analysis workflows for an essay
    
    Parameters
    ----------
    - (workflow_counts): Dictionary with counts of each analysis workflow for each case with the same assay
    '''
    
    # get all the analysis workflows
    analysis_workflows = []
    for case_id in workflow_counts:
        for workflow in workflow_counts[case_id]['analysis']:
            analysis_workflows.append(workflow)
    analysis_workflows = sorted(list(set(analysis_workflows)))    
    
    return analysis_workflows


def list_case_analysis_status(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with anlysis status of each case. Analysis status is
    True if at least 1 template for a case is complete
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case analysis extracted from the analysis review database
    '''

    D = {}
    for case_id in case_data:
        for d in case_data[case_id]:
            if case_id in D:
                D[case_id].append(d['valid'])
            else:
                D[case_id] = [d['valid']]

    for case_id in D:
        D[case_id] = bool(sum(D[case_id]))

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
                        'fileimport', 'import_fastq', 'bwamem', 'bwamem2', 'star_lane_level') AND Workflow_Inputs.case_id = ?;", (case,)).fetchall()
    conn.close()

    D = {}
    for i in data:
        limskey = i['limskey']
        library = i['library']
        wfrun_id = i['wfrun_id']
        
        assert wfrun_id not in D
        D[wfrun_id] = {'limskey': limskey, 'library': library}
       
    return D



def most_recent_analysis_workflow(case_data, creation_dates):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with a list of the most recent workflow for each analysis template of each case
       
    Parameters
    ----------
    - case_data (list): Dictionary with template information for each case
    - creation_dates (dict): Dictionary with creation dates of each workflow
    '''
        
    D = {}
    
    for case_id in case_data:
        most_recent = []
        for template in case_data[case_id]:
            if template['template']:
                L = []
                for i in ['Analysis', 'Data']:
                    for j in template['template'][i]:
                        for d in template['template'][i][j]:
                            workflow_id = d['workflow_id']
                            L.append(creation_dates[workflow_id])
                L.sort()
                try:
                    date = time.strftime('%Y-%m-%d', time.localtime(int(L[-1])))
                except:
                    date = 'NA'
                most_recent.append(date)
            else:
                most_recent = 'NA'
        D[case_id] = most_recent
        
    return D



def get_analysis_workflow_name(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with the worklow names of all workflow ids 
    in each template for each case
    
    Parameters
    ----------
    - case_data (list): Dictionary with template information for each case
    '''
    
    D = {}
    
    for case_id in case_data:
        for template in case_data[case_id]:
            if template['template']:
                for i in ['Analysis', 'Data']:
                    for j in template['template'][i]:
                        for d in template['template'][i][j]:
                            workflow_id = d['workflow_id']
                            workflow_name = d['workflow_name']
                            D[workflow_id] = workflow_name

    return D    


def get_case_workflow_samples(database, case_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary of workflow ids and list of corresponding samples for a single case
    
    Parameters
    ----------
    - database (str): Path to the database
    - case_id (str): Case of interest
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT Workflow_Inputs.wfrun_id, Libraries.donor_id, Libraries.tissue_type, \
                        Libraries.tissue_origin, Libraries.library_type, Libraries.group_id FROM \
                        Workflow_Inputs JOIN Libraries WHERE Workflow_Inputs.library=Libraries.library \
                        AND Libraries.case_id = ?;", (case_id,)).fetchall()     
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
    (dict) -> dict
    
    Returns a dictionary with a list of list of mising workflows for each template of each case
    
    Parameters
    ----------
    - case_data (list): Dictionary with analysis data organized by case 
    '''
    
    D = {}
    
    for case_id in case_data:
        missing = []
        for template in case_data[case_id]:
            if template['template']:
                L = []
                for workflow in template['template']['Analysis']:
                    if len(template['template']['Analysis'][workflow]) == 0:
                        L.append(workflow)
                L = list(set(L))                    
                missing.append(L) 
        D[case_id] = missing
    
    return D
    
    
 
    
    
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



def get_pipeline_deliverables(selection):
    '''
    (str) -> dict
    
    Returns a dictionary with the file extensions or file endings for each workflow
    for which output files are released as part of the standard WGS package 
    or as part of the MOH deliverables
    
    Parameters
    ----------
     None
    '''
    
    if selection == 'standard':
        deliverables = {'bammergepreprocessing': ['.bai', '.bam'],
                        'star': ['.bai', '.bam'],
                        'rsem': ['.genes.results', '.isoforms.results', '.transcript.bam'],
                        'hrdetect': ['.signatures.json'],
                        'mavis': ['.tab', '.zip'],
                        'starfusion': ['.tsv'],
                        'arriba': ['.tsv', '.fusions.pdf'],
                        'msisensor': ['.msi', '.msi_germline', '.msi_somatic', '.msi.booted'],
                        'purple': ['.purple_alternates.zip',
                                   'purple.cnv.gene.tsv',
                                   '.purple.cnv.somatic.tsv',
                                   '.solPrimary.purple.zip',
                                   '.purple.purity.range.tsv',
                                   '.purple.purity.tsv',       
                                   '.purple.qc',
                                   '.purple.segment.tsv',
                                   '.purple.somatic.vcf.gz',
                                   '.purple.somatic.vcf.gz.tbi',
                                   '.purple.purity.qc'],
                        'gridds': ['.purple.sv.vcf.gz',
                                   '.purple.sv.vcf.gz.tbi'],
                        'varianteffectpredictor': ['.mutect2.filtered.vep.vcf.gz',
                                                   '.mutect2.filtered.vep.vcf.gz.tbi',
                                                   '.mutect2.filtered.maf.gz'],
                        'delly': ['.somatic_filtered.delly.merged.vcf.gz',
                                  '.somatic_filtered.delly.merged.vcf.gz.tbi'],
                        'sequenza': ['.zip', 'summary.pdf', 'alternative_solutions.json'],
                        'bcl2fastq': ['.fastq.gz'],
                        'fileimportforanalysis': ['.fastq.gz'],
                        'fileimport': ['.fastq.gz'],
                        'import_fastq': ['.fastq.gz']}
    
    elif selection == 'MOH_pipeline':
        deliverables = {'bammergepreprocessing': ['.bai', '.bam'],
                        'star_call_ready': ['.bai', '.bam'],
                        'rsem': ['.genes.results', '.isoforms.results'],
                        'hrdetect': ['.signatures.json'],
                        'mavis': ['.tab', '.zip'],
                        'starfusion': ['.tsv'], 
                        'arriba': ['.tsv'],
                        'msisensor': ['.msi',
                                      '.msi_germline',
                                      '.msi_somatic',
                                      '.msi.booted'],
                        'purple': ['.purple_alternates.zip',
                                   'purple.cnv.gene.tsv',
                                   '.purple.cnv.somatic.tsv',
                                   '.solPrimary.purple.zip',
                                   '.purple.purity.range.tsv',
                                   '.purple.purity.tsv',       
                                   '.purple.qc',
                                   '.purple.segment.tsv',
                                   '.purple.somatic.vcf.gz',
                                   '.purple.somatic.vcf.gz.tbi',
                                   '.purple.purity.qc'],
                        'gridds': ['.purple.sv.vcf.gz',
                                   '.purple.sv.vcf.gz.tbi'],
                        'varianteffectpredictor': ['.mutect2.filtered.vep.vcf.gz',
                                                   '.mutect2.filtered.vep.vcf.gz.tbi',
                                                   '.mutect2.filtered.maf.gz'],
                        'delly': ['.somatic_filtered.delly.merged.vcf.gz',
                                  '.somatic_filtered.delly.merged.vcf.gz.tbi'],
                        'varscan': ['somatic.snp.vcf',
                                    'somatic.snp',
                                    'somatic.copynumber.filtered'],
                        'sequenza': ['.zip', 'alternative_solutions.json'],
                        'haplotypecaller': ['.g.vcf.gz',
                                            '.g.vcf.gz.tbi'],
                        'bcl2fastq': ['.fastq.gz'],
                        'fileimportforanalysis': ['.fastq.gz'],
                        'fileimport': ['.fastq.gz'],
                        'import_fastq': ['.fastq.gz']}
        
    else:
        deliverables = {}
    
    return deliverables



def get_cbioportal_deliverables():
    '''
    (None) -> dict
    
    Returns a dictionary with the file extensions or file endings for each workflow
    for which output files are released to cbioportal 
    
    Parameters
    ----------
     None
    '''

    deliverables = {'varianteffectpredictor': ['.maf.gz'],
                    'sequenza': ['.zip'],
                    'mavis': ['mavis_summary.tab'],
                    'rsem': ['.genes.results'],
                    'purple': ['.purple.purity.tsv',
                               '.purple.cnv.somatic.tsv']}
                    
    return deliverables



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
    - deliverables (None | dict): None or dictionary with file extensions of standard deliverables or MOH deliverables
    '''
        
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
   
    D = {}
    
    for case_id in case_data:
        for template in case_data[case_id]:
            for i in ['Analysis', 'Data']:
                for j in template['template'][i]:
                    for d in template['template'][i][j]:
                        workflow_id = d['workflow_id']
                        workflow_name = d['workflow_name']
                        # check that workflow is selected
                        if workflow_id in selected_workflows and selected_workflows[workflow_id]:
                            # get the workflow output files
                            outputfiles = workflow_outputfiles[workflow_id]
        
                            # check that only workflows in standard eliverables are used
                            if deliverables:
                                key = workflow_name.split('_')[0].lower()
                                if key in deliverables:
                                    # map all file endings of deliverables with files
                                    groups = list(itertools.product(outputfiles, deliverables[key]))
                                    # determine which files are part of the deliverables
                                    F = list(map(G, groups))
                                    L = [groups[k][0] for k in range(len(F)) if F[k]]
                                    if L:
                                        if case_id not in D:
                                            D[case_id] = {}
                                        if workflow_name not in D[case_id]:
                                            D[case_id][workflow_name] = {}
                                        if workflow_id not in D[case_id][workflow_name]:
                                            D[case_id][workflow_name][workflow_id] = L
                                        else:
                                            D[case_id][workflow_name][workflow_id].extend(L)
                            else:
                                if case_id not in D:
                                    D[case_id] = {}
                                if workflow_name not in D[case_id]:
                                    D[case_id][workflow_name] = {}
                                if workflow_id not in D[case_id][workflow_name]:
                                    D[case_id][workflow_name][workflow_id] = outputfiles
                                else:
                                    D[case_id][workflow_name][workflow_id].extend(outputfiles)
                                    
                            if workflow_name in D[case_id] and workflow_id in D[case_id][workflow_name]:
                                D[case_id][workflow_name][workflow_id] = sorted(list(set(D[case_id][workflow_name][workflow_id])))    
    
    
    return D



def create_case_analysis_json(case_data, selected_workflows, workflow_outputfiles, selection):
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
                       Values: standard or all  or MOH_pipeline
    '''
    
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
            
    # get the deliverables
    deliverables = get_pipeline_deliverables(selection)
    
    D = {}
            
    for case_id in case_data:
        for template in case_data[case_id]:
            for i in ['Analysis', 'Data']:
                for j in template['template'][i]:
                    for d in template['template'][i][j]:
                        workflow_id = d['workflow_id']
                        workflow_name = d['workflow_name']
                        # check that workflow is selected
                        if workflow_id in selected_workflows and selected_workflows[workflow_id]:
                            # get the workflow output files
                            outputfiles = workflow_outputfiles[workflow_id]
            
                            # check that only workflows in standard eliverables are used
                            if deliverables:
                                key = workflow_name.split('_')[0].lower()
                                if key in deliverables:
                                    # map all file endings of deliverables with files
                                    groups = list(itertools.product(outputfiles, deliverables[key]))
                                    # determine which files are part of the deliverables
                                    F = list(map(G, groups))
                                    L = [groups[k][0] for k in range(len(F)) if F[k]]
                                    if L:
                                        if case_id not in D:
                                            D[case_id] = {}
                                        if workflow_name not in D[case_id]:
                                            D[case_id][workflow_name] = {}
                                        if workflow_id not in D[case_id][workflow_name]:
                                            D[case_id][workflow_name][workflow_id] = L
                                        else:
                                            D[case_id][workflow_name][workflow_id].extend(L)
                            else:
                                if case_id not in D:
                                    D[case_id] = {}
                                if workflow_name not in D[case_id]:
                                    D[case_id][workflow_name] = {}
                                if workflow_id not in D[case_id][workflow_name]:
                                    D[case_id][workflow_name][workflow_id] = outputfiles
                                else:
                                    D[case_id][workflow_name][workflow_id].extend(outputfiles)
                            
                            if case_id in D:
                                if workflow_name in D[case_id] and workflow_id in D[case_id][workflow_name]:
                                    D[case_id][workflow_name][workflow_id] = sorted(list(set(D[case_id][workflow_name][workflow_id])))  
                       
    return D


def collect_tumor_sample_template(template):
    '''
    (dict) -> str
    
    Returns the tumor sample from the template
          
    Parameters
    ----------
    - template (dict): Dictionary with case template storing analysis data
    '''

    samples = {}
    for i in template['template']['Samples']:
        if template['template']['Samples'][i]['negate_tissue_type']:
            library_type = template['template']['Samples'][i]['library_type']
            sample = template['template']['Samples'][i]['sample']
            samples[library_type] = sample
    
    if len(samples) == 1:
        sample = list(samples.values())[0]
    else:
        assert 'WG' in samples
        sample = samples['WG']

    return sample    
 
    
    
 
def create_cbioportal_json(case_data, selected_workflows, workflow_outputfiles, segmentation):
    '''
    (dict, dict, dict, str)
    
    Returns a dictionary with information required for cbioportal upload
    for donors in a projct 
    
    Parameters
    ----------
    - case_data (dict): Dictionary with analysis templates for cases in a project
    - selected_workflows (dict): Dictionary with selected status of each workflow in project
    - workflow_outputfiles (dict): Dictionary with outputfiles for each workflow run
    - segmentation (str): Indicates if segmentation data comes from the sequenza or purple workflow.
                          Valid values: sequenza, purple
    '''
    
    # create a lambda to evaluate the deliverable files
    # x is a pair of (file, file_ending)
    G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
    
    cbioportal_workflows = ['varianteffectpredictor', 'rsem', 'mavis']
    # add the segmentation workflow
    cbioportal_workflows.insert(1, segmentation)
    deliverables = get_cbioportal_deliverables()    

    D = {}
    
    for case_id in case_data:
        for template in case_data[case_id]:
            donor = template['donor']
            try:
                sample = collect_tumor_sample_template(template)
            except:
                sample = ''
            if sample:
                for workflow in template['template']['Analysis']:
                    # get the generic workflow name                      
                    key = workflow.split('_')[0].lower()
                    if key in cbioportal_workflows:
                        if template['template']['Analysis'][workflow]:
                            for d in template['template']['Analysis'][workflow]:
                                workflow_id = d['workflow_id']
                                # check that workflow is selected
                                if workflow_id in selected_workflows and selected_workflows[workflow_id]:
                                    # get the workflow output files
                                    outputfiles = workflow_outputfiles[workflow_id]
                                    # map all file endings of deliverables with files
                                    groups = list(itertools.product(outputfiles, deliverables[key]))
                                    # determine which files are part of the deliverables
                                    F = list(map(G, groups))
                                    L = [groups[k][0] for k in range(len(F)) if F[k]]
                                    if L:
                                        if case_id not in D:
                                            D[case_id] = {}
                                        if donor not in D[case_id]:
                                            D[case_id][donor] = {}
                                        if sample not in D[case_id][donor]:
                                            D[case_id][donor][sample] = {}
                                        if key == 'purple':
                                            for i in L:
                                                if 'cnv' in i:
                                                    cnvfile = i
                                                elif 'purity' in i:
                                                    purityfile = i
                                            assert cnvfile and purityfile            
                                            D[case_id][donor][sample][workflow] = {'cnv': cnvfile, 'purity': purityfile}
                                        else:
                                            D[case_id][donor][sample][workflow] = L[0]

    return D
    
    
 

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
    
    for case_id in case_data:
        status = 0
        for template in case_data[case_id]:
            if template['template']:
                for i in ['Analysis', 'Data']:
                    for j in template['template'][i]:
                        for d in template['template'][i][j]:
                            if d['workflow_id'] in selected_workflows:
                                if selected_workflows[d['workflow_id']]:
                                    status = 1
                                    break
                for i in template['template']['Anchors']:
                    for j in template['template']['Anchors'][i]:
                        if j in selected_workflows:
                            if selected_workflows[j]:
                                status = 1
                                break
           
        D[case_id] = status
    
    return D



def identify_deliverables(project_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with boolean identicating if analysis pipeline and 
    cbioportal data should be released
       
    Parameters
    ----------
    - project_info (dict): Dictionary with project information
    '''
    
    deliverables = project_info['deliverables'].split(',')
        
    D ={'pipeline': False, 'cbioportal': False}
    for i in deliverables:
        if 'pipeline' in i.lower():
            D['pipeline'] = True
        if 'cbioportal' in i.lower():
            D['cbioportal'] = True
    
    return D



def get_workflow_names(database, case_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary of mapping each workflow id of a case to its name
        
    Parameters
    ----------
    - database (str): Path to the database
    - case_id (str): Case of interest
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT wfrun_id, wf FROM Workflows WHERE case_id = ?;", (case_id,)).fetchall()     
    conn.close()
    
    D = {}
    for i in data:
        workflow_id = i['wfrun_id']
        workflow_name = i['wf']
        D[workflow_id] = workflow_name
    
    return  D


def list_template_workflows(template):
   '''
   (dict) -> list
    
   Returns all the workflow identifiers of a single template of a case
    
   Parameters
   ----------
   - template (dict): Dictionary with case template storing analysis data 
   '''

   L = []

   for i in ['Analysis', 'Data']:
       if template['template']:
           for j in template['template'][i]:
               for d in template['template'][i][j]:
                   workflow_id = d['workflow_id']
                   L.append(workflow_id)
   L = list(set(L))
   
   return L
               
    
    
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
        line=dict(width=0.5, color='#888'),
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
        line_width=2))
    
    
    # color the nodes based on the number of connection
    node_adjacencies = [len(list(G.neighbors(node))) for node in G.nodes()]
    node_trace.marker.color = node_adjacencies
    
    # to change the size of the marker based on the number of connection
    #node_trace.marker.size = node_adjacencies
    
    # label the nodes with the workflow names
    node_text = [str(node) for node in G.nodes()]
    node_text = [workflow_names[i] for i in node_text]
    node_trace.text = node_text
    
    # generate figure
    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title='Workflow connections',
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    return fig

