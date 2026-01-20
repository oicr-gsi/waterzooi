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

   
# def get_amount_data(case_id, database, workflow_table='Workflows'):
#     '''
#     (str, str, str) -> dict
    
#     Returns a dictionary with the amount of data (ie, lane count) for each workflow in a case
        
#     Parameters
#     ----------
#     - case_id (str): Case identifier
#     - database (str): Path to the sqlite database
#     - workflow_table (str): Name of the table containing the workflow information in database
#     '''
    
#     conn = connect_to_db(database)
#     query = "SELECT DISTINCT {0}.lane_count, {0}.wfrun_id FROM {0} WHERE {0}.case_id = ?;".format(workflow_table)
#     data = conn.execute(query, (case_id,)).fetchall()
#     conn.close()

#     counts = {}
#     for i in data:
#         counts[i['wfrun_id']] = i['lane_count']
    
#     return counts


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
    for case in cases:
        samples = []
        for d in cases[case]:
            if d['template']:
                for i in d['template']['Samples']:
                    samples.append(d['template']['Samples'][i]['sample']) 
        D[case] = list(set(samples))
        
    return D
    

# def get_case_analysis_workflows(cases):
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
#             D[case] = []
#         for d in cases[case]:
#             callready = []
#             if 'Anchors' in d['template']:
#                 for i in d['template']['Anchors']:
#                     callready.append(d['template']['Anchors'][i]['workflows'])
#             callready = list(set(callready))
#             downstream = []
#             if 'Analysis' in d['template']:
#                 for i in d['template']['Analysis']:
#                     for j in d['template']['Analysis'][i]:
#                         downstream.append(j['workflow_id'])
#             downstream = list(set(downstream))
#             sequencing = {}
#             if 'Data' in d['template']:
#                 for i in d['template']['Data']['Sequencing']['workflows']:
#                     workflow_name = i['workflow_name']
#                     workflow_id = i['workflow_id']
#                     if workflow_name in sequencing:
#                         sequencing[workflow_name].append(workflow_id)
#                     else:
#                         sequencing[workflow_name] = [workflow_id]
#                 sequencing[workflow_name] = list(set(sequencing[workflow_name]))
#             analysis = {}
#             if 'Analysis' in d['template']:
#                 for i in d['template']['Analysis']:
#                     if i not in analysis:
#                         analysis[i] = []
#                     for j in d['template']['Analysis'][i]:
#                         analysis[i].append(j['workflow_id'])
#                 analysis[i] = list(set(analysis[i]))   
#             alignments = {}
#             if 'Data' in d['template']:
#                 for i in d['template']['Data']:
#                     if i != 'Sequencing':
#                         for k in d['template']['Data'][i]['workflows']:
#                             wfname = k['workflow_name']
#                             workflow_id = k['workflow_id']   
#                             if wfname not in sequencing:
#                                 if wfname not in alignments:
#                                     alignments[wfname] = []
#                                 alignments[wfname].append(workflow_id)
#                                 alignments[wfname] = list(set(alignments[wfname]))
                        
#             D[case].append({'callready': callready, 'downstream': downstream,
#                             'sequencing': sequencing, 'analysis': analysis,
#                             'alignments': alignments})
        
#     return D



def organize_analysis_workflows(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with workflow ids organized in sections for each template and each case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with analyses for each case
    '''
        
    D = {}
    
    for case_id in case_data:
        for template in case_data[case_id]:
            # track sequencing workflows
            seq = []            
            # get sequencing workflows
            sequencing = []
            for d in template['template']['Data']['Sequencing']:
                sequencing.append([d['workflow_id'], d['workflow_name']])
                seq.append(d['workflow_id'])
            # get alignment workflows
            alignments = []
            for i in template['template']['Data']:
                if i != 'Sequencing':
                    for d in template['template']['Data'][i]:
                        if d['workflow_id'] not in seq:
                            alignments.append([d['workflow_id'], d['workflow_name']])    
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
    
    # {"Samples": {"TumourWT": {"library_type": "WT",
    #                           "tissue_type": "R",
    #                           "negate_tissue_type": true,
    #                           "sample": "PANX_1779_Pa_P 101-060-017_LCM_WT"},
    #              "TumourWG": {"library_type": "WG", "tissue_type": "R",
    #                           "negate_tissue_type": true,
    #                           "sample": "PANX_1779_Pa_P 101-060-017_LCM_WG"},
                 
    #              "NormalWG": {"library_type": "WG", "tissue_type": "R",
    #                           "negate_tissue_type": false,
    #                           "sample": "PANX_1779_Ly_R 101-060-017_BuffyCoat_WG"}},
    #  "Data": {"Sequencing": [
    #      {"workflow_id": "vidarr:clinical/run/30f3f85fb20fe586aa30939647762aad947caf2d30d388207e4a765cad9ac6e6",
    #       "workflow_name": "bcl2fastq",
    #       "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG":
    #                    {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}],
    #           "inputs": ["NA"]},
    #          {"workflow_id": "vidarr:clinical/run/084876172db77249210ed8715466a12bb4ad2039c3028e67e43c0d8c995f2882",
    #           "workflow_name": "bcl2fastq",
    #           "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG":
    #                        {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}],
    #               "inputs": ["NA"]},
    #              {"workflow_id": "vidarr:clinical/run/d7392df73793867d2e0c18052854f64684461170d4800dd71894381ff5f80a69", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/a55acb4f73ae1ee0da34a39bb168fedb111ffb6bf6552a390f3c05b2ea4b6856", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/51146126538f9c98dfef7cb438106ae10737bd99af53f3d9edf51fcbf13232a8", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/ff0ef41835ee9d891530a9e40f7d855119c44d4e00fe892cf6a93d6a938f258b", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/195ed99ab1377ed439fc019f751cb1d19204b66a208033ce1278fa5e46a78ae3", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/2db568631e17ea96d9514e21e8ac25f52fda2288dd2128046d621b4fb5d76637", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/2bc701dc18e47bd12a8bfe5fe75fb3ca15570cd04221e8cdb187eff70dd8fe8e", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/d23b9536b8c3b52deeb05fd5f00211fab89afc4010cf8ea7df3562f73cd878ee", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/4b78dc4a7622b74283d6ea9b977863d462c660dcf6d1359c957138f7f07586fe", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/76479108cdda0584eba4f02f1cba55f7302511604d3ac95611f9ff58c8cfb3c0", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/ac6d7117e0f1f2dbb4b5217c535a12e36aabf59c55eb38eb052a4e2a0cba0a58", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}], "TumourWTalign": [{"workflow_id": "vidarr:clinical/run/d7392df73793867d2e0c18052854f64684461170d4800dd71894381ff5f80a69", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/76479108cdda0584eba4f02f1cba55f7302511604d3ac95611f9ff58c8cfb3c0", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}], "TumourWGalign": [{"workflow_id": "vidarr:clinical/run/30f3f85fb20fe586aa30939647762aad947caf2d30d388207e4a765cad9ac6e6", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/084876172db77249210ed8715466a12bb4ad2039c3028e67e43c0d8c995f2882", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/ff0ef41835ee9d891530a9e40f7d855119c44d4e00fe892cf6a93d6a938f258b", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/195ed99ab1377ed439fc019f751cb1d19204b66a208033ce1278fa5e46a78ae3", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/2db568631e17ea96d9514e21e8ac25f52fda2288dd2128046d621b4fb5d76637", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/d23b9536b8c3b52deeb05fd5f00211fab89afc4010cf8ea7df3562f73cd878ee", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/ac6d7117e0f1f2dbb4b5217c535a12e36aabf59c55eb38eb052a4e2a0cba0a58", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["NA"]}], "NormalWGalign": [{"workflow_id": "vidarr:clinical/run/a55acb4f73ae1ee0da34a39bb168fedb111ffb6bf6552a390f3c05b2ea4b6856", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/51146126538f9c98dfef7cb438106ae10737bd99af53f3d9edf51fcbf13232a8", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/2bc701dc18e47bd12a8bfe5fe75fb3ca15570cd04221e8cdb187eff70dd8fe8e", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}, {"workflow_id": "vidarr:clinical/run/4b78dc4a7622b74283d6ea9b977863d462c660dcf6d1359c957138f7f07586fe", "workflow_name": "bcl2fastq", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["NA"]}]},
         
    #      "Analysis": {"mutect2_matched":
    #                   [{"workflow_id": "vidarr:clinical/run/6b4ee54dc04f8c3fc3e481fcbe796da63446223609ca6cde01edec68cea05f9f",
    #                     "workflow_name": "mutect2_matched",
    #                     "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6", "vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c"]}], "variantEffectPredictor_matched": [{"workflow_id": "vidarr:clinical/run/05da2a0d37a1f416ecf25d3c08b63305fce4fa1ec262bd425654708adeec0488", "workflow_name": "variantEffectPredictor_matched", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/6b4ee54dc04f8c3fc3e481fcbe796da63446223609ca6cde01edec68cea05f9f"]}], "delly_matched": [{"workflow_id": "vidarr:clinical/run/744f84138c3542a09f8bf3565825778e8aa97e764818fb836c9f6eacd21f47ee", "workflow_name": "delly_matched", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6", "vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c"]}], "gridss": [{"workflow_id": "vidarr:clinical/run/7ef5fc4660a363e54639716b1332ed0ba7faee2f67ab402be4f81a8ef04f89e6", "workflow_name": "gridss", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6", "vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c"]}], "purple": [{"workflow_id": "vidarr:clinical/run/c5759642f1d1565722272d09807f89388119f0767cafaae294c1b0f5f433b333", "workflow_name": "purple", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6", "vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c", "vidarr:clinical/run/6b4ee54dc04f8c3fc3e481fcbe796da63446223609ca6cde01edec68cea05f9f", "vidarr:clinical/run/7ef5fc4660a363e54639716b1332ed0ba7faee2f67ab402be4f81a8ef04f89e6"]}], "bamMergePreprocessing_by_sample": [{"workflow_id": "vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6", "workflow_name": "bamMergePreprocessing_by_sample", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/df4859eb6fc3df248ba6c5dc96e5afd5126a64d478f91eae4d21480d96009c11", "vidarr:clinical/run/5cd133979e5b58540ab9276430ca5f23ce6e5416e177dd16ebd1cbaf966338eb", "vidarr:clinical/run/83ea2a6d548ac21d75339bc914dab72102c926570f456699dfc04391785449fd", "vidarr:clinical/run/c05529dcdb18f287128194286661f9578e47e5328a45bc00f9917b37cf7ccf5b", "vidarr:clinical/run/4c99bccb7dbe70524360515034978e356c771dc268b0112e4ed629655e2d887e", "vidarr:clinical/run/e0e5693f225c659a05643be53882f0a80e1da4bd87c84503d24e7e729d0b7562", "vidarr:clinical/run/a8a5392caf65f120120d6dd593e27d1bb548cc195866953f2da090a853434e77"]}, {"workflow_id": "vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c", "workflow_name": "bamMergePreprocessing_by_sample", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["vidarr:clinical/run/7282751b1816edb9b63e993d49869d59a743d797df25fb3cbb10e62cfe7571d2", "vidarr:clinical/run/2ad0644177bff4219c2e2dd54be33c6f768a75be3c09a9f40cf29c6b9e210471", "vidarr:clinical/run/71d572925ba9568bef3169ee70be2c08bea14f4d34a4d7a311868a48964478af", "vidarr:clinical/run/14a404d908ccacf7d2930fb61f01f8c0294e92832c9d398f83cbc3f10e1f9aba"]}], "bwaMem": [{"workflow_id": "vidarr:clinical/run/2ad0644177bff4219c2e2dd54be33c6f768a75be3c09a9f40cf29c6b9e210471", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["vidarr:clinical/run/2bc701dc18e47bd12a8bfe5fe75fb3ca15570cd04221e8cdb187eff70dd8fe8e"]}, {"workflow_id": "vidarr:clinical/run/e0e5693f225c659a05643be53882f0a80e1da4bd87c84503d24e7e729d0b7562", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/ff0ef41835ee9d891530a9e40f7d855119c44d4e00fe892cf6a93d6a938f258b"]}, {"workflow_id": "vidarr:clinical/run/7282751b1816edb9b63e993d49869d59a743d797df25fb3cbb10e62cfe7571d2", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["vidarr:clinical/run/a55acb4f73ae1ee0da34a39bb168fedb111ffb6bf6552a390f3c05b2ea4b6856"]}, {"workflow_id": "vidarr:clinical/run/14a404d908ccacf7d2930fb61f01f8c0294e92832c9d398f83cbc3f10e1f9aba", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["vidarr:clinical/run/51146126538f9c98dfef7cb438106ae10737bd99af53f3d9edf51fcbf13232a8"]}, {"workflow_id": "vidarr:clinical/run/4c99bccb7dbe70524360515034978e356c771dc268b0112e4ed629655e2d887e", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/d23b9536b8c3b52deeb05fd5f00211fab89afc4010cf8ea7df3562f73cd878ee"]}, {"workflow_id": "vidarr:clinical/run/c05529dcdb18f287128194286661f9578e47e5328a45bc00f9917b37cf7ccf5b", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/2db568631e17ea96d9514e21e8ac25f52fda2288dd2128046d621b4fb5d76637"]}, {"workflow_id": "vidarr:clinical/run/71d572925ba9568bef3169ee70be2c08bea14f4d34a4d7a311868a48964478af", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["vidarr:clinical/run/4b78dc4a7622b74283d6ea9b977863d462c660dcf6d1359c957138f7f07586fe"]}, {"workflow_id": "vidarr:clinical/run/83ea2a6d548ac21d75339bc914dab72102c926570f456699dfc04391785449fd", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/195ed99ab1377ed439fc019f751cb1d19204b66a208033ce1278fa5e46a78ae3"]}, {"workflow_id": "vidarr:clinical/run/df4859eb6fc3df248ba6c5dc96e5afd5126a64d478f91eae4d21480d96009c11", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/084876172db77249210ed8715466a12bb4ad2039c3028e67e43c0d8c995f2882"]}, {"workflow_id": "vidarr:clinical/run/a8a5392caf65f120120d6dd593e27d1bb548cc195866953f2da090a853434e77", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/ac6d7117e0f1f2dbb4b5217c535a12e36aabf59c55eb38eb052a4e2a0cba0a58"]}, {"workflow_id": "vidarr:clinical/run/5cd133979e5b58540ab9276430ca5f23ce6e5416e177dd16ebd1cbaf966338eb", "workflow_name": "bwaMem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/30f3f85fb20fe586aa30939647762aad947caf2d30d388207e4a765cad9ac6e6"]}], "rsem": [{"workflow_id": "vidarr:clinical/run/83c406fc6b2e0b6696d2138b203bff7cfb71699660146db4c750f2e5c63ab7fe", "workflow_name": "rsem", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/5c6335a63df24f423ee8f83d029bb7506ca7bcd5b7f0bc01275e986ae12b1e8a"]}], "star_call_ready": [{"workflow_id": "vidarr:clinical/run/5c6335a63df24f423ee8f83d029bb7506ca7bcd5b7f0bc01275e986ae12b1e8a", "workflow_name": "star_call_ready", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/76479108cdda0584eba4f02f1cba55f7302511604d3ac95611f9ff58c8cfb3c0", "vidarr:clinical/run/d7392df73793867d2e0c18052854f64684461170d4800dd71894381ff5f80a69"]}], "starfusion": [{"workflow_id": "vidarr:clinical/run/e95dcc79de60179f3c1b07f4035588f7ae670c505727af78feb02377af8c7432", "workflow_name": "starfusion", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/5c6335a63df24f423ee8f83d029bb7506ca7bcd5b7f0bc01275e986ae12b1e8a"]}], "arriba": [{"workflow_id": "vidarr:clinical/run/dcc7b431e83f4483a10184ad84137136ad32ebcbb30486a637a80a526b0652d6", "workflow_name": "arriba", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/5c6335a63df24f423ee8f83d029bb7506ca7bcd5b7f0bc01275e986ae12b1e8a"]}], "msisensor": [{"workflow_id": "vidarr:clinical/run/8f82920f937f42b8e0782765a45c7867ebea05d749ac8e49fa2d7f6b1f1c3a23", "workflow_name": "msisensor", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6", "vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c"]}], "star_lane_level": [{"workflow_id": "vidarr:clinical/run/54cbe0e438e2731fdf66e121ba274a8a33ea9608b6444e91bbf72802cbd9e21f", "workflow_name": "star_lane_level", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/d7392df73793867d2e0c18052854f64684461170d4800dd71894381ff5f80a69"]}, {"workflow_id": "vidarr:clinical/run/d021c3ac7b2fa0f7ae95f7b8f81528369fa54edd0dba3c5c58bd717c7d98e10f", "workflow_name": "star_lane_level", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/76479108cdda0584eba4f02f1cba55f7302511604d3ac95611f9ff58c8cfb3c0"]}], "mavis": [{"workflow_id": "vidarr:clinical/run/751861f203c7ca0cdf7dd12a6d8e4690bd4716af0beaf238a1e5521116e2c441", "workflow_name": "mavis", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WT": {"library_type": "WT", "tissue_type": "R", "negate_tissue_type": true}}, {"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/5c6335a63df24f423ee8f83d029bb7506ca7bcd5b7f0bc01275e986ae12b1e8a", "vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6", "vidarr:clinical/run/744f84138c3542a09f8bf3565825778e8aa97e764818fb836c9f6eacd21f47ee", "vidarr:clinical/run/dcc7b431e83f4483a10184ad84137136ad32ebcbb30486a637a80a526b0652d6", "vidarr:clinical/run/e95dcc79de60179f3c1b07f4035588f7ae670c505727af78feb02377af8c7432"]}], "hrDetect": [{"workflow_id": "vidarr:clinical/run/ef569288ce7b0b1d2b5343608f5580f3c0ce7b72b183ba8a71118cedb630f2e2", "workflow_name": "hrDetect", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}, {"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/6b4ee54dc04f8c3fc3e481fcbe796da63446223609ca6cde01edec68cea05f9f", "vidarr:clinical/run/c5759642f1d1565722272d09807f89388119f0767cafaae294c1b0f5f433b333"]}], "haplotypeCaller": [{"workflow_id": "vidarr:clinical/run/e91c429b0ee42bc8589155f8b447d2ae4f6a37aa40b87995d3845fc402b077f8", "workflow_name": "haplotypeCaller", "samples": [{"PANX_1779_Pa_P 101-060-017_LCM_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": true}}], "inputs": ["vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6"]}, {"workflow_id": "vidarr:clinical/run/f4567203383a504b23431a273faab11a86f9aa1d377ba0a0a0c56536803b6d1b", "workflow_name": "haplotypeCaller", "samples": [{"PANX_1779_Ly_R 101-060-017_BuffyCoat_WG": {"library_type": "WG", "tissue_type": "R", "negate_tissue_type": false}}], "inputs": ["vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c"]}]}, "Anchors": {"TumourWG": ["vidarr:clinical/run/b3d3f31abfa9a91ad43b5101909043de8926b99b92b5e89c4bd89fa6d80586f6"], "TumourWT": ["vidarr:clinical/run/5c6335a63df24f423ee8f83d029bb7506ca7bcd5b7f0bc01275e986ae12b1e8a"], "NormalWG": ["vidarr:clinical/run/7e329e4697f49970fef560b677c45682b33315470d1580c7c17a8bccebeac04c"]}}
   
    
        
        
    for case_id in case_data:
        if case_id not in D:
            D[case_id] = {'callready': [], 'downstream': [], 'analysis': {}}
        for template in case_data[case_id]:
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
        D[case_id] = most_recent
        
    return D


# def get_analysis_workflow_name(analysis):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with the name of workflows for each workflow id of
#     an analysis group of a single case
    
#     Parameters
#     ----------
#     - analysis (dict): Dictionary with workflow ids of analysis workflows of a single case
#     '''
    
#     D = {}
    
#     for workflow_name in analysis:
#         for workflow_id in analysis[workflow_name]:
#             assert workflow_id not in D
#             D[workflow_id] = workflow_name
    
#     return D
    
    
# def get_analysis_workflow_name(case_data):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with the worklow names of all workflow ids 
#     in each template for each case
    
#     Parameters
#     ----------
#     - case_data (list): Dictionary with template information for each case
#     '''
    
#     D = {}
    
#     for case_id in case_data:
#         for template in case_data[case_id]:
#             k = {}
#             for i in ['Analysis', 'Data']:
#                 for j in template['template'][i]:
#                     for d in template['template'][i][j]:
#                         workflow_id = d['workflow_id']
#                         workflow_name = d['workflow_name']
#                         k[workflow_id] = workflow_name
#             if case_id not in D:
#                 D[case_id] = []
#             D[case_id].append(k)
        
#     return D    
    

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
            L = []
            for i in ['Analysis', 'Data']:
                for workflow in template['template'][i]:
                    if len(template['template'][i][workflow]) == 0:
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
                    'sequenza': ['.zip', 'summary.pdf', 'alternative_solutions.json'],
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


# def create_analysis_json(case_data, selected_workflows, workflow_outputfiles, deliverables=None):
#     '''
#     (dict, dict, dict, dict, dict, dict, None | dict)
    
#     Returns a dictionary with workflow information for a given block (ie, sample pair)
#     and anchor bmpp parent workflow
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with analysis templates for all cases in a project
#     - selected_workflows (dict): Dictionary with selected status of each workflow in project
#     - workflow_outputfiles (dict): Dictionary with outputfiles for each workflow run
#     - deliverables (None | dict): None or dictionary with file extensions of standard deliverables
#     '''
        
#     # create a lambda to evaluate the deliverable files
#     # x is a pair of (file, file_ending)
#     G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
   
#     D = {}
        
#     for case in case_data:
#         for template in case_data[case]:
#             # make a list of workflows:
#             callready = template['callready']
#             workflows = {}
#             for i in ['sequencing', 'analysis', 'alignments']:
#                 for workflow in template[i]:
#                     if template[i][workflow]:
#                         workflows[workflow] = template[i][workflow]
        
#             # check that analysis workflows are selected
#             # do not include call ready workflows because they can be shared across templates
#             wfs = []
#             for i in workflows.values():
#                 wfs.extend(i)
            
#             analysis = [i for i in wfs if i not in callready]   
#             if any(map(lambda x: x in selected_workflows, analysis)):
#                 for workflow in workflows:
#                     for wfrunid in workflows[workflow]:
#                         # check if workflow is selected
#                         if selected_workflows[wfrunid]:
#                             # get workflow output files
#                             outputfiles = workflow_outputfiles[wfrunid]                        
                                        
#                             # check that only workflows in standard eliverables are used
#                             if deliverables:
#                                 key = workflow.split('_')[0].lower()
#                                 if key in deliverables:
#                                     # map all file endings of deliverables with files
#                                     groups = list(itertools.product(outputfiles, deliverables[key]))
#                                     # determine which files are part of the deliverables
#                                     F = list(map(G, groups))
#                                     L = [groups[k][0] for k in range(len(F)) if F[k]]
#                                     if L:
#                                         if case not in D:
#                                             D[case] = {}
#                                         if workflow in D[case]:
#                                             D[case][workflow].extend(L)
#                                         else:
#                                             D[case][workflow] = L
#                             else:
#                                 if case not in D:
#                                     D[case] = {}
#                                 if workflow in D[case]:
#                                     D[case][workflow].extend(outputfiles)
#                                 else:
#                                     D[case][workflow] = outputfiles
                            
#                             D[case][workflow] = sorted(list(set(D[case][workflow])))  
    
    
#     return D
                    


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





# def create_case_analysis_json(case, case_data, selected_workflows, workflow_outputfiles, selection):
#     '''
#     (str, list, dict, dict, str)
    
#     Returns a dictionary with workflow information for a given block (ie, sample pair)
#     and anchor parent workflow (bmpp or star)
    
#     Parameters
#     ----------
#     - case (str): Case unique identifier
#     - case_data (list): List of analysis templates for a single case in a project
#     - selected_workflows (dict): Dictionary with selected status of each workflow in project
#     - workflow_outputfiles (dict): Dictionary with outputfiles for each workflow run
#     - selection (str): Include files from all selected workflows or files from the standard deliverables
#                        Values: standard or all
#     '''
    
#     # create a lambda to evaluate the deliverable files
#     # x is a pair of (file, file_ending)
#     G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
            
#     # get the deliverables
#     if selection == 'standard':
#         deliverables = get_pipeline_standard_deliverables()
#     elif selection == 'all':
#         deliverables = {}
    
#     D = {}
            
#     for template in case_data:
#         # make a list of workflows:
#         callready = template['callready']
#         workflows = {}
#         for i in ['sequencing', 'analysis', 'alignments']:
#             for workflow in template[i]:
#                 if template[i][workflow]:
#                     workflows[workflow] = template[i][workflow]
            
#         # check that analysis workflows are selected
#         # do not include call ready workflows because they can be shared across templates
#         wfs = []
#         for i in workflows.values():
#             wfs.extend(i)
                
#         analysis = [i for i in wfs if i not in callready]   
#         if any(map(lambda x: x in selected_workflows, analysis)):
#             for workflow in workflows:
#                 for wfrunid in workflows[workflow]:
#                     # check if workflow is selected
#                     if selected_workflows[wfrunid]:
#                         # get workflow output files
#                         outputfiles = workflow_outputfiles[wfrunid]                        
                                            
#                         # check that only workflows in standard eliverables are used
#                         if deliverables:
#                             key = workflow.split('_')[0].lower()
#                             if key in deliverables:
#                                 # map all file endings of deliverables with files
#                                 groups = list(itertools.product(outputfiles, deliverables[key]))
#                                 # determine which files are part of the deliverables
#                                 F = list(map(G, groups))
#                                 L = [groups[k][0] for k in range(len(F)) if F[k]]
#                                 if L:
#                                     if case not in D:
#                                         D[case] = {}
#                                     if workflow in D[case]:
#                                         D[case][workflow].extend(L)
#                                     else:
#                                         D[case][workflow] = L
#                         else:
#                             if case not in D:
#                                 D[case] = {}
#                             if workflow in D[case]:
#                                 D[case][workflow].extend(outputfiles)
#                             else:
#                                 D[case][workflow] = outputfiles
                        
#                         D[case][workflow] = sorted(list(set(D[case][workflow])))  
                       
        
        
#     return D
    


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








    
    
# def create_cbioportal_json(case_data, selected_workflows, workflow_outputfiles, segmentation):
#     '''
#     (dict, dict, dict, str)
    
#     Returns a dictionary with information required for cbioportal upload
#     for donors in a projct 
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with analysis templates for cases in a project
#     - selected_workflows (dict): Dictionary with selected status of each workflow in project
#     - workflow_outputfiles (dict): Dictionary with outputfiles for each workflow run
#     - segmentation (str): Indicates if segmentation data comes from the sequenza or purple workflow.
#                           Valid values: sequenza, purple
#     '''
    
#     # create a lambda to evaluate the deliverable files
#     # x is a pair of (file, file_ending)
#     G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
    
#     cbioportal_workflows = ['varianteffectpredictor', 'rsem', 'mavis']
#     # add the segmentation workflow
#     cbioportal_workflows.insert(1, segmentation)
#     deliverables = get_cbioportal_deliverables()    
    
#     D = {}
    
#     for case in case_data:
#         for template in case_data[case]:
#             donor = template['donor']
#             samples = {}
#             for i in template['template']['Samples']:
#                 if template['template']['Samples'][i]['tissue_type'] != 'R':
#                     library_type = template['template']['Samples'][i]['library_type']
#                     sample = template['template']['Samples'][i]['sample']
#                     samples[library_type] = sample
#             if len(samples) == 1:
#                 samples = list(samples.values())[0]
#             else:
#                 assert 'WG' in samples
#                 sample = samples['WG']
        
#             for workflow in template['template']['Analysis']:
#                 # get the generic workflow name                      
#                 key = workflow.split('_')[0].lower()
#                 if key in cbioportal_workflows:
#                     if template['template']['Analysis'][workflow]:
#                         for d in template['template']['Analysis'][workflow]:
#                             workflow_id = d['workflow_id']
#                             if workflow_id in selected_workflows and selected_workflows[workflow_id]:
#                                 # get workflow output files
#                                 outputfiles = workflow_outputfiles[workflow_id]                        
#                                 # map all file endings of deliverables with files
#                                 groups = list(itertools.product(outputfiles, deliverables[key]))
#                                 # determine which files are part of the deliverables
#                                 F = list(map(G, groups))
#                                 L = [groups[k][0] for k in range(len(F)) if F[k]]
#                                 if L:
#                                     if donor not in D:
#                                         D[donor] = {}
#                                     if sample not in D[donor]:
#                                         D[donor][sample] = {}
#                                     if key == 'purple':
#                                         for i in L:
#                                             if 'cnv' in i:
#                                                 cnvfile = i
#                                             elif 'purity' in i:
#                                                 purityfile = i
#                                         assert cnvfile and purityfile            
#                                         D[donor][sample][workflow] = {'cnv': cnvfile, 'purity': purityfile}
#                                     else:
#                                         D[donor][sample][workflow] = L[0]

#     return D
    
 
    
 
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
            sample = collect_tumor_sample_template(template)
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
                                    if donor not in D:
                                        D[donor] = {}
                                    if sample not in D[donor]:
                                        D[donor][sample] = {}
                                    if key == 'purple':
                                        for i in L:
                                            if 'cnv' in i:
                                                cnvfile = i
                                            elif 'purity' in i:
                                                purityfile = i
                                        assert cnvfile and purityfile            
                                        D[donor][sample][workflow] = {'cnv': cnvfile, 'purity': purityfile}
                                    else:
                                        D[donor][sample][workflow] = L[0]

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




# def get_review_status(case_data, selected_workflows):
#     '''
#     (dict, dict) -> dict
    
#     Returns a dictionary with review status for each case in project
#     A case is considered to be reviewed if any worklow has been selected
    
#     Parameters
#     ----------
#     - case_data (dict): Dictionary with analysis template for each case
#     - selected_workflows (dict): Dictionary with selection status of each workflow in a project
#     '''
    
#     D = {}
    
#     for case in case_data:
#         status = 0
#         for template in case_data[case]:
#             for i in ['callready', 'downstream']:
#                 for wfrun_id in template[i]:
#                     if selected_workflows[wfrun_id]:
#                         status = 1
#                         break
#             for i in ['sequencing', 'alignments', 'analysis']:
#                 for workflow in template[i]:
#                     for wfrun_id in template[i][workflow]:
#                         if selected_workflows[wfrun_id]:
#                             status = 1
#                             break
#         D[case] = status
    
#     return D
    


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


# def list_template_workflows(template):
#     '''
#     (dict) -> list
    
#     Returns all the workflow identifiers of a single template of a case
    
#     Parameters
#     ----------
#     - template (dict): Dictionary with analysis data for a single template of a case
#     '''

#     L = []

#     if 'Data' in template['template']:
#         for i in template['template']['Data']:
#             for d in template['template']['Data'][i]['workflows']:
#                 L.append(d['workflow_id'])
#     if 'Analysis' in template['template']:
#         for workflow in template['template']['Analysis']:
#             for d in template['template']['Analysis'][workflow]:
#                 L.append(d['workflow_id'])
#     L = list(set(L))
    
#     return L
    


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

