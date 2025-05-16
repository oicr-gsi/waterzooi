# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:35:40 2023

@author: rjovelin
"""

import os
import itertools
from utilities import connect_to_db, get_children_workflows, remove_non_analysis_workflows, get_donors

 


# def get_WT_call_ready_cases(project_name, platform, database, library_type='WT'):
#     '''
#     (str, str, str, str) -> dict

#     Returns a dictionary with samples and libraries and bmpp and downstream workflow ids for each case in a project,
#     restricting data to specified platform and library type

#     Parameters
#     ----------
#     - project_name (str): Name of the project
#     - platform (str): Name of sequencing platform.
#                       Accepted values: novaseq, nextseq, hiseq, miseq
#     - database (str): Path to the sqlite database
#     - library_type (str): 2 letters-code indicating the type of library                   
#     '''

#     # get all the samples for project name 
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT DISTINCT Libraries.library, Libraries.case_id, Libraries.project_id, \
#                           Libraries.ext_id, Libraries.group_id, Libraries.library_type, \
#                           Libraries.tissue_type, Libraries.tissue_origin, \
#                           Workflows.wf, Workflows.wfrun_id, Workflow_Inputs.platform \
#                           from Workflow_Inputs JOIN Libraries JOIN Workflows \
#                           WHERE Libraries.project_id = '{0}' AND Workflow_Inputs.project_id = '{0}' \
#                           AND Workflows.project_id = '{0}' AND Workflow_Inputs.wfrun_id = Workflows.wfrun_id \
#                           AND Workflow_Inputs.library = Libraries.library AND Libraries.library_type = '{1}';".format(project_name, library_type)).fetchall()
#     conn.close()

#     cases = {}
#     for i in data:
#         # select star data sequenced on novaseq
#         if platform in i['platform'].lower():
#             if 'star_call_ready' in i['wf'].lower():
#                 if i['case_id'] not in cases:
#                     cases[i['case_id']] = {'project': i['project_id'], 'samples': [], 'libraries': [], 'star': []}
#                 cases[i['case_id']]['star'].append(i['wfrun_id'])
#                 sample = '_'.join([i['case_id'], i['tissue_type'], i['tissue_origin'], i['library_type'], i['group_id']]) 
#                 cases[i['case_id']]['samples'].append(sample)
#                 cases[i['case_id']]['libraries'].append(i['library'])
            
#     # get parent-children workflow relationships
#     parents = get_children_workflows(project_name, database)

#     # find the bmpp downstream workflows
#     for sample in cases:
#         downstream = []
#         for star in cases[sample]['star']:
#             if star in parents:
#                 # get the star downstream workflows
#                 children = parents[star]
#                 # removed any non-analysis workflow
#                 children = remove_non_analysis_workflows(children)
#                 # list all downtream workflows
#                 downstream.extend([i['children_id'] for i in children])
#                 # get the downstream workflows of downstream workflows
#                 # remove non-analysis workflows
#                 for workflow in downstream:
#                     if workflow in parents:
#                         L = remove_non_analysis_workflows(parents[workflow])
#                         downstream.extend([i['children_id'] for i in L])
#         cases[sample]['downstream'] = list(set(downstream)) 
        
#     return cases



# def get_WT_standard_deliverables():
#     '''
#     (None) -> dict
    
#     Returns a dictionary with the file extensions or file endings for each workflow
#     for which output files are released as part of the standard WT package
    
#     Parameters
#     ----------
#      None
#     '''
    
#     deliverables = {'star': ['.bai', '.bam'],
#                     'mavis': ['.zip', '.tab'],
#                     'rsem': ['.genes.results', '.isoforms.results', '.transcript.bam'],
#                     'starfusion': ['.tsv'],
#                     'arriba': ['.tsv']}
       
#     return deliverables


# def create_WT_project_block_json(project_name, database, blocks, block_status, selected_workflows, workflow_names, deliverables=None):
#     '''
#     (str, str, dict, dict, dict, dict, None | dict)
    
#     Returns a dictionary with workflow information for a given block (ie, sample pair)
#     and anchor star parent workflow
    
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
#         for sample in blocks[case]:
#             # check the selection status of the block
#             if block_status[case][sample] not in ['ready', 'review']:
#                 # block already reviewed and workflows selected
#                 anchor_wf = block_status[case][sample]
                                
#                 for workflow in blocks[case][sample][anchor_wf]['workflows']:
                    
#                     workflow = os.path.basename(workflow)
                    
#                     # get workflow name and version
#                     workflow_name = workflow_names[workflow][0]
#                     workflow_version = workflow_names[workflow][1]
                    
#                     # check workflow status
#                     if selected_workflows[workflow]:
                                                
#                         # get workflow output files
#                         # needed to sort outputs by sample pairs or by sample for call-ready workflows
#                         # even if all files are recorded
#                         outputfiles = workflow_outputfiles[workflow]
                        
#                         # check that only workflows in standard WGS deliverables are used
#                         if deliverables:
#                             # keep track of the files to be released                                            
#                             L = []
#                             key = workflow_name.split('_')[0].lower()
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
#                                         sample_id = j.replace('.', ';')
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



# def create_WT_block_json(database, project_name, case, blocks, sample, anchor_workflow, workflow_names, selected_workflows, selection):
#     '''
#     (str, str, dict, str, str, dict, dict, str)
    
#     Returns a dictionary with workflow information for a given block (ie, sample pair)
#     and anchor parent workflow (bmpp or star)
    
#     Parameters
#     ----------
#     - database (str): Path to the sqlite database
#     - project_name (str): Name of project of interest
#     - case (str): Donor identifier 
#     - blocks (dict): Dictionary with block information
#     - sample (str): Tumor sample
#     - anchor_workflow (str): bamMergePreprocessing parent workflow(s) or star_call_ready parent workflow
#     - workflow_names (dict): Dictionary with workflow name and version for each workflow in project
#     - selected_workflows (dict): Dictionary with selected status of each workflow in project
#     - selection (str): Include files from all selected workflows or files from the standard deliverables
#                        Values: standard or all
#     '''
    
#     libraries = map_limskey_to_library(project_name, database, table='Workflow_Inputs')
#     sample_names = map_library_to_sample(project_name, database, table = 'Libraries')
#     donors = map_library_to_case(project_name, database, table = 'Libraries')
#     workflow_outputfiles = get_workflow_output(project_name, database, libraries, sample_names, donors, 'Files')
    
    
#     # create a lambda to evaluate the deliverable files
#     # x is a pair of (file, file_ending)
#     G = lambda x: x[1] in x[0] and x[0][x[0].rindex(x[1]):] == x[1]
        
    
#     # get the deliverables
#     if selection == 'standard':
#         deliverables = get_WT_standard_deliverables()
#     elif selection == 'all':
#         deliverables = {}
    
#     # organize the workflows by block and samples
#     D = {}
#     # get the workflow ids for that block
#     for i in blocks[sample]:
#         if i['anchor_wf'] == anchor_workflow:
#             D[sample] = map(lambda x: x.strip(), i['workflows'].split(';'))
    
#     block_data = {}
#     for sample in D:
#         for workflow_id in D[sample]:
#             # check if workflow is selected
#             if selected_workflows[workflow_id]:
#                 # get workflow name and version
#                 workflow_name = workflow_names[workflow_id][0]
#                 workflow_version = workflow_names[workflow_id][1]
                
#                 # get workflow output files
#                 # needed to sort outputs by sample pairs or by sample for call-ready workflows
#                 # even if all files are recorded
#                 outputfiles = workflow_outputfiles[workflow_id]
                                
#                 # check that only workflows in standard WGS deliverables are used
#                 if deliverables:
#                     key = workflow_name.split('_')[0].lower()
#                     if key in deliverables:
#                         for j in outputfiles:
#                             # list all deliverable files
#                             L = []
#                             # gather all file paths for workflow and sample(s)
#                             files = [i[0] for i in outputfiles[j]]
#                             # map all file endings of deliverables with files
#                             groups = list(itertools.product(files, deliverables[key]))
#                             # determine which files are part of the deliverables
#                             F = list(map(G, groups))
#                             L = [groups[k][0] for k in range(len(F)) if F[k]]
                            
#                             if L:
#                                 sample_id = '.'.join(j.split(';'))
#                                 if case not in block_data:
#                                     block_data[case] = {}
#                                 if sample_id not in block_data[case]:
#                                     block_data[case][sample_id] = {}
#                                 if workflow_name not in block_data[case][sample_id]:
#                                     block_data[case][sample_id][workflow_name] = []
                                
#                                 d = {'workflow_id': workflow_id,
#                                      'workflow_version': workflow_version,
#                                      'files': L}
#                                 if d not in block_data[case][sample_id][workflow_name]:
#                                     block_data[case][sample_id][workflow_name].append(d)
                                    
#                 else:
#                     for j in outputfiles:
#                         sample_id = '.'.join(j.split(';'))
#                         d =  {'workflow_id': workflow_id, 
#                               'workflow_version': workflow_version,
#                               'files': [i[0] for i in outputfiles[j]]}
#                         if case not in block_data:
#                             block_data[case] = {}
#                         if sample_id not in block_data[case]:
#                             block_data[case][sample_id] = {}
#                         if workflow_name not in block_data[case][sample_id]:
#                             block_data[case][sample_id][workflow_name] = []
#                         if d not in block_data[case][sample_id][workflow_name]:
#                             block_data[case][sample_id][workflow_name].append(d)
                    
#     return block_data                
