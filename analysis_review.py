# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 16:36:23 2024

@author: rjovelin
"""

import sqlite3
import json
import argparse
import os
import hashlib
import itertools
from collections import deque
from generate_assays import generate_templates, list_qc_workflows, extract_assay_workflows
from db_helper import connect_to_db, define_columns, initiate_db, insert_multiple_records, \
    delete_unique_record, delete_multiple_records
from data_helper import load_data, clean_up_data, clean_up_workflows, is_case_info_complete
from commons import get_cases_md5sum, find_sequencing_attributes, convert_to_bool, \
    list_case_workflows, get_donor_name, compute_md5, case_to_update    
    



def collect_sample_workflows(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all workflows for all samples of a given case
        
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case   
    '''
    
    D = {}
    
    for d in case_data['workflow_runs']:
        wfrun_id = d['wfrunid']
        workflow = d['wf']
        limskeys = d['limsIds'].split(',')
        sequencing_attributes = find_sequencing_attributes(limskeys, case_data)
        for limskey in sequencing_attributes:
            tissue_type = sequencing_attributes[limskey]['tissue_type']
            #tissue_origin = sequencing_attributes[limskey]['tissue_origin']
            library_type = sequencing_attributes[limskey]['library_type']
            #group_id = sequencing_attributes[limskey]['group_id'] 
            sample = sequencing_attributes[limskey]['sample']
                    
            if sample not in D:
                D[sample] = {'tissue_type': tissue_type, 'library_type': library_type, 'workflows': [{'workflow': workflow, 'wfrun_id': wfrun_id}]}
            else:
                if {'workflow': workflow, 'wfrun_id': wfrun_id} not in D[sample]['workflows']:
                    D[sample]['workflows'].append({'workflow': workflow, 'wfrun_id': wfrun_id})            
        
    return D


    
def group_anchor_by_sample_type(anchor_samples):
    '''
    (dict) -> dict
    
    Returns a dictionary of anchor workflow ids for each sample type
        
    Parameters
    ----------
    - anchor_samples (dict): Dictionary with anchor workflow ids and their corresponding samples 
    '''
    
    D = {}
    
    for workflow in anchor_samples:
        for sample_type in anchor_samples[workflow]:
            if sample_type in D:
                D[sample_type].append(workflow)
            else:
                D[sample_type] = [workflow]
            D[sample_type] = list(set(D[sample_type]))
    
    return D


 
def group_anchor_workflows(anchor_samples):
    '''
    (dict) -> list
    
    Returns a list of tuples each containing a combination of anchor workflow ids
      
    Parameters
    ----------
    - anchor_samples (dict): Dictionary with anchor workflow ids and their corresponding samples 
    '''

    # group anchors bysample type
    D = group_anchor_by_sample_type(anchor_samples)
    # list all anchors according to the assay
    L = [D[i] for i in D]
    # find all combinations of anchors
    C = list(itertools.product(*L))    
 
    return C
 
   
def breadth_first_search(graph, start_node):
    '''
    (dict, str) -> list
    
    Returns a list of all the workflow ids connected found by traversing
    the graph starting at start_node
        
    Parameters
    ----------
    - graph (dict): Adjency list of workflow relationships
    - start_node (str): workflow_id used to start traversing the graph
    '''
    
    visited = set()
    queue = deque([start_node])
    traversal_order = []

    while queue:
        node = queue.popleft() # Dequeue from the left
        if node not in visited:
            traversal_order.append(node)
            visited.add(node)
            # Enqueue unvisited neighbors
            if node in graph:
                for neighbor in graph[node]:
                    if neighbor not in visited:
                        queue.append(neighbor)
    return traversal_order
 

def exclude_anchors(groups, focus_group):
    '''
    (list, tuple) -> list
    
    Returns a list of anchor workflow ids from each group in groups
    that are not in the focus griup 
        
    Parameters
    ----------
    - groups (list): List of tuples, each containing a combination of anchor workflow ids
    - focus_group (tuple): Tuple of anchor workflow ids     
    '''
    
    L = []
    for i in groups:
        if i != focus_group:
            for j in i:
                L.append(j)
    exclude = [i for i in L if i not in focus_group]
    
    return exclude
    
    
def exclude_workflows(groups, focus_group, parent_to_children_workflows):
    '''
    (list, tuple, dict) -> list

    Returns a list of anchor workflow ids in groups not present in the focus group
    including their immediate downstream workflows 

    Parameters
    ----------
    - groups (list): List of tuples, each containing a combination of anchor workflow ids
    - focus_group (tuple): Tuple of anchor workflow ids     
    - parent_to_children_workflows (dict): Dictionary with workflow ids and their immediate downstream workflows
    '''

    # exclude anchor workflows not in the focus group
    anchors = exclude_anchors(groups, focus_group)
    
    # exclude anchor workflows and their immediate downstream workflows
    excluded = []
    
    for i in parent_to_children_workflows:
        if i in anchors:
            excluded.append(i)
            excluded.extend(parent_to_children_workflows[i])    
    
    excluded = list(set(excluded))
    
    return excluded


def create_adjency_list(parent_to_children_workflows, excluded):
    '''
    (dict, list) -> dict    
    
    Returns a graph representation of the connected workflow ids without the excluded 
    ids
        
    Parameters
    ----------
    - parent_to_children_workflows (dict): Dictionary with workflow ids and their immediate downstream workflows
    - excluded (list): List of workflow ids to exclude from the adjency list 
    '''
    
    # create a dict of parent to children, removing parents without defined children workflows
    D = {}
    
    for i in parent_to_children_workflows:
        if 'NA' not in parent_to_children_workflows[i] and i != 'NA':
             if i not in excluded:
                 children = []
                 for j in parent_to_children_workflows[i]:
                     if j not in excluded and j != 'NA':
                         children.append(j)
                 if children:
                     D[i] = children
    

    return D




def identify_fastq_workflows(workflow_info, fastq_workflows):
    '''
    (dict, list) -> list
    
    Returns a list of workflow ids of fastq-generating workflows
        
    Parameters
    ----------
    - workflow_info (dict): Dictionary with workflow names mapped to workflow ids
    - fastq_workflows (list): List of fastq-generating workflow names
    '''
        
    sequences = [i for i in workflow_info if workflow_info[i] in fastq_workflows]

    return sequences


   
def find_related_workflows(groups, group, parent_to_children_workflows, workflow_info, fastq_workflows):
    '''
    (list, tuple, dict, dict, list) -> list
    
    Returns a list of all connected workflows from top (fastq-generating workflows)
    to the most downstream analysis workflows, traversing the call ready workflows
    in group
        
    Parameters
    ----------
    - groups (list): List of tuples, each containing a combination of anchor workflow ids
    - focus_group (tuple): Tuple of anchor workflow ids     
    - parent_to_children_workflows (dict): Dictionary with workflow ids and their immediate downstream workflows
    - workflow_info (dict): Dictionary with workflow names mapped to workflow ids
    - fastq_workflows (list): List of fastq-generating workflow names
    '''
    
    L = []
    
    # list al workflow ids of fastq-generating workflows
    sequences = identify_fastq_workflows(workflow_info, fastq_workflows)
    
    # exclude anchor worfklows and their immediate downstream workflows
    # that are not in group
    excluded = exclude_workflows(groups, group, parent_to_children_workflows)
    # create an adjency list of all the workflows without the excluded worfklows
    matrix = create_adjency_list(parent_to_children_workflows, excluded)
    # find all connected workflows starting at each sequence workflow
    for workflow_id in sequences:
        connected = breadth_first_search(matrix, workflow_id)
        # check that anchor workflow in group are among the connected workflows    
        keep = [i in connected for i in group]
        if any(keep):
            L.extend(connected)
    
    L = list(set(L))  
    
    return L
    
    

def get_downstream_workflows(parent_workflows):
    '''
    (dict) -> dict
    
    Returns a dictionary with list of downstream workflow ids for each workflow id
    
    Parameters
    ----------
    - parent_workflows (dict): Dictionary with child-parents workflow relationship
    '''
    
    D = {}
        
    for i in parent_workflows:
        for j in parent_workflows[i]:
            if j in D:
                D[j].append(i)
            else:
                D[j] = [i]
    return D
    
    

def convert_assay_to_template(assay):
    '''
    (dict) -> dict

    Returns a dictionary for collecting analysis data from a specific assay     
    
    Paramaters
    ----------
    - assay (dict): Dictionary with assay specifying expected anchors, samples and workflows
    '''
        
    template = {}
    
    for i in assay:
        template[i] = {}
        if i == 'Samples':
            for j in assay[i]:
                template[i][j] = {}    
                for k in assay[i][j]:
                    template[i][j][k] = '' 
                template[i][j]['sample'] = ''                 
        elif i == 'Analysis':
            for j in assay[i]:
                template[i][j] = []
        elif i == 'Data':
            for j in assay[i]:
                template[i][j] = []
        
    return template
    
    
    
def extract_workflow_information(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with information about all workflows for a case donor
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data   
    '''    

    D = {}

    for d in case_data['workflow_runs']:
        name = d['wf']
        wfrun = d['wfrunid']
        D[wfrun] = name
    
    return D
        
 
       
def identify_sample_type(tissue_type, library_source, negate_tissue_type):
    '''
    (str, str, bool) -> str
    
    Returns the sample type (eg NormalWG, TumourWT) based on the library source
    and tissue origin
        
    Parameters
    ----------
    - tissue_type (str): The tissue type, indicating normal or tumor origin
    - library_source (str): The code of the library
    - negate_tissue_type (bool): True if the tissue type indicates a tumor origin
    '''
        
    sample_name = ''
    if tissue_type == 'R' and negate_tissue_type == False:
        # normal sample
        sample_name = 'Normal' + library_source 
    elif tissue_type == 'R' and negate_tissue_type:
        # tumour sample
        sample_name = 'Tumour' + library_source
    
    if sample_name == '':
        print(tissue_type, library_source, negate_tissue_type)
    assert sample_name != ''
      
    
    return sample_name



def map_samples_to_workflows(samples_workflows):
    '''
    (dict) -> dict  
    
    Returns a dictionary with samples mapped to each workflow id in a case
            
    Parameters
    ----------
    - samples_workflows (dict): Dictionary with all the workflow ids of each sample in a case
    '''
    
    D = {}
    
    for sample in samples_workflows:
        tissue_type = samples_workflows[sample]['tissue_type']
        library_type = samples_workflows[sample]['library_type']    
        if tissue_type == 'R':
            negate_tissue_type = False
        else:
            negate_tissue_type = True
            tissue_type = 'R'       
        for d in samples_workflows[sample]['workflows']:
            workflow_id = d['wfrun_id']
            if workflow_id not in D:
                D[workflow_id] = {}
            if sample not in D[workflow_id]:
                D[workflow_id][sample] = {'library_type': library_type, 'tissue_type': tissue_type, 'negate_tissue_type': negate_tissue_type}
            else:
                assert D[workflow_id][sample]['library_type'] == library_type
                assert D[workflow_id][sample]['tissue_type'] == tissue_type
                assert D[workflow_id][sample]['negate_tissue_type'] == negate_tissue_type
    
    return D
    


def add_samples_to_template(connected_workflows, workflows_to_samples, template):
    '''
    (list, dict, dict) -> None
    
    Adds in place the samples information to the template
       
    Parameters
    ----------
    - connected_workflows (list): List of connected workflows for a group of call-ready workflows
    - workflows_to_samples (dict): Dictionary with samples mapped to each workflow id in a case
    - template (dict): Dictionary with data to collect
    '''
    
    for workflow_id in connected_workflows:
        for sample in workflows_to_samples[workflow_id]:
            library_type = workflows_to_samples[workflow_id][sample]['library_type']
            tissue_type = workflows_to_samples[workflow_id][sample]['tissue_type']
            negate_tissue_type = workflows_to_samples[workflow_id][sample]['negate_tissue_type']
            sample_type = identify_sample_type(tissue_type, library_type, negate_tissue_type)
            k = {'library_type': library_type,
                 'tissue_type': tissue_type,
                 'negate_tissue_type': negate_tissue_type,
                 'sample': sample}
                    
            for i in k:
                if bool(template['Samples'][sample_type][i]):
                    assert template['Samples'][sample_type][i] == k[i]
                else:
                    template['Samples'][sample_type][i] = k[i]
    
    


def add_analysis_to_template(connected_workflows, workflow_info, workflows_to_samples, child_to_parents_workflows, template):
    '''
    (list, dict, dict, dict, dict) -> None
            
    Adds in place the analysis information to the template
       
    Parameters
    ----------
    - connected_workflows (list): List of connected workflows for a group of call-ready workflows
    - workflow_info (dict): Dictionary with workflow name mapped to worfklow id
    - workflows_to_samples (dict): Dictionary with samples mapped to each workflow id in a case
    - child_to_parents_workflows (dict): Dictionary with parent worfklows mapped to worfklow ids
    - template (dict): Dictionary with data to collect
    '''
        
    for workflow_id in connected_workflows:
        workflow_name = workflow_info[workflow_id]
        samples = [{sample: workflows_to_samples[workflow_id][sample]} for sample in workflows_to_samples[workflow_id]]
        inputs = child_to_parents_workflows[workflow_id]
            
        d = {'workflow_id': workflow_id,
             'workflow_name': workflow_name,
             'samples': samples,
             'inputs': inputs}
        
        if workflow_name in template['Analysis']:
            # do not record identical workflows
            if d not in template['Analysis'][workflow_name]:
                template['Analysis'][workflow_name].append(d)
        

    
    
def add_data_to_template(connected_workflows, workflow_info, workflows_to_samples, child_to_parents_workflows, fastq_workflows, template): 
    '''
    (dict, dict, dict, dict) -> None
    
    Adds in place the sequencing and alignment data to the template
    
    Parameters
    ----------
    - connected_workflows (list): List of connected workflows for a group of call-ready workflows
    - workflow_info (dict): Dictionary with workflow name mapped to worfklow id
    - workflows_to_samples (dict): Dictionary with samples mapped to each workflow id in a case
    - child_to_parents_workflows (dict): Dictionary with parent worfklows mapped to worfklow ids
    - fastq_workflows (list): List of fastq-generating workflow names
    - template (dict): Dictionary with data to collect
    '''
        
    data_workflows = fastq_workflows + ['bwamem']
    
    for workflow_id in connected_workflows:
        workflow_name = workflow_info[workflow_id]
        samples = [{sample: workflows_to_samples[workflow_id][sample]} for sample in workflows_to_samples[workflow_id]]
        inputs = child_to_parents_workflows[workflow_id]
        
        d = {'workflow_id': workflow_id,
             'workflow_name': workflow_name,
             'samples': samples,
             'inputs': inputs}
        
        if workflow_name in fastq_workflows:
            if d not in template['Data']['Sequencing']:
                template['Data']['Sequencing'].append(d)
                
        if workflow_name in data_workflows:
            for k in samples:
                for sample in k:
                    library_type = k[sample]['library_type']
                    tissue_type = k[sample]['tissue_type']
                    negate_tissue_type = k[sample]['negate_tissue_type']
                    sample_type = identify_sample_type(tissue_type, library_type, negate_tissue_type)
                    sample_type = sample_type + 'align'       
            
                    if d not in template['Data'][sample_type]:
                        template['Data'][sample_type].append(d)
            
    
def add_anchors_to_template(connected_workflows, workflow_info, workflows_to_samples, template):
    '''
    (dict, dict, dict) -> None
            
    Adds in place the anchor workflow information to the template
       
    Parameters
    ----------
    - connected_workflows (list): List of connected workflows for a group of call-ready workflows
    - workflow_info (dict): Dictionary with workflow name mapped to worfklow id
    - workflows_to_samples (dict): Dictionary with samples mapped to each workflow id in a case
    - template (dict): Dictionary with data to collect
    '''
        
    for workflow_id in connected_workflows:
        workflow_name = workflow_info[workflow_id]
        if 'bammergepreprocessing' in workflow_name.lower() or 'star_call_ready' in workflow_name.lower():
            for sample in workflows_to_samples[workflow_id]:
                library_type = workflows_to_samples[workflow_id][sample]['library_type']
                tissue_type = workflows_to_samples[workflow_id][sample]['tissue_type']
                negate_tissue_type = workflows_to_samples[workflow_id][sample]['negate_tissue_type']    
                sample_type = identify_sample_type(tissue_type, library_type, negate_tissue_type)
                if sample_type not in template['Anchors']:
                    template['Anchors'][sample_type] = []
                if workflow_id not in template['Anchors'][sample_type]:
                    template['Anchors'][sample_type].append(workflow_id)
        

def fill_group_template(assay, connected_workflows, workflow_info, workflows_to_samples, child_to_parents_workflows, fastq_workflows):
    '''
    (dict, dict, dict, dict, dict) -> list
    
    Returns a list of templates with analysis data for each block
    
    Parameters
    ----------
    - assay (dict): Dictionary with assay information
    - connected_workflows (list): List of connected workflows for a group of call-ready workflows
    - workflow_info (dict): Dictionary with workflow name mapped to worfklow id
    - workflows_to_samples (dict): Dictionary with samples mapped to each workflow id in a case
    - child_to_parents_workflows (dict): Dictionary with parent worfklows mapped to worfklow ids
    - fastq_workflows (list): List of fastq-generating workflow names
    '''

    # create template from the assay
    template = convert_assay_to_template(assay)
    # add samples section to template
    add_samples_to_template(connected_workflows, workflows_to_samples, template)
    # add the analysis section
    add_analysis_to_template(connected_workflows, workflow_info, workflows_to_samples, child_to_parents_workflows, template)    
    # add alignment data
    add_data_to_template(connected_workflows, workflow_info, workflows_to_samples, child_to_parents_workflows, fastq_workflows, template)    
    # add anchor workflows
    add_anchors_to_template(connected_workflows, workflow_info, workflows_to_samples, template)
    
    return template


def evaluate_template_samples(template, assay):
    '''
    (dict, dict) -> (bool, str)
    
    Returns a boolean indicating if the template has missing samples, 
    and the corresponding error message
    
    Parameters
    ----------
    - template (dict): Dictionary collecting data for the given assay
    - assay (dict): Dictionary with assay information
    '''
    
    complete = True
    error = ''
    
    # evaluate expected samples
    for i in template['Samples']:
        for j in template['Samples'][i]:
            if template['Samples'][i][j] == '':
                complete = False
    if complete == False:
        error = 'missing samples'

    return complete, error


def evaluate_template_raw_data(template, assay):
    '''
    (dict, dict) -> (bool, str)
    
    Returns a boolean indicating if the raw data (sequences and alignments)
    is mmissing from the template, and the corresponding error message
    
    Parameters
    ----------
    - template (dict): Dictionary collecting data for the given assay
    - assay (dict): Dictionary with assay information
    '''
    
    complete = True
    error = ''

    # evaluate sequencing and alignments
    for i in template['Data']:
        if len(template['Data'][i]) == 0:
            complete = False
    if complete == False:
        error = 'missing raw data'        

    return complete, error




def evaluate_missing_analyses(template, assay):
    '''
    (dict, dict) -> (bool, list)
    
    Returns a boolean indicating if expected analyses workflows are missing
    from the template, and the corresponding error message
    
    Parameters
    ----------
    - template (dict): Dictionary collecting data for the given assay
    - assay (dict): Dictionary with assay information
    '''

    complete = True
    error = ''
    # evaluate analyses
    # check if each expected workflow has data
    missing = []
    for workflow in template['Analysis']:
        if len(template['Analysis'][workflow]) == 0:
            complete = False
            missing.append(workflow)
    if complete == False:
        missing = sorted(list(set(missing)))
        error = 'missing {0} workflows ({1})'.format(len(missing), ':'.join(missing))        

    return complete, error


def evaluate_extra_analyses(template, assay):
    '''
    (dict, dict) -> (bool, str)
    
    Returns a boolean indicating if there are multiple instances of the expected
    analysis workflows, and the corresponding error message
    
    Parameters
    ----------
    - template (dict): Dictionary collecting data for the given assay
    - assay (dict): Dictionary with assay information
    '''

    complete = True
    error = []

    # check if each expected workflow has a single record
    # skip call ready workflow because they may have a record each sample

    for workflow in template['Analysis']:
        
        # skipp bwamem and star lane level which can have multiple workflows
        if workflow.lower() not in ['bwamem', 'star_lane_level']:
            samples = {}
            extra = 0
            for d in template['Analysis'][workflow]:
                workflow_id = d['workflow_id']
                for k in d['samples']:
                    for sample in k:
                        if sample not in samples:
                            samples[sample] = []
                        samples[sample].append(workflow_id)
                        samples[sample] = list(set(samples[sample]))        
            for sample in samples:
                if len(samples[sample]) > 1:
                    extra += len(samples[sample]) - 1
        
            if extra:
                err = '{0}: {1} extra interations'.format(workflow, extra)
                error.append(err)
                complete = False
    
    error = ';'.join(error)
    
    return complete, error


def evaluate_assay(template, assay):
    '''
    (dict, dict) -> (bool, str)
    
    Returns a boolean indicating if the template holds all the expected data, 
    and the error message if template is incomplete
    
    Parameters
    ----------
    - template (dict): Dictionary collecting data for the given assay
    - assay (dict): Dictionary with assay information
    '''
    
    # records the completeness and error messages for the different template sections 
    complete = []
    error = []
    
    # evaluate samples
    cpl, err = evaluate_template_samples(template, assay)    
    complete.append(cpl)
    error.append(err)
    
    # evaluate raw data (sequences + alignments)
    cpl, err = evaluate_template_raw_data(template, assay)
    complete.append(cpl)
    error.append(err)
    
    # evaluate missing analyses workflows
    cpl, err = evaluate_missing_analyses(template, assay)
    complete.append(cpl)
    error.append(err)
    
    # # evaluate extra workflow occurences
    # cpl, err = evaluate_extra_analyses(template, assay)
    # complete.append(cpl)
    # error.append(err)
        
    # evaluate all sections
    while '' in error:
        error.remove('')
    error = ';'.join(error)
    complete = all(complete)
    
    return complete, error
    
    
def extract_anchor_samples(assay):
    '''
    (dict) -> dict
    
    Returns a dictionary with expected sample information for each anchor workflow
    
    Parameters
    ----------
    - assay (dict): Dictionary with assay information
    '''    
        
    D = {}
    
    for sample in assay['Anchors']:
        workflow = assay['Anchors'][sample]['workflows']
        if workflow not in D:
            D[workflow] = {}
        D[workflow][sample] = assay['Samples'][sample]    
    
    return D
    
    
    
def get_anchor_samples(samples_workflows, anchor_workflows):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with samples for each anchor workflow
    
    Parameters
    ----------
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case 
    - anchor_workflows (dict): Dictionary with expected sample information for each anchor workflow
    '''
    
    D = {}

    for sample in samples_workflows:
        tissue_type = samples_workflows[sample]['tissue_type']
        library_type = samples_workflows[sample]['library_type']
        for i in samples_workflows[sample]['workflows']:
            workflow_name = i['workflow']
            workflow_id = i['wfrun_id']
            if workflow_name in anchor_workflows:
                # find the sample name
                for j in anchor_workflows[workflow_name]:
                    if library_type == anchor_workflows[workflow_name][j]['library_type']:
                        if (tissue_type == anchor_workflows[workflow_name][j]['tissue_type'] and \
                            anchor_workflows[workflow_name][j]['negate_tissue_type'] == False) or \
                            (tissue_type != anchor_workflows[workflow_name][j]['tissue_type'] and \
                             anchor_workflows[workflow_name][j]['negate_tissue_type']):
                                # normal sample
                                sample_name = j
                                if  workflow_id not in D:
                                    D[workflow_id] = {}
                                if sample_name not in D[workflow_id]:
                                    D[workflow_id][sample_name] = []
                                D[workflow_id][sample_name].append(sample)
            
    return D                
                
 
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
            D[workflow].extend(children_workflows)
        else:
            D[workflow] = children_workflows
        D[workflow] = list(set(D[workflow]))
        # record the parent-child relationships of each parent and current workflow 
        for workflow_run in parent_workflows:
            if workflow_run in D:
                D[workflow_run].append(workflow)
            else:
                D[workflow_run] = [workflow]
            D[workflow_run] = list(set(D[workflow_run]))

    return D


def is_moh_case(case_data):
    '''
    (dict) -> bool    
    
    Returns True if any of the projects associated with the case are MOH
    
    Parameters
    ----------
    - case_data (str): Dictionary with case information
    '''
    
    L = []
    
    for d in case_data['project_info']:
        L.append('MOH Full Pipeline' in d['deliverables'])
    
    return any(L)             


def get_moh_assay(assay_name):
    '''
    (str) -> dict
    
    Returns a dictionary with the WGS or WGTS MOH assay 
    
    Parameters
    ----------
    - assay_name (str): Name of the MOH assay
    '''
        
    if 'WGTS' in assay_name:
        assays = {'Samples': {'TumourWT': {'library_type': 'WT',
                                            'tissue_type': 'R',
                                            'negate_tissue_type': True},
                               'TumourWG': {'library_type': 'WG',
                                            'tissue_type': 'R',
                                            'negate_tissue_type': True},
                               'NormalWG': {'library_type': 'WG',
                                            'tissue_type': 'R',
                                            'negate_tissue_type': False}},
                   'Data': {'Sequencing': {'workflows': ['bcl2fastq'], 'inputs': []},
                            'TumourWTalign': {'samples': ['TumourWT'],
                                              'inputs': [],
                                              'data': ['Sequencing'],
                                              'workflows': ['star_lane_level']},
                            'TumourWGalign': {'samples': ['TumourWG'],
                                              'inputs': [],
                                              'data': ['Sequencing'],
                                              'workflows': ['bwaMem']},
                            'NormalWGalign': {'samples': ['NormalWG'],
                                              'inputs': [],
                                              'data': ['Sequencing'],
                                              'workflows': ['bwaMem']}},
                   'Analysis':
                       {'mutect2_matched': {'samples': ['NormalWG', 'TumourWG'],
                                            'inputs': ['bamMergePreprocessing_by_sample']},
                       'variantEffectPredictor_matched': {'samples': ['NormalWG', 'TumourWG'],
                                                          'inputs': ['mutect2_matched']},
                       'delly_matched': {'samples': ['NormalWG', 'TumourWG'],
                                         'inputs': ['bamMergePreprocessing_by_sample']}, 
                       'gridss': {'samples': ['NormalWG', 'TumourWG'],
                                  'inputs': ['bamMergePreprocessing_by_sample']},
                       'purple': {'samples': ['NormalWG', 'TumourWG'],
                                  'inputs': ['mutect2_matched', 'bamMergePreprocessing_by_sample', 'gridss']},
                       'bamMergePreprocessing_by_sample': {'samples': ['NormalWG', 'TumourWG'],
                                                           'inputs': ['bwaMem']},
                       'bwaMem': {'samples': ['NormalWG', 'TumourWG'],
                                  'inputs': ['bcl2fastq']},
                       'rsem': {'samples': ['TumourWT'],
                                'inputs': ['star_call_ready']},
                       'star_call_ready': {'samples': ['TumourWT'],
                                           'inputs': ['bcl2fastq']},
                       'starfusion': {'samples': ['TumourWT'],
                                      'inputs': ['star_call_ready']},
                       'arriba': {'samples': ['TumourWT'],
                                  'inputs': ['star_call_ready']},
                       'msisensor': {'samples': ['NormalWG', 'TumourWG'],
                                     'inputs': ['bamMergePreprocessing_by_sample']},
                       'star_lane_level': {'samples': ['TumourWT'],
                                           'inputs': ['bcl2fastq']},
                       'mavis': {'samples': ['TumourWT', 'NormalWG', 'TumourWG'],
                                 'inputs': ['bamMergePreprocessing_by_sample',
                                            'star_call_ready', 'starfusion',
                                            'arriba', 'delly_matched']},
                       'hrDetect': {'samples': ['NormalWG', 'TumourWG'],
                                    'inputs': ['purple', 'mutect2_matched']},
                       'haplotypeCaller': {'samples': ['NormalWG', 'TumourWG'],
                                           'inputs': ['bamMergePreprocessing_by_sample']}},
                   'Anchors': {'TumourWT': {'workflows': 'star_call_ready'},
                               'TumourWG': {'workflows': 'bamMergePreprocessing_by_sample'},
                               'NormalWG': {'workflows': 'bamMergePreprocessing_by_sample'}}}
      

    elif 'WGS' in assay_name:
        assays = {'Samples': {'TumourWG': {'library_type': 'WG',
                                           'tissue_type': 'R',
                                           'negate_tissue_type': True},
                              'NormalWG': {'library_type': 'WG',
                                           'tissue_type': 'R',
                                           'negate_tissue_type': False}},
                  
                  'Data': {'Sequencing': {'workflows': ['bcl2fastq'], 'inputs': []},
                           'TumourWGalign': {'samples': ['TumourWG'],
                                             'inputs': [],
                                             'data': ['Sequencing'],
                                             'workflows': ['bwaMem']},
                           'NormalWGalign': {'samples': ['NormalWG'],
                                             'inputs': [],
                                             'data': ['Sequencing'],
                                             'workflows': ['bwaMem']}},
                  'Analysis':
                      {'bwaMem': {'samples': ['NormalWG', 'TumourWG'],
                                 'inputs': ['bcl2fastq']},
                       'bamMergePreprocessing_by_sample': {'samples': ['NormalWG', 'TumourWG'],
                                                           'inputs': ['bwaMem']},
                       'mutect2_matched': {'samples': ['NormalWG', 'TumourWG'],
                                            'inputs': ['bamMergePreprocessing_by_sample']},
                       'variantEffectPredictor_matched': {'samples': ['NormalWG', 'TumourWG'],
                                                          'inputs': ['mutect2_matched']},
                       'delly_matched': {'samples': ['NormalWG', 'TumourWG'],
                                         'inputs': ['bamMergePreprocessing_by_sample']}, 
                       'gridss': {'samples': ['NormalWG', 'TumourWG'],
                                  'inputs': ['bamMergePreprocessing_by_sample']},
                       'purple': {'samples': ['NormalWG', 'TumourWG'],
                                  'inputs': ['mutect2_matched', 'bamMergePreprocessing_by_sample', 'gridss']},
                       'msisensor': {'samples': ['NormalWG', 'TumourWG'],
                                     'inputs': ['bamMergePreprocessing_by_sample']},
                       'mavis': {'samples': ['NormalWG', 'TumourWG'],
                                 'inputs': ['bwaMem', 'delly_matched']}},
                       'hrDetect': {'samples': ['NormalWG', 'TumourWG'],
                                    'inputs': ['purple', 'mutect2_matched']},
                       'haplotypeCaller': {'samples': ['NormalWG', 'TumourWG'],
                                           'inputs': ['bamMergePreprocessing_by_sample']},
                  'Anchors': {'TumourWG': {'workflows': 'bamMergePreprocessing_by_sample'},
                              'NormalWG': {'workflows': 'bamMergePreprocessing_by_sample'}}}

    return assays
    
    

def generate_cache(provenance_data_file, assay_config_file, waterzooi_database, pinery, database, table='templates'):
    '''
    (str, str, str, str, str, str) -> None 
    
    Generates sqlite database with templates and review for all projects and cases in the
    provenance data file
    
    Parameters
    ----------
    - provenance_data_file (str): Path to the file with production data extracted from Shesmu
    - assay_config_file (str): Path to the assay config file
    - database (str): Path to the waterzooi database
    - database (str): Path to the sqlite database
    - pinery (str): URL to Pinery assay endpoint
    - table (str): Table in database storing the analysis data
    '''
    
    
    # defined fastq generating workflows
    fastq_workflows = ['bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
    
    # list QC workflows
    assay_configurations = extract_assay_workflows(assay_config_file)
    # make a list of QC workflows
    qc_workflows = list_qc_workflows(assay_configurations)
      
    # create database if file doesn't exist
    if os.path.isfile(database) == False:
        initiate_db(database, 'analysis_review', ['templates'])
    print('initiated database')    
    
    # collect the recorded md5sums of the donor data from the database
    recorded_md5sums = get_cases_md5sum(database, table = 'templates')
    print('pulled md5sums from database')
    
    # load data from file
    provenance_data = load_data(provenance_data_file)
    print('loaded data')
    
    # generate assays
    assays = generate_templates(waterzooi_database, assay_config_file, pinery)
        
    # track all cases in production
    P = []
          
    for case_data in provenance_data:
        case_id = case_data['case']
        P.append(case_id)
        assay = {}
        # check that no information is missing
        if is_case_info_complete(case_data):
            # record all case templates
            L = []
            # remove workflows that do not belong to the case
            case_data = clean_up_workflows(case_data)
            # compute the md5sum of the case info
            md5sum = compute_md5(case_data)
            # determine if case needs to be updated
            donor = get_donor_name(case_data)
            assay_name = case_data['assay']        
            
            if case_to_update(recorded_md5sums, case_id, md5sum):
                # open connection to database
                conn = connect_to_db(database)
                # delete case info from table
                delete_unique_record(case_id, conn, database, 'templates', 'case_id')
                conn.close()
                                
                            
                ## temporary hack to use manually defined MOH assays
                if is_moh_case(case_data):
                    assay = get_moh_assay(assay_name)
                else:
                    if assay_name in assays:
                        assay = assays[assay_name]
                        
                if assay:
                    # get all the workflow information
                    workflow_info = extract_workflow_information(case_data)
                    # get the anchor workflows and their expected samples
                    anchor_workflows = extract_anchor_samples(assay)
                    # get workflows of all samples for the case
                    samples_workflows = collect_sample_workflows(case_data)
                    # map samples to each workflow
                    workflows_to_samples = map_samples_to_workflows(samples_workflows)
                    # find the children of each workflow
                    parent_to_children_workflows = collect_workflow_relationships(case_data)
                    # find the parents of each workflow
                    child_to_parents_workflows = get_downstream_workflows(parent_to_children_workflows)
                    anchor_samples = get_anchor_samples(samples_workflows, anchor_workflows)
                                        
                    # # make groups of samples (can be 1, 2 or more samples depending on the assay)       
                    # groups = group_samples(anchor_samples, assay)
                    
                    # group anchor workflows
                    groups = group_anchor_workflows(anchor_samples)
                    
                    # fill the templates for each group
                    templates = []
                    for group in groups:
                        connected = find_related_workflows(groups, group, parent_to_children_workflows, workflow_info, fastq_workflows)
                        # remove QC workflows
                        connected = [i for i in connected if workflow_info[i] not in qc_workflows]
                        template = fill_group_template(assay, connected, workflow_info, workflows_to_samples, child_to_parents_workflows, fastq_workflows)
                        templates.append(template)
                    
                    
                    # check if the case has multiple projects
                    for i in case_data['project_info']:
                        project_id = i['project']
                        # evaluate template
                        for template in templates:
                            valid, error = evaluate_assay(template, assay)
                            L.append([case_id, donor, project_id, assay_name, json.dumps(template), str(int(valid)), error, md5sum])
                                     
                else:
                    for i in case_data['project_info']:
                        project_id = i['project']
                        L.append([case_id, donor, project_id, assay_name, json.dumps({}), str(0), 'no_assay', md5sum])
    
            if L:
                # connect to database
                conn = connect_to_db(database)
                # insert records
                insert_multiple_records(L, conn, database, 'templates', define_columns('analysis_review')['templates']['names'])
                # close database
                conn.close()
                
                
    
    
    # delete data for donors not in the provenance report
    if P:
        # make a list of cases in database that are not in production
        conn = connect_to_db(database)
        data = conn.execute('SELECT case_id FROM templates').fetchall()
        all_cases = [i['case_id'] for i in data]
        to_remove = [i for i in all_cases if i not in P]
        if to_remove:
            delete_multiple_records(to_remove, conn, database, 'templates', 'case_id')
        conn.close()



generate_cache('provenance_reporter.json', 'enabled_workflows.json', 'waterzooi_db_case.db', 'http://pinery.gsi.oicr.on.ca/assays', 'analysis_review_case.db', table='templates')






# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(prog = 'analysis_review.py', description='Script to generate the analyais review cache')
#     parser.add_argument('-p', '--provenance', dest = 'provenance', help = 'Path to the provenance reporter data json', required=True)
#     parser.add_argument('-db', '--database', dest='database', help='Path to the analysis review database', required=True)    
#     parser.add_argument('-cq', '--calcontaqc', dest = 'calcontaqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest', help = 'Path to the merged rnaseq calculateContamination database. Default is /scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest')


    
#     # get arguments from the command line
#     args = parser.parse_args()
#     # generate sqlite cache
#     generate_cache(args.provenance, args.database)    
    