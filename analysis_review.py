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
import copy
from generate_assays import generate_templates
from db_helper import connect_to_db, define_columns, initiate_db, insert_multiple_records, \
    delete_unique_record, delete_multiple_records
from data_helper import load_data, clean_up_data, clean_up_workflows, is_case_info_complete
from commons import get_cases_md5sum, find_sequencing_attributes, convert_to_bool, \
    list_case_workflows, get_donor_name, compute_md5, case_to_update    
    





    
def is_project_active(case_data):
    '''
    (dict) -> bool
    
    Returns True if the project is active and False otherwise
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data    
    '''
    
    active = [convert_to_bool(i['isActive']) for i in case_data['project_info']]
    return all(active)






def collect_sample_workflows(case_data):
    '''
    (dict, str) -> dict
    
    Returns a dictionary with all workflows for all samples of a given case
        
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case   
    '''
    
    D = {}
    
    #case = case_data['case']
    #donor = get_donor_name(case_data)
 
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



def get_assay_samples_groups(assay):
    '''
    (dict) -> list
    
    Return a list of sample types extracted from the Anchors section of the assay,
    grouped by workflow
    
    Parameters
    ----------
    - assay (dict): Dictionary with assay 
    '''
    
    # get the expected group of samples based on the assay workflows
    
    L = []
    d = {}
    for sample in assay['Anchors']:
        workflow = assay['Anchors'][sample]['workflows']
        if workflow in d:
            d[workflow].append(sample)
        else:
            d[workflow] = [sample]
    
    for workflow in d:
        if len(d[workflow]) > 1:
            L.append(d[workflow])
        else:
            L.extend(d[workflow])
       
    return L
    
    
def map_samples_to_sample_type(anchor_samples):
    '''
    (dict) -> dict
    
    Returns a dictionary of sample type for each sample
    
    Parameters
    ----------
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    '''
        
    D = {}
    for workflow in anchor_samples:
        for sample_type in anchor_samples[workflow]:
            for sample in anchor_samples[workflow][sample_type]:
                D[sample] = sample_type
    return D



def group_samples_by_sample_type(anchor_samples):
    '''
    (dict) -> dict
    
    Returns a dictionary with lists of samples for each sample type
    
    Parameters
    ----------
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    '''
    
    # consider not taking the reverse but instead listing all the samples, even if duplicates
    # shopuld correspond to the number of workflows, even if same sample for multiple worklows
    
    D = {}
    for workflow in anchor_samples:
        for sample_type in anchor_samples[workflow]:
            for sample in anchor_samples[workflow][sample_type]:
                if sample_type in D:
                    D[sample_type].append(sample)
                else:
                    D[sample_type] = [sample]

    return D
    
    
def group_samples(anchor_samples, assay):
    '''
    (dict, dict) -> list
    
    Returns a list of samples processed together in an assay
    
    Parameters
    ----------
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    - assay (dict): Dictionary with assay template
    '''

    # group samples by sample type
    D = group_samples_by_sample_type(anchor_samples)
    # list all samples according to the assay
    L = [D[i] for i in D]
    # find all combinations of samples
    C = list(itertools.product(*L))
    
    # group samples according to assay
    S = map_samples_to_sample_type(anchor_samples)
    A = get_assay_samples_groups(assay)
    groups = [] 
    for i in range(len(C)):
        groups.append(copy.deepcopy(A))
    for i in range(len(C)):
        for j in range(len(groups[i])):
            if type(groups[i][j]) == list:
                for k in range(len(groups[i][j])):
                    for s in C[i]:
                        if groups[i][j][k] == S[s]:
                            groups[i][j][k] = s
            else:
               for s in C[i]:
                   if groups[i][j] == S[s]:
                       groups[i][j] = s

    return groups                 
    


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
    
    
def identify_analysis_workflows(D, parent_workflows, children_workflows, anchors):
    '''
    (dict, dict, dict, list) -> list
    
    Returns a dictionary with analysis workflows grouped by sample pairs and blocks 
    (ie, originating from common call-ready workflows)
    
    Parameters
    ----------
    - D (dict): Dictionary with workflows of each normal, tumour sample pair
    - parent_worfklows (dict): Dictionary with parent workflows of each workflow for the normal, tumour sample pairs
    - children_workflows (dict): Dictionary with list of downstream workflow ids for each workflow id
    - anchors (list): List of anchor workflow run ids                           
    '''
    
    # track blocks
    blocks = {}

    # record all parent workflows
    L = []

    # get the analysis workflows immediately downstream of the anchor workflows
    for samples in D:
        for d in D[samples]:
            # get the workflow id
            wfrun_id = d['wfrun_id']
            # get the parent workflows
            parent = parent_workflows[wfrun_id]
            # check that parent workflows are in anchor list
            if parent and all([i in anchors for i in parent]):
                # get the anchor workflow
                parent_workflow = '.'.join(sorted(parent))
                if samples not in blocks:
                    blocks[samples] = {}
                if parent_workflow not in blocks[samples]:
                    blocks[samples][parent_workflow] = []
                # record workflow
                blocks[samples][parent_workflow].append(wfrun_id)
                # record anchor workflows
                #blocks[samples][parent_workflow].extend(parent_workflow.split('.'))
                blocks[samples][parent_workflow] = list(set(blocks[samples][parent_workflow]))
                L.append(wfrun_id)
    
    # find the downstream workflows
    for samples in blocks:
        for parent_workflow in blocks[samples]:
            downstream = []
            # get the downstream workflows of the workflow
            for workflow in blocks[samples][parent_workflow]:
                if workflow in children_workflows:
                    downstream.extend(children_workflows[workflow])
            blocks[samples][parent_workflow].extend(downstream)
                
    # record the anchor workflows and remove undefined workflows (when workflow have no parent or no child)
    for samples in blocks:
        for parent_workflow in blocks[samples]:
            blocks[samples][parent_workflow].extend(parent_workflow.split('.'))
            blocks[samples][parent_workflow] = list(set(blocks[samples][parent_workflow]))
            if 'NA' in blocks[samples][parent_workflow]:
                blocks[samples][parent_workflow].remove('NA') 
       
    return blocks                



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
                template[i][j] = {"workflows":[]}
        
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
        
 

def add_samples_to_template(assay, block, samples_workflows, template):
    '''
    (dict, dict, dict, dict) -> None
    
    Adds in place the samples information to the template
       
    Parameters
    ----------
    - assay (dict): Dictionary with assay information
    - block (dict): Dictionary with analyses organized in blocks of related workflows
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case
    - template (dict): Dictionary collecting data for the given assay
    '''
    
    samples = []
    for i in block:
        samples.extend(list(map(lambda x: x.strip(), i.split('|'))))
    # add samples information
    for i in samples:
        for j in assay['Samples']:
            if samples_workflows[i]['library_type'] == assay['Samples'][j]['library_type']:
                if (assay['Samples'][j]['negate_tissue_type'] == False and \
                    samples_workflows[i]['tissue_type'] == assay['Samples'][j]['tissue_type']) or \
                   (assay['Samples'][j]['negate_tissue_type'] and \
                    samples_workflows[i]['tissue_type'] != assay['Samples'][j]['tissue_type']):
                        template['Samples'][j]['sample'] = i
                        template['Samples'][j]['tissue_type'] = samples_workflows[i]['tissue_type']
                        template['Samples'][j]['library_type'] = samples_workflows[i]['library_type']
                        template['Samples'][j]['negate_tissue_type'] = assay['Samples'][j]['negate_tissue_type']
                        
                        

def add_analysis_to_template(block, workflows_information, parent_workflows, samples_workflows, template):
    '''
    (dict, dict, dict, dict, dict) -> None
            
    Adds in place the analysis information to the template
       
    Parameters
    ----------
    - block (dict): Dictionary with analyses organized in blocks of related workflows
    - workflows_information (dict): Information about each workflow for a given case
    - parent_workflows (dict): Dictionary with child-parents workflow relationship
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case 
    - template (dict): Dictionary collecting data for the given assay
    '''
    
    for sample_pair in block:
        for anchor_workflow in block[sample_pair]:
            for workflow_id in block[sample_pair][anchor_workflow]:
                name = workflows_information[workflow_id]
                # get the parent workflows
                inputs = parent_workflows[workflow_id]
                samples = [i for i in samples_workflows for j in samples_workflows[i]['workflows'] if workflow_id == j['wfrun_id']]
        
                d = {'workflow_id': workflow_id,
                     'workflow_name': name,
                     'samples': samples,
                     'inputs': inputs}
            
                if name in template['Analysis']:
                    # do not record identical workflows
                    if d not in template['Analysis'][name]:
                        template['Analysis'][name].append(d)



def get_alignment(template):
    '''
    (dict) -> dict
    
    Returns a dictionary with sample type with expected alignments for each sample in template
        
    Parameters
    ----------
    - template (str): Dictionary of template assay collecting data
    '''
        
    D = {}
    
    for i in template['Samples']:
        sample = template['Samples'][i]['sample']
        #assert sample not in D
        D[sample] = i + 'align'
    
    return D
    
    
def add_alignment_to_template(workflows_information, parent_workflows, samples_workflows, template):
    '''
    (dict, dict, dict, dict) -> None
    
    Adds in place the alignment data to the template
    
    Parameters
    ----------
    - workflows_information (dict): Information about each workflow for a given case
    - parent_workflows (dict): Dictionary with child-parents workflow relationship
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case 
    - template (dict): Dictionary collecting data for the given assay
    '''
    
    workflows = ['bcl2fastq', 'fileimportforanalysis', 'fileimport',
                       'import_fastq', 'bwamem']
    
    alignments = get_alignment(template)

    for workflow in template['Analysis']:
        for d in template['Analysis'][workflow]:
            workflow_id = d['workflow_id']
            # get the parent workflows
            inputs = parent_workflows[workflow_id]
            for wfrun_id in inputs:
                workflow_name = workflows_information[wfrun_id]
                if workflow_name.lower() in workflows:
                    sample = [i for i in samples_workflows for j in samples_workflows[i]['workflows'] if wfrun_id == j['wfrun_id']]
                    assert len(sample) == 1
                    sample = sample[0]
                    if sample in alignments:
                        alignment = alignments[sample]
                        b = {'workflow_id': wfrun_id,
                             'workflow_name': workflow_name,
                             'sample': sample}
                    
                        # add bwamem workflow information to the correct sample placeholder
                        if b not in template['Data'][alignment]['workflows']:
                                template['Data'][alignment]['workflows'].append(b)

                    
                    # for i in template['Data']:
                    #     if alignment in i.lower() and b not in template['Data'][i]['workflows']:
                    #         template['Data'][i]['workflows'].append(b)
                 
             

def add_anchors_to_template(assay, anchor_samples, template):
    '''
    (dict, dict, dict) -> None
            
    Adds in place the analysis information to the template
       
    Parameters
    ----------
    - assay (dict): Dictionary with assay information
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    - template (dict): Dictionary collecting data for the given assay
    '''
        
    for workflow in template['Analysis']:
        for sample in assay['Anchors']:
            if assay['Anchors'][sample]['workflows'] == workflow:
                for d in template['Analysis'][workflow]:
                    workflow_id = d['workflow_id']
                    if workflow_id in anchor_samples:
                        if sample in anchor_samples[workflow_id]:
                            for i in d['samples']:
                                if i in anchor_samples[workflow_id][sample]:
                                    if sample not in template['Anchors']:
                                        template['Anchors'][sample] = {}
                                    template['Anchors'][sample]['workflows'] = d['workflow_id']
    
    

def add_sequencing_to_template(workflows_information, parent_workflows, samples_workflows, template):
    '''
    (dict, dict, dict, dict) -> None
    
    Adds in place the sequencing data to the template
    
    Parameters
    ----------
    - workflows_information (dict): Information about each workflow for a given case
    - parent_workflows (dict): Dictionary with child-parents workflow relationship
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case 
    - template (dict): Dictionary collecting data for the given assay
    '''
    
    # make a list of bwamem workflows
    L = []
    for i in template['Data']:
        if i != 'Sequencing':
            for j in template['Data'][i]['workflows']:
                L.append(j['workflow_id'])
    L = list(set(L))
  
    # add sequencing data
    fastq_workflows = ['bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
    for wfrun_id in L:
        fastq_inputs = parent_workflows[wfrun_id]
        for i in fastq_inputs:
            if i != 'NA' and workflows_information[i].lower() in fastq_workflows:
                sample = [k for k in samples_workflows for m in samples_workflows[k]['workflows'] if i == m['wfrun_id']]
                assert len(sample) == 1
                sample = sample[0]
                f = {'workflow_id': i,
                     'workflow_name': workflows_information[i],
                     'sample': sample}
            elif i == 'NA' and workflows_information[wfrun_id].lower() in fastq_workflows:
                sample = [k for k in samples_workflows for m in samples_workflows[k]['workflows'] if wfrun_id == m['wfrun_id']]
                assert len(sample) == 1
                sample = sample[0]
                f = {'workflow_id': wfrun_id,
                     'workflow_name': workflows_information[wfrun_id],
                     'sample': sample}
                
            if f not in template['Data']['Sequencing']['workflows']:
                template['Data']['Sequencing']['workflows'].append(f)



def merge_group(groups):
    '''
    (dict) -> list
      
    Merge the blocks by sample and anchor
        
    Parameters
    ----------
    - groups (dict): Dictionary with analyses organized in blocks of related workflows
    '''
    
    # merge all dictionaries for each group  
    L = []
    for i in range(len(groups)):
        D = {}
        for j in range(len(groups[i])):
            for sample in groups[i][j]:
                assert sample not in D
                D[sample] = {}
                for anchor in groups[i][j][sample]:
                    assert anchor not in D[sample]
                    D[sample][anchor] = groups[i][j][sample][anchor]
        L.append(D)   
    
    return L

    

def count_workflows(blocks):
    '''
    (dict) -> dict
    
    Returns a dictionary with counts of each workflow in the analysis block
        
    Parameters
    ----------
    - block (dict): Dictionary with analyses organized in blocks of related workflows
    '''
    
    # identify workflows present in each subblocks
    counts = {}
    for sample in blocks:
        for anchor in blocks[sample]:
            for workflow in blocks[sample][anchor]:
                if workflow in counts:
                    counts[workflow] += 1
                else:
                   counts[workflow] = 1
    
    return counts



def clean_up_blocks(merged_blocks, workflow_information):
    '''
    (dict, dict) -> dict
    
    Removes in place the workflows from blocks that are parent/child workflows
    and processed with other workflows but that don't belong to the focus block
    
    Parameters
    ----------
    - merged_blocks (dict): Analysis blocks merged by sample and anchor
    - workflow_information (dict): Information about each workflow for a given donor
    '''
        
    for i in range(len(merged_blocks)):
        if len(merged_blocks[i]) > 1:
            # copy the block, because the same object may belong to multiple blocks
            # cleaning the blocks will remove unwanted workflows
            merged_blocks[i] = copy.deepcopy(merged_blocks[i])
            # count all the workflows
            counts = count_workflows(merged_blocks[i])
            # get the names of the workflows found across the sublocks
            names = []
            for j in counts:
                if counts[j] > 1:
                    names.append(workflow_information[j])
            # identify the workflows that should be found multiple times
            # but that are actually unique to a subblock. these are child and/or parent workflows
            # that are processed with other shared workflows but that do not belong to this block
            # these workflows should be removed
            to_remove = []
            for j in counts:
                if workflow_information[j] in names and counts[j] == 1:
                    to_remove.append(j)
            for j in to_remove:
                for sample in merged_blocks[i]:
                    for anchor in merged_blocks[i][sample]:
                        if j in merged_blocks[i][sample][anchor]:
                            merged_blocks[i][sample][anchor].remove(j)
                            # print('remove {0} with id {1} from block:'.format(workflow_information[j]['name'], j))
                            # print('sample {0}'.format(sample))
                            # print('anchor {0}'.format(anchor))
                            # print('---')



def get_library_type(sample, samples_workflows):
    '''
    (str, dict) -> str 
    
    Returns the library types associated with the sample(s)
        
    Parameters
    ----------
    - sample (str): Sample (or sample pair) associated with a block
    - samples_workflows (dict): Dictionary with workflow information for each sample of a donor 
    '''
            
    # sample may be a single sample or sample pair
    samples = list(map(lambda x: x.strip(), sample.split('|')))
    library_types = list(set([samples_workflows[i]['library_type'] for i in samples]))
    library_types = ';'.join(library_types)
    
    return library_types
    
   

def group_blocks(blocks, samples_workflows):
    '''
    (dict, dict) -> list
    
    Returns a dictionary with analyses organized in blocks of related workflows
        
    Parameters
    ----------
    - blocks (dict): Dictionary with workflows for a sample pair or sample
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case 
    '''
    
    # group blocks by library type of the samples
    T = {}
    for i in range(len(blocks)):
        sample = list(blocks[i].keys())
        assert len(sample) == 1
        sample = sample[0]
        # get the library types
        library_type = get_library_type(sample, samples_workflows)
        for anchor in list(blocks[i][sample].keys()):
            d = {}
            d[sample] = {}
            d[sample][anchor] = blocks[i][sample][anchor]
            if library_type in T:
                T[library_type].append(d)
            else:
                T[library_type] = [d]
            
    # create lists of lists of blocks for each library type
    L = [T[library_type] for library_type in T]
    # get all combinations of dictionaries to group the workflows by sample and anchor
    C = list(itertools.product(*L))
    
    return C
    
    
    
def fill_template(assay, block, samples_workflows, workflows_information, child_to_parents_workflows, anchor_samples):
    '''
    (dict, dict, dict, dict, dict, dict) -> dict
     
    Returns a template dictionary containing all the samples, sequencing, alignments,
    and analyses workflows for a case based on an expected assay
    
    Parameters
    ----------
    - assay (dict): Dictionary with assay information
    - block (dict): Dictionary with analyses organized in blocks of related workflows
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case 
    - workflows_information (dict): Information about each workflow for a given case
    - child_to_parents_workflows (dict): Dictionary with child-parents workflow relationship
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    '''
        
    template = convert_assay_to_template(assay)
    # add samples section to template
    add_samples_to_template(assay, block, samples_workflows, template)
    # add the analysis section
    add_analysis_to_template(block, workflows_information, child_to_parents_workflows, samples_workflows, template)
    # add alignment data
    add_alignment_to_template(workflows_information, child_to_parents_workflows, samples_workflows, template)
    # add sequencing data
    add_sequencing_to_template(workflows_information, child_to_parents_workflows, samples_workflows, template)
    # add anchor workflows
    add_anchors_to_template(assay, anchor_samples, template)

    return template



def get_sample_templates(merged_blocks, assay, samples_workflows, workflows_information, child_to_parents_workflows, anchor_samples):
    '''
    (dict, dict, dict, dict, dict) -> list
    
    Returns a list of templates with analysis data for each block
    
    Parameters
    ----------
    - merged_blocks (dict): Analysis blocks merged by sample and anchor
    - assay (dict): Dictionary with assay information
    - samples_workflows (dict): Dictionary with workflow information for each sample of a case 
    - workflows_information (dict): Information about each workflow for a given case
    - child_to_parents_workflows (dict): Dictionary with child-parents workflow relationship
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    '''

    # make a list of templates
    L = []
    for block in merged_blocks:
        L.append(fill_template(assay, block, samples_workflows, workflows_information, child_to_parents_workflows, anchor_samples))
        
    return L




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
        if len(template['Data'][i]['workflows']) == 0:
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
    error = ''

    # check if each expected workflow has a single record
    # skip call ready workflow because they may have a record each sample

    # get the anchow workflow names from the assay
    extra = 0
    anchor_workflows = extract_anchor_samples(assay)
    for workflow in template['Analysis']:
        if workflow not in anchor_workflows:
            if len(template['Analysis'][workflow]) > 1:
                complete = False
                extra += 1
    if complete == False:
        error = '{0} extra workflows'.format(extra)        

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
    
    # evaluate extra workflow occurences
    cpl, err = evaluate_extra_analyses(template, assay)
    complete.append(cpl)
    error.append(err)
        
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
                
     
        

     
        
     
        
     

     
        
     
        
     
        
     
        
     
        
     
        
     

def check_anchor_workflow(anchor_workflows, anchor_samples, workflow_info):
    '''
    (dict, dict, dict) -> bool
    
    Returns True if the workflows for each sample are anchor workflows
       
    Parameters
    ----------
    - anchor_workflows (dict): Dictionary with expected sample information for each anchor workflow
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    - workflow_info (dict): Information about each workflow for a given case
    '''
    
    L = [workflow_info[workflow]['name'] in anchor_workflows for workflow in anchor_samples]
       
    if L:
        return all(L)
    else:
        return False
    
   
    
def check_anchor_samples(anchor_workflows, anchor_samples):
   '''
   (dict, dict) -> bool

   Returns True if all the expected samples in the assay are defined
     
   Parameters
   ----------
   - anchor_workflows (dict): Dictionary with expected sample information for each anchor workflow
   - anchor_samples (dict): Dictionary with samples for each anchor worfklow
   '''

   expected_samples = []
   collected_samples = []
   for anchor_workflow in anchor_workflows:
       # make a list of sorted samples
       expected_samples.extend([i for i in anchor_workflows[anchor_workflow]])
   expected_samples = sorted(list(set(expected_samples)))
    

   for i in anchor_samples:
       for j in anchor_samples[i]:
           collected_samples.append(j)
   collected_samples = sorted(list(set(collected_samples)))
   

   samples = []
   for i in anchor_samples:
       for j in anchor_samples[i]:
           samples.extend(anchor_samples[i][j])
   samples.sort()     
   samples = list(set(samples))
   
   return expected_samples == collected_samples
        
    
    
    
def find_sample_groups_with_common_workflows(samples_workflows, group):
    '''
    (dict, list) -> dict
    
    Returns a dictionary of samples and workflows common to all samples in each group
    (ie, samples that have been processed together)
        
    Parameters
    ----------    
    - samples_workflows (dict): Dictionary of workflow information for each sample of a given case
    - group (list): List of samples
    '''

    D = {} 

    j = ' | '.join(sorted(group))
    # get the common workflows for each samples
    L = []
    for i in group:
        if i in samples_workflows:
            for d in samples_workflows[i]['workflows']:
                L.append(d)
    # keep the workflows common to all samples in the group
    K = []
    for d in L:
        if L.count(d) == len(group):
            if d not in K:
                K.append(d)
    # record the workflows for the sample group if they exist
    if K:
        D[j] = K

    return D
    
    
   
def identify_blocks(group, samples_workflows, parent_workflows, children_workflows, anchor_samples):
    '''    
    (list, dict, dict, dict, dict) -> list
    
    Returns a list of common workflows for paired samples or worklfows for single sample
       
    Parameters
    ----------
    - group (list): Samples processed together in an assay (single sample or sample pair)
    - samples_workflows (dict): Dictionary of workflow information for each sample of a given case
    - parent_workflows (dict): Dictionary with child-parents workflow relationship
    - children_workflows (dict): Dictionary with list of downstream workflow ids for each workflow id
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    '''
    
    L = []
    
    sample_group = []
    for i in range(len(group)):
        if type(group[i]) != list:
            sample_group.append([group[i]])
        else:
            sample_group.append(group[i])
    
    for samples in sample_group:
       # check if paired sample or single sample
       D = find_sample_groups_with_common_workflows(samples_workflows, samples)
       if D:
           blocks = identify_analysis_workflows(D, parent_workflows, children_workflows, list(anchor_samples.keys()))
           L.append(blocks)
      
    return L


def get_sample_groups(groups, samples_workflows, child_to_parents_workflows, parent_to_children_workflows, anchor_samples):
    '''
    (list, dict, dict, dict, dict) -> list
    
    Returns a list of analysis blocks for each sample pair/sample in the sample groups
    
    Parameters
    ----------
    - groups (list): List of samples processed together in an assay (single sample or sample pair)
    - samples_workflows (dict): Dictionary of workflow information for each sample of a given case
    - child_to_parents_workflows (dict): Dictionary with child-parents workflow relationship
    - parent_to_children_workflows (dict): Dictionary with list of downstream workflow ids for each workflow id
    - anchor_samples (dict): Dictionary with samples for each anchor worfklow
    '''

    L = []
    # for each group of samples, find common workflows for paired samples and/or workflows for single samples
    for group in groups:
        blocks = identify_blocks(group, samples_workflows, child_to_parents_workflows, parent_to_children_workflows, anchor_samples)
        if blocks:
            for d in blocks:
                if d and d not in L:
                    L.append(d)
    return L




def get_assay_workflows(assays):
    '''
    (dict) -> dict
    
    Returns a dictionary with all the workflows expected from the assay
    
    Parameters
    ----------
    - assays (dict): Dictionnary with all the assay types
    '''

    D = {}

    for assay in assays:
        anchors = list(set([assays[assay]['Anchors'][i]['workflows'] for i in assays[assay]['Anchors']]))
        workflows = list(set([i for i in assays[assay]['Analysis']]))

        D[assay] = {'anchors': anchors, 'workflows': workflows}

    return D



def get_case_workflows(workflow_info):
    '''
    (dict) -> set

    Returns a set of workflows for a given case    
    
    Parameters
    ----------
    - workflow_info (dict): Dictionary with workflow information for a given case
    '''
    
    S = set([workflow_info[i]['name'] for i in workflow_info])
    
    return S


def select_assay(assays, workflow_info):
    '''
    (dict, dict) -> str
    
    Returns the assay name matching the observed analysis workflows
    
    Parameters
    ----------
    - assays (dict): Dictionary with assays
    - workflow_info (dict): Dictionary with workflow information for a given case
    '''

    assay_wf = get_assay_workflows(assays)
    case_wf = get_case_workflows(workflow_info)

    assay = ''

    for i in assay_wf:
        if (set(assay_wf[i]['anchors']).intersection(case_wf) == set(assay_wf[i]['anchors'])) and \
            (set(assay_wf[i]['workflows']).intersection(case_wf) == set(assay_wf[i]['workflows'])):
                assay = i

    return assay




# assays = {
#   "WG.1":{
#     "Samples":{
#         "TumorWG":{"library_type":"WG",
#             "tissue_type":["P","M","X"]},
#         "MatchedNormalWG":{"library_type":"WG",
#                             "tissue_type":"R"}
#     },
#     "Data":{
#         "Sequencing":{"workflows":["bcl2fastq"],"inputs":[]},
#         "TumorWGalign":{"workflows":["bwaMem"],"samples":["TumorWG"],"inputs":[],"data":["Sequencing"]},
#         "MatchedNormalWGalign":{"workflows":["bwaMem"],"samples":["MatchedNormalWG"],"inputs":[],"data":["Sequencing"]}
#     },
#     "Analysis":{
#         "bamMergePreprocessing": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bwaMem"],"data":["wgalignT","wgalignN"]},
#         "mutect2_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing"]},
#         "sequenza": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["varscan"]},
#         "variantEffectPredictor_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["mutect2_matched"]},
#         "delly_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing"]},
#         "varscan": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing"]},
#         "strelka2_somatic": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing"]},
#         "variantMerging": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["strelka2_somatic"]},
#         "mavis": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["delly_matched", "bamMergePreprocessing"]}
#     },
#     "Anchors":{"TumorWG": {"workflows": "bamMergePreprocessing"},
#                "MatchedNormalWG": {"workflows": "bamMergePreprocessing"}}
#   },
    
#   "WG.2":{
#     "Samples":{
#         "TumorWG":{"library_type":"WG",
#             "tissue_type":["P","M","X"]},
#         "MatchedNormalWG":{"library_type":"WG",
#                             "tissue_type":"R"}
#     },
#     "Data":{
#         "Sequencing":{"workflows":["bcl2fastq"],"inputs":[]},
#         "TumorWGalign":{"workflows":["bwaMem"],"samples":["TumorWG"],"inputs":[],"data":["Sequencing"]},
#         "MatchedNormalWGalign":{"workflows":["bwaMem"],"samples":["MatchedNormalWG"],"inputs":[],"data":["Sequencing"]}
#     },
#     "Analysis":{
#         "bamMergePreprocessing_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bwaMem"],"data":["wgalignT","wgalignN"]},
#         "mutect2_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "sequenza_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["varscan_by_tumor_group"]},
#         "variantEffectPredictor_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["mutect2_matched_by_tumor_group"]},
#         "delly_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "varscan_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "strelka2_somatic_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "mavis": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group", "star_call_ready", "starfusion", "delly_matched_by_tumor_group"]},
#         "msisensor": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]}
#     },
#     "Anchors":{"TumorWG": {"workflows": "bamMergePreprocessing_by_tumor_group"},
#                "MatchedNormalWG": {"workflows": "bamMergePreprocessing_by_tumor_group"}}
#   },
  
  
#   "WG.3":{
#       "Samples":{
#           "TumorWG":{"library_type":"WG",
#               "tissue_type":["P","M","X"]},
#           "MatchedNormalWG":{"library_type":"WG",
#                               "tissue_type":"R"}
#       },
#       "Data":{
#           "Sequencing":{"workflows":["bcl2fastq"],"inputs":[]},
#           "TumorWGalign":{"workflows":["bwaMem"],"samples":["TumorWG"],"inputs":[],"data":["Sequencing"]},
#           "MatchedNormalWGalign":{"workflows":["bwaMem"],"samples":["MatchedNormalWG"],"inputs":[],"data":["Sequencing"]}
#       },
#       "Analysis":{
#           "bamMergePreprocessing_by_sample": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bwaMem"],"data":["wgalignT","wgalignN"]},
#           "mutect2_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#           "sequenza": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["varscan"]},
#           "variantEffectPredictor_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["mutect2_matched"]},
#           "delly_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#           "varscan": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#           "strelka2_somatic":{"workflows":["strelka2_somatic"],"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#           "variantMerging":{"workflows":["variantMerging"],"samples":["TumorWG","MatchedNormalWG"],"inputs":["strelka2_somatic", "mutect2_matched"]},
#           "mavis":{"workflows":["mavis"],"samples":["TumorWG","MatchedNormalWG"],"inputs":["delly_matched", "bamMergePreprocessing_by_sample"]}
#       },
#       "Anchors":{"TumorWG": {"workflows": "bamMergePreprocessing_by_sample"},
#                  "MatchedNormalWG": {"workflows": "bamMergePreprocessing_by_sample"}}
#     },

#   "WGTS.1":{
#     "Samples":{
#         "TumorWG":{
#             "library_type":"WG",
#             "tissue_type":["P","M","X"]
#         },
#         "MatchedNormalWG":{
#             "library_type":"WG",
#             "tissue_type":"R"
#         },
#         "TumorWT":{
#             "library_type":"WT",
#             "tissue_type":["P","M","X"]
#         }
#     },
#     "Data":{
#         "Sequencing":{"workflows":["bcl2fastq"],"inputs":[]},
#         "TumorWGalign":{"workflows":["bwaMem"],"samples":["TumorWG"],"inputs":[],"data":["Sequencing"]},
#         "MatchedNormalWGalign":{"workflows":["bwaMem"],"samples":["MatchedNormalWG"],"inputs":[],"data":["Sequencing"]},
#         "TumorWTalign":{"workflows":["star"],"samples":["TumorWT"],"inputs":[],"data":["Sequencing"]}
#     },
#     "Analysis":{
#         "bamMergePreprocessing_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bwaMem"],"data":["wgalignT","wgalignN"]},
#         "star_call_ready": {"samples":["TumorWT"],"inputs":["bcl2fastq"],"data":["wtalignT"]},
#         "mutect2_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "variantEffectPredictor_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["mutect2_matched_by_tumor_group"]},
#         "varscan_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "sequenza_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["varscan_by_tumor_group"]},
#         "delly_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "mavis": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group", "star_call_ready", "starfusion", "delly_matched_by_tumor_group"]},
#         "rsem": {"samples":["TumorWT"],"inputs":["star_call_ready"]},
#         "starfusion":{"samples":["TumorWT"],"inputs":["star_call_ready"]}
#     },
#     "Anchors":{"TumorWG": {"workflows": "bamMergePreprocessing_by_tumor_group"},
#                "MatchedNormalWG": {"workflows": "bamMergePreprocessing_by_tumor_group"},
#                "TumorWT": {"workflows": "star_call_ready"}}
#     },

#   "WGTS.2":{
#     "Samples":{
#         "TumorWG":{
#             "library_type":"WG",
#             "tissue_type":["P","M","X"]
#         },
#         "MatchedNormalWG":{
#             "library_type":"WG",
#             "tissue_type":"R"
#         },
#         "TumorWT":{
#             "library_type":"WT",
#             "tissue_type":["P","M","X"]
#         }
#     },
#     "Data":{
#         "Sequencing":{"workflows":["bcl2fastq"],"inputs":[]},
#         "TumorWGalign":{"workflows":["bwaMem"],"samples":["TumorWG"],"inputs":[],"data":["Sequencing"]},
#         "MatchedNormalWGalign":{"workflows":["bwaMem"],"samples":["MatchedNormalWG"],"inputs":[],"data":["Sequencing"]},
#         "TumorWTalign":{"workflows":["star"],"samples":["TumorWT"],"inputs":[],"data":["Sequencing"]}
#     },
#     "Analysis":{
#         "bamMergePreprocessing_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bwaMem"],"data":["wgalignT","wgalignN"]},
#         "star_call_ready": {"samples":["TumorWT"],"inputs":["bcl2fastq"],"data":["wtalignT"]},
#         "mutect2_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "variantEffectPredictor_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["mutect2_matched_by_tumor_group"]},
#         "varscan_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "sequenza_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["varscan_by_tumor_group"]},
#         "delly_matched_by_tumor_group": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_tumor_group"]},
#         "mavis": {"samples":["TumorWG","MatchedNormalWG","TumorWT"],"inputs":["bamMergePreprocessing_by_tumor_group", "star_call_ready", "starfusion", "delly_matched_by_tumor_group"]},
#         "rsem": {"samples":["TumorWT"],"inputs":["star_call_ready"]},
#         "starfusion":{"samples":["TumorWT"],"inputs":["star_call_ready"]},
#         "arriba":{"samples":["TumorWT"],"inputs":["star_call_ready"]}
#     },
#     "Anchors":{"TumorWG": {"workflows": "bamMergePreprocessing_by_tumor_group"},
#                "MatchedNormalWG": {"workflows": "bamMergePreprocessing_by_tumor_group"},
#                "TumorWT": {"workflows": "star_call_ready"}}
        
        
#   }
# }





def get_column_names(database, table):
    '''
    (str, str) -> list
    
    Returns a list of column headers from the table in database
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - table (str): Table in database
    '''
        
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    cur.execute('select * from {0}'.format(table))
    column_names = [i[0] for i in cur.description]
    conn.close()

    return column_names



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
                       'haplotypeCaller': {'samples': ['NormalWG'],
                                           'inputs': ['bwaMem']}},
                   'Anchors': {'TumourWT': {'workflows': 'star_call_ready'},
                               'TumourWG': {'workflows': 'bamMergePreprocessing_by_sample'},
                               'NormalWG': {'workflows': 'bamMergePreprocessing_by_sample'}}}



# hrdetect	Normal WG,Tumour WG	sample	calls.hrd	ALL
# haplotypeCaller	Normal WG	sample	calls.germline.mutations	ALL        
        
        
        
        

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
                       'mavis': {'samples': ['TumourWT', 'NormalWG', 'TumourWG'],
                                 'inputs': ['bwaMem', 'delly_matched']}},
                       'hrDetect': {'samples': ['NormalWG', 'TumourWG'],
                                    'inputs': ['purple', 'mutect2_matched']},
                       'haplotypeCaller': {'samples': ['NormalWG'],
                                           'inputs': ['bwaMem']},
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
                    # find the children of each workflow
                    parent_to_children_workflows = collect_workflow_relationships(case_data)
                    # find the parents of each workflow
                    child_to_parents_workflows = get_downstream_workflows(parent_to_children_workflows)
                    anchor_samples = get_anchor_samples(samples_workflows, anchor_workflows)
                    # make groups of samples (can be 1, 2 or more samples depending on the assay)       
                    groups = group_samples(anchor_samples, assay)
                    sample_groups = get_sample_groups(groups, samples_workflows, child_to_parents_workflows, parent_to_children_workflows, anchor_samples)
                    blocks = group_blocks(sample_groups, samples_workflows)    
                    merged_blocks = merge_group(blocks)
                    clean_up_blocks(merged_blocks, workflow_info)
                    templates = get_sample_templates(merged_blocks, assay, samples_workflows, workflow_info, child_to_parents_workflows, anchor_samples)
                
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



# assays = {
#   "RUO WGS - 80XT/30XN_1.0":{
#     "Samples":{
#         "TumorWG":{"library_type":"WG",
#             "tissue_type":["P","M","X", "n"]},
#         "MatchedNormalWG":{"library_type":"WG",
#                             "tissue_type":"R"}
#     },
#     "Data":{
#         "Sequencing":{"workflows":["bcl2fastq"],"inputs":[]},
#         "TumorWGalign":{"workflows":["bwaMem"],"samples":["TumorWG"],"inputs":[],"data":["Sequencing"]},
#         "MatchedNormalWGalign":{"workflows":["bwaMem"],"samples":["MatchedNormalWG"],"inputs":[],"data":["Sequencing"]}
#     },
#     "Analysis":{
#         "bamMergePreprocessing_by_sample": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bwaMem"],"data":["wgalignT","wgalignN"]}, 
#         "strelka2_somatic": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "variantEffectPredictor_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["mutect2_matched"]},
#         "mutect2_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "delly_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "varscan": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "sequenza": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["varscan"]},
#         "msisensor": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "gridss_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "purple": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample","mutect2_matched", "gridss_matched"]},
#         "hrDetect": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["purple","mutect2_matched"]}
#     },
#     "Anchors":{"TumorWG": {"workflows": "bamMergePreprocessing_by_sample"},
#                "MatchedNormalWG": {"workflows": "bamMergePreprocessing_by_sample"}}
#   },
       
#   "RUO WGTS - 80XT/30XN_1.0":{
#     "Samples":{
#         "TumorWG":{
#             "library_type":"WG",
#             "tissue_type":["P","M","X","n"]
#         },
#         "MatchedNormalWG":{
#             "library_type":"WG",
#             "tissue_type":"R"
#         },
#         "TumorWT":{
#             "library_type":"WT",
#             "tissue_type":["P","M","X","n"]
#         }
#     },
#     "Data":{
#         "Sequencing":{"workflows":["bcl2fastq"],"inputs":[]},
#         "TumorWGalign":{"workflows":["bwaMem"],"samples":["TumorWG"],"inputs":[],"data":["Sequencing"]},
#         "MatchedNormalWGalign":{"workflows":["bwaMem"],"samples":["MatchedNormalWG"],"inputs":[],"data":["Sequencing"]},
#         "TumorWTalign":{"workflows":["star_lane_level"],"samples":["TumorWT"],"inputs":[],"data":["Sequencing"]}
#     },
#     "Analysis":{
#         "bamMergePreprocessing_by_sample": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bwaMem"],"data":["wgalignT","wgalignN"]}, 
#         "star_call_ready": {"samples":["TumorWT"],"inputs":["bcl2fastq"],"data":["wtalignT"]},
#         "rsem": {"samples":["TumorWT"],"inputs":["star_call_ready"]},
#         "starfusion":{"samples":["TumorWT"],"inputs":["star_call_ready"]},
#         "strelka2_somatic": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "msisensor": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "varscan": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "sequenza": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["varscan"]},
#         "gridss_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "mutect2_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "delly_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["bamMergePreprocessing_by_sample"]},
#         "arriba": {"samples":["TumorWT"],"inputs":["star_call_ready"]},
#         "variantEffectPredictor_matched": {"samples":["TumorWG","MatchedNormalWG"],"inputs":["mutect2_matched"]},
#         "mavis": {"samples":["TumorWG","MatchedNormalWG", "TumorWT"],"inputs":["star_call_ready",
#                                                                                "delly_matched",
#                                                                                "bamMergePreprocessing_by_sample",
#                                                                                "arriba", "starfusion"]}
#   },
#     "Anchors":{"TumorWG": {"workflows": "bamMergePreprocessing_by_sample"},
#                "MatchedNormalWG": {"workflows": "bamMergePreprocessing_by_sample"},
#                "TumorWT": {"workflows": "star_call_ready"}}
#     }
# }



#generate_cache('provenance_reporter.json', 'enabled_workflows.json', 'waterzooi_db_case_new.db', 'http://pinery.gsi.oicr.on.ca/assays', 'analysis_review_case_new.db', table='templates')

#generate_cache('provenance_reporter.json', 'enabled_workflows.json', 'waterzooi_db_case_new.db', 'http://pinery.gsi.oicr.on.ca/assays', 'analysis_review_case_MOH.db', table='templates')

generate_cache('provenance_reporter.json', 'enabled_workflows.json', 'waterzooi_db_case_MOH.db', 'http://pinery.gsi.oicr.on.ca/assays', 'analysis_review_case_MOH.db', table='templates')



# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(prog = 'analysis_review.py', description='Script to generate the analyais review cache')
#     parser.add_argument('-p', '--provenance', dest = 'provenance', help = 'Path to the provenance reporter data json', required=True)
#     parser.add_argument('-db', '--database', dest='database', help='Path to the analysis review database', required=True)    
#     parser.add_argument('-cq', '--calcontaqc', dest = 'calcontaqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest', help = 'Path to the merged rnaseq calculateContamination database. Default is /scratch2/groups/gsi/production/qcetl_v1/calculatecontamination/latest')


    
#     # get arguments from the command line
#     args = parser.parse_args()
#     # generate sqlite cache
#     generate_cache(args.provenance, args.database)    
    