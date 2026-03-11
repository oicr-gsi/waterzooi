# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 16:01:53 2025

@author: rjovelin
"""

import json
import requests
import sqlite3




def extract_assays(pinery):
    '''
    (str) -> list
    
    Returns a list of assays from Pinery
    
    Parameters
    ----------
    - pinery (str): URL to Pinery assay endpoint
    '''
    
    # get the expected samples from the assays in Pinery
    headers = {'accept': 'application/json',}
    response = requests.get(pinery, headers=headers)
    
    if response.status_code == 200:
        L = response.json()
    else:
        L = []
    
    return L
    


def extract_assay_workflows(assay_config_file):
    '''
    (str) -> dict   
    
    Returns a dictionary with expected workflows for each assay in the assay config
    
    Parameters
    ----------
    - assay_config_file (str): Path to the assay config file
    '''
    
    infile = open(assay_config_file)
    workflows = json.load(infile)
    infile.close()
    
    D = {}
    
    for assay in workflows['values']:
        for version in workflows['values'][assay]['versions']:
            if assay not in D:
                D[assay] = {}
            assert version not in D[assay]
            D[assay][version] = list(workflows['values'][assay]['versions'][version]['workflows'].keys())
    
    return D



def record_assay_samples(assay):
    '''
    (dict) -> dict
    
    Returns a dictionary with the sample information of a given assay
    
    Parameters
    ----------
    - assay (dict): Dictionary with information for a specific assay pulled from Pinery 
    '''
    
    D = {}
    
    if 'tests' in assay:
        # check the assay samples:
        for d in assay['tests']:
            library_source = d['library_source_template_type']
            if 'tissue_type' not in d:
                ## assume this is a tumor only assay
                tissue_type = 'R'
                negate_tissue_type = True
            else:
                tissue_type = d['tissue_type']
                negate_tissue_type = d['negate_tissue_type']
            if tissue_type == 'R' and negate_tissue_type == False:
                # normal sample
                sample_name = 'Normal' + library_source 
                sample_tissue = tissue_type
            elif tissue_type == 'R' and negate_tissue_type:
                # tumour sample
                sample_name = 'Tumour' + library_source
                sample_tissue = 'R'
            # record sample information   
            D[sample_name] = {'library_type': library_source,
                              'tissue_type': sample_tissue,
                              'negate_tissue_type': negate_tissue_type}
    
    return D    
    
    
    
    
def list_qc_workflows(assay_configs):
    '''
    (dict) -> list
    
    Returns a list of QC workflows 
    
    Parameters
    ----------
    - assay_configs (dict): Dictionary with the expected workflows for each assay
    '''

    L = []
    
    for assay in assay_configs:
        for version in assay_configs[assay]:
            for workflow in assay_configs[assay][version]:
                if 'qc' in workflow.lower():
                    L.append(workflow)
                elif 'callability' in workflow:
                    L.append(workflow)
                elif 'metrics' in workflow.lower():
                    L.append(workflow)
                elif 'contamination' in workflow.lower():
                    L.append(workflow)           
                elif 'collector' in workflow.lower():
                    L.append(workflow)
                elif 'bcl2barcode' in workflow.lower():
                    L.append(workflow)
                elif 'tmbanalysis' in workflow.lower():
                    L.append(workflow)
                    
    return L


            
      

def is_sequencing_workflow(workflow):
    '''
    (str) -> bool
    
    Returns True if workflow is a fastq generating workflow
    
    Parameters
    ----------
    - workflow (str): Workflow name
    '''
    
    sequencing_workflows = ['bcl2fastq', 'casava', 'fileimport', 'fileimportforanalysis', 'import_fastq']

    return any([workflow == i or i in workflow for i in sequencing_workflows])

    
    
def record_sequencing_assay_data(assay_config, qc_workflows):
    '''
    (list, list) -> list
    
    Returns a list of sequencing workflows from the assay
    
    Parameters
    ----------
    - assay_config (list): List of workflows for a given assay
    - qc_workflows (list): List of QC workflows across all assays
    '''
  
    L = []  
  
    for workflow in assay_config:
        if is_sequencing_workflow(workflow):
            assert workflow not in qc_workflows
            L.append(workflow)
                
    return L
    


def is_lane_level_workflow(workflow, library_type):
    '''
    (str, str) -> bool
    
    Returns True if workflow is a lane level alignment workflow
    for the specific library type
    
    Parameters
    ----------
    - workflow (str): Workflow name
    - library_type (str): Library design 
    '''

    # star is the lane level aligner for transcriptomics
    if library_type == 'WT':
        return 'star_lane_level' in workflow.lower()
    else:
        # other library types use versions of bwamem
        return 'bwa' in workflow.lower()
      
   
    
def record_alignment_assay_data(samples, assay_config):
    '''
    (dict, list) -> dict
    
    Returns a dictionary with the lane level alignments information of an assay
    
    Parameters
    ----------
    - samples (dict): Dictionary with the sample information of a given assay
    - assay_config (list): List of workflows for a given assay
    '''    

    D = {}

    for sample in samples:
        k = sample + 'align'
        library_type = samples[sample]['library_type']
        D[k] = {"samples":[sample],"inputs":[],"data":["Sequencing"]}
        D[k]["workflows"] = [workflow for workflow in assay_config if is_lane_level_workflow(workflow, library_type)]
                
    return D


def record_anchors_assay_data(samples, assay_config):
    '''
    (dict, list) -> dict
    
    Returns a dictionary with the expected call ready workflows (anchor workflows) for each sample
    
    Parameters
    ----------
    - samples (dict): Dictionary with the sample information of a given assay
    - assay_config (list): List of workflows for a given assay
    '''

    D = {}
    
    for sample in samples:
        if samples[sample]['library_type'] == 'WT':
            call_ready_workflow = [workflow for workflow in assay_config if 'star_call_ready' in workflow.lower()]
        else:
            call_ready_workflow = [workflow for workflow in assay_config if 'bammergepreprocessing' in workflow.lower()]
        call_ready_workflow = list(set(call_ready_workflow))
        if call_ready_workflow:
            assert len(call_ready_workflow) == 1
            call_ready_workflow = call_ready_workflow[0]
        else:
            return {}
        
        D[sample] = {"workflows": call_ready_workflow}           

    return D            




def get_workflow_relationships(database):
    '''
    (str) -> dict
    
    Returns a dictionary with child and parent workflows
    
    Parameters
    ----------
    - database (str): Path to the waterzooi database
    '''
    
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    data = conn.execute("SELECT DISTINCT parents_id, children_id FROM Parents").fetchall() 
    conn.close()

    D = {}
    
    for i in data:
        parent = i['parents_id']
        child = i['children_id']
        if child in D:
            D[child].append(parent)
        else:
            D[child] = [parent]

    return D



def get_workflow_tissue_types(database):
    '''
    (str) -> dict
    
    Returns a dictionary with the sample inputs (library and tissue type) of each workflow
    
    Parameters
    ----------
    - database (str): Path to the waterzooi database
    '''
    
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    data = conn.execute("SELECT Workflows.wfrun_id, Libraries.library_type, Libraries.tissue_type FROM Workflows JOIN Workflow_Inputs \
                        JOIN Libraries WHERE Workflows.wfrun_id = Workflow_Inputs.wfrun_id \
                        AND Workflow_Inputs.library = Libraries.library").fetchall() 
    conn.close()

    D = {}
    
    for i in data:
        if i['wfrun_id'] not  in D:
            D[i['wfrun_id']] = []
        if i['tissue_type'] == 'R':
            d = {'library_type': i['library_type'], 'negate_tissue_type': False}
        else:
            d = {'library_type': i['library_type'], 'negate_tissue_type': True}
        if d not in D[i['wfrun_id']]:
            D[i['wfrun_id']].append(d)
        
    return D




def get_workflow_names(database):
    '''
    (str) -> dict
    
    
    Parameters
    ----------
    - database (str): Path to the waterzooi database
    '''

    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    data = conn.execute("SELECT wfrun_id, wf FROM Workflows").fetchall() 
    conn.close()

    D = {}
    
    for i in data:
        if i['wfrun_id'] in D:
            assert D[i['wfrun_id']] == i['wf']
        else:
            D[i['wfrun_id']] = i['wf']
            
    return D



def is_analysis_worfklow(workflow):
    '''
    (str) -> bool
    
    Returns True if the workflow is an expected analysis workflow (ie, not a QC workflow)
    
    Parameters
    ----------
    - workflow (str): Name of the workflow
    '''
        
    if any(['qc' in workflow.lower(), 'callability' in workflow.lower(),
            'metrics' in workflow.lower(), 'contamination' in workflow.lower(),
            'collector' in workflow.lower(), 'bcl2barcode' in workflow.lower(),
            'tmbanalysis' in workflow.lower()]):
        return False
    else:
        return True
            

def get_workflow_inputs(database):
    '''
    (str) -> dict
    
    Returns a dictionary with the workflow inputs (library types and workflow parents)
    for each worklow and its version extracted from the production data in the database
       
    Parameters
    ----------
    - database (str): Path to the waterzooi database
    '''
    
    workflow_relationships = get_workflow_relationships(database)
    workflow_names = get_workflow_names(database)
    workflow_library_types = get_workflow_tissue_types(database)    
        
    D = {}
    
    for child in workflow_relationships:
        if child != 'NA':
            name = workflow_names[child]
            if is_analysis_worfklow(name):
                # get the parents workflows
                parents_m = [workflow_names[i] for i in workflow_relationships[child] if i != 'NA' and is_analysis_worfklow(workflow_names[i])]
                # eliminate doublons
                seen = []
                parents = sorted([i for i in parents_m if i not in seen and not seen.append(i)])
                library_type = sorted(workflow_library_types[child], key=lambda x: (x['library_type'], x['negate_tissue_type']))
                d = {'library_types': library_type, 'parents': parents}  
            
                if name not in D:
                   D[name] = []
                if d not in D[name]:
                   D[name].append(d)
    
    return D




def parent_workflows_in_assay(assay_config, parents):
    '''
    (dict, list) -> bool
    
    Returns True is all the parent worklows (workflow and version) 
    of a given workflow are defined in the assay
        
    Parameters
    ----------
    - assay_config (dict): Dictionary with expected workflows for an assay
    - parents (list): List of dictionaries with parent workflow name and version
    '''
    
    L = []
    
    # check that parent workflows are in the assay
    for i in parents:
        L.append(i in assay_config)
    
    
    return all(L)
    
    


def record_analysis_assay_data(assay_config):
    '''
    (list)     
    
    Returns a dictionary with each analysis workflow from assay_config
    
    Parameters
    ----------
    - assay_config (list): Dictionary with expected assay workflows    
    '''
    
    D = {}
    
    for workflow in assay_config:
        if is_analysis_worfklow(workflow) and is_sequencing_workflow(workflow) == False:
            if workflow not in D:
                D[workflow] = {'samples': [], 'inputs': []}    
            
    return D 


def create_assay_template(assay_name, assay_version, qc_workflows, assay_config, assay):
    '''
    (str, str, list, list, dict) -> dict
    
    Returns a template representing the expected data for a given assay
        
    Paramaters
    ----------
    - assay_name (str): Name of the assay
    - assay_version (str): Version of the assay
    - qc_workflows (list): List of QC workflows
    - assay_config (list): List with expected assay workflows    
    - assay (dict): Dictionary with assay information from Pinery 
    '''
       
    D = {'Samples':{},
         'Data': {},
         'Analysis': {},
         'Anchors': {}}
    
    
    # record the samples
    samples = record_assay_samples(assay)
    D['Samples'] = samples
    # record the sequencing data
    sequencing = record_sequencing_assay_data(assay_config, qc_workflows)
    D["Data"]["Sequencing"] = {"workflows":sequencing,"inputs":[]}
    # record lane level data
    D["Data"].update(record_alignment_assay_data(samples, assay_config))
    # record anchor workflows (call read workflows)
    D["Anchors"] = record_anchors_assay_data(samples, assay_config)
    # record analysis workflows
    D['Analysis'] = record_analysis_assay_data(assay_config)
    
    for i in D:
        if len(D[i]) == 0:
            print(i, 'is empty')
            return {}
    
    return D



def generate_templates(assay_configurations, qc_workflows, pinery='http://pinery.gsi.oicr.on.ca/assays'):
    '''
    (dict, list, str) -> dict


    Returns a dictionary with the templates of all assays defined in pinery
    provided these templates do not have missing information and that the assay are defined in the 
    assay_config_file

    Parameters
    ----------
    - assay_configurations (dict): Dictionary with workflows for each assay and version
    - qc_workflows (list): List of QC workflkows
    - pinery (str): URL to Pinery assay endpoint
    '''

    assays = extract_assays(pinery)
    print('extracted assays')
    

    D = {}
    
    for assay in assays:
        assay_name = assay['name']
        assay_version = assay['version']
        print(assay_name, assay_version)
        if assay_name in assay_configurations:
            if assay_version in assay_configurations[assay_name]:
                assay_config = assay_configurations[assay_name][assay_version]
                template = create_assay_template(assay_name, assay_version, qc_workflows, assay_config, assay)
                
            else:
                print('version {0} not in assay_config'.format(assay_version))
                template = {}
        else:
            print('name {0} not in assay_config'.format(assay_name))
            template = {}
        
        if template:
            name = assay_name + '_' + assay_version
            assert name not in D
            D[name] = template
        
        print('-----')
        
            
    return D

