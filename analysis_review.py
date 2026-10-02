# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 16:36:23 2024

@author: rjovelin
"""

import json
import argparse
import os
import itertools


from commons import load_data, is_case_info_complete, get_cases_md5sum, find_sequencing_attributes, get_donor_name, \
    compute_md5, case_to_update, connect_to_db, define_columns, initiate_db, insert_multiple_records, \
    delete_unique_record, delete_multiple_records
        

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


def is_signoff_complete(case_data):
    '''
    (dict) -> bool

    Returns True if signoff is complete for all the lims Ids that pass QC
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''

    
    # evaluate all lims ids --> signoff is indicated by complete sequencing status   
    seq = json.loads(case_data['case_info']['sequencing'])
    # evaluate only lims ids for sequencing
    complete = []
    for i in seq:
        if i['type'] == 'FULL_DEPTH_SEQUENCING':
            complete.append(i['complete'])
    
    return all(complete)        


def map_lims_to_tests(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all sequencing lims ids passing QC for
    each tests in the case. 
    Assumption: sign off is complete (ie. sequencing status complete for all limds)
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''
    
    D = {}
        
    seq = json.loads(case_data['case_info']['sequencing'])
    # evaluate only lims ids for sequencing
    for i in seq:
        if i['type'] == 'FULL_DEPTH_SEQUENCING':
            # assumes signoff is complete
            assert i['complete']
            test = i['test']
            for j in i['limsIds']:
                if j['qcFailed'] == False:
                    limsid = j['id']
                    if test not in D:
                        D[test] = [limsid]
                    else:
                        D[test].append(limsid)
                        D[test].sort()
    return D


def map_tests_to_samples(case_data, tests):
    '''
    (dict, dict) -> dict

    Returns a dictionary matching all samples for each test
    Assumption: sign off is complete (ie. sequencing status complete for all limds)
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    - tests (dict): Dictionary with lims ids for each tests
    '''

    D = {}

    for test in tests:
        # find the sample corresponding to each lims id
        for i in case_data['sample_info']:
            limsid = i['limsId']
            sampleid = i['sampleId']
            if limsid in tests[test]:
                if test not in D:
                    D[test] = [sampleid]
                else:
                    D[test].append(sampleid)
                D[test] = sorted(list(set(D[test]))) 

    return D



def map_workflows_to_lims(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with lims ids for each workflow in a case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''
        
    D = {}
        
    for d in case_data['workflow_runs']:
        wfrunid = d['wfrunid']
        limsids = d['limsIds'].split(',')
        if wfrunid not in D:
            D[wfrunid] = limsids
        else:
            D[wfrunid].extend(limsids)
        D[wfrunid] = list(set(D[wfrunid]))

    return D




def map_expected_workflows_to_runid(workflow_info, case_workflows):
    '''
    (dict, list) -> bool
    
    Returns a dictionary with the workflow run ids, if they exist, for each expected
    workflow from the assay of a case
        
    Parameters
    ----------
    - workflow_info (dict): Dictionary mapping worfklow run ids to workflow names
    - case_workflows (list): List of expected workflows from the assay
    '''
    
    
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
    gridss_workflows = ['gridss_matched', 'gridss'] 
    bwa_workflows = ['bwaMem', 'bwamem2']  
       
    
    # map the workflow run ids of production workflows to the expected workflows
    D = {}
    for workflow in case_workflows:
        # collect run ids of all expected workflows
        D[workflow] = []
        for wfrunid in workflow_info:
            # sequencing workflows may be defined as bcl2fastq in the assay but the actual
            # sequencing workflow may differ
            # check if a run id exists for an alternative sequencing workflow
            if workflow in sequencing_workflows:
                for key in sequencing_workflows:
                    if key == workflow_info[wfrunid]:
                        if key not in D:
                            D[key] = []
                        D[key].append(wfrunid)
            # gridss is labeled gridss_matched in the assay
            # but the actual name may differ based on different naming schemes
            # for research and clinical - check if a run id exists for any gridss workflow
            elif workflow in gridss_workflows:
                for key in gridss_workflows:
                    if key == workflow_info[wfrunid]:
                        if key not in D:
                            D[key] = []
                        D[key].append(wfrunid)
            # bwaMem may be indicated in the assay/pipeline but bwaMem and/or bwamem2
            # might be running in production
            elif workflow in bwa_workflows:
                for key in bwa_workflows:
                    if key == workflow_info[wfrunid]:
                        if key not in D:
                            D[key] = []
                        D[key].append(wfrunid)
            else:
                if workflow_info[wfrunid] == workflow:
                    D[workflow].append(wfrunid)
        
    # remove empty gridds_matched, bcl2fastq and bwamem workflows if an alternative workflow was found
    for i in sequencing_workflows:
        if i != 'bcl2fastq' and i in D and len(D[i]) != 0 and len(D['bcl2fastq']) == 0:
            if 'bcl2fastq' in D:
                del D['bcl2fastq']
    if 'gridss' in D and len(D['gridss']) != 0 and 'gridss_matched' in D and len(D['gridss_matched']) == 0:
        del D['gridss_matched']
    if 'bwamem2' in D and len(D['bwamem2']) != 0 and 'bwaMem' in D and len(D['bwaMem']) == 0:
        del D['bwaMem']
    
            
    return D



def complete_expected_workflows(workflow_info, case_workflows):
    '''
    (dict, list) -> bool
    
    Returns True if all the expected workflows in case workflows have run in 
    in production and have assigned workflow run ids
    
    
    Parameters
    ----------
    - workflow_info (dict): Dictionary mapping worfklow run ids to workflow names
    - case_workflows (list): List of expected workflows from the assay
    '''

    # map the workflow run ids of production workflows to the expected workflows
    expected_workflows = map_expected_workflows_to_runid(workflow_info, case_workflows)

    complete = True
    for workflow in expected_workflows:
        if len(expected_workflows[workflow]) == 0:
            complete = False
    
    return complete


def identify_missing_workflows(workflow_info, case_workflows):
    '''
    (dict, list) -> bool
    
    Returns True if all the expected workflows in case workflows have run in 
    in production and have assigned workflow run ids
        
    Parameters
    ----------
    - workflow_info (dict): Dictionary mapping worfklow run ids to workflow names
    - case_workflows (list): List of expected workflows from the assay
    '''

    # map the workflow run ids of production workflows to the expected workflows
    expected_workflows = map_expected_workflows_to_runid(workflow_info, case_workflows)

    missing = []
    for workflow in expected_workflows:
        if len(expected_workflows[workflow]) == 0:
            missing.append(workflow)
    
    missing = list(set(missing))    
    
    return missing


def list_library_qualif_lims(case_data):
    '''
    (dict) -> list
    
    Returns a list of lim Ids used for library qualification
        
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''

    L = []
       
    seq = json.loads(case_data['case_info']['sequencing'])
    
    # evaluate only lims ids for sequencing
    for i in seq:
        if i['type'] == 'LIBRARY_QUALIFICATION':
            for j in i['limsIds']:
                L.append(j['id'])
            
    return L


def map_samples_to_lims(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with all lims id for each sample in a case    
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''


    # list lims used in library qualification
    qualif = list_library_qualif_lims(case_data)
    
    D = {}
        
    for i in case_data['sample_info']:
        sampleid = i['sampleId']
        limsid = i['limsId']
        # do not include lims for library qualification
        if limsid not in qualif:
            if sampleid not in D:
                D[sampleid] = [limsid]
            else:
                D[sampleid].append(limsid)
                D[sampleid].sort()
            
    return D


def sort_lims_by_samples(tests_samples, samples_lims):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with lists of lims for each sample, if multiple samples
    exist, for each test in a case
        
    Parameters
    ----------
    - test_samples (dict): Dictionart mapping tests with their samples
    - samples_lims (dict): Dictionary mapping samples with their lims ids
    '''
    
    D = {}
    
    for test in tests_samples:
        for sample in tests_samples[test]:
            limsids = samples_lims[sample]
            if test not in D:
                D[test] = [limsids]
            else:
                D[test].append(limsids)
     
    return D



def map_assay_test_to_test_case(assay_test, case_tests):
    '''
    (str, dict) -> str
    
    Returns the test from the case json corresponding to the test from the assay pipeline
    
    Parameters
    ----------
    - assay_test (str): Name of the test in the pipeline json
    - case_tests (dict): Dictionary with case tests from the case json
    '''
    
    
    L = []
    
    test_type, library_type = assay_test.split(':')
    
    for test in case_tests:
        if test_type == 'X':
            assert library_type == test
            L.append(test)
        else:
            if test.upper().startswith('T') or test.upper().startswith('N'):
                if test.startswith(test_type) and library_type in test:
                    L.append(test)
            else:
                if library_type == test:
                    L.append(test)
    
    L = list(set(L))
    
    assert len(L) == 1

    return L[0]




def get_assay_expected_workflows(pipeline_workflows, tests_samples, samples_lims):
    '''
    (dict, dict, dict) -> dict
    
    Returns a dictionary with expected workflows and corresponding lims according to the defined assays, pipeline
    and case information
    
    Parameters
    ----------
    - pipeline_workflows (dict): Dictionary with pipeline workflows
    - tests_samples (dict): Dictionary with all the samples mapping each test
    - samples_limns (dict): Dictionary with the lims mapping each assay 
    '''
       
    workflows = {}

    # loop over workflows in pipeline
    for workflow in pipeline_workflows:
        # get the expected tests for each workflow from the assay definition
        level = pipeline_workflows[workflow]['level']
        # collect the expected lims for each workflow depending on the case information
        if level == 'lane':
            # each lims id of each test should have a separate workflow run id
            tests = pipeline_workflows[workflow]['tests'].split(',')
            # collect the limsids for each test using case data
            for test in tests:
                # find the corresponding test in case data
                case_test = map_assay_test_to_test_case(test, tests_samples)
                # get the sample id - each test can have multiple samples
                for sampleid in tests_samples[case_test]:
                    # get the corresponding lims 
                    limsids = samples_lims[sampleid]
                    for lims in limsids:
                        if workflow in workflows:
                            workflows[workflow].append({'workflow': workflow,'test': [case_test], 'sampleid': sampleid, 'limsids': lims, 'parents': [], 'parent_workflows': []})
                        else:
                            workflows[workflow] = [{'workflow': workflow, 'test': [case_test], 'sampleid': sampleid, 'limsids': lims, 'parents': [], 'parent_workflows': []}]
        
        elif level == 'merge':
            
            if '|' in pipeline_workflows[workflow]['tests']:
                tests = pipeline_workflows[workflow]['tests'].split('|')
                # collect the limsids for each test using case data
                # get the expected combinations of lims for each combination of test samples
                L = []
                S = []
                
                case_tests = []
                
                for test in tests:
                    # find the corresponding test in case data
                    case_test = map_assay_test_to_test_case(test, tests_samples)
                    case_tests.append(case_test)
                    
                    l = []
                    s = []
                    # get the sample id - each test can have multiple samples
                    for sampleid in tests_samples[case_test]:
                        # get the corresponding lims 
                        limsids = samples_lims[sampleid]
                        l.append(limsids)
                        s.append(sampleid)
                    L.append(l)
                    S.append(s)
                
                combined_lims = list(itertools.product(*L))
                combined_samples = list(itertools.product(*S))
                
                
                # merge and sort each set of lims for each set of combined tests 
                for i in range(len(combined_lims)):
                    merged_lims = []
                    merged_samples = []
                    for j in combined_lims[i]:
                        merged_lims.extend(j)
                    for k in combined_samples[i]:
                        merged_samples.append(k)
                    
                    merged_lims = ','.join(sorted(merged_lims))
                    merged_samples = ','.join(sorted(merged_samples))
                               
                    if workflow in workflows:
                        workflows[workflow].append({'workflow': workflow, 'test': case_tests, 'sampleid': merged_samples, 'limsids': merged_lims, 'parents': [], 'parent_workflows': []})
                    else:
                        workflows[workflow] = [{'workflow': workflow, 'test': case_tests, 'sampleid': merged_samples, 'limsids': merged_lims, 'parents': [], 'parent_workflows': []}]
           
            elif ',' in pipeline_workflows[workflow]['tests'] or \
            (',' not in pipeline_workflows[workflow]['tests'] and '|' not in pipeline_workflows[workflow]['tests']):
                tests = pipeline_workflows[workflow]['tests'].split(',')  
                # collect the limsids for each test using case data
                for test in tests:
                    # find the corresponding test in case data
                    case_test = map_assay_test_to_test_case(test, tests_samples)
                    # get the sample id - each test can have multiple samples
                    for sampleid in tests_samples[case_test]:
                        # get the corresponding lims 
                        limsids = samples_lims[sampleid]
                        # the workflow has all the lims
                        limsids = ','.join(sorted(list(limsids)))
                        if workflow in workflows:
                            workflows[workflow].append({'workflow': workflow, 'test': [case_test], 'sampleid': sampleid, 'limsids': limsids, 'parents': [], 'parent_workflows': []})
                        else:
                            workflows[workflow] = [{'workflow': workflow, 'test': [case_test], 'sampleid': sampleid, 'limsids': limsids, 'parents': [], 'parent_workflows': []}]
                
            
    return workflows        



def get_production_workflows(samples_workflows, workflow_lims):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary mapping each workflow run id, their lims and sample to each workflow
    
    Parameters
    ----------
    - samples_workflows (dict): Dictionary mapping all lims for each sample
    - workflow_lims (dict): Dictionary mapping the lims to each workflow
    '''        
        
    # reorganize data: {workflow: {{'wfrunid': ,'limsids':, 'samples':}}}    
        
    D = {}
    
    for wfrunid in workflow_lims:
        lims = ','.join(sorted(workflow_lims[wfrunid]))
        samples = []
        names = []
        for sample in samples_workflows:
            for d in samples_workflows[sample]['workflows']:
                if d['wfrun_id'] == wfrunid:
                    samples.append(sample)
                    names.append(d['workflow'])
                    break
        names = list(set(names)) 
        assert len(names) == 1
        name = names[0]
        samples = ','.join(sorted(samples))
        if name not in D:
            D[name] = {}
        D[name][wfrunid] = {'limsids': lims, 'samples': samples}
            
    return D        



def is_incomplete_workflow_run(d):
    '''
    (dict) -> bool
    
    Returns True is any key in d is missing values (expect parents)
        
    Parameters
    ----------
    - d (dict): Dictionary with workflow run id information in case_analysis
    '''
    
    # exclude parents 
    vals = [d[i] for i in d.keys() if i != 'parents']
    return any(map(lambda x: x is None or len(x) == 0, vals))




def identify_workflows_with_missing_data(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> list

    Returns a list of workflows with missing data
            
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''

    missing = [workflow for workflow in cases_analysis if workflow not in expected_workflow_lims]
        
    # check if there are missing iterations
    for workflow in cases_analysis:
        if workflow in expected_workflow_lims:
            if len(cases_analysis[workflow]) < len(expected_workflow_lims[workflow]):
                missing.append(workflow)
                
    # check that all workflows have been identified
    for workflow in cases_analysis:
        for d in cases_analysis[workflow]:
            # analysis is incomplete if any workflow in assay has missing information
            if is_incomplete_workflow_run(d):
                missing.append(workflow)
                       
    missing = list(set(missing))

    return missing                    


def find_production_workflow(production_workflows, d):
    '''
    (dict, dict) -> dict    
    
    Returns a dictionary with prodcution data mapping the expected data for a specific workflow
                
    Parameters
    ----------
    - production_workflows (dict): Dictionary with case data extracted from the provenance reporter
    - d (dict): Dictionary with expected workflow information based on assay and case info
    '''
    
    # for sequencing workflows, the assay may indicate bcl2fastq but the 
    # sequencing workflows may be diferent if data is injected
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
    
    # gridss_matched is always indicated in the assays but the actual workflow
    # could be gridss or gridss_matched (same workflow but different names in research and clinical)
    gridss_workflows = ['gridss_matched', 'gridss']

    # bwaMem may be indicated in the assay but the actual workflow might be bwMem or bwamem2
    bwa_workflows = ['bwaMem', 'bwamem2']


    data = {'workflow': None, 'limsids': None, 'wfrunid': None, 'tests': None, 'samples': None, 'parents': []}
    workflow = d['workflow']
    expected_lims = d['limsids']
    expected_samples = d['sampleid']
    test = d['test']
    #  find the workflow in production with the expected limsids and samples
    if workflow in sequencing_workflows:
        # find the actual sequencing workflow as it may differ from assay
        for key in sequencing_workflows:
            if key in production_workflows:
                for wfrunid in production_workflows[key]:
                    limsids = production_workflows[key][wfrunid]['limsids']
                    samples = production_workflows[key][wfrunid]['samples']
                    if expected_samples == samples and expected_lims == limsids:
                        ### check that only 1 wfrunids match the requirement
                        assert data['wfrunid'] is None 
                        # update data collector
                        data['limsids'] = limsids
                        data['samples'] = samples
                        data['wfrunid'] = wfrunid
                        data['tests'] = test
                        data['workflow'] = key
    elif workflow in gridss_workflows:
        # find the gridss workflow as it may differ from assay
        for key in gridss_workflows:
            if key in production_workflows:
                for wfrunid in production_workflows[key]:
                    limsids = production_workflows[key][wfrunid]['limsids']
                    samples = production_workflows[key][wfrunid]['samples']
                    if expected_samples == samples and expected_lims == limsids:
                        ### check that only 1 wfrunids match the requirement
                        assert data['wfrunid'] is None 
                        # update data collector
                        data['limsids'] = limsids
                        data['samples'] = samples
                        data['wfrunid'] = wfrunid
                        data['tests'] = test
                        data['workflow'] = key
    elif workflow in bwa_workflows:
        # bwaMem and bwmem2 may both have been running in production
        # use workflow defined in pipeline if it exists
        # match alternative if expected workflow does not exist in production
        if workflow in production_workflows:
            for wfrunid in production_workflows[workflow]:
                limsids = production_workflows[workflow][wfrunid]['limsids']
                samples = production_workflows[workflow][wfrunid]['samples']
                if expected_samples == samples and expected_lims == limsids:
                    ### check that only 1 wfrunids match the requirement
                    assert data['wfrunid'] is None 
                    # update data collector
                    data['limsids'] = limsids
                    data['samples'] = samples
                    data['wfrunid'] = wfrunid
                    data['tests'] = test
                    data['workflow'] = workflow
        else:
            # find the bwa workflow as it may differ from assay
            for key in bwa_workflows:
                if key in production_workflows:
                    for wfrunid in production_workflows[key]:
                        limsids = production_workflows[key][wfrunid]['limsids']
                        samples = production_workflows[key][wfrunid]['samples']
                        if expected_samples == samples and expected_lims == limsids:
                            ### check that only 1 wfrunids match the requirement
                            assert data['wfrunid'] is None 
                            # update data collector
                            data['limsids'] = limsids
                            data['samples'] = samples
                            data['wfrunid'] = wfrunid
                            data['tests'] = test
                            data['workflow'] = key
    else:
        if workflow in production_workflows:
            for wfrunid in production_workflows[workflow]:
                limsids = production_workflows[workflow][wfrunid]['limsids']
                samples = production_workflows[workflow][wfrunid]['samples']
                if expected_samples == samples and expected_lims == limsids:
                    ### check that only 1 wfrunids match the requirement
                    assert data['wfrunid'] is None 
                    # update data collector
                    data['limsids'] = limsids
                    data['samples'] = samples
                    data['wfrunid'] = wfrunid
                    data['tests'] = test
                    data['workflow'] = workflow
     
    return data         


def map_expected_production_workflows(expected_workflow_lims, production_workflows):
    '''
    (dict, dict) -> dict    
    
    Returns a dictionary with prodcution data mapping the expected data from the assay ans case info
    with the production data available for a case in the provenance reporter
            
    Parameters
    ----------
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    - production_workflows (dict): Dictionary with case data extracted from the provenance reporter
    '''

    D = {}
        
    for workflow in expected_workflow_lims:
        D[workflow] = []
        for d in expected_workflow_lims[workflow]:
            data = find_production_workflow(production_workflows, d)
            D[workflow].append(data)
            
    return D            



def is_data_complete(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> bool    
    
    Returns True if each workflow in case_analysis has complete information
        
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''
        
    complete = True
        
    if cases_analysis.keys() != expected_workflow_lims.keys():
        complete = False
    
    for workflow in cases_analysis:
        if len(cases_analysis[workflow]) < len(expected_workflow_lims[workflow]):
            complete = False
    
    # check that all workflows have been identified
    for workflow in cases_analysis:
        for d in cases_analysis[workflow]:
            # analysis is incomplete if any workflow in assay has missing information
            if is_incomplete_workflow_run(d):
                complete = False
                
    return complete


def no_extra_data(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> bool    
    
    Returns True is each workflow in case_analysis have a single workflow run id
    matching the lims requirements
        
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''
    
    no_extra = True
        
    # check if there are extra workflows
    for workflow in cases_analysis:
        if len(cases_analysis[workflow]) > len(expected_workflow_lims[workflow]):
            no_extra = False
    
    return no_extra


def identify_extra_workflows(cases_analysis, expected_workflow_lims):
    '''
    (dict, dict) -> list    
    
    Returns a list of workflow with multiple run ids matching the lims requirements
    
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - expected_workflow_lims (dict): Dictionary with expected workflow and lims from assay and case info
    '''
    
    extra = []
       
    # check if there are extra workflows
    for workflow in cases_analysis:
        if len(cases_analysis[workflow]) > len(expected_workflow_lims[workflow]):
            extra.append(workflow)
    
    extra = list(set(extra))
    
    return extra


def reformat_pipeline_workflows(L):
    '''
    (list) -> dict
        
    Returns a dictionary with all expected workflows in a pipeline
    
    Parameters
    ----------
    - L (list): List of expected workflows for a given pipeline from pipelines.json
    '''
    
    D = {}
        
    for d in L:
        workflow = d['workflows']
        assert workflow not in D
        D[workflow] = d
    
    return D


def add_parent_workflows(case_analysis, parent_to_children_workflows):
    '''
    (dict, dict) -> dict
    
    Add the parent workflow run ids to each analysis workflow in case_analysis
    
    Parameters
    ----------
    - cases_analysis (dict): Dictionary with case production data
    - parent_to_children_workflows (dict): Dictionary with parent-children workflow relationships
    '''
    
    for workflow in case_analysis:
        for d in case_analysis[workflow]:
            wfrunid = d['wfrunid']
            for parent in parent_to_children_workflows:
                if wfrunid in parent_to_children_workflows[parent]:
                    d['parents'].append(parent)
            
    return case_analysis    







def review_data(provenance_data_file, assay_file, pipeline_file, database, table='templates'):
    '''
    (str, str, str, str, str) -> None 
    
    Generates sqlite database with templates and review for all projects and cases in the
    provenance data file
    
    Parameters
    ----------
    - provenance_data_file (str): Path to the file with production data extracted from Shesmu
    - assay_file (str): Path to the json file mapping assays and pipelines
    - pipeline_file (str): Path to the json file mapping workflows and pipelines
    - database (str): Path to the sqlite database
    - table (str): Table in database storing the analysis data
    '''
    
    #assays.json : list of assays/version and assigned pipelines/versions
    #pipelines.json : list of pipelines/versions with expected workflows and associated information
    
    
    infile = open(assay_file)
    assays = json.load(infile)
    infile.close()
    
    infile = open(pipeline_file)
    pipelines = json.load(infile)
    infile.close()
    
    # load production data
    provenance_data = load_data(provenance_data_file)
    print('loaded data')
    
    # create database if file doesn't exist
    if os.path.isfile(database) == False:
        initiate_db(database, 'analysis_review', ['templates'])
    print('initiated database')    
    
    # collect the recorded md5sums of the donor data from the database
    recorded_md5sums = get_cases_md5sum(database, table = 'templates')
    print('pulled md5sums from database')
    
    # track all cases in production
    P = []
          
    
    
    
    # make a list of problematic cases to explore later
    
    # exclude_cases = ['R5523_a141_GTNBP_0001_Bn_P',
    #                  'R5526_a120_BDWGTS_0198_Ut_M',
    #                  'R5526_a120_BIODIVA_0025_Om_M',
    #                  'R5526_a120_BIODIVA_0149_Om_M',
    #                  'R5526_a120_BIODIVA_0174_Ae_M',
    #                  'R5526_a120_BIODIVA_0200_So_M',
    #                  'R5526_a120_BIODIVA_0226_Ov_P',
    #                  'R5526_a120_BIODIVA_0286_nn_M',
    #                  'R5526_a120_BIODIVA_0392_Ov_P']
    
    exclude_cases = ['R5523_a141_GTNBP_0001_Bn_P']
    
    
    
    
    # 'R5523_a141_GTNBP_0001_Bn_P': multiple workflow runs with same lims
    # 'R5526_a120_BDWGTS_0198_Ut_M': tests labeled WG in case data, normal ? tumor?
    # 'R5526_a120_BIODIVA_0025_Om_M': tests labeled WG in case data, normal ? tumor?
    # 'R5526_a120_BIODIVA_0149_Om_M':  tests labeled WG in case data, normal ? tumor?
    
    
    
    for case_data in provenance_data:
        # record data to insert
        L = []
        case_id = case_data['case']
        
        if case_id in exclude_cases:
            continue
        
        print(case_id)
        
        
        P.append(case_id)
        # compute the md5sum of the case info
        md5sum = compute_md5(case_data)
        # check that case needs to be updated
        if case_to_update(recorded_md5sums, case_id, md5sum):
            # md5sums differ or case id not in analysis review cache
            # remove case id from database
            conn = connect_to_db(database)
            delete_unique_record(case_id, conn, database, 'templates', 'case_id')
            conn.close()
            
            # check that assay and version are defined in the assay config
            # get the assay name and version for the case
            assay_name = case_data['assay'].split('_')
            version = assay_name[-1]
            assay_version = 'v' + version
            assay_name = '_'.join(assay_name[:-1])
            
            # collect and evaluate data for each pipeline
            case_analysis = {}
            pipeline_analysis = {}            
                       
            # case data may be incomplete - check project and deliverables can be retrieved  
            try:
                project_ids = [case_data['project_info'][i]['project'] for i in range(len(case_data['project_info']))]
            except:
                project_ids = []
            
            try:
                donor = get_donor_name(case_data)
            except:
                donor = ''
        
        
            ### exclude biodiva and hbseq projects for now
            # issue witrh test --> WG
            
            
            if 'BIODIVA' in project_ids or 'HBSEQ' in project_ids:
                continue
            
            
        
        
            # check that case data is complete (all sections in the case dictionary are complete)
            if is_case_info_complete(case_data):
                
                
                # review analysis only if signoff is complete
                if is_signoff_complete(case_data):
                    
                                   
                    if assay_name in assays:
                        
                        if assay_version in assays[assay_name]:
                            
                            
                            # extract case data
                            # map tests to lims ids
                            tests_limsids = map_lims_to_tests(case_data)
                            # map tests to samples
                            tests_samples = map_tests_to_samples(case_data, tests_limsids)
                            # map samples to lims
                            samples_lims = map_samples_to_lims(case_data)
                            # sort lims ids by sample for each test (ie. if multiple samples per test exist)                            
                            tests_limsids = sort_lims_by_samples(tests_samples, samples_lims)
                            # collect all worfklows for the assay
                            # get all the workflow information
                            workflow_info = extract_workflow_information(case_data)
                            # map all samples to each workflow
                            # get workflows of all samples for the case
                            samples_workflows = collect_sample_workflows(case_data)
                            # map each workflow to its limsIds
                            workflow_lims = map_workflows_to_lims(case_data)
                            # find the parent-children worklow relationships
                            parent_to_children_workflows = collect_workflow_relationships(case_data)
                    
                            for pipeline_name in assays[assay_name][assay_version]:
                                
                                
                                
                                
                                
                                pipeline_version = assays[assay_name][assay_version][pipeline_name]
                                # get the pipeline expected worflows 
                                pipeline_workflows = reformat_pipeline_workflows(pipelines[pipeline_name][pipeline_version])
                    
                                
                                # did all expected workflows in config ran?
                                if complete_expected_workflows(workflow_info, pipeline_workflows):
                                    
                                    
                                    
                                    # get expected lims for each workflow based on the assay and the case
                                    expected_workflow_lims = get_assay_expected_workflows(pipeline_workflows, tests_samples, samples_lims)
                                    # get lims, samples and run ids for each workflow seen in production
                                    production_workflows = get_production_workflows(samples_workflows, workflow_lims)
                                    # did all the expected workflows ran for all tests (check lims)?
                                    pipeline_analysis = map_expected_production_workflows(expected_workflow_lims, production_workflows)
                                    # add parent workflows
                                    pipeline_analysis = add_parent_workflows(pipeline_analysis, parent_to_children_workflows)
                                                                    
                                    # check if missing data (workflows and parents)
                                    if is_data_complete(pipeline_analysis, expected_workflow_lims):
                                        
                                        
                                        # check if some workflows have extra iterations matching the required lims
                                        if no_extra_data(pipeline_analysis, expected_workflow_lims):
                                            # data passed all the checks
                                            valid = 1
                                            error = ''
                                            
                                        else:
                                            extra_workflows = identify_extra_workflows(pipeline_analysis, expected_workflow_lims)
                                            error = '[EXTRA WORKFLOWS]: Workflows have unexpected multiple runs {0}'.format(','.join(extra_workflows)) 
                                            valid = 0
                                        
                                    else:
                                        missing = identify_workflows_with_missing_data(pipeline_analysis, expected_workflow_lims)
                                        error = '[INCOMPLETE DATA]: Workflows are missing {0}'.format(','.join(sorted(list(set(missing)))))
                                        valid = 0
                               
                                else:
                                    # get the missing workflows
                                    missing_workflows = identify_missing_workflows(workflow_info, pipeline_workflows)
                                    error = '[MISSING WORKFLOWS]: missing {0}'.format(','.join(sorted(missing_workflows)))
                                    valid = 0
                                    pipeline_analysis = {}
                                
                                
                                case_analysis[pipeline_name] = {'pipeline_analysis': pipeline_analysis,
                                                                        'error': error,
                                                                        'valid': valid}
                        else:
                            error = '[ASSAY VERSION]: version not matching {0} in assay config'.format(assay_name)
                            valid = 0
                    else:
                        error = '[ASSAY]: assay {0} not in assay config'.format(assay_name)
                        valid = 0
                                       
                else:
                    error = '[INCOMPLETE SEQUENCING]: some tests have incomplete sequencing'
                    valid = 0
                    
            else:
                error = '[INCOMPLETE CASE]: case is missing some data'
                valid = 0
                
            # the case may be in multiple projects. record data for each project the case belongs to 
            if project_ids:
                for project in project_ids:
                    if case_analysis:
                        error = ';'.join(sorted([case_analysis[i]['error'] for i in case_analysis]))
                        valid = int(all([case_analysis[i]['valid'] for i in case_analysis]))
                        L.append([case_id, donor, project, assay_name, json.dumps(case_analysis), str(valid), error, md5sum])
                    else:
                        valid = 0
                        L.append([case_id, donor, project, assay_name, json.dumps({}), str(valid), error, md5sum])
            else:
                assert len(case_analysis) == 0
                valid = 0
                L.append([case_id, donor, 'NA', assay_name, json.dumps(case_analysis), str(valid), error, md5sum])
                            
            if L:
                conn = connect_to_db(database)
                insert_multiple_records(L, conn, database, 'templates', define_columns('analysis_review')['templates']['names'])
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

       
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = 'analysis_review.py', description='Script to generate the analyais review cache')
    parser.add_argument('-pv', '--provenance', dest = 'provenance', default = '/scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json',
                        help = 'Path to the provenance reporter data json. Default is /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json')
    parser.add_argument('-ad', '--analysis', dest='analysis_db', default = '/scratch2/groups/gsi/production/waterzooi/analysis_review_case.db', 
                        help='Path to the analysis review database')    
    parser.add_argument('-as', '--assays', dest='assay_file', help='Path to the json with assay definitions', required = True)    
    parser.add_argument('-pi', '--pipelines', dest='pipeline_file', help='Path to the json with pipeline definitions', required = True)    
    
    # get arguments from the command line
    args = parser.parse_args()
    # generate sqlite cache
    review_data(args.provenance, args.assay_file, args.pipeline_file, args.analysis_db)

       