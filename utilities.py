# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:04:02 2023

@author: rjovelin
"""

import sqlite3
import time
import string
import random
import requests
from db_helper import connect_to_db





def remove_non_analysis_workflows(L):
    '''
    (list) -> list
    
    Returns a list L of dictionaries with workflows, removing any non-analysis workflows
    
    Parameters
    ----------
    - L (list): List of dictionaries with workflow name and workflow ids
    '''
    
    non_analysis_workflows = ('wgsmetrics', 'insertsizemetrics', 'bamqc', 'calculatecontamination',
                              'callability', 'fastqc', 'crosscheckfingerprintscollector',
                              'fingerprintcollector', 'bwamem', 'bammergepreprocessing',
                              'ichorcna', 'tmbanalysis', 'casava', 'bcl2fastq',
                              'fileimportforanalysis', 'fileimport', 'import_fastq',
                              'dnaseqqc', 'hotspotfingerprintcollector', 'rnaseqqc')
    
    to_remove = [i for i in L if i['wf'].split('_')[0].lower() in non_analysis_workflows or i['wf'].lower() in non_analysis_workflows]
    for i in to_remove:
        L.remove(i)
    
    return L



def convert_epoch_time(epoch):
    '''
    (str) -> str
    
    Returns epoch time in readable format
    
    Parameters
    ----------
    - epoch (str)
    '''
    
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(epoch)))



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


def get_donors(project_name, database):
    '''
    (str, str) -> list
    
    Returns a list of donors for a given project
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
    
    # connect to db
    conn = connect_to_db(database)
    # extract library source
    data = conn.execute("SELECT DISTINCT case_id FROM Samples WHERE project_id = ?;", (project_name,)).fetchall()
    donors = list(set([i['case_id'] for i in  data]))
    conn.close()
    
    return donors



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
    
    
def get_case_md5sums(database, project_name):
    '''
    (str, str) -> dict
    
    Returns a dictionary with case, md5sums from the waterzooi database
    
    Parameters
    ----------
    - database (str): Path to the waterzooi sqlite database
    - project_name (str): Name of the project of interest
    '''
    
    # connect to db
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT case_id, md5 FROM Checksums WHERE project_id = ?;", (project_name,)).fetchall()
    conn.close()
    
    D = {}
    for i in data:
        case = i['case_id']
        md5sum = i['md5']
        D[case] = md5sum
    
    return D
    

def extract_case_signoff(case, nabu_key_file, nabu='https://nabu.gsi.oicr.on.ca/case'):
    '''
    (str, str, str) -> dict
    
    Returns a list of signoffs for that case
        
    Parameters
    ----------
    - case (str): Case identifier
    - nabu_key_file (str): File storing the nabu API key
    - nabu (str): URL to access the case in Nabu
    '''
    
    infile = open(nabu_key_file)
    nabu_key = infile.read().rstrip()
    infile.close()
    
    headers = {'accept': 'application/json',
               'X-API-KEY': nabu_key,}
    
    D = {}
    
    response = requests.get(nabu + '/{0}/sign-off'.format(case), headers=headers)
    if response.status_code == 200:
        for d in response.json():
            case = d['caseIdentifier']
            if case not in D:
                D[case] = {}
            step = d['signoffStepName']
            step = ' '.join(list(map(lambda x: x.lower().capitalize(), step.split('_'))))
            if step in D[case]:
                D[case][step].append(d)
            else:
                D[case][step] = [d]
    return D

    

def extract_nabu_signoff(cases, nabu_key_file, nabu='https://nabu.gsi.oicr.on.ca/case/sign-off'):
    '''
    (list, str, str) -> dict
    
    Returns a dictionary of signoffs for each case in cases
        
    Parameters
    ----------
    - cases (list): List of case identifiers
    - nabu_key_file (str): File storing the nabu API key
    - nabu (str): URL to access the signoffs in Nabu
    '''
    
    infile = open(nabu_key_file)
    nabu_key = infile.read().rstrip()
    infile.close()
    
    headers = {'accept': 'application/json',
               'X-API-KEY': nabu_key,}
    
    D = {}
    
    response = requests.get(nabu, headers=headers)
    if response.status_code == 200:
        for d in response.json():
            case_id = d['caseIdentifier']
            if case_id in cases:
                comment = d['comment']
                if comment and comment.startswith('G') and '-' in comment:
                    comment = comment.split('-')
                    c = ['-'.join([comment[0], comment[i]]) for i in range(1, len(comment))]
                else:
                    if comment:
                        c = [d['comment']]
                    else:
                        c = d['comment']
                d['comment'] = c
                if case_id not in D:
                    D[case_id] = {}
                step = d['signoffStepName']
                step = ' '.join(list(map(lambda x: x.lower().capitalize(), step.split('_'))))
                if step in D[case_id]:
                    D[case_id][step].append(d)
                else:
                    D[case_id][step] = [d]
    return D


def list_signoff_deliverables(signoffs):
    '''
    (dict) -> dict
    
    Returns a dictionary with the deliverables available for release signoff for each case 
    
    Parameters
    ----------
    - signoffs (dict): Dictionary with case signoffs
    '''
    
    D = {}
    
    for case in signoffs:
        D[case] = []
        if 'Release' in signoffs[case]:
            for d in signoffs[case]['Release']:
                D[case].append(d['deliverable'])
            D[case] = list(set(D[case]))
    
    return D
    

def remove_cases_with_no_approval_signoff(analysis_data, signoffs):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with analysis workflows for cases with release approval signoff
       
    Parameters
    ----------
    - analysis_data (dict): Dictionary with selected analysis workflows for each case 
    - signoffs (dict): Case signoffs extracted from Nabu
    '''
    
    to_remove = []
    for case_id in analysis_data:
        if case_id not in signoffs:
            to_remove.append(case_id)
        else:
            if 'Release Approval' not in signoffs[case_id]:
                to_remove.append(case_id)
            elif all([d['qcPassed'] for d in signoffs[case_id]['Release Approval']]) == False:
                to_remove.append(case_id)
    
    to_remove = list(set(to_remove))
    for case_id in to_remove:
        del analysis_data[case_id]

    return analysis_data



    
def remove_cases_with_competed_cbioportal_release(analysis_data, signoffs, deliverable):
    '''
    (dict, dict, str) -> list
    
    Returns a dictionary with analysis workflows for cases for which cbioportal signoff is not complete
        
    Parameters
    ----------
    - analysis_data (dict): Dictionary with selected analysis workflows for each case 
    - signoffs (dict): Case signoffs extracted from Nabu
    - deliverable (str): Selected option in waterzooi
    '''
    
    to_remove = []
    # remove cases for which release signoffs are complete    
    for case_id in analysis_data:
        if deliverable in ['sequenza', 'purple']:
            if 'Release' in signoffs[case_id]:
                for d in signoffs[case_id]['Release']:
                    if 'cbioportal' in d['deliverable'].lower() and d['qcPassed']:
                        # release complete 
                        to_remove.append(case_id)

    to_remove = list(set(to_remove))
    for case_id in to_remove:
         del analysis_data[case_id]        
    
    return analysis_data     



def remove_workflows_with_deliverable_signoff(analysis_data, signoffs, deliverable, deliverable_type):
    '''
    (dict, dict, str, str) -> dict
    
    Returns a dictionary with analysis workflows part of a deliverable that is not signed off for each case
      
    Parameters
    ----------
    - analysis_data (dict): Dictionary with selected analysis workflows for each case 
    - signoffs (dict): Case signoffs extracted from Nabu
    - deliverable (str): Selected option in waterzooi 
    - deliverable_type: type of deliverable :pipeline or fastq
    '''
    
    fastq_workflows = ['bcl2fastq', 'casava', 'fileimport', 'fileimportforanalysis', 'import_fastq']
    
    for case_id in analysis_data:
        # remove workflows
        remove_workflows = False
        if deliverable in ['all', 'selected', 'standard', 'MOH_pipeline']:
            if 'Release' in signoffs[case_id]:
                for d in signoffs[case_id]['Release']:
                    if deliverable_type in d['deliverable'].lower() and d['qcPassed']:
                        remove_workflows = True
                        break
        if remove_workflows:
            if deliverable_type == 'pipeline':
                L = [i for i in analysis_data[case_id] if i.lower() not in fastq_workflows]
            elif deliverable_type == 'fastq':
                L = [i for i in analysis_data[case_id] if i.lower() in fastq_workflows]
        else:
            L = []
        
        if L:
            for i in L:
                del analysis_data[case_id][i]
   
    remove_case = [case_id for case_id in analysis_data if len(analysis_data[case_id]) == 0]
    if remove_case:
        for case_id in remove_case:
            del analysis_data[case_id]
   
    return analysis_data                    
                        


def cbioportal_format(analysis_data):
    '''
    (dict) -> dict

    Returns a dictionary with cbioportal data keeping only donors and samples
    and removing case identifiers     
    
    Parameters
    ----------
    - analysis_data (dict): Dictionary with cbioportal data for one or more cases
    '''
    
    D = {}
    
    for case_id in analysis_data:
        for donor in analysis_data[case_id]:
            for sample in analysis_data[case_id][donor]:
                if donor not in D:
                    D[donor] = {}
                assert sample not in D[donor]
                D[donor][sample] = analysis_data[case_id][donor][sample]
    
    return D
    


def moh_format(analysis_data, donors):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with workflow information for a given block (ie, sample pair)
    and anchor bmpp parent workflow
    
    Parameters
    ----------
    - analysis_data (dict): Dictionary with analysis workflow data for one or more cases
    - donors (dict): Dictionary with donor id mapped to case id
    '''
        
    groups = {'purple': 'calls.copynumber',
              'sequenza': 'calls.copynumber',
              'varscan': 'calls.copynumber',
              'gridss': 'calls.copynumber',
              'rsem': 'calls.expression',
              'bammergepreprocessing': 'alignments_WG.callready',  
              'haplotypecaller': 'calls.germline.mutations',
              'starfusion': 'calls.fusions',
              'arriba': 'calls.fusions',
              'mavis': 'calls.structuralvariants',
              'delly': 'calls.structuralvariants',
              'star_call_ready': 'alignments_WT.callready',
              'varianteffectpredictor': 'calls.mutations',           
              'msisensor': 'calls.msi',
              'hrdetect': 'calls.hrd'}
              
                        
    D = {}
    
    for case_id in analysis_data:
        # get the donor
        donor = donors[case_id]
        for workflow in analysis_data[case_id]:
            name = workflow.split('_')[0].lower()
            assert name in groups
            group = groups[name]
            for workflow_id in analysis_data[case_id][workflow]:
                if donor not in D:
                    D[donor] = {}
                if group not in D[donor]:
                    D[donor][group] = []
                D[donor][group].extend(analysis_data[case_id][workflow][workflow_id])    
                    
    return D
    

def get_workflow_file_qc(database, case_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary with file qc status for each output file
    of every workflows of a given case
        
    Parameters
    ----------
    - database (str): Path to the waterzooi sqlite database
    - case_id (str): Case identifier
    '''
    
    # connect to db
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT wfrun_id, username, ticket, qcstatus FROM File_qc WHERE case_id = ?;", (case_id,)).fetchall()
    conn.close()
    
    D = {}

    for i in data:
        if i['wfrun_id'] not in D:
            D[i['wfrun_id']] = [i['qcstatus']]
        else:
            D[i['wfrun_id']].append(i['qcstatus'])
    
    return D    
    
    
def get_workflow_release_status(database, case_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary with the release status of each workflow id of a given case
    The release status is derived from the file qc status in Nabu of the workflow output files
    
    Parameters
    ----------
    - database (str): Path to the waterzooi sqlite database
    - case_id (str): Case identifier
    '''

    # get the file qc status for each output file of every workflows
    workflow_qc = get_workflow_file_qc(database, case_id)
    
    D = {}    

    for workflow_id in workflow_qc:
        if all(map(lambda x: x.isdigit(), workflow_qc[workflow_id])):
            if all(map(lambda x: int(x), workflow_qc[workflow_id])):
                D[workflow_id] = True
            elif any(map(lambda x: int(x), workflow_qc[workflow_id])):
                D[workflow_id] = True
            elif all(map(lambda x: int(x), workflow_qc[workflow_id])) == False:
                D[workflow_id] = False
        elif '1' in workflow_qc[workflow_id]:
            D[workflow_id] = True
        elif len(list(set(workflow_qc[workflow_id]))) == 1:
            D[workflow_id] = '?'
        
        
       
        
    return D    
         

    
def get_file_release_status(database, case_id):
    '''
    (str, str) -> dict
    
    Returns a dictionary with the release status of each file of a given case
    The release status is derived from the file qc status in Nabu 
    
    Parameters
    ----------
    - database (str): Path to the waterzooi sqlite database
    - case_id (str): Case identifier
    '''
    
    # connect to db
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT wfrun_id, filepath, username, ticket, qcstatus FROM File_qc WHERE case_id = ?;", (case_id,)).fetchall()
    conn.close()
    
    D = {}

    for i in data:
        assert i['filepath'] not in D
        if i['qcstatus'] == '1':
            D[i['filepath']] = True
        elif i['qcstatus'] == '0':
            D[i['filepath']] = False
        else:
            D[i['filepath']] = '?'

    return D    
    

def template_error_formatting(errors):
    '''
    (list) -> list
    
    Returns a list of error messages for each template of a single case formatted
    to display a list of elements
        
    Parameters
    ----------
    - errors (list): List with error messages for each template of a single case
    '''
    
    for i in range(len(errors)):
        # extract the list of workflows
        if '(' in errors[i]:
            wk = errors[i][errors[i].index('(')+1:errors[i].index(')')]
            # replace workflows  from error message
            message = errors[i].replace(wk, '').replace(' ()', '')
            wk = wk.split(':')
        else:
            message = errors[i]
            wk = ''
        errors[i] = {'message': message, 'workflows': wk}
        
        
        
    return errors
        
        
def case_error_formatting(errors):
    '''
    (str) -> dict
    
    Returns a dictionary with error message for the combined error message across
    templates for a single case
            
    Parameters
    ----------
    - errors (str): Combined error across templates of a single case
    '''
    
    # extract the list of workflows
    if '(' in errors:
        wk = errors[errors.index('(')+1:errors.index(')')]
        # replace workflows  from error message
        message = errors.replace(wk, '').replace(' ()', '')
        wk = wk.split(':')
    else:
        message = errors
        wk = ''
    errors = {'message': message, 'workflows': wk}
        
    return errors
   
        
   
    
        
        
        
    
    
    
    


