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



def connect_to_db(database):
    '''
    (str) -> sqlite3.Connection
    
    Returns a connection to SqLite database prov_report.db.
    This database contains information extracted from FPR
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    '''
    
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    return conn


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
            case = d['caseIdentifier']
            if case in cases:
                if case not in D:
                    D[case] = {}
                step = d['signoffStepName']
                step = ' '.join(list(map(lambda x: x.lower().capitalize(), step.split('_'))))
                if step in D[case]:
                    D[case][step].append(d)
                else:
                    D[case][step] = [d]
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
    
    
    
    
    