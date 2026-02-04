# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 21:11:22 2023

@author: rjovelin
"""


from utilities import connect_to_db, convert_epoch_time


def get_project_info(database, project_name=None):
    '''
    (str, str | None) -> list
    
    Returns a list with project information extracted from database for all projects 
    of for a single project if project_name is defined 
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - project_name (None | str): Project of interest
    '''
    
    # connect to db
    conn = connect_to_db(database)
    if project_name:
        # extract project info
        project = conn.execute('SELECT * FROM Projects WHERE project_id=?', (project_name,)).fetchall()
    else:
        project = conn.execute('SELECT * FROM Projects').fetchall()
    conn.close()
    
    return project


def get_cases(project_name, database):
    '''
    (str, str) -> list
    
    Returns a list of dictionaries with case information
    
    Paramaters
    -----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT case_id, assay, donor_id, ext_id, species, miso FROM Samples WHERE project_id = ?", (project_name,)).fetchall()
    conn.close()
    
    data = [dict(i) for i in data]
         
    return data



def extract_samples_libraries_per_case(project_name, database):
    '''
    (str, str) - > dict
    
    Returns a dictionary with samples sorted by tissue type and with libraries sorted by library type
        
    Parameters
    ----------
    - project_name (str): Name of project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    data = conn.execute("SELECT DISTINCT case_id, donor_id, tissue_type, tissue_origin, library_type, group_id, library FROM Libraries WHERE project_id = ?;", (project_name,)).fetchall()
    conn.close()

    D = {}
    
    for i in data:
        case = i['case_id']
        donor = i['donor_id']
        tissue_type = i['tissue_type']
        library_type = i['library_type']
        library = i['library']
        
        sample = '_'.join([donor, i['tissue_origin'] , tissue_type, library_type, i['group_id']])
        if tissue_type == 'R':
            tissue = 'normal'
        else:
            tissue = 'tumor'
        
        if case not in D:
            D[case] = {}
        if 'samples' not in D[case]:
            D[case]['samples'] = {}
        if 'libraries' not in D[case]:
            D[case]['libraries'] = {}
        
        if tissue not in D[case]['samples']:
            D[case]['samples'][tissue] = set()
        if library_type not in D[case]['libraries']:
            D[case]['libraries'][library_type] = set()
            
        D[case]['samples'][tissue].add(sample)
        D[case]['libraries'][library_type].add(library)
        
    return D            



def get_last_sequencing(project_name, database):
    '''
    (str, str) -> str
    
    Returns the date of the last sequencing for the project of interest
    
    Paramaters
    ----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    conn = connect_to_db(database)
    sequencing = conn.execute("SELECT DISTINCT Files.creation_date FROM Files JOIN Workflows \
                              WHERE Files.project_id = '{0}' AND Workflows.project_id = '{0}' \
                              AND Workflows.wfrun_id = Files.wfrun_id AND LOWER(Workflows.wf) in \
                              ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq');".format(project_name)).fetchall()
    conn.close()
    
    # get the most recent creation date of fastq generating workflows
    if sequencing:
        seq_dates = sorted([i['creation_date'] for i in sequencing])
        most_recent = seq_dates[-1]
    else:
        most_recent = 'NA'
        
    try:
        most_recent = convert_epoch_time(most_recent)    
        return most_recent
    except:
        return most_recent
    


def get_case_analysis_status(analysis_database, project_name=None):
    '''
    (str, str) -> dict
    
    Returns a dictionary with the analysis status of each case in a project if
    project name is specified or all projects otherwise.
    If a case has multiple analysis templates, it will return the status of a complete
    template if one exists
        
    Parameters
    ----------
    - analysis_database (str): Path to the database storing the analysis data
    - project_name (None str): Name of a specific project
    '''
    
    conn = connect_to_db(analysis_database)
    if project_name:
        data = conn.execute("SELECT project_id, case_id, valid FROM templates WHERE project_id = ?", (project_name,)).fetchall()
    else:
        data = conn.execute("SELECT project_id, case_id, valid FROM templates").fetchall()
    conn.close()
    
    D = {}
    for i in data:
        project = i['project_id']
        case = i['case_id']
        valid = i['valid']
        if project not in D:
            D[project] = {}
        if case in D[project]:
            D[project][case].append(int(valid))
        else:
            D[project][case] = [int(valid)]
    
    # return the status of the complete template if one exists
    for project in D:
        for case in D[project]:
            D[project][case] = sorted(D[project][case])
            D[project][case] = D[project][case][-1]
        
    return D


def count_completed_cases(analysis_status):
    '''
    (dict) -> dict
    
    Returns a dictionary with counts of cases with complete
    and incomplete analysis for each project     
       
    Parameters
    ----------
    - analysis_status (dict): Dictionary with analysis status of each case
                              for each project
    '''
    
    D = {}
    
    for project in analysis_status:
        # get the status of all cases
        status = list(analysis_status[project].values())
        complete = sum(status)
        incomplete = len(status) - complete
        
        D[project] = {'complete': complete, 'incomplete': incomplete}
    
    return D
    
    
def get_case_sequencing_status(database, project_name=None):
    '''
    (str, str) -> dict
    
    Returns a dictionary with the sequencing status of each case in a project if
    project name is specified or all projects otherwise.
            
    Parameters
    ----------
    - database (str): Path to the waterzooi database
    - project_name (None str): Name of a specific project
    '''
    
    conn = connect_to_db(database)
    if project_name:
        data = conn.execute("SELECT project_id, case_id, sequencing_status FROM Samples WHERE project_id = ?", (project_name,)).fetchall()
    else:
        data = conn.execute("SELECT project_id, case_id, sequencing_status FROM Samples").fetchall()
    conn.close()
    
    D = {}
    for i in data:
        project = i['project_id']
        case = i['case_id']
        sequencing_status = i['sequencing_status']
        if project not in D:
            D[project] = {}
        assert case not in D[project]
        D[project][case] = int(sequencing_status)
        
    return D



def count_complete_sequencing(sequencing_status):
    '''
    (dict) -> dict
    
    Returns a dictionary with counts of cases with complete and incomplete
    sequencing for each project     
       
    Parameters
    ----------
    - sequencing_status (dict): Dictionary with sequencing status of each case
                                for each project
    '''
    
    D = {}
    
    for project in sequencing_status:
        # get the status of all cases
        status = list(sequencing_status[project].values())
        complete = sum(status)
        incomplete = len(status) - complete
        
        D[project] = {'complete': complete, 'incomplete': incomplete}
    
    return D

