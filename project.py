# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 21:11:22 2023

@author: rjovelin
"""


from utilities import connect_to_db


def get_project_info(project_name, database):
    '''
    (str, str) -> list
    
    Returns a list with project information extracted from database for project_name 
    
    Parameters
    ----------
    - project_name (str): Project of interest
    - database (str): Path to the sqlite database
    '''
    
    # connect to db
    conn = connect_to_db(database)
    # extract project info
    project = conn.execute('SELECT * FROM Projects WHERE project_id=?', (project_name,)).fetchall()[0]
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
    data = conn.execute("SELECT DISTINCT case_id, donor_id, ext_id, species, sex, miso FROM Samples WHERE project_id = ?", (project_name,)).fetchall()
    data = [dict(i) for i in data]
    for i in range(len(data)):
        assay = data[i]['case_id'].split(':')[1]
        data[i]['assay'] = assay
      
    return data


# def get_sample_counts(project_name, database):
#     '''
#     (str, str) - > dict
    
#     Returns a dictionary with sample counts for each donor of a project of interest
    
#     Parameters
#     ----------
#     - project_name (str): Name of project of interest
#     - database (str): Path to the sqlite database
#     '''
    
#     conn = connect_to_db(database)
#     data = conn.execute("SELECT DISTINCT case_id, tissue_type, tissue_origin, library_type, group_id FROM Libraries WHERE project_id = ?;", (project_name,)).fetchall()
#     conn.close()

#     counts = {}
#     for i in data:
#         donor = i['case_id']
#         sample = '_'.join([i['case_id'], i['tissue_type'], i['tissue_origin'],
#                            i['library_type'], i['group_id']])
#         if i['tissue_type'] == 'R':
#             tissue = 'normal'
#         else:
#             tissue = 'tumor'
#         if donor not in counts:
#             counts[donor] = {}
#         if tissue not in counts[donor]:
#             counts[donor][tissue] = set()
#         counts[donor][tissue].add(sample)
        
#     for i in counts:
#         for j in ['normal', 'tumor']:
#             if j in counts[i]:
#                 counts[i][j] = len(counts[i][j])
#             else:
#                 counts[i][j] = 0
                    
#     return counts            


# def get_library_types(project_name, database):
#     '''
#     (str, str) -> list
    
#     Returns a list of different library types for a given project
    
#     Parameters
#     ----------
#     - project_name (str): Name of project of interest
#     - database (str): Path to the sqlite database
#     '''
    
#     # connect to db
#     conn = connect_to_db(database)
#     # extract library types
#     data = conn.execute("SELECT DISTINCT library_types FROM Projects WHERE project_id = ?;", (project_name,)).fetchall()
#     conn.close()
    
#     library_types = sorted(list(map(lambda x: x.strip(), data[0]['library_types'].split(','))))
   
#     return library_types


# def count_libraries(project_name, library_types, cases, database):
#     '''
#     (str, list, list, str) -> dict
    
#     Returns a dictionary with libraries for each library type and sample for a given project
       
#     Parameters
#     ----------
#     - project_name (str) Name of the project of interest
#     - library_types (list): List of library types recorded for project
#     - cases (list): List of dictionary with case information  
#     - database (str): Path to the sqlite database
#     '''
    
#     # connect to db
#     conn = connect_to_db(database)
#     # extract library source
#     data = conn.execute("SELECT DISTINCT case_id, library_type, library FROM Libraries WHERE project_id = '{0}';".format(project_name)).fetchall()
#     conn.close()
    
#     libraries= {}
    
#     # initiate the libraries dict
#     for i in cases:
#         libraries[i['case_id']] = {}
#         for j in library_types:
#             libraries[i['case_id']][j] = set()
    
#     # record libraries for each library type
#     for i in data:
#         libraries[i['case_id']][i['library_type']].add(i['library'])
    
#     return libraries



#######################

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




#########################

















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
    sequencing = conn.execute("SELECT DISTINCT Workflow_Inputs.run FROM Workflow_Inputs JOIN Files \
                              WHERE Workflow_Inputs.project_id = '{0}' AND Files.project_id = '{0}' \
                              AND Files.wfrun_id = Workflow_Inputs.wfrun_id \
                              AND LOWER(Files.workflow) in ('casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq');".format(project_name)).fetchall()
    conn.close()
    
    # get the most recent creation date of fastq generating workflows
    
    if sequencing:
        sequencing = list(set(sequencing))
        sequencing = [i['run'] for i in sequencing]
        sequencing = map(lambda x: x.split('_'), sequencing)
        seq_dates = [i for i in sequencing if any(list(map(lambda x: x.isdigit(), i)))]
        
        F = lambda y: list(map(lambda x: x.isdigit(), y)).index(True)
        date_index = list(map(lambda x: F(x), seq_dates))
        seq_dates = [seq_dates[i][date_index[i]] for i in range(len(date_index))]
        seq_date = sorted(list(set(seq_dates)))[-1]
        seq_date = '20' + str(seq_date)[:2] + '-' + str(seq_date)[2:4] + '-' + str(seq_date)[4:]
        
    else:
        seq_date = 'NA'
    return seq_date
