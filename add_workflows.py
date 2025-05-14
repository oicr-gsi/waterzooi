# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 15:36:41 2024

@author: rjovelin
"""


import os
import sqlite3
from utilities import connect_to_db
import argparse


def collect_workflows(db, project_name, table = 'Workflows'):
    '''
    (str, str) -> list
    
    Returns a list of workflows for project_name extracted from db database
    
    Parameters
    ----------
    - db (str): Path to the sqlite database
    - project_name (str): Name of project of interest
    - table (str): Table in db storing the workflow information. Default is Workflows
    '''
    
    conn = connect_to_db(db)
    data = conn.execute('SELECT wfrun_id FROM {0} WHERE {0}.project_id = "{1}"'.format(table, project_name)).fetchall()
    conn.close()
    
    L = [i['wfrun_id'] for i in data]
    
    return L
    

def collect_projects(main_db, table = 'Projects'):
    '''
    (str, str) -> list
    
    Returns a list of project names extracted from the Projects table in main database
    
    Parameters
    ----------
    - main_db (str): Path to the sqlite database storing production data
    - table (str): Table from main_db storing workflow information. Default is Projects 
    '''
    
    conn = connect_to_db(main_db)
    data = conn.execute('SELECT project_id FROM {0}'.format(table)).fetchall()
    conn.close()
    
    L = [i['project_id'] for i in data]
    
    return L


def add_missing_workflow(workflow_db, project_id, workflows, workflow_table='Workflows'):
    '''
    (str, str, list, str) -> None
    
    Add all workflows into workflow_table of worfklow_db. Set selected status to 0. 

    Parameters
    ----------
    - workflow_db (str): Path to the sqlite database storing Workflow information
    - project_id (str): Name of project of interest
    - workflows (list): List of workflow ids not alread recorded in the workflow_db
    - workflow_table : Table in db storing the workflow information. Default is Workflows
    '''

    column_names = ('wfrun_id', 'project_id', 'selected')
    
    conn = connect_to_db(workflow_db)
    for workflow_id in workflows:
        # insert data into table
        conn.execute('INSERT INTO {0} {1} VALUES {2}'.format(workflow_table, column_names, (workflow_id, project_id, 0)))
        conn.commit()
    conn.close()


def update_workflow_db(workflow_db, main_db, project_table, workflow_table):
    '''
    (str, str, str, str) -> None
    
    Update the workflow database with the workflows from the main database if these
    workflows are not already recorded
        
    Parameters
    ----------
    - workflow_db (str): Path to the sqlite database storing Workflow information
    - main_db (str): Path to the sqlite database storing production data
    - project_table (str): Table from main_db storing workflow information. Default is Projects 
    - workflow_table : Table in db storing the workflow information. Default is Workflows
    '''

    # make a list of project name
    projects = collect_projects(main_db, project_table)
    
    for project in projects:
        # collect the workflows from the main database
        new_workflows = collect_workflows(main_db, project, workflow_table)
        # collect the workflows from the workflow database
        recorded_workflows = collect_workflows(workflow_db, project, workflow_table)
        # make a list of workflows to be added
        workflows = list(set(new_workflows).difference(set(recorded_workflows)))
        # add the missing workflows for the project of focus
        add_missing_workflow(workflow_db, project, workflows, workflow_table)
        
        

def setup_database(database, table = 'Workflows'):
    '''
    (str, str) -> None
    
    Create a database with table Workflows
    
    Parameters
    ----------
    - database (str): Path to the sqlite database
    - table (str): Table in database. Defulat is Workflows
    '''
        
    # create dict to store column names for each table {table: [column names]}
    column_names = ['wfrun_id', 'project_id', 'selected']
    column_types = ['VARCHAR(572)', 'VARCHAR(128)', 'INT']
    
    # define table format including constraints    
    table_format = ', '.join(list(map(lambda x: ' '.join(x), list(zip(column_names, column_types)))))

    # connect to database
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    # create table
    cmd = 'CREATE TABLE {0} ({1})'.format(table, table_format)
    cur.execute(cmd)
    conn.commit()
    conn.close()





if os.path.isfile('workflows_case.db') == False:
    setup_database('workflows_case.db', 'Workflows')
update_workflow_db('workflows_case.db', 'waterzooi_db_case.db', 'Projects', 'Workflows')


# if __name__ == '__main__':

#     parser = argparse.ArgumentParser(prog = 'add_workflows.py', description='Updates the Workflows database with workflow ids from the main database', add_help=True)
#     parser.add_argument('-w', '--workflow_db', dest='workflow_db', help='Path to the database storing workflow information')
#     parser.add_argument('-m', '--main_db', dest = 'main_db', help = 'Path to the main database storing production information')
        
#     # get arguments from the command line
#     args = parser.parse_args()
#     # create database if doesn't exist
#     if os.path.isfile(args.workflow_db) == False:
#         setup_database(args.workflow_db, 'Workflows')
#     # add missing workflows
#     update_workflow_db(args.workflow_db, args.main_db, 'Projects', 'Workflows')
    
    