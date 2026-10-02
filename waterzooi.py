# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:32:40 2022

@author: rjovelin
"""

#import sqlite3
import os
import json
from flask import Flask, render_template, request, url_for, flash, redirect, Response, send_file
from werkzeug.exceptions import abort
import time
import pandas as pd


from waterzooi_helper import secret_key_generator, get_project_info, \
    get_project_level_deliverables, get_release_signoff, get_case_analysis_status, \
    count_completed_cases, extract_samples_libraries_per_case, collect_sequence_info, \
    get_fileqc, merge_qc_status_workflow, get_assays, get_platform_shortname, get_cases, \
    get_analysis_data, get_analysis_samples, get_analysis_workflows, error_formatting, \
    get_case_analysis_data, get_case_release_signoff, get_workflows_analysis_date, \
    most_recent_analysis_workflow, map_workflows_to_fileids, map_analysis_workflows, \
    organize_data, get_case_assay, get_case_parent_to_children_workflows, \
    get_case_children_to_parents_workflows, get_workflow_output_files, get_case_workflow_info, \
    get_input_sequences, add_workflow_qc_status, get_sequences_to_download, \
    get_data_release_signoff, get_data_release_approval_signoff, get_output_files, \
    get_workflow_outputs, prepare_analysis_json, prepare_cbioportal_json, count_cases, \
    plot_graph, plot_small_graph, get_last_sequencing, rename_case_id, create_graph_edges, \
    get_library_design    
        

import plotly.offline as pyo
import plotly.graph_objs as go
import plotly.io as pio


app = Flask(__name__)
#app.config['SECRET_KEY'] = secret_key_generator(10)
app.secret_key = secret_key_generator(10)



#database = 'waterzooi_db_case.db'
workflow_db = 'workflows_case.db'
#analysis_db = 'analysis_review_case.db'
nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'

#database = 'waterzooi_test_09092026.db'

database = 'waterzooi_test_10012026.db'



#analysis_db = 'analysis_review_test_09102026.db'
#analysis_db = 'analysis_review_test_09302026.db'
analysis_db = 'analysis_review_test_10022026.db'



database = 'waterzooi_test_09092026.db'
nabu_cache = 'nabu_cache.db'

workflow_deliv = 'workflow_deliverables.json'





@app.template_filter()
def find_workflow_id(generic_name, bmpp_children_workflows, library):
    '''
    (str, str, dict, str) -> str
    
    Returns the workflow id of a workflow that has generic name as substring and
    NA if no workflow has generic name as substring.
            
    Parameters
    ----------
    - generic_name (str): Generic workflow name, may be substring of workflow name in bmpp_children_workflows
    - bmpp_children_workflows (dict): Dictionary with downstream workflow information
    - library (str): Libraries of interest
    '''
    
    # make a list of downstream workflows for that bmpp run
    L = list(bmpp_children_workflows[library].keys())
    # create a same size list of generic name workflows
    workflows = [generic_name] * len(L)
    
    # define function to identify generic workflow as subtring of workflow name
    is_workflow = lambda x,y: x.split('_')[0].lower() == y.lower()
    
    # check if generic workflow is substring of bmpp children workflows 
    found = list(map(is_workflow, L, workflows))
    if any(found):
        return bmpp_children_workflows[library][L[found.index(True)]]['workflow_id']
    else:
        return 'NA'
    



@app.template_filter()
def readable_time(date):
    '''
    (str) -> str
    
    Returns epoch time in readable format
    
    Parameters
    ----------
    - date (str): Epoch time
    '''
    
    #return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(date)))
    return time.strftime('%Y-%m-%d', time.localtime(int(date)))


@app.template_filter()
def shorten_workflow_id(workflow_run_id):
    '''
    (str) -> str
    
    Shorten the workflow run id to 8 characters + trailing dots
             
    Parameters
    ----------
    - workflow_run_id (str): Workflow unique run identifier
    '''
    
    return workflow_run_id[:8] + '...'



@app.template_filter()
def format_identifier(identifier):
    '''
    (str) -> str
    
    Format case and assay for url
             
    Parameters
    ----------
    - identifier (str): case or assay identifier
    '''
    
    return identifier.replace('/', '+:+')




@app.template_filter()
def basename(workflow_run_id):
    '''
    (str) -> str
    
    Returns the basename of workflow id
             
    Parameters
    ----------
    - workflow_run_id (str): Workflow unique run identifier
    '''
    
    return os.path.basename(workflow_run_id)




@app.template_filter()
def format_created_time(created_time):
    '''
    (str) -> str
    
    Remove time in created time and keep only the date
                 
    Parameters
    ----------
    - created_time (str): Time a sample is created in the format Year-Month-DayTHour:Mn:SecZ
    '''
    
    return created_time[:created_time.index('T')]
    
    
@app.route('/')
def index():
    
    # extract project info and sort by project name
    projects = get_project_info(database)
    projects.sort(key = lambda d: d['project_id'])
    # get analysis status of each case in each project
    analysis_status = get_case_analysis_status(analysis_db)
    # count complete and incomplete cases
    analysis_counts = count_completed_cases(analysis_status)  
             
    return render_template('index.html',
                           projects=projects,
                           analysis_counts=analysis_counts)


@app.route('/<project_name>')
def project_page(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    # get case information
    cases = get_cases(project_name, database)
    # sort by case id
    cases = sorted(cases, key=lambda d: d['case_id']) 
    # extract signoff from the nabu cache
    signoffs = get_release_signoff(nabu_cache, project_name)
    # get the assays
    assays = get_assays(database, project_name)
    # get the samples and libraries for each case respectively sorted by tissue and library type
    samples_libraries = extract_samples_libraries_per_case(project_name, database)
    library_types = sorted(list(map(lambda x: x.strip(), project['library_types'].split(','))))
    seq_date = get_last_sequencing(project['project_id'], database)
    # get library definitions
    library_names = {i: get_library_design(i) for i in library_types}
    # get the analysis status of each case
    analysis_status = get_case_analysis_status(analysis_db, project_name)
    # count complete and incomplete cases
    analysis_counts = count_completed_cases(analysis_status)      
    
    return render_template('project.html',
                           project=project,
                           cases=cases,
                           assays=assays,
                           samples_libraries = samples_libraries,
                           seq_date=seq_date,
                           library_types = library_types,
                           library_names=library_names,
                           analysis_status=analysis_status,
                           analysis_counts=analysis_counts,
                           signoffs=signoffs
                           )
    

@app.route('/<project_name>/sequencing', methods = ['GET', 'POST'])
def sequencing(project_name):
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    # get sequence file information
    sequences = collect_sequence_info(project_name, database)       
    # get file qc for all files in project
    fileqc = get_fileqc(nabu_cache, project_name)
    # add release status of each workflow
    sequences = add_workflow_qc_status(sequences, fileqc)
    # get the assays
    assays = get_assays(database, project_name)
    # map the instrument short name to sequencing platform
    platform_names = get_platform_shortname(project_name, database)
 
    if request.method == 'POST':
        # get the user-selected sequencing platforms
        platforms = request.form.getlist('platform')
        # collect the sequencing information for the selected platforms        
        L = get_sequences_to_download(sequences, platform_names, platforms)
        # save to Excel file
        data = pd.DataFrame(L)
        outputfile = '{0}_libraries.xlsx'.format(project_name)
        data.to_excel(outputfile, index=False)
            
        return send_file(outputfile, as_attachment=True)

    else:
        return render_template('sequencing.html',
                               project=project,
                               sequences=sequences,
                               assays=assays,
                               platform_names=platform_names
                               )



@app.route('/<project_name>/<assay>', methods=['POST', 'GET'])
def analysis(project_name, assay):
    
    assay = assay.replace('+:+', '/')
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    # get analysis data
    analysis_data = get_analysis_data(analysis_db, project_name, assay)
    # sort cases id
    case_names = sorted(list(analysis_data.keys()))
          
    # get the samples from analysis data for each case
    samples = get_analysis_samples(analysis_data)
    
    # get the analysis workflows
    workflows = get_analysis_workflows(analysis_data)
    
    # extract signoff from the nabu cache
    signoffs = get_release_signoff(nabu_cache, project_name)
    # get project delievrables
    project_deliverables = project['deliverables'].split(',')
    # check if data release approval is signed off for each case
    data_release_approval = {case_id: get_data_release_approval_signoff(signoffs, case_id)[case_id] for case_id in case_names}
    # check if data release is already signed off
    data_release = {case_id: get_data_release_signoff(signoffs, project_deliverables, case_id, 'all')[case_id] for case_id in case_names}
    # check if cbioportal signoff exists or complete
    cbio_signoff = {case_id: get_data_release_signoff(signoffs, project_deliverables, case_id, 'cbioportal')[case_id] for case_id in case_names}
    # check if pipeline release exists or complete
    pipeline_signoff = {case_id: get_data_release_signoff(signoffs, project_deliverables, case_id, 'pipeline')[case_id] for case_id in case_names}
    
    # count cases
    complete, incomplete, complete_signedoff, complete_to_release, \
    pipeline_to_release, cbio_to_release = count_cases(analysis_data, data_release_approval,
                                                       data_release, pipeline_signoff,
                                                       cbio_signoff)
    
    # get the assays
    assays = get_assays(database, project_name)
       
    
    # add formatted error message
    for case_id in analysis_data:
        if analysis_data[case_id]['error']:
            err = error_formatting(analysis_data[case_id]['error'])
            analysis_data[case_id]['error_message'] = err
        
    if request.method == 'POST':
        deliverable = request.form.get('deliverable')
        
        # get the output files of each workflow for all cases 
        outputs = get_workflow_outputs(database, project_name)
        # keep only cases with complete data, data release appoval signoff and no release signoff
        analyses, workflow_outputs = {}, {}
        for case_id in analysis_data:
            if analysis_data[case_id]['valid'] and data_release_approval[case_id] and pipeline_signoff[case_id] == False:
                analyses[case_id] = analysis_data[case_id]
                workflow_outputs[case_id] = outputs[case_id]
    
        if deliverable == 'pipeline':
            # get pipeline deliverables
            infile = open(workflow_deliv)
            workflow_deliverables = json.load(infile)
            infile.close()
            
            # organize data for download
            downloadable_data = prepare_analysis_json(analyses, workflow_outputs, workflow_deliverables)
            
        else:
            # organize data for download
            downloadable_data = prepare_analysis_json(analyses, workflow_outputs)
                
        # send the json to outoutfile                    
        return Response(
            response=json.dumps(downloadable_data),
            mimetype="application/json",
            status=200,
            headers={"Content-disposition": "attachment; filename={0}.{1}.json".format(project_name, assay.replace(' ', '_'))})

    else:
        return render_template('assay.html',
                           project = project,
                           assays = assays,
                           current_assay = assay,
                           analysis_data = analysis_data,
                           project_deliverables = project_deliverables,
                           samples = samples,
                           workflows = workflows,
                           case_names = case_names,
                           signoffs=signoffs,
                           complete=complete,
                           incomplete=incomplete,
                           complete_signedoff=complete_signedoff,
                           complete_to_release=complete_to_release,
                           pipeline_to_release = pipeline_to_release,
                           cbio_to_release = cbio_to_release
                           )
                           

@app.route('/<project_name>/<assay>/<case_id>/',  methods=['POST', 'GET'])
def case_analysis(project_name, assay, case_id):
    
    assay = assay.replace('+:+', '/')
    case_id = case_id.replace('+:+', '/')
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    
    # get analysis data for case
    analysis_data = get_case_analysis_data(analysis_db, case_id, project_name, assay)

    # format error message
    if analysis_data[case_id]['error']:
        error = error_formatting(analysis_data[case_id]['error'])
    else:
        error = {}
        
    # extract signoff from the nabu cache
    signoffs = get_case_release_signoff(nabu_cache, case_id, project_name)
    
    # get the creation date of all workflows in each template
    creation_dates = get_workflows_analysis_date(case_id, project_name, database)
    # get the most recent workflow creation date
    most_recent = most_recent_analysis_workflow(analysis_data, case_id, creation_dates)
    
    # organize the data for easier sorting
    data, analysis_workflows, parent_to_children = organize_data(analysis_data, case_id)
    workflow_names = sorted(list(analysis_workflows.values()))
    workflow_runs = sorted(list(analysis_workflows.keys()))

    # get the release status for each workflow run id
    # map file swids to workflow runids
    workflow_outputs = map_workflows_to_fileids(case_id, project_name, database)
    # get the file qc of each file for the case
    fileqc = get_fileqc(nabu_cache, project_name, case_id)
    workflow_release_status = {wfrunid: merge_qc_status_workflow(workflow_outputs[wfrunid], fileqc) for wfrunid in workflow_outputs}

    # get the assays
    assays = get_assays(database, project_name)
        
    # check if analysis data validation
    valid = analysis_data[case_id]['valid']
        
    # get project delievrables
    project_deliverables = project['deliverables'].split(',')
    # get deliverables signoffs
    data_release_approval = get_data_release_approval_signoff(signoffs, case_id)
    # check if data release is already signed off
    data_release = get_data_release_signoff(signoffs, project_deliverables, case_id, 'all')
    # check if cbioportal signoff exists or complete
    cbio_signoff = get_data_release_signoff(signoffs, project_deliverables, case_id, 'cbioportal')
    # check if pipeline release exists or complete
    pipeline_signoff = get_data_release_signoff(signoffs, project_deliverables, case_id, 'pipeline')
    
    # create the graph edges
    edges = create_graph_edges(workflow_runs, parent_to_children)
    # create a figure
    fig = plot_graph(edges, analysis_workflows)
    # create the html plot
    plot_html = pyo.plot(fig, output_type='div', include_plotlyjs='cdn')
   
    # determine the number of columns in summary table
    table_cols = 7
    if 'fastq' in project['deliverables'].lower() or 'pipeline' in project['deliverables'].lower():
        table_cols += 1
    if 'cbioportal' in project['deliverables'].lower():
        table_cols += 1
       
   
    if request.method == 'POST':
        deliverable = request.form.get('deliverable')
                        
        # get the output files of each workflow for the case 
        workflow_outputfiles = get_workflow_outputs(database, project_name, case_id)
        # keep only cases with complete data, data release appoval signoff and no release signoff
        analyses, outputfiles = {}, {}
        if analysis_data[case_id]['valid'] and data_release_approval[case_id] and pipeline_signoff[case_id] == False:
            analyses[case_id] = analysis_data[case_id]
            outputfiles[case_id] = workflow_outputfiles[case_id]
                      
        if deliverable == 'pipeline':
            # get pipeline deliverables
            infile = open(workflow_deliv)
            workflow_deliverables = json.load(infile)
            infile.close()
            # organize data for download
            downloadable_data = prepare_analysis_json(analyses, outputfiles, workflow_deliverables)
        else:
            # organize data for download
            downloadable_data = prepare_analysis_json(analyses, outputfiles)
           
        # replace en dash in file name
        case_name = rename_case_id(case_id)
                   
        # send the json to outoutfile                    
        return Response(
            response=json.dumps(downloadable_data),
            mimetype="application/json",
            status=200,
            headers={"Content-disposition": "attachment; filename={0}.{1}.{2}.json".format(case_name, project_name, assay.replace(' ', '_'))})
        
    
    return render_template('case_assay.html',
                           project=project,
                           assays=assays,
                           assay=assay,
                           case_id=case_id,
                           most_recent=most_recent,
                           workflow_runs = workflow_runs,
                           workflow_release_status = workflow_release_status,
                           workflow_outputs = workflow_outputs,
                           data = data,
                           valid=valid,
                           error=error,
                           creation_dates=creation_dates,
                           signoffs=signoffs,
                           data_release_approval=data_release_approval,
                           data_release=data_release,
                           cbio_signoff=cbio_signoff,
                           pipeline_signoff=pipeline_signoff,
                           plot_html=plot_html,
                           table_cols=table_cols
                           )
                           
                           

@app.route('/workflow_view/<project_name>/<case_id>/<path:wfrunid>')
def show_workflow(project_name, case_id, wfrunid):
    
    case_id = case_id.replace('+:+', '/')
    wfrunid = wfrunid.replace('+:+', '/')
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    
    # get the file qc of each file for the case
    fileqc = get_fileqc(nabu_cache, project_name, case_id)
    # map file swids to workflow runids
    workflow_outputs = map_workflows_to_fileids(case_id, project_name, database)
    workflow_release_status = {wfrunid: merge_qc_status_workflow(workflow_outputs[wfrunid], fileqc) for wfrunid in workflow_outputs}

    # get workflow name and version
    workflow_info = get_case_workflow_info(database, case_id)
    workflow_name = workflow_info[wfrunid][0]
    workflow_version = workflow_info[wfrunid][1]
    
    # get the assay
    assay = get_case_assay(database, project_name, case_id)
    assay = '_'.join(assay.split('_')[:-1])
    
    # get the parent and children workflows
    parent_to_children = get_case_parent_to_children_workflows(database, case_id)
    child_to_parents = get_case_children_to_parents_workflows(parent_to_children)
    
    # get the output files
    outputfiles, files_to_swids = get_workflow_output_files(database, wfrunid)
    
    # get the input sequences
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis',
                            'fileimport', 'import_fastq']    
    if workflow_name.lower() not in sequencing_workflows:
        input_sequences = get_input_sequences(database, case_id, wfrunid)
    else:
        input_sequences = {}
    
    # create the graph edges
    case_workflows = {i: workflow_info[i][0] for i in workflow_info}    
    workflow_runs = sorted(list(case_workflows.keys()))
    
    # make a list of parents, children including workflow of interest
    workflow_runs = []
    parent_runs = []
    children_runs = []
    if wfrunid in parent_to_children and parent_to_children[wfrunid] != ['NA']:
        parent_runs.extend(parent_to_children[wfrunid])
    if wfrunid in child_to_parents and child_to_parents[wfrunid] != ['NA']:
        children_runs.extend(child_to_parents[wfrunid])
    
    workflow_runs.append(wfrunid)
    workflow_runs.extend(parent_runs)
    workflow_runs.extend(children_runs)
    workflow_runs = list(set(workflow_runs))
    
    edges = create_graph_edges(workflow_runs, parent_to_children)
    # create a figure
    fig = plot_small_graph(edges, case_workflows, wfrunid, parent_runs, children_runs)
    # create the html plot
    graph_div = pyo.plot(fig, auto_open=False, output_type='div', include_plotlyjs='cdn') 
    
    return render_template('workflow_info.html',
                       project=project,
                       case_id=case_id,
                       wfrunid=wfrunid,
                       workflow_name=workflow_name,
                       workflow_version=workflow_version,
                       workflow_release_status=workflow_release_status,
                       assay=assay,
                       fileqc=fileqc,
                       workflow_info=workflow_info,
                       child_to_parents=child_to_parents,
                       parent_to_children=parent_to_children,
                       outputfiles=outputfiles,
                       files_to_swids=files_to_swids,
                       input_sequences=input_sequences,
                       graph_div=graph_div
                       )


@app.route('/download_cases/<project_name>')
def download_cases_table(project_name):
    '''
    (str) -> None
    
    Download a table with project information in Excel format
    
    Parameters
    ----------
    - project_name (str): Name of project of interest
    '''
    
    # get case information
    cases = get_cases(project_name, database)
    # sort by case id
    cases = sorted(cases, key=lambda d: d['case_id']) 
    # get the samples and libraries for each case respectively sorted by tissue and library type
    samples_libraries = extract_samples_libraries_per_case(project_name, database)
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    library_types = sorted(list(map(lambda x: x.strip(), project['library_types'].split(','))))
    # get the analysis status of each case
    analysis_status = get_case_analysis_status(analysis_db, project_name)
     

    D = {}
    for i in cases:
        case_id = i['case_id']
        
        assay = '_'.join(i['assay'].split('_')[:-1])
        version = i['assay'].split('_')[-1]
        
        D[case_id] = {'donor': i['donor_id'], 'external_id': i['ext_id'], 'assay': assay, 'version': version}
        if analysis_status[project['project_id']][i['case_id']]:
            data_status = 'complete'
        else:
            data_status = 'incomplete'
        D[case_id]['analysis_status'] = data_status
        if 'normal' in samples_libraries[i['case_id']]['samples']:
            normal_count = len(samples_libraries[i['case_id']]['samples']['normal'])
        else:
            normal_count = 0
        D[case_id]['normal'] = normal_count
        if 'tumor' in samples_libraries[i['case_id']]['samples']:
            tumor_count = len(samples_libraries[i['case_id']]['samples']['tumor'])
        else:
            tumor_count = 0
        D[case_id]['tumor'] = tumor_count
        for library_type in library_types:
            if library_type in samples_libraries[i['case_id']]['libraries']:
                D[case_id][library_type] = len(samples_libraries[i['case_id']]['libraries'][library_type])
            else:
                D[case_id][library_type] = 0
        
        
    data = pd.DataFrame(D.values())
    data.to_excel('{0}_cases.xlsx'.format(project_name), index=False)
   
    return send_file("{0}_cases.xlsx".format(project_name), as_attachment=True)



@app.route('/download_analysis/<project_name>/<case_id>/<assay>')
def download_analysis_data(project_name, case_id, assay):
    
    assay = assay.replace('+:+', '/')
    case_id = case_id.replace('+:+', '/')
    
    # pull down analysis data
    analysis_data = get_case_analysis_data(analysis_db, case_id, project_name, assay)
    # get the output files of each workflow for the case 
    workflow_outputs = get_workflow_outputs(database, project_name, case_id)
    # organize data for download
    downloadable_data = prepare_analysis_json(analysis_data, workflow_outputs)
    
    # replace en dash in file name
    case_name = rename_case_id(case_id)    
    
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(downloadable_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}.{1}.{2}.json".format(case_name, project_name, assay.replace(' ', '_'))})


@app.route('/download_cbioportal/<project_name>/<case_id>/<assay>')
def download_cbioportal_data(project_name, case_id, assay):
 
    assay = assay.replace('+:+', '/')
    case_id = case_id.replace('+:+', '/')
        
    # strip version from assay
    # pull down analysis data
    analysis_data = get_case_analysis_data(analysis_db, case_id, project_name, assay)
    
    # get the output files of each workflow for the case 
    workflow_outputs = get_workflow_outputs(database, project_name, case_id)
    # organize data for cbioportal importer
    downloadable_data = prepare_cbioportal_json(analysis_data, workflow_outputs)
            
    # replace en dash in file name
    case_name = rename_case_id(case_id)   
    
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(downloadable_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}.{1}.{2}.cbioportal.json".format(case_name, project_name, assay)})



@app.route('/download_assay_cbioportal/<project_name>/<assay>')
def download_assay_cbioportal_data(project_name, assay):
 
    assay = assay.replace('+:+', '/')
        
    # pull down analysis data
    analyses = get_analysis_data(analysis_db, project_name, assay)
    # get the output files of each workflow for all cases 
    outputs = get_workflow_outputs(database, project_name)
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    # extract signoff from the nabu cache
    signoffs = get_release_signoff(nabu_cache, project_name)
    # get project delievrables
    project_deliverables = project['deliverables'].split(',')
    # check if data release approval is signed off for each case
    data_release_approval = {case_id: get_data_release_approval_signoff(signoffs, case_id)[case_id] for case_id in analyses}
    # check if cbioportal signoff exists or complete
    cbio_signoff = {case_id: get_data_release_signoff(signoffs, project_deliverables, case_id, 'cbioportal')[case_id] for case_id in analyses}
       
    # keep only cases with complete data, data release appoval signoff and no release signoff
    analysis_data, workflow_outputs = {}, {}
    for case_id in analyses:
        if analyses[case_id]['valid'] and data_release_approval[case_id] and cbio_signoff[case_id] == False:
            analysis_data[case_id] = analyses[case_id]
            workflow_outputs[case_id] = outputs[case_id]
        
    # organize data for cbioportal importer
    downloadable_data = prepare_cbioportal_json(analysis_data, workflow_outputs)
            
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(downloadable_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}.{1}.cbioportal.json".format(project_name, assay)})



# if __name__ == "__main__":
#     app.run(debug=True)


#if __name__ == "__main__":
#    #app.run(host='0.0.0.0', port='8080', debug=True)
#     app.run(host='5000')
#     #app.run(debug=True)