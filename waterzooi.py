# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:32:40 2022

@author: rjovelin
"""

#import sqlite3
import os
import json
from flask import Flask, render_template, request, url_for, flash, redirect, Response, send_file
#from werkzeug.exceptions import abort
import time
import pandas as pd
import matplotlib
matplotlib.use('agg')
from utilities import connect_to_db, get_library_design, secret_key_generator, get_case_md5sums, \
    extract_case_signoff, extract_nabu_signoff, list_signoff_deliverables, remove_cases_with_no_approval_signoff, \
    remove_cases_with_competed_cbioportal_release, remove_workflows_with_deliverable_signoff, \
    get_workflow_release_status, get_file_release_status 
from whole_genome import get_amount_data, get_workflows_analysis_date, get_workflow_file_count, \
    get_selected_workflows, update_wf_selection, get_input_sequences, get_cases_with_analysis,\
    get_case_analysis_samples, get_case_analysis_workflows, count_case_analysis_workflows,\
    most_recent_analysis_workflow, get_analysis_workflow_name, get_case_workflow_samples, \
    get_assays, get_missing_workflows, get_case_parent_to_children_workflows, \
    get_case_children_to_parents_workflows, get_case_workflow_info,\
    get_workflow_output_files, delete_cases_with_distinct_checksums,\
    map_donors_to_cases, list_assay_analysis_workflows, \
    get_sequencing_input, get_case_error_message, create_analysis_json, \
    get_workflow_outputfiles, get_pipeline_standard_deliverables,\
    create_case_analysis_json, get_review_status, identify_deliverables, \
    create_cbioportal_json, get_workflow_names, list_template_workflows, \
    create_graph_edges, plot_graph
from project import get_project_info, get_cases, get_last_sequencing, extract_samples_libraries_per_case, \
    get_case_analysis_status, count_completed_cases, get_case_sequencing_status, count_complete_sequencing
from sequencing import collect_sequence_info, get_platform_shortname

import plotly.offline as pyo
import plotly.graph_objs as go



app = Flask(__name__)
#app.config['SECRET_KEY'] = secret_key_generator(10)
app.secret_key = secret_key_generator(10)



# database = 'waterzooi_db_case.db'
# workflow_db = 'workflows_case.db'
# analysis_db = 'analysis_review_case.db'
# nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'


database = 'waterzooi_db_case_small.db'
workflow_db = 'workflows_case_small.db'
analysis_db = 'analysis_review_case_small.db'
nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'





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
    
    # analysis_database = 'analysis_review_case.db'
    # database = 'waterzooi_db_case.db'
    
    # extract project info
    projects = get_project_info(database)
    projects = sorted([(i['project_id'], i) for i in projects])
    projects = [i[1] for i in projects]
    
    # get analysis status of each case in each project
    analysis_status = get_case_analysis_status(analysis_db)
    # count complete and incomplete cases
    analysis_counts = count_completed_cases(analysis_status)  
    # get the sequencing status of each case
    sequencing_status = get_case_sequencing_status(database)
    # count cases with complete and incomplete sequencing
    sequencing_counts = count_complete_sequencing(sequencing_status)  
      
    return render_template('index.html',
                           projects=projects,
                           analysis_counts=analysis_counts,
                           sequencing_counts=sequencing_counts)


@app.route('/<project_name>')
def project_page(project_name):
    
    analysis_database = 'analysis_review_case.db'
    nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'
    
    
    #database = 'waterzooi.db'
    database = 'waterzooi_db_case.db'
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    # get case information
    cases = get_cases(project_name, database)
    # sort by case id
    cases = sorted(cases, key=lambda d: d['case_id']) 
    
    
    print(cases)
    
    # get signoffs
    case_names = [d['case_id'] for d in cases]
    signoffs = extract_nabu_signoff(case_names, nabu_key_file)
    # list the release deliverables for each case
    deliv = list_signoff_deliverables(signoffs)
    
    print(signoffs)
    
    # get the species
    species = ', '.join(sorted(list(set([i['species'] for i in cases]))))
    # get the assays
    assay_names = get_assays(database, project_name)
    assays = sorted(list(set(assay_names.split(','))))
    # get the samples and libraries for each case respectively sorted by tissue and library type
    samples_libraries = extract_samples_libraries_per_case(project_name, database)
    library_types = sorted(list(map(lambda x: x.strip(), project['library_types'].split(','))))
    seq_date = get_last_sequencing(project['project_id'], database)
    # get library definitions
    library_names = {i: get_library_design(i) for i in library_types}
    # get the analysis status of each case
    analysis_status = get_case_analysis_status(analysis_database, project_name)
    # count complete and incomplete cases
    analysis_counts = count_completed_cases(analysis_status)      
    # get the sequencing status of each case
    sequencing_status = get_case_sequencing_status(database, project_name)
    # count cases with complete and incomplete sequencing
    sequencing_counts = count_complete_sequencing(sequencing_status)
    
    return render_template('project.html', project=project, cases=cases,
                           assays=assays, assay_names = assay_names,
                           samples_libraries = samples_libraries,
                           seq_date=seq_date, species=species, 
                           library_types = library_types,
                           library_names=library_names,
                           analysis_status=analysis_status,
                           analysis_counts=analysis_counts,
                           sequencing_status=sequencing_status,
                           sequencing_counts=sequencing_counts,
                           signoffs=signoffs,
                           deliv=deliv
                           )
    

@app.route('/<project_name>/sequencing', methods = ['GET', 'POST'])
def sequencing(project_name):
    
    #database = 'waterzooi.db'
    # database = 'waterzooi_db_case.db'
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    # get sequence file information
    sequences = collect_sequence_info(project_name, database)       
    # get the assays
    assay_names = get_assays(database, project_name)
    assays = sorted(list(set(assay_names.split(','))))
    # map the instrument short name to sequencing platform
    platform_names = get_platform_shortname(project_name, database)
 
    if request.method == 'POST':
        
        platforms = request.form.getlist('platform')
                
        D = {}
        for i in sequences:
            d = {'Case': i['case'],
                 'Donor': i['donor'],
                 'DonorID': i['sample'],
                 'SampleID': i['group_id'],
                 'Sample': i['sample_id'],
                 'Description': i['group_description'],
                 'Library': i['library'],
                 'Library Type': i['library_type'],
                 'Tissue Origin': i['tissue_origin'],
                 'Tissue Type': i['tissue_type'],
                 'File Prefix': i['prefix']}
    
            # download all information if platforms are not selected
            if platforms:
                # check that platform is selected
                if platform_names[i['platform']] in platforms:
                    D[i['case']] = d     
            else:
                D[i['case']] = d
                    
        data = pd.DataFrame(D.values())
        outputfile = '{0}_libraries.xlsx'.format(project_name)
        data.to_excel(outputfile, index=False)
            
        return send_file(outputfile, as_attachment=True)

    else:
        return render_template('sequencing.html', project=project,
                               sequences=sequences, assays=assays,
                               platform_names=platform_names
                               )



@app.route('/<project_name>/<assay>', methods=['POST', 'GET'])
def analysis(project_name, assay):
    
    
    print('project_name',project_name)
    print('assay', assay)
    
    
    
    #database = 'waterzooi_db_test.db'
    # database = 'waterzooi_db_case.db'
    
    # workflow_db = 'workflows_case.db'
    # analysis_database = 'analysis_review_case.db'
    # nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'
    
    
    assay = assay.replace('+:+', '/')
        
    print('assay_modif', assay)
    
    
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    
    
    print('project', project)
    
    
    # get the deliverables
    deliverables = identify_deliverables(project)
    print('deliverables', deliverables)
    # get the cases with analysis data for that project and assay
    cases = get_cases_with_analysis(analysis_db, project_name, assay)
    
    print('cases', len(cases))
    print(cases.keys())
    
    # get signoffs
    signoffs = extract_nabu_signoff(cases, nabu_key_file)
    # list the release deliverables for each case
    deliv = list_signoff_deliverables(signoffs)


    
    
    print(project['deliverables'])
    
    
    
    # check that analysis is up to date with the waterzooi database
    md5sums = get_case_md5sums(database, project_name)
    # keep only cases with up to date data between resources
    delete_cases_with_distinct_checksums(cases, md5sums)
    # get the donor
    donors = map_donors_to_cases(cases)
    # get the assays
    assay_names = get_assays(database, project_name)
    assays = sorted(list(set(assay_names.split(','))))
    # get the samples analyzed in the assay for each case
    samples = get_case_analysis_samples(cases)
    case_names = sorted(list(cases.keys()))
    # re-organized the workflows
    case_data = get_case_analysis_workflows(cases)
    # count workflows
    workflow_counts = count_case_analysis_workflows(case_data)
    # get all the analysis workflows across each case of the same assay
    analysis_workflows = list_assay_analysis_workflows(case_data)
    # get the analysis status of each case
    analysis_status = get_case_analysis_status(analysis_db, project_name)
    if project_name in analysis_status:
        analysis_status = analysis_status[project_name]
    else:
        analysis_status = {}
    # get the combined error messages across template for each assay
    errors = get_case_error_message(cases)
    # get the sequencing status of each case
    sequencing_status = get_case_sequencing_status(database, project_name)
    # get the selected status of each workflows
    selected_workflows = get_selected_workflows(project_name, workflow_db, 'Workflows')    
    # get the review status of each case
    review_status = get_review_status(case_data, selected_workflows)
        
    if request.method == 'POST':
        deliverable = request.form.get('deliverable')
        # get the workflow output files
        workflow_outputfiles = get_workflow_outputfiles(database, project_name)
        if deliverable == 'selected':
            analysis_data = create_analysis_json(case_data, selected_workflows, workflow_outputfiles)
            filename = '{0}.pipeline.json'.format(project_name)
        elif deliverable == 'standard':
            standard_deliverables = get_pipeline_standard_deliverables()
            analysis_data = create_analysis_json(case_data, selected_workflows, workflow_outputfiles, standard_deliverables)
            filename = '{0}.pipeline.standard.json'.format(project_name)
        elif deliverable in ['purple', 'sequenza']:
            selected_workflows = get_selected_workflows(project_name, workflow_db, 'Workflows')
            # create json with workflow information for cbioportal importer
            analysis_data = create_cbioportal_json(cases, selected_workflows, workflow_outputfiles, deliverable)
            filename = '{0}.{1}.cbioportal.json'.format(project_name, assay)
        else:
            analysis_data = {}
            filename = '{0}.pipeline.json'.format(project_name)
        
        # keep only cases with proper signoff (completed release approval and deliverable not signed off)
        if analysis_data:
            analysis_data = remove_cases_with_no_approval_signoff(analysis_data, signoffs)
            if deliverable in ['selected', 'standard']:
                # remove workflows part of deliverables with complete signoff
                analysis_data = remove_workflows_with_deliverable_signoff(analysis_data, signoffs, deliverable, 'pipeline')
                analysis_data = remove_workflows_with_deliverable_signoff(analysis_data, signoffs, deliverable, 'fastq')
            elif deliverable in ['sequenza', 'purple']:
                # remove cases for which cbioportal release is signed off
                analysis_data = remove_cases_with_competed_cbioportal_release(analysis_data, signoffs, deliverable)
        
        return Response(
            response=json.dumps(analysis_data),
            mimetype="application/json",
            status=200,
            headers={"Content-disposition": "attachment; filename={0}".format(filename)})

    else:
        return render_template('assay.html',
                           project=project,
                           assays=assays,
                           current_assay = assay,
                           samples=samples,
                           cases=cases,
                           donors=donors,
                           case_names=case_names,
                           case_data=case_data,
                           analysis_workflows=analysis_workflows,
                           workflow_counts=workflow_counts,
                           analysis_status=analysis_status,
                           sequencing_status=sequencing_status,
                           errors=errors,
                           review_status=review_status,
                           deliverables=deliverables,
                           signoffs=signoffs,
                           deliv=deliv
                           )

@app.route('/<project_name>/<assay>/<case>/', methods = ['POST', 'GET'])
def case_analysis(project_name, assay, case):
    
    print('assay_case_analysis', assay)
    
    #### delete cases with distinct md5sums
    
    print(case)
    print(project_name)
    
    #database = 'waterzooi_db_test.db'
    # database = 'waterzooi_db_case.db'
    
    # workflow_db = 'workflows_case.db'
    # analysis_db = 'analysis_review_case.db'
    # nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'
    
    assay = assay.replace('+:+', '/')
    case = case.replace('+:+', '/')
    
    
    
    
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    deliverables = identify_deliverables(project)
    print(deliverables)
    
    
    # get the cases with analysis data for that project and assay
    cases = get_cases_with_analysis(analysis_db, project_name, assay)
    print('retrieved the cases')
    
    cases = {case: cases[case]}
    
    # get the case sign off
    case_signoffs = extract_case_signoff(case, nabu_key_file)
    # list the release deliverables for case
    deliv = list_signoff_deliverables(case_signoffs) 
       
    # get the release status of each workflow
    workflow_qc = get_workflow_release_status(database, case)
    
    
    print('get the workflow release status')    

    
    # check that analysis is up to date with the waterzooi database
    md5sums = get_case_md5sums(database, project_name)
    # keep only cases with up to date data between resources
    delete_cases_with_distinct_checksums(cases, md5sums)
    # get analysis data
    case_analysis = cases[case]
    case_data = get_case_analysis_workflows(cases) 
    case_data = case_data[case]
    # get the creation date of all workflows in each template
    creation_dates = get_workflows_analysis_date(project_name, database)
    # get the most recent creation date for each template
    most_recent = most_recent_analysis_workflow(case_data, creation_dates)
    # get the workflow names    
    workflow_names = [get_analysis_workflow_name(case_data[i]['analysis']) for i in range(len(case_data))]
    
    
    print('get workflow names')
    
    # get the sequencing status of the case
    sequencing_status = get_case_sequencing_status(database, project_name)
    sequencing_status = sequencing_status[project_name][case]
    # get the file count of each workflow in project
    file_counts = get_workflow_file_count(project_name, database)
    
    
    print('file counts')
    
    # get the amount of data for each workflow
    amount_data = get_amount_data(project_name, database)
    # get the samples corresponding to each worklow id
    samples = get_case_workflow_samples(database, case)
    
    print('got samples')
    
    # get the assays
    assay_names = get_assays(database, project_name)
    assays = sorted(list(set(assay_names.split(','))))
    
    print('got assays')
    
    
    # map limskeys and libraries to sequencing workflows
    seq_inputs = get_sequencing_input(database, case)
    
    print('got se inputs')
    
    
    # list expected workflows without workflow ids
    missing_workflows = get_missing_workflows(case_data)
    # get workflow relationships
    parent_to_children = get_case_parent_to_children_workflows(database, case)
    child_to_parents = get_case_children_to_parents_workflows(parent_to_children)
    
    print('got workflow relatinships')
    
    
    # extract selected status of each workflow
    selected = get_selected_workflows(project_name, workflow_db, 'Workflows')
    
    
    print('done all that')
    
    
    # map the workflow id to their names
    workflow_full_names = get_workflow_names(database, case)
    # plot the workflow relationships for each case template 
    figures = []
    for template in cases[case]:
        # list all the workflow ids of that template
        workflow_ids = list_template_workflows(template)
        # create the graph edges
        edges = create_graph_edges(workflow_ids, parent_to_children)
        # create a figure
        fig = plot_graph(edges, workflow_full_names)
        # create the html plot
        plot_html = pyo.plot(fig, output_type='div', include_plotlyjs='cdn')
        figures.append(plot_html)
   
    if request.method == 'POST':
        # get the list of checked workflows        
        selected_workflows = request.form.getlist('workflow')
        # make a list of workflows across all templates
        workflows = []
        for template in case_data:
            for i in ['sequencing', 'analysis', 'alignments']:
                for workflow in template[i]:
                    workflows.extend(template[i][workflow])
        update_wf_selection(workflows, selected_workflows, selected, workflow_db, 'Workflows')
        return redirect(url_for('case_analysis', project_name=project_name, assay=assay.replace('/', '+:+'), case=case))
    
    else:
        return render_template('case_assay.html',
                           project=project,
                           assays=assays,
                           assay=assay,
                           case=case,
                           case_data=case_data,
                           case_analysis=case_analysis,
                           most_recent=most_recent,
                           workflow_names=workflow_names,
                           file_counts=file_counts,
                           amount_data=amount_data,
                           creation_dates=creation_dates,
                           samples=samples,
                           selected=selected,
                           missing_workflows=missing_workflows,
                           child_to_parents=child_to_parents,
                           sequencing_status=sequencing_status,
                           seq_inputs=seq_inputs,
                           deliverables=deliverables,
                           figures=figures,
                           case_signoffs=case_signoffs,
                           deliv=deliv,
                           workflow_qc=workflow_qc
                           )






@app.route('/<project_name>/<assay>/<case>/<path:wfrunid>')
def show_workflow(project_name, assay, case, wfrunid):
    
    # database = 'waterzooi_db_case.db'
      
    assay = assay.replace('+:+', '/')
    case = case.replace('+:+', '/')
    wfrunid = wfrunid.replace('+:+', '/')
    
    
    # get the release status of each workflow
    workflow_qc = get_workflow_release_status(database, case)
    # get the file release status
    file_qc = get_file_release_status(database, case)
        
    
      
    # get the project info for project_name from db
    project = get_project_info(database, project_name)[0]
    # get the parent and children workflows
    parent_to_children = get_case_parent_to_children_workflows(database, case)
    child_to_parents = get_case_children_to_parents_workflows(parent_to_children)
    # get workflow name and version
    workflow_info = get_case_workflow_info(database, case)
    # get the output files
    outputfiles = get_workflow_output_files(database, wfrunid)
    # get the input sequences
    input_sequences = get_input_sequences(database, case, wfrunid)
        
    return render_template('workflow_info.html',
                       project=project,
                       case=case,
                       assay=assay,
                       workflow_info=workflow_info,
                       workflow_id=wfrunid,
                       child_to_parents=child_to_parents,
                       parent_to_children=parent_to_children,
                       outputfiles=outputfiles,
                       input_sequences=input_sequences,
                       workflow_qc=workflow_qc,
                       file_qc=file_qc
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
    
    # database = 'waterzooi_db_case.db'
    # analysis_database = 'analysis_review_case.db'
    
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
    # get the sequencing status of each case
    sequencing_status = get_case_sequencing_status(database, project_name)

    D = {}
    for i in cases:
        case = i['case_id']
        D[case] = {'donor': i['donor_id'], 'external_id': i['ext_id'], 'assay': i['assay']}
        if sequencing_status[project['project_id']][i['case_id']]:
            seq_status = 'complete'
        else:
            seq_status = 'incomplete'
        D[case]['sequencing_status'] = seq_status
        if analysis_status[project['project_id']][i['case_id']]:
            data_status = 'complete'
        else:
            data_status = 'incomplete'
        D[case]['analysis_status'] = data_status
        if 'normal' in samples_libraries[i['case_id']]['samples']:
            normal_count = len(samples_libraries[i['case_id']]['samples']['normal'])
        else:
            normal_count = 0
        D[case]['normal'] = normal_count
        if 'tumor' in samples_libraries[i['case_id']]['samples']:
            tumor_count = len(samples_libraries[i['case_id']]['samples']['tumor'])
        else:
            tumor_count = 0
        D[case]['tumor'] = tumor_count
        for library_type in library_types:
            if library_type in samples_libraries[i['case_id']]['libraries']:
                D[case][library_type] = len(samples_libraries[i['case_id']]['libraries'][library_type])
            else:
                D[case][library_type] = 0
        
        
    data = pd.DataFrame(D.values())
    data.to_excel('{0}_cases.xlsx'.format(project_name), index=False)
   
    return send_file("{0}_cases.xlsx".format(project_name), as_attachment=True)



@app.route('/download_analysis/<project_name>/<case>/<assay>/<selection>')
def download_analysis_data(project_name, case, assay, selection):
        
    # database = 'waterzooi_db_case.db'
    # workflow_db = 'workflows_case.db'
    # analysis_db = 'analysis_review_case.db'
    # nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'
    
    
    
    # get the case sign off
    case_signoffs = extract_case_signoff(case, nabu_key_file)
        
    assay = assay.replace('+:+', '/')
    case = case.replace('+:+', '/')
    
   
    # get the cases with analysis data for that project and assay
    cases = get_cases_with_analysis(analysis_db, project_name, assay)
    case_data = get_case_analysis_workflows(cases) 
    case_data = case_data[case]
    # extract selected status of each workflow
    selected_workflows = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the workflow output files
    workflow_outputfiles = get_workflow_outputfiles(database, project_name)
    # create json with workflow information for DARE
    analysis_data = create_case_analysis_json(case, case_data, selected_workflows, workflow_outputfiles, selection)
        
    # keep only data with proper signoff (completed release approval and deliverable not signed off)
    if analysis_data:
        analysis_data = remove_cases_with_no_approval_signoff(analysis_data, case_signoffs)
        # remove workflows part of deliverables with complete signoff
        analysis_data = remove_workflows_with_deliverable_signoff(analysis_data, case_signoffs, selection, 'pipeline')
        analysis_data = remove_workflows_with_deliverable_signoff(analysis_data, case_signoffs, selection, 'fastq')
        
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(analysis_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}.{1}.{2}.{3}.json".format(case, project_name, assay, selection)})


@app.route('/download_cbioportal/<project_name>/<case>/<assay>/<segmentation>')
def download_cbioportal_data(project_name, case, assay, segmentation):
        
    # database = 'waterzooi_db_case.db'
    # workflow_db = 'workflows_case.db'
    # analysis_db = 'analysis_review_case.db'
    # nabu_key_file = 'nabu-prod_qc-gate-etl_api-key'
    
    
    # get the case sign off
    case_signoffs = extract_case_signoff(case, nabu_key_file)
        
    assay = assay.replace('+:+', '/')
    case = case.replace('+:+', '/')
    
    # get the cases with analysis data for that project and assay
    cases = get_cases_with_analysis(analysis_db, project_name, assay)
    cases = {case: cases[case]}
    
    # extract selected status of each workflow
    selected_workflows = get_selected_workflows(project_name, workflow_db, 'Workflows')
    # get the workflow output files
    workflow_outputfiles = get_workflow_outputfiles(database, project_name)
    # create json with workflow information for cbioportal importer
    analysis_data = create_cbioportal_json(cases, selected_workflows, workflow_outputfiles, segmentation)
        
    # keep only data with proper signoff (completed release approval and deliverable not signed off)
    if analysis_data:
        analysis_data = remove_cases_with_no_approval_signoff(analysis_data, case_signoffs)
        # remove cases for which cbioportal release is signed off
        analysis_data = remove_cases_with_competed_cbioportal_release(analysis_data, case_signoffs, segmentation)
    
    # send the json to outoutfile                    
    return Response(
        response=json.dumps(analysis_data),
        mimetype="application/json",
        status=200,
        headers={"Content-disposition": "attachment; filename={0}.{1}.{2}.cbioportal.json".format(case, project_name, assay)})





# if __name__ == "__main__":
#     app.run(host='0.0.0.0')
    