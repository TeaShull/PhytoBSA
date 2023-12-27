#!/usr/bin/env python
import os
import sys
import multiprocessing
from datetime import datetime
from flask import Flask, render_template, request, session

modules_dir = os.path.join(os.getcwd(),'src','modules')
sys.path.append(modules_dir)

from thaleBSA_core import ThaleBSAParentFunctions, ThaleBSAUtilities
from config import (
    error_handler, SRC_DIR, LOG_DIR, OUTPUT_DIR, TEMPLATE_DIR, STATIC_DIR
)

#Initialize directory structure variables
src_dir = SRC_DIR
template_dir = TEMPLATE_DIR
static_dir = STATIC_DIR

# Initialize variables
reference_genome_name = ""
snpEff_db_name = ""
reference_genome_source = ""
threads_limit = 0
cleanup = True
known_snps = ""

# Get the number of threads available on the machine
available_threads = multiprocessing.cpu_count()
threads_limit = max(1, available_threads)  # Set half the available threads as the default limit

# Create instances of ThaleBSAUtilities
thale_bsa_parent_functions = ThaleBSAParentFunctions()
thale_bsa_utils = ThaleBSAUtilities()

#Check if user wants to just use the command line and variables.py instead of flask app. 
try: 
    if len(sys.argv) > 1:
        # Check if the first command line argument is set to 'cl'
        if sys.argv[1] == '-cl':
            error_handler('success', 'Command line argument is set to [-cl]. Running on command line.')
            error_handler('attempt', "Sourcing variables from variables.py")
            from variables import *
            
            experiment_dictionary = thale_bsa_utils.create_experiment_dictionary()

            thale_bsa_parent_functions.vcf_generation(
            experiment_dictionary, reference_genome_name, snpEff_db_name, 
            reference_genome_source, threads_limit, cleanup, known_snps
            )
            
            thale_bsa_parent_functions.bsa_analysis(experiment_dictionary)
            quit()
            else:
                error_handler('success', "Command line argument is not set to 'cl'. Starting flask app...")
        else:
            error_handler('success', "No command line arguments provided. Starting flask app...")
except Exception as e:
    error_handler('fail', 'Starting thaleBSA has failed: {e}')
    quit()

app = Flask(__name__, template_folder=template_dir, static_folder=static_dir)
app.secret_key = '1111'

@app.route('/', methods=['GET', 'POST'])
def index():
    global reference_genome_name, snpEff_db_name, reference_genome_source 
    global threads_limit, cleanup, known_snps


    if request.method == 'POST':
        reference_genome_name = request.form['reference_genome_name']
        snpEff_db_name = request.form['snpEff_db_name']
        reference_genome_source = request.form['reference_source']
        threads_limit = int(request.form['threads_limit'])
        cleanup = 'cleanup' in request.form
        known_snps = request.form['known_snps']


        # Check if the selected threads limit is equal to the maximum available threads
        if threads_limit == available_threads:
            warning_message = (
                "Are you sure you want to use all the resources available? "
                "This may cause performance issues."
            )

    return render_template(
        'index.html',
        #reference_genome_name=reference_genome_name,
        reference_genome_source=reference_genome_source,
        known_snps=known_snps,
        threads_limit=threads_limit,
        cleanup=cleanup,
        available_threads=available_threads
    )

@app.route('/run_create_experiment_dictionary', methods=['POST'])
def run_create_experiment_dictionary():
    error_handler('trigger', 'Create experiment dictionary triggered')
    try:
        experiment_dictionary = thale_bsa_utils.create_experiment_dictionary()  # Use the class method
        session['experiment_dictionary'] = experiment_dictionary
        return render_template(
            'index.html',
            message=experiment_dictionary,
            available_threads=available_threads
        )
    except:
        return render_template(
        'index.html',
        message='Flask trigger for building experiment dictionary seems to have failed.',
        available_threads=available_threads
        )

@app.route('/run_vcf_file_generation', methods=['POST'])
def run_vcf_file_generation():
    error_handler('trigger', 'VCF file generation triggered.')
    #try:
    experiment_dictionary = session.get('experiment_dictionary', {})
    thale_bsa_utils.vcf_file_generation(experiment_dictionary, reference_genome_name,
        reference_genome_source, threads_limit, cleanup, known_snps
    )  # Use the class method
    return render_template('index.html', 
        message="VCF file generation started.", 
        available_threads=available_threads
    )
    # except:
    #     return render_template(
    #         'index.html',
    #         message='Flask trigger for VCF file generation seems to have failed.',
    #         available_threads=available_threads
    #     )

@app.route('/run_data_analysis', methods=['POST'])
def run_data_analysis():
    error_handler('trigger', 'Run data analysis triggered')
    #try:
    experiment_dictionary = session.get('experiment_dictionary', {})
    thale_bsa_utils.data_analysis(experiment_dictionary)  # Use the class method
    return render_template('index.html', 
        message="Data analysis started.", 
        available_threads=available_threads
    )
    # #except Exception as e:
    #      return render_template('index.html', 
    #         message="Flask trigger for running analysis seems to have failed.", 
    #         available_threads=available_threads
    #     )

#Initialize Flask App
error_handler('attempt', 'Starting ThaleBSA Flask App')
try:
    if __name__ == '__main__':
        app.run(debug=True)
except Exception as e:
    error_handler('fail', 'Starting flask app has failed: {e}')