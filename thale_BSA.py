#!/usr/bin/env python
import os
import sys
import multiprocessing
from flask import Flask, render_template, request, session

modules_dir = os.path.join(os.getcwd(),'src','modules')
sys.path.append(modules_dir)

from thaleBSA_utilities import ThaleBSAUtilities
import config

#Initialize directory structure variables
src_dir = config.SRC_DIR
template_dir = config.TEMPLATE_DIR
static_dir = config.STATIC_DIR

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

# Create an instance of ThaleBSAUtilities
thale_bsa_utils = ThaleBSAUtilities(current_line_name)

#Check if user wants to just use the command line and variables.py instead of flask app. 
if len(sys.argv) > 1:
    # Check if the first command line argument is set to 'cl'
    if sys.argv[1] == 'cl':
        command_line = True
        print("Command line argument is set to 'cl'. Running on command line.")
        print("Sourcing variables from variables.py")
        from variables import *
        experiment_dictionary = thale_bsa_utils.create_experiment_dictionary()

        #thale_bsa_utils.vcf_file_generation(experiment_dictionary, 
        #reference_genome_name, snpEff_db_name, reference_genome_source, 
        #threads_limit, cleanup, known_snps
        #)
        thale_bsa_utils.data_analysis(experiment_dictionary, command_line)
        quit()
    else:
        command_line = False
        print("Command line argument is not set to 'cl'. Starting flask app...")
else:
    command_line = False
    # Your logic when no command line arguments are provided
    print("No command line arguments provided. Starting flask app...")

#Initialize Flask App
app = Flask(__name__, template_folder=template_dir, static_folder=static_dir)
app.secret_key = '1111'

@app.route('/', methods=['GET', 'POST'])
def index():
    global my_species, reference_genome_source, known_snps, threads_limit, cleanup

    if request.method == 'POST':
        my_species = request.form['species']
        reference_genome_source = request.form['reference_source']
        known_snps = request.form['known_snps']
        threads_limit = int(request.form['threads_limit'])
        cleanup = 'cleanup' in request.form
        snpEff_species_db = request.form['snpEff_species_db']

        # Check if the selected threads limit is equal to the maximum available threads
        if threads_limit == available_threads:
            warning_message = (
                "Are you sure you want to use all the resources available? "
                "This may cause performance issues."
            )

    return render_template(
        'index.html',
        my_species=my_species,
        reference_genome_source=reference_genome_source,
        known_snps=known_snps,
        threads_limit=threads_limit,
        cleanup=cleanup,
        available_threads=available_threads
    )

@app.route('/run_create_experiment_dictionary', methods=['POST'])
def run_create_experiment_dictionary():
    experiment_dictionary = thale_bsa_utils.create_experiment_dictionary()  # Use the class method
    session['experiment_dictionary'] = experiment_dictionary
    return render_template(
        'index.html',
        message=experiment_dictionary,
        available_threads=available_threads
    )

@app.route('/run_vcf_file_generation', methods=['POST'])
def run_vcf_file_generation():
    thale_bsa_utils.vcf_file_generation(experiment_dictionary, species,
        reference_genome_source, threads_limit, cleanup, known_snps
    )  # Use the class method
    return render_template('index.html', message="VCF file generation started.", available_threads=available_threads)

@app.route('/run_data_analysis', methods=['POST'])
def run_data_analysis():
    thale_bsa_utils.data_analysis()  # Use the class method
    return render_template('index.html', message="Data analysis started.", available_threads=available_threads)


if __name__ == '__main__':
    app.run(debug=True)