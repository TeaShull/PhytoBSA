#!/usr/bin/env python

from flask import Flask, render_template, request, session
from thaleBSA_utilities import ThaleBSAUtilities  # Import ThaleBSAUtilities

app = Flask(__name__)
app.secret_key = '1111'

# Initialize variables
my_species = ""
reference_genome_source = ""
known_snps = ""
threads_limit = 0
cleanup = False

# Get the number of threads available on the machine
import os
import multiprocessing

available_threads = multiprocessing.cpu_count()
threads_limit = max(1, available_threads // 2)  # Set half the available threads as the default limit

# Create an instance of ThaleBSAUtilities
thale_bsa_utils = ThaleBSAUtilities()

@app.route('/', methods=['GET', 'POST'])
def index():
    global my_species, reference_genome_source, known_snps, threads_limit, cleanup

    if request.method == 'POST':
        my_species = request.form['species']
        reference_genome_source = request.form['reference_source']
        known_snps = request.form['known_snps']
        threads_limit = int(request.form['threads_limit'])
        cleanup = 'cleanup' in request.form

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
    thale_bsa_utils.vcf_file_generation()  # Use the class method

@app.route('/run_data_analysis', methods=['POST'])
def run_data_analysis():
    thale_bsa_utils.data_analysis()  # Use the class method

if __name__ == '__main__':
    app.run(debug=True)