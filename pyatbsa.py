#!/usr/bin/env python

from flask import Flask, render_template, request, session

import subprocess

from pyatbsa_functions import *

app = Flask(__name__)
app.secret_key = '1111'

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run_create_experiment_dictionary', methods=['POST'])
def run_create_experiment_dictionary():
    experiment_dictionary = create_experiment_dictionary()
    session['experiment_dictionary'] = experiment_dictionary
    return render_template('index.html', message = experiment_dictionary)

@app.route('/run_vcf_file_generation', methods=['GET'])
def run_vcf_file_generation():
     vcf_file_generation()

@app.route('/run_data_analysis', methods=['GET'])
def run_data_analysis():
    data_analysis()

if __name__ == '__main__':
    app.run(debug=True)




