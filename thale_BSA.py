#!/usr/bin/env python
import os
import sys
import multiprocessing
from datetime import datetime
import argparse

modules_dir = os.path.join(os.getcwd(), 'src', 'modules')
sys.path.append(modules_dir)

from core import ThaleBSAParentFunctions
from utilities_general import FileUtilities
from utilities_logging import LogHandler

from config import (
    SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, TEMPLATE_DIR, STATIC_DIR
)

# Initialize core_log instance of LogHandler
'''
[NOTE]
All classes in this program must have a log passed to them upon initialization.

Every log is assigned a unique ID (ulid) upon initialization of the LogHandler class
ulid's are passed around in the LogHandler class instance.
the ulid is used to link all file outputs to its relevent log. 

Upon creating a new instance of a class, point the log output to an 
appropriate logger. 
If an appropriate logger doesn't exist in the log list, simply initialize one 
with an a good name, and add it to the list below to keep track. 

Current log list:
'core' - logs all parent functions, thale_bsa.py and flask front end. 
'vcf_gen' - logs all messages pertaining to parent_functions.vcf_generation
'analysis' - logs all messages peratining to parent_functions.bsa_analysis 
'''

def parse_command_line_arguments():
    # Argument parsing
    parser = argparse.ArgumentParser(description='PyAtBSA main command line script...')
    ## Command line inputs
    parser.add_argument('-cl', '--command_line', action='store_true', help='Run on the command line.')
    
    ## Analysis inputs
    parser.add_argument('-an', '--analysis', action='store_true', help='Run the analysis.')
    parser.add_argument('-n', '--line_name', type=str, help='name of the line you wish to analyze. Will be used to name output files.')
    parser.add_argument('-vt', '--vcf_table', type=str, help='path to the vcf table you wish to analyze.')

    reference_genome_name_help = f"""
        What is the name of your reference genome? this should be the base name of your fasta file. 
        example - Arabidopsis_thaliana.fa 
        reference_genome_name = Arabidopsis_thaliana""" 
    parser.add_argument('-rgn', '--reference_genome_name', default=None,
                        type=str, help=reference_genome_name_help)
    
    snpEff_species_db_help = ("What is the name of your snpEff database for your reference genome?")
    parser.add_argument('-ssdb', '--snpEff_species_db', default=None, type=str)
    
    parser.add_argument('-rgs', '--reference_genome_source', default=None, type=str)
    parser.add_argument('-t', '--threads_limit', default=None, type=str)
    parser.add_argument('-c','--cleanup', default=None, type=str)
    parser.add_argument('-ks','--known_snps', default=None, type=str)
    args = parser.parse_args()
    return args

def main():
    # initialize core log
    core_log = LogHandler('core')
    core_log.note(f'Core log begin. ulid: {core_log.ulid}')
    core_log.add_db_record()
    
    # parse command line arguments
    args = parse_command_line_arguments()
    line_name = args.line_name
    vcf_table = args.vcf_table
    reference_genome_name = args.reference_genome_name
    snpEff_species_db = args.snpEff_species_db
    reference_genome_source = args.reference_genome_source
    threads_limit = args.threads_limit
    cleanup = args.cleanup
    known_snps = args.known_snps
    
    # Create instances of ThaleBSAParentFunctions and FileUtilities. 
    parent_functions = ThaleBSAParentFunctions(core_log)
    file_utils = FileUtilities(core_log)

    # Check if user wants to the command line and variables.py instead of the Flask app.
    core_log.attempt('Parsing command line arguments...')
    try:
        # [If -cl arg] detected, attempt to run program in automatic mode. 
        if args.command_line:
            core_log.note('Command line argument is set to [-cl]. Running automatic command line operations.')
            experiment_dictionary = parent_functions.vcf_generation()
            parent_functions.bsa_analysis(experiment_dictionary)
            quit()
        else:
            core_log.note("Command line argument is not set.")

        # [if -an arg] accept line name and vcf table to run bsa_analysis 
        if args.analysis:
            core_log.note('Command line argument to run analysis detected.')

            core_log.attempt(f'Trying to create experiment_dictionary from arguments...')
            try:
                if line_name and vcf_table:
                    experiment_dictionary=file_utils.create_experiment_dictionary(
                        line_name, vcf_table
                    )
                    core_log.success(f'experiment_dictionary successfully created')
                else: 
                    core_log.fail('-ln and -vt not set. Aborting...')
            
            except Exception as e:
                core_log.fail('There was a failure trying to create experiment_dictionary from passed arguments:{e}')
            
            core_log.attempt(f'Attempting to begin analysis of {args.vcf_table}...')
            try: 
                parent_functions.bsa_analysis(experiment_dictionary)
                quit()
            except Exception as e:
                core_log.fail(f'There was an error while trying to start bsa_analysis:{e}')

        if not args.command_line and args.analysis:
            core_log.attempt('No command line arguments given. Starting Flask app....')

    except Exception as e:
        core_log.fail(f'Starting thaleBSA has failed: {e}')
        quit()

if __name__ == "__main__":
    main()
