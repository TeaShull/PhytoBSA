#!/usr/bin/env python
from settings.config import INPUT_DIR, OUTPUT_DIR, LOG_DIR
import os
import sys
import multiprocessing
from datetime import datetime
import argparse

from modules.core import ThaleBSAParentFunctions
from modules.utilities_general import FileUtilities, DictionaryUtilities
from modules.utilities_logging import LogHandler

def parse_program_arguments():
    parser = argparse.ArgumentParser(description='PyAtBSA main command line script...')

    ## Analysis inputs
    parser.add_argument('-an', '--analysis', action='store_true', help='Run the analysis.')
    parser.add_argument('-n', '--line_name', type=str, help='name of the line you wish to analyze. Will be used to name output files.')
    parser.add_argument('-vt', '--vcf_table', type=str, help='path to the vcf table you wish to analyze.')
    parser.add_argument('-st', '--segregation_type', type=str, help="Recessive(R) or Dominant(D)")

    ## VCF generation inputs
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
    return parser.parse_args()
def main():
    core_log = LogHandler('core')
    core_log.note(f'Core log begin. ulid: {core_log.ulid}')
    core_log.add_db_record()
    
    # Create instances of ThaleBSAParentFunctions and FileUtilities. 
    file_utils = FileUtilities(core_log)

    # parse command line arguments
    line_name = args.line_name
    vcf_table = args.vcf_table
    reference_genome_name = args.reference_genome_name
    snpEff_species_db = args.snpEff_species_db
    reference_genome_source = args.reference_genome_source
    threads_limit = args.threads_limit
    cleanup = args.cleanup
    known_snps = args.known_snps
    segregation_type = args.segregation_type
    

    # Check if user wants to the command line and variables.py instead of the Flask app.
    core_log.attempt('Parsing command line arguments...')
    try:
        # [if -an arg] accept line name and vcf table to run bsa_analysis 
        if args.analysis:
            core_log.note('Command line argument to run analysis detected.')
            core_log.attempt(f'Trying to create experiment_dictionary from arguments...')
            try:
                if line_name and vcf_table and segregation_type:
                    experiment_dictionary=parent_functions.populate_experiment_dictionary(
                        line_name, vcf_table, allele
                    )
                    core_log.success(f'experiment_dictionary successfully created')
                else: 
                    core_log.fail('-ln, -vt or -al not set. Aborting...')
            
            except Exception as e:
                core_log.fail('There was a failure trying to create experiment_dictionary from passed arguments:{e}')
            
            core_log.attempt(f'Attempting to begin analysis of {args.vcf_table}...')
            try: 
                parent_functions.bsa_analysis(experiment_dictionary)
                quit()
            except Exception as e:
                core_log.fail(f'There was an error while trying to start bsa_analysis:{e}')

        if not args.analysis:
            core_log.attempt('No command line arguments given. Running automatic command line operations.')
            dict_utils = DictionaryUtilities(core_log) #Initiate dict_utils
            dict_utils.experiment_detector() # Detect experiments in ./inputs
            dict_utils.check_vcf_gen_variables() #Source variables from settings/vcf_gen_variables.py if not passed by user
            dict_utils.generate_output_file_paths() #Generate neccissary output paths before running pipeline
            
            parent_functions = ThaleBSAParentFunctions(core_log, dict_utils)
            parent_functions.vcf_generation() #Generate VCF file to analyze
            parent_functions.bsa_analysis() #Analyze VCF file
            quit()

    except Exception as e:
        core_log.fail(f'Starting thaleBSA has failed: {e}')
        quit()

if __name__ == '__main__':
    sys.exit(main())
