#!/usr/bin/env python
from settings.config import INPUT_DIR, OUTPUT_DIR, LOG_DIR
import os
import sys
import multiprocessing
from datetime import datetime
import argparse

from modules.core import ThaleBSAParentFunctions
from modules.utilities_general import FileUtilities
from modules.utilities_logging import LogHandler

def parse_program_arguments():
    parser = argparse.ArgumentParser(description='PyAtBSA main command line script...')
    parser.add_argument('-cl', '--command_line', action='store_true', help='Run on the command line.')
    parser.add_argument('-an', '--analysis', action='store_true', help='Run the analysis.')
    parser.add_argument('-n', '--line_name', type=str, help='name of the line you wish to analyze. Will be used to name output files.')
    parser.add_argument('-vt', '--vcf_table', type=str, help='path to the vcf table you wish to analyze.')
    parser.add_argument('-st', '--segregation_type', type=str, help='Recessive(R) or Dominant(D)?')
    parser.add_argument('-rgn', '--reference_genome_name', default=None, type=str, help='What is the name of your reference genome? this should be the base name of your fasta file.')
    parser.add_argument('-ssdb', '--snpEff_species_db', default=None, type=str, help="What is the name of your snpEff database for your reference genome?")
    parser.add_argument('-rgs', '--reference_genome_source', default=None, type=str)
    parser.add_argument('-t', '--threads_limit', default=None, type=str)
    parser.add_argument('-c','--cleanup', default=None, type=str)
    parser.add_argument('-ks','--known_snps', default=None, type=str)
    return parser.parse_args()
def main():
    core_log = LogHandler('core')
    core_log.note(f'Core log begin. ulid: {core_log.ulid}')
    core_log.add_db_record()
    
    args = parse_program_arguments()
    parent_functions = ThaleBSAParentFunctions(core_log)
    file_utils = FileUtilities(core_log)
    try:
        if args.command_line:
            core_log.note('Command line argument is set to [-cl]. Running automatic command line operations.')
            experiment_dictionary = parent_functions.vcf_generation()
            parent_functions.bsa_analysis(experiment_dictionary)
            quit()
        if args.analysis:
            core_log.note('Command line argument to run analysis detected.')
            if args.line_name and args.vcf_table and args.segregation_type:
                experiment_dictionary=file_utils.create_experiment_dictionary(
                args.line_name, args.segregation_type, args.vcf_table
                )
                core_log.success(f'experiment_dictionary successfully created')
                parent_functions.bsa_analysis(experiment_dictionary)
                quit()
        if not args.command_line and args.analysis:
            core_log.attempt('No command line arguments given. Starting Flask app....')

    except Exception as e:
        core_log.fail(f'Starting thaleBSA has failed: {e}')
        quit()

if __name__ == '__main__':
    sys.exit(main())
