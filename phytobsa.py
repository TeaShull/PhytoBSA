#!/usr/bin/env python
import os
import sys
import argparse

from modules.utilities_variables import (
    Lines,AutomaticLineVariableDetector, VCFGenVariables, BSAVariables
)
from modules.utilities_logging import LogHandler
from modules.core_vcf_gen import VCFGenerator
from modules.core_bsa import BSA

# Initialize core_log instance of LogHandler
'''
All classes in this program must have a log passed to them upon initialization.

Every log is assigned a unique ID (ulid) upon init of LogHandler
ulid's are stored in LogHandler instance, and exported to the Lines data class. 
Ulids link file outputs to their log. 

Logging includes a log database, which will save run information. 

Current log list:
'core' - logs all parent functions, thale_bsa.py and flask front end. 
'vcf_gen' - logs all messages pertaining to parent_functions.vcf_generation
'analysis' - logs all messages peratining to parent_functions.bsa_analysis 
'''
def parse_program_arguments():
    # Main parser
    parser = argparse.ArgumentParser(description='PyAtBSA main command line script...')
    parser.add_argument('-a', '--automatic', action='store_true', help='Run PhytoBSA in automatic mode. Primary, default workflow.')
    subparsers = parser.add_subparsers(help='Sub-command help', dest='command')

    ## Analysis subparser
    parser_analysis = subparsers.add_parser('analysis', help='Run PhytoBSA analysis seperately.')
    ### REQUIRED
    parser_analysis.add_argument('-n', '--name', required=True, type=str, help='name of the line you wish to analyze. Will be used to name output files.')
    parser_analysis.add_argument('-vt', '--vcf_table', required=True, type=str, help='path to the vcf table you wish to analyze.')
    parser_analysis.add_argument('-st', '--segregation_type', required=True, type=str, help="Recessive(R) or Dominant(D)")
    ### OPTIONAL
    parser_analysis.add_argument('-ls', '--loess_span', type=float, default=0.3, help="Influences smoothing parameters.")
    parser_analysis.add_argument('-si', '--shuffle_iterations', type=int, default=1000, help="Iterations of bootstrapping during empirical cutoff calculations. Below 1000 can yeild inconsistant results")
    parser_analysis.add_argument('-sb', '--smooth_edges_bounds', type=int, default=15, help="Number of mirrored datapoints at chromosome edges to correct for loess edge bias. Increase if edge bias seems high")

    ## VCF generation subparser
    parser_vcf_gen = subparsers.add_parser('vcf_generator', help='Run PhytoBSA VCF Generator seperately')
    ### REQUIRED
    parser_vcf_gen.add_argument('-rgn', '--reference_genome_name', required=True, default=None, type=str, help='Reference genome name')
    parser_vcf_gen.add_argument('-wt', '--wt_input', required=True, default = None, type=str, help='Wild-type bulk fasta file(s)')
    parser_vcf_gen.add_argument('-mu', '--mu_input', required=True, default=None, type=str, help='mutant bulk fast files')
    ### OPTIONAL - the below arguments can be input via command line or can be sourced from ./settings/vcf_gen_variables.py
    parser_vcf_gen.add_argument('-rgs', '--reference_genome_source', default=None, type=str, help = 'Optional, if you wish the pipeline to download your reference from a url')
    parser_vcf_gen.add_argument('-ssdb', '--snpEff_species_db', default=None, type=str, help = 'The name of your snpEff database name.')
    parser_vcf_gen.add_argument('-ks','--known_snps', default=None, type=str, help = 'VCF file containing background SNPs. Helps improve output quality')
    parser_vcf_gen.add_argument('-t', '--threads_limit', default=None, type=str, help='Maximum threads you wish to use for analysis')
    parser_vcf_gen.add_argument('-p', '--call_variants_in_parallel', default=False, type=bool, help='Run gatk haplotype caller in parallel')
    parser_vcf_gen.add_argument('-c','--cleanup', default=None, type=bool, help='If true, intermediate files will de deleted. False for troubleshooting and archiving files.' )
    
    ## Log database subparser
    parser_log_db = subparsers.add_parser('logdb', help='Gather log information from log database')
    parser_log_db.add_argument('-vcf', "--vcf_ulid_log", default = None, type = str, help="Retrieve vcf log information based on vcf ulid input")
    parser_log_db.add_argument('-an', "--analysis_ulid_log", default = None, type = str, help="Retrieve analysis log information based on analysis ulid input")
    parser_log_db.add_argument('-name', "--line_name_log", default = None, type = str, help="Retrieve all information associated with the provided line name. ")    
    parser_log_db.add_argument('-core', '--core_ulid_log', default = None, type = str, help="Retrieve all information associated with the provided core log ulid")
    args = parser.parse_args()
    return args


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
        
    # Parse command line arguments
    args = parse_program_arguments()

    # initialize core log, unless trying to pull log information
    if args.command != 'logdb':
        core_log = LogHandler('core')
        core_log.note(f'Core log begin. ulid: {core_log.ulid}')
        core_log.add_db_record()
    
    # Determine which routine to run
    if args.command == 'analysis':
        line = Lines(core_log, args.name)
        
        line.usr_in_line_variables(args.vcf_table, args.segregation_type)
        bsa_vars = BSAVariables(core_log, 
            lines=line, 
            loess_span=args.loess_span, 
            smooth_edges_bounds=args.smooth_edges_bounds, 
            shuffle_iterations=args.shuffle_iterations
        )
    
        bsa = BSA(core_log, bsa_vars)    
        bsa.run_pipeline()
    
    elif args.command == 'vcf_generator':
        line = Lines(core_log, args.name)
        
        line.usr_in_line_variables(args.reference_genome_name, 
            args.wt_input, args.mu_input
        )
        vcf_gen_vars = VCFGenVariables(core_log, 
            lines=line,
            reference_genome_source=args.reference_genome_source, 
            snpEff_species_db=args.snpEff_species_db,
            reference_genome_name=args.reference_genome_source, 
            known_snps=args.known_snps, 
            threads_limit=args.threads_limit, 
            call_variants_in_parallel=args.call_variants_in_parallel,
            cleanup=args.cleanup
        )

        vcf_gen = VCFGenerator(core_log, vcf_gen_vars)
        vcf_gen.run_subprocess()

    elif args.automatic:
        auto_vars = AutomaticLineVariableDetector(core_log)
        auto_vars.automatic_line_variables()

        vcf_vars = VCFGenVariables(core_log, lines=auto_vars.lines)
        vcf_gen = VCFGenerator(core_log, vcf_vars)
        vcf_gen.run_subprocess()
        
        bsa_vars = BSAVariables(core_log, lines=vcf_gen.lines)
        bsa = BSA(core_log, bsa_vars)
        bsa.run_pipeline()

    elif args.command == 'logdb':
        logdb_utils = LogDbUtilites()
        if args.vcf_ulid_log:
            logdb_utils.print_vcf_log_data(args.vcf_ulid_log)
        if args.analysis_ulid_log:
            logdb_utils.print_analysis_log_data(args.analysis_ulid_log)
        if args.line_name_log:
            logdb_utils.print_line_name_data(args.line_name_log)
        if args.core_ulid_log:
            logdb_utils.print_core_ulid_data(args.core_ulid_log)

if __name__ == "__main__":
    main()
