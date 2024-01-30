#!/usr/bin/env python
import os
import sys
import argparse
import configparser

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
    parser_analysis.add_argument('-vt', '--vcf_table_path', required=True, type=str, help='path to the vcf table you wish to analyze.')
    parser_analysis.add_argument('-st', '--segregation_type', required=True, type=str, help="Recessive(R) or Dominant(D)")
    ### OPTIONAL
    parser_analysis.add_argument('-ls', '--loess_span', type=float, default=None, help="Influences smoothing parameters.")
    parser_analysis.add_argument('-si', '--shuffle_iterations', type=int, default=None, help="Iterations of bootstrapping during empirical cutoff calculations. Below 1000 can yeild inconsistant results")
    parser_analysis.add_argument('-sb', '--smooth_edges_bounds', type=int, default=None, help="Number of mirrored datapoints at chromosome edges to correct for loess edge bias. Increase if edge bias seems high")
    parser_analysis.add_argument('-fin', '--filter_indels', default=None, type=str, help="Filter out insertion-deletion mutations.")
    parser_analysis.add_argument('-fems', '--filter_ems', default=None, type=str, help="Filter results to only include mutations likely to arise from EMS treatment")
    parser_analysis.add_argument('-snpmsk', '--snpmask_path', default=None, type=str, help="Path to VCF file containing background snps.")

    ## VCF generation subparser
    parser_vcf_gen = subparsers.add_parser('vcf_generator', help='Run PhytoBSA VCF Generator seperately')
    ### REQUIRED
    parser_vcf_gen.add_argument('-wt', '--wt_input', required=True, default=None, type=str, help='Wild-type bulk fasta file(s)')
    parser_vcf_gen.add_argument('-mu', '--mu_input', required=True, default=None, type=str, help='mutant bulk fast files')
    ### OPTIONAL - the below arguments can be input via command line or sourced from ./settings/config.ini
    parser_vcf_gen.add_argument('-rgn', '--reference_genome_name', required=False, default=None, type=str, help='Reference genome name')
    parser_vcf_gen.add_argument('-ssdb', '--snpEff_species_db', default=None, type=str, help = 'The name of your snpEff database name.')
    parser_vcf_gen.add_argument('-rgs', '--reference_genome_source', default=None, type=str, help = 'Optional, if you wish the pipeline to download your reference from a url')
    parser_vcf_gen.add_argument('-ks','--known_snps', default=None, type=str, help = 'VCF file containing background SNPs. Helps improve output quality')
    parser_vcf_gen.add_argument('-t', '--threads_limit', default=None, type=str, help='Maximum threads you wish to use for analysis')
    parser_vcf_gen.add_argument('-p', '--call_variants_in_parallel', default=None, type=bool, help='Run gatk haplotype caller in parallel')
    parser_vcf_gen.add_argument('-c','--cleanup', default=None, type=bool, help='If true, intermediate files will de deleted. False for troubleshooting and archiving files.' )
    
    ## Log database subparser
    parser_log_db = subparsers.add_parser('logdb', help='Gather log information from log database')
    parser_log_db.add_argument('-vcf', "--vcf_ulid_log", default = None, type = str, help="Retrieve vcf log information based on vcf ulid input")
    parser_log_db.add_argument('-an', "--analysis_ulid_log", default = None, type = str, help="Retrieve analysis log information based on analysis ulid input")
    parser_log_db.add_argument('-name', "--line_name_log", default = None, type = str, help="Retrieve all information associated with the provided line name. ")    
    parser_log_db.add_argument('-core', '--core_ulid_log', default = None, type = str, help="Retrieve all information associated with the provided core log ulid")
    
    ## Settings subparser
    parser_settings = subparsers.add_parser('settings', help='Update default settings.')
    parser_settings.add_argument('--set_data_dir', default=None, type=str, help='set Data directory. This must be set for program to run')
    #vcf_gen default run settings. 
    parser_settings.add_argument('--set_reference_genome_name', required=False, default=None, type=str, help='Set default reference genome name')
    parser_settings.add_argument('--set_snpEff_species_db', default=None, type=str, help = 'Set default snpEff database name.')
    parser_settings.add_argument('--set_reference_genome_source', default=None, type=str, help = 'Set default reference genome source')
    parser_settings.add_argument('--set_known_snps', default=None, type=str, help = 'Set default VCF file containing background SNPs.')
    parser_settings.add_argument('--set_threads_limit', default=None, type=str, help='Set default maximum threads for analysis')
    parser_settings.add_argument('--set_call_variants_in_parallel', default=None, type=bool, help='Set default for running gatk haplotype caller in parallel')
    parser_settings.add_argument('--set_cleanup', default=None, type=bool, help='Set default for cleanup. If true, intermediate files will de deleted. False for troubleshooting and archiving files.' )
    #BSA default run settings.
    parser_settings.add_argument('--set_loess_span', type=float, help="Set default loess_span.")
    parser_settings.add_argument('--set_shuffle_iterations', type=int, help="Set default shuffle_iterations.")
    parser_settings.add_argument('--set_smooth_edges_bounds', type=int, help="Set default smooth_edges_bounds.")
    parser_settings.add_argument('--set_filter_indels', type=bool, help="Set default filter_indels.")
    parser_settings.add_argument('--set_filter_ems', type=bool, help="Set default filter_ems.")
    parser_settings.add_argument('--set_snpmask_path', type=str, help="Set default snp mask VCF path. VCF should contain known background SNPs")
    
    #parse args
    args = parser.parse_args()
    
    # If a parameter is not passed that exists in settings - it will be sourced. 
    # If settings parser is used (commands starting with set_,) config.ini will
    # be updated to reflect the new changes. 
    config = configparser.ConfigParser()
    
    config_ini = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'settings', 'config.ini')
    config.read(config_ini)
    
    if args.command == 'settings':
        for arg in vars(args):
            if arg.startswith('set_') and getattr(args, arg) is not None:
                    value = getattr(args, arg)
                    config_key = arg.replace('set_', '')
                    print(f'Default {config_key} is now set to {value}')
                    config.set('Settings', config_key, str(value))
        with open('config.ini', 'w') as configfile:
            config.write(configfile)
        quit()
    
    for arg in vars(args):
        if getattr(args, arg) is None and config.has_option('Settings', arg):
            setattr(args, arg, config.get('Settings', arg))
    return args
def main():
    
    # Parse command line arguments
    args = parse_program_arguments()
    if args.command == 'settings':
        print('Settings applied!')

    
    # Primary program
    from modules.utilities_variables import (
        Lines,AutomaticLineVariableDetector, VCFGenVariables, BSAVariables
    )
    from modules.utilities_logging import LogHandler
    from modules.core_vcf_gen import VCFGenerator
    from modules.core_bsa import BSA
    from modules.utilities_general import LogDbUtilites

    # initialize core log, unless trying to pull log information
    if args.command != 'logdb':
        core_log = LogHandler('core')
        core_log.note(f'Core log begin. ulid: {core_log.ulid}')
        core_log.add_db_record()
    
    # Determine which routine to run
    if args.command == 'analysis':
        line = Lines(core_log, args.name)
    
        line.usr_in_line_variables(args.vcf_table_path, arg.segregation_type)
        
        bsa_vars = BSAVariables(core_log, 
            lines=line, 
            loess_span=args.loess_span, 
            smooth_edges_bounds=args.smooth_edges_bounds, 
            shuffle_iterations=args.shuffle_iterations,
            filter_indels=args.filter_indels,
            filter_ems=args.filter_ems,
            snpmask_path=args.snpmask_path
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

        vcf_vars = VCFGenVariables(core_log, 
            lines=auto_vars.lines,
            reference_genome_source=args.reference_genome_source, 
            snpEff_species_db=args.snpEff_species_db,
            reference_genome_name=args.reference_genome_source, 
            known_snps=args.known_snps, 
            threads_limit=args.threads_limit, 
            call_variants_in_parallel=args.call_variants_in_parallel,
            cleanup=args.cleanup
        )
        vcf_gen = VCFGenerator(core_log, vcf_vars)
        vcf_gen.run_subprocess()
        
        bsa_vars = BSAVariables(core_log, 
            lines=vcf_gen.lines, 
            loess_span=args.loess_span, 
            smooth_edges_bounds=args.smooth_edges_bounds, 
            shuffle_iterations=args.shuffle_iterations,
            filter_indels=args.filter_indels,
            filter_ems=args.filter_ems,
            snpmask_path=args.snpmask_path
        )
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
