#!/usr/bin/env python
from modules.arg_parser import ArgumentParser
arg_parser = ArgumentParser()
args = arg_parser.args

from modules.utilities_variables import (
    Lines,AutomaticLineVariableDetector, VCFGenVariables, BSAVariables
)
from modules.utilities_logging import LogHandler
from modules.core_vcf_gen import VCFGenerator
from modules.core_bsa import BSA
from modules.utilities_general import LogDbUtilites


def main():
    # Primary program

    # initialize core log, unless trying to pull log information
    '''
    All classes in this program have LogHandler instance passed to them on init.

    Unique ID (ulid) are assigned upon init of LogHandler. 
    ulids link file outputs to their log. 

    ulid's are stored in LogHandler instances, and exported to the Lines data
    class.  

    LogHandler manages the log database to save run information. 

    log list:
    'core' - logs all parent functions, thale_bsa.py and flask front end. 
    'vcf_gen' - logs all messages pertaining to parent_functions.vcf_generation
    'analysis' - logs all messages peratining to parent_functions.bsa_analysis 
    '''
    
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
            reference_genome_name=args.reference_genome_name, 
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
