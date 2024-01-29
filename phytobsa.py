#!/usr/bin/env python
import os
import sys
import argparse

from modules.utilities_variables import (
    Lines, AutomaticVariables, VCFGenVariables, BSAVariables
)
from modules.utilities_logging import LogHandler
from modules.core_vcf_gen import VCFGenerator
from modules.core_bsa import BSA

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
    args = parser.parse_args()
    return args

def main():
    # initialize core log
    core_log = LogHandler('core')
    core_log.note(f'Core log begin. ulid: {core_log.ulid}')
    core_log.add_db_record()
    
    # Parse command line arguments
    args = parse_program_arguments()

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

if __name__ == "__main__":
    main()
