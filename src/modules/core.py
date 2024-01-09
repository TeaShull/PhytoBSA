import os
import subprocess
import fnmatch
from datetime import datetime
import pandas as pd
import re
from flask import session  # Keeping the imports for later use
from config import (
    LogHandler, INPUT_DIR, OUTPUT_DIR, MODULES_DIR, LOG_DIR, REFERENCE_DIR
)

from utilities_bsa_analysis import BSAAnalysisUtilities
from utilities_file import FileUtilities

class ThaleBSAParentFunctions:
    def __init__(self, logger):
        self.log = logger 

    def vcf_generation(self, experiment_dictionary,
                   reference_genome_name, snpEff_species_db,
                   reference_genome_source, threads_limit,
                   cleanup, known_snps)->dict:
        """Input: Experiment dictionary as well as paths and variables needed 
        to run VCFgen.sh. Subprocess VCFgen.sh takes raw reads(wild-type(wt) 
        and mutant(mu)), either single or paired end and generates VCF table 
        *.noknownsnps.table.
        
        Returns: updated experiment_dictionary containing the paths to 
        the noknownsnps.table."""
        core_ulid=experiment_dictionary['core_ulid']
        
        for key, value in experiment_dictionary.items():
            self.log.attempt(f"Generating VCF file for {key}...")
            self.log.delimiter(f"Shell [sh] VCF generator for {key} beginning...")
        
            try:
                current_line_name = key

                # generate log instance, add run info to sql db
                vcf_log = LogHandler(f'vcf_{current_line_name}')
                self.log.note(f'Logging for VCF Initialzed. Path: {vcf_log.log_path}')

                vcf_log.note(f'vcf_log initiated.')
                vcf_log.note(f'vcf_log ulid: {vcf_log.ulid}')
                experiment_dictionary[current_line_name]['vcf_ulid'] = vcf_log.ulid
                vcf_log.add_db_record(current_line_name, core_ulid)

                # Add output_path to experiment_dictionary. 
                output_name_prefix = f"{vcf_log.ulid}_-{current_line_name}"
                output_dir_path = os.path.join(OUTPUT_DIR, output_name_prefix)
                output_prefix = os.path.join(output_dir_path, output_name_prefix)
                experiment_dictionary[key]['output_dir_path'] = output_dir_path 

                # Add vcftable_path to experiment_dictionary.
                vcf_table_name = f"{output_name_prefix}.noknownsnps.table"
                vcf_table_path = os.path.join(
                    OUTPUT_DIR, current_line_name, vcf_table_name)
                experiment_dictionary[current_line_name]['vcf_table_path'] = vcf_table_path
                
                #Generate VCFgen.sh script path
                vcfgen_script_path = os.path.join(MODULES_DIR, 'VCFgen.sh')
                
                #Generate the knownSnps .vcf file path
                known_snps = os.path.join(REFERENCE_DIR, known_snps)
                # Retrieve allele and file input info from experiment_dictionary
                allele = value['allele']
                pairedness = value['pairedness']

                wt_input = ' '.join(value['wt'])
                wt_input = f'"{wt_input}"'
                self.log.note(f"wt_input:{wt_input}")

                mu_input = ' '.join(value['mu'])
                mu_input = f'"{mu_input}"'
                self.log.note(f"mu_input:{mu_input}")

                # Construct args to pass variable and experiment_dictionary to
                # VCFgen.sh.
                args = (
                    vcf_log.ulid,
                    current_line_name, 
                    allele,
                    INPUT_DIR,
                    wt_input,
                    mu_input,
                    pairedness,
                    output_dir_path,
                    output_prefix,
                    vcf_table_path,
                    REFERENCE_DIR,
                    reference_genome_name, 
                    snpEff_species_db,
                    reference_genome_source, 
                    known_snps,
                    threads_limit,
                    cleanup
                )

                # Construct the command for VCFgen.sh, passing the above variables
                cmd = f"{vcfgen_script_path} {' '.join(map(str, args))}"

                # Run vcfgen shell subprocess.
                process = subprocess.Popen(
                    cmd, cwd=MODULES_DIR, shell=True, stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, text=True
                )

                # Iterate over stdout from process andlog
                for line in process.stdout:
                    vcf_log.bash(line.strip())
                process.wait()
                
                self.log.note(f"VCF file generated for {current_line_name}.") 
                self.log.note("Log saved to {log_path}")
                self.log.note(f'VCF table path added to experiments_dictionary: {vcf_table_path}')
            
            except Exception as e:
                self.log.fail(f"Error while generating the VCF file for {current_line_name}: {e}")

        self.log.success("VCF file generation process complete")
        
        return experiment_dictionary

    def bsa_analysis(self, experiment_dictionary):
        '''
        Parent function for running BSA analysis. Given an experiment_dictionary
        containing the vcf_table_path and line_name, will output a list of 
        candidate mutations and 5 graphs, which will help narrow down regions of
        interest.  
        
        Dictionary Requirements:
        required 
        [key]: line_name
        [value] vcf_table_path, core_ulid, 

        optional
        [value] vcf_ulid
        
        vcf_table_path must be the path to a vcf table produced from VCFgen.sh. 
        support for VCF tables which are configured differently is on the roadmap. 
        
        '''
        self.log.attempt("Attempting to perform data analysis...")
        try:
            for key, value in experiment_dictionary.items():
                core_ulid = value['core_ulid']
                current_line_name = key
                vcf_ulid = value['vcf_ulid']
                vcf_table_path = value['vcf_table_path']
                # Configure an analysis logger for each line. 
                analysis_log = LogHandler(f'analysis_{current_line_name}')
                analysis_log.note('Analysis log initiated.') 
                analysis_log.note(f'VCF ulid:{vcf_ulid}')
                analysis_log.add_db_record(current_line_name, core_ulid, vcf_ulid)
                
                #FileUtilites instance that logs to analysis_log
                file_utils = FileUtilities(analysis_log)

                #Analysis operations. Loading VCF and feature production
                vcf_df = file_utils.load_vcf_table(vcf_table_path, current_line_name)
                
                bsa_analysis_utils = BSAAnalysisUtilities(current_line_name, vcf_ulid, analysis_log)
                vcf_df = bsa_analysis_utils.calculate_delta_snp_and_g_statistic(vcf_df)
                vcf_df = bsa_analysis_utils.drop_na_and_indels(vcf_df) 
                vcf_df = bsa_analysis_utils.loess_smoothing(vcf_df)
                vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff = (
                    bsa_analysis_utils.calculate_empirical_cutoffs(vcf_df)
                )
    
                #Saving and plotting outputs
                bsa_analysis_utils.sort_save_likely_candidates(vcf_df)
                bsa_analysis_utils.generate_plots(vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff)

            self.log.success("Data analysis complete")
        
        except Exception as e:
            self.log.fail(f"Error during data analysis: {e}")


