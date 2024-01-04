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

from analysis import AnalysisUtilities

class ThaleBSAParentFunctions:
    def __init__(self, logger):
        self.log = logger 

    def vcf_generation(self, experiment_dictionary,
                   reference_genome_name, snpEff_species_db,
                   reference_genome_source, threads_limit,
                   cleanup, known_snps):
        """Input raw reads, either single or paired end. 
        Returns: the VCF file is generated, and the function returns an updated
        experiment_dictionary now containing the paths to the noknownsnps.table
        files that were just generated. """
        for key, value in experiment_dictionary.items():
            self.log.attempt(f"Generating VCF file for {key}...")
            self.log.delimiter(f"Shell [sh] VCF generator for {key} beginning...")
        
            try:
                current_line_name = key

                # generate log instance
                vcf_log = LogHandler(f'vcf_{current_line_name}')
                self.log.note(f'Logging for VCF Initialzed. Path: {vcf_log.log_path}')
                
                experiment_dictionary[current_line_name]['vcf_ulid'] = vcf_log.ulid

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
                
                # Retrieve allele and file input info from experiment_dictionary
                allele = value['allele']

                wt_input = ' '.join(value['wt'])
                wt_input = f'"{wt_input}"'
                self.log.note(f"wt_input:{wt_input}")
                mu_input = ' '.join(value['mu'])
                mu_input = f'"{mu_input}"'

                # Construct args to pass variable and experiment_dictionary to
                # VCFgen.sh.
                args = (
                    vcf_log.ulid,
                    current_line_name, 
                    allele,
                    INPUT_DIR,
                    wt_input,
                    mu_input,
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
                self.log.note(f'VCF table path added to experiments_dictionary: {vcftable_path}')
            
            except Exception as e:
                self.log.fail(f"Error while generating the VCF file for {current_line_name}: {e}")

        self.log.success("VCF file generation process complete")
        
        return experiment_dictionary

    def bsa_analysis(self, experiment_dictionary):
        self.log.attempt("Attempting to perform data analysis...")
        
        try:
            for key, value in experiment_dictionary.items():
                current_line_name = key
                vcf_ulid = value['vcf_ulid']
                # Configure an analysis logger for each line. 
                analysis_log = LogHandler(f'analysis_{current_line_name}')

                #FileUtilites instance that logs to analysis_log
                file_utils = FileUtilities(analysis_log)

                #Analysis operations. Loading VCF and feature production
                vcf_df = file_utils.load_vcf_table(vcftable_path, current_line_name)
                
                analysis_utils = AnalysisUtilities(current_line_name, vcf_ulid, analysis_log)
                vcf_df = analysis_utils.calculate_delta_snp_and_g_statistic(vcf_df)
                vcf_df = analysis_utils.drop_na_and_indels(vcf_df) 
                vcf_df = analysis_utils.loess_smoothing(vcf_df)
                vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff = (
                    analysis_utils.calculate_empirical_cutoffs(vcf_df)
                )
                
                #Saving and plotting outputs
                analysis_utils.sort_save_likely_candidates(vcf_df)
                analysis_utils.generate_plots(vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff)

            self.log.success("Data analysis complete")
        
        except Exception as e:
            self.log.fail(f"Error during data analysis: {e}")

class FileUtilities:
    def __init__(self, logger):
        self.log = logger 

    def experiment_detector(self):
        self.log.attempt(f"Detecting experiment details in: {INPUT_DIR}...")
        try:
            expt_dict = {}

            # Iterate through the files in the directory
            for filename in os.listdir(INPUT_DIR):

                # Split the filename into parts based on dots
                parts = filename.split('.')
                self.log.attempt(f'Parsing {filename}')

                # Extract relevant information
                line_name = parts[0]
                allele = parts[1]
                
                if '_1' in allele:
                    allele = allele.rstrip('_1')
                if '_2' in allele:
                    allele = allele.rstrip('_2')

                segregation_type = parts[-3]
                pairedness = 'paired-end' if '_1' or '_2' in filename else 'single-read'

                # Use line_name as key for the dictionary
                key = line_name
                self.log.success(f"""{filename} parsed. 
                    key:{key}
                    allele:{allele}
                    segregation_type:{segregation_type}
                    pairedness:{pairedness}
                """)

                # Initialize or update the dictionary entry for the key
                if key not in expt_dict:
                    expt_dict[key] = {
                        'allele': allele,
                        'wt': [],
                        'mu': [],
                        'pairedness': pairedness
                    }

                # Add the file path to the appropriate list based on segregation_type
                file_path = os.path.join(INPUT_DIR, filename)

                if 'wt' in segregation_type:
                    wt_list = expt_dict[key]['wt']
                    wt_list.append(file_path)

                    # Sort the wt_list to ensure _1 and _2 files are in numeric order
                    expt_dict[key]['wt'] = sorted(wt_list, key=lambda x: int(x.split('_')[-1][0]))

                elif 'mu' in segregation_type:
                    mu_list = expt_dict[key]['mu']
                    mu_list.append(file_path)

                    # Sort the mu_list to ensure _1 and _2 files are in numeric order
                    expt_dict[key]['mu'] = sorted(mu_list, key=lambda x: int(x.split('_')[-1][0]))
            
            self.log.success(f'Experiment dictionary generated.')
            
            return expt_dict

        except Exception as e:
            self.log.fail(f"Error while detecting experiment details: {e}")
            
            return {}

    def load_vcf_table(self, current_line_table_path, current_line_name):
        """Load VCF table"""
        self.log.attempt(f"Attempting to load VCF table for line {current_line_name}")
        
        try:
            vcf_df = pd.read_csv(current_line_table_path, sep="\t")
            self.log.attempt(f"The VCF table for line {current_line_name} was successfully loaded.")
            return vcf_df
        
        except FileNotFoundError:
            self.log.fail(f"Error: File '{current_line_table_path}' not found.")
        
        except pd.errors.EmptyDataError:
            self.log.fail(f"Error: File '{current_line_table_path}' is empty.")
        
        except Exception as e:
            self.log.fail(f"An unexpected error occurred: {e}")

    def setup_output_directory(experiment_key):
        # Create the output directory path
        current_output_dir = os.path.join(OUT_DIR, experiment_key)
        self.log.attempt('Checking if output directory exists...')
        
        try: 
            # Check if the output directory exists, and create it if necessary
            if not os.path.exists(current_output_dir):
                self.log.attempt(f"Output directory does not exist. Creating: {current_output_dir}")
           
                os.makedirs(current_output_dir)
                self.log.success(f'Output directory created: {current_output_dir}')
            else:
                self.log.note(f"Output directory already exists: {current_output_dir}")
        
            return current_output_dir
        except Exception as e:
            self.log.fail('setting up output directory failed: {e}')

    def extract_ulid_from_filename(filename):
        ulid_pattern = re.compile(r'^[0-9A-HJKMNPQRSTVWXYZ]{26}')
        match = ulid_pattern.match(filename)
        if match:
            return match.group()
        else:
            return None

class ThaleBSASQLDB:
    """Handling retrieving from log database.
     IN PROGRESS...."""

    def __init__(self, logger, db_name="thale_bsa_sqldb.db"):
        self.conn = sqlite3.connect(db_name)
        self.log = logger 

    def get_vcf_data(self, analysis_id):
        """Retrieve the VCF data based on the analysis ID"""
        cursor = self.conn.execute('''
            SELECT line_name, vcf_id, vcf_log_path, analysis_log_path, run_date FROM thale_bsa_sqldb WHERE analysis_id = ?
        ''', (analysis_id,))
        result = cursor.fetchone()
        return result if result else None

    def close_connection(self):
        """Close the database connection"""
        self.conn.close()

    def print_paths_for_analysis_id(self, analysis_id):
        """Retrieve the paths based on the analysis ID"""
        cursor = self.conn.execute('''
            SELECT analysis_log_path, vcf_log_path FROM thale_bsa_sqldb WHERE analysis_id = ?
        ''', (analysis_id,))
        result = cursor.fetchone()

        if result:
            analysis_log_path, vcf_log_path = result
            print(f"Analysis ID: {analysis_id}")
            print(f"Analysis Log Path: {analysis_log_path}")
            print(f"VCF Log Path: {vcf_log_path}")
        else:
            print(f"No record found for Analysis ID '{analysis_id}'")



