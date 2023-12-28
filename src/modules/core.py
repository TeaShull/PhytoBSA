import os
import subprocess
import fnmatch
import logging
from datetime import datetime
import pandas as pd
from flask import session  # Keeping the imports for later use
from config import setup_logger, log_handler, INPUT_DIR, OUTPUT_DIR, MODULES_DIR, LOG_DIR
from analysis import AnalysisUtilities, PlotUtilities

class ThaleBSAParentFunctions:
    def __init__(self, logger):
        self.logger = logger 

    def vcf_generation(self, experiment_dictionary,
                       reference_genome_name, snpEff_db_name,
                       reference_genome_source, threads_limit,
                       cleanup, known_snps, vcf_table_uuid):

        log_handler(self.logger, 'info', 'Generating VCF files for experiments in dictionary')
        for key, value in experiment_dictionary.items():
            log_handler(self.logger, 'attempt', f"Generating VCF file for {key}...")
            try:
                # Construct cmd
                modules_dir = MODULES_DIR
                log_dir = LOG_DIR
                vcfgen_script_path = os.path.join(modules_dir, 'VCFgen.sh')
                args = (key, value['reads'], value['allele'],
                        reference_genome_name, snpEff_db_name,
                        reference_genome_source, threads_limit,
                        cleanup, known_snps, vcf_table_uuid
                        )

                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                log_name = f"{key}.VCF_generation.{timestamp}.{vcf_table_uuid}.log"
                log_path = os.path.join(log_dir, log_name)
                cmd = f"{vcfgen_script_path} {' '.join(map(str, args))} | tee {log_path}"

                subprocess.run(cmd, shell=True, text=True, check=True)

                log_handler(self.logger, 'success',
                              f"VCF file generated for {key}. Log saved to {log_path}"
                              )
            except Exception as e:
                log_handler('error',
                              f"Error while generating the VCF file for {key}: {e}"
                              )

        log_handler(self.logger, 'info', "VCF file generation process complete")

    def bsa_analysis(self, experiment_dictionary):
        
        log_handler(self.logger, 'info', "Attempting to perform data analysis...")
        try:
            for key, value in experiment_dictionary.items():
                current_line_name = key
                #Establishing VCF dataframe and variables
                # Configure the analysis logger for each line. 
                log_dir = LOG_DIR
                output_dir = OUTPUT_DIR
                timestamp = datetime.now().strftime("%Y.%m.%d_%H:%M")
                
                log_filename = f"{timestamp}_{current_line_name}_analysis.log"
                log_path = os.path.join(log_dir, log_filename)
                analysis_logger = setup_logger('analysis_logger', log_path)    
                
                file_utils = FileUtilities(analysis_logger)

                reads = value['reads']
                if value['allele'] == 'R':
                    allele = 'recessive'
                elif value['allele'] == 'D':
                    allele = 'dominant'
                print(f"analysis logger name? {analysis_logger.name}")
                log_handler(analysis_logger, 'note', f"{current_line_name} is labeled {allele}. Pairdness is {reads}")
                vcftable_name = f"{current_line_name}.noknownsnps.table"
                current_line_table_path = os.path.join(output_dir, current_line_name, vcftable_name)
                vcf_df = file_utils.load_vcf_table(current_line_table_path, current_line_name)

                #Analysis operations. sequential dataframe transformation
                analysis_utils = AnalysisUtilities(current_line_name, analysis_logger)
                vcf_df = analysis_utils.calculate_delta_snp_and_g_statistic(vcf_df)
                vcf_df = analysis_utils.drop_na_and_indels(vcf_df)
                vcf_df = analysis_utils.loess_smoothing(vcf_df)
                vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff = analysis_utils.calculate_empirical_cutoffs(vcf_df)
                
                #Saving and plotting outputs
                analysis_utils.sort_save_likely_candidates(vcf_df)
                plot_utils = PlotUtilities(current_line_name, analysis_logger)
                plot_utils.generate_plots(vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff)

            log_handler(self.logger, 'info', "Data analysis complete")
        except Exception as e:
            log_handler(self.logger, 'fail', f"Error during data analysis: {e}")

class FileUtilities:
    def __init__(self, logger):
        self.logger = logger 

    def detect_file_type(self, file):
        log_handler(self.logger, 'attempt', 
            f"Detecting if {file} is labeled recessive (R) or dominant (D)..."
        )
        try:
            if fnmatch.fnmatch(file, '*.R*'):
                label = "R"
            elif fnmatch.fnmatch(file, '*.D*'):
                label = "D"
            log_handler(self.logger, 'success', f"{file} is labeled {label}")
            return label
        except Exception as e:
            log_handler(self.logger, 'fail', f"Error while detecting file type: {e}")
            return None

    def create_experiment_dictionary(self):
        log_handler(self.logger, 'attempt', 
            "Creating a dictionary to store experiment details..."
            )
        try:
            input_dir = INPUT_DIR
            lines_dict = {}

            for file in os.listdir(input_dir):
                log_handler(self.logger, 'attempt', f"Parsing {file}...")
                try:
                    key = file.split(".")[0]
                    lines_dict[key] = lines_dict.get(
                        key, {'count': 0, 'allele': None}
                    )
                    lines_dict[key]['count'] += 1
                    lines_dict[key]['allele'] = self.detect_file_type(file)
                    log_handler(self.logger, 'success', 
                        f"{file} parsed and added to dictionary under key {key}"
                    )
                except Exception as e:
                    log_handler(self.logger, 'fail', 
                        f"Error parsing {file} for the experiment dictionary: {e}"
                    )

            # Convert counts and types to a readable format
            for key, value in lines_dict.items():
                log_handler(self.logger, 'attempt', "Making dictionary more readable")
                try:
                    if value['count'] == 4:
                        value['reads'] = "paired-end"
                    elif value['count'] == 2:
                        value['reads'] = "single-read"
                    log_handler(self.logger, 'success', "Dictionary is now more readable")
                except Exception as e:
                    log_handler(self.logger, 'fail', 
                        f"Error converting formatting counts and types: {e}"
                    )

            log_handler(self.logger, 'success', "Experiment Dictionary Created ")
            return lines_dict

        except Exception as e:
            log_handler(self.logger, 'fail', 
                f"Error while creating the experiment dictionary: {e}"
            )
            return {}

    def load_vcf_table(self, current_line_table_path, current_line_name):
        """Load VCF table"""
        log_handler(self.logger, 'attempt', 
            f"Attempting to load VCF table for line {current_line_name}"
        )
        try:
            vcf_df = pd.read_csv(current_line_table_path, sep="\t")
            log_handler(self.logger, 'attempt', 
                f"The VCF table for line {current_line_name} was successfully loaded."
            )
            return vcf_df
        except FileNotFoundError:
            log_handler(self.logger, 'fail', 
                f"Error: File '{current_line_table_path}' not found."
            )
        except pd.errors.EmptyDataError:
            log_handler(self.logger, 'fail', 
                f"Error: File '{current_line_table_path}' is empty."
            )
        except Exception as e:
            log_handler(self.logger, 'fail', f"An unexpected error occurred: {e}")

class GeneralUtilities:
    """General Utilities"""
    def __init__(self, logger):
        self.logger = logger 

    def int32_to_id(n):
        """returns a human readable id for each integer fed. 
            the id is reproducable and can create half a billion ids"""
        if n==0: return "0"
        chars="0123456789ACEFHJKLMNPRTUVWXY"
        length=len(chars)
        result=""
        remain=n
        while remain>0:
            pos = remain % length
            remain = remain // length
            result = chars[pos] + result
        return result
  

    def extract_vcf_id(file_name):
        # Get the base name of the file without the extension
        base_name, extension = os.path.splitext(file_name)

        # Split the base name by hyphen
        parts = base_name.split('-')

        # Check if there is at least one hyphen-separated part
        if len(parts) > 1:
            # Assume the vcf_id is the last part
            vcf_id = parts[-1]
            return vcf_id
        else:
            # If no hyphen-separated parts, or only one part, return None
            return None

class ThaleBSASQLDB:
    """Handling the log database, so that analyses can be associated with 
    their respective VCF file runs. IN PROGRESS...."""

    def __init__(self, logger, db_name="thale_bsa_sqldb.db"):
        self.conn = sqlite3.connect(db_name)
        self.logger = logger 
        self.create_tables()

        # Initialize the last used IDs
        self.last_vcf_id = self.get_last_id("vcf_id")
        self.last_analysis_id = self.get_last_id("analysis_id")

    def create_tables(self):
        """Create the VCF table with additional fields"""
        self.conn.execute('''
            CREATE TABLE IF NOT EXISTS thale_bsa_sqldb (
                analysis_id TEXT PRIMARY KEY,
                line_name TEXT,
                vcf_id INTEGER,
                vcf_log_path TEXT,
                analysis_log_path TEXT,
                run_date TEXT,
                vcf_file_id TEXT
            )
        ''')
        self.conn.commit()

    def get_last_id(self, id_type):
        """Retrieve the last used ID from the database"""
        cursor = self.conn.execute(f'''
            SELECT MAX({id_type}) FROM thale_bsa_sqldb
        ''')
        result = cursor.fetchone()
        return result[0] if result and result[0] is not None else 0

    def generate_vcf_id(self):
        """Increment and return the next vcf_id"""
        self.last_vcf_id += 1
        return self.last_vcf_id

    def generate_analysis_id(self):
        """Increment and return the next analysis_id"""
        self.last_analysis_id += 1
        return self.last_analysis_id

    def add_record(self, line_name, vcf_log_path, analysis_log_path, run_date):
        """Generate new IDs"""
        vcf_id = self.generate_vcf_id()
        analysis_id = self.generate_analysis_id()

        """Add a new record to the database"""
        self.conn.execute('''
            INSERT INTO thale_bsa_sqldb (analysis_id, line_name, vcf_id, vcf_log_path, analysis_log_path, run_date, vcf_file_id)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (analysis_id, line_name, vcf_id, vcf_log_path, analysis_log_path, run_date, f"{analysis_id}_{vcf_id}"))
        self.conn.commit()

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



