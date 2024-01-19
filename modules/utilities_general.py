from settings.config import INPUT_DIR, OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR

import os
import pandas as pd
import re
import sqlite3

class FileUtilities:
    """
    Utility class for handling file operations and inputs. 
    
    """
    def __init__(self, logger):
        self.log = logger 

    def experiment_detector(self)->dict:
        '''
        Detects potential BSA experiments in the inputs folder if files are 
        named according to the conventions outlined in the README
        
        For paired-end:
        <line_name>.<R or D>_<read number>.<wt or mu>.fq.gz  
        
        For single-read:
        <line_name>.<R or D>.<wt or mu>.fq.gz

        Args: None

        Returns: experiment_dictionary containing the detected information. If
        successful, all information needed to generate vcf files and run analysis
        will have been parsed. 
        '''
        self.log.attempt(f"Detecting experiment details in: {INPUT_DIR}")
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
                    expt_dict[key]['wt'] = sorted(
                        wt_list, key=lambda x: int(x.split('_')[-1][0])
                    )
                elif 'mu' in segregation_type:
                    mu_list = expt_dict[key]['mu']
                    mu_list.append(file_path)
                    # Sort the mu_list to ensure _1 and _2 files are in numeric order
                    expt_dict[key]['mu'] = sorted(
                        mu_list, key=lambda x: int(x.split('_')[-1][0])
                    )
            
            self.log.success(f'Experiment dictionary generated.')
            return expt_dict

        except Exception as e:
            self.log.fail(f"Error while detecting experiment details: {e}")
            return {}

    def check_vcf_gen_variables(self, experiment_dict: dict)-> dict:
        """
        Checks if the user has provided all the necessary variables
        experiment_dictionary will be populated with the contents of 
        vcf_gen_variables, if the user has not assigned any variables. 

        Please update this list if you add any new variables to vcf_gen_variables
        Args:
            experiment_dict:
                key:
                    reference_genome_name
                    snpEff_species_db
                    reference_genome_source
                    threads_limit
                    cleanup
                    known_snps
                    call_variants_in_parallel
        Returns:
        Variables sourced from settings.vcf_gen_variables module if they are missing
        True if variables are all accounted for
        """
        try:
            import settings.vcf_gen_variables as var
            for line, details in experiment_dict.items():
                for key, value in details.items():
                    if value:
                        self.log.note(f'Variable set by user|{key}:{value}')
                    if value is None:
                        try:
                            self.log.note(f'{key} variable not set. Sourcing from settings/vcf_gen_variables...')
                            sub_dict[key] = var.__dict__[key]
                            self.log.note(f'Variable assigned|{key}:{sub_dict[key]}')
                        except Exception as e:
                            self.log.fail(f'Aborting. Setting variable {key} failed: {e}')
            return experiment_dict
        except Exception as e:
            self.log.fail(f'There was an error while checking if variables for VCFgen.sh have been assigned: {e}')

    def load_vcf_table(self, current_line_table_path, current_line_name)->pd.DataFrame:
        """
        Loads VCF table into a pandas dataframe.
        
        Args:  
        current_line_table_path(str) - path to the vcf table to be loaded into df
        current_line_name(str) - name of the line associated with the vcf table

        Returns: 
        Pandas dataframe containing the information loaded from current_line_table_path
        """

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

    def save_experiment_details(self, experiment_dictionary):
        '''
        Saves experiment_dictionary information into a human-readable file, for
        easy veiwing. Allows easy access to information without searching through 
        logs. 

        Args: 
        experiment_dictionary(dict) - experiment dictionary that contains experiment info. 

        Returns: 
        Saves a run info .txt file in the current output directory. 
        '''

        self.log.attempt("Attempting to save run information for a quick overview of runconditions")
        try:
            for key, value in experiment_dictionary.items():
                info_filename= f"{self.log.ulid}_-{key}_experiment_details.txt"
                with open(info_filename, "w") as file:
                    file.write(f"Key: {key}\n")
                    file.write(f"Allele: {value['allele']}\n")
                    file.write(f"Pairedness: {value['pairedness']}\n")
                    file.write(f"WT Files:\n")
                    for wt_file in value['wt']:
                        file.write(f"- {wt_file}\n")
                    file.write(f"Mu Files:\n")
                    for mu_file in value['mu']:
                        file.write(f"- {mu_file}\n")
                    file.write("\n")
                    file.write(f"vcf log path: {value['vcf_log_path']}")
                    file.write(f"analysis log path: {value['analysis_log_path']}")
            self.log.success("Experiment details saved successfully.")
        
        except Exception as e:
            self.log.fail(f"Error saving experiment details: {e}")

    def setup_directory(self, directory):
        '''
        Checks if the directory path given exists, and creates it if it doesn't.

        Args:
        directory(str) - path you wish to create if it doesn't exist

        Returns: 
        None. Creates directory if it doesn't exist. 
        '''
        self.log.attempt('Checking if directory exists...')
        try: 
            # Check if the output directory exists, and create it if necessary
            if not os.path.exists(directory):
                self.log.attempt(f"Directory does not exist. Creating: {output_dir}")
                os.makedirs(directory)
                self.log.success(f'Directory created: {directory}')
            else:
                self.log.note(f"Directory already exists: {directory}")
        
        except Exception as e:
            self.log.fail(f'setting up directory failed: {e}')

    
    def create_experiment_dictionary(self, line, **kwargs):
        approved_inputs = [
            'allele',
            'pairedness',
            'reference_genome_name',
            'reference_genome_source',
            'cleanup',
            'vcf_table_path',
            'threads_limit',
            'wt_input',
            'mu_input',
            'known_snps_path',
            'call_variants_in_parallel',
            'snpEff_species_db',
            'core_ulid',
            'vcf_ulid',
            'analysis_ulid'
        ]

        path_variables =[
            'vcf_table_path',
            'known_snps_path'
        ]
        
        details = {}

        if not line:
            self.log.fail(f'line name can not be empty - it is required to create an experiment dictionary! Aborting')

        for key, value in kwargs.items():
            if key in approved_inputs and path_variables:
                if os.path.exists(value):
                    details[key] = value  # Key-value pair in the subdictionary
                else:
                    self.log.note(f"path:{value} does not exist. Checking for it in {INPUT_DIR}..")
                    input_path = os.path.join(INPUT_DIR, value)
                    if os.path.exists(input_path):
                        self.log.note(f" path:{input_path} found! Assigning path value to dictionary")
                    else: 
                        self.log.fail(f'{value} not found in {INPUT_DIR} or as a hard coded path. Aborting')
            elif key in approved_inputs:
                details[key] = value
            else: 
                self.log.fail(f"{key} is not in the approved_inputs for experiment dictionary creation. Dictionary or arguments may need updating.")

        experiment_dict = {
            line : details
        }

    return experiment_dict

    def _extract_ulid_from_file_path(self, file_path):
        ulid_pattern = re.compile(r'[0-9A-HJKMNPQRSTVWXYZ]{26}')
        match = ulid_pattern.search(file_path)
        if match:
            return match.group()
        else:
            return None

    def generate_vcf_file_paths(self, current_line_name, vcf_log, known_snps):
        """
        Generate file paths based on the given parameters.

        Args:
        current_line_name (str): The name used as a key in the experiment_dictionary.
        vcf_log: Some object with an 'ulid' attribute.
        known_snps (str): The name of the known_snps file.

        Returns:
        Tuple containing the generated file paths.
        """
        # Add output_path to experiment_dictionary. 
        output_name_prefix = f"{vcf_log.ulid}_-{current_line_name}"
        output_dir_path = os.path.join(OUTPUT_DIR, output_name_prefix)
        output_prefix = os.path.join(output_dir_path, output_name_prefix)
        
        # Add vcftable_path to experiment_dictionary.
        vcf_table_name = f"{output_name_prefix}.noknownsnps.table"
        vcf_table_path = os.path.join(OUTPUT_DIR, vcf_table_name)
        
        # Generate VCFgen.sh script path
        vcfgen_script_path = os.path.join(MODULES_DIR, 'subprocess_VCFgen.sh')
        
        # Generate the knownSnps .vcf file path
        known_snps_path = os.path.join(REFERENCE_DIR, known_snps)
        
        return output_dir_path, output_prefix, vcf_table_path, vcfgen_script_path, known_snps_path

class ThaleBSASQLDB:
    """
    Handling retrieving and entry from database.
     IN PROGRESS....
     """
    def __init__(self, logger, db_name="thale_bsa_sqldb.db"):
        self.conn = sqlite3.connect(db_name)
        self.log = logger 

    def get_vcf_data(self, analysis_id):
        """Retrieve the VCF data based on the analysis ID"""
        cursor = self.conn.execute('''
            SELECT line_name, vcf_id, vcf_log_path, analysis_log_path, 
                run_date FROM thale_bsa_sqldb WHERE analysis_id = ?
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



