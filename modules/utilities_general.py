from settings.config import INPUT_DIR, OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR

import os
import pandas as pd
import re
import sqlite3

class DictionaryUtilities:
    '''
    Utility class: Line dictionary creation and storage. Line_dict is a nested
    dictionary, with all details under the primary key of "line"
    
    All processes which alter line_dict occur within this class. 

    [NOTE] on the line_dict usage and structure
    I recognize that I could simply store the details as class variables, but I
    haven't found an nice way (i'm sure it exists) to iterate over "lines" for 
    bulk experiment runs. this structure allows the program to iterate over as 
    many "line_names" as are provided, facilitating fully automated bulk 
    processing of experiments. 


        [EXAMPLE line_dict]
        line_dict = {
            'line_name': {
                'segregation_type': '...',
                'wt_input': ['...', '...'],
                'mu_input': ['...', '...'],
                'pairedness': '...',
                ...
            }
        }
        '''
    def __init__(self, logger):
        self.log = logger

        self.line_dict={}

        self.approved_inputs = [
            'segregation_type',
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
        self.in_path_variables =[
            'vcf_table_path',
        ]

        self.ref_path_variables =[
            'known_snps_path'
        ]

    def experiment_detector(self) -> dict:
        '''
        Detects potential BSA experiments in the inputs folder if files are 
        named according to the conventions outlined in the README

        For paired-end:
        <line_name>.<R or D>_<read number>.<wt or mu>.fq.gz  
        
        For single-read:
        <line_name>.<R or D>.<wt or mu>.fq.gz

        Args: None

        Returns: line_dict containing the detected information. If
        successful, all information needed to generate vcf files and run analysis
        will have been parsed. 
        '''
        self.log.attempt(f"Detecting experiment details in: {INPUT_DIR}")
        try:
            for filename in os.listdir(INPUT_DIR):
                parts = filename.split('.')
                self.log.attempt(f'Parsing {filename}')
                line_name = parts[0]
                segregation_type = parts[1]
                
                # Processing segregation_type
                if '_1' in segregation_type:
                    segregation_type = segregation_type.rstrip('_1')
                if '_2' in segregation_type:
                    segregation_type = segregation_type.rstrip('_2')
                
                bulk_type = parts[-3]
                pairedness = 'paired-end' if '_1' or '_2' in filename else 'single-read'

                self.log.success(f"""{filename} parsed. 
                line_name:{line_name}
                segregation_type:{segregation_type}
                bulk_type:{bulk_type}
                pairedness:{pairedness}
                """)

                if line_name not in self.line_dict:
                    self.line_dict[line_name] = {
                        'segregation_type': segregation_type,
                        'wt_input': [],
                        'mu_input': [],
                        'pairedness': pairedness
                    }
                
                file_path = os.path.join(INPUT_DIR, filename)

                if bulk_type == 'wt':
                    self.line_dict[line_name]['wt_input'].append(file_path)
                elif bulk_type == 'mu':
                    self.line_dict[line_name]['mu_input'].append(file_path)
            
            for line_name, details in self.line_dict.items():
                details['wt_input'] = '" "'.join(sorted(details['wt_input']))
                details['mu_input'] = '" "'.join(sorted(details['mu_input']))

            self.log.success(f'Experiment dictionary generated.')

        except Exception as e:
            self.log.fail(f"Error while detecting experiment details: {e}")

    def populate_line_dict(self, line, **kwargs)-> dict:
        '''
        Organizes user inputs into a dictionary. 
        If the input needs to be a path, it checks if it exists. 
        If it doesn't exist, the function will check INPUT_DIR or REFERENCE_DIR 
        (depending if the path variable name is in "self.in_path_variables" or 
        "self.ref_path_variables" for the file and correct the pathing.

        input: line (which is the primary key for the dict, and is required)
        and user arguments passed as variables.

        output: line_dict, which organizes all experiment details
        under the line name as the primary key.

        Approved inputs and path variables need to be organized into the lists
        at the top of the function. Any time changes to inputs to the major
        functions of this program are made, organize them into this dictionary. 
        This way, all functions can flexably be assigned user inputs using an 
        easy to pass object. 
        '''
        details = {}
        self.log.attempt(f'Attempting to organize user inputs into an line dictionary')
        file_utils = self.FileUtilities(self.log)
        try:
            if not line:
                self.log.fail(f'line name can not be empty - it is required to create a line dictionary! Aborting')
            
            for key, details in kwargs.items():
                if key in self.approved_inputs and self.in_path_variables:
                    path = fileutils.process_path(INPUT_DIR, details[key])
                    details[key] = path
                
                elif key in self.approved_inputs and self.ref_path_variables:
                    path = fileutils.process_path(REFERENCE_DIR, details[key])
                
                elif key in self.approved_inputs:
                    details[key] = details
                    self.log.note(f'Arg to line dictionary| {key}:{details[key]}')
                
                else: 
                    self.log.fail(f"{key} is not in the approved_inputs for line dictionary creation. Dictionary or arguments may need updating.")

            if 'vcf_table_path' in kwargs.items() and 'vcf_ulid' not in kwargs.items():
                self.log.note('vcf table path provided, but vcf_ulid was not. Trying to extract ulid from file name...')
                details['vcf_ulid'] = file_utils.extract_ulid_from_file_path(vcf_table_path)

            self.line_dict = {
                line : details
            }
            self.log.success('User inputs successfully organized into line_dict.')
            return experiment_dict
            
        except Exception as e:
            self.log.fail(f'Error creating the line dictionary from user inputs {e}.')

    def check_vcf_gen_variables(self)-> 
        """
        Line dictionary will be populated with the contents of vcf_gen_variables, 
        if the user has not assigned any variables. 

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

    def generate_output_file_paths(self):
        """
        Generate file paths based on the given parameters.

        Args:
        current_line_name (str): The name used as a key in the self.line_dict.
        vcf_log: Some object with an 'ulid' attribute.
        known_snps (str): The name of the known_snps file.

        Returns:
        Tuple containing the generated file paths.
        """
        # Add output_path to self.line_dict. 
        for line, value in self.line_dict:
            output_name_prefix = f"{self.log.ulid}_-{line}"
            output_dir_path = self.process_path(OUTPUT_DIR, output_name_prefix)
            output_prefix = os.path.join(output_dir_path, output_name_prefix)
            value['output_dir_path'] = output_dir_path
            value['output_prefix'] = output_prefix
        
        # Add vcftable_path to self.line_dict.
        if value['vcf_table_path'] is None:
            vcf_table_name = f"{output_name_prefix}.noknownsnps.table"
            value['vcf_table_path'] = self.process_path(OUTPUT_DIR, vcf_table_name)
   
    def save_experiment_details(self):
        '''
        Saves self.line_dict information into a human-readable file, for
        easy veiwing. Allows easy access to information without searching through 
        logs. 

        Args: 
        self.line_dict(dict) - experiment dictionary that contains experiment info. 

        Returns: 
        Saves a run info .txt file in the current output directory. 
        '''

        self.log.attempt("Attempting to save run information for a quick overview of runconditions")
        try:
            for key, value in self.line_dict.items():
                info_filename= f"{self.log.ulid}_-{key}_experiment_details.txt"
        
                with open(info_filename, "w") as file:
                    file.write(f"Key: {key}\n")
                    file.write(f"Allele: {value['allele']}\n")
                    file.write(f"Pairedness: {value['pairedness']}\n")
                    file.write(f"WT Files:\n")
        
                    for wt_file in value['wt_input']:
                        file.write(f"- {wt_file}\n")
                    file.write(f"Mu Files:\n")
        
                    for mu_file in value['mu_input']:
                        file.write(f"- {mu_file}\n")
                    file.write("\n")
                    file.write(f"vcf log path: {value['vcf_log_path']}")
                    file.write(f"analysis log path: {value['analysis_log_path']}")
        
            self.log.success("Experiment details saved successfully.")
        
        except Exception as e:
            self.log.fail(f"Error saving experiment details: {e}")

class FileUtilities:
    # Utility class: Handling file operations and inputs. 
    def __init__(self, logger):
        self.log = logger 

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
                self.log.attempt(f"Directory does not exist. Creating: {directory}")
                os.makedirs(directory)
                self.log.success(f'Directory created: {directory}')
            else:
                self.log.note(f"Directory already exists: {directory}")
        
        except Exception as e:
            self.log.fail(f'setting up directory failed: {e}')

    def process_path(self, directory: str, path: str) -> str:
        if os.path.exists(path):
            self.log.note(f'Path found and assigned: {path}')
            return path
        else:
            self.log.note(f"path:{path} does not exist. Checking for it in {directory}..")
            dir_path = os.path.join(directory, path)
            if os.path.exists(input_path):
                self.log.note(f" path:{dir_path} found! Assigning path value.")
                return dir_path
            else:
                self.log.fail(f'path not found as {dir_path} or the hard coded ({path}). Aborting')
                return None

    def extract_ulid_from_file_path(self, file_path):
        ulid_pattern = re.compile(r'[0-9A-HJKMNPQRSTVWXYZ]{26}')
        match = ulid_pattern.search(file_path)
        if match:
            return match.group()
        else:
            return None

class LogDbUtilites:
    # Utility class: retrieving log information from the log database
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

        if result
            analysis_log_path, vcf_log_pat = result
            print(f"Analysis ID: {analysis_id}")
            print(f"Analysis Log Path: {analysis_log_path}")
            print(f"VCF Log Path: {vcf_log_path}")
        ele:
            print(f"No record found for Analysis ID '{analysis_id}'")



