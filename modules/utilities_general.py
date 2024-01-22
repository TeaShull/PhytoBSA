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

    def experiment_detector(self) -> dict:
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

                key = line_name
                self.log.success(f"""{filename} parsed. 
                key:{key}
                allele:{allele}
                segregation_type:{segregation_type}
                pairedness:{pairedness}
            """)
                if key not in expt_dict:
                    expt_dict[key] = {
                        'allele': allele,
                        'wt_input': [],
                        'mu_input': [],
                        'pairedness': pairedness
                    }
                file_path = os.path.join(INPUT_DIR, filename)

                if 'wt' in segregation_type or 'mu' in segregation_type:
                    bulk_list = expt_dict[key].get(segregation_type, [])
                    bulk_list.append(file_path)
                    expt_dict[key][segregation_type + "_input"] = sorted(
                        bulk_list, key=lambda x: int(x.split('_')[-1][0])
                    )
                    input_files = ' '.join(expt_dict[key][segregation_type + "_input"])
                    expt_dict[key][segregation_type + "_input"] = f'"{input_files}"'
                    del expt_dict[key][segregation_type]

            self.log.success(f'Experiment dictionary generated.')
            return expt_dict

            except Exception as e:s
                self.log.fail(f"Error while detecting experiment details: {e}")
                return {}


        except Exception as e:
            self.log.fail(f"Error while detecting experiment details: {e}")
            return {}

    def check_vcfgen_variables(self, reference_genome_name=None, snpEff_species_db=None, reference_genome_source=None, threads_limit=None, cleanup=None, known_snps=None):
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
            for line_name, value in experiment_dictionary.items():
                info_filename= f"{self.log.ulid}_-{line_name}_experiment_details.txt"
                with open(info_filename, "w") as file:
                    file.write(f"Key: {line_name}\n")
                    file.write(f"Allele: {value['segregation_type']}\n")
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

    def generate_output_file_paths(self, experiment_dictionary):
        """
        Generate file paths based on the given parameters.

        Args:
        current_line_name (str): The name used as a line_name in the experiment_dictionary.
        vcf_log: Some object with an 'ulid' attribute.
        known_snps (str): The name of the known_snps file.

        Returns:
        Tuple containing the generated file paths.
        """
        # Add output_path to experiment_dictionary. 
        for line, value in experiment_dictionary:
            output_name_prefix = f"{self.log.ulid}_-{line}"
            output_dir_path = self.process_path(OUTPUT_DIR, output_name_prefix)
            output_prefix = os.path.join(output_dir_path, output_name_prefix)
            value['output_dir_path'] = output_dir_path
            value['output_prefix'] = output_prefix
        
        # Add vcftable_path to experiment_dictionary.
        if value['vcf_table_path'] is None:
            vcf_table_name = f"{output_name_prefix}.noknownsnps.table"
            value['vcf_table_path'] = self.process_path(OUTPUT_DIR, vcf_table_name)
       
        return experiment_dictionary

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



