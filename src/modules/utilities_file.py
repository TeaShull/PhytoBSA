import os
import pandas as pd
import re
import sqlite3
from config import INPUT_DIR

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

    def save_experiment_details(self, experiment_dictionary):
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

    def setup_directory(self, output_dir):
        # Create the output directory path
        self.log.attempt('Checking if output directory exists...')
        try: 
            # Check if the output directory exists, and create it if necessary
            if not os.path.exists(output_dir):
                self.log.attempt(f"Directory does not exist. Creating: {output_dir}")
           
                os.makedirs(output_dir)
                self.log.success(f'Directory created: {output_dir}')
            else:
                self.log.note(f"Directory already exists: {output_dir}")
        except Exception as e:
            self.log.fail(f'setting up directory failed: {e}')


    def create_experiment_dictionary(self, line_name, vcf_table)->dict:
        '''
        Input: line_name and vcf_table.

        Returns: 
        experiment_dictionary:
            line_name
            vcf_table
            vcf_table_path
            core ulid
            vcf table ulid (if the file is labeled),
        '''
        
        self.log.attempt('Parsing arguments to create experiment_dictionary')
        try:
            self.log.note(f'line name parsed: {line_name}')
            self.log.note(f'vcf table parsed: {vcf_table}')
            
            core_ulid = self.log.ulid
            self.log.note(f'core_ulid:{core_ulid}')
                    
            vcf_table_path = os.path.join(INPUT_DIR, vcf_table)
            self.log.note(f'vcf table path created: {vcf_table_path}')

            vcf_ulid = self._extract_ulid_from_file_path(vcf_table_path)
            self.log.note(f'vcf ulid parsed:{vcf_ulid}')
                
        except Exception as e:
            self.log.fail(f'There was an error parsing the line name and vcf table path: {e}')
            quit()

        self.log.attempt(f'Checking if {vcf_table_path} exists...')
        try:
            if os.path.exists(vcf_table_path):
                self.log.success(f'Path exists: {vcf_table_path}')
                experiment_dictionary = {}
                experiment_dictionary[line_name] = {
                    'vcf_table_path': vcf_table_path,
                    'vcf_ulid': vcf_ulid,
                    'core_ulid' : core_ulid
                }
                return experiment_dictionary 
            else: 
                self.log.fail(f'vcf table path [{vcf_table_path}] does not exist.')
                quit()
        
        except Exception as e:
            self.log.fail(f'There was an error during experiment_dictionary creation:{e}')

    def _extract_ulid_from_file_path(self, file_path):
        ulid_pattern = re.compile(r'[0-9A-HJKMNPQRSTVWXYZ]{26}')
        match = ulid_pattern.search(file_path)
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


