import os
import pandas as pd
import re
import sqlite3
from config import ExperimentDictionary, INPUT_DIR

class FileUtilities:
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
        self.log.attempt(f"Detecting experiment details in: {INPUT_DIR}...")
        try:
            expt_dict = ExperimentDictionary()

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

    def check_vcfgen_variables(self, reference_genome_name, snpEff_species_db, reference_genome_source, threads_limit, cleanup, known_snps)
        # Check if the user has provided all the necessary variables
        self.log.attempt('Checking if runtime variables for VCFgen.sh subprocess are assigned...')
        try:
            if (
                reference_genome_name is None or
                snpEff_species_db is None or
                reference_genome_source is None or
                threads_limit is None or
                cleanup is None or
                known_snps is None
            ):
                # User hasn't provided all the variables, print an error message or handle accordingly
                self.log.fail("Error: Not all required variables have been provided.")
                return False
            else:
                # All variables are provided, continue with the rest of your code
                self.log.success("All variables are were provided. Well done! Proceeding...")
                return True
        except Exception as e:
            self.log.fail(f'There was an error while checking if variables for VCFgen.sh have been assigned:{e}')

    def load_vcf_table(
        self, current_line_table_path, current_line_name
    )->pd.DataFrame:
    
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


    def create_experiment_dictionary(self, line_name, vcf_table)->dict:
        '''

        Creates an experiment dictionary from line_name and vcf_table input. 
        Used to create experiment dictionaries when automatic experiment 
        detection is not initiated. 

        Args:
        line_name(str) 
        vcf_table(str)

        Returns: 
        experiment_dictionary(dict)-
            line_name:
                [line_name][vcf_table]:
                [line_name][vcf_table_path]:
                [line_name][core ulid]:
                [line_name][vcf table ulid]: (if the file is labeled),
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
                experiment_dictionary = ExperimentDictionary()
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
    """Handling retrieving and entry from database.
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


