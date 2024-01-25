from settings.config import INPUT_DIR, OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR

import os
import pandas as pd
import re
import sqlite3

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
        self.log.attempt(f'Checking if directory exists:{directory}')
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
            if os.path.exists(dir_path):
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

        if result:
            analysis_log_path, vcf_log_pat = result
            print(f"Analysis ID: {analysis_id}")
            print(f"Analysis Log Path: {analysis_log_path}")
            print(f"VCF Log Path: {vcf_log_path}")
        else:
            print(f"No record found for Analysis ID '{analysis_id}'")



