from settings.paths import (INPUT_DIR, 
    OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR, LOG_DATABASE_NAME, LOG_DATABASE_PATH
)

import os
import re
import sqlite3


class FileUtilities:
    def __init__(self, logger):
        self.log = logger 

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
                self.log.note(f"Directory does not exist. Creating: {directory}")
                os.makedirs(directory)
                self.log.success(f'Directory created: {directory}')
            else:
                self.log.success(f"Directory already exists: {directory}")
        
        except Exception as e:
            self.log.fail(f'setting up directory failed: {e}')

    def process_path(self, directories: list, path: str) -> str:
        if os.path.exists(path):
            self.log.note(f'Path found and assigned: {path}')
            return path
        else:
            for directory in directories:
                self.log.note(f"path:{path} does not exist. Checking for it in {directory}..")
                dir_path = os.path.join(directory, path)
                if os.path.exists(dir_path):
                    self.log.note(f" path:{dir_path} found! Assigning path value.")
                    return dir_path

            self.log.fail(f'Path not found in any of the provided directories or the hard coded ({path}). Aborting')
            return None

    def extract_ulid_from_file_path(self, file_path):
        ulid_pattern = re.compile(r'[0-9A-HJKMNPQRSTVWXYZ]{26}')
        match = ulid_pattern.search(file_path)
        if match:
            return match.group()
        else:
            return None


class LogDbUtilites:
    def __init__(self):
        self.conn = sqlite3.connect(LOG_DATABASE_PATH)

    def print_analysis_log_data(self, ulid):
        """Retrieve the paths based on the analysis ID"""
        cursor = self.conn.execute('''
        SELECT * FROM analysis
            WHERE analysis_ulid = ? OR vcf_ulid = ? OR core_ulid = ?
        ''', (ulid, ulid, ulid))
        result = cursor.fetchone()

        if result:
            print(f"analysis_ulid: {result[0]}")
            print(f"analysis_log_path: {result[1]}")
            print(f"analysis_timestamp: {result[2]}")
            print(f"name: {result[3]}")
            print(f"core_ulid: {result[4]}")
            print(f"vcf_ulid: {result[5]}")
            print(f"ratio_cutoff: {result[6]}")
            print(f"loess_span: {result[7]}")
            print(f"smooth_edges_bounds: {result[8]}")
            print(f"filter_indels: {result[9]}")
            print(f"filter_ems: {result[10]}")
            print(f"snpmask_path: {result[11]}")
        else:
            print(f"No database entry found for {ulid}")

    def print_vcf_log_data(self, ulid):
        cursor = self.conn.execute('''
        SELECT * FROM vcf 
            WHERE vcf_ulid = ? OR core_ulid = ?
        ''', (ulid, ulid)) 
        result = cursor.fetchone()
        if result:
            print(f"vcf_ulid: {result[0]}")
            print(f"vcf_log_path: {result[1]}")
            print(f"vcf_timestamp: {result[2]}")
            print(f"name: {result[3]}")
            print(f"core_ulid: {result[4]}")
            print(f"reference_genome_path: {result[5]}")
            print(f"snpeff_species_db: {result[6]}")
            print(f"reference_genome_source: {result[7]}")
            print(f"threads_limit: {result[8]}")
        else:
            print(f"No database entry found for {ulid}")

    def print_line_name_data(self, line_name):
        """Retrieve all entries based on the line name"""
        cursor = self.conn.execute('''
        SELECT vcf.*, analysis.* 
            FROM vcf 
            LEFT JOIN analysis 
            ON vcf.name = analysis.name 
            WHERE vcf.name = ?
        ''', (line_name,))
        results = cursor.fetchall()

        if results:
            for result in results:
                print("\nVCF Data:")
                if result[0] is not None:
                    print(f"vcf_ulid: {result[0]}")
                    print(f"vcf_log_path: {result[1]}")
                    print(f"vcf_timestamp: {result[2]}")
                    print(f"name: {result[3]}")
                    print(f"core_ulid: {result[4]}")
                    print(f"reference_genome_path: {result[5]}")
                    print(f"snpeff_species_db: {result[6]}")
                    print(f"reference_genome_source: {result[7]}")
                    print(f"threads_limit: {result[8]}")
                else:
                    print("No VCF data found for this line name.")

                print("\nAnalysis Data:")
                if result[9] is not None:
                    print(f"analysis_ulid: {result[9]}")
                    print(f"analysis_log_path: {result[10]}")
                    print(f"analysis_timestamp: {result[11]}")
                    print(f"name: {result[12]}")
                    print(f"core_ulid: {result[13]}")
                    print(f"vcf_ulid: {result[14]}")
                    print(f"ratio_cutoff: {result[15]}")
                    print(f"loess_span: {result[16]}")
                    print(f"smooth_edges_bounds: {result[17]}")
                    print(f"filter_indels: {result[18]}")
                    print(f"filter_ems: {result[19]}")
                    print(f"snpmask_path: {result[20]}")
                else:
                    print("No analysis data found for this line name.")
                print("\n")  # for separating different entries
        else:
            print("No results found for this line name.")

    def print_core_ulid_data(self, core_ulid):
        """Retrieve all entries based on the core ulid"""
        cursor = self.conn.execute('''
        SELECT core.*, vcf.*, analysis.* 
            FROM core 
            LEFT JOIN vcf 
            ON core.core_ulid = vcf.core_ulid 
            LEFT JOIN analysis 
            ON core.core_ulid = analysis.core_ulid 
            WHERE core.core_ulid = ?
        ''', (core_ulid,))
        results = cursor.fetchall()

        if results:
            for result in results:
                print("\nCore Data:")
                print(f"core_ulid: {result[0]}")
                print(f"core_log_path: {result[1]}")
                print(f"core_timestamp: {result[2]}")

                if result[3] is not None:
                    print("\nVCF Data:")
                    print(f"vcf_ulid: {result[3]}")
                    print(f"vcf_log_path: {result[4]}")
                    print(f"vcf_timestamp: {result[5]}")
                    print(f"name: {result[6]}")
                    print(f"core_ulid: {result[7]}")
                    print(f"reference_genome_path: {result[8]}")
                    print(f"snpeff_species_db: {result[9]}")
                    print(f"reference_genome_source: {result[10]}")
                    print(f"threads_limit: {result[11]}")

                if result[12] is not None:
                    print("\nAnalysis Data:")
                    print(f"analysis_ulid: {result[12]}")
                    print(f"analysis_log_path: {result[13]}")
                    print(f"analysis_timestamp: {result[14]}")
                    print(f"name: {result[15]}")
                    print(f"core_ulid:{result[16]}")
                    print(f"vcf_ulid: {result[17]}")
                    print(f"ratio_cutoff: {result[18]}")
                    print(f"loess_span: {result[19]}")
                    print(f"smooth_edges_bounds: {result[20]}")
                    print(f"filter_indels: {result[21]}")
                    print(f"filter_ems: {result[22]}")
                    print(f"snpmask_path: {result[23]}")

                print("\n")  # for separating different entries
        else:
            print(f"No database entries found for core ULID: {core_ulid}")


