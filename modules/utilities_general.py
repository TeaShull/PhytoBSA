from settings.paths import (INPUT_DIR, 
    OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR, LOG_DATABASE_NAME, LOG_DATABASE_PATH
)

import os
import re
import sqlite3

class FileUtilities:
    # Utility class: Handling file operations and inputs. 
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
                self.log.attempt(f"Directory does not exist. Creating: {directory}")
                os.makedirs(directory)
                self.log.success(f'Directory created: {directory}')
            else:
                self.log.note(f"Directory already exists: {directory}")
        
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
    # Utility class: retrieving log information from the log database
    def __init__(self):
        self.conn = sqlite3.connect(LOG_DATABASE_PATH)

    def print_analysis_log_data(self, ulid):
        """Retrieve the paths based on the analysis ID"""
        cursor = self.conn.execute('''
        SELECT analysis_ulid, line_name, core_ulid, 
            vcf_ulid, analysis_log_path, analysis_timestamp 
            FROM analysis
            WHERE analysis_ulid = ? OR vcf_ulid = ? OR core_ulid = ?
        ''', (ulid, ulid, ulid))
        result = cursor.fetchone()

        if result:
            print(f"analysis_ulid: {result[0]}")
            print(f"line_name: {result[1]}")
            print(f"core_ulid: {result[2]}")
            print(f"vcf_ulid: {result[3]}")
            print(f"analysis_log_path: {result[4]}")
            print(f"Analysis_timestamp: {result[5]}")
        else:
            print(f"No database entry found for {ulid}")

    def print_vcf_log_data(self, ulid):
        cursor = self.conn.execute('''
        SELECT vcf_ulid, line_name, core_ulid, vcf_log_path, vcf_timestamp 
            FROM vcf 
            WHERE vcf_ulid = ? OR core_ulid = ?
        ''', (ulid, ulid)) 
        result = cursor.fetchone()
        if result:
            print(f"VCF ulid: {result[0]}")
            print(f"Line Name: {result[1]}")
            print(f"Core ulid: {result[2]}")
            print(f"VCF Log Path: {result[3]}")
            print(f"VCF Timestamp: {result[4]}")
        else:
            print(f"No database entry found for {ulid}")

    def print_line_name_data(self, line_name):
        """Retrieve all entries based on the line name"""
        cursor = self.conn.execute('''
        SELECT vcf.vcf_ulid, vcf.vcf_log_path, vcf.vcf_timestamp, 
            analysis.analysis_ulid, analysis.line_name, analysis.core_ulid, 
            analysis.analysis_log_path, analysis.analysis_timestamp 
            FROM vcf 
            INNER JOIN analysis 
            ON vcf.line_name = analysis.line_name 
            WHERE vcf.line_name = ?
        ''', (line_name,))
        results = cursor.fetchall()

        if results:
            for result in results:
                print(f"VCF ULID: {result[0]}")
                print(f"VCF Log Path: {result[1]}")
                print(f"VCF Timestamp: {result[2]}")
                print(f"Analysis ULID: {result[3]}")
                print(f"Line Name: {result[4]}")
                print(f"Core ULID: {result[5]}")
                print(f"Analysis Log Path: {result[6]}")
                print(f"Analysis Timestamp: {result[7]}")
                print("\n")  # for separating different entries
        else:
            print("No results found for this line name.")

    def print_core_ulid_data(self, core_ulid):
        """Retrieve all entries based on the core ulid"""
        cursor = self.conn.execute('''
        SELECT vcf.vcf_ulid, vcf.vcf_log_path, vcf.vcf_timestamp, 
            analysis.analysis_ulid, analysis.line_name, analysis.core_ulid, 
            analysis.analysis_log_path, analysis.analysis_timestamp,
            core.core_log_path, core.core_timestamp
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
                print(f"Core ULID: {result[5]}")
                print(f"Core Log Path: {result[8]}")
                print(f"Core Timestamp: {result[9]}")

                if result[0] is not None:
                    print(f"VCF ULID: {result[0]}")
                    print(f"VCF Log Path: {result[1]}")
                    print(f"VCF Timestamp: {result[2]}")

                if result[3] is not None:
                    print(f"Analysis ULID: {result[3]}")
                    print(f"Line Name: {result[4]}")
                    print(f"Analysis Log Path: {result[6]}")
                    print(f"Analysis Timestamp: {result[7]}")

                print("\n")  # for separating different entries
        else:
            print(f"No database entries found for core ULID: {core_ulid}")


