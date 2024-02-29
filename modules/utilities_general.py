from settings.globals import (INPUT_DIR, 
    OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR, LOG_DATABASE_NAME, LOG_DATABASE_PATH
)

import os
import re
import sqlite3

import gzip
import shutil
import urllib.request
from urllib.parse import urlparse


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
            self.log.note(f"Path not hard coded. Traversing directories: {directories}")
            for directory in directories:
                for root, dirs, files in os.walk(directory):
                    self.log.note(f"Checking for path:{path} in {root}..")
                    
                    if path in files:
                        dir_path = os.path.join(root, path)
                        self.log.note(f"Path:{dir_path} found! Assigning path value.")
                        
                        return dir_path
            
            self.log.error(f'Path not found in any of the provided directories or the hard coded ({path}).')
            
            return None

    def extract_ulid_from_file_path(self, file_path):
        ulid_pattern = re.compile(r'[0-9A-HJKMNPQRSTVWXYZ]{26}')
        match = ulid_pattern.search(file_path)
        if match:
            return match.group()
        else:
            return None

    def _write_lines_class_attrs(self, lines_instance, file):
        for attr in dir(lines_instance):
            if not attr.startswith("__"):
                value = getattr(lines_instance, attr)
                if not callable(value):
                    file.write(f"  {attr}: {value}\n")

    def write_instance_vars_to_file(self, instances, filename):
        with open(filename, 'w') as f:
            for instance in instances:
                f.write(f"Instance of {instance.__class__.__name__}:\n")
                for attr in dir(instance):
                    if not attr.startswith("__"):
                        value = getattr(instance, attr)
                        if not callable(value):
                            f.write(f"  {attr}: {value}\n")
                f.write("\n")

    def _unzip_file(self, file_path: str) -> str:
        self.log.attempt(f"Attempting to unzip {file_path}")
        with gzip.open(file_path, 'rb') as f_in:
            with open(file_path[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        file_path = file_path[:-3]
        self.log.success(f"File unzipped to {file_path}")

        return file_path

    def _download_file(self, file_path: str, file_source: str) -> str:
        self.log.attempt(f"File doesn't exist, sourcing from: {file_source}")
        parsed_url = urlparse(file_source)
        source_extension = os.path.splitext(parsed_url.path)[1]
        
        if file_path.endswith(source_extension) is False:
            self.log.note(f"Source extension {source_extension} and file name extension don't match. Attempting to fix...")
            file_path += source_extension
        
        self.log.attempt(f"Downloading from URL: {file_source}")
        urllib.request.urlretrieve(file_source, file_path)
        self.log.success(f"File downloaded and saved at {file_path}")

        return file_path 

    def parse_file(self, file_path: str, file_source: str, destination: str)-> str:
        self.log.attempt("Attempting to parse file...")
        try:
            if not os.path.isfile(file_path) and file_source:
                file_name = os.path.basename(file_path)  # Extract the filename from the path
                file_path = os.path.join(destination, file_name)  # Reconstruct the path with the destination directory
                file_path = self._download_file(file_path, file_source)
            
            if os.path.isfile(file_path) and file_path.endswith('.gz'):
                file_path = self._unzip_file(file_path)

            if file_path:
                return file_path
            else:
                self.log.error('There was an error while parsing snpmask. Make sure the snpmask file exists or that the snpmask URL in references.db is valid.')
                self.log.error('You may need to reconfigure the entry for your reference in references.db. Sse the functions found in ./refdb_manager -h to do this')
        
        except Exception as e:
            self.log.error("Parsing file failed.")
            self.log.error(e)
            
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

    def get_line_name_data(self, line_name):
        """Retrieve all entries based on the line name"""
        cursor = self.conn.execute('''
        SELECT vcf.*, analysis.* 
            FROM vcf 
            LEFT JOIN analysis 
            ON vcf.name = analysis.name 
            WHERE vcf.name = ?
        ''', (line_name,))
        results = cursor.fetchall()
        return results

    def print_line_name_data(self, line_name):
        """Print all entries based on the line name"""
        results = self.get_line_name_data(line_name)
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

class RefDbUtilities:
    def __init__(self, logger, ref_name):
        self.log = logger
        self.ref_name = ref_name
        self.conn = None
        self.cursor = None
        self.db_path = 'references.db'
        self.open_connection()

    def open_connection(self):
        self.log.attempt("Opening reference database connection...")
        try:
            self.conn = sqlite3.connect(self.db_path)
            self.cursor = self.conn.cursor()
        except Exception as e:
            self.log.note(f"Failed to connect to the database: {e}")
            raise

    def fetch_ref_var(self, field):
        self.log.attempt(f"Attempting to retrieve {field} for {self.ref_name}")
        try:
            self.cursor.execute(f"""
                SELECT {field}
                FROM RefVariables
                WHERE reference_name = ?
            """, (self.ref_name,))
            
            row = self.cursor.fetchone()
            if row is None:
                return None
            
            self.log.success(f"{field} for {self.ref_name} retrieved!")

            return row[0]  # Return the single field value
        
        except Exception as e:
            self.log.fail(f"Failed to retrieve {field}: {e}")
            raise

    def _close_connection(self):
        if self.conn:
            self.conn.close()