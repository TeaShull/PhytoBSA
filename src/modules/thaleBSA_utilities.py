import os
import subprocess
import fnmatch
import logging
import uuid
from datetime import datetime
from flask import session
from config import error_handler, INPUT_DIR, MODULES_DIR, LOG_DIR

#General Utilities Class
class ThaleBSAUtilities:

    @staticmethod
    def detect_file_type(file):
        error_handler('attempt', 
            f"Detecting if {file} is labeled recessive (R) or dominant (D)..."
        )
        try:
            if fnmatch.fnmatch(file, '*.R*'):
                label = "R"
            elif fnmatch.fnmatch(file, '*.D*'):
                label = "D"
            error_handler('success', f"{file} is labeled {label}")
            return label
        except Exception as e:
            error_handler('fail', f"Error while detecting file type: {e}")
            return None

    def extract_uuid(file_path):
        error_handler('attempt', "attempting to retrieve UUID from vcf table...")
        try:
            with open(file_path, 'r') as file:
                # Read the first line from the file
                first_line = file.readline().strip()

                # Check if the line contains '# UUID:'
                if '# UUID:' in first_line:
                    # Extract characters directly after '# UUID:'
                    uuid_value = first_line.split('# UUID:', 1)[1].strip()
                    return uuid_value
                else:
                    error_handler('fail', "No UUID found in the first line of the table.")
        except FileNotFoundError:
            error_handler('fail', f"File not found: {file_path}")
        except Exception as e:
            error_handler('fail', f"Error: {e}")

    def create_experiment_dictionary(self):
        error_handler('attempt', "Creating a dictionary to store experiment details...")
        try:
            input_dir = INPUT_DIR
            lines_dict = {}

            for file in os.listdir(input_dir):
                error_handler('attempt', f"Parsing {file}...")
                try:
                    key = file.split(".")[0]
                    lines_dict[key] = lines_dict.get(
                        key, {'count': 0, 'allele': None}
                    )
                    lines_dict[key]['count'] += 1
                    lines_dict[key]['allele'] = self.detect_file_type(file)
                    error_handler('success', 
                        f"{file} parsed and added to dictionary under key {key}"
                    )
                except Exception as e:
                    error_handler('fail', 
                        f"Error parsing {file} for the experiment dictionary: {e}"
                    )

            # Convert counts and types to a readable format
            for key, value in lines_dict.items():
                error_handler('attempt', "Making dictionary more readable")
                try:
                    if value['count'] == 4:
                        value['reads'] = "paired-end"
                    elif value['count'] == 2:
                        value['reads'] = "single-read"
                    error_handler('success', "Dictionary is now more readable")
                except Exception as e:
                    error_handler('fail', f"Error converting formatting counts and types: {e}")

            error_handler('success', "Experiment Dictionary Created ")
            return lines_dict

        except Exception as e:
            error_handler('fail', f"Error while creating the experiment dictionary: {e}")
            return {}

    def vcf_file_generation(self, experiment_dictionary,
                             reference_genome_name, snpEff_db_name,
                             reference_genome_source, threads_limit,
                             cleanup, known_snps
                             ):

        error_handler('attempt', 'Generating VCF files for experiments in dictionary')
        for key, value in experiment_dictionary.items():
            print(f"Generating VCF file for {key}...")
            try:
                # Construct cmd
                modules_dir = MODULES_DIR
                log_dir = LOG_DIR
                vcf_table_uuid = str(uuid.uuid4())
                vcfgen_script_path = os.path.join(modules_dir, 'VCFgen.sh')
                args = (key, value['reads'], value['allele'],
                        reference_genome_name, snpEff_db_name,
                        reference_genome_source, threads_limit,
                        cleanup, known_snps, vcf_table_uuid
                        )
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                log_name = f"{key}.VCF_file_generation.{timestamp}.{vcf_table_uuid}.log"
                log_path = os.path.join(log_dir, log_name)
                cmd = f"{vcfgen_script_path} {' '.join(map(str, args))} | tee {log_path}"

                # Run cmd
                subprocess.run(cmd, shell=True, text=True, check=True)

                error_handler('success', f"VCF file generated for {key}. Log saved to {log_path}")
            except Exception as e:
                error_handler('fail', f"Error while generating the VCF file for {key}: {e}")

        error_handler('success', "VCF file generation process complete")


    def data_analysis(self, experiment_dictionary):
        error_handler('attempt', "Attempting to perform data analysis...")
        try:
            modules_dir = MODULES_DIR
            log_dir = LOG_DIR
            analysis_uuid = str(uuid.uuid4())
            analysis_script = os.path.join(modules_dir, 'analysis.py')
            for key in experiment_dictionary:
                # Construct cmd
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                log_name = f"{key}.Data_Analysis_{timestamp}.{analysis_uuid}.log"
                log_path = os.path.join(log_dir, log_name)
                cmd = f'python {analysis_script} {key} {analysis_uuid} | tee {log_path}'
                # Run cmd
                subprocess.run(cmd, shell=True, text=True, check=True)

            error_handler('success', "Data analysis complete")
        except subprocess.CalledProcessError as e:
            error_handler('fail', f"Error during data analysis. Command returned non-zero exit code. Output: {e.output}")
        except Exception as e:
            error_handler('fail', f"Error during data analysis: {e}")