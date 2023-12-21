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
    def __init__(self, vcf_uuid, analysis_uuid):
        self.vcf_uuid
        self.analysis_uuid

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
    try:
        with open(file_path, 'r') as file:
            # Read the first line from the file
            first_line = file.readline().strip()

            # Check if the line starts with '#UUID'
            if first_line.startswith('#UUID'):
                # Extract characters directly after '#UUID'
                uuid_value = first_line[len('#UUID'):].strip()
                print(f"UUID value: {uuid_value}")
            else:
                print("No matching pattern found in the first line.")
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error: {e}")

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
                vcfgen_script_path = os.path.join(modules_dir, 'VCFgen.sh')
                args = (key, value['reads'], value['allele'],
                        reference_genome_name, snpEff_db_name,
                        reference_genome_source, threads_limit,
                        cleanup, known_snps
                        )
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                log_name = f"{key}.VCF_file_generation.{timestamp}.log"
                cmd = f"{vcfgen_script_path} {' '.join(map(str, args))} | tee {log_name}"

                # Run cmd
                subprocess.run(cmd, shell=True, text=True, check=True)

                error_handler('success', f"VCF file generated for {key}. Log saved to {log_name}")
            except Exception as e:
                error_handler('fail', f"Error while generating the VCF file for {key}: {e}")

        error_handler('success', "VCF file generation process complete")


    def data_analysis(self, experiment_dictionary, command_line):
        error_handler('attempt', "Attempting to perform data analysis...")
        try:
            modules_dir = MODULES_DIR
            analysis_script = os.path.join(modules_dir, 'analysis.py')
            for key in experiment_dictionary:
                # Construct cmd
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                log_name = f"{key}.Data_Analysis_{timestamp}.log"
                cmd = f'python {analysis_script} {key} | tee {log_name}'

                # Run cmd
                subprocess.run(cmd, shell=True, text=True, check=True)

            error_handler('success', "Data analysis complete")
        except subprocess.CalledProcessError as e:
            error_handler('fail', f"Error during data analysis. Command returned non-zero exit code. Output: {e.output}")
        except Exception as e:
            error_handler('fail', f"Error during data analysis: {e}")