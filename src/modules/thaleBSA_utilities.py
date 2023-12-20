import os
import subprocess
import fnmatch
import logging
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
        try:
            for key, value in experiment_dictionary.items():
                try:
                    self.vcf_logger.info(f"Generating VCF file for {key}...")
                    modules_dir = MODULES_DIR
                    vcfgen_script_path = os.path.join(modules_dir, 'VCFgen.sh')
                    args = (key, value['reads'], value['allele'],
                            reference_genome_name, snpEff_db_name,
                            reference_genome_source, threads_limit,
                            cleanup, known_snps
                            )
                    cmd = [vcfgen_script_path, *args]

                    # Redirect stdout to the logging system
                    process = subprocess.run(
                        cmd, text=True, stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT, check=True
                    )
                    
                    # Save the output to a log file named with current_line_name and timestamp
                    current_line_name = key.replace(".", "_")
                    timestamp_output = datetime.now().strftime("%Y%m%d_%H%M%S")
                    log_dir = LOG_DIR
                    log_name = f'vcf_generation_{current_line_name}_{timestamp_output}.log'
                    log_path = os.path.join(log_dir, log_name)
                    with open(log_path, 'w') as log_file:
                        log_file.write(process.stdout)
                    
                    # Log the event with current_line_name and timestamp
                    error_handler('success', f"VCF file generated for {key}. log saved to {log_name_output}")
                except Exception as e:
                    self.vcf_logger.error(f"Error while generating the VCF file for {key}: {e}")

            error_handler('success', "VCF file generation process complete")

        except Exception as e:
            error_handler('success', f"Error during VCF file generation: {e}")

    def data_analysis(self, experiment_dictionary, command_line):
        error_handler('attempt', "Attempting to perform data analysis...")
        try:
            modules_dir = MODULES_DIR
            analysis_script = os.path.join(modules_dir, 'analysis.py')
            for key in experiment_dictionary:
                cmd = ['python', analysis_script, key]

                # Redirect stdout to the logging system
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                data_analysis_log_name = f'{key}.Data_Analysis_{timestamp}.log'
                log_dir = LOG_DIR
                log_file_path = os.path.join(log_dir, data_analysis_log_name)
                
                process = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, check=True
                )

                output = process.stdout

                for line in process.stdout.splitlines():
                # Print the line to stdout
                    print(line)

                # Log the line
                    logging.info(line)

                    with open(log_file_path, 'w') as log_file:
                        log_file.write(output)

                # Log the event
                error_handler('success', 
                    f"Script execution complete. Output saved to {log_file_path}"
                    )

            error_handler('success', "Data analysis complete")
        except Exception as e:
            error_handler('fail', f"Error during data analysis: {e}")

