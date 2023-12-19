#!/usr/bin/env python

import os
import subprocess
import fnmatch
from flask import session
import config

class ThaleBSAUtilities:

    @staticmethod
    def detect_file_type(file):
        """Detect if a file is labeled as Recessive (R) or Dominant (D)."""
        print(f"(Attempt) Detecting if {file} is labeled recessive (R) or dominant (D)...")
        try:
            if fnmatch.fnmatch(file, '*.R*'):
                label = "R"
            elif fnmatch.fnmatch(file, '*.D*'):
                label = "D"
            print(f"(Success) {file} is labeled {label}")
            return label
        except Exception as e:
            print(f"(Fail) Error while detecting file type: {e}")
            return None

    def create_experiment_dictionary(self):
        """Create a dictionary to store experiment details."""
        print(f"(Attempt) Creating a dictionary to store experiment details...")
        try:
            input_dir = config.INPUT_DIR
            lines_dict = {}

            for file in os.listdir(input_dir):
                print(f"(Attempt) Parsing {file}...")
                try:
                    key = file.split(".")[0]
                    lines_dict[key] = lines_dict.get(
                        key, {'count': 0, 'allele': None}
                    )
                    lines_dict[key]['count'] += 1
                    lines_dict[key]['allele'] = self.detect_file_type(file)
                    print(f"(Success) {file} parsed and added to dictionary under key {key}")
                except Exception as e:
                    print(f"Error parsing {file} for the experiment dictionary: {e}")
                num_keys = len(lines_dict)
                print(f"(Success) Dictionary created. {num_keys} experiments detected. ")
            # Convert counts and types to a readable format
            for key, value in lines_dict.items():
                print(f"(Attempt) Making dictionary more readable")
                try:
                    if value['count'] == 4:
                        value['reads'] = "paired-end"
                    elif value['count'] == 2:
                        value['reads'] = "single-read"
                    print(f"(Success) Dictionary is now more readable")
                except Exception as e:
                    print(f"(Fail) Error converting formatting counts and types: {e}")
            return lines_dict

            print(f"(Success) Experiment Dictionary Created ")
        except Exception as e:
            print(f"(Fail) Error while creating the experiment dictionary: {e}")
            return {}

    def vcf_file_generation(self, experiment_dictionary,
                             reference_genome_name, snpEff_db_name, 
                             reference_genome_source, threads_limit, 
                             cleanup, known_snps
                             ):
        """Generate VCF files based on the experiment details."""
        print("(Attempt) Generating VCF files for experiments...")
        try:
            for key, value in experiment_dictionary.items():
                try:
                    print(f" (Attempt) Generating VCF file for {key}...")
                    vcfgen_script_path = os.path.join(os.getcwd(), 'VCFgen.sh')
                    args = (key, value['reads'], value['allele'],
                             reference_genome_name, snpEff_db_name, 
                             reference_genome_source, threads_limit, 
                             cleanup, known_snps
                             )
                    cmd = [vcfgen_script_path, *args]
                    subprocess.run(cmd, text=True)
                    print(f"(Success) VCF file generated for {key}")
                except Exception as e:
                    print(f"(Fail) Error while generating the VCF file for {key}: {e}")
            print("(Success) VCF file generation process complete")

        except Exception as e:
            print(f"(Fail) Error during VCF file generation: {e}")

    def data_analysis(self, experiment_dictionary, command_line):
        """Perform data analysis based on the experiment details."""
        print(f" (Attempt) Attempting to perform data analysis...")
        try:
            if not command_line:
                experiment_dictionary = session.get('experiment_dictionary', {})
            modules_dir = config.MODULES_DIR
            analysis_script = os.path.join(modules_dir, 'analysis.py')
            for key in experiment_dictionary:
                cmd = ['python', analysis_script, key]
                subprocess.run(cmd, text=True)
            print(f"(Success) Data analysis complete")
        except Exception as e:
            print(f"(Fail) Error during data analysis: {e}")


# ####### SNP mask generation ########
# # Create snp mask (not ready yet)
# # cmd = ['snp_mask.sh']
# # subprocess.run(cmd text=True)
# # files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]
