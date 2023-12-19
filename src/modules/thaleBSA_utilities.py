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
        print(f"Detecting if {file} is labeled recessive (R) or dominant (D)...")
        try:
            if fnmatch.fnmatch(file, '*.R*'):
                label = "R"
            elif fnmatch.fnmatch(file, '*.D*'):
                label = "D"
            print(f"{file} is labeled {label}")
            return label
        except Exception as e:
            print(f"Error while detecting file type: {e}")
            return None

    def create_experiment_dictionary(self):
        """Create a dictionary to store experiment details."""
        print("Creating a dictionary to store experiment details...")
        try:
            input_dir = config.INPUT_DIR
            lines_dict = {}

            for file in os.listdir(input_dir):
                print(f"Parsing {file}...")
                try:
                    key = file.split(".")[0]
                    self.lines_dict[key] = self.lines_dict.get(
                        key, {'count': 0, 'allele': None}
                    )
                    self.lines_dict[key]['count'] += 1
                    self.lines_dict[key]['allele'] = self.detect_file_type(file)
                    print(f"{file} parsed and added to dictionary under key {key}")
                except Exception as e:
                    print(f"Error parsing {file} for the experiment dictionary: {e}")

            # Convert counts and types to a readable format
            for key, value in self.lines_dict.items():
                print(f"Making dictionary more readable")
                try:
                    if value['count'] == 4:
                        value['reads'] = "paired-end"
                    elif value['count'] == 2:
                        value['reads'] = "single-read"
                except Exception as e:
                    print(f"Error converting formatting counts and types: {e}")
            return self.lines_dict

        except Exception as e:
            print(f"Error while creating the experiment dictionary: {e}")
            return {}

    def vcf_file_generation(self, experiment_dictionary,
                             reference_genome_name, snpEff_db_name, 
                             reference_genome_source, threads_limit, 
                             cleanup=True, known_snps
                             ):
        """Generate VCF files based on the experiment details."""
        print("Attempting to generate VCF files for experiments...")
        try:
            for key, value in experiment_dictionary.items():
                try:
                    print(f"Generating VCF file for {key}...")
                    vcfgen_script_path = os.path.join(os.getcwd(), 'VCFgen.sh')
                    args = (key, value['reads'], value['allele'],
                             reference_genome_name, snpEff_db_name, 
                             reference_genome_source, threads_limit, 
                             cleanup, known_snps
                             )
                    cmd = [vcfgen_script_path, *args]
                    subprocess.run(cmd, text=True)
                except Exception as e:
                    print(f"Error while generating the VCF file for {key}: {e}")
            print("VCF file generation process complete")

        except Exception as e:
            print(f"Error during VCF file generation: {e}")

    def data_analysis(self, experiment_dictionary, command_line=False):
        """Perform data analysis based on the experiment details."""
        print(f"Attempting to perform data analysis...")
        try:
            if not command_line:
                experiment_dictionary = session.get('experiment_dictionary', {})
            modules_dir = config.MODULES_DIR
            analysis_script = os.path.join(modules_dir, 'analysis.py')
            for key in experiment_dictionary:
                cmd = ['python', analysis_script, key]
                subprocess.run(cmd, text=True)
        except Exception as e:
            print(f"Error during data analysis: {e}")


# ####### SNP mask generation ########
# # Create snp mask (not ready yet)
# # cmd = ['snp_mask.sh']
# # subprocess.run(cmd text=True)
# # files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]
