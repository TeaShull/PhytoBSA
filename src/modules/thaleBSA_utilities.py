#!/usr/bin/env python

import os
import subprocess
import fnmatch
from flask import session
import config

class ThaleBSAUtilities:
    def __init__(self):
        self.lines_dict = {}

    @staticmethod
    def detect_file_type(file):
        """Detect if a file is labeled as Recessive (R) or Dominant (D)."""
        if fnmatch.fnmatch(file, '*.R*'):
            return "R"
        elif fnmatch.fnmatch(file, '*.D*'):
            return "D"
        return None

    def create_experiment_dictionary(self):
        """Create a dictionary to store experiment details."""
        input_dir = config.INPUT_DIR

        for file in os.listdir(input_dir):
            key = file.split(".")[0]
            self.lines_dict[key] = self.lines_dict.get(key, {'count': 0, 'allele': None})
            self.lines_dict[key]['count'] += 1
            self.lines_dict[key]['allele'] = self.detect_file_type(file)

        # Convert counts and types to a readable format
        for key, value in self.lines_dict.items():
            if value['count'] == 4:
                value['reads'] = "paired-end"
            elif value['count'] == 2:
                value['reads'] = "single-read"

        return self.lines_dict

    def vcf_file_generation(self):
        """Generate VCF files based on the experiment details."""
        experiment_dictionary = session.get('experiment_dictionary', {})
        for key, value in experiment_dictionary.items():
            vcfgen_script_path = os.path.join(os.getcwd(), 'VCFgen.sh')
            cmd = [vcfgen_script_path, key, value['reads'], value['allele']]
            subprocess.run(cmd, text=True)

    def data_analysis(self):
        """Perform data analysis based on the experiment details."""
        experiment_dictionary = session.get('experiment_dictionary', {})
        modules_dir = config.MODULES_DIR
        analysis_script = os.path.join(modules_dir,'analysis.py')
        for key in experiment_dictionary:
            cmd = ['python', analysis_script, key]
            subprocess.run(cmd, text=True)


# ####### SNP mask generation ########
# # Create snp mask (not ready yet)
# #cmd = ['snp_mask.sh']
# #subprocess.run(cmd text=True)
# #files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]
