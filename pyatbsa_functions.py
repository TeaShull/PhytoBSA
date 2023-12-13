#!/usr/bin/env python

import os
import subprocess
import fnmatch

from flask import session

def detect_file_type(file):
    """Detect if a file is labeled as Recessive (R) or Dominant (D)."""
    if fnmatch.fnmatch(file, '*.R*'):
        return "R"
    elif fnmatch.fnmatch(file, '*.D*'):
        return "D"
    return None

def create_experiment_dictionary():
    """Create a dictionary to store experiment details."""
    lines_dict = {}
    input_dir = './input/'

    for file in os.listdir(input_dir):
        key = file.split(".")[0]
        lines_dict[key] = lines_dict.get(key, {'count': 0, 'allele': None})
        lines_dict[key]['count'] += 1
        lines_dict[key]['allele'] = detect_file_type(file)

    # Convert counts and types to a readable format
    for key, value in lines_dict.items():
        if value['count'] == 4:
            value['reads'] = "paired-end"
        elif value['count'] == 2:
            value['reads'] = "single-read"

    return lines_dict

def vcf_file_generation(experiment_dictionary):
    """Generate VCF files based on the experiment details."""
    for key, value in experiment_dictionary.items():
        cmd = ['./code/VCFgen.sh', key, value['reads'], value['allele']]
        subprocess.run(cmd, text=True)

def data_analysis(experiment_dictionary):
    """Perform data analysis based on the experiment details."""
    for key in experiment_dictionary:
        cmd = ['python', './code/analysis.py', key]
        subprocess.run(cmd, text=True)

if __name__ == "__main__":
    experiment_dict = create_experiment_dictionary()
    session['experiment_dictionary'] = experiment_dict

    vcf_file_generation(experiment_dict)
    data_analysis(experiment_dict)


# ####### SNP mask generation ########
# # Create snp mask (not ready yet)
# #cmd = ['snp_mask.sh']
# #subprocess.run(cmd text=True)
# #files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]
