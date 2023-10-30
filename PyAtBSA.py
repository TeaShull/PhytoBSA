#!/usr/bin/env python

import sys
import os
import argparse
import subprocess

import glob
import fnmatch

# Print delimeter
print("""
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Welcome to PyAtBSA. Checking files for formating
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
""")

# Check if files are formatted properly, and if paired end or not
lines_dict = dict()

for file in os.listdir('./input/'):
	if ((fnmatch.fnmatch(file, '*_1*')) 
		and (fnmatch.fnmatch(file, '*wt.fq.gz*')) 
		or (fnmatch.fnmatch(file, '*mu.fq.gz*'))): 
		line = file.split("_")[0]
		lines_dict[line] = lines_dict.get(line, 0) + 1

	elif ((fnmatch.fnmatch(file, '*_2*')) 
		and (fnmatch.fnmatch(file, '*wt.fq.gz*')) 
		or (fnmatch.fnmatch(file, '*mu.fq.gz*'))):
		line = file.split("_")[0]
		lines_dict[line] = lines_dict.get(line, 0) + 1

	elif ((fnmatch.fnmatch(file, '*wt.fq.gz*')) 
		or (fnmatch.fnmatch(file, '*mu.fq.gz*')) 
		and not (fnmatch.fnmatch(file, '*_2*')) 
		and not (fnmatch.fnmatch(file, '*_1*'))):
		line = file.split(".")[0]
		lines_dict[line] = lines_dict.get(line, 0) + 1
	else:
		print("""
		check that your inputs are named properly, 
		or perhaps if there are spurious(or no) files
		in input. 
		
		paired-end file names must be formatted as follows: 
		<line>_1.wt.fq.gz 
		<line>_2.wt.fq.gz
		<line>_1.mu.fq.gz 
		<line>_2.mu.fq.gz  
		
		unpaired-end file names must be formatted as follows: 
		<line>.wt.fq.gz 
		<line>.mu.fq.gz
		""")
		quit()

#check for multiple lines to genotype. Relevant for SNP masking later
# if single_test >= 4 or pair_test >= 8:
# 	multiple_lines = True

##Read file names
files_fq = [os.path.basename(x) for x in glob.glob('./input/*')]


####### VCF file generation ########
vcf_generation_script = './code/VCFgen.sh'

reads = str()

for key in lines_dict:
	if lines_dict[key] == 4:
		reads = "paired-end"
	elif lines_dict[key] == 2:
		reads = "single-read"

	cmd = ['vcf_generation_script', key, reads]
	subprocess.run(cmd, text=True)

####### SNP mask generation ########
# Create snp mask (not ready yet)
#cmd = ['snp_mask.sh']
#subprocess.run(cmd text=True)
#files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]


####### Data analysis ########
analysis_script = './code/analysis.py'

for key in lines_dict:
	cmd = ['python', analysis_script, key]
	subprocess.run(cmd, text=True)
