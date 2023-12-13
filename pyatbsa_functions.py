#!/usr/bin/env python

import sys
import os
import argparse
import subprocess

import glob
import fnmatch

from flask import session

def create_experiment_dictionary():
	lines_dict = dict()
	input_dir = './input/'
	#Check if files are paired-end or single-read
	for file in os.listdir(input_dir):
		if ((fnmatch.fnmatch(file, '*_1*')) 
			and (fnmatch.fnmatch(file, '*wt.fq.gz*')) 
			or (fnmatch.fnmatch(file, '*mu.fq.gz*'))): 
			line = file.split(".")[0]
			lines_dict[line] = lines_dict.get(line, 0) + 1

		elif ((fnmatch.fnmatch(file, '*_2*')) 
			and (fnmatch.fnmatch(file, '*wt.fq.gz*')) 
			or (fnmatch.fnmatch(file, '*mu.fq.gz*'))):
			line = file.split(".")[0]
			lines_dict[line] = lines_dict.get(line, 0) + 1

		elif ((fnmatch.fnmatch(file, '*wt.fq.gz*')) 
			or (fnmatch.fnmatch(file, '*mu.fq.gz*')) 
			and not (fnmatch.fnmatch(file, '*_2*')) 
			and not (fnmatch.fnmatch(file, '*_1*'))):
			line = file.split(".")[0]
			lines_dict[line] = lines_dict.get(line, 1) + 1
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

	#convert dict int values to list, and replace pairedness 
	#quantification with a readable string
	
	for key, value in lines_dict.items():
		if lines_dict[key] == 4:
			lines_dict[key] = ["paired-end"]
		elif lines_dict[key] == 2:
			lines_dict[key] = ["single-read"]
		
	#Detect if files are labeled Recessive or Dominant. Add 
	#designation to the list. 
	marked_keys = {}
	
	for file in os.listdir(input_dir):
		for key in lines_dict:
			key_match = f"*{key}*"
			if ((fnmatch.fnmatch(file, key_match))
				and (fnmatch.fnmatch(file, '*.R*'))
				and key not in marked_keys):

				lines_dict[key].append("R")
				marked_keys[key] = True

			elif ((fnmatch.fnmatch(file, key_match))
				and (fnmatch.fnmatch(file, '*.D*'))
				and key not in marked_keys):
				
				lines_dict[key].append("D")
				marked_keys[key] = True

	return(lines_dict)

def vcf_file_generation():
	experiment_dictionary = session.get('experiment_dictionary', {})
	
	for key in experiment_dictionary:
		reads = experiment_dictionary[key][0]
		allele = experiment_dictionary[key][1]
		cmd = ['./code/VCFgen.sh', key, reads, allele]
		subprocess.run(cmd, text=True)

def data_analysis():
	experiment_dictionary = session.get('experiment_dictionary', {})
	
	for key in experiment_dictionary:
	  cmd = ['python', ./code/analysis.py, key]
	  subprocess.run(cmd, text=True)

# ####### SNP mask generation ########
# # Create snp mask (not ready yet)
# #cmd = ['snp_mask.sh']
# #subprocess.run(cmd text=True)
# #files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]
