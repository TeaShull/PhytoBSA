#!/usr/bin/env python

import os
import sys
import time
import datetime
import argparse
import csv
import glob
import subprocess
import platform
import shlex
import pandas
import numpy
import matplotlib.pyplot as plt
import fnmatch

#welcome, print delimeter with delightful ascii DNA "art"
print("""
    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
	                    Welcome
	         -. .-.   .-. .-.   .-. .-.   .  
	         ||\|||\ /|||\|||\ /|||\|||\ /|
	         |/ \|||\|||/ \|||\|||/ \|||\||
	         ~   `-~ `-`   `-~ `-`   `-~ `-
    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
	          Checking files for formating
    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
	""")

#Check if files are formatted properly, and if paired end or not
single_test = 0
pair_test = 0

for file in os.listdir('./input/'):
	if (fnmatch.fnmatch(file, '*_1*')) and (fnmatch.fnmatch(file, '*wt.fq.gz*')) or (fnmatch.fnmatch(file, '*mu.fq.gz*')): 
		pair_test += 1
	elif (fnmatch.fnmatch(file, '*_2*')) and (fnmatch.fnmatch(file, '*wt.fq.gz*')) or (fnmatch.fnmatch(file, '*mu.fq.gz*')):
		pair_test += 1
	elif (fnmatch.fnmatch(file, '*wt.fq.gz*')) or (fnmatch.fnmatch(file, '*mu.fq.gz*')) and not (fnmatch.fnmatch(file, '*_2*')) or (fnmatch.fnmatch(file, '*_1*')):
		single_test += 1
	else:
		print("""
		check that your inputs are named properly, or perhaps if there are spurious(or no) files in 
		input. 
		
		file names must be in <line>_1.wt.fq.gz <line>_2.wt.fq.gz <line>_1.mu.fq.gz <line>_2.mu.fq.gz  
		format for paired end reads and <line>.wt.fq.gz <line>.mu.fq.gz for unpaired. 

		Currently,there is no support for batch producing VCFs for mixes of single-read and 
		paired-end files. This script won't do a good job detecting both. 
		""")
		quit()

if single_test >= 2:
    paired = 'single-read'

if pair_test >= 3:
	paired = 'paired-end'

#check for multiple lines to genotype. Relevant for SNP masking later
if single_test >= 4 or pair_test >= 8:
	multiple_lines = True

# #Read file names
files_fq = [os.path.basename(x) for x in glob.glob('./input/*')]

#create a list from the line names
lines = []
if paired == 'paired-end':
	for f in files_fq:
		lines.append(f.split("_")[0])
elif paired == 'single-read':
	for f in files_fq:
		lines.append(f.split(".wt")[0])

#remove duplicate entrys from paired reads and mu/wt designations using a dictionary
lines_dict = list(dict.fromkeys(lines))


#Iterate through detected files and produce VCFs
for key in lines_dict:
	subprocess.call(['./code/VCFgen.sh', key, paired])

for key in lines_dict:
	count_vcfs
	subprocess.call(['./code/VCFgen.sh', key, paired])


# files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]