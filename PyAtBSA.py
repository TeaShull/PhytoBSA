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
import pandas as pd
import numpy as np
import statsmodels.api as sm
import plotly.express as px
import plotly.graph_objects as go
from scipy.interpolate import interp1d
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

#for key in lines_dict:
#	count_vcfs
#	subprocess.call(['./code/VCFgen.sh', key, paired])


# files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]


#######data analysis########
print("""
    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
	         -. .-.   .-. .-.   .-. .-.   .  
	         ||\|||\ /|||\|||\ /|||\|||\ /|
	         |/ \|||\|||/ \|||\|||/ \|||\||
	         ~   `-~ `-`   `-~ `-`   `-~ `-
    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
	 Producing plots and identifying putative causal mutations
    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
	""")

# lowess = sm.nonparametric.lowess

# for key in lines_dict:
# 	vcftable = "./output/" + key + ".ems.table"
# 	plotname = "./output/" + key + ".rough_BSA_map_" + ".jpeg"
# 	print(vcftable)
# 	df = pd.read_csv(vcftable, sep="\t", header=0)
# 	print(df)
# 	df['ratio'] = ((df.wt_ref)/(df.wt_ref+df.wt_alt))-((df.mu_ref)/(df.mu_ref+df.mu_alt))
# #make figure, add scatter plot values and facetes
# 	fig = px.scatter(df, x=df['pos'], y=df['ratio'], facet_col="chr",
# 		opacity=0.8, color_discrete_sequence=['blue'])
# 	fig.update_layout(dict(plot_bgcolor = 'white'))
# 	fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')
# 	fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
# 	fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')
# 	fig.update_yaxes(showline=True, linewidth=1, linecolor='black')
# 	fig.update_layout(title=dict(text="Bulk segrigant linkage map, unsmoothed",
# 		font=dict(color='black')))
# 	fig.update_traces(marker=dict(size=3))


# #add LOWESS line
# 	chr_facets=df["chr"].unique()
# 	for i in chr_facets:
# 		X=df[df['chr']==i]['pos'].values
# 		print(X) 
# 		Y=df[df['chr']==i]['ratio'].values
# 		print(Y)
# 		y_hat = lowess(Y,X, frac=1/5)
# 		fig.add_traces(go.Scatter(x=y_hat[:,0], y=y_hat[:1], name = 'Lowess', line=dict(color='red'), facet_col=i))
# 		fig.write_image("plotname")
# #add gstatistic
