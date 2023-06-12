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
<<<<<<< HEAD
from scipy.interpolate import interp1d
=======
import plotly.io as pio
import scipy
from scipy.signal import find_peaks
>>>>>>> ee0054a (add plot production to .py script, add archive management and cleanup to vcfgen.sh script)
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
<<<<<<< HEAD


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
=======


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

for key in lines_dict:
	vcftable = "./output/archive/" + key/ + key + ".ems.table"
	plotname = "./output/archive/" + key/ + key + ".BSA_linkage_map" + ".png"
	finaltablename = "./output/archive/" + key/ + key + ".complete_table" + ".tsv"
	df = pd.read_csv(vcftable, sep="\t", header=0)

	#calculate ratio, add to dataframe and calculate mean ratio for peak threshold later. 
	df['ratio'] = ((df.wt_ref)/(df.wt_ref+df.wt_alt))-((df.mu_ref)/(df.mu_ref+df.mu_alt))
	ratio_mean = df["ratio"].mean()

	#fit LOWESS curve for each chromosome, add fitted values to dataframe and identify peaks
	chr_facets=df["chr"].unique()
	df_list = []
	lowess = sm.nonparametric.lowess

	for i in chr_facets:
    
    	df_chr = df[df['chr']==i].copy()
    
    	X=df_chr['pos'].values
    
    	Y=df_chr['ratio'].values
    
    	y_hat = lowess(Y,X, frac=0.29)[:,1]
    
    	df_chr['yhat'] = y_hat
    
    	df_list.append(df_chr)
    
    	signal = df_chr['yhat'].to_numpy().flatten()
    
    	peaks = scipy.signal.find_peaks(signal, height=ratio_mean, 
                            threshold=None, 
                            distance=None, 
                            prominence=None, 
                            width=None, 
                            wlen=None, 
                            rel_height=0.5, 
                            plateau_size=None)
    
    	h = peaks[1]['peak_heights']
    
    	if len(h) > 0:
        	h.sort()
        	max = h[0]
        	min = h[-1]
    	else:
        	max = 0
        	min = 0
    
    	df_chr['peak'] = [1 if (np.isclose(max, x) or (max < x)) and (np.isclose(min, x) or (min > x))
                      	else 0 for x in df_chr['yhat']]
    
    df = pd.concat(df_list)

	df_peaks = df.loc[df['peak'] == 1]

	df.to_csv('file_name.tsv', sep='\t')

	chr_facets_p=df_peaks["chr"].unique()

	fig = px.scatter(df, x=df['pos'], y=df['ratio'],
    	facet_col="chr",
    	opacity=0.8,
   		color_discrete_sequence=['goldenrod'],
    	trendline="lowess",
    	trendline_options=dict(frac=0.29),
    	trendline_color_override="blue")
            
	fig.add_trace(go.Scatter(x=[1],
    	y=[1],
    	mode='lines',
    	name='Lowess Fitted Ratio',
    	line=dict(color="blue")))
            
	fig.update_layout(dict(plot_bgcolor = 'white'))
	fig.update_xaxes(matches=None)
	fig.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True))
	fig.for_each_xaxis(lambda x: x.update(title = ''))
	fig.for_each_yaxis(lambda y: y.update(title = ''))

	fig.add_annotation(
    	showarrow=False,
    	xanchor='center',
    	xref='paper', 
    	x=0.5, 
    	yref='paper',
    	y=-0.12,
    	text='Position'
	)

	fig.add_annotation(
    	showarrow=False,
    	xanchor='center',
    	xref='paper', 
    	x=-0.065, 
    	yanchor='middle',
    	yref='paper',
    	y=0.6,
    	textangle=-90,
    	text='Ratio'
	)


	fig.update_xaxes(showgrid=True, 
	    gridwidth=0.5, 
	    gridcolor='lightgrey')

	fig.update_xaxes(showline=True, 
	    linewidth=1, 
	    linecolor='black')

	fig.update_xaxes(rangemode="tozero")

	fig.update_yaxes(showgrid=True, 
	    gridwidth=0.5, 
	    gridcolor='lightgrey')

	fig.update_yaxes(showline=True, 
	    linewidth=1, 
	    linecolor='black')

	fig.update_yaxes(range=[0, 0.8])

	fig.update_layout(title=dict(text="Linkage map",
	    font=dict(color='black')))

	fig.update_traces(marker=dict(size=1))

	fig['data'][0]['showlegend'] = True
	fig['data'][0]['name'] = 'Polymorphism ratio'

	max_pos=df_peaks['pos'].max()
	min_pos=df_peaks['pos'].min()
	max_yhat_pos=df_peaks.loc[df_peaks['yhat'] == df_peaks['yhat'].max()]['pos'].values[0]

	for i in chr_facets_p:
	    fig.add_vrect(x0=min_pos, x1=max_pos, col=i,
	        fillcolor="red", opacity=0.2, line_width=0)
	    fig.add_vline(x=max_yhat_pos, 
	                  line_dash="dot", 
	                  col=i, 
	                  line_width=1)
	    
	fig.add_trace(go.Scatter(x=[0,0], 
	    y=[0,0], 
	    mode='lines', 
	    line=dict(color='black', width=1, dash='dot'),
	    name='Identified Peak',))

	pio.write_image(fig, plotname, scale=6, width=1080, height=1080)
	df_peaks.to_csv(finaltablename, sep='\t')
	
# #add gstatistic


>>>>>>> ee0054a (add plot production to .py script, add archive management and cleanup to vcfgen.sh script)
