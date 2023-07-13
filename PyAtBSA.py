#!/usr/bin/env python

import os
import sys
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
import plotly.io as pio
import scipy
from scipy.signal import find_peaks
import fnmatch
import io

#Print delimeter
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

#Check if files are formatted properly, and if paired end or not
pair_test = 0

#probably a more straightforward and safer way of doing this using fq headers...
for file in os.listdir('./input/'):
	if (fnmatch.fnmatch(file, '*_1*')) and (fnmatch.fnmatch(file, '*wt.fq.gz*')) or (fnmatch.fnmatch(file, '*mu.fq.gz*')): 
		pair_test += 1
	elif (fnmatch.fnmatch(file, '*_2*')) and (fnmatch.fnmatch(file, '*wt.fq.gz*')) or (fnmatch.fnmatch(file, '*mu.fq.gz*')):
		pair_test += 1
	elif (fnmatch.fnmatch(file, '*wt.fq.gz*')) or (fnmatch.fnmatch(file, '*mu.fq.gz*')) and not (fnmatch.fnmatch(file, '*_2*')) or (fnmatch.fnmatch(file, '*_1*')):
		pair_test += 3
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

if pair_test%2 == 0:
	paired = 'paired-end'

if pair_test%3 == 0:
	paired = 'single-read'

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

#Create snp mask (not ready yet)
# for key in lines_dict:
# 	count_vcfs
# 	subprocess.call(['./code/VCFgen.sh', key, paired])

# files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]


#######data analysis########
print("""
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Producing plots and identifying putative causal mutations
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
""")

######establish functions, exc######
def gStatistic_Array(o1, o3, o2, o4):
    # Calculate G-statistc using numpy array input
    np.seterr(all='ignore')

    e1 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o1+o3)/(o1+o2+o3+o4), 0)
    e2 = np.where(o1+o2+o3+o4!=0, (o1+o2)*(o2+o4)/(o1+o2+o3+o4), 0)
    e3 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o1+o3)/(o1+o2+o3+o4), 0)
    e4 = np.where(o1+o2+o3+o4!=0, (o3+o4)*(o2+o4)/(o1+o2+o3+o4), 0)

    llr1 = np.where(o1/e1>0, 2*o1*np.log(o1/e1), 0.0)
    llr2 = np.where(o2/e2>0, 2*o2*np.log(o2/e2), 0.0)
    llr3 = np.where(o3/e3>0, 2*o3*np.log(o3/e3), 0.0)
    llr4 = np.where(o4/e4>0, 2*o4*np.log(o4/e4), 0.0)

    return np.where(e1*e2*e3*e4==0, 0.0, llr1+llr2+llr3+llr4)

lowess = sm.nonparametric.lowess

#######Establish variables######

for key in lines_dict:
	vcftable = os.path.join("./output/" + key + "/" + key + ".noknownsnps.table")
	lowess_plotname = os.path.join("./output/" + key + "/" + key + ".BSA_linkage_map" + ".png")
	gstat_plotname = os.path.join("./output/" + key + "/" + key + ".BSA_g-stat" + ".png")
	finalcompletetablename = os.path.join("./output/" + key + "/" + key + ".complete_table" + ".tsv")
	finallikelycandidatestablename = os.path.join("./output/" + key + "/" + key + ".likely_candidates" + ".tsv")

	#read in ems table

	df = pd.read_csv(vcftable, sep="\t")

	#Calculate delta SNP ratio
	df['ratio'] = ((df.wt_ref)/(df.wt_ref+df.wt_alt))-((df.mu_ref)/(df.mu_ref+df.mu_alt))

	#Drop indels, leave only SNPs
	df[
	    df["ref"].apply(lambda x: len(x) == 1)
	    & df["alt"].apply(lambda x: len(x) == 1)
	]

	#Calculate ratio_mean for peak identification lower limit
	ratio_mean = df["ratio"].mean()

	#df.drop(df[df['ratio'] <= 0.3].index, inplace = True) # drop ratios below the defined value. Sometimes necissary

	#Drop NAs. Infrequently will arise while div/0 seemingly
	df.dropna(axis=0, how='any', subset="ratio", inplace=True)

	#######lowess smoothing by chromosome######
	#at some point, try not going chromosome by chromosome. May help with Lowess edge bias? 
	chr_facets=df["chr"].unique()
	df_list = []

	lowess = sm.nonparametric.lowess
	lowess_span=0.275

	for i in chr_facets:
	    
	    # Fit LOWESS ratio ~ position
	    df_chr = df[df['chr']==i].copy()

	    X=df_chr['pos'].values
	    
	    Y=df_chr['ratio'].values
	    
	    y_hat = lowess(Y,X, frac=lowess_span)[:,1]
	    
	    df_chr['yhat'] = y_hat # fitted values
	    
	    df_list.append(df_chr)
	    
	    # Find peaks in the LOWESS fitted curve
	    signal = df_chr['yhat'].to_numpy().flatten()
	    
	    peaks = scipy.signal.find_peaks(
	        signal, 
	        height=0.5, 
	        threshold=None, 
	        distance=None, 
	        prominence=None, 
	        width=None,
	        wlen=None, 
	        rel_height=None,
	        plateau_size=None)
	    
	    if len(peaks[0]) > 0:
	    
	        pi = peaks[0] 
	        h = peaks[1]['peak_heights']

	        # Get index of the tallest peak
	        pi_max = (pi[np.argmax(h)]).flatten()

	        peak_widths = scipy.signal.peak_widths(
	            signal, 
	            peaks=pi_max, 
	            rel_height=0.5, 
	            prominence_data=None, 
	            wlen=None)

	        peak_bounds = peak_widths[2:4]

	        peak_bound_min = np.round(peak_bounds[0])
	        peak_bound_max = np.round(peak_bounds[-1])

	        peak_bound_min_pos = df_chr['pos'].iloc[peak_bound_min].item()
	        peak_bound_max_pos = df_chr['pos'].iloc[peak_bound_max].item()

	    # label peaks
	        if len(peak_bounds[0]) != 0:
	            df_chr['peak'] = [1 if x >= peak_bound_min_pos and x <= peak_bound_max_pos else 0 for x in df_chr['pos']]
	    
	    else:
	        
	        df_chr['peak'] = 0
	        
	df = pd.concat(df_list)

	df_peaks = df.loc[df['peak'] == 1]

	#Calculate G-statistic to glean more info about potential causal SNPs
	df['G_S'] = gStatistic_Array(df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'])

	#calculate 98th percentile of G-statistics
	GS_98p = np.percentile(df['G_S'], 98)

	df['G_S_98p'] = [1 if (np.isclose(x, GS_98p) or (x > GS_98p))
	                else 0 for x in df['G_S']]

	dfLC11=df.loc[(df['peak'] == 1) & (df['G_S_98p'] == 1)]
	dfLC01=df.loc[(df['peak'] == 0) & (df['G_S_98p'] == 1)]
	dfLC10=df.loc[(df['peak'] == 1) & (df['G_S_98p'] == 0)]

	df_likely_cands = pd.concat([dfLC11, dfLC01, dfLC10])

	#output final table and candidates

	df_likely_cands_sorted = df_likely_cands.sort_values(by = ['G_S', 'ratio', 'peak'], ascending = [False, False, False], na_position = 'first')
	
	df_likely_cands_sorted.to_csv(finallikelycandidatestablename, sep='\t', index=False)
	df.to_csv(finalcompletetablename, sep='\t',index=False)

	chr_facets_p=df_peaks["chr"].unique()

	######Generate lowess-ratio graph###### (at some point, perhaps considering using matplotlib. Plotly seems convoluted)
	chr_facets_p=df_peaks["chr"].unique()

	#establish fig
	fig = px.scatter(
	    df,
	    x=df['pos'],
	    y=df['ratio'],
	    facet_col="chr",
	    opacity=0.8,
	    color_discrete_sequence=['goldenrod'],
	    trendline="lowess",
	    trendline_options=dict(frac=lowess_span),
	    trendline_color_override="blue",
	)
	            
	#format fig
	fig.update_traces(
	    marker=dict(size=3)
	)

	fig.update_layout(
	    dict(plot_bgcolor = 'white')
	)

	fig.update_xaxes(
	    matches=None
	)

	fig.for_each_xaxis(
	    lambda xaxis: xaxis.update(showticklabels=True)
	)

	fig.for_each_xaxis(
	    lambda xaxis: fig.update_xaxes(tickangle=0)
	)

	fig.for_each_xaxis(
	    lambda x: x.update(title = '')
	)

	fig.for_each_yaxis(
	    lambda y: y.update(title = '')
	)


	fig.add_annotation(
	    showarrow=False,
	    xanchor='center',
	    xref='paper', 
	    x=0.5, 
	    yref='paper',
	    y=-0.1,
	    text='Position'
	)

	fig.add_annotation(
	    showarrow=False,
	    xanchor='center',
	    xref='paper', 
	    x=-0.8, 
	    yanchor='middle',
	    yref='paper',
	    y=0.6,
	    textangle=-90,
	    text='Ratio'
	)


	fig.update_xaxes(
	    showgrid=True, 
	    gridwidth=0.5, 
	    gridcolor='lightgrey'
	)

	fig.update_xaxes(
	    showline=True, 
	    linewidth=1, 
	    linecolor='black'
	)

	fig.update_xaxes(
	    rangemode="tozero"
	)

	fig.update_yaxes(
	    showgrid=True, 
	    gridwidth=0.5, 
	    gridcolor='lightgrey'
	)

	fig.update_yaxes(
	    showline=True, 
	    linewidth=1, 
	    linecolor='black'
	)

	fig.update_yaxes(
	    range=[0, 1]
	)

	fig.update_layout(
	    title=dict(
	        text="Linkage map",font=dict(color='black')
	    )
	)


	fig['data'][0]['showlegend'] = True
	fig['data'][0]['name'] = 'Polymorphism ratio'

	max_pos=df_peaks['pos'].max()
	min_pos=df_peaks['pos'].min()

	        
	fig.add_trace(
	go.Scatter(
	    x=[2],
	    y=[2],
	    mode="markers",
	    name="G-stat > 98 percentile",
	    line=go.scatter.Line(color="red"),
	    marker=dict(size=4),
	    showlegend=True), row=1, col=1
	)


	fig.add_trace(
	    go.Scatter(
	        x=[1],
	        y=[1],
	        mode='lines',
	        name='Lowess Fitted Ratio',
	        line=dict(color="blue")))

	fig.add_trace(
	    go.Scatter(
	        x=[0,0], 
	        y=[0,0], 
	        mode='lines', 
	        line=dict(color='black', width=1, dash='dot'),
	        name='Identified Peak',
	    )
	)

	#add g-stat > 98 percentile markers
	for i in chr_facets:
	    df_chr = df[df['chr']==i].copy()
	    outliers = df_chr.loc[df_chr['G_S_98p'] == 1 , ["gene","pos","ratio"]]
	    fig.add_trace(
	    go.Scatter(
	        x=outliers['pos'],
	        y=outliers['ratio'],
	        mode="markers",
	        name="G-stat > 98 percentile",
	        line=go.scatter.Line(color="red"),
	        marker=dict(size=4),
	        showlegend=False), row=1, col=i
	    )
	    
	#add peak markers
	for i in chr_facets_p:
	    df_peaks_chr=df_peaks[df_peaks['chr']==i]
	    df_chr=df[df['chr'] == i]
	    max_pos=df_peaks_chr['pos'].max()
	    min_pos=df_peaks_chr['pos'].min()

	    if max_pos > 0:
	        max_yhat_pos=df_peaks_chr.loc[df_peaks_chr['yhat'] == df_peaks_chr['yhat'].max()]['pos'].values[0]
	        fig.add_vrect(
	            x0=min_pos, 
	            x1=max_pos, 
	            col=i,
	            fillcolor="red", 
	            opacity=0.2, 
	            line_width=0
	        )
	        
	        fig.add_vline(
	            x=max_yhat_pos, 
	            line_dash="dot", 
	            col=i, 
	            line_width=1
	       )

	pio.write_image(fig, lowess_plotname, scale=6, width=1080, height=540)

	#Generate G-statistic graph
	#establish g-stat lowess fig
	fig = px.scatter(df, x=df['pos'], y=df['G_S'],
	    facet_col="chr",
	    opacity=0.8,
	    color_discrete_sequence=['goldenrod'],
	    trendline="lowess",
	    trendline_options=dict(frac=lowess_span),
	    trendline_color_override="blue")

	#format g-stat lowess fig
	fig.update_traces(
	    marker=dict(size=3)
	)

	fig.update_layout(
	    dict(plot_bgcolor = 'white')
	)

	fig.update_xaxes(
	    matches=None
	)

	fig.for_each_xaxis(
	    lambda xaxis: xaxis.update(showticklabels=True)
	)

	fig.for_each_xaxis(
	    lambda xaxis: fig.update_xaxes(tickangle=0)
	)

	fig.for_each_xaxis(
	    lambda x: x.update(title = '')
	)

	fig.for_each_yaxis(
	    lambda y: y.update(title = '')
	)


	fig.add_annotation(
	    showarrow=False,
	    xanchor='center',
	    xref='paper', 
	    x=0.5, 
	    yref='paper',
	    y=-0.1,
	    text='Position'
	)

	fig.add_annotation(
	    showarrow=False,
	    xanchor='center',
	    xref='paper', 
	    x=-0.8, 
	    yanchor='middle',
	    yref='paper',
	    y=0.6,
	    textangle=-90,
	    text='Ratio'
	)


	fig.update_xaxes(
	    showgrid=True, 
	    gridwidth=0.5, 
	    gridcolor='lightgrey'
	)

	fig.update_xaxes(
	    showline=True, 
	    linewidth=1, 
	    linecolor='black'
	)

	fig.update_xaxes(
	    rangemode="tozero"
	)

	fig.update_yaxes(
	    showgrid=True, 
	    gridwidth=0.5, 
	    gridcolor='lightgrey'
	)

	fig.update_yaxes(
	    showline=True, 
	    linewidth=1, 
	    linecolor='black'
	)

	fig.update_layout(
	    title=dict(
	        text="G-statistic Lowess",font=dict(color='black')
	    )
	)


	fig.add_annotation(
	    showarrow=False,
	    xanchor='center',
	    xref='paper', 
	    x=0.5, 
	    yref='paper',
	    y=-0.1,
	    text='Position'
	)

	fig.add_annotation(
	    showarrow=False,
	    xanchor='center',
	    xref='paper', 
	    x=-0.8, 
	    yanchor='middle',
	    yref='paper',
	    y=0.6,
	    textangle=-90,
	    text='G-statistic'
	)

	fig.add_trace(
	    go.Scatter(
	    x=[1],
	    y=[1],
	    mode='lines',
	    name='Lowess Fitted G-stat',
	    line=dict(color="blue")
	    )
	)


	fig['data'][0]['showlegend'] = True
	fig['data'][0]['name'] = 'G-statistic'

	max_pos=df_peaks['pos'].max()
	min_pos=df_peaks['pos'].min()

	    
	fig.add_trace(go.Scatter(x=[0,0], 
	    y=[0,0], 
	    mode='lines', 
	    line=dict(color='black', width=1, dash='dot'),
	    name='Identified Peak',))

	pio.write_image(fig, gstat_plotname, scale=6, width=1080, height=540)

	print(f'''
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Results for {key} generated. 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
''')