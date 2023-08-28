#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import platform

import csv
import glob
import shlex
import fnmatch
import io

import pandas as pd
import numpy as np

from plotnine import (
    ggplot,
    aes,
    facet_grid,
    geom_point,
    geom_line,
    geom_hline,
    ggtitle,
    scale_x_continuous,
    xlab,
    scale_y_continuous,
    ylab,
    theme,
    element_line,
    element_text,
    theme_linedraw
)

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
# if single_test >= 4 or pair_test >= 8:
# 	multiple_lines = True

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

files_vcf = [os.path.basename(x) for x in glob.glob('./VCFs/*')]


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

# Establish ratio and g-stat functions
def deltaSNParray(wtr, wta, mur, mua):
    #Calculate delta-SNP ratio
    return ((wtr)/(wtr+wta))-((mur)/(mur+mua))

def gStatistic_Array(o1, o3, o2, o4):
    # Calculate G-statistic using numpy array input
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


# Establish variables

for key in lines_dict:

	# Read in ems table, calculate metrics, format
	vcftable = os.path.join("./output/" + key + "/" + key + ".noknownsnps.table")
	df = pd.read_csv(vcftable, sep="\t")

	## Calculate delta SNP ratio
	df['ratio'] = deltaSNParray(df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'])

	## Calculate G-statistic to glean more info about potential causal SNPs
	df['G_S'] = gStatistic_Array(df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'])

	## Drop indels, leave only SNPs
	df[
	    df["ref"].apply(lambda x: len(x) == 1)
	    & df["alt"].apply(lambda x: len(x) == 1)
	]

	## Drop NAs. Infrequently will arise while div/0 seemingly
	df.dropna(axis=0, how='any', subset="ratio", inplace=True)

	# lowess Smoothing of ratio and G-statistic by chromosome
	chr_facets=df["chr"].unique()
	df_list = []

	## Set lowess parameters
	lowess = sm.nonparametric.lowess
	lowess_span=0.3

	for i in chr_facets:
	    # Establish copies of data by chromosome, and smooth edges by inverting duplicated data
	    df_chr = df[df['chr']==i].copy()

	    print(df_chr)

	    ## establish inverted dataframes

	    positions = df_chr['pos'].to_numpy()
	    deltas=[]
	    # Produce deltas
	    for i, pos in enumerate(positions):
	        if i == 0:
	            deltas.append(pos)
	        if i > 0:
	            deltas.append(positions[i] - positions[i-1])
	    
	    ## Create negative and positive delta extensions
	    deltas_pos_inv = deltas[::-1][-15:-1]
	    deltas_neg_inv = deltas[::-1][1:15]

	    ## Extend deltas, making "mirrored" delta positions at both ends of chr
	    deltas_mirrored_ends=[]
	    deltas_mirrored_ends.extend(deltas_pos_inv+deltas+deltas_neg_inv)
	    
	    ## Create new pseudo positions for lowess smoothing by cumulative addition of deltas
	    psuedo_pos = []
	    for i, pos in enumerate(deltas_mirrored_ends):
	        if i == 0:
	            psuedo_pos.append(0)
	        if i > 0:
	            psuedo_pos.append(pos+psuedo_pos[i-1])

	    ## Create "mirrored" data from both ends of chr
	    df_chr_inv_neg = df_chr[::-1].iloc[-15:-1]
	    df_chr_inv_pos = df_chr[::-1].iloc[1:15]

	    ## Concat dataframe
	    df_chr_smooth_list = [df_chr_inv_neg, df_chr, df_chr_inv_pos]

	    df_chr = pd.concat(df_chr_smooth_list, ignore_index=False)
	            
	    ## Add psuedo positions to dataframe
	    df_chr['psuedo_pos'] = psuedo_pos

	    ## Fit LOWESS ratio and G-statistic ~ position
	    X=df_chr['psuedo_pos'].values

	    ## Fit ratio 
	    ratio_Y=df_chr['ratio'].values
	    df_chr['ratio_yhat'] = lowess(ratio_Y,X, frac=lowess_span)[:,1]
	   
	    ## Fit G-Statistic
	    G_S_Y = df_chr['G_S'].values
	    df_chr['G_S_yhat'] = lowess(G_S_Y,X, frac=lowess_span)[:,1]
	    
	    ## Produce ratio-scaled G-statistic and fit
	    df_chr['RS_G'] = df_chr['G_S']*df_chr['ratio']
	    
	    RS_G_Y = df_chr['RS_G'].values
	    df_chr['RS_G_yhat'] = lowess(RS_G_Y,X, frac=lowess_span)[:,1]

	    ## remove edge smoothing (mirrored) data
	    df_chr = df_chr[14:-14]
	    df_chr.drop(axis=1, inplace=True, columns='psuedo_pos')
	    df_list.append(df_chr)
	    print(df_chr)

	df = pd.concat(df_list)

	# Randomize data 1000 times, calculate ratio, GS and fitted values to emperically derive cutoffs
	pos = df[['pos']].copy()
	dfShwt = df[['pos', 'wt_ref', 'wt_alt']].copy()
	dfShmu = df[['mu_ref', 'mu_alt']].copy()

	## Establish lists
	smGstatAll = []
	smRatioAll = []
	RS_GAll = []
	smRS_G_yhatAll = []
	smRatio_yhatAll = []

	## Break the phenotype / read count / position association. iter and smooth 1000 times, collect metrics
	for i in range(1000):
	    dfShPos = pos.sample(frac=1)
	    dfShwt = dfShwt.sample(frac=1)
	    dfShmu = dfShmu.sample(frac=1)

	    smPos = dfShPos['pos'].to_numpy()
	    sm_wt_ref = dfShwt['wt_ref'].to_numpy()
	    sm_wt_alt = dfShwt['wt_alt'].to_numpy()
	    sm_mu_ref = dfShmu['mu_ref'].to_numpy()
	    sm_mu_alt = dfShmu['mu_alt'].to_numpy()
	    
	    #Calculate G-stats per iteration,collect results
	    smGstat = gStatistic_Array(sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt)
	    smGstatAll.extend(smGstat)

	    #Calculate snpRatio per iteration,collect results
	    smRatio = deltaSNParray(sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt)
	    smRatioAll.extend(smRatio)

	    #Calculate ratio-scaled G-stat per iteration,collect results
	    smRS_G = smRatio*smGstat
	    RS_GAll.extend(smRS_G)

	    #Lowess smooth per iteration, collect results
	    smRS_G_yhatAll.extend(lowess(smRS_G,smPos, frac=lowess_span)[:,1])

	# Collect emperical cutoffs from 1000x iter, identify candidates
	G_S_95p = np.percentile(smGstatAll, 95)
	RS_G_95p = np.percentile(RS_GAll, 95)
	RS_G_Y_99p = np.percentile(smRS_G_yhatAll, 99.99)
	ratio_Y_99p = np.percentile(smRatio_yhatAll, 99.99)

	df['G_S_001p'] = [1 if (np.isclose(x, G_S_95p) or (x > G_S_95p))
	                else 0 for x in df['G_S']]

	df['RS_G_001p'] = [1 if (np.isclose(x, G_S_95p) or (x > G_S_95p))
	                else 0 for x in df['G_S']]

	df['RS_G_yhat_001p'] = [1 if (np.isclose(x, G_S_95p) or (x > G_S_95p))
	                else 0 for x in df['G_S']]

	df_likely_cands = df.loc[df['RS_G_yhat_001p'] == 1]

	## Use emperically derived 99.99 percentile cutoff of fitted ratio-scaled g-stat to identify likely candidates
	df_likely_cands_sorted = df_likely_cands.sort_values(
	    by = ['G_S','RS_G_yhat'], 
	    ascending = [False, False], 
	    na_position = 'first'
	    )

	# Save complete table
	finalcompletetablename = os.path.join("./output/" + key + "/" + key + ".complete_table" + ".tsv")
	df.to_csv(finalcompletetablename, sep='\t', index=False)

	## save likely candidate list	
	finallikelycandidatestablename = os.path.join("./output/" + key + "/" + key + ".likely_candidates" + ".tsv")
	df_likely_cands.to_csv(finallikelycandidatestablename, sep='\t', index=False)

	# Generate Plots
	## Establish plot names
	GS_plotname = os.path.join("./output/" + key + "/" + key + ".GS" + ".png")
	GS_yhat_plotname = os.path.join("./output/" + key + "/" + key + ".GS_yhat" + ".png")

	RS_GS_plotname = os.path.join("./output/" + key + "/" + key + ".RS_GS" + ".png")
	RS_GS_yhat_plotname = os.path.join("./output/" + key + "/" + key + ".RS_GS" + ".png")

	deltaSNP = os.path.join("./output/" + key + "/" + key + ".deltaSNP" + ".png")
	deltaSNP_yhat_plotname = os.path.join("./output/" + key + "/" + key + ".deltaSNP_yhat" + ".png")

	## G-stat plot
	chart = ggplot(df, aes('pos_mb', y='G_S'))
	title = ggtitle("G-stastic")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("G-statistic")
	plot = (chart 
	    + points 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + theme(panel_spacing=0.025)
	)
	plot.save(filename = GS_plotname, height=6, width=8, units = 'in', dpi=1000)

	#G-stat yhat plot
	chart = ggplot(df, aes('pos_mb', y='G_S_yhat'))
	title = ggtitle("Lowess smoothed G-stastic")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Fitted G-statistic")
	plot = (chart 
	    + points
	    + lines
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + theme(panel_spacing=0.025)
	)
	plot.save(filename = GS_yhat_plotname, height=6, width=8, units = 'in', dpi=1000)
	
	## Ratio-scaled G-Stat plot
	chart = ggplot(df, aes('pos_mb', y='RS_G'))
	title = ggtitle("Ratio-scaled G statistic")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Ratio-scaled G-statistic")

	plot = (chart 
	    + points 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	)
	plot.save(filename = RS_GS_plotname, height=6, width=8, units = 'in', dpi=1000)

	## Fitted Ratio-Scaled G-stat plot
	df['pos_mb'] = df['pos']*0.000001
	chart = ggplot(df, aes('pos_mb', y='RS_G_yhat'))
	points = geom_point(color='goldenrod', size=0.8)
	lines = geom_line(color='blue')
	themes = theme_linedraw()
	facets = facet_grid('. ~ chr', space='free_x', scales='free_x')
	title = ggtitle("Lowess smoothed ratio-scaled G statistic")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Fitted Ratio-scaled G-statistic")
	spacing = theme(panel_spacing=0.025)
	cutoff = geom_hline(yintercept = RS_G_Y_99p,color='red',linetype="dashed", size=0.3) # add one horizonal line
	plot = (chart 
	    + points 
	    + lines 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	    + cutoff
	)
	plot.save(filename = RS_GS_yhat_plotname, height=6, width=8, units = 'in', dpi=1000)


	## Delta SNP Ratio plot
	chart = ggplot(df, aes('pos_mb', y='ratio'))
	title = ggtitle("Delta SNP ratio")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Ratio")

	plot = (chart 
	    + points 
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	)
	plot.save(filename = deltaSNP, height=6, width=8, units = 'in', dpi=1000)

	#Fitted delta SNP ratio plot
	chart = ggplot(df, aes('pos_mb', y='ratio_yhat'))
	title = ggtitle("Delta SNP ratio")
	axis_x =  xlab("Position (Mb)")
	axis_y = ylab("Ratio")

	plot = (chart 
	    + points
	    + lines
	    + themes 
	    + facets 
	    + title 
	    + axis_x
	    + axis_y
	    + spacing
	)
	plot.save(filename = deltaSNP_yhat_plotname, height=6, width=8, units = 'in', dpi=1000)

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