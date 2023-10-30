import sys
import os
import csv
import warnings

import pandas as pd
import numpy as np
import statsmodels.api as sm

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
    theme_linedraw,
    annotate
)

from analysis_functions import *

# Assign current line name to variable
current_line_name = sys.argv[1]

print(f"""
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Producing plots and putative causal mutations for {current_line_name}
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
""")

# Establish variables
# Read in ems table, calculate metrics, format
vcftable = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".noknownsnps.table"
)

df = pd.read_csv(vcftable, sep="\t")

## Calculate delta SNP ratio
df['ratio'] = delta_snp_array(
	df['wt_ref'], 
	df['wt_alt'],
	df['mu_ref'],
	df['mu_alt']
)

## Calculate G-statistic to glean more info about potential causal SNPs
df['G_S'] = g_statistic_array(
	df['wt_ref'], 
	df['wt_alt'], 
	df['mu_ref'], 
	df['mu_alt']
)

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
    # Establish chr copies and smooth edges by inverting duplicated data
    df_chr = df[df['chr']==i].copy()

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

    ## Extend deltas, making "mirrored" delta positions at ends of chr
    deltas_mirrored_ends=[]
    deltas_mirrored_ends.extend(deltas_pos_inv+deltas+deltas_neg_inv)
    
    ## Create new pseudo positions for lowess smoothing 
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

df = pd.concat(df_list)

print("""
Calculating empirical cutoff for LOESS ratio-scaled g-stats 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
""")

# Randomize data 1000 times, calculate metrics to derive cutoffs
## subset position, wt and mu
dfShPos = df[['pos']].copy()
dfShwt = df[['wt_ref', 'wt_alt']].copy()
dfShmu = df[['mu_ref', 'mu_alt']].copy()

## Break the phenotype / read count / position association.
## Iterate and smooth 1000 times 
#Calculate empirical cutoffs for gstatistics, 
gs_cutoff, rsg_cutoff, rsg_y_cutoff = empircal_cutoff(dfShPos, dfShwt, dfShmu)

df['G_S_05p'] = [1 if (np.isclose(x, gs_cutoff) or (x > gs_cutoff))
                else 0 for x in df['G_S']]

df['RS_G_05p'] = [1 if (np.isclose(x, rsg_cutoff) or (x > rsg_cutoff))
                else 0 for x in df['RS_G']]

df['RS_G_yhat_01p'] = [1 if (np.isclose(x, rsg_y_cutoff) or (x > rsg_y_cutoff))
                else 0 for x in df['RS_G_yhat']]

df_likely_cands = df.loc[df['RS_G_yhat_01p'] == 1]

## Use 99.99 percentile cutoff to identify likely candidates
df_likely_cands_sorted = df_likely_cands.sort_values(
    by = ['G_S','RS_G_yhat'], 
    ascending = [False, False], 
    na_position = 'first'
    )

# Save complete table
finalcompletetablename = os.path.join(
	"./output/" 
	+ current_line_name + "/" 
	+ current_line_name 
	+ ".complete_table" 
	+ ".tsv"
)

df.to_csv(finalcompletetablename, sep='\t', index=False)

## save likely candidate list	
finallikelycandidatestablename = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".likely_candidates" 
	+ ".tsv"
)

df_likely_cands.to_csv(
	finallikelycandidatestablename, 
	sep='\t', 
	index=False
)

# Generate Plots
## Establish plot names
gs_plotname = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".gs" 
	+ ".png"
)

gs_yhat_plotname = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".gs_yhat" 
	+ ".png"
)

rs_gs_yhat_plotname = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".rsgs_yhat" 
	+ ".png"
)

delta_snp_plotname = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".delta_snp" 
	+ ".png"
)

delta_snp_yhat_plotname = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".delta_snp_yhat" 
	+ ".png"
)

rs_gs_plotname = os.path.join(
	"./output/" 
	+ current_line_name 
	+ "/" 
	+ current_line_name 
	+ ".rsgs" 
	+ ".png"
)

## Establish aesthics and general save settings
warnings.filterwarnings( "ignore", module = "plotnine\..*" )
df['pos_mb'] = df['pos']*0.000001

points = geom_point(color='goldenrod', size=0.8)
lines = geom_line(color='blue')
themes = theme_linedraw()
facets = facet_grid('. ~ chr', space='free_x', scales='free_x')
spacing = theme(panel_spacing=0.025)

## G-stat plot
chart = ggplot(df, aes('pos_mb', y='G_S'))
title = ggtitle("G-stastic")
axis_x =  xlab("Position (Mb)")
axis_y = ylab("G-statistic")

# Establish plot object
plot = (chart
	+ points
    + themes 
    + facets 
    + title 
    + axis_x
    + axis_y
    + spacing
)

# Save plot
plot.save(
	filename = gs_plotname, 
	height=6, 
	width=8, 
	units = 'in', 
	dpi=500
)

#G-stat yhat plot
chart = ggplot(df, aes('pos_mb', y='G_S_yhat'))
title = ggtitle("Lowess smoothed G-stastic")
axis_x =  xlab("Position (Mb)")
axis_y = ylab("Fitted G-statistic")
cutoff = geom_hline(
	yintercept = gs_cutoff, 
	color = 'red', 
	linetype = "dashed", 
	size = 0.3
)


# Establish plot object
plot = (chart
	+ points
    + themes 
    + facets 
    + title 
    + axis_x
    + axis_y
    + spacing

    + lines
    + cutoff
)

plot.save(
	filename = gs_yhat_plotname, 
	height=6, 
	width=8, 
	units = 'in', 
	dpi=500
)

## Ratio-scaled G-Stat plot
chart = ggplot(df, aes('pos_mb', y='RS_G'))
title = ggtitle("Ratio-scaled G statistic")
axis_x =  xlab("Position (Mb)")
axis_y = ylab("Ratio-scaled G-statistic")
cutoff = geom_hline(
	yintercept = rsg_cutoff, 
	color = 'red', 
	linetype = "dashed", 
	size = 0.3
)

# Establish plot object
plot = (chart
	+ points
    + themes 
    + facets 
    + title 
    + axis_x
    + axis_y
    + spacing

    + cutoff
)

plot.save(
	filename = rs_gs_plotname, 
	height=6, 
	width=8, 
	units = 'in', 
	dpi=500
)

## Delta SNP Ratio plot
chart = ggplot(df, aes('pos_mb', y='ratio'))
title = ggtitle("Delta SNP ratio")
axis_x =  xlab("Position (Mb)")
axis_y = ylab("Ratio")

# Establish plot object
plot = (chart
	+ points
    + themes 
    + facets 
    + title 
    + axis_x
    + axis_y
    + spacing
)

plot.save(
	filename = delta_snp_plotname, 
	height=6, 
	width=8, 
	units = 'in', 
	dpi=500
)


## Fitted delta SNP ratio plot
chart = ggplot(df, aes('pos_mb', y='ratio_yhat'))
title = ggtitle("Fitted Delta SNP ratio")
axis_x =  xlab("Position (Mb)")
axis_y = ylab("Fitted delta SNP ratio")

# Establish plot object
plot = (chart
	+ points
    + themes 
    + facets 
    + title 
    + axis_x
    + axis_y
    + spacing

    + lines
)

plot.save(
	filename = delta_snp_yhat_plotname, 
	height=6, 
	width=8, 
	units = 'in', 
	dpi=500
)


## Fitted Ratio-Scaled G-stat plot
chart = ggplot(df, aes('pos_mb', y='RS_G_yhat'))
title = ggtitle("Lowess smoothed ratio-scaled G statistic")
axis_x =  xlab("Position (Mb)")
axis_y = ylab("Fitted Ratio-scaled G-statistic")
cutoff = geom_hline(
	yintercept = rsg_y_cutoff, 
	color = 'red', 
	linetype = "dashed", 
	size = 0.3
)

# Establish plot object

plot = (chart
	+ points
    + themes 
    + facets 
    + title 
    + axis_x
    + axis_y
    + spacing
    
    + lines
    + cutoff

)

# Save plot
plot.save(
	filename = rs_gs_plotname, 
	height=6, 
	width=8, 
	units = 'in', 
	dpi=500
)


print(f'''
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Results for {current_line_name} generated. 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
''')