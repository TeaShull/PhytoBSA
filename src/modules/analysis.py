#!/usr/bin/env python

import sys
import os
import pandas as pd
import numpy as np
from plotnine import (
    ggplot, aes, geom_point, geom_line, theme_linedraw,
    facet_grid, theme, ggtitle, xlab, ylab, geom_hline
)
import statsmodels as sm
import config 
from analysis_utilities import AnalysisUtilities

# Set up directories
src_dir = config.SRC_DIR
base_dir = config.BASE_DIR
input_dir = config.INPUT_DIR
output_dir = config.OUTPUT_DIR
log_dir = config.LOG_DIR

# Initialize AnalysisUtilities instance
if __name__ == "__main__":
    try:
        current_line_name = sys.argv[1]
    except IndexError:
        print("Error parsing current_line_name as arg passed from analysis.py")
        print(f"Aborting analysis of {current_line_name}")
        quit()

    print({current_line_name})
    analysis_utils = AnalysisUtilities()
    current_line_out_dir = os.path.join(output_dir, current_line_name)

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

    # Load VCF table
    vcftable = os.path.join(current_line_out_dir, f"{current_line_name}.noknownsnps.table")

    try:
        df = pd.read_csv(vcftable, sep="\t")
        print(f"The VCF table for line {current_line_name} was successfully loaded.")
    except FileNotFoundError:
        print(f"Error: File '{vcftable}' not found.")
        analysis_utils.abort_analysis(current_line_name)
    except pd.errors.EmptyDataError:
        print(f"Error: File '{vcftable}' is empty.")
        analysis_utils.abort_analysis(current_line_name)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        analysis_utils.abort_analysis(current_line_name)

    # Calculate delta SNP ratio and G-statistic
    try:
        df['ratio'] = analysis_utils.delta_snp_array(df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'])
        df['G_S'] = analysis_utils.g_statistic_array(df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'])
        print("Calculation of delta-SNP ratios and G-statistics were successful.")
    except Exception as e:
        print(f"An error occurred during calculation: {e}")
        analysis_utils.abort_analysis(current_line_name)

    # Clean up data by dropping indels and NAs
    try:
        df = df[(df["ref"].apply(lambda x: len(x) == 1)) & (df["alt"].apply(lambda x: len(x) == 1))]
        df.dropna(axis=0, how='any', subset=["ratio"], inplace=True)
        print("Indels dropped, and NaN values cleaned successfully.")
    except Exception as e:
        print(f"An error occurred during data processing: {e}")
        analysis_utils.abort_analysis(current_line_name)

    # LOESS smoothing of ratio and G-stat by chromosome
    try:
        chr_facets = df["chr"].unique()
        df_list = []
        lowess = sm.api.nonparametric.lowess
        lowess_span = 0.3

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
            df_chr['GS_yhat'] = lowess(G_S_Y,X, frac=lowess_span)[:,1]
            
            ## Produce ratio-scaled G-statistic and fit
            df_chr['RS_G'] = df_chr['G_S']*df_chr['ratio']
            
            RS_G_Y = df_chr['RS_G'].values
            df_chr['RS_G_yhat'] = lowess(RS_G_Y,X, frac=lowess_span)[:,1]

            ## remove edge smoothing (mirrored) data
            df_chr = df_chr[14:-14]
            df_chr.drop(axis=1, inplace=True, columns='psuedo_pos')
            df_list.append(df_chr)

        df = pd.concat(df_list)

        print("LOESS smoothing calculations successful.")
    except Exception as e:
        print(f"An error occurred during LOESS smoothing: {e}")
        analysis_utils.abort_analysis(current_line_name)
    print("""
Calculating empirical cutoff for LOESS ratio-scaled g-stats 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
    """)

    # Calculate empirical cutoffs
    try:
        dfShPos = df[['pos']].copy()
        dfShwt = df[['wt_ref', 'wt_alt']].copy()
        dfShmu = df[['mu_ref', 'mu_alt']].copy()
        gs_cutoff, rsg_cutoff, rsg_y_cutoff = analysis_utils.empirical_cutoff(dfShPos, dfShwt, dfShmu)

        df['G_S_05p'] = [1 if (np.isclose(x, gs_cutoff) or (x > gs_cutoff)) else 0 for x in df['G_S']]
        df['RS_G_05p'] = [1 if (np.isclose(x, rsg_cutoff) or (x > rsg_cutoff)) else 0 for x in df['RS_G']]
        df['RS_G_yhat_01p'] = [1 if (np.isclose(x, rsg_y_cutoff) or (x > rsg_y_cutoff)) else 0 for x in df['RS_G_yhat']]
        df_likely_cands = df.loc[df['RS_G_yhat_01p'] == 1]
        
        print(f"G-statistic cutoff = {gs_cutoff}.")
        print(f"Ratio-scaled G-statistic cutoff = {rsg_cutoff}.")
        print(f"LOESS smoothed Ratio-scaled G-statistic cutoff = {rsg_y_cutoff}.")
        print("Empirical cutoff via randomization calculation completed.")
    except Exception as e:
        print(f"An error occurred during cutoff calculations: {e}")
        analysis_utils.abort_analysis(current_line_name)
    
    try:
        # Identify likely candidates using G-stat and smoothed ratio-scaled G-stat
        likely_cands_sorted = df_likely_cands.sort_values(
            by=['G_S', 'RS_G_yhat'],
            ascending=[False, False],
            na_position='first'
        )

        # Save DataFrames to CSV files
        results_table_name = os.path.join(
            current_line_out_dir, f"{current_line_name}_results_table.tsv"
        )

        candidates_table_name = os.path.join(
            current_line_out_dir, f"{current_line_name}_candidates_table.tsv"
        )

        df.to_csv(results_table_name, sep='\t', index=False)
        likely_cands_sorted.to_csv(candidates_table_name, sep='\t', index=False)

        print("Results and candidates tables generated.")
    except Exception as e:
        print(f"An error occurred during table generation: {e}")
        analysis_utils.abort_analysis(current_line_name)
    print("""
Generating plots 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
    """)

    # List of plots. Articulate the plot scenarios using the following format:
    # ('y_column', 'title_text', 'ylab_text', cutoff_value=None, lines=False)
    plot_scenarios = [
        ('G_S', 'G-statistic', 'G-statistic', None, False),
        ('GS_yhat', 'Lowess smoothed G-statistic', 'Fitted G-statistic', gs_cutoff, True),
        ('RS_G', 'Ratio-scaled G statistic', 'Ratio-scaled G-statistic', rsg_cutoff, False),
        ('ratio', 'Delta SNP ratio', 'Ratio', None, False),
        ('ratio_yhat', 'Fitted Delta SNP ratio', 'Fitted delta SNP ratio', None, True),
        ('RS_G_yhat', 'Lowess smoothed ratio-scaled G statistic', 'Fitted Ratio-scaled G-statistic', rsg_y_cutoff, True),
    ]

    try:
        analysis_utils.current_line_name = current_line_name
        for plot_scenario in plot_scenarios:
            analysis_utils.plot_data(df, *plot_scenario)

        print(f'''
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
Results for {current_line_name} generated. 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
        ''')
    except Exception as e:
        print(f"An error occurred during plot generation: {e}")
        analysis_utils.abort_analysis(current_line_name)