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
    analysis_utils = AnalysisUtilities(current_line_name)
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
    except pd.errors.EmptyDataError:
        print(f"Error: File '{vcftable}' is empty.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    # Calculate delta SNP ratio and G-statistic
    print(f" Attempting to calculate delta-SNP ratios and G-statistics.")
    try:
        suppress = False
        df['ratio'] = analysis_utils.delta_snp_array(
            df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'], suppress
        )
        
        df['G_S'] = analysis_utils.g_statistic_array(
            df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'], suppress
        )
        
        print("Calculation of delta-SNP ratios and G-statistics were successful.")
    except Exception as e:
        print(f"An error occurred during calculation: {e}")

    # Clean up data by dropping indels and NAs
    analysis_utils.drop_NA_and_indels(df)

    # LOESS smoothing of ratio and G-stat by chromosome
    try:
        lowess_span = 0.3
        smooth_edges_bounds = 15

        df = analysis_utils.smooth_chr_facets(df, lowess_span, smooth_edges_bounds)
        print("LOESS smoothing calculations successful.")
    except Exception as e:
        print(f"An error occurred during LOESS smoothing: {e}")
    
    print("""
Calculating empirical cutoff for LOESS ratio-scaled g-stats 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
    """)

    # Calculate empirical cutoffs
    try:
        dfShPos = df[['pos']].copy()
        dfShwt = df[['wt_ref', 'wt_alt']].copy()
        dfShmu = df[['mu_ref', 'mu_alt']].copy()
        iterations = 1000
        lowess_span = 0.3

        gs_cutoff, rsg_cutoff, rsg_y_cutoff = analysis_utils.empirical_cutoff(
            dfShPos, dfShwt, dfShmu, iterations, lowess_span
        )

        df['G_S_05p'] = [1 if (np.isclose(x, gs_cutoff) 
            or (x > gs_cutoff)) else 0 for x in df['G_S']
        ]
        df['RS_G_05p'] = [1 if (np.isclose(x, rsg_cutoff) 
            or (x > rsg_cutoff)) else 0 for x in df['RS_G']
        ]
        df['RS_G_yhat_01p'] = [1 if (np.isclose(x, rsg_y_cutoff) 
            or (x > rsg_y_cutoff)) else 0 for x in df['RS_G_yhat']
        ]
      
        df_likely_cands = df.loc[df['RS_G_yhat_01p'] == 1]
        
        print(f"G-statistic cutoff = {gs_cutoff}.")
        print(f"Ratio-scaled G-statistic cutoff = {rsg_cutoff}.")
        print(f"LOESS smoothed Ratio-scaled G-statistic cutoff = {rsg_y_cutoff}.")
        print("Empirical cutoff via randomization calculation completed.")
    except Exception as e:
        print(f"An error occurred during cutoff calculations: {e}")
    
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

    print(f"Attempting to produce and save plots...")
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
        print(f"An error occurred while producing and saving plots: {e}")
