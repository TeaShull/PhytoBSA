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
from config import (
    error_handler, print_delimiter, 
    BASE_DIR, SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR
)
from analysis_utilities import AnalysisUtilities

# Set up directories
src_dir = SRC_DIR
base_dir = BASE_DIR
input_dir = INPUT_DIR
output_dir = OUTPUT_DIR
log_dir = LOG_DIR

# Initialize AnalysisUtilities instance and read in VCF table produced in VCFgen.sh
if __name__ == "__main__":
    try:
        current_line_name = sys.argv[1]
    except IndexError:
        error_handler('fail','Error parsing current_line_name as arg passed from analysis.py')
        print(f"Aborting analysis of {current_line_name}")
        quit()

    print_delimiter(f"Beginning analysis of {current_line_name}. Loading VCF tables.")
    analysis_utils = AnalysisUtilities(current_line_name)
    current_line_out_dir = os.path.join(output_dir, current_line_name)

    # Load VCF table
    vcftable_name = f"{current_line_name}.noknownsnps.table"
    vcftable_path = os.path.join(current_line_out_dir, vcftable_name)
    error_handler('attempt', f"Attempting to load VCF table for line {current_line_name}")
    try:
        df = pd.read_csv(vcftable_path, sep="\t")
        error_handler('attempt', f"The VCF table for line {current_line_name} was successfully loaded.")
    except FileNotFoundError:
        error_handler('fail', f"Error: File '{vcftable_path}' not found.")
    except pd.errors.EmptyDataError:
        error_handler('fail', f"Error: File '{vcftable_path}' is empty.")
    except Exception as e:
        error_handler('fail', f"An unexpected error occurred: {e}")

    # Calculate delta SNP ratio and G-statistic
    print_delimiter(f"Calculating Segregation Statistics for {current_line_name}.")
    error_handler('attempt', f"Initialize calculations for delta-SNP ratios and G-statistics for {current_line_name}")
    try:
        suppress = False
        df['ratio'] = analysis_utils.delta_snp_array(
            df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'], suppress
        )
        
        df['G_S'] = analysis_utils.g_statistic_array(
            df['wt_ref'], df['wt_alt'], df['mu_ref'], df['mu_alt'], suppress
        )
        
        error_handler('success', "Calculation of delta-SNP ratios and G-statistics were successful.")
    except Exception as e:
        error_handler('fail', f"An error occurred during calculation: {e}")

    # Clean up data by dropping indels and NAs
    analysis_utils.drop_NA_and_indels(df)

    # LOESS smoothing of ratio and G-stat by chromosome
    error_handler('attempt', "Initialize LOESS smoothing calculations.")
    try:
        lowess_span = 0.3
        smooth_edges_bounds = 15

        df = analysis_utils.smooth_chr_facets(df, lowess_span, smooth_edges_bounds)
        error_handler('success', "LOESS smoothing calculations successful.")
    except Exception as e:
        error_handler('fail', f"An error occurred during LOESS smoothing: {e}")
    
    print_delimiter('Calculating empirical cutoff for LOESS ratio-scaled g-stats')

    # Calculate empirical cutoffs
    error_handler('attempt', "Initialize calculation of empirical cutoffs")
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
        
        error_handler('success', f"G-statistic cutoff = {gs_cutoff}.")
        error_handler('success', f"Ratio-scaled G-statistic cutoff = {rsg_cutoff}.")
        error_handler('success', f"LOESS smoothed Ratio-scaled G-statistic cutoff = {rsg_y_cutoff}.")
        error_handler('success', f"Empirical cutoff via randomization for {current_line_name} completed.")
    except Exception as e:
        print(f"An error occurred during cutoff calculations: {e}")
    
    error_handler('attempt', 'Initialize the identification of likely candidates')
    try:
        # Identify likely candidates using G-stat and smoothed ratio-scaled G-stat
        likely_cands_sorted = df_likely_cands.sort_values(
            by=['G_S', 'RS_G_yhat'],
            ascending=[False, False],
            na_position='first'
        )

        # Save DataFrames to CSV files
        results_table_name =f"{current_line_name}_results_table.tsv"
        results_table_path = os.path.join(
            current_line_out_dir, results_table_name
        )

        candidates_table_name = f"{current_line_name}_candidates_table.tsv"
        candidates_table_path = os.path.join(
            current_line_out_dir, candidates_table_name
        )

        df.to_csv(results_table_path, sep='\t', index=False)
        likely_cands_sorted.to_csv(candidates_table_path, sep='\t', index=False)
        
        error_handler('success', f"Results and candidates tables for {current_line_name} generated.")
    except Exception as e:
        error_handler('fail', f"An error occurred during {current_line_name} table generation: {e}")
    
    print_delimiter('Generating plots')

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

    error_handler('attempt', f"Attempting to produce and save plots for {current_line_name}...")
    try:
        analysis_utils.current_line_name = current_line_name
        for plot_scenario in plot_scenarios:
            analysis_utils.plot_data(df, *plot_scenario)
        print_delimiter(f"Results for {current_line_name} generated.")
    except Exception as e:
        error_handler('fail', f"An error occurred while producing and saving plots: {e}")
