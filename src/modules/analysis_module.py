import os
import pandas as pd
import numpy as np
from plotnine import (
    ggplot, aes, geom_point, geom_line, theme_linedraw,
    facet_grid, theme, ggtitle, xlab, ylab, geom_hline
)
import statsmodels as sm

from config import (error_handler, print_delimiter, OUTPUT_DIR)

from analysis_utilities import AnalysisUtilities

# All functions in this module are used sequentually by bsa_analysis.ThaleBSAParentFunctions()


def load_vcf_table(current_line_table_path, current_line_name):
    """Load VCF table"""
    error_handler('attempt', 
        f"Attempting to load VCF table for line {current_line_name}"
    )
    try:
        vcf_df = pd.read_csv(current_line_table_path, sep="\t")
        error_handler('attempt', 
            f"The VCF table for line {current_line_name} was successfully loaded."
        )
        return vcf_df
    except FileNotFoundError:
        error_handler('fail', 
            f"Error: File '{current_line_table_path}' not found."
        )
    except pd.errors.EmptyDataError:
        error_handler('fail', 
            f"Error: File '{current_line_table_path}' is empty."
        )
    except Exception as e:
        error_handler('fail', f"An unexpected error occurred: {e}")

def calculate_delta_snp_and_g_statistic(vcf_df, current_line_name):
    analysis_utils = AnalysisUtilities(current_line_name)
    """Calculate delta SNP ratio and G-statistic"""
    error_handler('attempt', f"Initialize calculations for delta-SNP ratios and G-statistics")
    try:
        suppress = False
        vcf_df['ratio'] = analysis_utils.delta_snp_array(
            vcf_df['wt_ref'], vcf_df['wt_alt'], 
            vcf_df['mu_ref'], vcf_df['mu_alt'], 
            suppress
        )
        
        vcf_df['G_S'] = analysis_utils.g_statistic_array(
            vcf_df['wt_ref'], vcf_df['wt_alt'], 
            vcf_df['mu_ref'], vcf_df['mu_alt'], 
            suppress
        )

        vcf_df['RS_G'] = vcf_df['ratio'] * vcf_df['G_S']
        return vcf_df
        error_handler('success', "Calculation of delta-SNP ratios and G-statistics were successful.")
    except Exception as e:
        error_handler('fail', f"An error occurred during calculation: {e}")

def drop_na_and_indels(vcf_df, current_line_name):
    error_handler('attempt', 'Removing NAs and indels')
    try:
        # Use .loc for assignment to avoid the warning
        vcf_df.loc[:, "ref"] = vcf_df["ref"].apply(
            lambda x: x if len(x) == 1 else np.nan
        )
        vcf_df.loc[:, "alt"] = vcf_df["alt"].apply(
            lambda x: x if len(x) == 1 else np.nan
        )
        vcf_df.dropna(axis=0, how='any', subset=["ratio"], inplace=True)
        error_handler('success',
            f'Indels dropped, and NaN values for {current_line_name} cleaned successfully.'
        )
        return vcf_df
    except Exception as e:
        error_handler('fail',
            f"An error occurred during data processing: {e}"
        )
        return None


def calculate_empirical_cutoffs(vcf_df, current_line_name):
    analysis_utils = AnalysisUtilities(current_line_name)
    """Calculate empirical cutoffs"""
    error_handler('attempt', "Initialize calculation of empirical cutoffs")
    try:
        vcf_df_position = vcf_df[['pos']].copy()
        vcf_df_wt = vcf_df[['wt_ref', 'wt_alt']].copy()
        vcf_df_mu = vcf_df[['mu_ref', 'mu_alt']].copy()
        iterations = 1000
        lowess_span = 0.3

        gs_cutoff, rsg_cutoff, rsg_y_cutoff = analysis_utils.empirical_cutoff(
            vcf_df_position, vcf_df_wt, vcf_df_mu, iterations, lowess_span
        )

        vcf_df['G_S_05p'] = [1 if (np.isclose(x, gs_cutoff) 
            or (x > gs_cutoff)) else 0 for x in vcf_df['G_S']
        ]
        vcf_df['RS_G_05p'] = [1 if (np.isclose(x, rsg_cutoff) 
            or (x > rsg_cutoff)) else 0 for x in vcf_df['RS_G']
        ]
        vcf_df['RS_G_yhat_01p'] = [1 if (np.isclose(x, rsg_y_cutoff) 
            or (x > rsg_y_cutoff)) else 0 for x in vcf_df['RS_G_yhat']
        ]
        return vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff

        error_handler('success', f"G-statistic cutoff = {gs_cutoff}.")
        error_handler('success', 
            f"Ratio-scaled G-statistic cutoff = {rsg_cutoff}."
        )
        error_handler('success', 
            f"LOESS smoothed Ratio-scaled G-statistic cutoff = {rsg_y_cutoff}."
        )
        error_handler('success', 
            f"Empirical cutoff via randomization for {current_line_name} completed.")
        
    except Exception as e:
        print(f"An error occurred during cutoff calculations: {e}")
    
def sort_save_likely_candidates(vcf_df, current_line_name):
    """Identify likely candidates"""
    output_dir = OUTPUT_DIR
    error_handler('attempt', 'Initialize the identification of likely candidates')
    try:
        # Identify likely candidates using G-stat and smoothed ratio-scaled G-stat
        vcf_df_likely_cands = vcf_df.loc[vcf_df['RS_G_yhat_01p'] == 1]
        likely_cands_sorted = vcf_df_likely_cands.sort_values(
            by=['G_S', 'RS_G_yhat'],
            ascending=[False, False],
            na_position='first'
        )

        # Save DataFrames to CSV files
        results_table_name =f"{current_line_name}_results_table.tsv"
        results_table_path = os.path.join(
            output_dir, results_table_name
        )
        vcf_df.to_csv(results_table_path, sep='\t', index=False)

        candidates_table_name = f"{current_line_name}_candidates_table.tsv"
        candidates_table_path = os.path.join(
            output_dir, candidates_table_name
        )
        likely_cands_sorted.to_csv(candidates_table_path, sep='\t', index=False)
        
        error_handler('success', 
            f"Results and candidates tables for {current_line_name} generated."
        )
    except Exception as e:
        error_handler('fail', 
            f"An error occurred during {current_line_name} table generation: {e}"
        )
    
def generate_plots(vcf_df, current_line_name, gs_cutoff, rsg_cutoff, rsg_y_cutoff):
    """Generate and save plots. Plot scenarios are below"""
    analysis_utils = AnalysisUtilities(current_line_name)
    # Plot scenarios format:
    # ('y_column', 'title_text', 'ylab_text', cutoff_value=None, lines=False)
    plot_scenarios = [
        ('G_S', 'G-statistic', 'G-sttistic', None, False),
        ('GS_yhat', 'Lowess smoothed G-statistic', 'Fitted G-statistic', 
            gs_cutoff, True
        ),
        ('RS_G', 'Ratio-scaled G statistic', 'Ratio-scaled G-statistic',
            rsg_cutoff, False
         ),
        ('ratio', 'Delta SNP ratio', 'Ratio', None, False),
        ('ratio_yhat', 'Fitted Delta SNP ratio', 'Fitted delta SNP ratio',
            None, True
         ),
        ('RS_G_yhat', 'Lowess smoothed ratio-scaled G statistic', 
            'Fitted Ratio-scaled G-statistic', rsg_y_cutoff, True
        ),
    ]

    error_handler('attempt', 
        f"Attempting to produce and save plots for {current_line_name}..."
    )
    try:
        for plot_scenario in plot_scenarios:
            analysis_utils.plot_data(vcf_df, *plot_scenario)
        print_delimiter(f"Results for {current_line_name} generated.")
    except Exception as e:
        error_handler('fail', 
            f"An error occurred while producing and saving plots: {e}"
        )
