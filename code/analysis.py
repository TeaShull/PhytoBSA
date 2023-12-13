import sys
import os

import pandas as pd
import statsmodels as sm
from plotnine import (
    ggplot, aes, geom_point, geom_line, theme_linedraw,
    facet_grid, theme, ggtitle, xlab, ylab, geom_hline
)

from analysis_functions import *

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

vcftable = os.path.join(
    .,
    output,
    current_line_name,
    f"{current_line_name}.noknownsnps.table"
)

# Read in vcf table generated by VCFgen.sh
try:
    df = pd.read_csv(vcftable, sep="\t")
    print(f"{current_line_name} VCF table successfully loaded.")
except FileNotFoundError:
    print(f"Error: File '{vcftable}' not found.")
    
except pd.errors.EmptyDataError:
    print(f"Error: File '{vcftable}' is empty.")

except Exception as e:
    print(f"An unexpected error occurred: {e}")

# Calculate statistics
try:
    # Calculate delta SNP ratio
    df['ratio'] = delta_snp_array(
        df['wt_ref'], 
        df['wt_alt'],
        df['mu_ref'],
        df['mu_alt']
    )

    # Calculate G-statistic
    df['G_S'] = g_statistic_array(
        df['wt_ref'], 
        df['wt_alt'], 
        df['mu_ref'], 
        df['mu_alt']
    )

    print("Calculation of delta SNP ratio and G-statistic successful.")

except Exception as e:
    print(f"An error occurred during calculation: {e}")

# Clean up data by dropping indels and NAs
try:
    df = df[
        df["ref"].apply(lambda x: len(x) == 1)
        & df["alt"].apply(lambda x: len(x) == 1)
    ]

    df.dropna(axis=0, how='any', subset=["ratio"], inplace=True)
    print("Indels dropped, and NaN values cleaned successfully.")

except Exception as e:
    print(f"An error occurred during data processing: {e}")

"""
lowess Smoothing of ratio and G-stat by chromosome. This helps extract a 
signal from the noise
"""
try:
    chr_facets = df["chr"].unique()
    df_list = []

    lowess = sm.nonparametric.lowess
    lowess_span = 0.3

    for i in chr_facets:
        # Correct for Loess edge bias - extend chr edges via inversion of 15 rows.
        df_chr = df[df['chr'] == i].copy()

        positions = df_chr['pos'].to_numpy()
        deltas = []

        for i, pos in enumerate(positions):
            if i == 0:
                deltas.append(pos)
            if i > 0:
                deltas.append(positions[i] - positions[i-1])

        deltas_pos_inv = deltas[::-1][-15:-1]
        deltas_neg_inv = deltas[::-1][1:15]

        deltas_mirrored_ends = []
        deltas_mirrored_ends.extend(deltas_pos_inv + deltas + deltas_neg_inv)

        pseudo_pos = [0] + [
            pos + pseudo_pos[i-1]
            for i, pos in enumerate(deltas_mirrored_ends) if i > 0
        ]

        df_chr_inv_neg = df_chr[::-1].iloc[-15:-1]
        df_chr_inv_pos = df_chr[::-1].iloc[1:15]

        df_chr_smooth_list = [df_chr_inv_neg, df_chr, df_chr_inv_pos]

        df_chr = pd.concat(df_chr_smooth_list, ignore_index=False)

        df_chr['pseudo_pos'] = pseudo_pos

        # Fit LOWESS ratio and G-statistic ~ position
        x = df_chr['pseudo_pos'].values

        ratio_y = df_chr['ratio'].values
        df_chr['ratio_yhat'] = lowess(ratio_y, x, frac=lowess_span)[:, 1]

        # Fit G-Statistic
        gs_y = df_chr['G_S'].values
        df_chr['GS_yhat'] = lowess(gs_y, x, frac=lowess_span)[:, 1]

        # Produce ratio-scaled G-statistic and fit
        df_chr['RS_G'] = df_chr['G_S']*df_chr['ratio']

        ratio_scaled_g_stat = df_chr['RS_G'].values
        df_chr['RS_G_yhat'] = lowess(ratio_scaled_g_stat, x, frac=lowess_span)[:, 1]

        # remove edge smoothing data
        df_chr = df_chr[14:-14]
        df_chr.drop(columns='pseudo_pos', inplace=True)
        df_list.append(df_chr)

    df = pd.concat(df_list)
    print("LOESS smoothing calculations successful.")

except Exception as e:
    print(f"An error occurred during LOESS smoothing: {e}")

print("""
Calculating empirical cutoff for LOESS ratio-scaled g-stats 
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<
""")

    """
    Breaking the phenotype / read count / position association.
    Iterate and loess smooth 1000 times to Calculate empirical cutoffs for 
    gstatistics, ratio_scaled g-statistics and smoothed ratio_scaled g-statistics. 
    """
try:
    # subset position, wt, and mu
    dfShPos = df[['pos']].copy()
    dfShwt = df[['wt_ref', 'wt_alt']].copy()
    dfShmu = df[['mu_ref', 'mu_alt']].copy()

    gs_cutoff, rsg_cutoff, rsg_y_cutoff = empircal_cutoff(dfShPos, dfShwt, dfShmu)

	df['G_S_05p'] = [
	    1 if (np.isclose(x, gs_cutoff) or (x > gs_cutoff)) else 0
	    for x in df['G_S']
	]

	df['RS_G_05p'] = [
	    1 if (np.isclose(x, rsg_cutoff) or (x > rsg_cutoff)) else 0
	    for x in df['RS_G']
	]

	df['RS_G_yhat_01p'] = [
	    1 if (np.isclose(x, rsg_y_cutoff) or (x > rsg_y_cutoff)) else 0
	    for x in df['RS_G_yhat']
	]

	df_likely_cands = df.loc[df['RS_G_yhat_01p'] == 1]

    print("Randomization and metric calculation completed.")

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
        os.getcwd(), 
        "output",
        current_line_name,
        f"{current_line_name}_results_table.tsv"
    )

    candidates_table_name = os.path.join(
        os.getcwd(), 
        "output",
        current_line_name,
        f"{current_line_name}_candidates_table.tsv"
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

""" 
List of plots. Articulate the plot scenarios using the following format:
('y_column', 'title_text', 'ylab_text', cutoff_value=None, lines=False)
"""
plot_scenarios = [
    ('G_S', 'G-statistic', 'G-statistic', None, False),
    ('GS_yhat', 'Lowess smoothed G-statistic', 'Fitted G-statistic', gs_cutoff, True),
    ('RS_G', 'Ratio-scaled G statistic', 'Ratio-scaled G-statistic', rsg_cutoff, False),
    ('ratio', 'Delta SNP ratio', 'Ratio', None, False),
    ('ratio_yhat', 'Fitted Delta SNP ratio', 'Fitted delta SNP ratio', None, True),
    ('RS_G_yhat', 'Lowess smoothed ratio-scaled G statistic', 'Fitted Ratio-scaled G-statistic', rsg_y_cutoff, True),
]

try:
    for plot_scenario in plot_scenarios:
        plot_data(df, *plot_scenario)

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