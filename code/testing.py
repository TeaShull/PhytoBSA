import sys
import os
import warnings
import pandas as pd
import statsmodels.api as sm
from plotnine import (
    ggplot, aes, geom_point, geom_line, theme_linedraw,
    facet_grid, theme, ggtitle, xlab, ylab, geom_hline
)
from analysis_functions import *

# Assign current line name to variable
current_line_name = sys.argv[1]

# ...

# Read in ems table, calculate metrics, format
vcftable = os.path.join(
    "./output/",
    current_line_name,
    current_line_name + ".noknownsnps.table"
)

df = pd.read_csv(vcftable, sep="\t")

# ...

# Drop indels, leave only SNPs
df = df[
    (df["ref"].apply(lambda x: len(x) == 1))
    & (df["alt"].apply(lambda x: len(x) == 1))
]

# Drop NAs. Infrequently will arise while div/0 seemingly
df.dropna(axis=0, how='any', subset=["ratio"], inplace=True)

# ...

# lowess Smoothing of ratio and G-statistic by chromosome
chr_facets = df["chr"].unique()
df_list = []

# ...

# Calculate empirical cutoffs for gstatistics
gs_cutoff, rsg_cutoff, rsg_y_cutoff = empircal_cutoff(dfShPos, dfShwt, dfShmu)

# ...

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

# ...

# Generate Plots
warnings.filterwarnings("ignore", module="plotnine\..*")

# Define a function for plotting
def plot_data(df, y_column, title_text, ylab_text, cutoff_value=None, lines=False):
    chart = ggplot(df, aes('pos_mb', y=y_column))
    title = ggtitle(title_text)
    axis_x = xlab("Position (Mb)")
    axis_y = ylab(ylab_text)

    plot = (chart + geom_point(color='goldenrod', size=0.8) +
            theme_linedraw() + facet_grid('. ~ chr', space='free_x', scales='free_x') +
            title + axis_x + axis_y + theme(panel_spacing=0.025))

    if cutoff_value is not None:
        cutoff = geom_hline(yintercept=cutoff_value, color='red', linetype="dashed", size=0.3)
        plot += cutoff

    if lines:
        plot += geom_line(color='blue')

    # Save plot
    plot.save(filename=f"./output/{y_column.lower()}_plotname.png", height=6, width=8, units='in', dpi=500)

# List of scenarios
scenarios = [
    ('G_S', 'G-statistic', 'G-statistic', None, False),
    ('GS_yhat', 'Lowess smoothed G-statistic', 'Fitted G-statistic', gs_cutoff, True),
    ('RS_G', 'Ratio-scaled G statistic', 'Ratio-scaled G-statistic', rsg_cutoff, False),
    ('ratio', 'Delta SNP ratio', 'Ratio', None, False),
    ('ratio_yhat', 'Fitted Delta SNP ratio', 'Fitted delta SNP ratio', None, True),
    ('RS_G_yhat', 'Lowess smoothed ratio-scaled G statistic', 'Fitted Ratio-scaled G-statistic', rsg_y_cutoff, True),
]

for scenario in scenarios:
    plot_data(df, *scenario)

# ...