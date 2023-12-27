import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
from plotnine import (
    ggplot, aes, geom_point, geom_line, theme_linedraw,
    facet_grid, theme, ggtitle, xlab, ylab, geom_hline
)
from config import (
    error_handler, BASE_DIR, SRC_DIR,
    INPUT_DIR, MODULES_DIR, OUTPUT_DIR
)

class AnalysisUtilities:
    LOWESS_SPAN = 0.3

    def __init__(self, current_line_name):
        self.current_line_name = current_line_name

    def drop_NA_and_indels(self, df):
        error_handler('attempt', 'Removing NAs and indels')
        try:
            # Use .loc for assignment to avoid the warning
            df.loc[:, "ref"] = df["ref"].apply(
                lambda x: x if len(x) == 1 else np.nan
            )
            df.loc[:, "alt"] = df["alt"].apply(
                lambda x: x if len(x) == 1 else np.nan
            )
            df.dropna(axis=0, how='any', subset=["ratio"], inplace=True)
            error_handler('success',
                'Indels dropped, and NaN values cleaned successfully.'
            )
            return df
        except Exception as e:
            error_handler('fail',
                f"An error occurred during data processing: {e}"
            )
            return None

    def delta_snp_array(self, wtr, wta, mur, mua, suppress):
        if not suppress:
            error_handler('attempt',
                f"Calculate delta-SNP ratios for {self.current_line_name}..."
            )
        try:
            result = ((wtr) / (wtr + wta)) - ((mur) / (mur + mua))
            if not suppress:
                error_handler('success',
                    f"Delta-SNP calculation successful for {self.current_line_name}"
                )
            return result
        except Exception as e:
            error_handler('fail',
                f"Error in delta_snp_array for {self.current_line_name}: {e}"
            )
            return None

    def g_statistic_array(self, o1, o3, o2, o4, suppress):
        if not suppress:
            error_handler('attempt',
                f"Calculate G-statistics for {self.current_line_name}...."
            )
        try:
            np.seterr(all='ignore')

            zero_mask = o1 + o2 + o3 + o4 != 0
            denominator = o1 + o2 + o3 + o4

            e1 = np.where(zero_mask, (o1 + o2) * (o1 + o3) / (denominator), 0)
            e2 = np.where(zero_mask, (o1 + o2) * (o2 + o4) / (denominator), 0)
            e3 = np.where(zero_mask, (o3 + o4) * (o1 + o3) / (denominator), 0)
            e4 = np.where(zero_mask, (o3 + o4) * (o2 + o4) / (denominator), 0)

            llr1 = np.where(o1 / e1 > 0, 2 * o1 * np.log(o1 / e1), 0.0)
            llr2 = np.where(o2 / e2 > 0, 2 * o2 * np.log(o2 / e2), 0.0)
            llr3 = np.where(o3 / e3 > 0, 2 * o3 * np.log(o3 / e3), 0.0)
            llr4 = np.where(o4 / e4 > 0, 2 * o4 * np.log(o4 / e4), 0.0)

            result = np.where(
                e1 * e2 * e3 * e4 == 0, 0.0, llr1 + llr2 + llr3 + llr4
            )

            if not suppress:
                error_handler('success',
                    f"G-statistic calculation complete for {self.current_line_name}"
                )
            return result
        except Exception as e:
            error_handler('fail',
                f"Error in g_statistic_array for {self.current_line_name}: {e}"
            )
            return None

    def smooth_chr_facets(self, df, lowess_span, smooth_edges_bounds):
        error_handler('attempt', 'Smoothing chromosome facets')
        lowess_function = sm.nonparametric.lowess
        df_list = []

        def smooth_single_chr(df_chr, chr):
            error_handler('attempt',
                f"LOESS of chr:{chr} for {self.current_line_name}..."
            )
            psuedo_pos = []
            positions = df_chr['pos'].to_numpy()
            deltas = ([pos - positions[i - 1]
                       if i > 0 else pos for i, pos in enumerate(positions)]
                      )
            deltas_pos_inv = deltas[::-1][-smooth_edges_bounds:-1]
            deltas_neg_inv = deltas[::-1][1:smooth_edges_bounds]
            deltas_mirrored_ends = deltas_pos_inv + deltas + deltas_neg_inv

            psuedo_pos = []
            for i, pos in enumerate(deltas_mirrored_ends):
                if i == 0:
                    psuedo_pos.append(0)
                if i > 0:
                    psuedo_pos.append(pos + psuedo_pos[i - 1])

            df_chr_inv_neg = df_chr[::-1].iloc[-smooth_edges_bounds:-1]
            df_chr_inv_pos = df_chr[::-1].iloc[1:smooth_edges_bounds]
            df_chr_smooth_list = [df_chr_inv_neg, df_chr, df_chr_inv_pos]
            df_chr = pd.concat(df_chr_smooth_list, ignore_index=False)

            df_chr['pseudo_pos'] = psuedo_pos
            X = df_chr['pseudo_pos'].values

            ratio_Y = df_chr['ratio'].values
            df_chr['ratio_yhat'] = (
                lowess_function(ratio_Y, X, frac=lowess_span)[:, 1]
            )

            G_S_Y = df_chr['G_S'].values
            df_chr['GS_yhat'] = (
                lowess_function(G_S_Y, X, frac=lowess_span)[:, 1]
            )

            df_chr['RS_G'] = df_chr['G_S'] * df_chr['ratio']
            RS_G_Y = df_chr['RS_G'].values
            df_chr['RS_G_yhat'] = (
                lowess_function(RS_G_Y, X, frac=lowess_span)[:, 1]
            )

            df_chr = df_chr[smooth_edges_bounds:-smooth_edges_bounds].drop(
                columns='pseudo_pos'
            )

            error_handler('success',
                f"LOESS of chr:{chr} for {self.current_line_name} complete"
            )
            return df_chr

        chr_facets = df["chr"].unique()

        for i in chr_facets:
            df_chr = df[df['chr'] == i]
            result = smooth_single_chr(df_chr, i)
            if result is not None:
                df_list.append(result)

        return pd.concat(df_list)

    def empirical_cutoff(self, vcf_df_position, vcf_df_wt,
                         vcf_df_mu, shuffle_iterations, lowess_span):
        error_handler('attempt',
            f"Calculate empirical cutoff for {self.current_line_name}..."
        )
        try:
            lowess = sm.nonparametric.lowess
            smGstatAll, smRatioAll, RS_GAll, smRS_G_yhatAll = [], [], [], []
            suppress = True
            for _ in range(shuffle_iterations):
                dfShPos = vcf_df_position.sample(frac=1)
                dfShwt = vcf_df_wt.sample(frac=1)
                dfShmu = vcf_df_mu.sample(frac=1)

                smPos = dfShPos['pos'].to_numpy()
                sm_wt_ref = dfShwt['wt_ref'].to_numpy()
                sm_wt_alt = dfShwt['wt_alt'].to_numpy()
                sm_mu_ref = dfShmu['mu_ref'].to_numpy()
                sm_mu_alt = dfShmu['mu_alt'].to_numpy()

                smGstat = self.g_statistic_array(
                    sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, suppress
                )
                smGstatAll.extend(smGstat)

                smRatio = self.delta_snp_array(
                    sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, suppress
                )
                smRatioAll.extend(smRatio)

                smRS_G = smRatio * smGstat
                RS_GAll.extend(smRS_G)

                smRS_G_yhatAll.extend(lowess(
                    smRS_G, smPos, frac=lowess_span)[:, 1]
                )

            G_S_95p = np.percentile(smGstatAll, 95)
            RS_G_95p = np.percentile(RS_GAll, 95)
            RS_G_Y_99p = np.percentile(smRS_G_yhatAll, 99.99)

            result = G_S_95p, RS_G_95p, RS_G_Y_99p
            error_handler('success',
                f"Empirical cutoff calculation completed for {self.current_line_name}"
            )
            return result
        except Exception as e:
            error_handler('fail',
                f"Error in empirical_cutoff for {self.current_line_name}: {e}"
            )
            return None, None, None

    def plot_data(self, df, y_column, title_text,
                  ylab_text, cutoff_value=None, lines=False
                  ):
        warnings.filterwarnings("ignore", module="plotnine\..*")
        error_handler('attempt',
            f"Plot data and save plots for {self.current_line_name}..."
        )

        try:
            mb_conversion_constant = 0.000001
            df['pos_mb'] = df['pos'] * mb_conversion_constant
            chart = ggplot(df, aes('pos_mb', y=y_column))
            title = ggtitle(title_text)
            axis_x = xlab("Position (Mb)")
            axis_y = ylab(ylab_text)

            if cutoff_value is not None:
                cutoff = geom_hline(yintercept=cutoff_value, color='red',
                                    linetype="dashed", size=0.3
                                    )
                plot = (chart
                        + geom_point(color='goldenrod', size=0.8)
                        + theme_linedraw()
                        + facet_grid('. ~ chr', space='free_x', scales='free_x')
                        + title
                        + axis_x
                        + axis_y
                        + theme(panel_spacing=0.025)
                        + cutoff
                        )
            else:
                plot = (chart
                        + geom_point(color='goldenrod', size=0.8)
                        + theme_linedraw()
                        + facet_grid('. ~ chr', space='free_x', scales='free_x')
                        + title
                        + axis_x
                        + axis_y
                        + theme(panel_spacing=0.025)
                        )

            if lines:
                plot += geom_line(color='blue')

            # Save plot
            output_dir = OUTPUT_DIR 
            plot_name = f"{self.current_line_name}_{y_column.lower()}.png"
            file_path_name = os.path.join(
                output_dir, self.current_line_name, plot_name
            )
            plot.save(
                filename=file_path_name,
                height=6,
                width=8,
                units='in',
                dpi=500
            )

            error_handler('success', f"Plot saved {file_path_name}")
        except Exception as e:
            error_handler('fail',
                f"Plotting data failed for {self.current_line_name}: {e}"
            )



# Assuming OUTPUT_DIR is defined somewhere
# OUTPUT_DIR = ...

# Example usage
# Replace with your actual data and parameters
# analysis = AnalysisUtilities("YourCurrentLineName")
# df = pd.DataFrame(...)  # Your DataFrame
# analysis.drop_NA_and_indels(df)
# analysis.plot_data(df, y_column='your_column', title_text='Your Title', ylab_text='Your Y-axis Label')

