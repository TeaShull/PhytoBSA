import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from plotnine import ggplot, aes, geom_point, geom_line, theme_linedraw, facet_grid, theme, ggtitle, xlab, ylab, geom_hline
import config

class AnalysisUtilities:
    LOWESS_SPAN = 0.3

    def __init__(self, current_line_name):
        self.current_line_name = current_line_name

    def drop_NA_and_indels(self, df):
        print(f"(Attempt) Removing NAs and indels")
        try:
            df = df[(df["ref"].apply(lambda x: len(x) == 1)) & (df["alt"].apply(lambda x: len(x) == 1))]
            df.dropna(axis=0, how='any', subset=["ratio"], inplace=True)
            print("(Success)Indels dropped, and NaN values cleaned successfully.")
            return df
        except Exception as e:
            print(f"An error occurred during data processing: {e}")
            return None


    def delta_snp_array(self, wtr, wta, mur, mua, suppress):
        """Calculate delta-SNP ratio."""
        try:
            if suppress == False:
                print(f"(Attempt) calculate delta-SNP ratios for {self.current_line_name}...")
            result = ((wtr) / (wtr + wta)) - ((mur) / (mur + mua))
            if suppress == False:
                print(f"(Success) Delta-SNP ratio calculation successful for {self.current_line_name}")
            return result
        except Exception as e:
            print(f"(Fail) Error in delta_snp_array for {self.current_line_name}: {e}")
            return None

    def g_statistic_array(self, o1, o3, o2, o4, suppress):
        """Calculate G-statistic using numpy array input."""
        if suppress == False:
            print(f"(Attempt) calculate G-statistics for {self.current_line_name}....")
        try:
            np.seterr(all='ignore')

            nonzero_mask = o1 + o2 + o3 + o4 != 0

            e1 = np.where(nonzero_mask, (o1 + o2) * (o1 + o3) / (o1 + o2 + o3 + o4), 0)
            e2 = np.where(nonzero_mask, (o1 + o2) * (o2 + o4) / (o1 + o2 + o3 + o4), 0)
            e3 = np.where(nonzero_mask, (o3 + o4) * (o1 + o3) / (o1 + o2 + o3 + o4), 0)
            e4 = np.where(nonzero_mask, (o3 + o4) * (o2 + o4) / (o1 + o2 + o3 + o4), 0)

            llr1 = np.where(o1 / e1 > 0, 2 * o1 * np.log(o1 / e1), 0.0)
            llr2 = np.where(o2 / e2 > 0, 2 * o2 * np.log(o2 / e2), 0.0)
            llr3 = np.where(o3 / e3 > 0, 2 * o3 * np.log(o3 / e3), 0.0)
            llr4 = np.where(o4 / e4 > 0, 2 * o4 * np.log(o4 / e4), 0.0)

            result = np.where(e1 * e2 * e3 * e4 == 0, 0.0, llr1 + llr2 + llr3 + llr4)
            if suppress == False:
                print(f"(Success) G-statistic calculation completed for {self.current_line_name}")
            return result
        except Exception as e:
            print(f"(Fail) Error in g_statistic_array for {self.current_line_name}: {e}")
            return None

    def smooth_chr_facets(self, df, lowess_span, smooth_edges_bounds):
        lowess_function = sm.nonparametric.lowess
        df_list = []

        def smooth_single_chr(df_chr):
            psuedo_pos = []
            positions = df_chr['pos'].to_numpy()
            deltas = [pos - positions[i - 1] if i > 0 else pos for i, pos in enumerate(positions)]
            deltas_pos_inv = deltas[::-1][-smooth_edges_bounds:-1]
            deltas_neg_inv = deltas[::-1][1:smooth_edges_bounds]
            deltas_mirrored_ends = deltas_pos_inv + deltas + deltas_neg_inv
            
            psuedo_pos = []
            for i, pos in enumerate(deltas_mirrored_ends):
                if i == 0:
                    psuedo_pos.append(0)
                if i > 0:
                    psuedo_pos.append(pos+psuedo_pos[i-1])

            df_chr_inv_neg = df_chr[::-1].iloc[-smooth_edges_bounds:-1]
            df_chr_inv_pos = df_chr[::-1].iloc[1:smooth_edges_bounds]
            df_chr_smooth_list = [df_chr_inv_neg, df_chr, df_chr_inv_pos]
            df_chr = pd.concat(df_chr_smooth_list, ignore_index=False)

            df_chr['pseudo_pos'] = psuedo_pos
            X = df_chr['pseudo_pos'].values

            ratio_Y = df_chr['ratio'].values
            df_chr['ratio_yhat'] = lowess_function(ratio_Y, X, frac=lowess_span)[:, 1]

            G_S_Y = df_chr['G_S'].values
            df_chr['GS_yhat'] = lowess_function(G_S_Y, X, frac=lowess_span)[:, 1]

            df_chr['RS_G'] = df_chr['G_S'] * df_chr['ratio']
            RS_G_Y = df_chr['RS_G'].values
            df_chr['RS_G_yhat'] = lowess_function(RS_G_Y, X, frac=lowess_span)[:, 1]

            df_chr = df_chr[smooth_edges_bounds:-smooth_edges_bounds].drop(columns='pseudo_pos')
            return df_chr

        chr_facets = df["chr"].unique()

        for i in chr_facets:
            df_chr = df[df['chr'] == i].copy()
            df_list.append(smooth_single_chr(df_chr))

        return pd.concat(df_list)

    def empirical_cutoff(self, posin, wt, mu, shuffle_iterations, lowess_span):
        """Calculate empirical cutoff."""
        print(f"(Attempt) Calculate empirical cutoff for {self.current_line_name}...")
        try:
            lowess = sm.nonparametric.lowess
            smGstatAll, smRatioAll, RS_GAll, smRS_G_yhatAll = [], [], [], []
            suppress = True
            for _ in range(shuffle_iterations):
                dfShPos = posin.sample(frac=1)
                dfShwt = wt.sample(frac=1)
                dfShmu = mu.sample(frac=1)

                smPos = dfShPos['pos'].to_numpy()
                sm_wt_ref = dfShwt['wt_ref'].to_numpy()
                sm_wt_alt = dfShwt['wt_alt'].to_numpy()
                sm_mu_ref = dfShmu['mu_ref'].to_numpy()
                sm_mu_alt = dfShmu['mu_alt'].to_numpy()

                smGstat = self.g_statistic_array(sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, suppress)
                smGstatAll.extend(smGstat)

                smRatio = self.delta_snp_array(sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, suppress)
                smRatioAll.extend(smRatio)

                smRS_G = smRatio * smGstat
                RS_GAll.extend(smRS_G)

                smRS_G_yhatAll.extend(lowess(smRS_G, smPos, frac=lowess_span)[:, 1])

            G_S_95p = np.percentile(smGstatAll, 95)
            RS_G_95p = np.percentile(RS_GAll, 95)
            RS_G_Y_99p = np.percentile(smRS_G_yhatAll, 99.99)

            result = G_S_95p, RS_G_95p, RS_G_Y_99p
            print(f"(Success) Empirical cutoff calculation completed for {self.current_line_name}")
            return result
        except Exception as e:
            print(f"(Fail) Error in empirical_cutoff for {self.current_line_name}: {e}")
            return None, None, None

    def plot_data(self, df, y_column, title_text, ylab_text, cutoff_value=None, lines=False):
        """Plot data."""
        print(f"(Attempt) Plot data and save plots for {self.current_line_name}...")
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
            output_dir = config.OUTPUT_DIR
            plot_name = f"{self.current_line_name}_{y_column.lower()}.png"
            file_path_name = os.path.join(output_dir, self.current_line_name, plot_name)
            plot.save(filename=file_path_name, 
                height=6, 
                width=8, 
                units = 'in', 
                dpi=500
            )

            print(f"(Success) Plots saved for {self.current_line_name}")
        except Exception as e:
            print(f"(Fail) Plotting data failed for {self.current_line_name}: {e}")