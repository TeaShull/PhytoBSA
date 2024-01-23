from settings.config import (BASE_DIR, INPUT_DIR, MODULES_DIR, OUTPUT_DIR)

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
from plotnine import (
    ggplot, aes, geom_point, geom_line, theme_linedraw,
    facet_grid, theme, ggtitle, xlab, ylab, geom_hline
)
from modules.utilities_general import FileUtilities

"""
Module for the bsa_analysis parent function. Save for a few file utilities used 
in the parent function, the core of the read-depth analysis between the 
wild-type and mutant bulks are stored here.
"""

class BSAAnalysisUtilities:
    def __init__(self, current_line_name, vcf_df, vcf_ulid, logger):
        self.current_line_name = current_line_name
        self.log = logger
        self.vcf_ulid = vcf_ulid
        self.analysis_out_prefix = f'{self.log.ulid}_-{current_line_name}_analysis'

        if self.vcf_ulid:
            self.analysis_out_path = os.path.join(
                OUTPUT_DIR, 
                f'{self.vcf_ulid}_-{current_line_name}',
                self.analysis_out_prefix
            )
        else:
            self.analysis_out_path = os.path.join(
                OUTPUT_DIR,
                self.analysis_out_prefix
            )        
        #Make the analysis_out_path directory if it does not exist
        file_utils=FileUtilities(self.log)
        file_utils.setup_directory(self.analysis_out_path)

        #Initialize loess function parameters
        self.loess_function = sm.nonparametric.lowess
        self.loess_span = 0.3
        self.smooth_edges_bounds = 15
        self.shuffle_iterations = 1000

    def drop_indels(self):
        """
        Drops insertion/deletions from VCF dataframe.
        
        Args: 
        vcf_df(pd.DataFrame)
        VCF dataframe
        
        Returns: 
        VCF dataframe with no indels
        """
        self.vcf_df.loc[~(self.vcf_df["ref"].str.len() > 1) 
            & ~(self.vcf_df["alt"].str.len() > 1)]

    def drop_na(self):
        """
        Drops rows with NaN values from VCF dataframe.
        
        Args: 
        vcf_df(pd.DataFrame)
        VCF dataframe
        
        Produces: 
        VCF dataframe with no NaN values
        """
        self.vcf_df.dropna(axis=0, how='any', subset=["ratio"])

    def filter_genotypes(self, segregation_type: str):
        """
        Filter genotypes in the 'mu:wt_GTpred' column of a DataFrame based on 
        the specified allele.
        Args:
            allele (str): The allele value to filter the genotypes. Filters: 
                          R = Recessive seg '1/1:0/1', '0/1:0/0', '0/1:0/1'.
                          D = Dominant seg '0/1:0/0', '1/1:0/0', '0/1:0/1'.
            df (pd.DataFrame): The input DataFrame containing the genotypes.

        Returns:
            pandas.DataFrame: A filtered DataFrame containing only the rows with 
            matching genotypes.

            [IMPORTANT NOTE]
            #0/1:0/1 is included because of occasianal leaky genotypying by GATK 
            haplotype caller. Nearly 100% of negative delta SNP values arise from 
            0/1:0/1 situations. To retain information without losing data that
            may help fit GAM or LOESS, we will retain 0/1:0/1 loci and instead 
            cut the negative values out after calculating the delta allele
            see: self.drop_genos_with_negative_ratios
        """
        self.log.attempt('Attempting to filter genotypes based on segrigation pattern...')
        try:
            if segregation_type == 'R':
                self.log.note('Filtering genotypes based an a recessive segrigation pattern')
                seg_filter = ['1/1:0/1', '0/1:0/0', '0/1:0/1']
            elif segregation_type == 'D':
                self.log.note('Filtering genotypes based an a dominant segrigation pattern')
                seg_filter = ['0/1:0/0', '1/1:0/0', '0/1:0/1']  
            else: 
                self.log.fail(f'Allele type:{segregation_type} is not a valid selection! Aborting.')
            
            self.vcf_df[self.vcf_df['mu:wt_GTpred'].isin(seg_filter)]
            self.log.success('Genotypes filtured based on segrigation pattern')
        
        except Exception as e:
            self.log.fail(f'There was an error while filtering genotypes:{e}')        

    def filter_ems_mutations(self):
        """
        Filter mutations likely to be from EMS for analysis and return a filtered DataFrame.

        Args:
            df (pd.DataFrame): The input DataFrame containing the mutations.

        Returns:
            pd.DataFrame: A filtered DataFrame containing the mutations.
        """
        self.vcf_df[(self.vcf_df['ref'].isin(['G', 'C', 'A', 'T'])) & (self.vcf_df['alt'].isin(['A', 'T', 'G', 'C']))]

    def drop_genos_with_negative_ratios(self):
        '''
        Removes those genotypes that give rise to negative delta SNP ratios.
        These genotypes are nearly always the result of 0/1:0/1 situations, 
        which exist only because GATK haplotype caller sometimes doesn't genotype 
        everything perfectly and data that is useful for fitting LOESS or GAM 
        models get cleaned out if 0/1:0/1 genotypes aren't included. 
        '''
        self.log.attempt('Trying to remove Genotypes that produce negative delta SNP ratios')
        try: 
            self.vcf_df[self.vcf_df['ratio'] >= 0]
            self.log.success('Genotypes that produce negative delta SNP ratios removed.')
        except Exception as e:
            self.log.fail(f'There was an error removing genotypes that produce nagative delta snp ratios:{e}')
    
    def calculate_delta_snp_and_g_statistic(self):
        """
        Calculate delta SNP ratio and G-statistic
        
        Args: 
        vcf_df 
        pd.DataFrame VCF with no NA and indels.
        
        Returns: 
        pd.DataFrame VCF with delta-snp and g-stat calculated
        """
        
        self.log.attempt(f"Initialize calculations for delta-SNP ratios and G-statistics")
        try:
            suppress = False
            self.vcf_df['ratio'] = self._delta_snp_array(
                self.vcf_df['wt_ref'], self.vcf_df['wt_alt'], 
                self.vcf_df['mu_ref'], self.vcf_df['mu_alt'],
                suppress
            )
            self.vcf_df['G_S'] = self._g_statistic_array(
                self.vcf_df['wt_ref'], self.vcf_df['wt_alt'], 
                self.vcf_df['mu_ref'], self.vcf_df['mu_alt'],
                suppress
            )
            self.vcf_df['RS_G'] = self.vcf_df['ratio'] * self.vcf_df['G_S']
            self.log.success("Calculation of delta-SNP ratios and G-statistics were successful.")
            return vcf_df
        except Exception as e:
            self.log.fail( f"An error occurred during calculation: {e}")

    def _delta_snp_array(self, wtr: np.ndarray, wta: np.ndarray, mur: np.ndarray, mua: np.ndarray, suppress: bool)-> np.ndarray:
        """
        Calculates delta SNP feature, which quantifies divergence in 
            read depths between the two bulks.

        Args: 
            wtr, wta, mur, mua (numpy array)
            Read depths of wt reference and alt reads
            and the read depths of mutant reference and alt reads 
            
            suppress(boolian)
            so that the calculate_empirical_cutoffs function does not print 
            excessively into the log during the 1000x iteration of shuffling.

        Returns: Delta-snp calculation, which is a quantification of 
            allelic segregation at each polymorphic site.
        """ 

        if not suppress:
            self.log.attempt(f"Calculate delta-SNP ratios for {self.current_line_name}...")
        try:
            result = ((wtr) / (wtr + wta)) - ((mur) / (mur + mua))
            if not suppress:
                self.log.success(f"Delta-SNP calculation successful for {self.current_line_name}")
            return result
        
        except Exception as e:
            self.log.fail(f"Error in delta_snp_array for {self.current_line_name}: {e}")
            return None

    def _g_statistic_array(self, wtr, wta, mur, mua, suppress)->np.ndarray:
        """
        Calculates g-statistic feature, which is a more statistically driven 
        approach to calculating read-depth divergence from expected values. 
        Chi square ish.
        """ 
        if not suppress:
            self.log.attempt(f"Calculate G-statistics for {self.current_line_name}....")
        try:
            np.seterr(all='ignore')

            zero_mask = wtr + mur + wta + mua != 0
            denominator = wtr + mur + wta + mua

            e1 = np.where(zero_mask, (wtr + mur) * (wtr + wta) / (denominator), 0)
            e2 = np.where(zero_mask, (wtr + mur) * (mur + mua) / (denominator), 0)
            e3 = np.where(zero_mask, (wta + mua) * (wtr + wta) / (denominator), 0)
            e4 = np.where(zero_mask, (wta + mua) * (mur + mua) / (denominator), 0)

            llr1 = np.where(wtr / e1 > 0, 2 * wtr * np.log(wtr / e1), 0.0)
            llr2 = np.where(mur / e2 > 0, 2 * mur * np.log(mur / e2), 0.0)
            llr3 = np.where(wta / e3 > 0, 2 * wta * np.log(wta / e3), 0.0)
            llr4 = np.where(mua / e4 > 0, 2 * mua * np.log(mua / e4), 0.0)

            result = np.where(
                e1 * e2 * e3 * e4 == 0, 0.0, llr1 + llr2 + llr3 + llr4
            )
            if not suppress:
                self.log.success(f"G-statistic calculation complete for {self.current_line_name}"
                )
            return result

        except Exception as e:
            self.log.fail(f"Error in g_statistic_array for {self.current_line_name}: {e}")
            return None

    def loess_smoothing(self)->pd.DataFrame:
        """
        LOESS smoothing of ratio and G-stat by chromosome
        
        Input: Cleaned dataframe with delta SNPs and G-stats calculated
        
        Returns: Dataframe containing LOESS fitted values for ratio, g-stat and 
        ratio-scaled g-stat
        """
        self.log.attempt("Initialize LOESS smoothing calculations.")
        self.log.attempt(f"span: {self.loess_span}, Edge bias correction: {self.smooth_edges_bounds}") 
        try:

            vcf_df = self._smooth_chr_facets(vcf_df)
            self.log.success("LOESS smoothing calculations successful.")
            return vcf_df
        
        except Exception as e:
            self.log.fail( f"An error occurred during LOESS smoothing: {e}")
    
    def _smooth_chr_facets(self, df)->pd.DataFrame:
        """
        Internal Function for smoothing chromosome facets using LOESS. 
        Uses function "smooth_single_chr" to interate over chromosomes as facets
        to generate LOESS smoothed values for g-statistics and delta-SNP feature
        
        Input: vcf_df
        
        output: vcf_df updated with gs, ratio, and gs-ratio yhat values
        """        
        df_list = []
        chr_facets = df["chr"].unique()

        def _smooth_single_chr(df_chr, chr):
            """
            Input: df_chr chunk, extends the data 15 data values in each
            direction (to mitigate LOESS edge bias), fits smoothed values and 
            subsequently removes the extended data. 
            Returns: df with fitted values included.
            """
            self.log.attempt(f"LOESS of chr:{chr} for {self.current_line_name}...")

            positions = df_chr['pos'].to_numpy()
            
            deltas = ([pos - positions[i - 1]
                       if i > 0 else pos for i, pos in enumerate(positions)]
                      )
            deltas_pos_inv = deltas[::-1][-self.smooth_edges_bounds:-1]
            deltas_neg_inv = deltas[::-1][1:self.smooth_edges_bounds]
            deltas_mirrored_ends = deltas_pos_inv + deltas + deltas_neg_inv

            psuedo_pos = []
            for i, pos in enumerate(deltas_mirrored_ends):
                if i == 0:
                    psuedo_pos.append(0)
                if i > 0:
                    psuedo_pos.append(pos + psuedo_pos[i - 1])

            df_chr_inv_neg = df_chr[::-1].iloc[-self.smooth_edges_bounds:-1]
            df_chr_inv_pos = df_chr[::-1].iloc[1:self.smooth_edges_bounds]
            df_chr_smooth_list = [df_chr_inv_neg, df_chr, df_chr_inv_pos]
            df_chr = pd.concat(df_chr_smooth_list, ignore_index=False)

            df_chr['pseudo_pos'] = psuedo_pos
            X = df_chr['pseudo_pos'].values

            ratio_Y = df_chr['ratio'].values
            df_chr['ratio_yhat'] = (
                self.loess_function(ratio_Y, X, frac=self.loess_span)[:, 1]
            )

            G_S_Y = df_chr['G_S'].values
            df_chr['GS_yhat'] = (
                self.loess_function(G_S_Y, X, frac=self.loess_span)[:, 1]
            )

            df_chr['RS_G'] = df_chr['G_S'] * df_chr['ratio']
            RS_G_Y = df_chr['RS_G'].values
            df_chr['RS_G_yhat'] = (
                self.loess_function(RS_G_Y, X, frac=self.loess_span)[:, 1]
            )

            df_chr = df_chr[self.smooth_edges_bounds:-self.smooth_edges_bounds].drop(
                columns='pseudo_pos'
            )

            self.log.success(f"LOESS of chr:{chr} for {self.current_line_name} complete"
            )
            return df_chr

        self.log.attempt('Smoothing chromosome facets')
        try:
            for i in chr_facets:
                df_chr = df[df['chr'] == i]
                result = _smooth_single_chr(df_chr, i)
                if result is not None:
                    df_list.append(result)
            self.log.success('Chromosome facets LOESS smoothed')
            self.vcf_df =  pd.concat(df_list)
        except Exception as e:
            self.log.fail(f'There was an error during LOESS smoothing of chromosome facets:{e}')

    def calculate_empirical_cutoffs(self):
        """
        Calculate empirical cutoffs.
        
        Input: processed VCF dataframe.  
        
        produces:self.gs_cutoff, self.rsg_cutoff, self.rsg_y_cutoff 
        """
        self.log.attempt("Initialize calculation of empirical cutoffs")
        self.log.note(f"breaking geno/pheno association. iterations:{self.shuffle_iterations}, LOESS span:{self.loess_span}")
        
        try:
            vcf_df_position = self.vcf_df[['pos']].copy()
            vcf_df_wt = self.vcf_df[['wt_ref', 'wt_alt']].copy()
            vcf_df_mu = self.vcf_df[['mu_ref', 'mu_alt']].copy()

            self.gs_cutoff, self.rsg_cutoff, self.rsg_y_cutoff = self._empirical_cutoff(
                vcf_df_position, vcf_df_wt, vcf_df_mu 
            )

            self.vcf_df['G_S_05p'] = [1 if (np.isclose(x, self.gs_cutoff) 
                or (x > self.gs_cutoff)) else 0 for x in self.vcf_df['G_S']
            ]
            self.vcf_df['RS_G_05p'] = [1 if (np.isclose(x, self.rsg_cutoff) 
                or (x > self.rsg_cutoff)) else 0 for x in self.vcf_df['RS_G']
            ]
            self.vcf_df['RS_G_yhat_01p'] = [1 if (np.isclose(x, self.rsg_y_cutoff) 
                or (x > self.rsg_y_cutoff)) else 0 for x in self.vcf_df['RS_G_yhat']
            ]

            self.log.success(f"G-statistic cutoff = {self.gs_cutoff}.")
            self.log.success(f"Ratio-scaled G-statistic cutoff = {self.rsg_cutoff}.")
            self.log.success(f"LOESS smoothed Ratio-scaled G-statistic cutoff = {self.rsg_y_cutoff}.")
            self.log.success(f"Empirical cutoff via randomization for {self.current_line_name} completed.")

        except Exception as e:
            self.log.fail(f"An error occurred during cutoff calculations: {e}")

    def _empirical_cutoff(self, vcf_df_position, vcf_df_wt, vcf_df_mu):
        """
        randomizes the input read_depths, breaking the position/feature link.
        this allows the generation of a large dataset which has no linkage 
        information, establishing an empirical distribution of potenial 
        delta-snps and g-statistics. Given the data provided, we can then set 
        reasonable cutoff values for the fitted values
        
        There is probably a less computationally intensive statistical framework 
        for doing this, especially for the g-statistics....
        
        input: various arrays pulled from vcf_df. 
            vcf_df_position(array) - genome positions
            vcf_df_wt(array) - wt read depth
            vcf_df_mu(array) - mu read depth

            shuffle_iterations(int) - how many iterations to perform. 
                [More = More consistancy, more compute 
                Less = less consistancy, less compute]
        """
        self.log.attempt(f"Calculate empirical cutoff for {self.current_line_name}...")
        
        try:
            lowess = sm.nonparametric.lowess
            smGstatAll, smRatioAll, RS_GAll, smRS_G_yhatAll = [], [], [], []
            suppress = True
            for _ in range(self.shuffle_iterations):
                dfShPos = vcf_df_position.sample(frac=self.loess_span)
                dfShwt = vcf_df_wt.sample(frac=self.loess_span)
                dfShmu = vcf_df_mu.sample(frac=self.loess_span)

                smPos = dfShPos['pos'].to_numpy()
                sm_wt_ref = dfShwt['wt_ref'].to_numpy()
                sm_wt_alt = dfShwt['wt_alt'].to_numpy()
                sm_mu_ref = dfShmu['mu_ref'].to_numpy()
                sm_mu_alt = dfShmu['mu_alt'].to_numpy()

                smGstat = self._g_statistic_array(
                    sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, suppress
                )
                smGstatAll.extend(smGstat)

                smRatio = self._delta_snp_array(
                    sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, suppress
                )
                smRatioAll.extend(smRatio)
                
                smRS_G = smRatio * smGstat
                RS_GAll.extend(smRS_G)

                smRS_G_yhatAll.extend(lowess(
                    smRS_G, smPos, frac=self.loess_span)[:, 1]
                )
            G_S_95p = np.percentile(smGstatAll, 95)
            RS_G_95p = np.percentile(RS_GAll, 95)
            RS_G_Y_99p = np.percentile(smRS_G_yhatAll, 99.99)
            result = G_S_95p, RS_G_95p, RS_G_Y_99p

            self.log.success(f"Empirical cutoff calculation completed for {self.current_line_name}")
            return result
        
        except Exception as e:
            self.log.fail(f"Error in empirical_cutoff for {self.current_line_name}: {e}")
            return None, None, None        

    def sort_save_likely_candidates(self):
        """
        Identify likely candidates
        """
        self.log.attempt('Initialize the identification of likely candidates')
        self.log.note(f'associated VCF table ulid: {self.vcf_ulid}')

        try:
            # Identify likely candidates using G-stat and smoothed ratio-scaled G-stat
            vcf_df_likely_cands = vcf_df[
            (vcf_df['RS_G_yhat_01p'] == 1) |
            (vcf_df['G_S_05p'] == 1) |
            (vcf_df['RS_G_05p'] == 1)
            ]
            #sort
            likely_cands_sorted = vcf_df_likely_cands.sort_values(
                by=['G_S', 'RS_G_yhat'],
                ascending=[False, False],
                na_position='first'
            )
            # Save DataFrames to CSV files
            results_table_name =f"{self.analysis_out_prefix}_results_table.tsv"
            results_table_path = os.path.join(
                self.analysis_out_path, results_table_name
            )
            vcf_df.to_csv(results_table_path, sep='\t', index=False)

            candidates_table_name = f"{self.analysis_out_prefix}_candidates_table.tsv"
            candidates_table_path = os.path.join(
                self.analysis_out_path, candidates_table_name
            )
            likely_cands_sorted.to_csv(candidates_table_path, sep='\t', index=False)
            
            self.log.success(f"Results and candidates tables for {self.current_line_name} generated.")

        except Exception as e:
            self.log.fail(f"An error occurred during {self.current_line_name} table generation: {e}")

    def generate_plots(self):
        """
        Generate and save plots. Plot scenarios are below
        Plot scenarios format:
        ('y_column', 'title_text', 'ylab_text', cutoff_value=None, lines=False)
        """
        plot_scenarios = [
            ('G_S', 'G-statistic', 'G-statistic', self.gs_cutoff, False),
            ('GS_yhat', 'Lowess smoothed G-statistic', 'Fitted G-statistic', None, True),
            ('RS_G', 'Ratio-scaled G statistic', 'Ratio-scaled G-statistic', self.rsg_cutoff, False),
            ('ratio', 'Delta SNP ratio', 'Ratio', None, False),
            ('ratio_yhat', 'Fitted Delta SNP ratio', 'Fitted delta SNP ratio', None, True),
            ('RS_G_yhat', 'Lowess smoothed ratio-scaled G statistic', 'Fitted Ratio-scaled G-statistic', self.rsg_y_cutoff, True),
        ]

        self.log.attempt(f"Attempting to produce and save plots for {self.current_line_name}...")
        
        try:
            for plot_scenario in plot_scenarios:
                self._plot_data(*plot_scenario)
            
            self.log.delimiter(f"Results for {self.current_line_name} generated.")

        except Exception as e:
            self.log.fail(f"An error occurred while producing and saving plots: {e}")
        
    def _plot_data(self, y_column, title_text, ylab_text, cutoff_value=None, lines=False):
        """
        Generate and save plots.
        """
        warnings.filterwarnings("ignore", module="plotnine\..*")

        self.log.attempt(f"Plot data and save plots for {self.current_line_name}...")
        try:
            mb_conversion_constant = 0.000001
            self.self.vcf_df['pos_mb'] = self.self.vcf_df['pos'] * mb_conversion_constant
            chart = ggplot(self.vcf_df, aes('pos_mb', y=y_column))
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
            plot_name = f"{self.analysis_out_prefix}_{y_column.lower()}.png"
            file_path_name = os.path.join(self.analysis_out_path, plot_name)
            plot.save(
                filename=file_path_name,
                height=6,
                width=8,
                units='in',
                dpi=500
            )
            self.log.success(f"Plot saved {file_path_name}")
        
        except Exception as e:
            self.log.fail(f"Plotting data failed for {self.current_line_name}: {e}")

