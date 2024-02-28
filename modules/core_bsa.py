import os
import bisect

import numpy as np
import pandas as pd
import statsmodels.api as sm

from plotnine import *
from itertools import groupby
from multiprocessing import Pool

from modules.utilities_logging import LogHandler
from modules.utilities_general import FileUtilities
from settings.globals import THREADS_LIMIT

"""
Core module for bsa analysis
Input variable class: BSA_variables from utilities_variables module.
The read-depth analysis between the wild-type and mutant bulks are stored here.

Most logging and error handling is done outside of the BSA class. 
"""

class BSA:
    def __init__(self, logger, bsa_vars):
        self.log = logger 
        self.bsa_vars = bsa_vars

    def _initialize_log(self, line):
        self.log.delimiter(f'Initializing BSA pipeline log for {line.name}')
        bsa_log = LogHandler(f'analysis_{line.name}')
        line.analysis_ulid = bsa_log.ulid 
        
        if not line.vcf_ulid:#store associated VCF ulid for records and file pathing
            file_utils = FileUtilities(self.log)
            line.vcf_ulid = file_utils.extract_ulid_from_file_path(line.vcf_table_path)
        
        bsa_log.add_db_record(
            name = line.name, 
            core_ulid = self.log.ulid, 
            vcf_ulid = line.vcf_ulid, 
            ratio_cutoff = self.bsa_vars.ratio_cutoff, 
            loess_span = self.bsa_vars.loess_span, 
            smooth_edges_bounds = self.bsa_vars.smooth_edges_bounds, 
            filter_indels = self.bsa_vars.filter_indels, 
            filter_ems = self.bsa_vars.filter_ems, 
            snpmask_path = self.bsa_vars.snpmask_path,
            segregation_type = line.segregation_type
        )

        return bsa_log

    def _load_data(self, line, bsa_log):

        line.vcf_df = self.bsa_vars.load_vcf_table(line.vcf_table_path)

    def _filter_data(self, line, bsa_log):
        #Filter data based on boolean values set in Check config.ini or argparse
        data_filter = DataFiltering(bsa_log, line.name)
        line.vcf_df = data_filter.filter_genotypes(line.segregation_type, line.vcf_df)
        
        if self.bsa_vars.filter_indels and line.segregation_type != 'QTL': 
            line.vcf_df = data_filter.drop_indels(line.vcf_df)
        
        if self.bsa_vars.filter_ems and line.segregation_type != 'QTL': 
            line.vcf_df = data_filter.filter_ems_mutations(line.vcf_df)
        
        if self.bsa_vars.snpmask_path and self.bsa_vars.mask_snps: 
            line.snpmask_df = self.bsa_vars.load_snpmask(
                self.bsa_vars.snpmask_path, self.bsa_vars.snpmask_url)
            line.vcf_df = data_filter.mask_known_snps(line.snpmask_df, line.vcf_df)

    def _produce_features(self, line, bsa_log):
        # Produce SNP ratio and G-statistic
        feature_prod = FeatureProduction(bsa_log, line.name) #Init FeatureProduction class
        line.vcf_df = feature_prod.calculate_delta_snp_and_g_statistic(line.vcf_df)
        
        #some additional filtering for better results...
        data_filter = DataFiltering(bsa_log, line.name)
        line.vcf_df = data_filter.drop_na(line.vcf_df)
        line.vcf_df = data_filter.drop_genos_below_ratio_cutoff(line.vcf_df, self.bsa_vars.ratio_cutoff)
        
        line.vcf_df.reset_index(drop=True)
        
        #Produced fitted G and fitted ratio to improve signal/noise
        line.vcf_df = feature_prod.fit_model(line.vcf_df, sm.nonparametric.lowess, 
            self.bsa_vars.loess_span, self.bsa_vars.smooth_edges_bounds
        )

    def _produce_null_models(self, line, bsa_log):
        #Init Feature Production Class
        feature_prod = FeatureProduction(bsa_log, line.name)
        
        #Generate H0* distrubtions. Locus-wise for fitted, global otherwise 
        null_models = feature_prod.calculate_null_models(
            line.vcf_df, sm.nonparametric.lowess, self.bsa_vars.loess_span, 
            self.bsa_vars.shuffle_iterations, self.bsa_vars.ratio_cutoff
        )
        #Aggregate values from unsmoothed features.
        null_models = feature_prod.aggregate_unsmoothed_values(null_models)
        
        #Label the dataframe with percentiles from the null models. 
        feature_prod.label_df_with_percentiles(line.vcf_df, null_models)
        
        #Trim out extra data from edge bias adjustment
        line.vcf_df, null_models = feature_prod.remove_extra_data(line.vcf_df, 
            null_models
        )

        return null_models

    def _save_and_plot_outputs(self, line, null_models, bsa_log):
        #Generate output prefix
        line.analysis_out_prefix = self.bsa_vars.gen_bsa_out_prefix(line.name, 
            line.analysis_ulid, line.vcf_ulid
        )
        
        #init TableAndPlots class
        table_and_plots = TableAndPlots(
            bsa_log,
            line.name,
            line.analysis_out_prefix
        )

        #Plot Null models 
        table_and_plots.plot_null_histograms(null_models)
        
        #Identify likely candidates based on H0* (No phenotype/genotype link)
        table_and_plots.process_and_save_candidates(line.vcf_df)
        table_and_plots.generate_plots(line.vcf_df)
        
        #Save run parameters for posterity
        file_utils = FileUtilities(self.log)
        out_vars = f"{line.analysis_out_prefix}.bsa_vars.txt"
        class_instances = [self.bsa_vars, line]
        file_utils.write_instance_vars_to_file(class_instances, out_vars)

    def __call__(self):
        #Run BSA pipeline
        for line in self.bsa_vars.lines:
            try:
                bsa_log = self._initialize_log(line)
                self._load_data(line, bsa_log)
                self._filter_data(line, bsa_log)
                self._produce_features(line, bsa_log)
                
                null_models = self._produce_null_models(line, bsa_log)
                self._save_and_plot_outputs(line, null_models, bsa_log)
            except Exception as e:
                self.log.error(f"Running BSA for line {line.name} failed unexpectedly. {e}")
                self.log.error("Continuing to attempt analysis of other lines...")
                continue

class DataFiltering:
    def __init__ (self, logger, name):
        self.log = logger
        self.name = name
    
    def drop_indels(self, vcf_df: pd.DataFrame)-> pd.DataFrame:
        """
        Drops insertion/deletions from VCF dataframe.
        
        Args: 
        vcf_df(pd.DataFrame)
        VCF dataframe
        
        Returns: 
        VCF dataframe with no indels
        """
        self.log.attempt('Attempting to drop indels...')
        try:
            self.log.note(f'Input dataframe length: {len(vcf_df)}')
            vcf_df = vcf_df.loc[~(vcf_df["ref"].str.len() > 1) 
                & ~(vcf_df["alt"].str.len() > 1)
            ]
            self.log.note(f'Filtered dataframe length: {len(vcf_df)}')
            
            self.log.success("Indels dropped")
            
            return vcf_df
        
        except AttributeError:
            self.log.error("'ref' and 'alt' columns should only contain strings. VCF may not be properly formatted. Aborting...")
        
        except KeyError:
            self.log.error("'ref' or 'alt' column not found in the DataFrame. Please ensure they exist.")
    
    def drop_na(self, vcf_df: pd.DataFrame)-> pd.DataFrame:
        """
        Drops rows with NaN values from VCF dataframe.
        
        Args: 
            vcf_df(pd.DataFrame)
            VCF dataframe
        
        Produces: 
            VCF dataframe with no NaN values
        """

        self.log.attempt('Attempting to drop NaN values...')
        try:
            self.log.note(f'Input dataframe length: {len(vcf_df)}')
            vcf_df = vcf_df.dropna(axis=0, how='any', subset=["ratio"])
            self.log.note(f'Filtered dataframe length: {len(vcf_df)}')
        
            self.log.success('NaN values dropped')
        
            return vcf_df
        
        except Exception as e:
            self.log.error(f'There was an error dropping NaN values{e}')     

    def filter_genotypes(self, segregation_type: str, vcf_df: pd.DataFrame)-> pd.DataFrame:
        """
        Filter genotypes in the 'mu:wt_GTpred' column of a DataFrame based on 
        the specified allele. 1/1:0/0 is included because in some situations, 
        read depths
        
        Args:
            segregation_type (str): The allele value to filter the genotypes. 
                Filters: 
                        R = Recessive seg '1/1:0/1', '0/1:0/1', 1/1:0/0.
                        D = Dominant seg '0/1:0/0', '0/1:0/1', 1/1:0/0.
            
            vcf_df (pd.DataFrame): The input DataFrame containing the genotypes.

        Returns:
            pd.DataFrame: Filtered DataFrame containing only the rows with 
            matching genotypes.

        [EXTRA INFO]
        Nearly 100% of negative delta SNP values arise from 
        0/1:0/1 situations. To retain information without losing data that
        may help fit GAM or LOESS, we will retain 0/1:0/1 loci and instead 
        cut the negative values out after calculating the delta allele
        see: self.drop_genos_with_negative_ratios
        
        """
        self.log.attempt('Attempting to filter genotypes based on segregation pattern...')
        try:
            if segregation_type == 'R':
                self.log.note('Filtering genotypes based an a recessive segregation pattern')
                seg_filter = ['1/1:0/1', '0/1:0/1']
            
            elif segregation_type == 'D':
                self.log.note('Filtering genotypes based an a dominant segregation pattern')
                seg_filter = ['0/1:0/0', '0/1:0/1','1/1:0/0']  
            
            elif segregation_type =='QTL':
                self.log.note('Segregation type is QTL. Skipping filtering...')
                
                return vcf_df

            else: 
                self.log.error(f'Segregation type value is not valid or not assigned. Continuing with broadest filter...')
                self.log.warning("Data should be mostly OK, but some spurious mutations won't be filtered. Null model generation will be more resource intensive")
                seg_filter = ['0/1:0/0', '0/1:0/1', '1/1:0/1', '1/1:0/0'] 
            
            try:
                self.log.note(f'Input dataframe length: {len(vcf_df)}')
                vcf_df = vcf_df[vcf_df['mu:wt_GTpred'].isin(seg_filter)]
                self.log.note(f'Filtered dataframe length: {len(vcf_df)}')
                
                self.log.success('Genotypes filtured based on segregation pattern')
                
                return vcf_df

            except KeyError as e:
                self.log.note('Key error. VCF dataframe should have the following headers: ')
                self.log.note('chrom pos ref alt gene snpEffect snpVariant snpImpact mu:wt_GTpred mu_ref mu_alt wt_ref wt_alt')
                self.log.error(f"Dataframe doesn't contain {e} column. Aborting...")
            
        
        except Exception as e:
            self.log.error(f'There was an error while filtering genotypes:{e}')        

    def filter_ems_mutations(self, vcf_df: pd.DataFrame)-> pd.DataFrame:
        """
        Filter mutations likely to be from EMS for analysis and return a filtered DataFrame.

        Args:
            vcf_df (pd.DataFrame): The input DataFrame containing the mutations.

        Returns:
            pd.DataFrame: A filtered DataFrame containing the mutations.
        """

        self.log.attempt('Filtering varients to only include those likely to arise through EMS exposure...')
        ems_snps = [('G', 'A'), ('C', 'T'), ('A', 'G'), ('T', 'C')]
        self.log.note(f'Input dataframe length: {len(vcf_df)}')
        
        
        vcf_df = vcf_df[vcf_df[['ref', 'alt']].apply(tuple, axis=1).isin(ems_snps)]
        
        self.log.note(f'Filtered dataframe length: {len(vcf_df)}')
        self.log.success('Varients filtered.')
        
        return vcf_df

    def drop_genos_below_ratio_cutoff(self, vcf_df: pd.DataFrame, ratio_cutoff)-> pd.DataFrame:
        '''
        Removes those genotypes that give rise to negative delta SNP ratios.
        args:
            vcf_df: pd.DataFrame - VCF dataframe
        
        Returns:
            vcf_df: pd.DataFrame - Filtered dataframe with no negative delta-snp
            values
        '''

        self.log.attempt(f'Removing Genotypes that produce delta SNP ratios below cutoff:{ratio_cutoff}')
        try: 
            self.log.note(f'Input dataframe length: {len(vcf_df)}')

            vcf_df = vcf_df[vcf_df['ratio'] >= ratio_cutoff]
            
            self.log.success(f'Genotypes that produce delta SNP ratios below cutoff ({ratio_cutoff}) removed.')
            self.log.note(f'Filtered dataframe length: {len(vcf_df)}')

            return vcf_df

        except Exception as e:
            self.log.error(f'There was an error removing genotypes that produce nagative delta snp ratios:{e}')

    def mask_known_snps(self, snpmask_df: pd.DataFrame, vcf_df: pd.DataFrame) -> pd.DataFrame:
        '''
        This fuction applys the SNP mask. If the user has provided a collection
        of background snps in a suitable format, (headers aren't case sensitive
        but they must exists between the two VCFs) those background snps will
        be removed before analysis proceeds. This leads a much cleaner output, 
        and is particularly useful if aligning EMS mutants in an background that
        diverges from the reference genome, as there will be many spurious 
        background variants that are mostly irrelevant to your analysis. 

        args: 
            snpmask_df: pd.DataFrame - df containing background snps
            vcf_df: pd.DataFrame - VCF dataframe
        returns:
            vcf_df: pd.DataFrame - Filtered dataframe without background snps
        '''
        self.log.attempt("Filtering known snps from provided snpmask file..")
        try:        
            # Convert column names to lower case
            snpmask_df.columns = snpmask_df.columns.str.lower()
            vcf_df.columns = vcf_df.columns.str.lower()

            # Create a set from the 'chrom', 'pos', 'ref', 'alt' columns
            known_snps_set = set(zip(snpmask_df['chrom'], snpmask_df['pos'], 
                                    snpmask_df['ref'], snpmask_df['alt']))

            # Filter the vcf_df to only include rows not in the known_snps_set
            self.log.note(f'Input dataframe length: {len(vcf_df)}')

            vcf_df = vcf_df[~vcf_df[['chrom', 'pos', 'ref', 'alt']]
                                .apply(tuple, axis=1).isin(known_snps_set)]
 
            self.log.note(f'Filtered dataframe length: {len(vcf_df)}')

            return vcf_df

        except Exception as e:
            self.log.error(f"There was an error while filtering known snps. {e}")


class FeatureProduction:    
    def __init__(self, logger, name):
        self.log = logger
        self.name = name
    
    @staticmethod
    def _delta_snp_array(wtr: np.ndarray, wta:np.ndarray, mur: np.ndarray, mua: np.ndarray)-> np.ndarray:
        """
        Calculates delta SNP feature, which quantifies divergence in 
            read depths between the two bulks. More data driven than g-stat
            method. 

        Args: 
            wtr, wta, mur, mua (numpy array)
            Read depths of wt reference and alt reads
            and the read depths of mutant reference and alt reads 

        Returns: Delta-snp calculation, which is a quantification of 
            allelic segregation at each polymorphic site.
        """ 

        return ((wtr) / (wtr + wta)) - ((mur) / (mur + mua))
    
    @staticmethod
    def _g_statistic_array(wtr: np.ndarray, wta: np.ndarray, mur: np.ndarray, mua: np.ndarray)->np.ndarray:
        """
        Calculates g-statistic feature, which is a Chi squared-ish approach to 
        calculating read-depth divergence from expected values. 
        """ 

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

        return np.where(e1 * e2 * e3 * e4 == 0, 0.0, llr1 + llr2 + llr3 + llr4)

    def calculate_delta_snp_and_g_statistic(self, vcf_df: pd.DataFrame)-> pd.DataFrame:
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
            # calculate delta snps
            
            wt_ref = vcf_df['wt_ref'].values
            wt_alt = vcf_df['wt_alt'].values
            mu_ref = vcf_df['mu_ref'].values
            mu_alt = vcf_df['mu_alt'].values
            
            vcf_df['ratio'] = FeatureProduction._delta_snp_array(
                wt_ref, wt_alt, mu_ref, mu_alt
            )
            vcf_df['G_S'] = FeatureProduction._g_statistic_array(
                   wt_ref, wt_alt, mu_ref, mu_alt
            )
            vcf_df['RS_G'] = vcf_df['ratio'].values * vcf_df['G_S'].values

            self.log.success("Calculation of delta-SNP ratios and G-statistics was successful.")
            return vcf_df
        
        except Exception as e:
            self.log.error(f"An error occurred during calculation: {e}")

    def _create_mirrored_data(self, df_chrom, smooth_edges_bounds):
        self.log.attempt(f'Mirroring {smooth_edges_bounds} datapoints at chromosome ends to correct for loess edge bias...')
        try:
            # Get the positions
            positions = df_chrom['pos'].to_numpy()

            # Calculate the deltas
            deltas = np.diff(positions, prepend=0)

            # Calculate mirrored ends
            deltas_pos_inv = deltas[1:smooth_edges_bounds][::-1]
            deltas_neg_inv = deltas[-smooth_edges_bounds:-1][::-1]
            deltas_mirrored_ends = np.concatenate([deltas_pos_inv, deltas, deltas_neg_inv])

            # Calculate pseudo_pos
            pseudo_pos = np.cumsum(deltas_mirrored_ends)

            # Mirror the data
            df_chrom_inv_neg = df_chrom.iloc[1:smooth_edges_bounds][::-1]
            df_chrom_inv_pos = df_chrom.iloc[-smooth_edges_bounds:-1][::-1]
            df_chrom = pd.concat([df_chrom_inv_neg, df_chrom, df_chrom_inv_pos], ignore_index=False)

            # Add pseudo_pos to the DataFrame
            df_chrom['pseudo_pos'] = pseudo_pos
            
            self.log.success(f"Mirrored data created!")
            return df_chrom

        except Exception as e:
            self.log.warning(f"Error in creating mirrored data: {e}")
            return None
        
    def _fit_values(self, df_chrom, smoothing_function, loess_span):

            X = df_chrom['pseudo_pos'].values
                        
            for col in ['ratio', 'G_S', 'RS_G']:
                Y = df_chrom[col].values
                df_chrom[f'{col}_yhat'] = smoothing_function(Y, X, frac=loess_span)[:, 1]
            
            return df_chrom

    def _fit_single_chrom(self, df_chrom, chrom, smoothing_function, loess_span, smooth_edges_bounds):
        self.log.attempt(f"Fitting values in chromosome:{chrom}")
        try:
            
            df_chrom = self._create_mirrored_data(df_chrom, smooth_edges_bounds)
            df_chrom = self._fit_values(df_chrom, smoothing_function, loess_span)

            # Sort by 'chrom' and 'pos'
            df_chrom = df_chrom.sort_values(by=['chrom', 'pseudo_pos'])
            
            self.log.success(f"Values fit in chromosome: {chrom}") 
            return df_chrom

        except Exception as e:
            self.log.warning(f"There was an error while smoothing chrom:{chrom}:{e}")

        return None

    def _fit_chrom_facets(self, vcf_df, smoothing_function, loess_span, smooth_edges_bounds):
        df_list = []
        chrom_facets = vcf_df["chrom"].unique()

        for chrom in chrom_facets:
            df_chrom = vcf_df[vcf_df['chrom'] == chrom]
            result = self._fit_single_chrom(df_chrom, chrom, smoothing_function, loess_span, smooth_edges_bounds)
            
            if result is not None:
                df_list.append(result)

        return pd.concat(df_list).reset_index(drop=True)
    
    def fit_model(self, vcf_df: pd.DataFrame, smoothing_function, loess_span: float, smooth_edges_bounds: int)->pd.DataFrame:
        """
        LOESS smoothing of ratio and G-stat by chromosome
        
        Input: Cleaned dataframe with delta SNPs and G-stats calculated
        
        Returns: Dataframe containing LOESS fitted values for ratio, g-stat and 
        ratio-scaled g-stat
        """

        self.log.attempt("Initialize LOESS smoothing calculations.")
        self.log.attempt(f"span: {loess_span}, Edge bias correction: {smooth_edges_bounds}") 
        try:
            vcf_df = self._fit_chrom_facets(vcf_df, smoothing_function, loess_span, smooth_edges_bounds)
        
            self.log.success("LOESS smoothing calculations successful.")    
            
            return vcf_df
        
        except Exception as e:
            self.log.error( f"An error occurred during LOESS smoothing: {e}")
    
    @staticmethod
    def _null_models(smChr, smPseudoPos, sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, smoothing_function, loess_span: float, ratio_cutoff: float):
            '''
            This process is parallelized, so everything is pretty sanatized.

            shuffle - calculate statistics - mask positions below cutoff - sort 
            - fit - return 
            '''
            
            #Break the phenotype/genotype association by shuffling read depth data
            np.random.shuffle(sm_wt_ref)
            np.random.shuffle(sm_wt_alt)
            np.random.shuffle(sm_mu_ref)
            np.random.shuffle(sm_mu_alt)
            
            #calc delta snp
            smRatio = FeatureProduction._delta_snp_array(
                sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt
            )
            #calc g-stat
            smGstat = FeatureProduction._g_statistic_array(
                sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt
            )
            # calc ratio scaled g-stat
            smRS_G = smRatio * smGstat
            
            # Mask those values and positions below the delta snp ratio cutoff 
            # (makes null model fit trimmed data)
            mask = smRatio >= ratio_cutoff
            smPseudoPos = smPseudoPos[mask] 
            smRatio = smRatio[mask]
            smGstat = smGstat[mask]
            smChr = smChr[mask]
            smRS_G = smRS_G[mask]

            # Sort arrays by chromosome
            sorted_indices = np.argsort(smChr)
            smChr = smChr[sorted_indices]
            smPseudoPos = smPseudoPos[sorted_indices]
            smRatio = smRatio[sorted_indices]
            smGstat = smGstat[sorted_indices]
            smRS_G = smRS_G[sorted_indices]

            # Calculate smoothed values for each chromosome separately
            smRatio_y = []
            smGstat_y = []
            smRS_G_y = []
            for chrom, indices in groupby(enumerate(smChr), key=lambda x: x[1]):
                indices = [i for i, _ in indices]
                smRatio_y.extend(smoothing_function(smRatio[indices], smPseudoPos[indices], frac=loess_span)[:, 1])
                smGstat_y.extend(smoothing_function(smGstat[indices], smPseudoPos[indices], frac=loess_span)[:, 1])
                smRS_G_y.extend(smoothing_function(smRS_G[indices], smPseudoPos[indices], frac=loess_span)[:, 1])

            return smChr, smPseudoPos, smRatio, np.array(smRatio_y), smGstat, np.array(smGstat_y), smRS_G, np.array(smRS_G_y)

    def _initialize_array(self, smChr, smPos, smPseudoPos, shuffle_iterations, list_name):
        '''
        init structured numpy arrays in memory. Saves some overhead
        '''
        self.log.attempt(f"Initializing null model structured array: {list_name}")
        try: 
            dtype = [('chrom', 'U10'), ('pos', int), ('pseudo_pos', int), ('value', float, shuffle_iterations)]
            lst = np.zeros(len(smChr), dtype=dtype)
            lst['chrom'] = smChr
            lst['pos'] = smPos
            lst['pseudo_pos'] = smPseudoPos

            self.log.success(f"Null model structured array initialized: {list_name}")
            
            return lst
        
        except Exception as e:
            self.log.error(f"Initializing structured array for {list_name} null model failed: {e}")

    def calculate_null_models(self, vcf_df: pd.DataFrame, smoothing_function, loess_span: float, shuffle_iterations: int, ratio_cutoff: float)->tuple:
        '''
        Process that generates the null models. This is fairly resource hungry, 
        so it is run in parallel. For each position, shuffle_iterations(n) 
        bootstrapped values are calculated. 

        [...]
        This function needs work but is tricky to refactor because 
        staticmethods are needed to play nice with parallel processing.
        '''
        self.log.attempt('Bootstrapping to generate null models...')
        try:

            smChr = vcf_df['chrom'].to_numpy()
            smPos = vcf_df['pos'].to_numpy()
            smPseudoPos = vcf_df['pseudo_pos'].to_numpy()
            sm_wt_ref = vcf_df['wt_ref'].to_numpy()
            sm_wt_alt = vcf_df['wt_alt'].to_numpy()
            sm_mu_ref = vcf_df['mu_ref'].to_numpy()
            sm_mu_alt = vcf_df['mu_alt'].to_numpy()

            args = [smChr, smPos, smPseudoPos, shuffle_iterations]

            sm_ratio_array = self._initialize_array(*args, 'sm_ratio')
            sm_ratio_y_array = self._initialize_array(*args, 'sm_ratio_y')
            sm_g_stat_array = self._initialize_array(*args, 'sm_g_stat')
            sm_g_stat_y_array = self._initialize_array(*args, 'sm_g_stat_y')
            sm_ratio_scaled_g_array = self._initialize_array(*args, 'sm_ratio_scaled_g')
            sm_ratio_scaled_g_y_array = self._initialize_array(*args, 'sm_ratio_scaled_g_y')

            position_counts = {}
            position_indices = {}
            unique_index = 0
            for chrom, pseudoPos in zip(smChr, smPseudoPos):
                if (chrom, pseudoPos) not in position_counts:
                    position_counts[(chrom, pseudoPos)] = 0
                    position_indices[(chrom, pseudoPos)] = unique_index
                    unique_index += 1
            
            self.log.note(f"Distributing bootstrapping calculations to number of cores:{THREADS_LIMIT}")
            self.log.note(f"Bootstrapped vales per position:{shuffle_iterations}")
            total_values = len(smChr)*shuffle_iterations
            self.log.note(f"Total values to be generated:{total_values}")

            args = [(smChr, smPseudoPos, sm_wt_ref, sm_wt_alt, sm_mu_ref, sm_mu_alt, smoothing_function, 
                loess_span, ratio_cutoff) for _ in range(THREADS_LIMIT)]

            iteration = 0
            total_values_added = 0
            with Pool() as pool:
                while not all(count == shuffle_iterations for count in position_counts.values()):
                    results = pool.starmap(FeatureProduction._null_models, args)
                    iteration += 1
                    self.log.note(f"Iteration: {iteration}")

                    self.log.note(f"Progress: {total_values_added}/{total_values}")
                    
                    for result in results:
                        for chrom, pseudo_pos, ratio, ratio_y, gstat, gstat_y, rs_g, rs_g_y in zip(*result):
                            try:
                                key = (chrom, pseudo_pos)
                                if position_counts[key] < shuffle_iterations:
                                    index = position_counts[key]
                                    array_index = position_indices[key]
                                    sm_ratio_array['value'][array_index, index] = ratio
                                    sm_ratio_y_array['value'][array_index, index] = ratio_y
                                    sm_g_stat_array['value'][array_index, index] = gstat
                                    sm_g_stat_y_array['value'][array_index, index] = gstat_y
                                    sm_ratio_scaled_g_array['value'][array_index, index] = rs_g
                                    sm_ratio_scaled_g_y_array['value'][array_index, index] = rs_g_y

                                    position_counts[key] += 1
                                    total_values_added += 1
                            except Exception as e:
                                self.log.warning(f"Counts not added at key:{key}, index:{index}, array index:{array_index}")
                                continue

            self.log.success('Bootstrapping complete!')

            return sm_ratio_array, sm_ratio_y_array, sm_g_stat_array, sm_g_stat_y_array, sm_ratio_scaled_g_array, sm_ratio_scaled_g_y_array

        except Exception as e:
            self.log.error(f'Bootstrapping to generate empirical cutoffs failed:{e}')

            return None, None, None, None, None, None, None

    def aggregate_unsmoothed_values(self, null_models):
        '''
        Aggregates unsmooted values into a list. This is done because there is no
        position or chromosome factor for unsmoothed values. Therefore, the null
        models for unsmoothed values are aggregated. 
        '''
       # Unpack null_models
        sm_ratio_array, sm_ratio_y_array, sm_g_stat_array, sm_g_stat_y_array, sm_ratio_scaled_g_array, sm_ratio_scaled_g_y_array = null_models


        # Define a list of structured arrays
        arrays = [sm_ratio_array, sm_g_stat_array, sm_ratio_scaled_g_array]

        # Initialize an empty list for each structured array
        aggregated_values_lsts = [[] for _ in arrays]

        # Iterate over the structured arrays and the corresponding empty lists
        for array, values in zip(arrays, aggregated_values_lsts):
            # Concatenate all values from the structured array into a single list
            values.extend(array['value'].flatten())

        sm_ratio_lst, sm_g_stat_lst, sm_ratio_scaled_g_lst = aggregated_values_lsts
        
        return sm_ratio_lst, sm_ratio_y_array, sm_g_stat_lst, sm_g_stat_y_array, sm_ratio_scaled_g_lst, sm_ratio_scaled_g_y_array

    def _calculate_percentile(self, value, sorted_array):
        '''
        Quick function to calc percentile. Built in libraries for this were
        painfully slow
        '''
        sorted_list = (
            sorted_array.tolist() 
                if isinstance(sorted_array, np.ndarray) 
                else sorted_array
        )
        
        idx = bisect.bisect_left(sorted_list, value)
        percentile = idx / len(sorted_list)
        
        return percentile

    def _assign_percentiles(self, vcf_df, columns_and_lists):
        for column, lst in columns_and_lists:
            self.log.note(f"assigning {column} percentiles based on null model")
            sorted_lst = sorted(lst)
            vcf_df[column + '_percentile'] = vcf_df[column].apply(
                self._calculate_percentile, sorted_array=sorted_lst
            )

    def _assign_yhat_percentiles(self, vcf_df, columns_and_arrays_yhat):
        for column, lst in columns_and_arrays_yhat:
            self.log.note(f"assigning {column} percentiles based on null model")
            column_array = vcf_df[column].to_numpy()

            for i, value in enumerate(column_array):
                if i < len(lst):
                    values_array = sorted(lst[i]['value'])
                    self._calculate_and_assign_percentiles(
                        i, column, value, values_array, vcf_df
                    )
                else:
                    self.log.warning(f"Skipping index {i} due to missing data")

    def _calculate_and_assign_percentiles(self, i, column, value, values_array, vcf_df):
        vcf_df.at[i, column + '_percentile'] = self._calculate_percentile(
            value, values_array
        )
        vcf_df.at[i, column + '_null_5'] = np.percentile(values_array, 5)
        vcf_df.at[i, column + '_null_25'] = np.percentile(values_array, 25)
        vcf_df.at[i, column + '_null_50'] = np.percentile(values_array, 50)
        vcf_df.at[i, column + '_null_75'] = np.percentile(values_array, 75)
        vcf_df.at[i, column + '_null_95'] = np.percentile(values_array, 95)
    
    def label_df_with_percentiles(self, vcf_df: pd.DataFrame, null_models)->pd.DataFrame:
        '''
        Labels the dataframe with the percetiles needed to both generate the 
        likely candidates list, and to create the null model ribbons on the 
        fitted graphs
        '''
        self.log.attempt(f"Labeling dataframe with percentiles based on null models...")
        try:
            (sm_ratio_lst, sm_ratio_y_array, sm_g_stat_lst, sm_g_stat_y_array, 
                sm_ratio_scaled_g_lst, sm_ratio_scaled_g_y_array
            ) = null_models

            columns_and_arrays_yhat = [('ratio_yhat', sm_ratio_y_array), 
                                    ('G_S_yhat', sm_g_stat_y_array), 
                                    ('RS_G_yhat', sm_ratio_scaled_g_y_array)]

            columns_and_lists = [('ratio', sm_ratio_lst), 
                                ('G_S', sm_g_stat_lst), 
                                ('RS_G', sm_ratio_scaled_g_lst)]

            self._assign_percentiles(vcf_df, columns_and_lists)
            self._assign_yhat_percentiles(vcf_df, columns_and_arrays_yhat)
            
            self.log.success(f"Dataframe labeled with percentiles")
            return vcf_df

        except Exception as e:
            self.log.error(f"Labeling dataframe failed:{e}")

    def remove_extra_data(self, vcf_df, null_models):
        '''
        Removes the extra data created during the edge bias correction process. 
        This data is retained until now so that the null models include the 
        correction process as well
        '''
        self.log.attempt(f"Attempting to remove psuedo positions and mirrored data from dataframe and stuctured array")
        try:
            vcf_df = vcf_df.drop(columns='pseudo_pos')
            vcf_df = vcf_df.dropna(subset=['ratio_yhat', 'G_S_yhat', 'RS_G_yhat'])
            vcf_df = vcf_df.drop_duplicates(subset=['chrom', 'pos'])

            (sm_ratio_lst, sm_ratio_y_array, sm_g_stat_lst, sm_g_stat_y_array, 
                sm_ratio_scaled_g_lst, sm_ratio_scaled_g_y_array
            ) = null_models

            # For each structured array in null_models, drop 'pseudo_pos' and remove duplicates
            new_null_models = []
            for arr in null_models:
                if isinstance(arr, np.ndarray) and 'pseudo_pos' in arr.dtype.names:
                    new_arr = np.copy(arr[['chrom', 'pos', 'value']])  # create a new array with only 'chrom', 'pos', 'values' fields
                    new_arr = np.unique(new_arr, axis=0)  # remove duplicates
                    new_null_models.append(new_arr)

                else:
                    new_null_models.append(arr)

            self.log.success("Mirrored data removed!")
            return vcf_df, new_null_models
        
        except Exception as e:
            self.log.error(f"There was an error while removing mirrored data and psuedo positions:{e}")


class TableAndPlots:
    def __init__(self, logger, name, analysis_out_prefix):
        self.log = logger
        self.name = name
        self.analysis_out_prefix = analysis_out_prefix
        
        self.impact_mapping = {
            'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1
        }
   
    def _theme(self):
        return theme(
            #general
            plot_background = element_rect(fill = 'white'),
            panel_background = element_rect(fill = 'white'),
            panel_border = element_rect(color = 'grey', size = 0.5, alpha=0.6),
            panel_grid = element_line(size = 0.3, color = 'grey', alpha=0.5),
            panel_grid_minor = element_line(size = 0.5, color = 'grey', 
                alpha=0.5, linetype = 'dotted'
            ),
            
            # Set the title text properties
            title = element_text(color = 'black', size = 10, weight = 'bold'),
            
            # Set the axis text properties
            axis_text = element_text(color = 'black', size = 8),
            axis_ticks_major = element_line(color='grey', size = 0.1, alpha=0.5),
            
            # Set the legend text properties
            legend_text = element_text(color = 'black', size = 8),
            legend_title = element_text(color = 'black', size = 10, 
                weight = 'bold'
            ),

            # Set the strip properties
            strip_background = element_rect(fill = 'white', color = 'grey',
                size = 0.2
            ),
            strip_text = element_text(color = 'grey', size = 8, weight = 'bold'),        
        )

    def _create_histograms(self, dataframes):
        return {
            'sm_ratio_df': ["Ratio Frequency Distribution", f"{self.analysis_out_prefix}-NULL-ratio.png"],
            'sm_g_stat_df': ["G-statistic Frequency Distribution", f"{self.analysis_out_prefix}-NULL_g_stat.png"],
            'sm_ratio_scaled_g_df': ["Ratio-scaled G-statistic Frequency Distribution", f"{self.analysis_out_prefix}-NULL_ratio_scaled_g.png"],
            'sm_ratio_y_df': ["LOESS smoothed Ratio Frequency Distribution", f"{self.analysis_out_prefix}-NULL_ratio_y.png"],
            'sm_g_stat_y_df': ["LOESS smoothed G-statistic Frequency Distribution", f"{self.analysis_out_prefix}-NULL-g_stat_y.png"],
            'sm_ratio_scaled_g_y_df': ["LOESS smoothed Ratio-scaled G-statistic Frequency Distribution", f"{self.analysis_out_prefix}-NULL-ratio_scaled_g_y.png"]
        }

    def _sample_and_explode_df(self, df):
        num_samples = min(20, len(df))
        sampled_df = df.sample(n=num_samples, random_state=1)
        
        exploded_df = sampled_df.explode('value')
        exploded_df['value'] = pd.to_numeric(exploded_df['value'])
        cutoffs = exploded_df.groupby(['chrom', 'pos'])['value'].quantile(0.95).reset_index()
        cutoffs.rename(columns={'value': 'cutoff'}, inplace=True)
        exploded_df = pd.merge(exploded_df, cutoffs, on=['chrom', 'pos'])
        
        return exploded_df

    def _plot_histogram(self, data, title, filename, faceted):
        
        if len(data) > 100_000 and not faceted:
            self.log.note(f"Data for histogram exceeds 100,000. Subsampling null model histogram down to 100,000 for faster plotting....")
            data = data.sample(n=100_000, random_state=1)

        plot = (ggplot(data, aes(x='value')) +
                geom_histogram(fill='#69b3a2', color="#000000", alpha=0.7) +
                ggtitle(title) +
                self._theme())

        if faceted:
            plot = (plot +
                    facet_wrap('~chrom + pos') +
                    geom_vline(data=aes(xintercept='cutoff'), color="red", linetype="dashed", size=0.3) +
                    theme(panel_spacing=0.01) 
                    )
        else:
            cutoff = np.percentile(data['value'], 95)
            plot = (plot + 
                geom_vline(xintercept=cutoff, color="red", linetype="dashed", size=0.5)
            )
        
        # Save the plot
        plot.save(filename, dpi=600)

    def plot_null_histograms(self, null_models):
        '''
        Function that handles plotting the null model histograms. The histograms
        for the null models that are lists (The ones without "y" in their name)
        are plotted as a single histogram. The calculation of their indiviual stats
        are not relient upon position, so they just have a single global null model

        Those with "y" it there name are structured arrays with chrom pos value, 
        and each position contains a null model, because the fitting procress causes
        closer positioned snps to have biased results. To visualize a sample of the
        position wise models, a subsample of 20 are taken from the array and graphed. 

        You will notice that there is more variance in those positions with close neighbors, 
        which is also visiable in the null model percentile ribbons in the main graphs.  
        '''
        # Unpack null_models
        sm_ratio_lst, sm_ratio_y_array, sm_g_stat_lst, sm_g_stat_y_array, sm_ratio_scaled_g_lst, sm_ratio_scaled_g_y_array = null_models

        # Create dataframes
        dataframes = {
            'sm_ratio_df': pd.DataFrame(sm_ratio_lst, columns=['value']),
            'sm_g_stat_df': pd.DataFrame(sm_g_stat_lst, columns=['value']),
            'sm_ratio_scaled_g_df': pd.DataFrame(sm_ratio_scaled_g_lst, columns=['value']),
            'sm_ratio_y_df': pd.DataFrame(sm_ratio_y_array.tolist(), columns=['chrom', 'pos', 'value']),
            'sm_g_stat_y_df': pd.DataFrame(sm_g_stat_y_array.tolist(), columns=['chrom', 'pos', 'value']),
            'sm_ratio_scaled_g_y_df': pd.DataFrame(sm_ratio_scaled_g_y_array.tolist(), columns=['chrom', 'pos', 'value'])
        }

        # Create histograms
        histograms = self._create_histograms(dataframes)

        for df_name, histogram in histograms.items():
            if '_y_' in df_name:
                exploded_df = self._sample_and_explode_df(dataframes[df_name])
                self._plot_histogram(exploded_df, histogram[0], histogram[1], faceted=True)
            else:
                self._plot_histogram(dataframes[df_name], histogram[0], histogram[1], faceted=False)
    
    def _sort_candidates(self, df):
        '''
        Sorts likely candidates df by impact rank and percentile rankings.
        '''
        df['impact_rank'] = df['snpimpact'].apply(
            lambda x: max(self.impact_mapping.get(i, 0) for i in x.split(':'))
        )
        sorted_df = df.sort_values(
            by=['impact_rank', 
                'RS_G_yhat_percentile', 
                'G_S_yhat_percentile', 
                'ratio_yhat_percentile', 
                'RS_G_percentile', 
                'G_S_percentile', 
                'ratio_percentile'], 
            ascending=False
        )

        return sorted_df.drop('impact_rank', axis=1)

    def _get_likely_candidates(self, vcf_df):
        '''
        Creates a new df of likely candidates. Those that are above the 95 
        percentile of their null models are considered. 
        '''
        conditions = (
            (vcf_df['ratio_percentile'] > 0.95) |
            (vcf_df['ratio_yhat_percentile'] > 0.95) |
            (vcf_df['G_S_percentile'] > 0.95) |
            (vcf_df['G_S_yhat_percentile'] > 0.95) |
            (vcf_df['RS_G_percentile'] > 0.95) |
            (vcf_df['RS_G_yhat_percentile'] > 0.95) 
        )

        likely_candidates = vcf_df[conditions].copy()

        if likely_candidates.empty:
            self.log.nocandidates()
        
        return likely_candidates

    def _save_to_csv(self, df, filename):
        '''
        Save the input dataframe to file
        '''
        df.to_csv(f"{self.analysis_out_prefix}_{filename}.csv", index=False)

    def process_and_save_candidates(self, vcf_df):
        '''
        Retrieve candidates, sort them and save them. 
        '''
        
        likely_candidates = self._sort_candidates(
            self._get_likely_candidates(vcf_df)
        )
        self._save_to_csv(vcf_df, 'all')
        self._save_to_csv(likely_candidates, 'likely_candidates')

    def _create_plot(self, vcf_df, y_column, title_text, ylab_text, cutoff_value=None, lines=False):
        '''
        Where the majority of the main plot specifications are done. _theme
        contains the general theming of the plots. 
        
        Those plots with y in the title(the fitted values) will contain a nested 
        ribbon showing the 5-25-50-75-95 percentiles, so that it is easy to see 
        the distrubtion of the data if there was no genotype/phenotype link.
        '''
        try:
            mb_conversion_constant = 0.000001
            vcf_df['pos_mb'] = vcf_df['pos'] * mb_conversion_constant
            chart = ggplot(vcf_df, aes('pos_mb', y=y_column))
            title = ggtitle(title_text)
            axis_x = xlab("Position (Mb)")
            axis_y = ylab(ylab_text)

            plot = (chart
                    + geom_point(color='goldenrod', size=0.8)
                    + self._theme()
                    + facet_grid('. ~ chrom', space='free_x', scales='free_x')
                    + title
                    + axis_x
                    + axis_y
                    + theme(panel_spacing=0.025)
            )
            
            if cutoff_value is not None:
                cutoff = geom_hline(yintercept=cutoff_value, color='red',
                                    linetype="dashed", size=0.3)
                plot += cutoff

            cutoff_column = y_column + '_null_50'
            if cutoff_column in vcf_df.columns:

                # Calculate the middle position for each 'chrom'
                mid_pos_df = vcf_df.groupby('chrom')['pos_mb'].apply(lambda x: (x.max() + x.min()) / 2).reset_index()
                mid_pos_df.columns = ['chrom', 'mid_pos']

                # Calculate the average of the 50th percentile for each 'chrom'
                avg_50_percentile_df = vcf_df.groupby('chrom')[cutoff_column].mean().reset_index()
                avg_50_percentile_df.columns = ['chrom', 'avg_50_percentile']

                # Merge middle position and average 50th percentile data frames
                label_df = pd.merge(mid_pos_df, avg_50_percentile_df, on='chrom')


                plot += geom_ribbon(aes(ymin=y_column +'_null_5', ymax=y_column +'_null_95'), fill='gray', alpha=0.2)
                plot += geom_ribbon(aes(ymin=y_column +'_null_25', ymax=y_column +'_null_75'), fill='gray', alpha=0.2)
                plot += geom_line(aes(y=y_column +'_null_50'), size=0.4, color='black', alpha=0.2)
                plot += geom_line(aes(y=y_column +'_null_50'), size=0.4, color='gray', alpha=0.2)
                plot += geom_text(data=label_df, mapping=aes(x='mid_pos', y='avg_50_percentile', label="r'$H_{0}^{*}$'"), size=6, ha='center', color='purple', alpha=1)
                # Add "null" annotation to the average 50th percentile position in each facet

            plot += geom_point(color='goldenrod', size=0.9)

            if lines:
                plot += geom_line(color='blue', size=0.7, alpha=0.85)
                plot += geom_line(color='blue', size=0.7, alpha=0.9)

            return plot

        except Exception as e:
            self.log.error(f"Plot creation failed for {self.name}, column {y_column}: {e}") 
            
            return None
    
    def _create_and_save_plot(self, df, y_column, title_text, ylab_text, cutoff=None, lines=False):
        '''
        organizes the creation and saving of plots. An intermedary between 
        generate_plots and _create_plot
        '''
        plot = self._create_plot(
            df, y_column, title_text, ylab_text, cutoff, lines
        )
        
        if plot:
            plot_path = f"{self.analysis_out_prefix}_{y_column.lower()}.png"
            plot.save(filename=plot_path, height=6, width=8, units='in', dpi=500)

    def generate_plots(self, vcf_df):
        plot_scenarios = [
            ('G_S', 'G-statistic', 'G-statistic'),
            ('G_S_yhat', 'Fitted G-statistic', 'Fitted G-statistic', None, True),
            ('RS_G', 'Ratio-scaled G statistic', 'Ratio-scaled G-statistic'),
            ('ratio', 'Delta SNP ratio', 'Ratio'),
            ('ratio_yhat', 'Fitted Delta SNP ratio', 'Fitted delta SNP ratio', None, True),
            ('RS_G_yhat', 'Fitted ratio-scaled G statistic', 'Fitted Ratio-scaled G-statistic', None, True),
            ('ratio_percentile', 'Delta SNP ratio percentile', 'Ratio percentile', 0.95),
            ('G_S_percentile', 'G-statistic percentile', 'G-statistic percentile', 0.95),
            ('RS_G_percentile', 'Ratio-scaled G statistic percentile', 
             'Ratio-scaled G-statistic percentile', 0.95),
        ]

        for scenario in plot_scenarios:
            self._create_and_save_plot(vcf_df, *scenario)