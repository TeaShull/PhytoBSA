from settings.config import INPUT_DIR, REFERENCE_DIR, MODULES_DIR, VCF_GEN_SCRIPT

import subprocess

from modules.utilities_bsa_analysis import BSAAnalysisUtilities
from modules.utilities_general import FileUtilities
from modules.utilities_logging import LogHandler

class ThaleBSAParentFunctions:
    def __init__(self, logger, dict_utils):
        self.log = logger
        self.dict_utils = dict_utils

    def vcf_generation(self):
        """
        Input: Experiment dictionary as well as paths and variables needed 
        to run VCFgen.sh. Subprocess VCFgen.sh takes raw reads(wild-type(wt) 
        and mutant(mu)), either single or paired end and generates VCF table 
        *.noknownsnps.table.

        Nonetypes - if nothing is provided, line_dict will be 
        detected automatically from ./input and variables will be sourced from 
        the variables.py module.
        
        Returns: updated line_dict containing the paths to 
        the noknownsnps.table.
        """

        # Generate the knownSnps .vcf file path
        self.log.note('Beginning VCF generation process for line_dict')
        try:
            for line, details in self.dict_utils.line_dict.items():
                # Print variables contained in line_dict to core_log
                self.log.delimiter('Experiment Dictionary passed to vcf_generation function')
                self.log.print(f'Line:{line}')
                for details_key, value in details.items():
                    self.log.print(f"{details_key}: {value}")
                self.log.print(' ')
                '''
                UNPACK REQUIRED VARIABLES: VCF_GEN_SCRIPT(subprocess_VCFgen.sh)
                '''
                wt_input = details['wt_input']
                mu_input = details['mu_input']
                pairedness = details['pairedness']
                output_dir_path = details['output_dir_path']
                output_prefix = details['output_prefix']
                vcf_table_path = details['vcf_table_path']
                reference_genome_name = details['reference_genome_name']
                snpEff_species_db = details['snpEff_species_db']
                reference_genome_source = details['reference_genome_source']
                known_snps_path = details['known_snps_path']
                threads_limit = details['threads_limit']
                call_variants_in_parallel = details['call_variants_in_parallel']
                cleanup = details['cleanup']
                
                # generate log instance, add run info to sql db
                vcf_log = LogHandler(f'vcf_{line}')
                details['vcf_ulid'] = vcf_log.ulid
                vcf_log.add_db_record(line, self.log.ulid)
                self.log.note(f'Logging for VCF subprocess Initialized and ulid added to log database. Path: {vcf_log.log_path}')
                
                # Construct args to pass variables to VCF_GEN_SCRIPT(subprocess_VCFgen.sh)
                args = (
                    vcf_log.ulid,
                    line, 
                    INPUT_DIR,
                    wt_input,
                    mu_input,
                    pairedness,
                    output_dir_path,
                    output_prefix,
                    vcf_table_path,
                    REFERENCE_DIR,
                    reference_genome_name, 
                    snpEff_species_db,
                    reference_genome_source, 
                    known_snps_path,
                    threads_limit,
                    call_variants_in_parallel,
                    cleanup
                )

                # Construct the command for VCFgen.sh, passing the above variables
                cmd = f"{VCF_GEN_SCRIPT} {' '.join(map(str, args))}"
                # Run vcfgen shell subprocess.
                process = subprocess.Popen(
                    cmd, cwd=MODULES_DIR, shell=True, stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, text=True
                )
                # Iterate over stdout from process and log
                for line in process.stdout:
                    vcf_log.bash(line.strip())
                process.wait()
                self.log.note(f"VCF file generated for {line}.") 
                self.log.note(f"Log saved to {log_path}")
                self.log.note(f'VCF table path added to experiments_dictionary: {vcf_table_path}')
                
            self.log.success("VCF file generation process complete")

        except Exception as e:
            self.log.fail(f"Error while generating the VCF file for {line}: {e}")

    def bsa_analysis(self):
        '''
        Parent function for running BSA analysis. Given an line_dict
        containing the vcf_table_path, line_name and segregation_type, it will 
        output a list of candidate mutations and assotiated graphs, which will 
        help narrow down regions of interest.  
        
        Dictionary Requirements:
        required 
        [key]: line_name
        [details] vcf_table_path, core_ulid, 
        [details] segregation_type (what segrigation pattern? R or D?)

        optional
        [details] vcf_ulid
        
        vcf_table_path must be the path to a vcf table produced from VCFgen.sh. 
        support for VCF tables which are configured differently is on the roadmap. 
        
        '''
        self.log.attempt("Attempting to perform data analysis...")
        try:
            # Run analysis
            for line, details in self.dict_utils.line_dict.items():
                # Print variables contained in line_dict to core_log
                self.log.delimiter('Experiment Dictionary passed to bsa_analysis function:')
                self.log.print(f'Line:{line}')
                for details_key, value in details.items():
                    self.log.print(f"   {details_key}: {value}")
                self.log.print(' ')
                '''
                UNPACK REQUIRED VARIABLES: BSA analysis
                '''
                vcf_table_path = details['vcf_table_path']
                vcf_ulid = details['vcf_ulid']
                segregation_type = details['segregation_type']


                # Configure an analysis logger for each line.
                analysis_log = LogHandler(f'analysis_{line}')
                self.log.note(f'Analysis log initialized. Path: {analysis_log.log_path}')
                analysis_log.add_db_record(line, self.ulid, vcf_ulid)

                # FileUtilites instance that logs to analysis_log
                file_utils = FileUtilities(analysis_log)

                '''
                Initiate bsa analysis utilities. vcf_df is passed along and 
                updated in the BSAAnalysisUtilities class. It goes in, and 
                really doesn't need to come out, except for clarity's sake. 
                ''' 
                vcf_df = file_utils.load_vcf_table(vcf_table_path, line) #load vcf_df
                bsa_analysis_utils = BSAAnalysisUtilities(
                    line, vcf_df, vcf_ulid, analysis_log
                )
                ## Filter genotypes based on segregation pattern
                bsa_analysis_utils.filter_genotypes(segregation_type)
                
                ## data cleaning and orginization
                bsa_analysis_utils.filter_ems_mutations()
                bsa_analysis_utils.drop_indels()
                
                ## Feature production
                bsa_analysis_utils.calculate_delta_snp_and_g_statistic()
                bsa_analysis_utils.drop_na()
                bsa_analysis_utils.loess_smoothing()
                bsa_analysis_utils.calculate_empirical_cutoffs()
                
                ## Saving and plotting outputs
                bsa_analysis_utils.sort_save_likely_candidates()
                bsa_analysis_utils.generate_plots()
            self.log.success("Data analysis complete")
        
        except Exception as e:
            self.log.fail(f"Error during data analysis: {e}")
