import subprocess
from settings.config import (INPUT_DIR, REFERENCE_DIR, MODULES_DIR, 
                            VCF_GEN_SCRIPT)
from modules.utilities_bsa_analysis import BSAAnalysisUtilities
from modules.utilities_general import FileUtilities
from modules.utilities_logging import LogHandler


def get_details_data(details):
    return (details.get('wt_input'), details.get('mu_input'), 
           details.get('pairedness'), details.get('output_dir_path'), 
           details.get('output_prefix'), details.get('vcf_table_path'), 
           details.get('reference_genome_name'))

def get_args(ulid, line, wt, mu, paired, output_path, output_prefix, vcf_path, 
           ref_genome_name, species_db, genome_source, snps_path, limit, is_parallel, 
           is_cleanup):
    return (ulid, line, INPUT_DIR, wt, mu, paired, output_path, 
           output_prefix, vcf_path, REFERENCE_DIR, ref_genome_name, 
           species_db, genome_source, snps_path, limit, is_parallel,
           is_cleanup)

class ThaleBSAParentFunctions:
    def __init__(self, logger, dict_utils):
        self.log = logger
        self.dict_utils = dict_utils

    def vcf_generation(self):
        self.log.note('Beginning VCF generation process')
        try:
            for line, details in self.dict_utils.line_dict.items():
                self.log.delimiter('Experiment Dictionary passed')
                self.log.print(' ')
                wt, mu, pairedness, output_dir_path, output_prefix, vcf_table_path, ref_genome_name = get_details_data(details)
                vcf_log = LogHandler(f'vcf_{line}')
                details['vcf_ulid'] = vcf_log.ulid
                vcf_log.add_db_record(line, self.log.ulid)
                
                args = get_args(vcf_log.ulid, line, wt, mu, pairedness, 
                             output_dir_path, output_prefix, vcf_table_path, 
                             ref_genome_name, snpEff_species_db,reference_genome_source,
                             known_snps_path, threads_limit, call_variants_in_parallel,
                             cleanup)
                cmd = f"{VCF_GEN_SCRIPT} {' '.join(map(str, args))}"
                self.run_command(cmd, line)
            self.log.success("VCF file generation process complete")

        except Exception as e:
            self.log.fail(f"Error in VCF file gen: {line}: {e}")

    def run_command(self, cmd, line):
        process = subprocess.Popen(cmd, cwd=MODULES_DIR, shell=True, 
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
                            text=True)
        vcf_log = LogHandler(f'vcf_{line}')
        for line in process.stdout: 
            vcf_log.bash(line.strip())
        process.wait()
        self.log.note(f"VCF file generated for {line}") 

    def bsa_analysis(self):
        self.log.attempt("Attempting to perform data analysis...")
        try:
            for line, details in self.dict_utils.line_dict.items():
                self.log.delimiter('Experiment Dictionary pass')
                self.log.print(' ')
                self.run_analysis(line, details)
            self.log.success("Data analysis complete")
        
        except Exception as e:
            self.log.fail(f"Error during data analysis: {e}")
  
    def run_analysis(self, line, details):
        vcf_table_path = details.get('vcf_table_path')
        vcf_ulid = details.get('vcf_ulid')
        segregation_type = details.get('segregation_type')

         ## Filter genotypes based on segregation pattern
        bsa_analysis_utils = self.init_analysis_utils(line)
        bsa_analysis_utils.filter_genotypes(segregation_type)

        # FileUtilites instance that logs to analysis_log
        file_utils = self.inti_analysis_log_and_add_record(line, vcf_ulid)
        vcf_df = file_utils.load_vcf_table(vcf_table_path, line) #load vcf_df
        vcf_df = self.update_bsa_analysis_utils_vcf_df(bsa_analysis_utils, vcf_df, vcf_ulid)

        ## data cleaning and orginization
        vcf_df = bsa_analysis_utils.do_filtering(vcf_df)

        ## Feature production
        vcf_df = self.do_vcf_processing(bsa_analysis_utils, vcf_df)

        ## Saving and plotting outputs
         bsa_analysis_utils.do_outputs()
        
    def init_analysis_utils(self, line):
        analysis_log = LogHandler(f'analysis_{line}')
        self.log.note(f'Analysis log initialized')
        return BSAAnalysisUtilities(line, None, None, analysis_log)
        
    def inti_analysis_log_and_add_record(line, vcf_ulid):
        # Initialize an analysis logger for each line
        analysis_log = LogHandler(f'analysis_{line}')
        analysis_log.add_db_record(line, self.log.ulid, vcf_ulid)
        return FileUtilities(analysis_log)

    def update_bsa_analysis_utils_vcf_df(bsa_analysis_utils, analysis_log, vcf_table_path, line):
        # Initiate bsa analysis utilities and add vcf_df
        bsa_analysis_utils = BSAAnalysisUtilities(
            line, None, None, analysis_log
        )
        vcf_df = file_utils.load_vcf_table(vcf_table_path, line) #load vcf_df
        bsa_analysis_utils.vcf_df = vcf_df
        return bsa_analysis_utils

    def do_filtering(bsa_analysis_utils, vcf_df):
        bsa_analysis_utils.filter_ems_mutations()
        bsa_analysis_utils.drop_indels()
        return vcf_df

    def do_vcf_processing(bsa_analysis_utils, vcf_df):
        bsa_analysis_utils.calculate_delta_snp_and_g_statistic()
        bsa_analysis_utils.drop_na()
        bsa_analysis_utils.loess_smoothing()
        bsa_analysis_utils.calculate_empirical_cutoffs()
        return vcf_df
