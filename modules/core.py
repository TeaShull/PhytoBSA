from settings.config import INPUT_DIR, REFERENCE_DIR, MODULES_DIR

import subprocess

from modules.utilities_bsa_analysis import BSAAnalysisUtilities
from modules.utilities_general import FileUtilities
from modules.utilities_logging import LogHandler

class ThaleBSAParentFunctions:
    def __init__(self, logger):
        self.log = logger 

    def vcf_generation(self, experiment_dictionary=None,
                   reference_genome_name=None, snpEff_species_db=None,
                   reference_genome_source=None, threads_limit=None,
                   cleanup=None, known_snps=None)->dict:
        """
        Input: Experiment dictionary as well as paths and variables needed 
        to run VCFgen.sh. Subprocess VCFgen.sh takes raw reads(wild-type(wt) 
        and mutant(mu)), either single or paired end and generates VCF table 
        *.noknownsnps.table.

        Nonetypes - if nothing is provided, experiment_dictionary will be 
        detected automatically from ./input and variables will be sourced from 
        the variables.py module.
        
        Returns: updated experiment_dictionary containing the paths to 
        the noknownsnps.table.
        """
        core_ulid=self.log.ulid
        try:    
            self.log.attempt('Checking if experiment dictionary exists...')
            if not experiment_dictionary:
                self.log.warning('Experiment details undefined. Auto generating from files in ./input folder')
                file_utils = FileUtilities(self.log)
                experiment_dictionary = file_utils.experiment_detector()
            else:
                self.log.note('Experiment details provided. Checking other variables...')

            check_variables=file_utils.check_vcfgen_variables(
                reference_genome_name, 
                snpEff_species_db, 
                reference_genome_source, 
                threads_limit, 
                cleanup, 
                known_snps
            ) 
            if check_variables is not None:
                (
                    reference_genome_name, 
                    snpEff_species_db, 
                    reference_genome_source, 
                    threads_limit, 
                    cleanup, 
                    known_snps
                ) = check_variables

        except Exception as e:
            self.log.fail(f'Parsing variables for subprocess_VCFgen.sh failed:{e}')
        
        self.log.note('Beginning VCF generation process for experiment_dictionary')
        try:
            for key, value in experiment_dictionary.items():
                self.log.attempt(f"Generating VCF file for {key}...")
                self.log.delimiter(f"Shell [sh] VCF generator for {key} beginning...")

                current_line_name = key
                # generate log instance, add run info to sql db
                vcf_log = LogHandler(f'vcf_{current_line_name}')

                self.log.note(f'Logging for VCF subprocess Initialized. Path: {vcf_log.log_path}')
                experiment_dictionary[current_line_name]['vcf_ulid'] = vcf_log.ulid
                vcf_log.add_db_record(current_line_name, core_ulid)
                
                #Generate file paths needed for vcf generation
                file_utils = FileUtilities(vcf_log)
                (
                    output_dir_path,
                    output_prefix,
                    vcf_table_path, 
                    vcfgen_script_path, 
                    known_snps_path
                ) = file_utils.generate_vcf_file_paths(current_line_name, vcf_log, known_snps)
               
                #save vcf_table_path to the experiment dictionary for downstream use
                value['vcf_table_path'] = vcf_table_path

                # Retrieve allele and file input info from experiment_dictionary
                allele = value['allele'] #Recessive or dominant?
                pairedness = value['pairedness'] #Paired-end or single?
                # Pull input files from dictionary
                wt_input = ' '.join(value['wt'])
                wt_input = f'"{wt_input}"' #Wild-type bulk input files
                self.log.note(f"wt_input:{wt_input}")
                mu_input = ' '.join(value['mu'])
                mu_input = f'"{mu_input}"'#Mutant bulk input files
                self.log.note(f"mu_input:{mu_input}")

                # Construct args to pass variables to VCFgen.sh.
                args = (
                    vcf_log.ulid,
                    current_line_name, 
                    allele,
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
                    cleanup
                )
                # Construct the command for VCFgen.sh, passing the above variables
                print(vcfgen_script_path)
                cmd = f"{vcfgen_script_path} {' '.join(map(str, args))}"
                # Run vcfgen shell subprocess.
                process = subprocess.Popen(
                    cmd, cwd=MODULES_DIR, shell=True, stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, text=True
                )
                # Iterate over stdout from process andlog
                for line in process.stdout:
                    vcf_log.bash(line.strip())
                process.wait()
                self.log.note(f"VCF file generated for {current_line_name}.") 
                self.log.note(f"Log saved to {log_path}")
                self.log.note(f'VCF table path added to experiments_dictionary: {vcf_table_path}')
                
            self.log.success("VCF file generation process complete")
            return experiment_dictionary

        except Exception as e:
            self.log.fail(f"Error while generating the VCF file for {current_line_name}: {e}")


    def bsa_analysis(self, experiment_dictionary):
        '''
        Parent function for running BSA analysis. Given an experiment_dictionary
        containing the vcf_table_path and line_name, will output a list of 
        candidate mutations and 5 graphs, which will help narrow down regions of
        interest.  
        
        Dictionary Requirements:
        required 
        [key]: line_name
        [value] vcf_table_path, core_ulid, 

        optional
        [value] vcf_ulid
        
        vcf_table_path must be the path to a vcf table produced from VCFgen.sh. 
        support for VCF tables which are configured differently is on the roadmap. 
        
        '''
        self.log.attempt("Attempting to perform data analysis...")
        
        try:
            # Print experiment_dictionary information passed to function
            for key, value in experiment_dictionary.items():
                self.log.delimiter('Experiment Dictionary passed to bsa_analysis function')
                self.log.print(f"Key: {key}")
                for inner_key, inner_value in value.items():
                    self.log.print(f"{inner_key}: {inner_value}")
                self.log.print(' ')
                core_ulid = self.log.ulid
            
            # Run analysis
            for key, value in experiment_dictionary.items():

                core_ulid = self.log.ulid

                current_line_name = key
                vcf_ulid = value['vcf_ulid'] 
                vcf_table_path = value['vcf_table_path']
                
                # Configure an analysis logger for each line.
                analysis_log = LogHandler(f'analysis_{current_line_name}')

                self.log.note(f'Analysis log initialized. Path: {analysis_log.log_path}')
                analysis_log.add_db_record(current_line_name, core_ulid, vcf_ulid)
                
                #FileUtilites instance that logs to analysis_log
                file_utils = FileUtilities(analysis_log)

                #Analysis operations. Loading VCF and producing features
                vcf_df = file_utils.load_vcf_table(vcf_table_path, current_line_name)
                
                bsa_analysis_utils = BSAAnalysisUtilities(
                    current_line_name, vcf_ulid, analysis_log
                )
                vcf_df = bsa_analysis_utils.calculate_delta_snp_and_g_statistic(
                    vcf_df
                )
                vcf_df = bsa_analysis_utils.drop_na_and_indels(vcf_df) 
                vcf_df = bsa_analysis_utils.loess_smoothing(vcf_df)
                
                vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff = (
                    bsa_analysis_utils.calculate_empirical_cutoffs(vcf_df)
                )
    
                #Saving and plotting outputs
                bsa_analysis_utils.sort_save_likely_candidates(vcf_df)
                
                bsa_analysis_utils.generate_plots(
                    vcf_df, gs_cutoff, rsg_cutoff, rsg_y_cutoff
                )

            self.log.success("Data analysis complete")
        
        except Exception as e:
            self.log.fail(f"Error during data analysis: {e}")
