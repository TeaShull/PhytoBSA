from settings.config import INPUT_DIR, REFERENCE_DIR, MODULES_DIR, VCF_GEN_SCRIPT

import subprocess

from modules.utilities_bsa_analysis import BSAAnalysisUtilities
from modules.utilities_general import FileUtilities
from modules.utilities_logging import LogHandler

class ThaleBSAParentFunctions:
    def __init__(self, logger):
        self.log = logger

    def populate_experiment_dictionary(self, line, **kwargs)-> dict:
        '''
        Organizes user inputs into a dictionary. 
        If the input needs to be a path, it checks if it exists. If it doesn't
        the function will check INPUT_DIR for the file and correct the pathing.

        input: line (which is the primary key for the dict, and is required)
        and user arguments passed as variables.

        output: experiment_dictionary, which organizes all experiment details
        under the line name as the primary key.

        Approved inputs and path variables need to be organized into the lists
        at the top of the function. Any time changes to inputs to the major
        functions of this program are made, organize them into this dictionary. 
        This way, all functions can flexably be assigned user inputs using an 
        easy to pass object. 
        '''
        approved_inputs = [
            'segregation_type',
            'pairedness',
            'reference_genome_name',
            'reference_genome_source',
            'cleanup',
            'vcf_table_path',
            'threads_limit',
            'wt_input',
            'mu_input',
            'known_snps_path',
            'call_variants_in_parallel',
            'snpEff_species_db',
            'core_ulid',
            'vcf_ulid',
            'analysis_ulid'
        ]
        in_path_variables =[
            'vcf_table_path',
        ]

        ref_path_variables =[
            'known_snps_path'
        ]
        details = {}

        self.log.attempt(f'Attempting to organize user inputs into an experiment dictionary')
        file_utils = FileUtilities(self.log)
        try:
            if not line:
                self.log.fail(f'line name can not be empty - it is required to create an experiment dictionary! Aborting')
            for key, details in kwargs.items():
                if key in approved_inputs and in_path_variables:
                    path = fileutils.process_path(INPUT_DIR, details[key])
                    details[key] = path
                elif key in approved_inputs and ref_path_variables:
                    path = fileutils.process_path(REFERENCE_DIR, details[key])
                elif key in approved_inputs:
                    details[key] = details
                    self.log.note(f'Arg to experiment dictionary| {key}:{details[key]}')
                else: 
                    self.log.fail(f"{key} is not in the approved_inputs for experiment dictionary creation. Dictionary or arguments may need updating.")

            if 'vcf_table_path' in kwargs.items() and 'vcf_ulid' not in kwargs.items():
                self.log.note('vcf table path provided, but vcf_ulid was not. Trying to extract ulid from file name...')
                details['vcf_ulid'] = file_utils.extract_ulid_from_file_path(vcf_table_path)

            experiment_dict = {
                line : details
            }
            self.log.success('User inputs successfully organized into experiment_dictionary.')
            return experiment_dict
            
        except Exception as e:
            self.log.fail(f'Error creating the experiment dictionary from user inputs {e}.')

    def vcf_generation(self, experiment_dictionary=None)->dict:
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

        # Generate the knownSnps .vcf file path
        self.log.note('Beginning VCF generation process for experiment_dictionary')
        try:
        experiment_dictionary = file_utils.generate_output_file_paths(experiment_dictionary)
            for line, details in experiment_dictionary.items():
                # Print variables contained in experiment_dictionary to core_log
                self.log.delimiter('Experiment Dictionary passed to vcf_generation function')
                self.log.print(f'Line:{line}')
                for details_key, value in details.items():
                    self.log.print(f"{details_key}: {value}")
                self.log.print(' ')
                '''
                Unpack required variables for [line]
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
                
                #Generate file paths needed for vcf generation
                file_utils = FileUtilities(vcf_log)
                # Construct args to pass variables to VCFgen.sh.
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
                print(vcfgen_script_path)
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
            return experiment_dictionary

        except Exception as e:
            self.log.fail(f"Error while generating the VCF file for {line}: {e}")

    def bsa_analysis(self, experiment_dictionary):
        '''
        Parent function for running BSA analysis. Given an experiment_dictionary
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
            for line, details in experiment_dictionary.items():
                # Print variables contained in experiment_dictionary to core_log
                self.log.delimiter('Experiment Dictionary passed to bsa_analysis function:')
                self.log.print(f'Line:{line}')
                for details_key, value in details.items():
                    self.log.print(f"   {details_key}: {value}")
                self.log.print(' ')
                '''
                Unpack required variables to run analysis
                '''
                vcf_table_path = details['vcf_table_path']
                vcf_ulid = details['vcf_ulid']
                segregation_type = details['segregation_type']


                # Configure an analysis logger for each line.
                analysis_log = LogHandler(f'analysis_{line}')
                self.log.note(f'Analysis log initialized. Path: {analysis_log.log_path}')
                analysis_log.add_db_record(
                    line, self.ulid, vcf_ulid
                )

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
            
                # Push vcf_df out of the class, for posterity
                vcf_df = basa_analysis_utils.vcf_df
                
                ## Saving and plotting outputs
                bsa_analysis_utils.sort_save_likely_candidates()
                bsa_analysis_utils.generate_plots()
            self.log.success("Data analysis complete")
        
        except Exception as e:
            self.log.fail(f"Error during data analysis: {e}")
