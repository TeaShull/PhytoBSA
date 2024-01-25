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