import subprocess

from modules.utilities_logging import LogHandler
from settings.config import MODULES_DIR

class VCFGenerator:
    def __init__(self, logger, vcf_gen_vars):
        self.vcf_gen_vars = vcf_gen_vars
        self.log = logger

    def run_subprocess(self):
        """
        Input: Subprocess VCFgen.sh takes raw reads(wild-type(wt) 
        and mutant(mu)), either single or paired end and generates VCF table 
        *.noknownsnps.table.

        Nonetypes - if nothing is provided, line_dict will be 
        detected automatically from ./input and variables will be sourced from 
        the variables.py module.
        
        Returns: updated line_dict containing the paths to 
        the noknownsnps.table.
        """
        # Generate the knownSnps .vcf file path
        for line in self.vcf_gen_vars.lines:
            self.log.note(f'Initializing vcf_generation subprocess log for {line.name}')
            vcf_log = LogHandler(f'vcf_{line.name}')
            line.vcf_ulid = vcf_log.ulid
            vcf_log.add_db_record(line.name, line.vcf_ulid)

            #Generate output paths
            vcf_out_paths = self.vcf_gen_vars.gen_vcf_output_paths(
                line.name, line.vcf_ulid
            )
            (
                line.vcf_output_dir, 
                line.vcf_output_prefix, 
                line.vcf_table_path
            ) = vcf_out_paths

            #Generate line.vcf_gen_cmd
            line.vcf_gen_cmd = self.vcf_gen_vars.make_vcfgen_command(line)
            # Run vcfgen shell subprocess.
            process = subprocess.Popen(
                line.vcf_gen_cmd, 
                cwd=MODULES_DIR, 
                shell=True, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                text=True
            )
            # Iterate over stdout from process and log
            for line in process.stdout:
                vcf_log.bash(line.strip())
            process.wait()

            vcf_log.note(f"VCF file generated for {line}.") 
            vcf_log.success("VCF file generation process complete")