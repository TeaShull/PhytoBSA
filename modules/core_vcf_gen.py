import subprocess

from modules.utilities_logging import LogHandler

class VCFGenerator
    def __init__(self, vcf_gen_vars, logger):
        self.log = logger
        self.vcf_gen_vars = vcf_gen_vars

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
        self.log.note('Beginning VCF generation process for line_dict')
        

        for line in self.vcf_gen_vars.lines
            self.log = LogHandler(f'vcf_{line}')
            line.vcf_ulid = self.log.ulid
            vcf_log.add_db_record(line, line.vcf_ulid)

            #Generate output paths
            vcf_out_paths = gen_vcf_output_paths(name, line.vcf_ulid)
            (
                line.vcf_output_dir_path, 
                line.vcf_output_prefix, 
                line.vcf_table_path
            ) = vcf_out_paths

            #Generate line.vcf_gen_cmd
            line.vcf_gen_cmd = make_vcfgen_command(line)
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

            self.log.note(f"VCF file generated for {line}.") 
            self.log.success("VCF file generation process complete")

        except Exception as e:
            self.log.fail(f"Error while generating the VCF file for {line}: {e}")