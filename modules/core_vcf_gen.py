import os
import subprocess
import re

from Bio import SeqIO
import glob
import pandas as pd

from modules.utilities_general import FileUtilities
from modules.utilities_logging import LogHandler
from settings.globals import MODULES_DIR, REFERENCE_DIR, THREADS_LIMIT

"""
Core module for VCF generation subprocess.
Input variable class: vcf_vars from utilities_variables module.

This module processes bulk segregant fasta files, and formats them for BSA. 
"""


class VCFGenerator:
    
    def __init__(self, logger, vcf_vars):
        self.vcf_vars = vcf_vars
        self.log = logger #pass core_log from main phytobsa script

    def _create_chromosomal_fasta(self, input_file: str, output_file: str, *patterns: list):
        if not os.path.isfile(output_file):
            self.log.attempt(f"Creating {output_file} with only chromosomal DNA...")
            self.log.note(f"Omitting contigs containing {', '.join(patterns)}")
            with open(input_file, "r") as in_handle, open(output_file, "w") as out_handle:
                for record in SeqIO.parse(in_handle, "fasta"):
                    if not any(re.search(pattern, record.description) 
                        for pattern in patterns
                    ):
                        SeqIO.write(record, out_handle, "fasta")
        else:
            self.log.note(f"Chromosomal fasta exists: {output_file}")
            self.log.note("Proceeding...")
    
    def _cleanup_files(self, output_dir_path: str, cleanup_filetypes: list):
        for file_type in cleanup_filetypes:
            if file_type == '*.table':
                self.log.warning("*.table not allowed in cleanup_filetypes. This is the point of running VCF_gen in the first place. Continuing...")
                continue
            files = glob.glob(os.path.join(output_dir_path, file_type))
            for file in files:
                if os.path.isfile(file):
                    os.remove(file)

    def __call__(self):
        """
        Input: VCFGenVariables class instance, populated with the variables needed
        to run the process. 
        
        Required variables can be seen in core_variables.VCFVariables

        the make_vcfgen_command returns the cmd passed to subprocess_VCFgen.sh 
        This command needs to be kept in order, as the variable parsing in the
        shell script is ordered.
        """
        
        #Generate chrs paths
        self.vcf_vars.gen_reference_chrs_paths()
        
        # Parse reference genome path. If exists - proceed. 
        # If compressed, decompress. 
        # If not there, download from reference_genome_source and decompress
        file_utils = FileUtilities(self.log)
        self.vcf_vars.reference_genome_path = file_utils.parse_file(
            self.vcf_vars.reference_genome_path, 
            self.vcf_vars.reference_genome_source,
            REFERENCE_DIR
        )
        # Create a new fasta file, with unneeded contigs filtered based on 
        # user defaults. (mitochondria, exc)
        self._create_chromosomal_fasta(self.vcf_vars.reference_genome_path, 
            self.vcf_vars.reference_chrs_fa_path, *self.vcf_vars.omit_chrs_patterns
        )

        for line in self.vcf_vars.lines:
            self.log.delimiter(f'Initializing vcf_generation subprocess log for {line.name}')
            try:
                #Initialize vcf log for the current line name
                vcf_log = LogHandler(f'vcf_{line.name}')
                line.vcf_ulid = vcf_log.ulid
                vcf_log.add_db_record(
                    name=line.name, 
                    core_ulid=self.log.ulid,
                    reference_genome_path=self.vcf_vars.reference_genome_path,
                    snpeff_species_db=self.vcf_vars.snpeff_species_db,
                    reference_genome_source=self.vcf_vars.reference_genome_source,
                    omit_chrs_patterns=self.vcf_vars.omit_chrs_patterns,
                    threads_limit=THREADS_LIMIT
                )
                
                #Generate output paths for process
                vcf_out_paths = self.vcf_vars.gen_vcf_output_paths(
                    line.name, line.vcf_ulid
                )
                (
                    line.vcf_output_dir, 
                    line.vcf_output_prefix, 
                    line.vcf_table_path,
                    line.snpeff_report_path,
                    line.snpeff_out_path, 
                    line.snpsift_out_path
                ) = vcf_out_paths

                #Generate line.vcf_gen_cmd
                line.vcf_gen_cmd = self.vcf_vars.make_vcfgen_command(line)

                #Run vcfgen shell subprocess.
                process = subprocess.Popen(
                    line.vcf_gen_cmd, 
                    cwd=MODULES_DIR, 
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, 
                    text=True
                )
                # Iterate over stdout from process and log
                for stdout_line in process.stdout:
                    vcf_log.bash(stdout_line.strip())
                process.wait()

                vcf_log.note(f"VCF file generated for {line.name}.") 
                vcf_log.success("VCF file generation process complete")
                
                # Formatting VCF output to make it dataframe friendly
                vcf_format = VCFFormat(line.snpsift_out_path, line.vcf_table_path)
                vcf_format.format_fields()
                vcf_format.remove_complex_genotypes()
                
                #Cleanup files if cleanup = True
                if self.vcf_vars.cleanup:
                    self._cleanup_files(line.vcf_output_dir, self.vcf_vars.cleanup_filetypes)
            
            except Exception as e:
                self.log.error(f"There was an error while running vcf_generation process for {line.name}. Continuing to see if other lines can be run...")
                continue
            

class VCFFormat:
    
    def __init__(self, snpsift_out_path, vcf_table_path):
        self.snpsift_out_path = snpsift_out_path
        self.vcf_table_path = vcf_table_path
        self.vcf_df = None

        self.header_map = {
            "CHROM": "chrom",
            "POS": "pos",
            "REF": "ref",
            "ALT": "alt",
            "ANN[*].GENE": "gene",
            "ANN[*].EFFECT": "snpEffect",
            "ANN[*].HGVS_P": "snpVariant",
            "ANN[*].IMPACT": "snpImpact",
            "GEN[*].GT": "mu:wt_GTpred",
            "GEN[mu].AD": ["mu_ref", "mu_alt"],
            "GEN[wt].AD": ["wt_ref", "wt_alt"]
        }

    def format_fields(self, snpsift_out_path, vcf_table_path):
        with open(snpsift_out_path, 'r') as f_in, open(vcf_table_path, 'w') as f_out:
            header_line = next(f_in) #grab header line
            headers = header_line.split('\t') # articulate each header
            output_headers = [self.header_map.get(header, header) for header in headers] #map each header to self.header_map to create output_headers
            f_out.write("\t".join(output_headers) + "\n") #write the output headers to file

            for line in f_in:
                fields = line.split('\t') 
                new_fields = [] #new fields list
                #create header, field tuples and interate through, splitting the allele depth columns (AD)
                for header, field in zip(headers, fields): #headers is static - from "articulate each header", fields from for line in f_in. creates dynamic tuples of input headers and feilds
                    if "AD" in header and ',' in field:
                        new_fields.extend(field.split(','))
                    else:
                        new_fields.append(field)
                f_out.write('\t'.join(new_fields) + "\n")

    def remove_complex_genotypes(self):
        complex_geno_tablename = f"{self.vcf_table_path}.complex_genos"
        with open(self.vcf_table_path, 'r') as f:
            lines = f.readlines()

        with open(self.vcf_table_path, 'w') as f, open(complex_geno_tablename, 'w') as f_too_complex:
            for line in lines:
                if len(line.split('\t')) == 13:
                    f.write(line)
                else:
                    f_too_complex.write(line)
    