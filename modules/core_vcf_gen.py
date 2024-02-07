import os
import subprocess
import urllib.request
from urllib.parse import urlparse
import gzip
import shutil
import re

from Bio import SeqIO
import pandas as pd

from modules.utilities_logging import LogHandler
from settings.paths import MODULES_DIR, REFERENCE_DIR

"""
Core module for VCF generation subprocess.
Input variable class: vcf_vars from utilities_variables module.

This module processes bulk segregant fasta files, and formats them for BSA. 
"""


class VCFGenerator:
    def __init__(self, logger, vcf_vars):
        self.vcf_vars = vcf_vars
        self.log = logger #pass core_log from main phytobsa script


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
        
        # Parse reference genome path. If exists - proceed. 
        # If compressed, decompress. 
        # If not there, download from reference_genome_source and decompress
        self.vcf_vars.reference_genome_path = self._parse_reference_genome(
            self.vcf_vars.reference_genome_path, 
            self.vcf_vars.reference_genome_source
        )
        # Create a new fasta file, with unneeded contigs filtered based on 
        # user defaults. (mitochondria, exc)
        self._create_chromosomal_fasta(self.vcf_vars.reference_genome_path, 
            self.vcf_vars.reference_chrs_fa_path, *self.vcf_vars.omit_chrs_patterns
        )

        for line in self.vcf_vars.lines:
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
            
            
            #Initialize vcf log for the current line name
            self.log.delimiter(f'Initializing vcf_generation subprocess log for {line.name}')
            vcf_log = LogHandler(f'vcf_{line.name}')
            line.vcf_ulid = vcf_log.ulid
            vcf_log.add_db_record(line.name, line.vcf_ulid)

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
            vcf_format.remove_repetitive_nan() #vcf_table_path file created here
            vcf_format.format_fields()
            vcf_format.remove_complex_genotypes()
            vcf_format.add_header() #pd.DataFrame created
            vcf_format.remove_repetitive_snpeff_labels()#pd.DataFrame saved as csv 
            
            #Cleanup files if cleanup = True
            if self.vcf_vars.cleanup:
                self._cleanup_files(line.vcf_output_dir)
            

    def _parse_reference_genome(self, ref_genome_path: str, ref_genome_source: str)-> str:
        self.log.attempt("Attempting to parse reference genome...")
        try:
            if not os.path.isfile(ref_genome_path) and ref_genome_source:
                self.log.attempt(f"Reference genome doesn't exist, sourcing from: {ref_genome_source}")
                parsed_url = urlparse(ref_genome_source)
                source_extension = os.path.splitext(parsed_url.path)[1]
                
                if ref_genome_path.endswith(source_extension) is False:
                    self.log.note(f"Source extension {source_extension} and reference genome name extension don't match. Attempting to fix...")
                    ref_genome_path += source_extension
                
                self.log.attempt(f"Downloading from URL: {ref_genome_source}")
                urllib.request.urlretrieve(ref_genome_source, ref_genome_path)
                self.log.success(f"Reference genome downloaded and saved at {ref_genome_path}")
            
            if os.path.isfile(ref_genome_path) and ref_genome_path.endswith('.gz'):
                self.log.attempt(f"Attempting to unzip {ref_genome_path}")
                with gzip.open(ref_genome_path, 'rb') as f_in:
                    with open(ref_genome_path[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                ref_genome_path = ref_genome_path[:-3]
                self.log.success(f"File unzipped to {ref_genome_path}")

            return ref_genome_path
        
        except Exception as e:
            self.log.error("Parsing reference genome failed.")
            self.log.error(e)
            return None


    def _create_chromosomal_fasta(self, input_file: str, output_file: str, *patterns: list):
        if not os.path.isfile(output_file):
            self.log.attempt(f"Creating {output_file} with only chromosomal DNA...")
            print(f"Omitting contigs containing {', '.join(patterns)}")
            with open(input_file, "r") as in_handle, open(output_file, "w") as out_handle:
                for record in SeqIO.parse(in_handle, "fasta"):
                    if not any(re.search(pattern, record.id) for pattern in patterns):
                        SeqIO.write(record, out_handle, "fasta")
        else:
            self.log.note(f"Chromosomal fasta exists: {output_file}")
            self.log.note("Proceeding...")
    

    def _cleanup_files(self, output_dir_path: str, cleanup_filetypes: bool):
        for file_type in cleanup_filetypes:
            if file_type == '*.table':
                print("*.table not allowed in cleanup_filetypes. This is the point of running VCF_gen in the first place. Continuing...")
                continue
            files = glob.glob(os.path.join(output_dir_path, file_type))
            for file in files:
                if os.path.isfile(file):
                    os.remove(file)


class VCFFormat:
    def __init__(self, snpsift_out_path, vcf_table_path):
        self.snpsift_out_path = snpsift_out_path
        self.vcf_table_path = vcf_table_path
        self.vcf_df = None
    

    def remove_repetitive_nan(self):
        with open(self.snpsift_out_path, 'r') as f_in, open(self.vcf_table_path, 'w') as f_out:
            for line in f_in:
                f_out.write(line.replace('NaN:', ''))


    def format_fields(self):
        with open(self.vcf_table_path, 'r') as f:
            lines = f.readlines()

        with open(self.vcf_table_path, 'w') as f:
            for line in lines:
                fields = line.split('\t')
                for i in [8, 9, 10]:
                    fields[i] = fields[i].replace(',', '\t')
                f.write('\t'.join(fields))


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


    def add_header(self):
        self.vcf_df = pd.read_csv(self.vcf_table_path, sep='\t', header=None)

        # Add the header
        self.vcf_df.columns = [
            'chrom', 
            'pos', 
            'ref', 
            'alt', 
            'gene', 
            'snpEffect', 
            'snpVariant', 
            'snpImpact', 
            'mu:wt_GTpred', 
            'mu_ref', 
            'mu_alt', 
            'wt_ref', 
            'wt_alt'
        ]
    
    
    def remove_repetitive_snpeff_labels(self):

        columns_to_process = ["gene", "snpEffect", "snpVariant", "snpImpact"]
        # Iterate over snpEff columns and use "set" to get unique labels only
        for column in columns_to_process:
            self.vcf_df[column] = (
                self.vcf_df[column].astype(str)
                                .str.split(":")
                                .apply(lambda x: ":".join(set(x)))
            )

        # Write the DataFrame back to the file
        self.vcf_df.to_csv(self.vcf_table_path, sep='\t', index=False)