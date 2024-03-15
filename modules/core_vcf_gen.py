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
        
    def __call__(self):
        self.vcf_vars.gen_reference_chrs_paths()
        self._parse_reference_genome_path()
        self._create_chromosomal_fasta(self.vcf_vars.reference_genome_path, 
            self.vcf_vars.reference_chrs_fa_path, *self.vcf_vars.omit_chrs_patterns
        )

        for line in self.vcf_vars.lines:
            self.log.delimiter(f'Initializing vcf_generation subprocess log for {line.name}')
            try:
                vcf_log = self._initialize_vcf_log(line)
                self._generate_output_paths(line, vcf_log)
                self._run_vcfgen_subprocess(line, vcf_log)
                self._format_vcf_output(line, vcf_log)
                if self.vcf_vars.cleanup:
                    self._cleanup_files(line.vcf_output_dir, self.vcf_vars.cleanup_filetypes, vcf_log)
            except Exception as e:
                self.log.error(f"There was an error while running vcf_generation process for {line.name}:{e}.")
                self.log.error("Continuing to see if other lines can be run...")
                continue

    def _parse_reference_genome_path(self):
        file_utils = FileUtilities(self.log)
        self.vcf_vars.reference_genome_path = file_utils.parse_file(
            self.vcf_vars.reference_genome_path, 
            self.vcf_vars.reference_genome_source,
            REFERENCE_DIR
        )

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

    def _initialize_vcf_log(self, line):
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

        self.vcf_vars.log = vcf_log # redirect logging in vcf_vars to the vcf log

        return vcf_log

    def _generate_output_paths(self, line, vcf_log):
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

    def _run_vcfgen_subprocess(self, line, vcf_log):
        line.vcf_gen_cmd = self.vcf_vars.make_vcfgen_command(line)
        process = subprocess.Popen(
            line.vcf_gen_cmd, 
            cwd=MODULES_DIR, 
            shell=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True
        )
        vcf_log = LogHandler(f'vcf_{line.name}')
        for stdout_line in process.stdout:
            vcf_log.bash(stdout_line.strip())
        process.wait()
        vcf_log.note(f"VCF file generated for {line.name}.") 
        vcf_log.success("VCF file generation process complete")

    def _format_vcf_output(self, line, vcf_log):
        try:
            self.log.attempt("Attempting to format vcf table...")
            vcf_format = VCFFormat(line.snpsift_out_path, line.vcf_table_path, vcf_log)
            vcf_format.format_fields()
            vcf_format.remove_complex_genotypes()
        except Exception as e:
            self.log.error(f"There was an error while formatting vcf table:{e}")

    def _cleanup_files(self, output_dir_path: str, cleanup_filetypes: list):
        for file_type in cleanup_filetypes:
            if file_type == '*.table':
                self.log.warning("*.table not allowed in cleanup_filetypes. This is the point of running VCF_gen in the first place. Continuing...")
                continue
            files = glob.glob(os.path.join(output_dir_path, file_type))
            for file in files:
                if os.path.isfile(file):
                    os.remove(file)


class VCFFormat:
    
    def __init__(self, snpsift_out_path, vcf_table_path, logger):
        self.snpsift_out_path = snpsift_out_path
        self.vcf_table_path = vcf_table_path
        self.vcf_df = None
        self.log = logger

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

    def format_fields(self):
        self.log.attempt(f"Attempting to format vcf table:{self.vcf_table_path}")
        try:
            with open(self.snpsift_out_path, 'r') as f_in, open(self.vcf_table_path, 'w') as f_out:
                header_line = next(f_in).strip()
                headers = header_line.split('\t')
                output_headers = []
                for header in headers:
                    mapped_header = self.header_map.get(header, header)
                    if isinstance(mapped_header, list):
                        output_headers.extend(mapped_header)
                    else:
                        output_headers.append(mapped_header)
                f_out.write("\t".join(output_headers) + "\n") 

                for line in f_in:
                    fields = line.strip().split('\t') 
                    new_fields = []
                    for header, field in zip(headers, fields):
                        if "AD" in header and ',' in field:
                            new_fields.extend(field.split(','))
                        else:
                            new_fields.append(field)
                    f_out.write('\t'.join(new_fields) + "\n")

            self.log.success(f"VCF table formatted successfully!")

        except IOError as e:
            self.log.error(f"IOError occurred: {e}")
        
        except Exception as e:
            self.log.error(f"There was an error while formatting the feilds for vcf file:{e}")
    
    def remove_complex_genotypes(self):
        self.log.attempt("Removing complex genotypes (can be found in .complex_genos file)")
        try:
            complex_geno_tablename = f"{self.vcf_table_path}.complex_genos"
            with open(self.vcf_table_path, 'r') as f:
                lines = f.readlines()

            with open(self.vcf_table_path, 'w') as f, open(complex_geno_tablename, 'w') as f_too_complex:
                for line in lines:
                    if len(line.split('\t')) == 13:
                        f.write(line)
                    else:
                        f_too_complex.write(line)
        except Exception as e:
            self.log.error(F"Removing complex genotypes failed for an unknown reason: {e}")