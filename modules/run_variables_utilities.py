from settings.config import INPUT_DIR, OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR
import settings.vcf_gen_variables

import os
import pandas as pd
import re
import sqlite3

class Lines:
    def __init__(self, name):
        self.name = name
        self.segregation_type = None
        self.wt_input = []
        self.mu_input = []
        self.pairedness = None
        self.vcf_table_path = None
        self.output_dir_path
        self.output_prefix
        self.core_ulid = None
        self.vcf_ulid = None
        self.analysis_ulid = None

        self.in_path_variables=[
           'vcf_table_path',
        ]
        self.ref_path_variables=[
           'known_snps_path'
        ]
class RuntimeVariables:
    def __init__(self, logger):
        self.log = logger

        self.lines=[]        
        
        self.reference_genome_name = None
        self.reference_genome_source = None
        self.cleanup = None
        self.threads_limit = None
        self.known_snps_path = None
        self.call_variants_in_parallel = None
        self.snpEff_species_db = None

    def _process_file(self, filename):
        self.log.attempt(f'Parsing {filename}')
        for filename in os.listdir(INPUT_DIR):
            line_name, segregation_type, bulk_type, pairedness = self._parse_filename(filename)
            
            self.log.success(f"""{filename} parsed. 
            line_name:{line_name}
            segregation_type:{seg_type}
            bulk_type:{bulk_type}
            pairedness:{pairedness}
            """)

            line = self._get_line(line_name)
            
            line.segregation_type = segregation_type
            line.pairedness = pairedness
            file_path = os.path.join(INPUT_DIR, filename)
            
            self._append_file_path(bulk_type, line, file_path)

    def _parse_filename(self, filename):

        parts = filename.split('.')
        self.log.attempt(f'Parsing {filename}')
        line_name = parts[0]
        segregation_type = parts[1]

        if '_1' in segregation_type:
            segregation_type = segregation_type.rstrip('_1')
        if '_2' in segregation_type:
            segregation_type = segregation_type.rstrip('_2')

        bulk_type = parts[-3]
        pairedness = 'paired-end' if '_1' or '_2' in filename else 'single-read'

        return line_name, segregation_type, bulk_type, pairedness

    def _get_line(self, line_name):
        for l in self.lines:
            if l.name == line_name:
                return l
                
        line = Lines(line_name)
        self.lines.append(line)
        return line

    def _append_file_path(self, bulk_type, line, file_path):
        if bulk_type == 'wt':
            line.wt_input.append(file_path)
        elif bulk_type == 'mu':
            line.mu_input.append(file_path)

    def _sort_file_paths(self):
        for line in self.lines:
            sorted_wt_inputs = " ".join(sorted(line.wt_input))
            sorted_mu_inputs = " ".join(sorted(line.mu_input))
            line.wt_input = sorted_wt_inputs
            line.mu_input = sorted_mu_inputs

    def autodetect_line_variables(self):
        self.log.attempt(f"Detecting experiment details in: {INPUT_DIR}")
        try:
            for filename in os.listdir(INPUT_DIR):
                self._process_file(filename)

            self._sort_file_paths()

            self.log.success(f'Line and run details generated.')
        except Exception as e:
            self.log.fail(f"Error while detecting line and run details: {e}")

    def _process_input(self, key, details):
        file_utils = self.FileUtilities(self.log)
        if key in self.in_path_variables:
            return file_utils.process_path(INPUT_DIR, details)
        elif key in self.ref_path_variables:
            return file_utils.process_path(REFERENCE_DIR, details)
        else:
            return details

    def usr_in_line_variables(self, line_name, **kwargs):
        self.log.attempt('Attempting to organize inputs into Line data class...')
       
        try:
            line = Lines(line_name)
            for key, details in kwargs.items():
                processed_input = self._process_input(key, details)
                if hasattr(line, key):
                    setattr(line, key, processed_input)
                self.log.note(f'Arg to line dictionary| {key}:{details}')

            if 'vcf_table_path' in kwargs and 'vcf_ulid' not in kwargs:
                msg = 'vcf table path provided, vcf_ulid not. Extracting ulid...'
                self.log.note(msg)
                file_utils = self.FileUtilities(self.log)
                line.vcf_ulid = file_utils.extract_ulid_from_file_path(
                                 kwargs['vcf_table_path'])

            self.lines.append(line)
            self.log.success('Inputs successfully organized into line_dict.')
            
        except Exception as e:
            self.log.fail(f'Error creating line dictionary from inputs {e}.') 

    def _generate_output_dir_path(self, line):
        file_utils = FileUtilities(self.log)
        try:
            output_name_prefix = f"{self.log.ulid}_-{line.name}"
            line.output_dir_path = os.path.join(OUTPUT_DIR, output_name_prefix)
            file_utils.setup_directory(line.output_dir_path)
            self.log.success('Output directory path set.')
        except Exception as e:
            self.log.fail(f'Error producing output directory path: {e}')

    def _generate_output_prefix(self, line):
        try:
            line.output_prefix = os.path.join(line.output_dir_path, line.output_name_prefix)
            self.log.success('Output prefix set.')
        except Exception as e:
            self.log.fail(f'Error producing output prefix: {e}')

    def _check_and_add_vcf_table_path(self, line):
        try:
            if line.vcf_table_path is None:
                vcf_table_name = f"{line.output_prefix}.noknownsnps.table"
                line.vcf_table_path = os.path.join(OUTPUT_DIR, vcf_table_name)
                self.log.success('VCF table path set.')
            else:
                self.log.success('VCF table path already set, proceeding.')
        except Exception as e:
            self.log.fail(f'Error setting vcf_table_path: {e}')

    def generate_output_file_paths(self):
        for line in self.lines:
            self._generate_output_dir_path(line)
            self._generate_output_prefix(line)
            self._check_and_add_vcf_table_path(line)

    def apply_settings(self):
        for attribute in vars(self).keys():
            if getattr(self, attribute) is None and hasattr(vcf_gen_variables, attribute):
                setattr(self, attribute, getattr(vcf_gen_variables, attribute))
