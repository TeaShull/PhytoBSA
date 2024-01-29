from settings.config import (INPUT_DIR, OUTPUT_DIR, MODULES_DIR, REFERENCE_DIR, 
     VCF_GEN_SCRIPT)
from settings import vcf_gen_variables
from modules.utilities_general import FileUtilities

import os
import pandas as pd
class Lines:
    '''
    This class contains information and methods pertaining to variables
    needed for both vcf generation and bsa pipeline. The variables contained 
    here are dynamic, and many are used in both processes. I attempt to make
    clear what variables are needed for each process.  

    The variables of this class can be populated by automatic_line_variables
    in the AutomaticLineVariableDetector class or the usr_in_line_variables 
    within this current class. 
    '''
    
    __slots__ = [
        'log', 'name', 'segregation_type', 'vcf_table_path', 
        'mu_input', 'wt_input', 'pairedness', 'vcf_gen_cmd', 
        'vcf_output_prefix', 'vcf_output_dir', 'vcf_ulid', 'vcf_df', 
        'analysis_out_prefix', 'analysis_out_path', 'gs_cutoff', 
        'rsg_cutoff', 'rsg_y_cutoff', 'analysis_ulid', 
        'in_path_variables', 'ref_path_variables', 'snpeff_dir', 
        'snpeff_out_filename'
    ]

    def __init__(self, name, logger):
        self.log = logger

        # BSA and VCF_gen variables
        self.name = name
        self.vcf_table_path = None 

        # VCF_gen variables 
        self.mu_input = []
        self.wt_input = []
        self.pairedness = None
        self.vcf_gen_cmd = None
        self.vcf_output_prefix = None
        self.vcf_output_dir = None
        self.vcf_ulid = None
        self.snpeff_dir = None
        self.snpeff_out_filename = None


        # BSA variables
        self.segregation_type = None
        self.vcf_df = None
        self.analysis_out_prefix = None
        self.analysis_out_path = None
        self.gs_cutoff = None
        self.rsg_cutoff = None
        self.rsg_y_cutoff = None
        self.analysis_ulid = None

        '''
        Variables that require path parsing. Variables in in_path_variables
        will be double checked - first if they are hard-coded, or if they 
        are in ./data/inputs. Same deal with ref_path_variables, but for 
        the ./data/references folder.
        '''  
        self.in_path_variables=[
           'vcf_table_path',
           'mu_input',
           'wt_input'
        ]
        self.ref_path_variables=[
           'known_snps_path'
        ]

    def _process_input(self, key, details):
        file_utils = self.FileUtilities(self.log)
        if key in self.in_path_variables or key in self.ref_path_variables:
            paths = details.split()
            dir = INPUT_DIR if key in self.in_path_variables else REFERENCE_DIR
            processed_paths = [file_utils.process_path(dir, path) for path in paths]
            return ' '.join(processed_paths)
            return details
    
    def usr_in_line_variables(self, **kwargs):
        self.log.attempt('Attempting to organize inputs into Line data class...')
        
        try:
            self.lines = Lines(name)
            for key, details in kwargs.items():
                try:
                    processed_input = self._process_input(key, details)
                except KeyError as ke:
                    self.log.fail(f'Key not found: {ke}')
                    continue
                 
                try:
                    if hasattr(line, key):
                        setattr(line, key, processed_input)
                except AttributeError as ae:
                    self.log.fail(f'Attribute error: {ae}')
        
                self.log.note(f'User input added: |{name}|{key}:{processed_input} ')
        except Exception as e:
            self.log.fail(f'There was an error processing user inputs:{e}')

class AutomaticLineVariableDetector:
    def __init__(self, logger):
        self.log = logger
        self.lines = []

    def _parse_filename(self, filename):
        parts = filename.split('.')
        self.log.attempt(f'Parsing {filename}')
        name = parts[0]
        segregation_type = parts[1]

        if '_1' in segregation_type:
            segregation_type = segregation_type.rstrip('_1')
        if '_2' in segregation_type:
            segregation_type = segregation_type.rstrip('_2')

        bulk_type = parts[-3]
        pairedness = 'paired-end' if '_1' or '_2' in filename else 'single-read'

        return name, segregation_type, bulk_type, pairedness

    def _get_line(self, name):
        for l in self.lines:
            if l.name == name:
                return l
                
        line = Lines(name, self.log)
        self.lines.append(line)
        return line

    def _append_file_path(self, bulk_type, line, file_path):
        if bulk_type == 'wt':
            line.wt_input.append(file_path)
        elif bulk_type == 'mu':
            line.mu_input.append(file_path)

    def _process_file(self, filename):
        self.log.attempt(f'Parsing {filename}')
        details = self._parse_filename(filename)
        name, segregation_type, bulk_type, pairedness = details
        
        self.log.success(f"""{filename} parsed. 
        name:{name}
        segregation_type:{segregation_type}
        bulk_type:{bulk_type}
        pairedness:{pairedness}
        """)

        line = self._get_line(name)
        
        line.segregation_type = segregation_type
        line.pairedness = pairedness
        file_path = os.path.join(INPUT_DIR, filename)
        
        self._append_file_path(bulk_type, line, file_path)
    
    def _sort_file_paths(self):
        for line in self.lines:
            sorted_wt_inputs = "' '".join(sorted(line.wt_input))
            sorted_mu_inputs = "' '".join(sorted(line.mu_input))
            line.wt_input = sorted_wt_inputs
            line.mu_input = sorted_mu_inputs

    def automatic_line_variables(self):
        self.log.attempt(f"Detecting experiment details in: {INPUT_DIR}")
        try:
            for filename in os.listdir(INPUT_DIR):
                self._process_file(filename)

            self._sort_file_paths()

            self.log.success(f'Line and run details generated.')
        except Exception as e:
            self.log.fail(f"Error while detecting line and run details: {e}")

class VCFGenVariables:
    def __init__(self, logger,
                 lines, 
                 reference_genome_name=None, 
                 snpEff_species_db=None,
                 reference_genome_source=None, 
                 known_snps_path=None, 
                 threads_limit=None, 
                 call_variants_in_parallel=None, 
                 cleanup=None 
                 ):
        
        self.log = logger
        self.lines=lines

        self.reference_genome_name = reference_genome_name
        self.reference_genome_source = reference_genome_source
        self.cleanup = cleanup
        self.threads_limit = threads_limit
        self.known_snps_path = known_snps_path
        self.call_variants_in_parallel = call_variants_in_parallel
        self.snpEff_species_db = snpEff_species_db

    def make_vcfgen_command(self, line):
        self._apply_settings() #if None, source from settings.vcf_gen_variables
        rgp = os.path.join(REFERENCE_DIR, self.reference_genome_name)
        reference_genome_prefix = rgp
        args = (
            line.vcf_ulid,
            line.name, 
            INPUT_DIR,
            line.wt_input,
            line.mu_input,
            line.pairedness,
            line.vcf_output_dir,
            line.vcf_output_prefix,
            line.vcf_table_path,
            line.snpeff_dir,
            line.snpeff_out_filename,
            REFERENCE_DIR,
            self.reference_genome_source, 
            reference_genome_prefix, 
            self.snpEff_species_db,
            self.known_snps_path,
            self.threads_limit,
            self.call_variants_in_parallel,
            self.cleanup
        )
        cmd = f"{VCF_GEN_SCRIPT} {' '.join(map(str, args))}"
        return cmd

    def _apply_settings(self):
        for attribute in vars(self).keys():
            if getattr(self, attribute) is None and hasattr(vcf_gen_variables, attribute):
                setattr(self, attribute, getattr(vcf_gen_variables, attribute))

    def gen_vcf_output_paths(self, name, vcf_ulid):
        output_name = f"{vcf_ulid}-_{name}"
        vcf_output_dir_path = os.path.join(OUTPUT_DIR, output_name)
        vcf_output_prefix = os.path.join(vcf_output_dir_path, output_name) 
        vcf_table_path= f"{vcf_output_prefix}.noknownsnps.table"
        
        snpeff_dir= os.path.join(OUTPUT_DIR, 'snpEff')
        snpeff_out_filename=os.path.join(snpeff_dir, output_name)
        
        return (vcf_output_dir_path, vcf_output_prefix, 
            vcf_table_path, snpeff_dir, snpeff_out_filename
        )

class BSAVariables:
    def __init__(self, logger,
             lines, 
             loess_span, 
             smooth_edges_bounds, 
             shuffle_iterations):
        
        self.log = logger
        self.lines = lines

        self.loess_span = loess_span
        self.smooth_edges_bounds = smooth_edges_bounds 
        self.shuffle_iterations = shuffle_iterations

    def load_vcf_table(self, vcf_table_path)->pd.DataFrame:
        """
        Loads VCF table into a pandas dataframe.
        
        Args:  
        current_line_table_path(str) - path to the vcf table to be loaded into df
        current_line_name(str) - name of the line associated with the vcf table

        Returns: 
        Pandas dataframe containing the information loaded from current_line_table_path
        """
        self.log.attempt(f"Attempting to load VCF table for line {current_line_name}")
        try:
            vcf_df = pd.read_csv(current_line_table_path, sep="\t")
            self.log.attempt(f"The VCF table for line {current_line_name} was successfully loaded.")
            return vcf_df
        
        except FileNotFoundError:
            self.log.fail(f"Error: File '{current_line_table_path}' not found.")
        
        except pd.errors.EmptyDataError:
            self.log.fail(f"Error: File '{current_line_table_path}' is empty.")
        
        except Exception as e:
            self.log.fail(f"An unexpected error occurred: {e}")

    def gen_bsa_out_prefix(self, name, ulid, vcf_ulid): #called in core_bsa.py
        analysis_out_prefix = f'{ulid}_-{name}'
        if vcf_ulid:
            out_path = os.path.join(
                OUTPUT_DIR, 
                f'{vcf_ulid}_-{name}',
                analysis_out_prefix
            )
        else:
            out_path = os.path.join(
                OUTPUT_DIR,
                analysis_out_prefix
            )
        file_utils = FileUtilities(self.log)
        file_utils.setup_directory(out_path)
        return out_path