import os
import pandas as pd
from pathlib import Path

from settings.paths import (INPUT_DIR, OUTPUT_DIR, REFERENCE_DIR, MODULES_DIR, 
    VCF_GEN_SCRIPT
)
from modules.utilities_general import FileUtilities

class Lines:
    '''
    This class contains information and methods pertaining to variables
    needed for both vcf generation and bsa pipeline. The variables contained 
    here are dynamic, and many are used in both processes.

    The variables of this class can be populated by automatic_line_variables
    in the AutomaticLineVariableDetector class or the usr_in_line_variables 
    within this current class. 
    '''
    
    __slots__ = ['log', 'name', 'segregation_type', 'vcf_table_path', 
        'mu_input', 'wt_input', 'pairedness', 'vcf_gen_cmd', 
        'vcf_output_prefix', 'vcf_output_dir', 'vcf_ulid', 'vcf_df', 
        'analysis_out_prefix', 'analysis_out_path', 'gs_cutoff', 
        'rsg_cutoff', 'rsg_y_cutoff', 'analysis_ulid', 
        'in_path_variables', 'ref_path_variables', 'snpeff_dir', 
        'snpeff_out_filename', 'snpeff_out_path', 'snpsift_out_path'
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
        self.snpeff_out_path = None
        self.snpsift_out_path = None


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

    def _process_input(self, key, details):
        file_utils = self.FileUtilities(self.log)
        if key in self.in_path_variables or key in self.ref_path_variables:
            paths = details.split()
            dir = INPUT_DIR if key in self.in_path_variables else REFERENCE_DIR
            processed_paths = [file_utils.process_path(path, dir) for path in paths]
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
    __slots__ = ['log', 'lines', 'reference_genome_path', 
        'reference_genome_source', 'omit_chrs_patterns', 'snpEff_species_db', 
        'threads_limit', 'call_variants_in_parallel', 'cleanup', 
        'cleanup_filetypes', 'reference_chrs_fa_path', 
        'reference_chrs_dict_path'
    ]
    
    def __init__(self, logger,
        lines, 
        reference_genome_path, 
        reference_genome_source,
        omit_chrs_patterns,
        snpEff_species_db,
        threads_limit, 
        call_variants_in_parallel, 
        cleanup, 
        cleanup_filetypes
        ):
        self.log = logger
        
        self.lines = lines # can be single Line instance or list to iterate 
        self.reference_genome_path = reference_genome_path
        self.reference_genome_source = reference_genome_source
        self.omit_chrs_patterns = omit_chrs_patterns
        self.snpEff_species_db = snpEff_species_db
        self.threads_limit = threads_limit
        self.call_variants_in_parallel = call_variants_in_parallel
        self.cleanup = cleanup
        self.cleanup_filetypes = cleanup_filetypes
        
        # Generated by _gen_reference_chrs_paths
        self.reference_chrs_fa_path = None 
        self.reference_chrs_dict_path = None
        #Generate chrs paths
        self._gen_reference_chrs_paths()
    
    def make_vcfgen_command(self, line):
        args = (
            line.vcf_ulid, #
            line.name, 
            line.wt_input,
            line.mu_input,
            line.pairedness,
            line.vcf_output_prefix,
            line.snpeff_dir,
            line.snpeff_out_path,
            line.snpsift_out_path,
            self.reference_chrs_fa_path, #made in modules.core_vcf_gen
            self.reference_chrs_dict_path, #made in modules.core_vcf_gen
            self.snpEff_species_db,
            self.threads_limit,
            self.call_variants_in_parallel,
            self.cleanup
        )
        cmd = f"{VCF_GEN_SCRIPT} {' '.join(map(str, args))}"
        return cmd

    def gen_vcf_output_paths(self, name, vcf_ulid):
        self.log.attempt(f"Generating output paths and directories for line {name} with ulid:{vcf_ulid}")
        try:
            file_utils = FileUtilities(self.log)
            
            output_name = f"{vcf_ulid}-_{name}"
            vcf_output_dir_path = os.path.join(OUTPUT_DIR, output_name)
            file_utils.setup_directory(vcf_output_dir_path)
            
            vcf_output_prefix = os.path.join(vcf_output_dir_path, output_name) 
            vcf_table_path= f"{vcf_output_prefix}.table"
            snpsift_out_path=f"{vcf_output_prefix}.snpsift"
            
            snpeff_dir= os.path.join(OUTPUT_DIR, 'snpEff')
            snpeff_out_path=os.path.join(snpeff_dir, output_name)
            file_utils.setup_directory(snpeff_out_path)
            self.log.success("File paths and directories created sucessfully")
            return (vcf_output_dir_path, vcf_output_prefix, 
                vcf_table_path, snpeff_dir, snpeff_out_path, snpsift_out_path
            )
        except Exception as e:
            self.log.error(f"Error generating VCF output paths: {e}")
            raise
    def _get_ref_name(self, file_path):
        self.log.attempt('Attempting to retrieve reference base name from inputs...')
        try: 
            # Get the file name without the path
            file_name = os.path.basename(file_path)
            self.log.note(f"File name:{file_name}")
            # Remove the file extension(s)
            base_name = Path(file_name).stem
            if base_name.endswith('.fa'):
                base_name = Path(base_name).stem
            if base_name.endswith('.fasta'):
                base_name = Path(base_name).stem
            self.log.success(f"Base name: {base_name}")
        except Exception as e:
            self.log.fail(f"There as an error retrieving reference base name from inputs...")
        return base_name
    
    def _check_ref_path(self, ref_genome_path):
        self.log.attempt(f"Checking reference genome path for completeness:{ref_genome_path}")
        try:
            # Make sure we weren't just fed a directory
            if os.path.isdir(ref_genome_path):
                raise ValueError(f"Expected file, got directory: '{ref_genome_path}'")
            #Check extensions
            ext = os.path.splitext(ref_genome_path)[1]
            gz_ext = os.path.splitext(os.path.splitext(ref_genome_path)[0])[1]
            if ext not in ['.fa', '.fasta', '.gz'] or (ext == '.gz' and gz_ext not in ['.fa', '.fasta']):
                raise ValueError(f"Invalid extension: '{gz_ext+ext}'. Expected '.fa', '.fasta', '.fa.gz' or '.fasta.gz'.")

            #Check if the file already exists as a hard-coded path
            if os.path.isfile(ref_genome_path):
                self.log.note(f"Hard-coded path exists: {ref_genome_path}")
                return ref_genome_path
            else:
            # If file doesn't exist already, the base name shall be used when 
            # downloading from ref_genome_source in core_vcf_gen.py
                if self.reference_genome_source:
                    self.log.note(f"Assuming reference genome path is not hard coded. using {ref_genome_path} as base reference genome name")
                    self.log.note(f"Later we will attempt to source reference from user provided URL: {self.reference_genome_source}")
                    ref_file_path = os.path.join(REFERENCE_DIR, ref_genome_path)
                    return ref_file_path
                else:
                    self.log.fail(f"Reference genome path {ref_genome_path} does not exist, and reference genome source (a URL to download the reference from) was not provided. Aborting...")
        except Exception as e:
            self.log.fail(f"There was a critical failure while checking reference genome path:{e}")

    def _gen_reference_chrs_paths(self):
        self.reference_genome_path = self._check_ref_path(
            self.reference_genome_path
        )
        # Retrieve the base name of the reference genome
        ref_name = self._get_ref_name(self.reference_genome_path)
        
        #Reference genome prefix for .chrs.fa file
        ref_genome_prefix = os.path.join(REFERENCE_DIR, ref_name)
        self.reference_chrs_dict_path = (f"{ref_genome_prefix}.chrs.dict")
        self.reference_chrs_fa_path = (f"{ref_genome_prefix}.chrs.fa")

class BSAVariables:
    __slots__ = ['log', 'lines', 'loess_span', 'smooth_edges_bounds', 
        'shuffle_iterations', 'filter_indels', 'filter_ems', 'snpmask_path', 'snpmask_df'
    ]
    
    def __init__(self, logger,
        lines, 
        loess_span, 
        smooth_edges_bounds, 
        shuffle_iterations,
        filter_indels,
        filter_ems, 
        snpmask_path
        ):
        self.log = logger
        
        self.lines = lines # can be single Line instance or list to iterate 
        self.loess_span = loess_span
        self.smooth_edges_bounds = smooth_edges_bounds 
        self.shuffle_iterations = shuffle_iterations
        self.filter_indels = filter_indels
        self.filter_ems = filter_ems
        self.snpmask_path = snpmask_path
        self.snpmask_df = None

    def load_vcf_table(self, vcf_table_path)->pd.DataFrame:
        """
        Loads VCF table into a pandas dataframe.
        
        Args:  
        vcf_table_path - path to the vcf table to be loaded into df
        
        Returns: 
        Pandas dataframe containing the information loaded from current_line_table_path
        """
        file_utils = FileUtilities(self.log)
        self.log.attempt(f"Attempting to load {vcf_table_path} for line {current_line_name}")
        try:
            directories = [INPUT_DIR, REFERENCE_DIR] 
            vcf_table_path = file_utils.process_path(vcf_table_path, directories)
            vcf_df = pd.read_csv(current_line_table_path, sep="\t", comment='#')
            self.log.success(f"vcf_df loaded: {vcf_table_path}")
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