import os
import vcf
import csv
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
        'analysis_out_prefix', 'gs_cutoff', 'rsg_cutoff', 'rsg_y_cutoff', 
        'analysis_ulid', 'in_path_variables', 'ref_path_variables', 
        'snpeff_dir', 'snpeff_out_path', 'snpsift_out_path', 
        'snpeff_report_path', 'rs_cutoff', 'snpmask_df'
    ]

    def __init__(self, logger, name):
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
        self.snpeff_report_path = None
        self.snpeff_out_path = None
        self.snpsift_out_path = None


        # BSA variables
        self.segregation_type = None
        self.vcf_df = None
        self.snpmask_df = None
        self.analysis_out_prefix = None
        self.gs_cutoff = None
        self.rs_cutoff = None
        self.rsg_cutoff = None
        self.rsg_y_cutoff = None
        self.analysis_ulid = None

        
        # Input Variables that require path parsing. _process_input checks if 
        # hard-coded or if in the /inputs folder  
        self.in_path_variables=[
           'vcf_table_path',
           'mu_input',
           'wt_input'
        ]


    def _process_input(self, key, details):
        file_utils = FileUtilities(self.log)
        if key in self.in_path_variables:
            paths = details.split() #for space seperated lists mu_input and wt_input
            dir = [INPUT_DIR]
            processed_paths = [file_utils.process_path(dir, path) for path in paths]
            return ' '.join(processed_paths)
        else:
            return details
    

    def usr_in_line_variables(self, **kwargs):
        self.log.attempt('Attempting to organize inputs into Line data class...')
        try:
            for key, details in kwargs.items():
                try:
                    processed_input = self._process_input(key, details)
                
                except KeyError as ke:
                    self.log.fail(f'Key not found: {ke}')
                    
                    continue
                 
                try:
                    if hasattr(self, key):
                        setattr(self, key, processed_input)
                        self.log.note(f'User input added: |{self.name}|{key}:{processed_input} ')

                except AttributeError as ae:
                    self.log.fail(f'Attribute error: {ae}')

        except Exception as e:
            self.log.fail(f'There was an error processing user inputs:{e}')


class AutomaticLineVariableDetector:
    def __init__(self, logger):
        self.log = logger
        self.lines = []


    def _parse_filename(self, filename):
        parts = filename.split('.')
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
        details = self._parse_filename(filename)
        name, segregation_type, bulk_type, pairedness = details
        
        self.log.note(f""" 
        -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                  Details Detected
        -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        name:{name}
        segregation_type:{segregation_type}
        bulk_type:{bulk_type}
        pairedness:{pairedness}
        -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=""")

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
                print(" ")
                self.log.attempt(f"Parsing {filename}....")
                self._process_file(filename)
                self.log.success(f"{filename} parsed!")
            
            self.log.note("Sorting file paths...")
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
            line.name, 
            line.wt_input,
            line.mu_input,
            line.pairedness,
            line.vcf_output_prefix,
            line.snpeff_report_path,
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
            self.log.note(f"Output path variable generated:{output_name}")
            
            vcf_output_dir_path = os.path.join(OUTPUT_DIR, output_name)
            self.log.note(f"Output path variable generated:{vcf_output_dir_path}")
            
            file_utils.setup_directory(vcf_output_dir_path) #setup vcf_gen output directory
            
            vcf_output_prefix = os.path.join(vcf_output_dir_path, output_name) 
            self.log.note(f"Output path variable generated:{vcf_output_prefix}")

            snpeff_out_path=f"{vcf_output_prefix}.se.vcf"
            self.log.note(f"Output path variable generated:{snpeff_out_path}")
            
            snpsift_out_path=f"{vcf_output_prefix}.snpsift"
            self.log.note(f"Output path variable generated:{snpsift_out_path}")
            
            vcf_table_path= f"{vcf_output_prefix}.table"
            self.log.note(f"Output path variable generated:{vcf_table_path}")
            
            snpeff_report_dir=os.path.join(vcf_output_dir_path, 'snpEff')
            self.log.note(f"Output path variable generated:{snpeff_report_dir}")
            
            file_utils.setup_directory(snpeff_report_dir) #setup snpeff_report directory
            
            snpeff_report_path = os.path.join(snpeff_report_dir, output_name)
            self.log.note(f"Output path variable generated:{snpeff_report_path}")

            self.log.success("File paths and directories created sucessfully")
            return (vcf_output_dir_path, vcf_output_prefix, 
                vcf_table_path, snpeff_report_path, snpeff_out_path, 
                snpsift_out_path
            )
        except Exception as e:
            self.log.fail(f"Error generating VCF output paths: {e}")

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
            self.log.fail(f"There as an error retrieving reference base name from inputs:{e}")
        return base_name
    

    def _check_ref_path(self, ref_genome_path):
        self.log.attempt(f"Checking reference genome path for completeness:{ref_genome_path}")

        # Raise error if we were fed a directory
        if os.path.isdir(ref_genome_path):
            self.log.fail(f"Expected file, got directory: '{ref_genome_path}'")
            return

        # Check extensions
        ext = os.path.splitext(ref_genome_path)[1]
        gz_ext = os.path.splitext(os.path.splitext(ref_genome_path)[0])[1]
        if ext not in ['.fa', '.fasta', '.gz'] or (ext == '.gz' and gz_ext not in ['.fa', '.fasta']):
            self.log.fail(f"Invalid extension: '{gz_ext+ext}'. Expected '.fa', '.fasta', '.fa.gz' or '.fasta.gz'")
            return

        # Check if the file already exists as a hard-coded path
        if os.path.isfile(ref_genome_path):
            self.log.note(f"Hard-coded path exists: {ref_genome_path}")
            return ref_genome_path

        # If file doesn't exist already, the base name shall be used when 
        # downloading from ref_genome_source in core_vcf_gen.py
        self.log.note(f"Assuming reference genome path is not hard coded. Using {ref_genome_path} as base reference genome name")
        ref_file_path = os.path.join(REFERENCE_DIR, ref_genome_path)

        if os.path.isfile(ref_file_path):
            self.log.note(f"Found {ref_genome_path} in the references directory. Proceeding...")
            return ref_file_path

        if not self.reference_genome_source:
            self.log.fail(f"Reference genome path {ref_genome_path} does not exist, and reference genome source (a URL to download the reference from) was not provided. Aborting...")
            return

        self.log.note(f"{ref_genome_path} is not the full path or found in ./references.") 
        self.log.note(f"Later we will attempt to source reference from user provided URL: {self.reference_genome_source}")



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
        'shuffle_iterations', 'filter_indels', 'filter_ems', 'snpmask_path', 
        'snpmask_df', 'ratio_cutoff'
    ]
    

    def __init__(self, logger,
        lines, 
        loess_span, 
        smooth_edges_bounds, 
        shuffle_iterations,
        filter_indels,
        filter_ems, 
        snpmask_path,
        ratio_cutoff
        ):
        self.log = logger
        
        self.lines = lines # can be single Line instance or list to iterate 
        self.loess_span = loess_span
        self.smooth_edges_bounds = smooth_edges_bounds 
        self.shuffle_iterations = shuffle_iterations
        self.filter_indels = filter_indels
        self.filter_ems = filter_ems
        self.snpmask_path = snpmask_path
        self.ratio_cutoff = ratio_cutoff
        self.snpmask_df = None


    def load_vcf_table(self, vcf_table_path)->pd.DataFrame:
        """
        Loads VCF table into a pandas dataframe.

        Args:  
        vcf_table_path - path to the vcf table to be loaded into df

        Returns: 
        Pandas dataframe containing the information loaded from vcf_table_path
        """
        file_utils = FileUtilities(self.log)
        self.log.attempt(f"Attempting to load {vcf_table_path}")
        try:
            directories = [INPUT_DIR] 
            vcf_table_path = file_utils.process_path(directories, vcf_table_path)
            df = pd.read_csv(vcf_table_path, sep="\t")
            self.log.note(f"{vcf_table_path} loaded. Proceeding...")
        
            return df
        
        except pd.errors.EmptyDataError:
            self.log.fail(f"Error: File '{vcf_table_path}' is empty.")
        
        except Exception as e:
            self.log.fail(f"An unexpected error occurred: {e}")


    def load_snpmask(self, snpmask_path)->pd.DataFrame:
        """
        Handles SNP mask and converts it into a pandas dataframe.

        Args:  
        snpmask_path - path to the vcf table to be loaded into df

        Returns: 
        Pandas dataframe containing the information loaded from vcf_table_path
        """
        file_utils = FileUtilities(self.log)
        self.log.attempt(f"Attempting to load {snpmask_path}")
        try:
            directories = [REFERENCE_DIR] 
            snpmask_path = file_utils.process_path(directories, snpmask_path)
            df = pd.read_csv(snpmask_path, sep="\t")
            
            return df

        except pd.errors.ParserError:
            self.log.note(f"snpmask file doesn't appear to be formatted for pandas. Formatting....")
            temp_table_path = f"{snpmask_path}.temp"
            backup_table_path = f"{snpmask_path}.backup"
            headers = ['CHROM','POS','REF','ALT' ]

            with open(snpmask_path, 'r') as file:
                vcf_reader = vcf.Reader(file, strict_whitespace=True)
                with open(temp_table_path, 'w', newline='') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(headers)
                    while True:
                        try:
                            record = next(vcf_reader)
                        except StopIteration:
                            break
                        except IndexError:
                            self.log.warning(f"Skipping record due to incorrect number of fields:{record}")
                            continue

                        try:
                            chrom = record.CHROM
                            pos = record.POS
                            ref = record.REF
                            alt = ','.join(str(a) for a in record.ALT) 
                            writer.writerow([chrom, pos, ref, alt])
                        except AttributeError:
                            self.log.warning(f"Skipping record due to missing field(s): {record}")

            # Save original file as backup
            os.rename(snpmask_path, backup_table_path)
            self.log.note(f"Backup saved:{backup_table_path}")

            # Replace original file with the new one
            os.rename(temp_table_path, snpmask_path)
            self.log.note(f"{snpmask_path} converted to friendlier format for masking. Original file can be found here: {backup_table_path}")

            # Load the newly written file into a DataFrame
            df = pd.read_csv(snpmask_path, sep="\t", names=headers)

            self.log.success(f"vcf_df loaded and saved: {snpmask_path}")

            return df

        except pd.errors.EmptyDataError:
            self.log.fail(f"Error: File '{snpmask_path}' is empty.")

        except Exception as e:
            self.log.fail(f"An unexpected error occurred: {e}")


    def gen_bsa_out_prefix(self, name, ulid, vcf_ulid): #called in core_bsa.py
        try:
            analysis_out_prefix = f'{ulid}-_{name}'
            if vcf_ulid:
                out_path = os.path.join(
                        OUTPUT_DIR, 
                        f'{vcf_ulid}-_{name}',
                        analysis_out_prefix, 
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
        
        except Exception as e:
            self.log.error(f"Error generating BSA output prefix: {e}")
            
            return None

