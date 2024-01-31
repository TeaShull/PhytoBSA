import os
import subprocess
import urllib.request
from urllib.parse import urlparse
import gzip
import shutil
import re

from Bio import SeqIO

from modules.utilities_logging import LogHandler
from settings.paths import MODULES_DIR

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
        self._download_reference_genome(self.vcf_gen_vars.reference_genome_path, 
            self.vcf_gen_vars.reference_genome_source
        )
        # Create a new fasta file, with unneeded contigs filtered based on user defaults. (mitochondria, exc)
        self._create_chromosomal_fasta(self.vcf_gen_vars.reference_genome_path, 
            self.reference_chrs_path, self.vcf_gen_vars.omit_chrs_patterns
        )

        # Generate the knownSnps .vcf file path
        for line in self.vcf_gen_vars.lines:
            self.log.note(f'Initializing vcf_generation subprocess log for {line.name}')
            vcf_log = LogHandler(f'vcf_{line.name}')
            line.vcf_ulid = vcf_log.ulid
            vcf_log.add_db_record(line.name, line.vcf_ulid)

            #Generate output paths for process
            vcf_out_paths = self.vcf_gen_vars.gen_vcf_output_paths(
                line.name, line.vcf_ulid
            )
            (
                line.vcf_output_dir, 
                line.vcf_output_prefix, 
                line.vcf_table_path,
                line.snpeff_dir,
                line.snpeff_out_path, 
                line.snpsift_out_path
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
            
            # Formatting VCF output to make it dataframe friendly
            vcf_format = VCFFormat(line.vcf_table_path)
            vcf_format.remove_repetitive_nan(line.snpsift_out_path)
            vcf_format.format_fields()
            vcf_format.remove_complex_genotypes()
            vcf_format.add_headers()
            
            #Cleanup files if cleanup = True
            if cleanup:
                self._cleanup_files(line.vcf_output_dir)
            
    def _download_reference_genome(ref_genome_path, ref_genome_source):
        # Parse the source URL to get the file extension
        parsed_url = urlparse(ref_genome_source)
        source_extension = os.path.splitext(parsed_url.path)[1]

        # Adjust the local file path if necessary
        local_extension = os.path.splitext(ref_genome_path)[1]
        if source_extension != local_extension:
            ref_genome_path += source_extension

        if not os.path.isfile(ref_genome_path):
            urllib.request.urlretrieve(ref_genome_source, ref_genome_path)

            # If the file is compressed, uncompress it
            if source_extension == '.gz':
                with gzip.open(ref_genome_path, 'rb') as f_in:
                    with open(ref_genome_path[:-3], 'wb') as f_out:  # Remove '.gz' from the end of the path
                        shutil.copyfileobj(f_in, f_out)

    def _create_chromosomal_fasta(input_file, output_file, *patterns):
        if not os.path.isfile(output_file):
            print(f"Creating {output_file} with only chromosomal DNA...")
            print(f"Omitting contigs containing {', '.join(patterns)}")
            with open(input_file, "r") as in_handle, open(output_file, "w") as out_handle:
                for record in SeqIO.parse(in_handle, "fasta"):
                    if not any(re.search(pattern, record.id) for pattern in patterns):
                        SeqIO.write(record, out_handle, "fasta")
        else:
            print(f"Chromosomal fasta exists: {output_file}")
            print("Proceeding...")
    
    def _cleanup_files(output_dir_path, cleanup_filetypes):
        for file_type in cleanup_filetypes:
            if file_type == '*.table':
                print("*.table not allowed in cleanup_filetypes. This is the point of running VCF_gen in the first place. Continuing...")
                continue
            files = glob.glob(os.path.join(output_dir_path, file_type))
            for file in files:
                if os.path.isfile(file):
                    os.remove(file)

class VCFFormat:
    def __init__(self, vcf_table_path):
        self.vcf_table_path = vcf_table_path
    
    def remove_repetitive_nan(self, input):
        with open(input, 'r') as f_in, open(self.vcf_table_path, 'w') as f_out:
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

    def add_headers(self):
        # Read the file into a pandas DataFrame
        df = pd.read_csv(self.vcf_table_path, sep='\t', header=None)
    
        # Add the header
        df.columns = [
            'chr', 
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
    
        # Write the DataFrame back to the file
        df.to_csv(self.vcf_table_path, sep='\t', index=False)
    

