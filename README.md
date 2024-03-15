- [PhytoBSA](#phytobsa)
  - [Installation](#installation)
    - [Environment installation](#environment-installation)
    - [Setup directories](#setup-directories)
  - [Usage](#usage)
    - [Variables](#variables)
      - [Reference Name](#reference-name)
    - [Running](#running)
    - [Command Line Arguments](#command-line-arguments)
      - [./phytobsa analysis](#phytobsa-analysis)
      - [./phytobsa vcf\_generator](#phytobsa-vcf_generator)
    - [Output](#output)
    - [Reference Database Manager](#reference-database-manager)
      - [Features](#features)
      - [Implementation Details](#implementation-details)
      - [Reference Form Configuration](#reference-form-configuration)
        - [Fields Explanation:](#fields-explanation)
        - [Example:](#example)
# PhytoBSA

This python program analyzes and visualizes bulk segregant analysis (BSA) data. It takes sequenced segrigant bulks as an input and outputs a list of likely casual polymorphisms underlying the phenotypic segrigation of the two bulks. Likely causal polymorphisms are identified using g-statistics and a comparison of the ratio of reference to non-reference reads in each bulk. Currently, this script is optimized for running Arabidopsis BSA experiments, but can handle other organisms with some tinkering.

For a simple explanation of the experimental design of BSA -> [link to article](https://doi.org/10.1104/pp.17.00415)

This pipeline uses:

  - Delta-allele calculation. The ratio of reference read depth to total read depth in each bulk are subtracted from one another, creating a ratio that indicates phenotypic linkage.
  
  - The G-statistic calculation described in the publication: [link to publication](https://doi.org/10.1186/s12859-020-3435-8)

## Installation
### Environment installation
Install and activate the conda environment from the environment.yml file in the ./conda folder. I highly recommend using mamba (https://mamba.readthedocs.io) to install this environment, as the environment is fairly complex and conda's environment solver is comparitivly very inefficient (conda sometimes freezes trying to resolve this environment).

`mamba env create --f ./conda/environment.yml`

`mamba activate pyatbsa`

or, if you are attempting to use conda (not recommended, but probably possible)

`conda env create -f ./conda/environment.yml`

`conda activate phytobsa`

### Setup directories
run:  
`./setup.py`

## Usage
Put the fq.gz files you want analyzed into the data/input folder. You can put 
multiple experiments in the folder and they will be analyzed. 

The files must be formatted as follows:  
  
  For paired-end  
  `<line_name>.<R or D>_<read number>.<wt or mu>.fq.gz`  
    example experiment:  
    "line.R_1.wt.fq.gz"  
    "line.R_2.wt.fq.gz"   
    "line.R_1.mu.fq.gz"   
    "line.R_2.mu.fq.gz"   
  
  for unpaired  
  `<line_name>.<R or D>.<wt or mu>.fq.gz`  
    example experiment:    
    "line.R.wt.fq.gz"  
    "line.R.mu.fq.gz"       

 ### Variables
 Edit ./settings/vcf_gen_variables.sh to alter the variables for VCF file generation. The script assumes recessive polymorphisms by default. If your mutation of interest is dominant, change the mutation variable accordingly. 

 If you are analyzing an organism different than Arabidopsis, this is the place to put your reference genome and background SNPs link. In order for SNPeff to run properly, your organism needs to be named as it is in the SNPeff database. This is currently clunky to configure. Easier support for other organisms is on the roadmap. 

 #### Reference Name

- `-r`, `--reference_name`: 
  - Name of the reference genome for the input sequences. This name maps to a reference genome, SnpEff library and a snp mask. 
  --New reference names can be added to the references database using the ./refdb_manager (see [Reference Database Manager](#reference-database-manager) section for more details)
  - Type: str
  - 
### Running
 
 `./phytobsa.py -a` 

### Command Line Arguments

#### ./phytobsa analysis 
- `-ls`, `--loess_span`: 
  - Influences smoothing parameters.
  - Type: float

- `-si`, `--shuffle_iterations`: 
  - Iterations of bootstrapping during empirical cutoff calculations. Below 1000 can yield inconsistent results.
  - Type: int

- `-sb`, `--smooth_edges_bounds`: 
  - Number of mirrored datapoints at chromosome edges to correct for loess edge bias. Increase if edge bias seems high.
  - Type: int

- `-fin`, `--filter_indels`: 
  - Filter out insertion-deletion mutations.
  - Type: str

- `-fems`, `--filter_ems`: 
  - Filter results to only include mutations likely to arise from EMS treatment.
  - Type: str

- `-rco`, `--ratio_cutoff`: 
  - Used to filter results based on a ratio cutoff number. Increase to 0.2 or 0.3 if there is a lot of noise at lower ratio bounds.
  - Type: float

- `-msk`, `--mask_snps`: 
  - Set to true if you have a snpmask file configured and would like to mask known SNPs in your analysis.
  - Type: bool

- `-cc`, `--critical_cutoff`: 
  - Set the critical cutoff value for what is considered a significant polymorphism.
  - Type: float

- `-m`, `--method`: 
  - Set the method of generating the null hypothesis. Either simulate or bootstrap.
  - Type: str

#### ./phytobsa vcf_generator

- `-p`, `--call_variants_in_parallel`: 
  - Run GATK haplotype caller in parallel.
  - Type: bool

- `-c`, `--cleanup`: 
  - If true, intermediate files will be deleted. False for troubleshooting and archiving files.
  - Type: bool

- `-cft`, `--cleanup_filetypes`: 
  - Filetypes to clean out after VCF generation is complete. Format should be ['file_suffix', exc]. Example: ['.tmp', '*.metrics']
  - Type: list

- `-ocp`, `--omit_chrs_patterns`: 
  - Header patterns to omit from reference chromosomes. Useful for removing >mt (mitochondrial) and other unneeded reference sequences.
  - Type: list


### Output

Will output 6 plots: Calculated ratios, G statistics, and ratio-scaled G statistics as well as their corresponding lowess smoothed graphs. Red dashed lines represent calculated empirical cutoffs for likely candidate genes. Causal SNPs nearly always appear clearly at the top of a the lowess smoothed ratio-scaled g statistic graph. 

Example: 474-3 in https://doi.org/10.1104/pp.17.00415. Pipeline correctly identifies early stop codon in *SHR*

SRA Runs - SRR5029628(474_3_wt); SRR5029636 (474_3_mut) 

Ratio Scaled G-statistics
![RS_G_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/7a73e741-4722-4a

### Reference Database Manager

To use the Reference Database Manager, follow these steps:

1. **Fill out Configuration**: Before using the manager, ensure that the `ref_form.ini` file is filled out with the necessary configuration information.

2. **Run the Script**: Execute the script `refdb_manager.py` from the command line.

3. **Command Options**:
   - `create`: Creates a new entry in the reference database using the information provided in the `ref_form.ini` configuration file.
   - `delete`: Deletes an existing entry from the reference database. Requires specifying the `reference_name` of the entry to delete.
   - `list`: Lists all entries currently stored in the reference database. Optionally, abbreviates long URLs for better readability.

#### Features

- **Database Management**: Utilizes SQLite to manage the reference database, ensuring data integrity and efficient querying.
- **Configuration Reader**: Reads configuration information from the `ref_form.ini` file, simplifying the process of adding new entries.
- **User-Friendly Interface**: Offers a clear command-line interface with intuitive commands for creating, deleting, and listing entries.

#### Implementation Details

The Reference Database Manager consists of several components:
- **ConfigReader**: Reads configuration information from the `ref_form.ini` file.
- **RefDbUtilities**: Handles interactions with the SQLite database, including opening/closing connections, creating tables, inserting/deleting entries, and listing entries.
- **Argument Parsing**: Uses the `argparse` module to parse command-line arguments and execute the appropriate actions based on user input.

This manager provides a convenient way to manage reference information, facilitating the organization and retrieval of essential data for the `phytobsa` pipeline.

#### Reference Form Configuration

The `ref_form.ini` file contains configuration settings for reference genomes used in the `phytobsa` pipeline. Each entry in the file corresponds to a specific reference genome and provides essential information for analysis.

##### Fields Explanation:

- `reference_name`: 
  - Description: This is the name used as the reference name when running `./phytobsa` processes. All other information in this form is retrieved based on this name.

- `reference_genome_path`: 
  - Description: The path or file name where the reference genome is saved. If it doesn't exist, specify the desired path/name for saving. If the genome already exists, provide the path/name of the existing file.
  
- `reference_genome_source`: 
  - Description: This field can be None or a URL containing the reference genome FASTA file. It indicates the source of the reference genome data.
  
- `snpeff_species_db`: 
  - Description: The SNPeff database that matches the reference genome. This database allows SNPs to be labeled with their likely impact on protein function, aiding in the identification of causal mutations.
  
- `snpmask_path`: 
  - Description: Path to a file containing known SNPs. The file should include at least the chromosome, position, reference allele, and alternate allele fields.
  
- `snpmask_url`: 
  - Description: If a SNPMask is available online, provide the link here. The file will be saved with the specified filename in the `./data/references/` directory.

##### Example:

```ini
[RefDB]
reference_name = Solanum_lycopersicum
reference_genome_path = Solanum_lycopersicum_SL2.fa
reference_genome_source = ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL2.50.dna.toplevel.fa.gz
snpeff_species_db = Solanum_lycopersicum
snpmask_path = Solanum_lycopersicum_SL2.snpmask.vcf
snpmask_url = ftp://ftp.ensemblgenomes.org/pub/plants/release-32/vcf/solanum_lycopersicum/solanum_lycopersicum.vcf.gz



