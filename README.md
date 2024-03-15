- [PhytoBSA](#phytobsa)
  - [Experimental Design of BSA](#experimental-design-of-bsa)
  - [Key Features](#key-features)
  - [Installation](#installation)
    - [Environment installation](#environment-installation)
    - [Setting up the Data Directory](#setting-up-the-data-directory)
  - [Usage](#usage)
- [Default Settings](#default-settings)
  - [General Settings](#general-settings)
  - [VCF Generation Default Settings](#vcf-generation-default-settings)
  - [BSA Default Settings](#bsa-default-settings)
      - [Reference Name](#reference-name)
- [Running](#running)
  - [Automatic Mode](#automatic-mode)
  - [./phytobsa analysis](#phytobsa-analysis)
  - [./phytobsa vcf\_generator](#phytobsa-vcf_generator)
  - [Output](#output)
- [Log Database Utilities](#log-database-utilities)
  - [Functionality](#functionality)
  - [Logging Functions](#logging-functions)
  - [Usage](#usage-1)
- [Log Database Utilities](#log-database-utilities-1)
  - [Functionality](#functionality-1)
  - [Logging Functions](#logging-functions-1)
  - [Usage](#usage-2)
- [Reference Database Manager](#reference-database-manager)
  - [Reference Form Configuration](#reference-form-configuration)
# PhytoBSA

PhytoBSA is a Python program designed for analyzing and visualizing bulk segregant analysis (BSA) data. It takes sequenced segregant bulks as input and outputs a list of likely causal polymorphisms underlying the phenotypic segregation of the two bulks. PhytoBSA has been extensivly tested in Arabidopsis EMS screen populations, and lightly
tested on Rice and Tomoto QTL analysis. 

## Experimental Design of BSA
For a simple explanation of the experimental design of BSA, refer to [this article](https://doi.org/10.1104/pp.17.00415).

## Key Features

- Automatic mode
  - Once the install is configured and your files are formatted, running PhytoBSA is as simple as running 
```bash
./phytobsa -a 
```

- SNP masking. 
  - SNP masking allows the inclusion of files that contain lists of background SNPS in VCF format, so that known snps are excluded from analysis. 
    - This feature is particularly useful if your lines are divergent from the reference genome. 

- Paralell Haplotype Calling
  - if activated, Haplotypes can be called on chunks of chromosomes, scaled to the CPU resources available. This dramatically increases runtime efficiency during this
    step of VCF generation  

- Delta-Allele Calculation
  - Calculates the delta-allele, which is the ratio of reference read depth to total read depth in each bulk. By subtracting these ratios from one another, it creates a value indicating phenotypic linkage.

- G-Statistic Calculation
  - Implements the G-statistic calculation described in the publication [here](https://doi.org/10.1186/s12859-020-3435-8). This statistic helps identify likely causal polymorphisms by comparing the ratio of reference to non-reference reads in each bulk.

- ULIDs for File and Analysis Identification
  - Utilizes ULIDs (Universally Unique Lexicographically Sortable Identifiers) to ensure that each generated file and analysis is uniquely identified. This allows the program to handle concurrent analyses pointing to the same output directory without conflicts.

- Robust Logging of Run Parameters
  - Incorporates robust logging of run parameters, making debugging and reproducibility of results easier to track. This logging system captures all relevant parameters used in each analysis, aiding in result interpretation and replication.

- Bayesian-Based Simulation for Critical Cutoff Values
  - Utilizes Bayesian-based simulation to produce analytically tractable and robust critical cutoff values for SNPs. This approach enhances the reliability of identifying significant polymorphisms.

## Installation
### Environment installation
Install and activate the conda environment from the environment.yml file in the ./conda folder. I highly recommend using mamba (https://mamba.readthedocs.io) to install this environment, as the environment is fairly complex and conda's environment solver is comparitivly very inefficient (conda sometimes freezes trying to resolve this environment).

`mamba env create --f ./conda/environment.yml`

`mamba activate pyatbsa`

or, if you are attempting to use conda (not recommended, but probably possible)

`conda env create -f ./conda/environment.yml`

`conda activate phytobsa`

### Setting up the Data Directory

Before running the program, it's essential to set up the data directory. The data directory serves as the storage location for various files required by the program, including reference databases and large genomic files. Adequate storage space should be available in this directory to accommodate these files. This is also the 
directory in which the outputs of the program will be stored. 

To set up the data directory:

**Configure Data Directory**:

`./phytobsa settings --set_data_dir <path-to-directory>`

Replace `<path-to-directory>` with the desired path for the data directory.

Upon first execultion of the program, the specified data directory 
  will be installed initialized.

Once the data directory is set up, the program will be able to access the required files and directories smoothly.


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

# Default Settings

PhytoBSA offers default settings that can be applied to streamline the analysis process. These default settings allow users to set preferred configurations for various parameters, ensuring consistency and reducing the need for manual configuration for each run. Below is a breakdown of the default settings available for configuration:

## General Settings
These settings are automatically applied if not explictly passed in any mode. 

- `--set_data_dir`: Set the data directory. This must be set for the program to run.
- `--set_threads_limit`: Set the threads limit for BSA and for VCF generation. If not set, threads will be detected, and threads -2 will be used.
- `--set_reference_name`: Set the name of the reference genome.

## VCF Generation Default Settings
These settings are automatically applied if not explicitly provided in automatic or VCF generation mode.

- `--set_call_variants_in_parallel`: Set default for running GATK Haplotype Caller in parallel.
- `--set_cleanup`: Set default for cleanup. If true, intermediate files will be deleted; false for troubleshooting and archiving files.
- `--set_cleanup_filetypes`: Set default for cleanup file types. An ordered list of globs for files to clear out after VCF generation process.
- `--set_omit_chrs_patterns`: Set defaults for filtering reference chromosome contigs. Useful for filtering non-genomic reference contigs to speed up VCF generation.

## BSA Default Settings
These settings are automatically applied if not explicitly passed to automatic or BSA mode.

- `--set_loess_span`: Set default Loess span.
- `--set_shuffle_iterations`: Set default shuffle iterations.
- `--set_smooth_edges_bounds`: Set default smooth edges bounds.
- `--set_filter_indels`: Set default filter indels.
- `--set_filter_ems`: Set default filter EMS.
- `--set_ratio_cutoff`: Set default ratio cutoff bound.
- `--set_mask_snps`: Set default mask SNPs boolean value.
- `--set_critical_cutoff`: Set default critical cutoff value.
- `--set_method`: Set the default method of generating the null hypothesis. Either 'simulate' or 'bootstrap'.

Users can apply these default settings using the `phytobsa settings` command with the corresponding options. This feature is particularly useful for users who primarily work with specific reference genomes, species, or analysis methodologies, as it eliminates the need for repetitive configuration adjustments.

 #### Reference Name

- `-r`, `--reference_name`: 
  - Name of the reference genome for the input sequences. This name maps to a reference genome, SnpEff library and a snp mask.
  --New reference names can be added to the references database using the ./refdb_manager (see [Reference Database Manager](#reference-database-manager) section for more details)
  - Type: str
  
# Running
 
## Automatic Mode
 Assuming the data directory is configured, phytobsa conda environment is activated and you files are properly formated, 
 all that is needed to run the analysis is the following: 

 `./phytobsa.py -a` 


## ./phytobsa analysis 
This command line argument allows the running of independant analysis.

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

## ./phytobsa vcf_generator
This command allows the generation of VCF files independantly of running the analysis. 
As with all other commands, you can set the default settings using ./phytobsa settings

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

## Output

The tool will generate a total of 6 plots, including calculated ratios, G statistics, and ratio-scaled G statistics, along with their corresponding lowess smoothed graphs.

Ratio-scaled G-statistics are calculated by multiplying the G-statistic with the ratio. This combined metric tends to provide more stable results compared to using either feature alone.

In the plots, you will notice nested gray ribbons that represent the null model, showing percentiles (1st, 25th, 50th, 75th, and 99th) at each position. These null models are generated through a Bayesian simulation process using bootstrapped reads.

During the simulation, the reference reads are distributed binomially, and a conjugate prior (beta) describes the frequency of reference alleles in each bulk sample. The allele frequencies of the reference alleles are simulated from the posterior distribution after updating with the bootstrapped values.

By using bootstrapping, the tool breaks the link between phenotypes and genotypes, while the Bayesian simulation introduces a soft constraint to the simulated values, guiding extreme values towards a reference allele frequency of 0.5.

Significant polymorphisms are identified based on their position above the critical cutoff percentile in the null model.

Example: 474-3 in https://doi.org/10.1104/pp.17.00415. Pipeline correctly identifies early stop codon in *SHR*

SRA Runs - SRR5029628(474_3_wt); SRR5029636 (474_3_mut) 

Ratios, Lowess smoothed  
<img src="https://github.com/TeaShull/PhytoBSA/assets/125574642/3158aac9-60da-4ce9-9e70-c099e80c1082" width="600">

G-statistics, Lowess smoothed  
<img src="https://github.com/TeaShull/PhytoBSA/assets/125574642/b3c8ca0b-a799-4252-9fc5-281edf9e5dfa" width="600">

Ratio Scale G-statistics, Lowess smoothed  
<img src="https://github.com/TeaShull/PhytoBSA/assets/125574642/709423e3-4313-4595-894a-3dc43ea89ee2" width="600">

Finally, a list of the likely candidates will be produced: 

# Log Database Utilities

The Log Database Utilities module provides functions to interact with a log database, allowing users to easily track and retrieve runtime parameters and associated information. This is crucial for ensuring reproducibility and comparability of results across different runs of analysis or processing tasks.

## Functionality

1. **Print Analysis Log Data**
   - Retrieves and prints information related to an analysis based on the analysis ID provided (ulid).
   - Command: `logdb -an ANALYSIS_ULID`

2. **Print VCF Log Data**
   - Retrieves and prints information related to a Variant Call Format (VCF) based on the VCF ID or core ID provided (ulid).
   - Command: `logdb -vcf VCF_ULID`

3. **Get Line Name Data**
   - Retrieves all entries related to a specific line name and returns the results as a list.
   - Command: `logdb -name LINE_NAME`

4. **Print Line Name Data**
   - Prints all entries related to a specific line name, including both VCF data and Analysis data.
   - Command: `logdb -name LINE_NAME`

5. **Print Core ID Data**
   - Retrieves and prints information related to a core ID, including core log data, VCF data, and Analysis data linked to that core ID.
   - Command: `logdb -core CORE_ULID`

## Logging Functions

1. **Create Tables**
   - Creates the necessary tables in the log database to store core, VCF, and analysis log data.

2. **Add Database Record**
   - Adds records to the log database based on the type of log (core, VCF, or analysis) and the provided parameters.

## Usage

Users can utilize the logdbutils functions to store, retrieve, and analyze runtime parameters and associated data in a structured manner. By logging this information, users can maintain a record of the processes and configurations used for each run, enabling reproducibility and comparison of results across different executions.

# Log Database Utilities

The Log Database Utilities module provides functions to interact with a log database, allowing users to easily track and retrieve runtime parameters and associated information. This is crucial for ensuring reproducibility and comparability of results across different runs of analysis or processing tasks.

## Functionality

1. **Print Analysis Log Data**
   - Retrieves and prints information related to an analysis based on the analysis ID provided (ulid).
   - Command: `logdb -an ANALYSIS_ULID`

2. **Print VCF Log Data**
   - Retrieves and prints information related to a Variant Call Format (VCF) based on the VCF ID or core ID provided (ulid).
   - Command: `logdb -vcf VCF_ULID`

3. **Get Line Name Data**
   - Retrieves all entries related to a specific line name and returns the results as a list.
   - Command: `logdb -name LINE_NAME`

4. **Print Line Name Data**
   - Prints all entries related to a specific line name, including both VCF data and Analysis data.
   - Command: `logdb -name LINE_NAME`

5. **Print Core ID Data**
   - Retrieves and prints information related to a core ID, including core log data, VCF data, and Analysis data linked to that core ID.
   - Command: `logdb -core CORE_ULID`

## Logging Functions

1. **Create Tables**
   - Creates the necessary tables in the log database to store core, VCF, and analysis log data.

2. **Add Database Record**
   - Adds records to the log database based on the type of log (core, VCF, or analysis) and the provided parameters.

## Usage

Users can utilize the logdbutils functions to store, retrieve, and analyze runtime parameters and associated data in a structured manner. By logging this information, users can maintain a record of the processes and configurations used for each run, enabling reproducibility and comparison of results across different executions.

# Reference Database Manager

To use the Reference Database Manager, follow these steps:

1. **Fill out Configuration**: Before using the manager, ensure that the `ref_form.ini` file is filled out with the necessary configuration information.

2. **Run the Script**: Execute the script `refdb_manager.py` from the command line.

3. **Command Options**:
   - `create`: Creates a new entry in the reference database using the information provided in the `ref_form.ini` configuration file.
   - `delete`: Deletes an existing entry from the reference database. Requires specifying the `reference_name` of the entry to delete.
   - `list`: Lists all entries currently stored in the reference database. Optionally, abbreviates long URLs for better readability.


## Reference Form Configuration

The `ref_form.ini` file contains configuration settings for reference genomes used in the `phytobsa` pipeline. Each entry in the file corresponds to a specific reference genome and provides essential information for analysis.

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

Example:

```ini
[RefDB]
reference_name = Solanum_lycopersicum
reference_genome_path = Solanum_lycopersicum_SL2.fa
reference_genome_source = ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL2.50.dna.toplevel.fa.gz
snpeff_species_db = Solanum_lycopersicum
snpmask_path = Solanum_lycopersicum_SL2.snpmask.vcf
snpmask_url = ftp://ftp.ensemblgenomes.org/pub/plants/release-32/vcf/solanum_lycopersicum/solanum_lycopersicum.vcf.gz



