
# PhytoBSA
PhytoBSA is a Python program designed for analyzing and visualizing bulk segregant analysis (BSA) data. It takes sequenced segregant bulks as input and outputs a list of likely causal polymorphisms underlying the phenotypic segregation of the two bulks. PhytoBSA has been extensively tested in Arabidopsis EMS screen populations, and lightly
tested on Rice and Tomato QTL analysis. 

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

- Parallel Haplotype Calling
  - if activated, Haplotypes can be called on chunks of chromosomes, scaled to the CPU resources available. This dramatically increases runtime efficiency during this
    step of VCF generation  

- Delta-Allele Calculation
  - Calculates the delta-allele, which is the ratio of reference read depth to total read depth in each bulk. By subtracting these ratios from one another, it creates a value indicating phenotypic linkage.

- G-Statistic Calculation
  - Implements the G-statistic calculation described in the publication [here](https://doi.org/10.1186/s12859-020-3435-8). This statistic helps identify likely causal polymorphisms by comparing the ratio of reference to non-reference reads in each bulk.

- Bayesian Simulation of Null Models
  - Utilizes Bayesian simulation to produce analytically tractable and robust critical cutoff values for SNPs. This approach enhances the reliability of identifying significant polymorphisms.

- ULIDs for File and Analysis Identification
  - Utilizes ULIDs (Universally Unique Lexicographically Sortable Identifiers) to ensure that each generated file and analysis is uniquely identified. This allows the program to handle concurrent analyses pointing to the same output directory without conflicts. Every process has a log, which is linked to the generated files through the ULID. 

- Robust Logging of Run Parameters
  - Incorporates robust logging of run parameters. This logging system captures all relevant parameters used in each analysis, aiding in result interpretation and replication.

## Output

The tool will generate a total of 9 plots, including calculated ratios, G statistics, and ratio-scaled G statistics, along with their corresponding percentiles and lowess smoothed graphs.

Ratio-scaled G-statistics are calculated by multiplying the G-statistic with the ratio. This combined metric tends to provide more stable results compared to using either feature alone.

### Identification of Significant Polymorphisms

In the generated plots, you will observe nested gray ribbons illustrating the null model, displaying percentiles (1st, 25th, 50th, 75th, and 99th) at each position. These null models are created through a Bayesian simulation process using bootstrapped reads.

The tool employs chromosome-wise random resampling with replacement to sever the link between phenotypes and genotypes. This random resampling, along with the Bayesian simulation, imposes a soft constraint on the simulated values, guiding values towards a reference allele frequency of 0.5. Additionally, the random draw from the posterior introduces extra variance in the in the simulated read depths. More details can be found in the footnote: [Null Model Simulation](#null-model-simulation)

Significant polymorphisms are identified based on their position above the critical cutoff percentile in the null model.

Example: 474-3 in https://doi.org/10.1104/pp.17.00415. Pipeline correctly identifies early stop codon in *SHR*

SRA Runs - SRR5029628(474_3_wt); SRR5029636 (474_3_mut) 

Ratios, Lowess smoothed  
<img src="https://github.com/TeaShull/PhytoBSA/assets/125574642/3158aac9-60da-4ce9-9e70-c099e80c1082" width="600">

G-statistics, Lowess smoothed  
<img src="https://github.com/TeaShull/PhytoBSA/assets/125574642/b3c8ca0b-a799-4252-9fc5-281edf9e5dfa" width="600">

Ratio Scale G-statistics, Lowess smoothed  
<img src="https://github.com/TeaShull/PhytoBSA/assets/125574642/709423e3-4313-4595-894a-3dc43ea89ee2" width="600">

Finally, a list of the likely candidates will be produced, filtered based on where each locus lies in the null model percentiles and SNP impact. 


# Table of Contents
- [Installation](#installation)
  - [Environment installation](#environment-installation)
  - [Setting up the Data Directory](#setting-up-the-data-directory)
- [Commands](#commands)
  - [./phytobsa -a (Automatic mode)](#phytobsa--a-automatic-mode)
  - [./phytobsa analysis](#phytobsa-analysis)
  - [./phytobsa vcf\_generator](#phytobsa-vcf_generator)
- [Default Settings](#default-settings)
  - [Set General Defaults](#set-general-defaults)
  - [Set VCF Generation Defaults](#set-vcf-generation-defaults)
  - [Set BSA Defaults](#set-bsa-defaults)
- [Log Database Utilities](#log-database-utilities)
- [Reference Database Manager](#reference-database-manager)
  - [Reference Form Configuration](#reference-form-configuration)
- [Footnotes](#footnotes)
  - [Null Model Simulation](#null-model-simulation)

# Installation
## Environment installation
Install and activate the conda environment from the environment.yml file in the ./conda folder. I highly recommend using mamba (https://mamba.readthedocs.io) to install this environment, as the environment is fairly complex and conda's environment solver is comparatively very inefficient (conda sometimes freezes trying to resolve this environment).

`mamba env create --f ./conda/environment.yml`

`mamba activate pyatbsa`

or, if you are attempting to use conda (not recommended, but probably possible)

`conda env create -f ./conda/environment.yml`

`conda activate phytobsa`

## Setting up the Data Directory

Before running the program, it's essential to set up the data directory. The data directory serves as the storage location for various files required by the program, including reference databases and large genomic files. Adequate storage space should be available in this directory to accommodate these files. This is also the 
directory in which the outputs of the program will be stored. 

To set up the data directory:

**Configure Data Directory**:

`./phytobsa settings --set_data_dir <path-to-directory>`

Replace `<path-to-directory>` with the desired path for the data directory.

Upon first execution of the program, the specified data directory 
  will be installed initialized.

Once the data directory is set up, the program will be able to access the required files and directories smoothly.


# Commands
## ./phytobsa -a (Automatic mode)
***File Formating***  
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

***Reference Configuration***  
Reference name must be configured. 
Currently, there are a few pre-configured references that will auto-populate genomes, snpmasks and the appropriate SnpEff library
- Arabidopsis_thaliana
- Zea_mays
- Oryza_sativa_Japanica
- Solanum_lycopersicum

An example of setting the Arabidopsis thaliana default via the command line- 
`./phytobsa settings --set_reference_name Arabidopsis_thaliana`

This will autopopulate the program with reference genome names, the URL to download the reference genome files, the appropriate SnpEff library to use, 
and a good general snpmask. 

For more information about configuring new genomes for the reference database manager: [Reference Database Manager](#reference-database-manager)

***Settings configuration***   
For automatic mode to run well, configure your default settings using ./phytobsa settings (see [Default Settings](#default-settings) or ./phytobsa settings -h)
Another option is to directly modify settings.ini


## ./phytobsa analysis 
This command line argument allows the running of analysis independently oh the automatic workflow.
This can make the fine tuning of analysis easier. 

***Required***  

If running analysis separately, these variables can't be set using the config. 
| Command              | ab    | Description                                                                                      | Type |
| -------------------- | ----- | ------------------------------------------------------------------------------------------------ | ---- |
| `--name`             | `-n`  | name of the line for output files.                                        | str  |
| `--vcf_table_path`   | `-vt` | path to the VCF table or file in input or output.                | str  |
| `--segregation_type` | `-st` | Specify the segregation type. (R), (D) or (QTL). | str  |

  
***Options*** (see [Default Settings](#default-settings) to alter default values)

| Command                 | ab      | Description                                                                               | Type  |
| ----------------------- | ------- | ----------------------------------------------------------------------------------------- | ----- |
| `--reference_name`      | `-r`    | Name of the reference genome. This name maps to a reference genome.                       | str   |
| `--loess_span`          | `-ls`   | Smoothing parameter.                                                                      | float |
| `--shuffle_iterations`  | `-si`   | Iterations of bootstrapping during null model generation.                         | int   |
| `--smooth_edges_bounds` | `-sb`   | Mirrored data at chrom edges to correct for loess edge bias. | int   |
| `--filter_indels`       | `-fin`  | Filter out insertion-deletion mutations.                                                  | str   |
| `--filter_ems`          | `-fems` | Filter results to only include mutations likely to arise from EMS treatment.              | str   |
| `--ratio_cutoff`        | `-rco`  | Used to filter results based on a ratio cutoff number.                                    | float |
| `--mask_snps`           | `-msk`  | Mask known SNPs in analysis.                                                              | bool  |
| `--critical_cutoff`     | `-cc`   | Set the critical cutoff value for significant polymorphism.                               | float |



## ./phytobsa vcf_generator  
This command allows the generation of VCF files independently of running the analysis. 
As with all other commands, you can set the default settings using ./phytobsa settings
***Required***  

| Command           | ab    | Description                                                                                                           | Type |
| ----------------- | ----- | --------------------------------------------------------------------------------------------------------------------- | ---- |
| `--name`      | `-n`  | Specify the name which will be used to name output files.                                                             | str  |
| `--wt_input` | `-wt` | Specify the path(s) to the wild-type bulk fasta file(s). Can be file in input or direct path. Format: "FILE_1 FILE_2" | str  |
| `--mu_input` | `-mu` | Specify the path(s) to the mutant bulk fasta file(s). Can be file in input or direct path. Format: "FILE_1 FILE_2"    | str  |


***Options*** (see [Default Settings](#default-settings) to alter default values)

| Command                       | ab     | Description                                                                                         | Type |
| ----------------------------- | ------ | --------------------------------------------------------------------------------------------------- | ---- |
| `--reference_name`            | `-r`   | Name of the reference genome. This name maps to a reference genome, SnpEff library, and a snp mask. | str  |
| `--call_variants_in_parallel` | `-p`   | Run GATK haplotype caller in parallel.                                                              | bool |
| `--cleanup`                   | `-c`   | Cleanup intermediate files?                                                                         | bool |
| `--cleanup_filetypes`         | `-cft` | Filetypes to clean out after VCF generation is complete. Example: ['.tmp', '*.metrics']             | list |
| `--omit_chrs_patterns`        | `-ocp` | Header patterns to omit from reference chromosomes. Removes unneeded reference sequences.           | list |



# Default Settings  

PhytoBSA offers default settings that can be applied to streamline the analysis process. These default settings allow users to set preferred configurations for various parameters, ensuring consistency and reducing the need for manual configuration for each run. Below is a breakdown of the default settings available for configuration:

*Note - you can also configure these settings directly in settings/config.ini, if you so wish*

## Set General Defaults  
These settings are automatically applied if not explicitly passed in any mode. 
| Command                | Description                                                                          | Possible Values                                       |
| ---------------------- | ------------------------------------------------------------------------------------ | ----------------------------------------------------- |
| `--set_reference_name` | Set the name of the reference genome.                                                | string [see refdb info](#reference-database-manager)) |
| `--set_data_dir`       | Set the data directory. This must be set for the program to run.                     | Path to data directory                                |
| `--set_threads_limit`  | Set the threads limit for BSA and for VCF generation. default is detected threads -2 | int                                                   |



## Set VCF Generation Defaults  
These settings are automatically applied if not explicitly provided in automatic or VCF generation mode.

- `--set_call_variants_in_parallel`: 
  - Set default for running GATK Haplotype Caller in parallel.

- `--set_cleanup`: 
  - Set default for cleanup. If true, intermediate files will be deleted; false for troubleshooting and archiving files.

- `--set_cleanup_filetypes`: 
  - Set default for cleanup file types. An ordered list of globs for files to clear out after VCF generation process.

- `--set_omit_chrs_patterns`: 
  - Set defaults for filtering reference chromosome contigs. Useful for filtering non-genomic reference contigs to speed up VCF generation.

## Set BSA Defaults  
These settings are automatically applied if not explicitly passed to automatic or BSA mode.  

| Command                     | Description                          | Possible Values                      |
| --------------------------- | ------------------------------------ | ------------------------------------ |
| `--set_loess_span`          | Set default Loess span.              | Float between 0 and 1                |
| `--set_shuffle_iterations`  | Set default shuffle iterations.      | Int (100-10000 recommended)          |
| `--set_smooth_edges_bounds` | Set default smooth edges bounds.     | int                                  |
| `--set_filter_indels`       | Set default filter indels.           | True False                           |
| `--set_filter_ems`          | Set default filter EMS.              | True False                           |
| `--set_ratio_cutoff`        | Set default ratio cutoff bound.      | Float -1 and 1 (0.1-0.3 recommended) |
| `--set_mask_snps`           | Set default mask SNPs boolean value. | True False                           |
| `--set_critical_cutoff`     | Set default critical cutoff value.   | Float (0.95-0.99 recommended)        |


Users can apply these default settings using the `phytobsa settings` command with the corresponding options. This feature is particularly useful for users who primarily work with specific reference genomes, species, or analysis methodologies, as it eliminates the need for repetitive configuration adjustments.

# Log Database Utilities

The Log Database Utilities module provides functions to interact with a log database, allowing users to easily track and retrieve runtime parameters and associated information. This is crucial for ensuring reproducibility and comparability of results across different runs of analysis or processing tasks.

| Command                     | Description                                                                    | Possible Values |
| --------------------------- | ------------------------------------------------------------------------------ | --------------- |
| `--print_analysis_log_data` | prints info related to an analysis based on analysis ULID (ulid).              | analysis ULID   |
| `--print_vcf_log_data`      | prints info related to a VCF process based on ULID.                            | vcf ULID        |
| `--print_line_name_data`    | Prints info related to a line name, including both VCF data and Analysis data. | line name       |
| `--print_core_id_data`      | Prints info related to a core ID based on core ULID.                           | core ULID       |

# Reference Database Manager

To use the Reference Database Manager, follow these steps:

1. **Fill out Configuration**: Before using the manager, ensure that the `ref_form.ini` file is filled out with the necessary configuration information.

2. **Run the Script**: Execute the script `refdb_manager.py` from the command line.

3. **Command Options**:
   - `create`: Creates a new entry in the reference database using the information provided in the `ref_form.ini` configuration file.
   - `delete`: Deletes an existing entry from the reference database. Requires specifying the `reference_name` of the entry to delete.
   - `list`: Lists all entries currently stored in the reference database. Optionally, `-ab` abbreviates long URLs for better readability.


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
```

# Footnotes

## Null Model Simulation
Null model generation begins by retrieving position and read depths in the wild type and mutant bulk for each chromosome. 

The reference and alternative reads within each chromosome are then resampled, to represent random segregation of the two 
bulks. This process breaks the link between phenotype and genotype. for each position, Coverage (C) is calculated, which is the sum of the reference and alternative reads at each position.  

A new set of reference read depths are produced from the locus-wise posterior distributions, which is produced by the following 
process.

1. **Binomial Distribution for Reference Reads:**

   Where R represents the count of reference allele occurrences in a locus in each bulk with coverage (C) reads, 
   and θ is *p*, the reference allele frequency in each bulk. 
   
   The distribution of reference reads (R) follows a binomial distribution with parameters C and *p*, where:
   
   R ~Binomial(C, *p*)

2. **Conjugate Prior Beta Distribution for Allele Frequencies:**

   The reference allele frequencies at each position are modeled using a beta distribution, serving as a conjugate prior for the binomial likelihood.

   If θ represents the reference allele frequency at each position, then θ follows a beta distribution with parameters α and β, where:
   
   *p* ~ Beta(α, β)
   
   Here, α and β are shape parameters of the beta distribution, characterizing the prior beliefs about the reference allele frequency.

3. **Update of Allele Frequencies with Bootstrapped Values:**  

   A conservative prior distribution of allele frequencies (Beta(2, 2)) undergoes an update from the bootstrapped data. The reference allele frequency *p* is then sampled from the posterior, 
   
   The product of *p* and coverage C produces the simulated reference read depth. The alternative read is then produced:
  
   simulated alternative read depth = simulated reference read depth - C

4. **Calculate Features from the simulated data**  
   G-statistics and delta-snp calculations are performed on the simulated arrays  

***Smoothing of values***   
   The Null model values are smoothed in each iteration of bootstrapping, producing a distribution of potential signals given the data. 
   From many iterations, locus-specific null models are produced, which effectively model the influence of near neighbors on the 
   smoothed values. 

***Unsmoothed values***  
  The Null model for unsmoothed values in each iteration are simply aggregated on
  a chromosome-by-chromosome basis, and used to produce percentile values for each
  value in the observed data. This allows simple filtering of those values over
  the critical cutoff. 

