- [PhytoBSA](#phytobsa)
  - [Installation](#installation)
    - [Environment installation](#environment-installation)
    - [Setup directories](#setup-directories)
  - [Usage](#usage)
    - [Variables](#variables)
    - [Running](#running)
    - [Command Line Arguments](#command-line-arguments)
      - [BSA Analysis Options](#bsa-analysis-options)
      - [VCF Generation Options](#vcf-generation-options)
      - [Reference Name](#reference-name)
    - [Output](#output)
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

### Running
 
 `./phytobsa.py -cl` 

### Command Line Arguments

#### BSA Analysis Options
- `-ls`, `--loess_span`: 
  - Type: float
  - Default: None
  - Description: Influences smoothing parameters.

- `-si`, `--shuffle_iterations`: 
  - Type: int
  - Default: None
  - Description: Iterations of bootstrapping during empirical cutoff calculations. Below 1000 can yield inconsistent results.

- `-sb`, `--smooth_edges_bounds`: 
  - Type: int
  - Default: None
  - Description: Number of mirrored datapoints at chromosome edges to correct for loess edge bias. Increase if edge bias seems high.

- `-fin`, `--filter_indels`: 
  - Type: str
  - Default: None
  - Description: Filter out insertion-deletion mutations.

- `-fems`, `--filter_ems`: 
  - Type: str
  - Default: None
  - Description: Filter results to only include mutations likely to arise from EMS treatment.

- `-rco`, `--ratio_cutoff`: 
  - Type: float
  - Default: None
  - Description: Used to filter results based on a ratio cutoff number. Increase to 0.2 or 0.3 if there is a lot of noise at lower ratio bounds.

- `-msk`, `--mask_snps`: 
  - Type: bool
  - Default: None
  - Description: Set to true if you have a snpmask file configured and would like to mask known SNPs in your analysis.

- `-cc`, `--critical_cutoff`: 
  - Type: float
  - Default: None
  - Description: Set the critical cutoff value for what is considered a significant polymorphism.

- `-m`, `--method`: 
  - Type: str
  - Default: None
  - Description: Set the method of generating the null hypothesis. Either simulate or bootstrap.

#### VCF Generation Options

- `-p`, `--call_variants_in_parallel`: 
  - Type: bool
  - Default: None
  - Description: Run GATK haplotype caller in parallel.

- `-c`, `--cleanup`: 
  - Type: bool
  - Default: None
  - Description: If true, intermediate files will be deleted. False for troubleshooting and archiving files.

- `-cft`, `--cleanup_filetypes`: 
  - Type: list
  - Default: None
  - Description: Filetypes to clean out after VCF generation is complete. Format should be ['file_suffix', exc]. Example: ['.tmp', '*.metrics']

- `-ocp`, `--omit_chrs_patterns`: 
  - Type: list
  - Default: None
  - Description: Header patterns to omit from reference chromosomes. Useful for removing >mt (mitochondrial) and other unneeded reference sequences.


#### Reference Name

- `-r`, `--reference_name`: 
  - Type: str
  - Required: False
  - Description: Name of the reference genome.


### Output

Will output 6 plots: Calculated ratios, G statistics, and ratio-scaled G statistics as well as their corresponding lowess smoothed graphs. Red dashed lines represent calculated empirical cutoffs for likely candidate genes. Causal SNPs nearly always appear clearly at the top of a the lowess smoothed ratio-scaled g statistic graph. 

Example: 474-3 in https://doi.org/10.1104/pp.17.00415. Pipeline correctly identifies early stop codon in *SHR*

SRA Runs - SRR5029628(474_3_wt); SRR5029636 (474_3_mut) 

Ratio Scaled G-statistics
![RS_G_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/7a73e741-4722-4a



