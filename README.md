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


## Output

Will output 6 plots: Calculated ratios, G statistics, and ratio-scaled G statistics as well as their corresponding lowess smoothed graphs. Red dashed lines represent calculated empirical cutoffs for likely candidate genes. Causal SNPs nearly always appear clearly at the top of a the lowess smoothed ratio-scaled g statistic graph. 

Example: 474-3 in https://doi.org/10.1104/pp.17.00415. Pipeline correctly identifies early stop codon in *SHR*

SRA Runs - SRR5029628(474_3_wt); SRR5029636 (474_3_mut) 

Ratio Scaled G-statistics
![RS_G_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/7a73e741-4722-4a1b-86be-1cc10b185535)

Ratio Scale G-statistics, Lowess smoothed
![RS_G_yhat_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/b7e8dd00-af16-42c7-a396-cad954f27de9)

G-statistics
![GS_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/0214466f-749e-4f3f-910b-50560840b647)

G-statistics, Lowess smoothed
![GS_yhat_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/3217b417-aac3-4e72-993c-71995421a01a)

Ratio
![SNP_ratio_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/21b42e24-5dbb-4c5b-b32d-ff7d2072b42e)

Ratio, Lowess smooted 
![ratio_yhat_4773](https://github.com/TeaShull/PyAtBSA/assets/125574642/074361b6-3c43-4d12-9976-1a7e530e2535)
finally, a table of the likely candidates will be created, which should contain your mutation of interest.
