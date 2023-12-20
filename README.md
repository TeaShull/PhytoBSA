# PyAtBSA
This python script analyzes and visualizes bulk segrigant analysis (BSA) data. It takes sequenced segrigant bulks as an input and outputs 
a list of likely casual polymorphisms underlying the phenotypic segrigation of the two bulks. Likely causal polymorphisms are identified using g-statistics and a
comparison of the ratio of reference to non-reference reads in each bulk. Currently, this script is optimized for running Arabidopsis BSA experiments, but can handle 
other organisms with some tinkering.   

For a simple explanation of the experimental design of BSA -> https://doi.org/10.1104/pp.17.00415 

My aim is to consolidate multiple techniques for identifiying causal EMS-induced SNPs behind *Arabidopsis thaliana* phenotypes 
using BSA. This pipeline uses:

  Delta-allele calculation. The ratio of reference read depth to total read depth in each bulk are subtracted from one another, creating a ratio that indicates phenotypic linkage.
  
  The G-statistic calculation described in the publication: https://doi.org/10.1186/s12859-020-3435-8
  
  A SNP mask generation script in order to help users mask background SNPs (Can be found in VCFsnpmask.sh). 
  This script will eventually be incorperated into the main  program, such that if a user supplies multiple lines in the same background, 
  the script will automatically generate a SNP mask out of the common snps between the generated VCFs. For now, users must manually produce their own SNP masks. 

## installation
Install and activate the conda environment from the pyatbsa_env.yml file. I highly recommend using mamba (https://mamba.readthedocs.io) to install this environment, as the environment is fairly complex and conda's environment solver is comparitivly very inefficient (conda sometimes freezes trying to resolve this environment). 

`mamba env update --name pyatbsa --file pyatbsa_env.yml`

`mamba activate pyatbsa`

or, if you are attempting to use conda (not recommended, but probably possible)

`conda env update --name pyatbsa --file pyatbsa_env.yml`

`conda activate pyatbsa`

## Usage
Put the fq.gz files you want analyzed into the /input folder. You can put multiple experiments in the folder and they will be analyzed. 
The files must be formatted as follows:

  paired-end reads - "line_1.wt.fq.gz" "line_2.wt.fq.gz" "line_1.mu.fq.gz" "line_2.mu.fq.gz"

  unpaired reads - "line.wt.fq.gz" "line.mu.fq.gz" 

 ### Variables
 Edit ./code/variables.sh to alter run conditions. The script assumes recessive polymorphisms by default. If your mutation of interest is dominant, change the mutation variable accordingly. 

 Alter the threads variable to half the number of threads available on your machine. 

 If you are analyzing an organism different than Arabidopsis, this is the place to put your reference genome and background SNPs link. In order for SNPeff to run properly, your organism needs to be named as it is in the SNPeff database. This is currently clunky to configure. Easier support for other organisms is on the roadmap. 

### Running
 
 `./PyAtBSA.py -cl` 


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



