# PyAtBSA
This python script analyzes and visualizes bulk segrigant analysis (BSA) data. It takes sequenced segrigant bulks as an input and outputs 
a list of likely casual alleles underlying the phenotypic segrigation of the two bulks. Likely causal alleles are identified using g-statistics and a
comparison of the ratio of reference to non-reference alleles in each bulk. Currently, this script is optimized for running Arabidopsis BSA experiments, but can handle 
other organisms with some tinkering.   

For a simple explanation of the experimental design of BSA -> https://doi.org/10.1104/pp.17.00415 

My aim is to consolidate multiple techniques for identifiying causal EMS-induced SNPs behind Arabidopsis Thaliana phenotypes 
using BSA. This pipeline uses:

  Delta-allele calculation. The ratio of reference read depth to total read depth in each bulk are subtracted from one another, creating a ratiometrically derived linkage map.
  
  The G-statistic calculation described in the publication: https://doi.org/10.1186/s12859-020-3435-8
  
  A SNP mask generation script in order to help users mask background SNPs (Can be found in VCFsnpmask.sh). 
  This script will eventually be incorperated into the main  program, such that if a user supplies multiple lines in the same background, 
  the script will automatically generate a SNP mask out of the common snps between the generated VCFs. For now, users must manually produce their own SNP masks. 

## installation
Install and activate the conda environment from the pyatbsa_conda_env_file.yaml file included in the /code directory. I highly recommend using mamba (https://mamba.readthedocs.io) to install this environment, as the environment is fairly complex and conda's environment solver is comparitivly very inefficient (conda sometimes freezes trying to resolve this environment). 

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

 `./PyAtBSA.py` 

## Output

Will output a linkage map with lowess smoothed ratios and SNPS in the top 98th percentile g-statistic highlighted in red
![4743 BSA_linkage_map](https://github.com/TeaShull/PyAtBSA/assets/125574642/e5979fc8-7310-4a00-acec-5bde5c8aaf5c)

Also, the script will output a graph of the calculated g-statistics, also lowess smoothed. 

![4743 BSA_g-stat](https://github.com/TeaShull/PyAtBSA/assets/125574642/2f1a7ff2-0eda-4258-bbbc-d6a2895997c4)


finally, a table of the likely candidates will be created, which should contain your mutation of interest. 



