# PyAtBSA
This is my attempt to write a python script to analyze and visualize bulk segrigant analysis (BSA) data, taking sequenced segrigant bulks as an input, and outputing 
a list of putative casual alleles behind the phenotypic segrigation of the two bulks. Putative causal alleles are identified using g-statistics and a
comparison of the ratio of reference to non-reference alleles in each bulk. Currently, this script is optimized for running Arabidopsis BSA experiments, but can handle 
other organisms with some tinkering.   

For a good, simple explanation of the underlying experimental design of BSA analysis -> https://doi.org/10.1104/pp.17.00415 

To use, install and activate the conda environment from the pyatbsa_conda_env_file.yaml file included in the /code directory. I highly recommend using mamba (https://mamba.readthedocs.io) to install this environment, as the environment is fairly complex and conda is comparitivly inefficient (conda sometimes freezes trying to resolve this environment). 

Put the fq.gz files you want analyzed into the /input folder. You can put multiple experiments in the folder and they will be analyzed. 
The files must be formatted as follows:

  paired-end reads - "line_1.wt.fq.gz" "line_2.wt.fq.gz" "line_1.mu.fq.gz" "line_2.mu.fq.gz"

  unpaired reads - "line.wt.fq.gz" "line.mu.fq.gz" 

  run PyAtBSA.py 

  
Conceptually, there is nothing new here. My hope is to consolidate multiple techniques for identifiying causal EMS-induced SNPs behind Arabidopsis Thaliana phenotypes 
using BSA. My implementation of this pipeline uses:

  the delta-allele calculation described in the publication: https://doi.org/10.1104/pp.17.00415 (also linked earlier for a quick explanation of BSA experimental design)
  
  The G-statistic calculation described in the publication: https://doi.org/10.1186/s12859-020-3435-8
  
  A SNP mask generation script in order to help users mask background SNPs (Can be found in VCFsnpmask.sh). 
  This script will eventually be incorperated into the main  program, such that if a user supplies multiple lines in the same background, 
  the script will automatically generate a SNP mask out of the common snps between the generated VCFs. For now, users must manually produce their own SNP masks. 

output: 
Will output a linkage map with lowess smoothed ratios and SNPS in the top 98th percentile g-statistic highlighted in red
![4743 BSA_linkage_map](https://github.com/TeaShull/PyAtBSA/assets/125574642/e5979fc8-7310-4a00-acec-5bde5c8aaf5c)

Also, the script will output a graph of the calculated g-statistics, also lowess smoothed. 

![4743 BSA_g-stat](https://github.com/TeaShull/PyAtBSA/assets/125574642/2f1a7ff2-0eda-4258-bbbc-d6a2895997c4)

finally, a table of the likely candidates will be created, which should contain your mutation of interest. 
