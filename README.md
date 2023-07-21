# PyAtBSA
This python script analyzes and visualizes bulk segrigant analysis (BSA) data. It takes sequenced segrigant bulks as an input and outputs 
a list of likely casual alleles underlying the phenotypic segrigation of the two bulks. Likely causal alleles are identified using g-statistics and a
comparison of the ratio of reference to non-reference alleles in each bulk. Currently, this script is optimized for running Arabidopsis BSA experiments, but can handle 
other organisms with some tinkering.   

For a simple explanation of the experimental design of BSA -> https://doi.org/10.1104/pp.17.00415 

##installation
Install and activate the conda environment from the pyatbsa_conda_env_file.yaml file included in the /code directory. I highly recommend using mamba (https://mamba.readthedocs.io) to install this environment, as the environment is fairly complex and conda is comparitivly inefficient (conda sometimes freezes trying to resolve this environment). 

'mamba env update --name pyatbsa --file pyatbsa_env.yml'

or, if you are attempting to use conda (not recommended, but probably possible)

'conda env update --name pyatbsa --file pyatbsa_env.yml'

##Usage
Put the fq.gz files you want analyzed into the /input folder. You can put multiple experiments in the folder and they will be analyzed. 
The files must be formatted as follows:

  paired-end reads - "line_1.wt.fq.gz" "line_2.wt.fq.gz" "line_1.mu.fq.gz" "line_2.mu.fq.gz"

  unpaired reads - "line.wt.fq.gz" "line.mu.fq.gz" 

 './PyAtBSA.py' 

##Output
My aim is to consolidate multiple techniques for identifiying causal EMS-induced SNPs behind Arabidopsis Thaliana phenotypes 
using BSA. This pipeline uses:

  the delta-allele calculation described in the publication: https://doi.org/10.1104/pp.17.00415 (also linked earlier for a quick explanation of BSA experimental design)
  
  The G-statistic calculation described in the publication: https://doi.org/10.1186/s12859-020-3435-8
  
  A SNP mask generation script in order to help users mask background SNPs (Can be found in VCFsnpmask.sh). 
  This script will eventually be incorperated into the main  program, such that if a user supplies multiple lines in the same background, 
  the script will automatically generate a SNP mask out of the common snps between the generated VCFs. For now, users must manually produce their own SNP masks. 





