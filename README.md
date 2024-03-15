[01HS1PNQVNJ32VQH6V7GYYT2Q0-_4743_likely_candidates.csv](https://github.com/TeaShull/PhytoBSA/files/14619674/01HS1PNQVNJ32VQH6V7GYYT2Q0-_4743_likely_candidates.csv)- [PhytoBSA](#phytobsa)
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
- `--set_snpmask_path`: Set default SNP mask VCF path. VCF should contain known background SNPs.
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

Will output 6 plots: Calculated ratios, G statistics, and ratio-scaled G statistics as well as their corresponding lowess smoothed graphs. Red dashed lines represent calculated empirical cutoffs for likely candidate genes. Causal SNPs nearly always appear clearly at the top of a the lowess smoothed ratio-scaled g statistic graph. 

Example: 474-3 in https://doi.org/10.1104/pp.17.00415. Pipeline correctly identifies early stop codon in *SHR*

SRA Runs - SRR5029628(474_3_wt); SRR5029636 (474_3_mut) 


Ratio Scale G-statistics, Lowess smoothed
![01HS1PNQVNJ32VQH6V7GYYT2Q0-_4743_ratio_yhat](https://github.com/TeaShull/PhytoBSA/assets/125574642/709423e3-4313-4595-894a-3dc43ea89ee2)


G-statistics, Lowess smoothed
![01HS1PNQVNJ32VQH6V7GYYT2Q0-_4743_g_s_yhat](https://github.com/TeaShull/PhytoBSA/assets/125574642/b3c8ca0b-a799-4252-9fc5-281edf9e5dfa)

Ratio, Lowess smooted 
![01HS1PNQVNJ32VQH6V7GYYT2Q0-_4743_ratio_yhat](https://github.com/TeaShull/PhytoBSA/assets/125574642/3158aac9-60da-4ce9-9e70-c099e80c1082)

Finally, a list of the likely candidates will be produced: 
[Uchrom,pos,ref,alt,gene,snpeffect,snpvariant,snpimpact,mu:wt_gtpred,mu_ref,mu_alt,wt_ref,wt_alt,ratio,G_S,RS_G,ratio_yhat,G_S_yhat,RS_G_yhat,ratio_percentile,G_S_percentile,RS_G_percentile,ratio_yhat_percentile,ratio_yhat_null_1,ratio_yhat_null_5,ratio_yhat_null_25,ratio_yhat_null_50,ratio_yhat_null_75,ratio_yhat_null_95,ratio_yhat_null_99,G_S_yhat_percentile,G_S_yhat_null_1,G_S_yhat_null_5,G_S_yhat_null_25,G_S_yhat_null_50,G_S_yhat_null_75,G_S_yhat_null_95,G_S_yhat_null_99,RS_G_yhat_percentile,RS_G_yhat_null_1,RS_G_yhat_null_5,RS_G_yhat_null_25,RS_G_yhat_null_50,RS_G_yhat_null_75,RS_G_yhat_null_95,RS_G_yhat_null_99
4,17692534,C,T,SHR:AT4G09525:AT4G37660:AT4G09515:AT4G09535:NAGS2:NAGS2,stop_gained:upstream_gene_variant:upstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant,p.Arg222*:NaN:NaN:NaN:NaN:NaN:NaN,HIGH:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,0,43,38,22,0.6333333333333333,56.768346969985956,35.95328641432444,0.630551221952931,58.42499647683181,34.486293495035774,0.833184090909091,0.8722443181818181,0.8644704545454546,1.0,0.1753453439237463,0.216401551210226,0.29398731654782656,0.33918376223823254,0.38503800660098003,0.43958436786988614,0.47522407464150185,1.0,8.423360630584567,10.041833926934688,13.903322522664505,17.03201924452976,19.908784251817274,24.99579301961561,28.951442210429857,1.0,1.6514398916474915,2.7577144519049357,4.934432286986327,6.721194406057113,8.526447613075478,11.878846023763652,14.261317244559898
4,18088198,C,T,AT4G38760:CYP18-3:AT4G09735,splice_acceptor_variant&intron_variant:upstream_gene_variant:downstream_gene_variant,NaN:NaN:NaN,HIGH:MODIFIER:MODIFIER,1/1:0/1,6,62,34,13,0.6351689612015019,52.58087004351695,33.397736604611836,0.6269975721699516,57.94526668260699,34.38384383647178,0.8355613636363637,0.8495943181818182,0.8505522727272727,1.0,0.16940188536320283,0.21435240590806054,0.294867041049284,0.3394862847804554,0.38510863653612226,0.43817960470463346,0.47844948217761857,1.0,8.425138136688421,10.112402618199116,14.088602750294669,17.128751953760872,20.24751667910326,25.328604533045333,29.418299202806057,1.0,1.7000036122267812,2.8518958678204083,4.974842912599637,6.822431923330111,8.689267896204916,11.995133635223747,14.418579305927466
4,16694522,C,T,AT4G35070:AT4G35070:AT4G35070:AT4G35080:AT4G35080:AT4G35080,missense_variant:5_prime_UTR_premature_start_codon_gain_variant:5_prime_UTR_variant:upstream_gene_variant:upstream_gene_variant:upstream_gene_variant,p.Pro12Leu:NaN:NaN:NaN:NaN:NaN,MODERATE:LOW:MODIFIER:MODIFIER:MODIFIER:MODIFIER,0/1:0/1,4,52,32,7,0.7490842490842491,60.54617118132328,45.354183174287954,0.5533176859330421,47.2981483312369,30.63543990812101,0.9313454545454546,0.8899784090909091,0.9077329545454546,0.9976,0.12932149303293555,0.1993289699179438,0.27790728981518875,0.33940125330609083,0.39085565572265946,0.4567145401237477,0.5068338550601705,1.0,5.971954369025904,8.733888945008374,12.777922131168648,16.05701956901632,20.471275505126776,25.74560767478886,28.23874415329986,1.0,1.4240159386553763,2.087855189011417,4.53551372832199,6.357183742064739,8.345875366219513,12.534811593452332,14.789821961014606
4,17095776,C,T,FPP5:FPP5:RPL8C:AT4G36110:AT4G36140:AT4G36140,missense_variant:missense_variant:upstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant,p.Ala201Thr:p.Ala201Thr:NaN:NaN:NaN:NaN,MODERATE:MODERATE:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,0,50,30,21,0.5882352941176471,53.77752970857016,31.63384100504127,0.6354010879939035,59.121349244565565,34.66340928328275,0.7896170454545455,0.8562125,0.8398511363636364,1.0,0.14646782672472322,0.20959349748085912,0.2839789357491832,0.34227620422575056,0.3899266570223838,0.45051957020687927,0.4909492117354908,1.0,7.377105011365602,9.29507486592882,13.04784110434131,16.337717877420467,20.34240944746218,25.061524457808577,27.109449389221464,1.0,1.583049776642864,2.3553801485214962,4.627966145262357,6.56022920252839,8.171456633686342,12.109762012287142,14.242830949291182
4,17423140,C,T,MAPKKK21:AT4G36945:AT4G36960:AT4G36960:NAPRT1,synonymous_variant:upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:downstream_gene_variant,p.Ser235Ser:NaN:NaN:NaN:NaN,LOW:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,0,51,31,16,0.6595744680851063,62.03377937197014,40.91589703257605,0.6348349947210825,59.036061274527725,34.63988475960665,0.8580829545454546,0.8964306818181819,0.8888659090909091,1.0,0.18035894315911102,0.21723644439074638,0.29264364740942134,0.3400283038401847,0.38634857301497183,0.43923143137943177,0.47376925940041736,1.0,8.165631876638452,10.164928216183853,13.854917290630691,16.883366251293605,19.957876746906024,24.834779967882024,28.500518073863923,1.0,1.610032952521331,2.7039991713020077,4.830857373618891,6.654941648968704,8.38988345975188,11.815389041383465,14.085074550302913
4,17626553,C,T,RNC1:AT4G37483:CYCB1-1:PER50:AT4G37480:AT4G37500,synonymous_variant:upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:downstream_gene_variant:downstream_gene_variant,p.Leu440Leu:NaN:NaN:NaN:NaN:NaN,LOW:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,1,67,29,16,0.6297385620915034,61.79239950102303,38.91305680995797,0.6316421867345468,58.57558993263733,34.521564524427085,0.8298136363636364,0.8952170454545455,0.8796454545454545,1.0,0.18009405110216298,0.21894777956561096,0.2921917378266127,0.33987130583454217,0.38586607288707586,0.4398695401996722,0.47376114587116946,1.0,8.257015423027356,10.159769789939633,13.846253953547984,16.888655579365043,19.91632985129993,24.87188688634488,28.574838918954367,1.0,1.6169152731225809,2.7021775285022653,4.864158953278368,6.656754769782987,8.403677637533246,11.825695851909133,14.114439258182506
4,17122658,C,T,AT4G36180:AT4G36190:AT4G36195:AT4G36195:AT4G36195:AT4G36170,synonymous_variant:upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:downstream_gene_variant,p.Glu347Glu:NaN:NaN:NaN:NaN:NaN,LOW:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,0,51,25,15,0.625,54.07315394252994,33.79572121408121,0.6350325984453542,59.065697929718525,34.64800582566625,0.8249477272727272,0.8575340909090909,0.8534159090909091,1.0,0.18042582302636928,0.21620056508036412,0.2930551272601705,0.3398434274983696,0.3863058523998443,0.4388444108061175,0.4738678048529425,1.0,8.12236135752919,10.168160410190715,13.827038587182454,16.89086820832858,19.919762457667016,24.820457891484764,28.45551093631098,1.0,1.6058543352508143,2.705094859437859,4.817201456561163,6.669277937210818,8.373316738537609,11.809173315804077,14.12192388211571
4,17765007,C,T,HAT22:AT4G09665:MYB87-AT4G09665,upstream_gene_variant:downstream_gene_variant:intergenic_region,NaN:NaN:NaN,MODIFIER:MODIFIER:MODIFIER,1/1:0/1,0,41,37,16,0.6981132075471698,61.10347080877011,42.65713999857536,0.6299059754670269,58.336935253135216,34.466444040743546,0.8916318181818181,0.8921954545454546,0.8968340909090909,1.0,0.17382137857641553,0.2168874546423929,0.2944365423873564,0.3394003404293937,0.38444650921925716,0.4397682444908039,0.4757295973434313,1.0,8.459380831186747,9.992878932014579,13.882945057751828,17.09193480848744,19.965876968362675,24.95890911808298,29.02326702814274,1.0,1.6623731742699877,2.804554812896439,4.924223379588684,6.769738846745954,8.538314744364467,11.896246628609514,14.307113577611753
4,18547798,C,T,AT4G40000:AT4G40000:RABA4B:SRK2F:AT4G40011:AT4G40000-SRK2F,upstream_gene_variant:upstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:intergenic_region,NaN:NaN:NaN:NaN:NaN:NaN,MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,3,56,39,15,0.6713747645951036,61.59553866902681,41.35369027402648,0.6256133497471504,57.759763863161005,34.34718967250326,0.8694102272727273,0.8945704545454546,0.8910897727272727,1.0,0.17276912268242367,0.2223210571752072,0.29264077604219757,0.33694834234240345,0.38448230584699,0.44106264143931356,0.4822846604588544,1.0,8.564487063704789,10.345317786184227,14.148011236722603,17.123045745332917,20.426738199459788,25.707229996567314,29.574719698811023,1.0,1.5826644049786578,2.8727124015359125,5.017300428909985,6.809027093761269,8.724557732638873,12.232521776704068,14.766288109378825
4,17510130,C,T,AT4G37210:AT4G37220:AT4G37210:HHO5:AT4G37190:HHO5:HHO5:AT4G37190:HCF164,upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:intron_variant,NaN:NaN:NaN:NaN:NaN:NaN:NaN:NaN:NaN,MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,0,56,32,22,0.5925925925925926,59.654694429579514,35.35093003234341,0.6324172152128651,58.68435712305183,34.54809942851402,0.7928909090909091,0.8861693181818182,0.861028409090909,1.0,0.18034072731679066,0.21751785980099742,0.292517721023766,0.34007845010801996,0.38637587294767134,0.4393364738221473,0.4738237339921542,1.0,8.180636350371781,10.164064361111471,13.862482064406061,16.8911554310538,19.96072124727005,24.840866618840813,28.512742456833504,1.0,1.6111664996717288,2.703700620760651,4.834556474767554,6.657869201894096,8.39035446522852,11.81708060596512,14.089687953815822
4,17044526,C,T,CSP1:CSP1:CSP1:AT4G36010:AT4G36032:AT4G36010:AT4G36032:ARO3:ATJ11,5_prime_UTR_variant:5_prime_UTR_variant:5_prime_UTR_variant:upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:downstream_gene_variant:downstream_gene_variant,NaN:NaN:NaN:NaN:NaN:NaN:NaN:NaN:NaN,MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,0,60,30,22,0.5769230769230769,59.31831183406111,34.2221029811891,0.5813754657449812,52.12476695089484,33.78620150958883,0.7762465909090909,0.8840886363636363,0.8559386363636363,1.0,0.14500224611447704,0.20924794297763005,0.28303818867097,0.3430257333803568,0.3896154963392666,0.4506016830649909,0.4924799728582878,1.0,7.239345772964364,9.1528511428539,12.996295216026232,16.367333383595405,20.281239135604345,25.087807147253045,27.12816349182707,1.0,1.5912438502437938,2.345403714399253,4.598204940553726,6.537776966699883,8.164734852482361,12.086478206131316,14.290776609931015
4,16374845,C,T,CYCLASE1:SEP1:ASK11:AT4G34215:PGDH1,upstream_gene_variant:upstream_gene_variant:upstream_gene_variant:downstream_gene_variant:intron_variant,NaN:NaN:NaN:NaN:NaN,MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,5,71,32,29,0.45880069025021575,38.539909373554835,17.68213702276772,0.5308435175648264,43.38867938918407,26.770822091866613,0.6352204545454545,0.7603647727272728,0.7236443181818182,0.9832,0.0604508718307373,0.1636615407495102,0.25753501124704875,0.3401834301548936,0.41252167157273956,0.4831489647012214,0.5567754938056807,0.9976,3.826516376438784,6.257652744755051,12.299337600672349,17.243903688975735,22.696315315496534,30.843824879036365,38.83559870282994,1.0,-0.3363807379723674,1.1365647392156524,3.8496953293181333,6.8931378420834255,10.067256963053605,14.978731199712076,20.797261243056283
4,17967896,C,T,AT4G38350:AT4G38350:AT4G38370:DTX45:DTX45:DTX45:DTX45:DTX45:LAZ1:LAZ1,upstream_gene_variant:upstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:intron_variant:intron_variant,NaN:NaN:NaN:NaN:NaN:NaN:NaN:NaN:NaN:NaN,MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,0/1:0/1,4,48,26,20,0.48829431438127086,29.540342661680434,14.424381366573053,0.6291747819589235,58.23780516715762,34.44472303970888,0.6715579545454545,0.6826306818181819,0.6827943181818181,1.0,0.1453316992974905,0.20871325042876057,0.28269062014603846,0.3425516956127684,0.38919078271928803,0.45055137571292925,0.4930438868246531,1.0,7.188764819861465,9.141185832893592,12.984325175374211,16.38273855659554,20.29312534758942,25.16276443176957,27.13818871372003,1.0,1.5943074931116956,2.3418800456020654,4.5984791626256705,6.517105879715328,8.174530197540005,12.110922431789836,14.308885591621898
4,15685807,C,T,AT4G32510:AKT5:SHM3:AT4G32510:SHM3,3_prime_UTR_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant:downstream_gene_variant,NaN:NaN:NaN:NaN:NaN,MODIFIER:MODIFIER:MODIFIER:MODIFIER:MODIFIER,1/1:0/1,2,35,31,27,0.48042870456663556,27.012482311478028,12.977571884032546,0.48173535830883424,35.49567692798892,21.19646457823179,0.6608045454545455,0.6566988636363637,0.6618522727272728,0.9223,0.0501943919156405,0.16517389441361555,0.2884596320300151,0.34986546599167373,0.4082616527333114,0.4963321995426252,0.5596650678183103,0.988,6.74375110806397,7.884287331152621,12.473559474499725,16.62249051408014,21.27705285556521,29.503611097895202,36.94801928195035,0.9904,0.813827584427881,1.8035905687449643,4.512074900956163,6.368622968581495,9.77465956739545,14.86223215240119,20.6082323784598
ploading 01HS1PNQVNJ32VQH6V7GYYT2Q0-_4743_likely_candidates.csv…]()

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



