#!/usr/bin/env bash
#Edit these variables if you are operating with a different reference genome. This has only been tested with arabidopsis so far, but I'm sure it can be used for other organisms. 

#species name you want used in the script
my_species='Arabidopsis_thaliana'

#is the mutation 'recessive' or 'dominant'?
mutation='recessive'

#reference genomea and genome features file (gff) to be sourced. You can download seperately and this will be skipped if you stick to the $my_species.fa and $my_species.gff format for the reference genome. 
reference_genome_source='ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'
gff_file_source='ftp://ftp.ensemblgenomes.org/pub/plants/release-56/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gff3.gz'

#set how many threads you would like the mapping processes to use. Please note that sometimes 2 processes are run simultaniously, so half the number of cores you wish to use. 
# for example, if you have 24 threads available, setting threads to 10 should be OK provided you aren't doing anything else CPU intensive.
threads='8'

# set this variable to False if you want all files output by this script to be saved in archive. This takes up quite some space, but good if you want to keep your sorted Bams and so on for other workflows.
cleanup=True



knownsnps="./references/Arabidopsis_thaliana.vcf"
