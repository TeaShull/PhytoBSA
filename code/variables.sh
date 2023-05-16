#!/usr/bin/env bash
#Edit these variables if you are operating with a different reference genome. This has only been tested with arabidopsis so far, but I'm sure it can be used for other organisms. 

#species name you want used in the script
my_species='Arabidopsis_thaliana'

#reference genome to be sourced. You can download seperately and this will be skipped if you stick to the $my_species.fa format for the reference genome. 
reference_genome_source='ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'

#set how many threads you would like the mapping processes to use. Please note that sometimes 2 processes are run simultaniously, so half the number of cores you wish to use. 
# for example, if you have 24 threads available, setting threads to 10 should be OK provided you aren't doing anything else CPU intensive.
threads='8'
