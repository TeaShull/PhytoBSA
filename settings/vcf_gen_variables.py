#!/usr/bin/env python

#Edit these variables for command line only runs.  

# What is the name of your reference genome? this should be the base name of your fasta file. 
# example - Arabidopsis_thaliana.fa 
# reference_genome_name = Arabidopsis_thaliana 
reference_genome_name = 'Arabidopsis_thaliana'

#What is the snpEff database you want to use? (run snpEff databases to see list. 'grep <search> | snpEff databases' to search)
snpEff_species_db = 'Arabidopsis_thaliana'

#link to reference genome download. 
reference_genome_source = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'

# How many threads do you wish to use? set this below the max threads on your machine. 
threads_limit = '20'

#cleanup intermedia files? 
cleanup = "False"

# vcf of known background snps. 
# Drop the file into the ./references folder or give a hard coded path.
# Makes your final data much more compact and clean. 
known_snps_path = 'Arabidopsis_thaliana.vcf'

# Run gatk haplotype caller on each chromosome in parallel? The will generally run
# faster, but is a beta feature and could be buggy
call_variants_in_parallel = "True"
