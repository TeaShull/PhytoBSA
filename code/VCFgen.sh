#!/usr/bin/env bash
#pull in variables
source ./code/variables.sh

echo	"    Preparing reference"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# #If reference genome is missing, make a new one 
# if ! [ -f ./references/$my_species.fa ]; then
#   curl -o ./references/$my_species.fa.gz $reference_genome_source
#   gzip -d ./references/$my_species.fa.gz
# fi

# awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' ./references/$my_species.fa > ./references/$my_species.chrs.fa
fa=./references/$my_species.chrs.fa

# #creating .fai and index files
# samtools faidx $fa
# bwa index -p ./references/$my_species.chrs.fa -a is $fa

echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"    Preparing VCF for $1"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo    "    reads seem to be $2"

# #set input files as paired-end or single-read
# if [ "$2" == "paired-end" ]; then
#   fa_in_wt="./input/$1_1.wt.fq.gz ./input/$1_2.wt.fq.gz"
#   fa_in_mu="./input/$1_1.mu.fq.gz ./input/$1_2.mu.fq.gz"
# fi

# if [ "$2" == "single-read" ]; then
#   fa_in_wt="./input/$1.wt.fq.gz"
#   fa_in_wt="./input/$1.mu.fq.gz"
# fi

# #mapping
# bwa mem -t $threads -M $fa $fa_in_wt > ./output/$1_wt.sam &
# bwa mem -t $threads -M $fa $fa_in_mu > ./output/$1_mu.sam

echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"    converting sam to bam"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# samtools view -bSh -@ $threads ./output/$1_mu.sam > ./output/$1_mu.bam &
# samtools view -bSh -@ $threads ./output/$1_wt.sam > ./output/$1_wt.bam
# wait

# #fix paired end
# if [ "$2" == "paired-end" ]; then
# 	samtools fixmate ./output/$1_mu.bam ./output/$1_mu.fix.bam &
# 	samtools fixmate ./output/$1_wt.bam ./output/$1_wt.fix.bam
# fi

echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"    Sorting by coordinate"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# #Coordinate sorting
# picard SortSam I=./output/$1_mu.fix.bam O=./output/$1_mu.sort.bam SORT_ORDER=coordinate &
# picard SortSam I=./output/$1_wt.fix.bam O=./output/$1_wt.sort.bam  SORT_ORDER=coordinate
# wait

##consider using sambamba for this step. I hear it is much faster, and this seems to take longer than it should
# picard MarkDuplicates I=./output/$1_mu.sort.bam O=./output/$1_mu.sort.md.bam METRICS_FILE=./output/$1_mu.matrics.txt ASSUME_SORTED=true &
# picard MarkDuplicates I=./output/$1_wt.sort.bam O=./output/$1_wt.sort.md.bam METRICS_FILE=./output/$1_wt.matrics.txt ASSUME_SORTED=true
# wait

# #add header for gatk
# picard AddOrReplaceReadGroups I=./output/$1_mu.sort.md.bam O=./output/$1_mu.sort.md.rg.bam RGLB=$1_mu RGPL=illumina RGSM=$1_mu RGPU=run1 SORT_ORDER=coordinate &
# picard AddOrReplaceReadGroups I=./output/$1_wt.sort.md.bam O=./output/$1_wt.sort.md.rg.bam RGLB=$1_wt RGPL=illumina RGSM=$1_wt RGPU=run1 SORT_ORDER=coordinate
# wait

# picard BuildBamIndex INPUT=./output/$1_mu.sort.md.rg.bam &
# picard BuildBamIndex INPUT=./output/$1_wt.sort.md.rg.bam
# wait

#GATK HC Variant calling
gatk3 HaplotypeCaller -R $fa -I ./output/$1_mu.sort.md.rg.bam -I ./output/$1_wt.sort.md.rg.bam -o ./output/$1.hc.vcf -minReadsPerAlignStart 7 -gt_mode DISCOVERY -out_mode EMIT_ALL_SITES -writeFullFormat -stand_call_conf 10 -nct 2 -variant_index_type 1AR -variant_index_parameter 128000 -allowPotentiallyMisencodedQuals #the last argument is necessary for old sequencing results where the quality scores do not match the HC restriction: https://www.biostars.org/p/94637/; I also tried --fix_misencoded_quality_scores -fixMisencodedQuals from the same link but I received an error message. "Bad input: while fixing mis-encoded base qualities we encountered a read that was correctly encoded; we cannot handle such a mixture of reads so unfortunately the BAM must be fixed with some other tool"

############prepering for R#########################
#Exclude indels from a VCF
#gatk -R $fa -T SelectVariants --variant ./output/$1.hc.vcf -o ./output/$1.selvars.vcf --selectTypeToInclude SNP

#now make it into a table
gatk3 VariantsToTable -R $fa -V ./output/$1.hc.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -o ./output/$1.table

#snpEff
$java -jar programs/snpEff/snpEff.jar -c programs/snpEff/snpEff.config $snpEffDB -s ./output/snpEff_summary.html ./output/$1.hc.vcf > ./output/$1.se.vcf