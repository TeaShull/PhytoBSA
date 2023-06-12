#!/usr/bin/env bash
#pull in variables
source ./code/variables.sh

echo	"    Preparing reference"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

if [ ! -d "./output/archive" ]; then
	mkdir ./output/archive
fi 

if [ ! -d "./output/archive/${1}" ]; then
	mkdir ./output/archive/${1}
fi

if [ ! -d "./output/archive/VCFs" ]; then
	mkdir ./output/archive/VCFs
fi

if [ ! -d "./output/archive/${1}/snpEff" ]; then
	mkdir ./output/archive/${1}/snpEff
fi

#If reference genome is missing, make a new one 
if ! [ -f ./references/$my_species.fa ]; then
  curl -o ./references/$my_species.fa.gz $reference_genome_source
  gzip -d ./references/$my_species.fa.gz
fi

awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' ./references/$my_species.fa > ./references/$my_species.chrs.fa
fa=./references/$my_species.chrs.fa

# #creating .fai and index files
# samtools faidx $fa
# bwa index -p ./references/$my_species.chrs.fa -a is $fa

#create dictory for gatk haplotype caller
<<<<<<< HEAD
picard CreateSequenceDictionary -R .fa -O ./references/$my_species.chrs.dict
=======
if [ ! -d "./references/$my_species.chrs.dict" ]; then
	picard CreateSequenceDictionary -R .fa -O ./references/$my_species.chrs.dict
fi
>>>>>>> ee0054a (add plot production to .py script, add archive management and cleanup to vcfgen.sh script)

echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"    Preparing VCF for $1"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo    "    reads seem to be $2"

#set input files as paired-end or single-read
if [ "$2" == "paired-end" ]; then
  fa_in_wt="./input/$1_1.wt.fq.gz ./input/$1_2.wt.fq.gz"
  fa_in_mu="./input/$1_1.mu.fq.gz ./input/$1_2.mu.fq.gz"
fi

if [ "$2" == "single-read" ]; then
  fa_in_wt="./input/$1.wt.fq.gz"
  fa_in_wt="./input/$1.mu.fq.gz"
fi

#mapping
bwa mem -t $threads -M $fa $fa_in_wt > ./output/$1_wt.sam &
bwa mem -t $threads -M $fa $fa_in_mu > ./output/$1_mu.sam

echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"    converting sam to bam"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

samtools view -bSh -@ $threads ./output/$1_mu.sam > ./output/$1_mu.bam &
samtools view -bSh -@ $threads ./output/$1_wt.sam > ./output/$1_wt.bam
wait

#fix paired end
if [ "$2" == "paired-end" ]; then
	samtools fixmate ./output/$1_mu.bam ./output/$1_mu.fix.bam &
	samtools fixmate ./output/$1_wt.bam ./output/$1_wt.fix.bam
fi

echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"    Sorting by coordinate"
echo	"    >=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

#Coordinate sorting
picard SortSam I=./output/$1_mu.fix.bam O=./output/$1_mu.sort.bam SORT_ORDER=coordinate &
picard SortSam I=./output/$1_wt.fix.bam O=./output/$1_wt.sort.bam  SORT_ORDER=coordinate
wait

#consider using sambamba for this step. I hear it is much faster, and this seems to take longer than it should
picard MarkDuplicates I=./output/$1_mu.sort.bam O=./output/$1_mu.sort.md.bam METRICS_FILE=./output/$1_mu.matrics.txt ASSUME_SORTED=true &
picard MarkDuplicates I=./output/$1_wt.sort.bam O=./output/$1_wt.sort.md.bam METRICS_FILE=./output/$1_wt.matrics.txt ASSUME_SORTED=true
wait

#add header for gatk
picard AddOrReplaceReadGroups I=./output/$1_mu.sort.md.bam O=./output/$1_mu.sort.md.rg.bam RGLB=$1_mu RGPL=illumina RGSM=$1_mu RGPU=run1 SORT_ORDER=coordinate &
picard AddOrReplaceReadGroups I=./output/$1_wt.sort.md.bam O=./output/$1_wt.sort.md.rg.bam RGLB=$1_wt RGPL=illumina RGSM=$1_wt RGPU=run1 SORT_ORDER=coordinate
wait

picard BuildBamIndex INPUT=./output/$1_mu.sort.md.rg.bam &
picard BuildBamIndex INPUT=./output/$1_wt.sort.md.rg.bam
wait

<<<<<<< HEAD
#GATK HC Variant calling
#gatk HaplotypeCaller -R $fa -I ./output/$1_mu.sort.md.rg.bam -I ./output/$1_wt.sort.md.rg.bam -O ./output/$1.hc.vcf -output-mode EMIT_ALL_CONFIDENT_SITES --native-pair-hmm-threads 20 #the last argument is necessary for old sequencing results where the quality scores do not match the HC restriction: https://www.biostars.org/p/94637/; I also tried --fix_misencoded_quality_scores -fixMisencodedQuals from the same link but I received an error message. "Bad input: while fixing mis-encoded base qualities we encountered a read that was correctly encoded; we cannot handle such a mixture of reads so unfortunately the BAM must be fixed with some other tool"

############prepering for R#########################
#Exclude indels from a VCF
#gatk -R $fa -T SelectVariants -V ./output/$1.hc.vcf -O ./output/$1.selvars.vcf --selectTypeToInclude SNP

#now make it into a table

# #snpEff
#snpEff $my_species -s ./output/snpEff_summary.html ./output/$1.hc.vcf > ./output/$1.se.vcf
=======
GATK HC Variant calling
gatk HaplotypeCaller -R $fa -I ./output/$1_mu.sort.md.rg.bam -I ./output/$1_wt.sort.md.rg.bam -O ./output/$1.hc.vcf -output-mode EMIT_ALL_CONFIDENT_SITES --native-pair-hmm-threads 20 #the last argument is necessary for old sequencing results where the quality scores do not match the HC restriction: https://www.biostars.org/p/94637/; I also tried --fix_misencoded_quality_scores -fixMisencodedQuals from the same link but I received an error message. "Bad input: while fixing mis-encoded base qualities we encountered a read that was correctly encoded; we cannot handle such a mixture of reads so unfortunately the BAM must be fixed with some other tool"

###########prepering for R#########################
Exclude indels from a VCF
gatk -R $fa -T SelectVariants -V ./output/$1.hc.vcf -O ./output/$1.selvars.vcf --selectTypeToInclude SNP

#now make it into a table

#snpEff
snpEff $my_species -s ./output/archive/${1}/snpEff/${1}_snpEff_summary.html ./output/$1.hc.vcf > ./output/$1.se.vcf
>>>>>>> ee0054a (add plot production to .py script, add archive management and cleanup to vcfgen.sh script)

#extract snpEFF data and variant information into a table, remove repetative NaN's and retain only those polymorphisms likely to arise from EMS.
SnpSift extractFields -s ":" -e "NaN" ./output/$1.se.vcf CHROM POS  REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].IMPACT" "GEN[*].GT" "GEN[$1_mu].AD" "GEN[$1_wt].AD" > ./output/$1.table

grep -e $'G\tA' -e $'C\tT' -e $'A\tG' -e $'T\tC' ./output/$1.table > ./output/$1.ems.table.tmp
sed -i 's/NaN://g' ./output/$1.ems.table.tmp

<<<<<<< HEAD
#if [ "$mutation" = 'recessive' ]; then 
	grep -F -e '1/1:0/1' -e '0/1:0/0' ./output/$1.ems.table.tmp > ./output/$1.ems.table
#else 
#grep -e $'0/1:0/0' -e '1/1:0/0' ./output/$1.ems.table.tmp > ./output/$1.ems.table
#fi

rm ./output/*.tmp
=======
if [ "$mutation" = 'recessive' ]; then 
	grep -F -e '1/1:0/1' -e '0/1:0/0' ./output/$1.ems.table.tmp > ./output/$1.ems.table
else 
grep -e $'0/1:0/0' -e '1/1:0/0' ./output/$1.ems.table.tmp > ./output/$1.ems.table
fi

>>>>>>> ee0054a (add plot production to .py script, add archive management and cleanup to vcfgen.sh script)
awk -i inplace -F'\t' -vOFS='\t' '{ gsub(",", "\t", $9) ; gsub(",", "\t", $10) ; gsub(",", "\t", $11) ; print }' ./output/$1.ems.table

#remove complex genotypes
awk -i inplace -F'\t' 'NF==13' ./output/$1.ems.table

#remove non-numeric chromasomes. This will get rid of chloroplastic and mitochondrial polymorphisms. 
awk -i inplace '$1 == ($1+0)' ./output/$1.ems.table

#add headers
<<<<<<< HEAD
sed -i '1s/^/'chr'\t'pos'\t'ref'\t'alt'\t'gene'\t'snpEffect'\t'snpVariant'\t'snpImpact'\t'mu:wt_GTpred'\t'wt_alt'\t'wt_ref'\t'mu_alt'\t'mu_ref'\n/' ./output/$1.ems.table
=======
sed -i '1s/^/'chr'\t'pos'\t'ref'\t'alt'\t'gene'\t'snpEffect'\t'snpVariant'\t'snpImpact'\t'mu:wt_GTpred'\t'wt_alt'\t'wt_ref'\t'mu_alt'\t'mu_ref'\n/' ./output/$1.ems.table


#clean up and organize. Change cleanup variable to "False", or comment out to disable.  
if [ cleanup=="True" ]; then
	rm ./output/*.tmp
	rm *.bam
	rm *.sam
	rm ${1}.table
	rm ${}1}.hc.vcf
fi

mv ./output/* ./output/archive/${1}
>>>>>>> ee0054a (add plot production to .py script, add archive management and cleanup to vcfgen.sh script)
