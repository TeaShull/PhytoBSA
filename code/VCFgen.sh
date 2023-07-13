#!/usr/bin/env bash
#pull in variables
source ./code/variables.sh

echo	"Preparing references and directory structure"
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

#make directory structure if it doesn't exist. 
if [ ! -d "./output/$1" ]; then
	mkdir ./output/$1
fi

if [ ! -d "./output/VCFs" ]; then
	mkdir ./output/VCFs
fi

if [ ! -d "./output/$1/snpEff" ]; then
	mkdir ./output/$1/snpEff
fi

# #get genomic features file if it doesn't exist (for later...)
# if ! [ -f ./references/$my_species.gff ]; then
#   curl -o ./references/$my_species.gff.gz $gff_file_source
#   gzip -d ./references/$my_species.gff.gz
# fi

#make reference genome if it doesn't exist 
if ! [ -f ./references/$my_species.fa ]; then
  curl -o ./references/$my_species.fa.gz $reference_genome_source
  gzip -d ./references/$my_species.fa.gz
fi

#make .chrs file if it doesn't exist, set reference variable
if ! [ -f ./references/$my_species.chrs.fa ]; then
	awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' ./references/$my_species.fa > ./references/$my_species.chrs.fa
fi

fa=./references/$my_species.chrs.fa

# #creating .fai and index files if they don't exist
if ! [ -f ./references/$my_species.chrs.fa.fai ]; then
	samtools faidx $fa
	bwa index -p ./references/$my_species.chrs.fa -a is $fa
fi

#create dictory for gatk haplotype caller if it doesn't exist
if [ ! -f "./references/$my_species.chrs.dict" ]; then
	picard CreateSequenceDictionary -R .fa -O ./references/$my_species.chrs.dict
fi

echo	"References prepared, preparing VCF for $1"
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo  "$1 reads seem to be $2"
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

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
bwa mem -t $threads -M $fa $fa_in_wt > ./output/$1/$1_wt.sam &
bwa mem -t $threads -M $fa $fa_in_mu > ./output/$1/$1_mu.sam

echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"Converting sam to bam"
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

samtools view -bSh -@ $threads ./output/$1/$1_mu.sam > ./output/$1/$1_mu.bam &
samtools view -bSh -@ $threads ./output/$1/$1_wt.sam > ./output/$1/$1_wt.bam
wait

#fix paired end
if [ "$2" == "paired-end" ]; then
	samtools fixmate ./output/$1/$1_mu.bam ./output/$1/$1_mu.fix.bam &
	samtools fixmate ./output/$1/$1_wt.bam ./output/$1/$1_wt.fix.bam
fi

echo	"Sorting by coordinate"
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

#Coordinate sorting
picard SortSam I=./output/$1/$1_mu.fix.bam O=./output/$1/$1_mu.sort.bam SORT_ORDER=coordinate &
picard SortSam I=./output/$1/$1_wt.fix.bam O=./output/$1/$1_wt.sort.bam  SORT_ORDER=coordinate
wait

#consider using sambamba for this step. I hear it is much faster, and this seems to take longer than it should
picard MarkDuplicates I=./output/$1/$1_mu.sort.bam O=./output/$1/$1_mu.sort.md.bam METRICS_FILE=./output/$1/$1_mu.matrics.txt ASSUME_SORTED=true &
picard MarkDuplicates I=./output/$1/$1_wt.sort.bam O=./output/$1/$1_wt.sort.md.bam METRICS_FILE=./output/$1/$1_wt.matrics.txt ASSUME_SORTED=true
wait

#add header for gatk
picard AddOrReplaceReadGroups I=./output/$1/$1_mu.sort.md.bam O=./output/$1/$1_mu.sort.md.rg.bam RGLB=$1_mu RGPL=illumina RGSM=$1_mu RGPU=run1 SORT_ORDER=coordinate &
picard AddOrReplaceReadGroups I=./output/$1/$1_wt.sort.md.bam O=./output/$1/$1_wt.sort.md.rg.bam RGLB=$1_wt RGPL=illumina RGSM=$1_wt RGPU=run1 SORT_ORDER=coordinate
wait

picard BuildBamIndex INPUT=./output/$1/$1_mu.sort.md.rg.bam O=./output/$1/${1}_mu.sort.md.rg.bai &
picard BuildBamIndex INPUT=./output/$1/$1_wt.sort.md.rg.bam O=./output/$1/${1}_wt.sort.md.rg.bai
wait

echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"Calling haplotypes. This may take awhile..."
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# #GATK HC Variant calling (I believe there is a way to parrelalize this by chr, and merge at the end. Improve later)
gatk HaplotypeCaller -R $fa -I ./output/$1/$1_mu.sort.md.rg.bam -I ./output/$1/$1_wt.sort.md.rg.bam -O ./output/$1/$1.hc.vcf -output-mode EMIT_ALL_CONFIDENT_SITES --native-pair-hmm-threads 20

snpEff, labeling snps with annotations and potential impact on gene function

snpEff $my_species -s ./output/$1/snpEff/${1}_snpEff_summary.html ./output/$1/$1.hc.vcf > ./output/$1/$1.se.vcf

echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo	"Haplotypes called and snps labeled. Cleaning data...."
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

#extract snpEFF data and variant information into a table, remove repetative NaN's and retain only those polymorphisms likely to arise from EMS.
SnpSift extractFields -s ":" -e "NaN" ./output/$1/$1.se.vcf CHROM POS  REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].IMPACT" "GEN[*].GT" "GEN[$1_mu].AD" "GEN[$1_wt].AD" > ./output/$1/$1.table

grep -e $'G\tA' -e $'C\tT' -e $'A\tG' -e $'T\tC' ./output/$1/$1.table > ./output/$1/$1.ems.table.tmp

sed -i 's/NaN://g' ./output/$1/$1.ems.table.tmp

#only tested on recessive mutations so far - need to get some verified dominant datasets from SRA to test.... This needs some attention, but works for now
if [ "$mutation" = 'recessive' ]; then 
	grep -F -e '1/1:0/1' -e '0/1:0/0' -e '1/1:0/0' ./output/$1/$1.ems.table.tmp > ./output/$1/$1.ems.table
else 
	grep -F -e $'0/1:0/0' -e '1/1:0/0' ./output/$1/$1.ems.table.tmp > ./output/$1/$1.ems.table
fi

awk -i inplace -F'\t' -vOFS='\t' '{ gsub(",", "\t", $9) ; gsub(",", "\t", $10) ; gsub(",", "\t", $11) ; print }' ./output/$1/$1.ems.table

#remove complex genotypes
awk -i inplace -F'\t' 'NF==13' ./output/$1/$1.ems.table

#remove non-numeric chromasomes. This will get rid of chloroplastic and mitochondrial polymorphisms. 
awk -i inplace '$1 == ($1+0)' ./output/$1/$1.ems.table

#remove known snps
awk 'FNR==NR{a[$1$2];next};!($1$2 in a) || $1~/#CHROM/' $knownsnps ./output/$1/$1.ems.table > ./output/$1/$1.noknownsnps.table

#add headers
sed -i '1s/^/'chr'\t'pos'\t'ref'\t'alt'\t'gene'\t'snpEffect'\t'snpVariant'\t'snpImpact'\t'mu:wt_GTpred'\t'mu_ref'\t'mu_alt'\t'wt_ref'\t'wt_alt'\n/' ./output/$1/$1.noknownsnps.table

#clean up and organize. Change cleanup variable to "False", or comment out to disable.  
if [ $cleanup ]; then
	rm ./output/$1/*.tmp
	rm ./output/$1/*.bam
	rm ./output/$1/*.sam
	rm ./output/$1/*.idx
	rm ./output/$1/*.bai
	rm ./output/$1/*.matrics.txt
	rm ./output/$1/${1}.table
	rm ./output/$1/${1}.ems.table
	rm ./output/$1/${1}.hc.vcf

fi

