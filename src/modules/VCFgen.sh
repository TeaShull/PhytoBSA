#!/usr/bin/env bash
# Organize variables passed from thale_bsa_utils.vcf_file_generation(args)
# args - 
# (key, value['reads'], value['allele'],
#     reference_genome_name, snpEff_db_name,
#     reference_genome_source, threads_limit,
#     cleanup, known_snps
# )

line_name = $1
pairedness = $2
allele_R_or_D = $3
reference_genome_name = $4
snpEff_db_name = $5
reference_genome_source = $6
threads_limit = $7
threads_halfed = ${threads_limit}/2
cleanup = $8
known_snps = $9

echo "Preparing references and directory structure"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# make directory structure if it doesn't exist.
if [ ! -d "./output/${line_name}" ]; then
    mkdir "./output/${line_name}"
fi

if [ ! -d "./output/VCFs" ]; then
    mkdir "./output/VCFs"
fi

if [ ! -d "./output/${line_name}/snpEff" ]; then
    mkdir "./output/${line_name}/snpEff"
fi

# make reference genome if it doesn't exist
if ! [ -f "./references/$reference_genome_name.fa" ]; then
    curl -o "./references/$reference_genome_name.fa.gz" "$reference_genome_source"
    gzip -d "./references/$reference_genome_name.fa.gz"
fi

# make .chrs file if it doesn't exist, set reference variable
if ! [ -f "./references/$reference_genome_name.chrs.fa" ]; then
    awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' "./references/$reference_genome_name.fa" > "./references/$reference_genome_name.chrs.fa"
fi

reference_genome="./references/$reference_genome_name.chrs.fa"

# creating .fai and index files if they don't exist
if ! [ -f "./references/$reference_genome_name.chrs.fa.fai" ]; then
    samtools faidx "$reference_genome"
    bwa index -p "./references/$reference_genome_name.chrs.fa" -a is "$reference_genome"
fi

# create dictionary for gatk haplotype caller if it doesn't exist
if [ ! -f "./references/$reference_genome_name.chrs.dict" ]; then
    picard CreateSequenceDictionary \
        -R "$reference_genome" \
        -O "./references/$reference_genome_name.chrs.dict"
fi

echo "References prepared, preparing VCF for ${line_name}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "${line_name} reads seem to be ${pairedness}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# set input files as paired-end or single-read
if [ "${pairedness}" == "paired-end" ]; then
    input_files_wt="./input/${line_name}.${allele_R_or_D}_1.wt.fq.gz ./input/${line_name}.${allele_R_or_D}_2.wt.fq.gz"
    input_files_mu="./input/${line_name}.${allele_R_or_D}_1.mu.fq.gz ./input/${line_name}.${allele_R_or_D}_2.mu.fq.gz"
fi

if [ "${pairedness}" == "single-read" ]; then
    input_files_wt="./input/${line_name}.${allele_R_or_D}.wt.fq.gz"
    input_files_mu="./input/${line_name}.${allele_R_or_D}.mu.fq.gz"
fi

# mapping
bwa mem \
    -t "$threads_halfed" \
    -M "$reference_genome" \
    -v 1 \
    $input_files_wt > "./output/${line_name}/${line_name}_wt.sam" &
bwa mem \
    -t "$threads_halfed" \
    -M "$reference_genome" \
    -v 1 \
    $input_files_mu > "./output/${line_name}/${line_name}_mu.sam"

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "Converting sam to bam"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

samtools view \
    -bSh \
    -@ "$threads_halfed" \
    "./output/${line_name}/${line_name}_mu.sam" > "./output/${line_name}/${line_name}_mu.bam" &
samtools view \
    -bSh \
    -@ "$threads_halfed" \
    "./output/${line_name}/${line_name}_wt.sam" > "./output/${line_name}/${line_name}_wt.bam"
wait

# fix paired end
if [ "${pairedness}" == "paired-end" ]; then
    samtools fixmate \
        "./output/${line_name}/${line_name}_mu.bam" "./output/${line_name}/${line_name}_mu.fix.bam" &
    samtools fixmate \
        "./output/${line_name}/${line_name}_wt.bam" "./output/${line_name}/${line_name}_wt.fix.bam"
fi

echo "Sorting by coordinate"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# Coordinate sorting
picard SortSam \
    I="./output/${line_name}/${line_name}_mu.fix.bam" \
    O="./output/${line_name}/${line_name}_mu.sort.bam" \
    SORT_ORDER=coordinate &
picard SortSam \
    I="./output/${line_name}/${line_name}_wt.fix.bam" \
    O="./output/${line_name}/${line_name}_wt.sort.bam" \
    SORT_ORDER=coordinate
wait

# Mark duplicates. Consider using sambamba for this step.
picard MarkDuplicates \
    I="./output/${line_name}/${line_name}_mu.sort.bam" \
    O="./output/${line_name}/${line_name}_mu.sort.md.bam" \
    METRICS_FILE="./output/${line_name}/${line_name}_mu.metrics.txt" \
    ASSUME_SORTED=true &
picard MarkDuplicates \
    I="./output/${line_name}/${line_name}_wt.sort.bam" \
    O="./output/${line_name}/${line_name}_wt.sort.md.bam" \
    METRICS_FILE="./output/${line_name}/${line_name}_wt.metrics.txt" \
    ASSUME_SORTED=true
wait

# add header for gatk
picard AddOrReplaceReadGroups \
    I="./output/${line_name}/${line_name}_mu.sort.md.bam" \
    O="./output/${line_name}/${line_name}_mu.sort.md.rg.bam" \
    RGLB="${line_name}_mu" \
    RGPL=illumina \
    RGSM="${line_name}_mu" \
    RGPU=run1 \
    SORT_ORDER=coordinate &
picard AddOrReplaceReadGroups \
    I="./output/${line_name}/${line_name}_wt.sort.md.bam" \
    O="./output/${line_name}/${line_name}_wt.sort.md.rg.bam" \
    RGLB="${line_name}_wt" \
    RGPL=illumina \
    RGSM="${line_name}_wt" \
    RGPU=run1 \
    SORT_ORDER=coordinate
wait

# build BAM index
picard BuildBamIndex \
    INPUT="./output/${line_name}/${line_name}_mu.sort.md.rg.bam" \
    O="./output/${line_name}/${1}_mu.sort.md.rg.bai" &
picard BuildBamIndex \
    INPUT="./output/${line_name}/${line_name}_wt.sort.md.rg.bam" \
    O="./output/${line_name}/${1}_wt.sort.md.rg.bai"
wait

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "Calling haplotypes. This may take awhile..."
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# GATK HC Variant calling
gatk HaplotypeCaller \
    -R "$reference_genome" \
    -I "./output/${line_name}/${line_name}_mu.sort.md.rg.bam" \
    -I "./output/${line_name}/${line_name}_wt.sort.md.rg.bam" \
    -O "./output/${line_name}/${line_name}.hc.vcf" \
    -output-mode EMIT_ALL_CONFIDENT_SITES \
    --native-pair-hmm-threads_limit $threads_limit

# snpEff, labeling snps with annotations and potential impact on gene function
snpEff "$reference_genome_name" \
    -s "./output/${line_name}/snpEff/${1}_snpEff_summary.html" \
    "./output/${line_name}/${line_name}.hc.vcf" > "./output/${line_name}/${line_name}.se.vcf"

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "Haplotypes called and snps labeled. Cleaning data."
echo	">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

#extract snpEFF data and variant information into a table, remove repetative NaN's and retain only those polymorphisms likely to arise from EMS.
SnpSift extractFields \
	-s ":" \
	-e "NaN" \
	./output/${line_name}/${line_name}.se.vcf \
	CHROM POS  REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].IMPACT" "GEN[*].GT" "GEN[${line_name}_mu].AD" "GEN[${line_name}_wt].AD" > ./output/${line_name}/${line_name}.table

grep -e $'G\tA' -e $'C\tT' -e $'A\tG' -e $'T\tC' ./output/${line_name}/${line_name}.table > ./output/${line_name}/${line_name}.ems.table.tmp

sed -i 's/NaN://g' ./output/${line_name}/${line_name}.ems.table.tmp

#[mu:wt] genotypes. Grab appropriate genotypes for analysis. 0/1:0/1 included in both analyses due to occasianal leaky genotyping by GATK HC. 
if [ "${allele_R_or_D}" = 'R' ]; then 
	grep -F -e '1/1:0/1' -e '0/1:0/0' -e '0/1:0/1' ./output/${line_name}/${line_name}.ems.table.tmp > ./output/${line_name}/${line_name}.ems.table
else 
	grep -F -e '0/1:0/0' -e '1/1:0/0' -e '0/1:0/1' ./output/${line_name}/${line_name}.ems.table.tmp > ./output/${line_name}/${line_name}.ems.table
fi

awk -i inplace -F'\t' -vOFS='\t' '{ gsub(",", "\t", $9) ; gsub(",", "\t", ${line_name}0) ; gsub(",", "\t", ${line_name}1) ; print }' ./output/${line_name}/${line_name}.ems.table

#remove complex genotypes
awk -i inplace -F'\t' 'NF==13' ./output/${line_name}/${line_name}.ems.table

#remove non-numeric chromasomes. This will get rid of chloroplastic and mitochondrial polymorphisms. 
awk -i inplace '${line_name} == (${line_name}+0)' ./output/${line_name}/${line_name}.ems.table

#remove known snps
awk 'FNR==NR{a[${line_name}${pairedness}];next};!(${line_name}${pairedness} in a) || ${line_name}~/#CHROM/' $knownsnps ./output/${line_name}/${line_name}.ems.table > ./output/${line_name}/${line_name}.noknownsnps.table

#add headers
sed -i '1s/^/'chr'\t'pos'\t'ref'\t'alt'\t'gene'\t'snpEffect'\t'snpVariant'\t'snpImpact'\t'mu:wt_GTpred'\t'mu_ref'\t'mu_alt'\t'wt_ref'\t'wt_alt'\n/' ./output/${line_name}/${line_name}.noknownsnps.table

#clean up and organize. Change cleanup variable to "False", or comment out to disable.  
if [ "$cleanup" = True ]; then
	rm ./output/${line_name}/*.tmp
	rm ./output/${line_name}/*.bam
	rm ./output/${line_name}/*.sam
	rm ./output/${line_name}/*.idx
	rm ./output/${line_name}/*.bai
	rm ./output/${line_name}/*.matrics.txt
	rm ./output/${line_name}/${1}.table
	rm ./output/${line_name}/${1}.ems.table
	rm ./output/${line_name}/${1}.hc.vcf

fi

