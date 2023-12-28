#!/usr/bin/env bash

# Organize variables passed from thale_bsa_utils.vcf_file_generation(args)
line_name="${1}"
pairedness="${2}"
allele_R_or_D="${3}"
reference_genome_name="${4}"
snpEff_db_name="${5}"
reference_genome_source="${6}"
threads_limit="${7}"
threads_halfed="$((threads_limit / 2))"
cleanup="${8}"
known_snps="${9}"
vcf_table_uuid="${10}"

# Initialize path variables
formatted_timestamp=$(date "+%Y.%m.%d ~%H:%M")
reference_dir="./references"
reference_genome_path="$reference_dir/$reference_genome_name.fa"
reference_chrs_path="$reference_dir/$reference_genome_name.chrs"
input_dir="./input"
input_name_prefix="${input_dir}/${line_name}.${allele_R_or_D}"
output_dir="./output/${line_name}"
output_file_prefix="$output_dir/${line_name}"
snpeff_dir="$output_dir/snpEff"

echo "Line Name: $line_name"
echo "Pairedness: $pairedness"
echo "Allele (R or D): $allele_R_or_D"
echo "Reference Genome: $reference_genome_name"
echo "SNPEff DB Name: $snpEff_db_name"
echo "Reference Genome Source: $reference_genome_source"
echo "Threads Limit: $threads_limit"
echo "Threads Halved: $threads_halfed"
echo "Cleanup: $cleanup"
echo "Known SNPs: $known_snps"
echo "VCF Table UUID: $vcf_table_uuid"
echo "Formatted Timestamp: $formatted_timestamp"
echo "Reference Genome Path: $reference_genome_path"
echo "Reference Chromosomes Path: $reference_chrs_path"
echo "Input Directory: $input_dir"
echo "Input Name Prefix: $input_name_prefix"
echo "Output Directory: $output_dir"
echo "Output File Prefix: $output_file_prefix"
echo "SNPEff Directory: $snpeff_dir"


echo "$formatted_timestamp Preparing references and directory structure"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# Create directory structure if it doesn't exist
mkdir -p "${output_dir}"
mkdir -p "${snpeff_dir}"
mkdir -p "${reference_dir}"

# Make reference genome if it doesn't exist
if ! [ -f "$reference_genome_path" ]; then
    curl -o "$reference_genome_path.gz" "$reference_genome_source" && \
        gzip -d "$reference_genome_path.gz"
fi

# Make .chrs file if it doesn't exist, set reference variable
if ! [ -f "$reference_chrs_path" ]; then
    awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' \
    "$reference_genome_path" > "$reference_chrs_path"
fi

# creating .fai and index files if they don't exist
if ! [ -f "${reference_chrs_path}.fa.fai" ]; then
    samtools faidx "${reference_chrs_path}"
    bwa index -p "${reference_chrs_path}.fa" -a is "${reference_genome_path}"
fi
${reference_chrs_path}.fa
echo "reference chrm path = ${reference_chrs_path}.fa.fai"
# create dictionary for gatk haplotype caller if it doesn't exist
if [ ! -f "${reference_chrs_path}.dict" ]; then
    picard CreateSequenceDictionary \
        -R "$reference_genome_path" \
        -O "${reference_chrs_path}.dict"
fi

echo "References prepared, preparing VCF for ${line_name}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp ${line_name} reads seem to be ${pairedness}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "Current Working Directory: $(pwd)"

# Set input files as paired-end or single-read
if [ "${pairedness}" == "paired-end" ]; then
    input_files_wt="${input_name_prefix}_1.wt.fq.gz \
        ${input_name_prefix}_2.wt.fq.gz"
    input_files_mu="${input_name_prefix}_1.mu.fq.gz \
        ${input_name_prefix}_2.mu.fq.gz"
fi

if [ "${pairedness}" == "single-read" ]; then
    input_files_wt="${input_name_prefix}.wt.fq.gz"
    input_files_mu="${input_name_prefix}.mu.fq.gz"
fi
echo "old fa : ${reference_chrs_path}.fa"
# Mapping
bwa mem \
    -t "$threads_halfed" \
    -M "${reference_chrs_path}.fa" \
    -v 1 $input_files_wt > "${output_file_prefix}_wt.sam" &
bwa mem \
    -t "$threads_halfed" \
    -M "${reference_chrs_path}.fa" \
    -v 1 $input_files_mu > "${output_file_prefix}_mu.sam"

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Converting sam to bam"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

samtools view \
    -bSh \
    -@ "$threads_halfed" \
    "${output_file_prefix}_mu.sam" > "${output_file_prefix}_mu.bam" &
samtools view \
    -bSh \
    -@ "$threads_halfed" \
    "${output_file_prefix}_wt.sam" > "${output_file_prefix}_wt.bam"
wait

# Fix paired-end
if [ "${pairedness}" == "paired-end" ]; then
    samtools fixmate "${output_file_prefix}_wt.bam" "${output_file_prefix}_wt.fix.bam" &
    samtools fixmate "${output_file_prefix}_mu.bam" "${output_file_prefix}_mu.fix.bam"
fi

echo "$formatted_timestamp Sorting by coordinate"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# Coordinate sorting
picard SortSam \
    I="${output_file_prefix}_mu.fix.bam" \
    O="${output_file_prefix}_mu.sort.bam" \
    SORT_ORDER=coordinate &
picard SortSam \
    I="${output_file_prefix}_wt.fix.bam" \
    O="${output_file_prefix}_wt.sort.bam" \
    SORT_ORDER=coordinate
wait

# Mark duplicates. Consider using sambamba for this step.
picard MarkDuplicates \
    I="${output_file_prefix}_mu.sort.bam" \
    O="${output_file_prefix}_mu.sort.md.bam" \
    METRICS_FILE="${output_file_prefix}_mu.metrics.txt" \
    ASSUME_SORTED=true &
picard MarkDuplicates \
    I="${output_file_prefix}_wt.sort.bam" \
    O="${output_file_prefix}_wt.sort.md.bam" \
    METRICS_FILE="${output_file_prefix}_wt.metrics.txt" \
    ASSUME_SORTED=true
wait

# Add header for GATK
picard AddOrReplaceReadGroups \
    I="${output_file_prefix}_mu.sort.md.bam" \
    O="${output_file_prefix}_mu.sort.md.rg.bam" \
    RGLB="${line_name}_mu" \
    RGPL=illumina RGSM="${line_name}_mu" \
    RGPU=run1 \
    SORT_ORDER=coordinate &
picard AddOrReplaceReadGroups \
    I="${output_file_prefix}_wt.sort.md.bam" \
    O="${output_file_prefix}_wt.sort.md.rg.bam" \
    RGLB="${line_name}_wt" \
    RGPL=illumina RGSM="${line_name}_wt" \
    RGPU=run1 SORT_ORDER=coordinate
wait

# Build BAM index
picard BuildBamIndex \
    INPUT="${output_file_prefix}_mu.sort.md.rg.bam" \
    O="${output_file_prefix}_mu.sort.md.rg.bai" &
picard BuildBamIndex \
    INPUT="${output_file_prefix}_wt.sort.md.rg.bam" \
    O="${output_file_prefix}_wt.sort.md.rg.bai"
wait

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Calling haplotypes. This may take awhile..."
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# GATK HC Variant calling
gatk HaplotypeCaller \
    -R "$reference_genome_path" \
    -I "${output_file_prefix}_mu.sort.md.rg.bam" \
    -I "${output_file_prefix}_wt.sort.md.rg.bam" \
    -O "${output_file_prefix}.hc.vcf" \
    -output-mode EMIT_ALL_CONFIDENT_SITES \
    --native-pair-hmm-threads_limit "$threads_limit"

# snpEff, labeling snps with annotations and potential impact on gene function
snpEff "$reference_genome_name" \
    -s "${output_dir}/snpEff/${1}_snpEff_summary.html" \
    "${output_file_prefix}.hc.vcf" > "${output_file_prefix}.se.vcf"

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Haplotypes called and snps labeled. Cleaning data."
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# Extract snpEFF data and variant information into a table, 
SnpSift extractFields \
    -s ":" \
    -e "NaN" "${output_file_prefix}.se.vcf" \
    CHROM POS REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].HGVS_P" \
    "ANN[*].IMPACT" "GEN[*].GT" "GEN[${line_name}_mu].AD" \
    "GEN[${line_name}_wt].AD" > "${output_file_prefix}.table"

# remove repetitive NaN's, retain EMS point mutations

sed -i 's/NaN://g' "${output_file_prefix}.ems.table.tmp"

grep \
    -e $'G\tA' \
    -e $'C\tT' \
    -e $'A\tG' \
    -e $'T\tC' \
    "${output_file_prefix}.table" > "${output_file_prefix}.ems.table.tmp"


# [mu:wt] genotypes. Grab appropriate genotypes for analysis. 
# 0/1:0/1 included in both analyses due to occasional leaky genotyping by GATK HC.
if [ "${allele_R_or_D}" = 'R' ]; then 
    grep \
        -F \
        -e '1/1:0/1' \
        -e '0/1:0/0' \
        -e '0/1:0/1' \
        "${output_file_prefix}.ems.table.tmp" > "${output_file_prefix}.ems.table"
else 
    grep \
        -F \
        -e '0/1:0/0' \
        -e '1/1:0/0' \
        -e '0/1:0/1' \
        "${output_file_prefix}.ems.table.tmp" > "${output_file_prefix}.ems.table"
fi

awk \
    -i inplace \
    -F'\t' \
    -vOFS='\t' \
    '{ gsub(",", "\t", $9) ; gsub(",", "\t", ${line_name}0) ; \
    gsub(",", "\t", ${line_name}1) ; print }' "${output_file_prefix}.ems.table"

# Remove complex genotypes
awk -i inplace -F'\t' 'NF==13' "${output_file_prefix}.ems.table"

# Get rid of chloroplastic and mitochondrial polymorphisms.
awk \
    -i inplace \
    '${line_name} == (${line_name}+0)' "${output_file_prefix}.ems.table"

noknownsnps_tablename="${output_file_prefix}.noknownsnps.table"
# Remove known snps
awk \
    'FNR==NR{a[${line_name}${pairedness}];next};!(${line_name}${pairedness} in a) \
    || ${line_name}~/#CHROM/' "$known_snps" "${output_file_prefix}.ems.table" \
    > "$noknownsnps_tablename"

# Add headers
sed \
    -i \
    '1s/^/'chr'\t'pos'\t'ref'\t'alt'\t'gene'\t'snpEffect'\t'snpVariant'\
    \t'snpImpact'\t'mu:wt_GTpred'\t'mu_ref'\t'mu_alt'\t'wt_ref'\t'wt_alt'\n/' \
    "$noknownsnps_tablename"

#add unique identifiers? 

# Clean up. Change cleanup variable to "False", or comment out to disable.
if [ "$cleanup" = "True" ]; then
    rm "$output_dir/${line_name}"/*.tmp
    rm "$output_dir/${line_name}"/*.bam
    rm "$output_dir/${line_name}"/*.sam
    rm "$output_dir/${line_name}"/*.idx
    rm "$output_dir/${line_name}"/*.bai
    rm "$output_dir/${line_name}"/*.matrics.txt
    rm "$output_dir/${line_name}/${1}.table"
    rm "$output_dir/${line_name}/${1}.ems.table"
    rm "$output_dir/${line_name}/${1}.hc.vcf"
fi
