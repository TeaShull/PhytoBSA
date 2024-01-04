#!/usr/bin/env bash
source VCFfunctions.sh

# Organize variables passed from thale_bsa_utils.vcf_file_generation(args)
vcf_ulid="${1}"
current_line_name="${2}"
allele="${3}"
input_dir="${4}"
wt_input="${5}"
mu_input="${6}"
output_dir_path="${7}"
output_prefix="${8}"
vcf_table_path="${9}"
reference_dir="${10}"
reference_genome_name="${11}"
snpEff_species_db="${12}"
reference_genome_source="${13}"
known_snps="${14}"
threads_limit="${15}"
cleanup="${16}"

# Generate other needed variables
formatted_timestamp=$(date "+%Y.%m.%d ~%H:%M")
reference_genome_path="$reference_dir/$reference_genome_name.fa"
reference_chrs_path="$reference_dir/$reference_genome_name.chrs"
reference_chrs_fa_path="$reference_dir/$reference_genome_name.chrs.fa"
snpeff_dir="$output_dir_path/snpEff"

#Final file output names
noknownsnps_tablename="${output_prefix}.noknownsnps.table"
ems_file_name="${output_prefix}.ems.table"

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Variables initiated."
echo "Passed to subprocess:"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "vcf_ulid: ${vcf_ulid}"
echo "current_line_name: ${current_line_name}"
echo "allele: ${allele}"
echo "input_dir: ${input_dir}"
echo "wt_input: ${wt_input}"
echo "mu_input: ${mu_input}"
echo "output_dir_path: ${output_dir_path}"
echo "output_prefix: ${output_prefix}"
echo "vcf_table_path: ${vcf_table_path}"
echo "reference_dir: ${reference_dir}"
echo "reference_genome_name: ${reference_genome_name}"
echo "snpEff_species_db: ${snpEff_species_db}"
echo "reference_genome_source: ${reference_genome_source}"
echo "known_snps: ${known_snps}"
echo "threads_limit: ${threads_limit}"
echo "cleanup: ${cleanup}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "Generated variables:"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "formatted_timestamp: ${formatted_timestamp}"
echo "reference_genome_path: ${reference_genome_path}"
echo "reference_chrs_path: ${reference_chrs_path}"
echo "reference_chrs_fa_path: ${reference_chrs_fa_path}"
echo "snpeff_dir: ${snpeff_dir}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "VCFgen output file names:"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "noknownsnps_tablename: ${noknownsnps_tablename}"
echo "ems_file_name: ${ems_file_name}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Preparing references and directory structure"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

#Create directory structure if it doesn't exist
mkdir -p "${output_dir}"
mkdir -p "${snpeff_dir}"
mkdir -p "${reference_dir}"

# Make reference genome if it doesn't exist
if ! [ -f "$reference_genome_path" ]; then
    curl -o "$reference_genome_path" "$reference_genome_source" && \
        gzip -d "$reference_genome_path.gz"
fi

# Make .chrs file if it doesn't exist, set reference variable
if ! [ -f "$reference_chrs_fa_path" ]; then
    awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' \
    $reference_genome_path > $reference_chrs_fa_path
fi

# creating .fai and index files if they don't exist
if ! [ -f "${reference_chrs_fa_path}.fai" ]; then
    samtools faidx "${reference_chrs_fa_path}"
    bwa index -p "${reference_chrs_fa_path}" -a is "${reference_genome_path}"
fi


# create dictionary for gatk haplotype caller if it doesn't exist
if [ ! -f "${reference_chrs_path}.dict" ]; then
    picard CreateSequenceDictionary \
        -R "${reference_chrs_fa_path}" \
        -O "${reference_chrs_path}.dict"
fi

echo "References prepared, preparing VCF for ${line_name}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp ${line_name} reads seem to be ${pairedness}"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "Current Working Directory: $(pwd)"
# Set input files as paired-end or single-read

echo "wt_input: $wt_input"
echo "mu_input: $mu_input"

# Mapping
bwa mem \
    -t $threads_halfed \
    -M "${reference_chrs_fa_path}" \
    $wt_input > "${output_prefix}_wt.sam" &
bwa mem \
    -t $threads_halfed \
    -M "${reference_chrs_fa_path}" \
    $mu_input > "${output_prefix}_mu.sam"
wait

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Converting sam to bam"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

samtools view \
    -bSh \
    -@ "$threads_halfed" \
    "${output_prefix}_mu.sam" > "${output_prefix}_mu.bam" &
samtools view \
    -bSh \
    -@ "$threads_halfed" \
    "${output_prefix}_wt.sam" > "${output_prefix}_wt.bam"
wait

# Fix paired-end
if [ "${pairedness}" == "paired-end" ]; then
    samtools fixmate "${output_prefix}_wt.bam" \
        "${output_prefix}_wt.fix.bam" &
    samtools fixmate "${output_prefix}_mu.bam" \
        "${output_prefix}_mu.fix.bam"
fi

echo "$formatted_timestamp Sorting by coordinate"
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# Coordinate sorting
picard SortSam \
    I="${output_prefix}_mu.fix.bam" \
    O="${output_prefix}_mu.sort.bam" \
    SORT_ORDER=coordinate &

picard SortSam \
    I="${output_prefix}_wt.fix.bam" \
    O="${output_prefix}_wt.sort.bam" \
    SORT_ORDER=coordinate

wait

# Mark duplicates. Consider using sambamba for this step.
picard MarkDuplicates \
    I="${output_prefix}_mu.sort.bam" \
    O="${output_prefix}_mu.sort.md.bam" \
    METRICS_FILE="${output_prefix}_mu.metrics.txt" \
    ASSUME_SORTED=true &
picard MarkDuplicates \
    I="${output_prefix}_wt.sort.bam" \
    O="${output_prefix}_wt.sort.md.bam" \
    METRICS_FILE="${output_prefix}_wt.metrics.txt" \
    ASSUME_SORTED=true
wait

# Add header for GATK
picard AddOrReplaceReadGroups \
    I="${output_prefix}_mu.sort.md.bam" \
    O="${output_prefix}_mu.sort.md.rg.bam" \
    RGLB="${line_name}_mu" \
    RGPL=illumina \
    RGSM="${line_name}_mu" \
    RGPU=run1 \
    SORT_ORDER=coordinate &
picard AddOrReplaceReadGroups \
    I="${output_prefix}_wt.sort.md.bam" \
    O="${output_prefix}_wt.sort.md.rg.bam" \
    RGLB="${line_name}_wt" \
    RGPL=illumina \
    RGSM="${line_name}_wt" \
    RGPU=run1 \
    SORT_ORDER=coordinate
wait

# Build BAM index
picard BuildBamIndex \
    INPUT="${output_prefix}_mu.sort.md.rg.bam" \
    O="${output_prefix}_mu.sort.md.rg.bai" &
picard BuildBamIndex \
    INPUT="${output_prefix}_wt.sort.md.rg.bam" \
    O="${output_prefix}_wt.sort.md.rg.bai"
wait

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Calling haplotypes. This may take awhile..."
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# GATK HC Variant calling
gatk HaplotypeCaller \
    -R "$reference_chrs_fa_path" \
    -I "${output_prefix}_mu.sort.md.rg.bam" \
    -I "${output_prefix}_wt.sort.md.rg.bam" \
    -O "${output_prefix}.hc.vcf" \
    -output-mode EMIT_ALL_CONFIDENT_SITES \
    --native-pair-hmm-threads "$threads_limit"

# snpEff, labeling snps with annotations and potential impact on gene function
snpEff "$snpEff_db_name" \
    -s "${output_dir}/snpEff/${1}_snpEff_summary.html" \
    "${output_prefix}.hc.vcf" > "${output_prefix}.se.vcf"

echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"
echo "$formatted_timestamp Haplotypes called and snps labeled. Cleaning data."
echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=<"

# Extract snpEFF data and variant information into a table, 
SnpSift extractFields \
    -s ":" \
    -e "NaN" "${output_prefix}.se.vcf" \
    CHROM POS REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].HGVS_P" \
    "ANN[*].IMPACT" "GEN[*].GT" "GEN[${line_name}_mu].AD" \
    "GEN[${line_name}_wt].AD" > "${output_prefix}.table.tmp"

# remove repetitive NaN's, retain EMS point mutations

sed -i 's/NaN://g' "${output_prefix}.table.tmp"

grep \
    -e $'G\tA' \
    -e $'C\tT' \
    -e $'A\tG' \
    -e $'T\tC' \
    "${output_prefix}.table.tmp" > "${output_prefix}.ems.table.tmp"

# [mu:wt] genotypes. Grab appropriate genotypes for analysis. 
# 0/1:0/1 included in both analyses due to occasional leaky genotyping by GATK HC.
if [ "${allele}" = 'R' ]; then 
    grep \
        -F \
        -e '1/1:0/1' \
        -e '0/1:0/0' \
        -e '0/1:0/1' \
        "${output_prefix}.ems.table.tmp" > "${ems_file_name}"
else 
    grep \
        -F \
        -e '0/1:0/0' \
        -e '1/1:0/0' \
        -e '0/1:0/1' \
        "${output_prefix}.ems.table.tmp" > "${ems_file_name}"
fi

awk \
    -i inplace \
    -F'\t' \
    -vOFS='\t' \
    '{ gsub(",", "\t", $9) ; \
    gsub(",", "\t", '${line_name}'0) ; \
    gsub(",", "\t", '${line_name}'1) ; print }' \
    "${ems_file_name}"


# Remove complex genotypes
awk \
    -i inplace \
    -F'\t' 'NF==13' \
    "${ems_file_name}"

# Get rid of chloroplastic and mitochondrial polymorphisms.
awk -i inplace \
    '${line_name} == (${line_name}+0)' \
    "${ems_file_name}"

# Remove known snps
awk \
    'FNR==NR{a[${line_name}${pairedness}];next};!(${line_name}${pairedness} in a) || ${line_name}~/#CHROM/' "$known_snps" "${output_prefix}.ems.table" \
    > "$noknownsnps_tablename"

# Add headers
sed \
    -i \
    '1s/^/'chr'\t'pos'\t'ref'\t'alt'\t'gene'\t'snpEffect'\t'snpVariant'\t'snpImpact'\t'mu:wt_GTpred'\t'mu_ref'\t'mu_alt'\t'wt_ref'\t'wt_alt'\n/' "$noknownsnps_tablename"


# Clean up. Change cleanup variable to "False", or comment out to disable.
if [ "$cleanup" = "True" ]; then
    rm "${output_dir_path}"/*.tmp
    rm "${output_dir_path}"/*.bam
    rm "${output_dir_path}"/*.sam
    rm "${output_dir_path}"/*.idx
    rm "${output_dir_path}"/*.bai
    rm "${output_dir_path}"/*.matrics.txt
    rm "${output_dir_path}.hc.vcf"
fi
