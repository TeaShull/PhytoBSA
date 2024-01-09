#!/usr/bin/env bash

# print a separator for better readability
print_separator() {
    echo ">=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<"
}

# Prints the variable info in an associative array of variables with a header.
# $1, Associative array of variables, $2 - "message"
print_variable_info() {
    # Ensure arguments passed is at least 2
    if [ "$#" -lt 2 ]; then
        echo "Error: Insufficient arguments for print_variable_info"
        return 1
    fi

    #assign message content and associative array name. 
    local -n array=$1
    local custom_message="$2"
    shift 2

    #Print variables
    print_separator
    echo "$formatted_timestamp $custom_message."
    print_separator

    for var in "${!array[@]}"; do
        echo "${var}: ${array[$var]}"
    done
}

assign_values() {
    # Check if the associative array is provided
    local -n array=$1
    if [ "$#" -eq 0 ]; then
        echo "Error: No associative array provided."
        return 1
    fi

    # Loop through the associative array
    for key in "${!array[@]}"; do
        # Assign the value to a variable with the same name as the key
        declare -g "${key}=${array[${key}]}"
    done
}


# Print a timestamp and message
print_message() {
    echo "$formatted_timestamp $1"
    print_separator
}

# Create directories if they don't exist
create_directories() {
    # Ensure at least one directory is provided
    if [ "$#" -eq 0 ]; then
        echo "Error: No directories provided for create_directories"
        return 1
    fi

    for dir in "${@}"; do
        mkdir -p "${dir}"
    done
}

# Download reference genome if it doesn't exist
download_reference_genome() {
    local reference_genome_path="$1"
    local reference_genome_source="$2"

    if [ ! -f "${reference_genome_path}" ]; then
        curl -o "${reference_genome_path}" "${reference_genome_source}" && gzip -d "${reference_genome_path}.gz"
    fi
}

# Create .chrs file if it doesn't exist
create_chrs_file() {
    local reference_genome_path="$1"
    local reference_chrs_fa_path="$2"

    if [ ! -f "${reference_chrs_fa_path}" ]; then
        awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' "${reference_genome_path}" > "${reference_chrs_fa_path}"
    fi
}

# Create .fai and index files if they don't exist
create_fai_and_index() {
    local reference_genome_path="$1"
    local reference_chrs_fa_path="$2"

    if [ ! -f "${reference_chrs_fa_path}.fai" ]; then
        samtools faidx "${reference_chrs_fa_path}"
        bwa index -p "${reference_chrs_fa_path}" -a is "${reference_genome_path}"
    fi
}

# Create a dictionary for GATK haplotype caller if it doesn't exist
create_sequence_dictionary() {
    local reference_chrs_fa_path="$1"
    local reference_chrs_path="$2"

    if [ ! -f "${reference_chrs_path}.dict" ]; then
        picard CreateSequenceDictionary -R "${reference_chrs_fa_path}" -O "${reference_chrs_path}.dict"
    fi
}

# remove known SNPs
remove_known_snps() {
    local known_snps_file="$1"
    local input_file="$2"
    local output_file="$3"
    
    awk -v line_name="${line_name}" -v pairedness="${pairedness}" \
        'FNR==NR{a[$1$2];next};!(($1$2) in a) || $1~/#CHROM/' \
        "$known_snps_file" "$input_file" > "$output_file"
}

#remove polymorphisms in mitochondria and chloroplast
remove_nongenomic_polymorphisms() {
    local line_name="$1"
    local file="$2"
    awk -i inplace '$1 == ($1+0)' "$file"
}


# Add neccissary headers to file. 
add_headers() {
    local file="$1"
    sed -i '1s/^/chr\tpos\tref\talt\tgene\tsnpEffect\tsnpVariant\tsnpImpact\tmu:wt_GTpred\tmu_ref\tmu_alt\twt_ref\twt_alt\n/' "$file"
}

extract_fields_snpSift() {
    local input_file="$1"
    local output_file="$2"
    local line_name="$3"
    local extract_fields="CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].HGVS_P ANN[*].IMPACT GEN[*].GT GEN[${line_name}_mu].AD GEN[${line_name}_wt].AD"
    SnpSift extractFields -s ":" -e "NaN" "${input_file}" $extract_fields > "${output_file}"
}

remove_repetitive_nan() {
    local file="$1"
    sed -i 's/NaN://g' "$file"
}

# Grab mutations likely to be from ems for analysis
filter_ems_mutations() {
    local input_file="$1"
    local output_file="$2"
    grep -e $'G\tA' -e $'C\tT' -e $'A\tG' -e $'T\tC' "$input_file" > "$output_file"
}

# Filter [mu:wt] genotypes based on the allele
filter_genotypes() {
    local allele="$1"
    local input_file="$2"
    local output_file="$3"
    
    if [ "$allele" = 'R' ]; then
        grep -F -e '1/1:0/1' -e '0/1:0/0' -e '0/1:0/1' "$input_file" > "$output_file"
    else
        grep -F -e '0/1:0/0' -e '1/1:0/0' -e '0/1:0/1' "$input_file" > "$output_file"
    fi
}

# Format fields in the EMS file
format_ems_file() {
    local file="$1"
    local line_name="$2"
    awk -i inplace -F'\t' -vOFS='\t' '{ gsub(",", "\t", $9) ; gsub(",", "\t", $10) ; gsub(",", "\t", $11) ; print }' "$file"
}

# Remove complex genotypes from the EMS file
remove_complex_genotypes() {
    local file="$1"
    awk -i inplace -F'\t' 'NF==13' "$file"
}

# Clean up temporary files
cleanup_files() {
    local output_dir_path="$1"
    local cleanup=$2
    if [ "${cleanup}" == True ]; then
        rm "${output_dir_path}"/*.tmp
        rm "${output_dir_path}"/*.bam
        rm "${output_dir_path}"/*.sam
        rm "${output_dir_path}"/*.idx
        rm "${output_dir_path}"/*.bai
        rm "${output_dir_path}"/*.matrics.txt
        rm "${output_prefix}.hc.vcf"
    fi
}
