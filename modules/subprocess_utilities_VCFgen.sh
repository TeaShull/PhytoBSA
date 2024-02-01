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
        wait
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
    print_separator
    echo "$formatted_timestamp $1"
    print_separator
}

# Create .fai and index files if they don't exist
create_fai_and_index() {
    local reference_chrs_fa_path="$2"

    if [ ! -f "${reference_chrs_fa_path}.fai" ]; then
        samtools faidx "${reference_chrs_fa_path}"
        bwa index -p "${reference_chrs_fa_path}" -a is "${reference_chrs_fa_path}"
    fi
}

# Create a dictionary for GATK haplotype caller if it doesn't exist
create_sequence_dictionary() {
    local reference_chrs_fa_path="$1"
    local reference_chrs_dict_path="$2"

    if [ ! -f "${reference_chrs_dict_path}" ]; then
        echo "Dictory doesn't exist yet - creating ${reference_chrs_dict_path}"
        picard CreateSequenceDictionary -R "${reference_chrs_fa_path}" -O "${reference_chrs_dict_path}"
    fi
}

split_and_call_haplotypes() {
    local output_prefix=$1
    local reference_chrs_fa_path=$2
    local picard_addorreplacereadgroups_output_wt=$3
    local picard_addorreplacereadgroups_output_mu=$4
    local threads_limit=$5
    local gatk_haplotypecaller_output=$6

    # Create an array of chromosome names (assuming the reference fasta file is indexed)
    mapfile -t chrs < <(awk '{print $1}' "${reference_chrs_fa_path}.fai")
    # Calculate threads per chromosome
    local threads_per_chr=$((threads_limit / ${#chrs[@]}))

    # Create an array to hold output file names
    declare -a output_files

    # Loop over chromosomes
    for chr in "${chrs[@]}"; do
        # Generate output file name
        local output_file="${output_prefix}_${chr}.vcf"
        output_files+=("$output_file")

        # Call HaplotypeCaller for each chromosome
        gatk HaplotypeCaller \
            -R "$reference_chrs_fa_path" \
            -I "$picard_addorreplacereadgroups_output_mu" \
            -I "$picard_addorreplacereadgroups_output_wt" \
            -O "$output_file" \
            -L "$chr" \
            -output-mode EMIT_ALL_CONFIDENT_SITES \
            --native-pair-hmm-threads "$threads_per_chr" &
    done

    # Wait for all background jobs to finish
    wait
    echo "Haplotypes called chromosome by chromosomes. Merging VCF files..."

    # Merge VCF files
    gatk MergeVcfs -I "${output_files[@]}" -O "$gatk_haplotypecaller_output"

    # Remove individual chromosome VCF files
    for output_file in "${output_files[@]}"; do
        rm -f "$output_file"
    done
}

extract_fields_snpSift() {
    local input_file="$1"
    local output_file="$2"
    local line_name="$3"
    local extract_fields="CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].HGVS_P ANN[*].IMPACT GEN[*].GT GEN[${line_name}_mu].AD GEN[${line_name}_wt].AD"
    SnpSift extractFields -s ":" -e "NaN" "$input_file" "$extract_fields" >"$output_file"
}
