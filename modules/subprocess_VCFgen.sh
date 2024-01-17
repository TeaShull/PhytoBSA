#!/usr/bin/env bash
source subprocess_utilities_VCFgen.sh

main() {
    declare -A passed_variables
    passed_variables["vcf_ulid"]=${1}
    passed_variables["current_line_name"]=${2}
    passed_variables["allele"]=${3}
    passed_variables["input_dir"]=${4}
    passed_variables["wt_input"]=${5}
    passed_variables["mu_input"]=${6}
    passed_variables["pairedness"]=${7}
    passed_variables["output_dir_path"]=${8}
    passed_variables["output_prefix"]=${9}
    passed_variables["vcf_table_path"]=${10}
    passed_variables["reference_dir"]=${11}
    passed_variables["reference_genome_name"]=${12}
    passed_variables["snpEff_species_db"]=${13}
    passed_variables["reference_genome_source"]=${14}
    passed_variables["known_snps_path"]=${15}
    passed_variables["threads_limit"]=${16}
    passed_variables["cleanup"]=${17}

    ## Printing all assigned variables for logging purposes
    print_variable_info passed_variables "Variables passed to VCFgen.sh"
    ## Assigning variables in array to their respective keys
    assign_values passed_variables

    ## Generate other needed variables
    declare -A generated_variables
    generated_variables["reference_genome_path"]=${reference_dir}/${reference_genome_name}.fa
    generated_variables["reference_chrs_path"]=${reference_dir}/${reference_genome_name}.chrs
    generated_variables["reference_chrs_fa_path"]=${reference_dir}/${reference_genome_name}.chrs.fa
    generated_variables["snpeff_dir"]=${output_dir_path}/snpEff
    generated_variables["snpeff_out_filename"]=${output_dir_path}/snpEff/${vcf_ulid}-_${current_line_name}
    generated_variables["ems_file_name"]="${output_prefix}.ems.table"

    ## Printing all generated variables for logging purposes
    print_variable_info generated_variables "Variables generated in VCFgen.sh"
    ## Assigning variables in array to their respective keys
    assign_values generated_variables

    ## Prepare references and directory structure
    print_message "Preparing references and directory structure"
    create_directories ${output_dir} ${snpeff_dir} ${reference_dir}
    
    download_reference_genome "${reference_genome_path}" "${reference_genome_source}"
    
    create_chrs_file "${reference_genome_path}" "${reference_chrs_fa_path}"
    create_fai_and_index ${reference_genome_path} "${reference_chrs_fa_path}"
    create_sequence_dictionary "${reference_chrs_fa_path}" "${reference_chrs_path}"

    echo "Refrences and directories prepared. Proceeding with mapping...."

    print_message "Mapping"
    # Align reads using BWA. A more modern aligner for this may be implemented 
    # sometime  
    bwa mem \
        -t "$threads_halfed" \
        -M "${reference_chrs_fa_path}" \
        $wt_input > "${output_prefix}_wt.sam" &
    
    bwa mem \
        -t "$threads_halfed" \
        -M "${reference_chrs_fa_path}" \
        $mu_input > "${output_prefix}_mu.sam"
    wait

    #Create binary alignment map for more effecient processing
    print_message "Converting sam to bam"
    samtools view \
        -bSh \
        -@ "$threads_halfed" \
        "${output_prefix}_mu.sam" > "${output_prefix}_mu.bam" &
    samtools view \
        -bSh \
        -@ "$threads_halfed" \
        "${output_prefix}_wt.sam" > "${output_prefix}_wt.bam"
    
    echo "..."
    wait

    # Fix paired-end
    # Ensures that mapping information between read pairs is accurate. 
    if [ "${pairedness}" == "paired-end" ]; then
        print_message "Reads are paired-end. Running samtools fixmate"
        samtools fixmate "${output_prefix}_wt.bam" "${output_prefix}_wt.fix.bam" &
        samtools fixmate "${output_prefix}_mu.bam" "${output_prefix}_mu.fix.bam"
    fi
    
    echo "..."
    wait

    #SortSam
    # Sorting ensures that reads are organized in genomic order. GATK haplotype
    # caller reassembles variant regions denovo - it needs regions to be ordered. 
    # Reassembly makes HC slow, but it is pretty accurate. 
    print_message "Sorting by coordinate"
    picard SortSam \
        -I "${output_prefix}_mu.fix.bam" \
        -O "${output_prefix}_mu.sort.bam" \
        -SORT_ORDER coordinate &
    picard SortSam \
        -I "${output_prefix}_wt.fix.bam" \
        -O "${output_prefix}_wt.sort.bam" \
        -SORT_ORDER coordinate
    wait

    # Mark duplicates. Reads with identical start and stop positions, and are 
    # formed during PCR amplification. Marking them allows accurate assesment 
    # of read depth, in that only unique reads at variants are counted.
    print_message "Marking duplicates"
    picard MarkDuplicates \
        -I "${output_prefix}_mu.sort.bam" \
        -O "${output_prefix}_mu.sort.md.bam" \
        -METRICS_FILE "${output_prefix}_mu.metrics.txt" \
        -ASSUME_SORTED true &
    picard MarkDuplicates \
        -I "${output_prefix}_wt.sort.bam" \
        -O "${output_prefix}_wt.sort.md.bam" \
        -METRICS_FILE "${output_prefix}_wt.metrics.txt" \
        -ASSUME_SORTED true
    wait

    # Format headers so BAMs can be fed through GATK haplotype caller
    print_message "Adding header for GATK"
    picard AddOrReplaceReadGroups \
        -I "${output_prefix}_mu.sort.md.bam" \
        -O "${output_prefix}_mu.sort.md.rg.bam" \
        -RGLB "${current_line_name}_mu" \
        -RGPL illumina \
        -RGSM "${current_line_name}_mu" \
        -RGPU run1 \
        -SORT_ORDER coordinate &
    picard AddOrReplaceReadGroups \
        -I "${output_prefix}_wt.sort.md.bam" \
        -O "${output_prefix}_wt.sort.md.rg.bam" \
        -RGLB "${current_line_name}_wt" \
        -RGPL illumina \
        -RGSM "${current_line_name}_wt" \
        -RGPU run1 \
        -SORT_ORDER coordinate
    wait

    print_message "Building BAM index"
    # Build BAM index
    # Needed to run haplotyple caller. Increases the speed of accessing and 
    # retrieving data within genomic regions during variant calling. Allows GATK HC
    # to skip directly to the region of interest.  
    picard BuildBamIndex \
        -INPUT "${output_prefix}_mu.sort.md.rg.bam" \
        -O "${output_prefix}_mu.sort.md.rg.bai" &
    picard BuildBamIndex \
        -INPUT "${output_prefix}_wt.sort.md.rg.bam" \
        -O "${output_prefix}_wt.sort.md.rg.bai"
    wait

    print_message "Calling haplotypes. This may take a while..."
    # GATK HC Variant calling
    # Haplotype caller looks for regions with varience and locally reconstructs
    # the region using the available reads, and calls variants for the region. 
    # Time consuming but accurate
    gatk HaplotypeCaller \
        -R "$reference_chrs_fa_path" \
        -I "${output_prefix}_mu.sort.md.rg.bam" \
        -I "${output_prefix}_wt.sort.md.rg.bam" \
        -O "${output_prefix}.hc.vcf" \
        -output-mode EMIT_ALL_CONFIDENT_SITES \
        --native-pair-hmm-threads "$threads_limit"

    print_message "SnpEff: Labeling SNPs with annotations and potential impact on gene function"
    
    # snpEff, labeling SNPs
    # snpEff labels the variants in haplotype caller with likely impact of variants 
    # on gene function (early stop/start codons, missense mutations, exc) using
    # databases assembled from annotated reference files 
    # (gff, transcriptomes and genomes). 
    snpEff "$snpEff_species_db" \
        -s "${snpeff_out_filename}" \
        "${output_prefix}.hc.vcf" > "${output_prefix}.se.vcf"
    wait
    print_message "Haplotypes called and SNPs labeled. Cleaning data."

    
    # Extracting SNPeff data and variant information into a table
    # built in data processing of snpEff labels... 
    extract_fields_snpSift \
        "${output_prefix}.se.vcf" \
        "${output_prefix}.table.tmp" \
        "${current_line_name}"

    # Clean up data table
    remove_repetitive_nan "${output_prefix}.table.tmp"
    filter_ems_mutations "${output_prefix}.table.tmp" "${output_prefix}.ems.table.tmp"
    filter_genotypes "${allele}" "${output_prefix}.ems.table.tmp" "${ems_file_name}"
    format_ems_file "${ems_file_name}" "${current_line_name}"
    remove_complex_genotypes "${ems_file_name}"

    # Get rid of chloroplastic and mitochondrial polymorphisms.
    remove_nongenomic_polymorphisms "$current_line_name" "${ems_file_name}"
    
    # Remove known SNPs
    remove_known_snps "$known_snps_path" "${output_prefix}.ems.table" "$vcf_table_path"
    
    # Add headers
    add_headers "$vcf_table_path"

    # Clean up temporary files
    cleanup_files "${output_dir_path}" "${cleanup}"
}

# Execute the main function
main "${@}"