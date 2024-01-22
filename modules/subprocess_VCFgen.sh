source subprocess_utilities_VCFgen.sh

main() {
    declare -A passed_variables
    passed_variables["vcf_ulid"]=${1}
    passed_variables["current_line_name"]=${2}
    passed_variables["input_dir"]=${3}
    passed_variables["wt_input"]=${4}
    passed_variables["mu_input"]=${5}
    passed_variables["pairedness"]=${6}
    passed_variables["output_dir_path"]=${7}
    passed_variables["output_prefix"]=${8}
    passed_variables["vcf_table_path"]=${9}
    passed_variables["reference_dir"]=${10}
    passed_variables["reference_genome_name"]=${11}
    passed_variables["snpEff_species_db"]=${12}
    passed_variables["reference_genome_source"]=${13}
    passed_variables["known_snps_path"]=${14}
    passed_variables["threads_limit"]=${15}
    passed_variables["call_variants_in_parallel"]=${16}
    passed_variables["cleanup"]=${17}

    ## Printing all assigned variables for logging purposes
    print_variable_info passed_variables "Variables passed to subprocess_VCFgen.sh"
    ## Assigning variables in array to their respective keys
    assign_values passed_variables

    ## Generate other needed variables 
    # ...get rid of the / paths! pass them with os.path.join for consistancy...
    declare -A generated_variables
    generated_variables["reference_genome_path"]="${reference_dir}/${reference_genome_name}.fa"
    generated_variables["reference_chrs_path"]="${reference_dir}/${reference_genome_name}.chrs"
    generated_variables["reference_chrs_fa_path"]="${reference_dir}/${reference_genome_name}.chrs.fa"
    generated_variables["snpeff_dir"]="${output_dir_path}/snpEff"
    generated_variables["snpeff_out_filename"]="${output_dir_path}/snpEff/${vcf_ulid}-_${current_line_name}"
    
    # Bulktupe is to be used to iterate through the wt and mu bulk files, 
    # to ensure easy naming consistancy... Hope this isn't just harder to read
    bulktype=("wt" "mu")     
    for bulk in "${bulktype[@]}"; do
        generated_variables["bwa_output_sam_${bulk}"]="${output_prefix}_${bulk}.sam"
        generated_variables["samtools_output_bam_${bulk}"]="${output_prefix}_${bulk}.bam"
        generated_variables["samtools_fixmate_output_${bulk}"]="${output_prefix}_${bulk}.fix.bam"
        generated_variables["samtools_sortsam_output_${bulk}"]="${output_prefix}_${bulk}.sort.bam"
        generated_variables["picard_markduplicates_output_${bulk}"]="${output_prefix}_${bulk}.sort.md.bam"
        generated_variables["picard_addorreplacereadgroups_output_${bulk}"]="${output_prefix}_${bulk}.sort.md.rg.bam"
        generated_variables["picard_buildbamindex_output_${bulk}"]="${output_prefix}_${bulk}.sort.md.rg.bai"
    done

    generated_variables["gatk_haplotypecaller_output"]="${output_prefix}.hc.vcf"
    generated_variables["snpeff_output"]="${output_prefix}.se.vcf"
    generated_variables["snpsift_output"]=${output_prefix}.snpsift.table.tmp
    generated_variables["tmp_table_file_name"]="${output_prefix}.table.tmp"

    ## Printing all generated variables for logging purposes
    print_variable_info generated_variables "Variables generated in subprocess_VCFgen.sh"
    ## Assigning variables in array to their respective keys
    assign_values generated_variables

    ## Prepare references and directory structure
    print_message "Preparing references and directory structure"
    create_directories ${output_dir} ${snpeff_dir} ${reference_dir}
    
    download_reference_genome "${reference_genome_path}" "${reference_genome_source}"
    
    create_chrs_file "${reference_genome_path}" "${reference_chrs_fa_path}"
    create_fai_and_index ${reference_genome_path} "${reference_chrs_fa_path}"
    create_sequence_dictionary "${reference_chrs_fa_path}" "${reference_chrs_path}"

    echo "References and directories prepared. Proceeding with mapping...."

       print_message "Mapping"
    # Align reads using BWA. A more modern aligner for this may be implemented 
    # sometime  
    for bulk in "${bulktype[@]}"; do
        bwa mem \
            -t "$threads_halfed" \
            -M "${reference_chrs_fa_path}" \
            "bwa_output_sam_${bulk}" > "bwa_output_sam_${bulk}"
    done

    # Create binary alignment map for more efficient processing
    print_message "Converting sam to bam"
    for bulk in "${bulktype[@]}"; do
        samtools view \
            -bSh \
            -@ "$threads_halfed" \
            "bwa_output_sam_${bulk}" > "samtools_output_bam_${bulk}" &
    done
    
    echo "..."
    wait

    # Fix paired-end
    # Ensures that mapping information between read pairs is accurate. 
    for bulk in "${bulktype[@]}"; do
        if [ "${pairedness}" == "paired-end" ]; then
            print_message "Reads are paired-end. Running samtools fixmate"
            samtools fixmate \
            "samtools_output_bam_${bulk}" \
            "samtools_fixmate_output_${bulk}" &
        fi
    done
    
    echo "..."
    wait

    # SortSam
    # Sorting ensures that reads are organized in genomic order. GATK haplotype
    # caller reassembles variant regions de novo - it needs regions to be ordered. 
    # Reassembly makes HC slow, but it is pretty accurate. 
    print_message "Sorting by coordinate"
    for bulk in "${bulktype[@]}"; do
        picard SortSam \
            -I "samtools_fixmate_output_${bulk}" \
            -O "samtools_sortsam_output_${bulk}" \
            -SORT_ORDER coordinate &
    done
    wait

    # Mark duplicates. Reads with identical start and stop positions, and are 
    # formed during PCR amplification. Marking them allows accurate assessment 
    # of read depth, in that only unique reads at variants are counted.
    print_message "Marking duplicates"
    for bulk in "${bulktype[@]}"; do
        picard MarkDuplicates \
            -I "samtools_sortsam_output_${bulk}" \
            -O "picard_markduplicates_output_${bulk}" \
            -METRICS_FILE "${output_prefix}_${bulk}.metrics.txt" \
            -ASSUME_SORTED true &
    done
    wait

    # Format headers so BAMs can be fed through GATK haplotype caller
    print_message "Adding header for GATK"
    for bulk in "${bulktype[@]}"; do
        picard AddOrReplaceReadGroups \
            -I "picard_markduplicates_output_${bulk}" \
            -O "picard_addorreplacereadgroups_output_${bulk}" \
            -RGLB "${current_line_name}_${bulk}" \
            -RGPL illumina \
            -RGSM "${current_line_name}_${bulk}" \
            -RGPU run1 \
            -SORT_ORDER coordinate &
    done
    wait

    print_message "Building BAM index"
    # Build BAM index
    # Needed to run haplotype caller. Increases the speed of accessing and 
    # retrieving data within genomic regions during variant calling. Allows GATK HC
    # to skip directly to the region of interest.  
    for bulk in "${bulktype[@]}"; do
        picard BuildBamIndex \
            -INPUT "picard_addorreplacereadgroups_output_${bulk}" \
            -O "picard_buildbamindex_output_${bulk}" &
    done
    wait

    print_message "Calling haplotypes. This may take a while..."

    # GATK HC Variant calling
    # Haplotype caller looks for regions with variance and locally reconstructs
    # the region using the available reads and calls variants for the region. 
    # Time-consuming but accurate
    
    if [ $call_variants_in_parallel = True ]; then
        split_and_call_haplotype $output_prefix $output_dir_path $reference_chrs_fa_path $threads_limit
    else
        gatk HaplotypeCaller \
            -R "$reference_chrs_fa_path" \
            -I "$picard_addorreplacereadgroups_output_mu" \
            -I "$picard_addorreplacereadgroups_output_wt" \
            -O "$gatk_haplotypecaller_output" \
            -output-mode EMIT_ALL_CONFIDENT_SITES \
            --native-pair-hmm-threads "$threads_limit"
    fi

    print_message "SnpEff: Labeling SNPs with annotations and potential impact on gene function"
    
    # snpEff, labeling SNPs
    # snpEff labels the variants in haplotype caller with likely impacts of variants 
    # on gene function (early stop/start codons, missense mutations, etc.) using
    # databases assembled from annotated reference files 
    # (gff, transcriptomes, and genomes). 
    snpEff "$snpEff_species_db" \
        -s "$snpeff_out_filename" \
        "$gatk_haplotypecaller_output" > "${snpeff_output}"
    wait
    print_message "Haplotypes called and SNPs labeled. Cleaning data."

    # Extracting SNPeff data and variant information into a table
    # built-in data processing of snpEff labels... 
    extract_fields_snpSift \
        "$snpeff_output" \
        "$snpsift_output" \
        "$current_line_name"
    
    # Clean up data table
    remove_repetitive_nan "$snpsift_output" "$tmp_table_file_name"
    format_fields "$tmp_table_file_name" "$current_line_name"
    remove_complex_genotypes "$tmp_table_file_name"
    
    # Remove known SNPs
    remove_known_snps "$known_snps_path" "$tmp_table_file_name" "$vcf_table_path"
    
    # Add headers
    add_headers "$vcf_table_path"

    # Clean up temporary files
    cleanup_files "$output_dir_path" "$cleanup"
}

# Execute the main function
main "${@}"
