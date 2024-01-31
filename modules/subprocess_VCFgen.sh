source ./subprocess_utilities_VCFgen.sh

main() {
    declare -A passed_variables=(
        ["vcf_ulid"]=${1}
        ["current_line_name"]=${2}
        ["wt_input"]=${3}
        ["mu_input"]=${4}
        ["pairedness"]=${5}
        ["output_prefix"]=${6}
        ["snpeff_dir"]=${7}
        ["snpeff_out_path"]=${8}
        ["snpsift_out_path"]=${9}
        ["reference_chrs_fa_path"]=${10}
        ["snpEff_species_db"]=${11}
        ["threads_limit"]=${12}
        ["call_variants_in_parallel"]=${13}
        ["cleanup"]=${14}
    )
    print_variable_info passed_variables "Variables passed to subprocess_VCFgen.sh"
    assign_values passed_variables #assign variables based on keys for easy access

    declare -A generated_variables=(
        ["threads_halfed"]=$((threads_limit / 2))
        ["bwa_output_sam_wt"]="${output_prefix}_wt.sam"
        ["bwa_output_sam_mu"]="${output_prefix}_mu.sam"
        ["samtools_output_bam_wt"]="${output_prefix}_wt.bam"
        ["samtools_output_bam_mu"]="${output_prefix}_mu.bam"
        ["samtools_fixmate_output_wt"]="${output_prefix}_wt.fix.bam"
        ["samtools_fixmate_output_mu"]="${output_prefix}_mu.fix.bam"
        ["samtools_sortsam_output_wt"]="${output_prefix}_wt.sort.bam"
        ["samtools_sortsam_output_mu"]="${output_prefix}_mu.sort.bam"
        ["picard_markduplicates_output_wt"]="${output_prefix}_wt.sort.md.bam"
        ["picard_markduplicates_output_mu"]="${output_prefix}_mu.sort.md.bam"
        ["picard_addorreplacereadgroups_output_wt"]="${output_prefix}_wt.sort.md.rg.bam"
        ["picard_addorreplacereadgroups_output_mu"]="${output_prefix}_mu.sort.md.rg.bam"
        ["picard_buildbamindex_output_wt"]="${output_prefix}_wt.sort.md.rg.bai"
        ["picard_buildbamindex_output_mu"]="${output_prefix}_mu.sort.md.rg.bai"
        ["gatk_haplotypecaller_output"]="${output_prefix}.hc.vcf"
        ["snpeff_output"]="${output_prefix}.se.vcf"
    )
    print_variable_info generated_variables "Variables generated in subprocess_VCFgen.sh"
    assign_values generated_variables

    create_fai_and_index "${reference_genome_path}" "${reference_chrs_fa_path}"
    create_sequence_dictionary "${reference_chrs_fa_path}"
    echo "References and directories prepared. Proceeding with mapping...."

    print_message "Mapping"
    # Align reads using BWA. A more modern aligner for this may be implemented 
    # sometime  
    bwa mem \
        -t $threads_halfed \
        -M $reference_chrs_fa_path \
        $wt_input > $bwa_output_sam_wt & 
    bwa mem \
        -t $threads_halfed \
        -M $reference_chrs_fa_path \
        $mu_input > $bwa_output_sam_mu

    echo "..."
    wait

    print_message "Converting sam to bam"
    #Create binary alignment map for more effecient processing
    samtools view \
        -bSh \
        -@ $threads_halfed \
        $bwa_output_sam_mu > $samtools_output_bam_mu & 
    samtools view \
        -bSh \
        -@ "$threads_halfed" \
        $bwa_output_sam_wt > $samtools_output_bam_wt

    echo "..."
    wait

    print_message "Fix paired-end reads"
    # Ensures that mapping information between read pairs is accurate. 
    if [ "${pairedness}" == "paired-end" ]; then
        samtools fixmate \
            $samtools_output_bam_mu \
            $samtools_fixmate_output_mu & 
        samtools fixmate \
            $samtools_output_bam_wt \
            $samtools_fixmate_output_wt
    fi

    echo "..."
    wait

    print_message "Sorting by coordinate"
    # Sorting ensures that reads are organized in genomic order. GATK haplotype
    # caller reassembles variant regions denovo - it needs regions to be ordered. 
    # Reassembly makes HC slow, but it is pretty accurate. 
    picard SortSam \
        -I $samtools_fixmate_output_mu \
        -O $samtools_sortsam_output_mu \
        -SORT_ORDER coordinate &
    picard SortSam \
        -I $samtools_fixmate_output_wt \
        -O $samtools_sortsam_output_wt \
        -SORT_ORDER coordinate
    wait

    print_message "Marking duplicates"
    # Mark duplicates. Reads with identical start and stop positions, and are 
    # formed during PCR amplification. Marking them allows accurate assesment 
    # of read depth, in that only unique reads at variants are counted.
    picard MarkDuplicates \
        -I $samtools_sortsam_output_mu \
        -O $picard_markduplicates_output_mu \
        -METRICS_FILE "${output_prefix}_mu.metrics.txt" \
        -ASSUME_SORTED true &
    picard MarkDuplicates \
        -I $samtools_sortsam_output_wt \
        -O $picard_markduplicates_output_wt \
        -METRICS_FILE "${output_prefix}_wt.metrics.txt" \
        -ASSUME_SORTED true
    wait

    print_message "Adding header for GATK"
    # Format headers so BAMs can be fed through GATK haplotype caller
    picard AddOrReplaceReadGroups \
        -I $picard_markduplicates_output_mu \
        -O $picard_addorreplacereadgroups_output_mu \
        -RGLB "${current_line_name}_mu" \
        -RGPL illumina \
        -RGSM "${current_line_name}_mu" \
        -RGPU run1 \
        -SORT_ORDER coordinate &
    picard AddOrReplaceReadGroups \
        -I $picard_markduplicates_output_wt \
        -O $picard_addorreplacereadgroups_output_wt \
        -RGLB "${current_line_name}_wt" \
        -RGPL illumina \
        -RGSM "${current_line_name}_wt" \
        -RGPU run1 \
        -SORT_ORDER coordinate
    wait

    print_message "Building BAM index"
    # Needed to run haplotyple caller. Increases the speed of accessing and 
    # retrieving data within genomic regions during variant calling. Allows GATK HC
    # to skip directly to the region of interest.  
    picard BuildBamIndex \
        -INPUT $picard_addorreplacereadgroups_output_mu \
        -O $picard_buildbamindex_output_mu &
    picard BuildBamIndex \
        -INPUT $picard_addorreplacereadgroups_output_wt \
        -O $picard_buildbamindex_output_wt
    wait

    print_message "Calling haplotypes. This may take a while..."
    # Haplotype caller looks for regions with varience and locally reconstructs
    # the region using the available reads, and calls variants for the region. 
    # Time consuming but accurate. Calling variants in parallel allows faster
    # processing time. 
    if [ $call_variants_in_parallel = True ]; then
        split_and_call_haplotype $output_prefix $output_dir_path $reference_chrs_fa_path $threads_limit
    else
        gatk HaplotypeCaller \
            -R $reference_chrs_fa_path \
            -I $picard_addorreplacereadgroups_output_mu \
            -I $picard_addorreplacereadgroups_output_wt \
            -O $gatk_haplotypecaller_output \
            -output-mode EMIT_ALL_CONFIDENT_SITES \
            --native-pair-hmm-threads "$threads_limit"
    fi

    print_message "SnpEff: Labeling SNPs with annotations and potential impact on gene function"
    # snpEff labels the variants in haplotype caller with likely impact of variants 
    # on gene function (early stop/start codons, missense mutations, exc) using
    # databases assembled from annotated reference files 
    # (gff, transcriptomes and genomes). 
    snpEff $snpEff_species_db \
        -s $snpeff_out_path \
        $gatk_haplotypecaller_output > $snpeff_output
    wait
    print_message "Haplotypes called and SNPs labeled. Cleaning data."

    # Extracting SNPeff data and variant information into a table
    # built in data processing of snpEff labels... 
    extract_fields_snpSift \
        $snpeff_output \
        $snpsift_out_path \
        $current_line_name
}

# Execute the main function
main "${@}"
