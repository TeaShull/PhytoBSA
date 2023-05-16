#!/usr/bin/env bash

#list your vcf files to be intersected
vcfs="$@"

#index and tabix the vcfs
for files in $vcfs
  do
  bgzip "$files.hc.vcf"
  tabix "$files.hc.vcf.gz"
  vcfs_filelist+="./VCFs/$files.hc.vcf.gz "
done

bcftools isec -n=+2 -c all -o ws-0_all_common_SNPs.vcf $vcfs_filelist

