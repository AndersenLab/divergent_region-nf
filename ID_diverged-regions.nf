#!/usr/bin/env nextflow 

date = new Date().format( 'yyyyMMdd' )


params.help                   = null
params.config                 = null
params.out_base               = null
params.cpu                    = "4"
params.snv_vcf                = params.snv_vcf
params.c_sizes                = params.c_sizes
params.g_regions              = params.g_regions
params.window_size            = params.window_size
params.sv_fraction            = params.sv_fraction
params.vcfanno_config         = params.vcfanno_config
params.nsv_svcov              = params.nsv_svcov
params.cov_svcov              = params.cov_svcov
params.percent_outlier        = params.percent_outlier


if (params.help) {
    log.info '''
    ____   _                                         __     ____                _                  _   __ ______
   / __ \\ (_)_   __ ___   _____ ____ _ ___   ____   / /_   / __ \\ ___   ____ _ (_)____   ____     / | / // ____/
  / / / // /| | / // _ \\ / ___// __ `// _ \\ / __ \\ / __/  / /_/ // _ \\ / __ `// // __ \\ / __ \\   /  |/ // /_    
 / /_/ // / | |/ //  __// /   / /_/ //  __// / / // /_   / _, _//  __// /_/ // // /_/ // / / /  / /|  // __/    
/_____//_/  |___/ \\___//_/    \\__, / \\___//_/ /_/ \\__/  /_/ |_| \\___/ \\__, //_/ \\____//_/ /_/  /_/ |_//_/       
                             /____/                                  /____/                                     
'''
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow ID_diverged-regions.nf --out_base Analysis"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--out_base             String                Name of folder to output results"
    log.info "--bamdir               String                Name of folder where bam files are"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info "--snv_vcf              FILE                 Location to the small variant VCF to use for analysis"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=2)"
    log.info "--window_size          INTEGER              Size of window for variant and depth counts(default=1000)"
    log.info "--percent_outlier      INTEGER              Percent of bins for outlier variant count calculation"
    log.info "--sv_fraction          INTEGER              Fraction of window to be overlapped by structural variant to be considered for masking"
    log.info "--nsv_svcov            INTEGER              Number of SVs to be present in window for SV+Cov mask"
    log.info "--cov_svcov            INTEGER              Minimum Coverage for SV+Cov mask"
    log.info "--email                STRING               email address for job notifications"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "BEDtools               v2.27.1"
    log.info "R                      vXX"
    log.info "R-tidyverse            vXX"
    log.info "R-data.table           vXX"
    log.info "--------------------------------------------------------"    
    exit 1
} else {
    log.info ""
    log.info "Small Variant VCF                       = ${params.snv_vcf}"
    log.info "cpu                                     = ${params.cpu}"
    log.info "Output folder                           = ${params.out}"
    log.info "Window Size                             = ${params.window_size}"
    log.info "Percent outlier                         = ${params.percent_outlier}"   
    log.info "Window SV fraction                      = ${params.sv_fraction}"
    log.info "SV Count SV+Cov                         = ${params.nsv_svcov}"
    log.info "Coverage SV+Cov                         = ${params.cov_svcov}"
    log.info ""
}

/*
~ ~ ~ > * Define Contigs 
*/
CONTIG_LIST = ["I", "II", "III", "IV", "V", "X"]
Channel.from(CONTIG_LIST)
       .into{contigs_counts;
             contigs_windows;
             contigs_other}

// initialize VCF channel and spread to different channels for analysis
small_vcf = Channel.fromPath(params.snv_vcf)
small_index = Channel.fromPath("${params.snv_vcf}" + ".tbi")

// initialize chromsome sizes 
chromosome_lengths = Channel.fromPath(params.c_sizes)

// initialize chromsome sizes 
genomic_regions = Channel.fromPath(params.g_regions)


// grab bam files for coverage analysis
bamdir = params.bamdir

bams = Channel.fromPath( bamdir+'*.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${bamdir}" }
              .map {  path -> [ path.name.replace(".bam",""), path ] }

bais = Channel.fromPath( bamdir+'*.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${bamdir}" }
              .map { path -> [ path.name.replace(".bam.bai",""), path ] }

bams
  .join(bais)
  .set { bams_bais }

// grab strain SV files - tuple: [SM, SV bed]
sv_bin = params.sv_bin

sm_sv = Channel.fromPath( sv_bin+'*_SnpEff.bed' )
              .ifEmpty { error "Cannot find any SV file in: ${sv_bin}" }
              .map {  path -> [ path.name.replace("_SnpEff.bed",""), path ] }


// move this below VCF annotation when we decide to switch analysis to using masked VCF.

small_vcf
  .spread(small_index)
  .into { smallvcf_to_names;
          smallvcf_counts;
          smallvcf_masking;
          smallvcf_ldprune;
          }

process extract_sample_names {

  executor 'local'

  input:
    set file(vcf), file(index) from smallvcf_to_names

  output:
    file("sample_names.csv") into sample_names

  """
    bcftools query -l ${vcf} > sample_names.csv
  """

}


process make_windows {

  executor 'local'

  input:
    file(chrlen) from chromosome_lengths

  output:
    file("windows.bed") into genomic_windows
    file("windows.bed") into genomic_windows_coverage
    file("windows.bed") into genomic_windows_ale

  """
    bedtools makewindows -g ${chrlen} -w ${params.window_size} > windows.bed
  """

}

sample_names
    .splitCsv()
    .map{ row -> [row[0]] }
    .spread(smallvcf_counts)
    .spread(genomic_windows)
    .set { samples_to_counts }

process count_variants {

  tag {"${SM}"}

  publishDir "${params.out}/Strain_Counts", mode: "copy", pattern: "*_variant_counts.txt"

  input:
    set val(SM), file(vcf), file(index), file(windows) from samples_to_counts

  output:
    set val(SM), file("${SM}_variant_counts.txt") into sample_counts
    set val(SM), file("${SM}.vcf.gz") into sample_vcf


  """
    bcftools view -s ${SM} ${vcf} |\
    bcftools filter -i 'GT="alt"' -Oz -o ${SM}.vcf.gz
    bedtools coverage -a ${windows} -b ${SM}.vcf.gz -counts > ${SM}_variant_counts.txt
  """

}

sample_counts
  .spread(genomic_regions)
  .set{sample_counts_to_outliers}


process find_outliers {

  tag {"${SM}"}

  publishDir "${params.out}/Strain_Outliers/Data", mode: 'copy', pattern: "*Outlier_Counts.tsv"
  publishDir "${params.out}/Strain_Outliers/Plots", mode: 'copy', pattern: "*Outlier_Counts.png"

  input:
    set val(SM), file(SM_counts), file(gregions) from sample_counts_to_outliers

  output:
    set val(SM), file("*Outlier_Counts.tsv") into sample_outliers
    file("*Outlier_Counts.png") into sample_outliers_plot


  """
    Rscript --vanilla `which Classify_Outliers.R` ${SM} ${gregions}
  """

}

bams_bais
  .spread(genomic_windows_coverage)
  .into{bam_bai_window;
       bam_bai_ale}

process calculate_coverage {

  tag {"${SM}"}

  publishDir "${params.out}/Strain_Coverage/", mode: 'copy', pattern: "*regions.bed.gz"
  publishDir "${params.out}/Strain_Coverage/", mode: 'copy', pattern: "*regions.bed.gz.csi"

  input:
    set val(SM), file(smbam), file(smbam_index), file(windows) from bam_bai_window

  output:
    set val(SM), file(windows), file("${SM}.regions.bed.gz"), file("${SM}.regions.bed.gz.csi") into sample_coverage


  """
    mosdepth ${SM} ${smbam} -b ${windows}
  """

}

sample_coverage
  .join(sample_outliers)
  .join(sm_sv)
  .set{outlier_coverage}

/*
process calculate_ALE {

  tag {"${SM}"}

  input:
    set val(SM), file(smbam), file(smbam_index), file(windows) from bam_bai_ale

  output:
    set val(SM), file(windows), file("${SM}_ALE.bed") into sample_ale


  """
    ${params.ale} ${smbam} ${params.reference} ${SM}.ale

    grep -v "#" ${SM}.ale |\\
    awk 'BEGIN {FS=" ";OFS="\\t"} {\$1=\$1} {print \$1, \$2, \$2, \$3, \$4+\$5+\$6+\$7}' |\\
    sed 's/^0/I/g'|\\
    sed 's/^1/II/g' |\\
    sed 's/^2/III/g' |\\
    sed 's/^3/IV/g' |\\
    sed 's/^4/V/g' |\\
    sed 's/^5/X/g' |\\
    sed 's/^6/MtDNA/g' |\\
    grep -v MtDNA |\\
    awk '\$1 != 0 && \$2 != 0 {print}' > ${SM}_ALE.bed
  """

}

process window_ALE {

  memory '64 GB'
  cpus 8


  tag {"${SM}"}

  publishDir "${params.out}/Strain_ALE/", mode: 'copy', pattern: "*_window_ale.tsv"

  input:
    set val(SM), file(windows), file(smale) from sample_ale

  output:
    set val(SM), file("${SM}_window_ale.tsv") into sample_window_ale


  """
    bedtools intersect -loj -a ${windows} -b ${smale} > temp_ALEwindow.bed

    Rscript --vanilla `which window_ALE.R` temp_ALEwindow.bed ${SM}
  """
}


sample_coverage
  .join(sample_outliers)
  .join(sample_window_ale)
  .set{outlier_coverage}

*/



process outlier_coverage {

  tag {"${SM}"}

  publishDir "${params.out}/Strain_Coverage_Outliers/", mode: 'copy', pattern: "*_Processed_Outliers.tsv"

  input:
    set val(SM), file(windows), file(covbed), file(covind), file(outliers), file(sv_bed) from outlier_coverage

  output:
    set val(SM), file(windows), file(covbed), file(covind), file(outliers), file(sv_bed), file("${SM}_Processed_Outliers.tsv") into sm_processed_outliers
    file("${SM}_Processed_Outliers.tsv") into sm_processed_outliers_to_join


  """
    grep -e DEL -e INS ${sv_bed} |\\
    awk '\$3-\$2 < 2e4 {print}' |\
    awk '\$4 > 1 || \$5 == "INS" {print}' |\
    awk '\$10 == "1/1"{print}' |\
    awk -F"\t" '!seen[\$1, \$2, \$3, \$5]++' |\
    cut -f -3 |\
    bedtools coverage -a ${windows} -b stdin > ${SM}_nSVs.bed

    Rscript --vanilla `which Coverage_Outlier.R` ${SM} ${outliers} ${covbed} ${SM}_nSVs.bed
  """

}

sm_processed_outliers_to_join
  .toSortedList()
  .set{combined_outliers}

process combine_strain_outliers {

  publishDir "${params.out}/Strain_Coverage_Outliers/", mode: 'copy', pattern: "Population_Outlier.tsv"
  publishDir "${params.out}/Strain_Coverage_Outliers/", mode: 'copy', pattern: "Count_Threshold.tsv"

  memory '64 GB'

  input:
    file combined_outliers

  output:    
    file("Population_Outlier.tsv") into pop_processed_outliers_to_join
    file("Count_Threshold.tsv") into variant_count_threshold

  """
  head_file=`ls *Processed_Outliers.tsv | head -1`

  head -1 \$head_file > Population_Outlier.tsv; tail -n +2 -q *_Processed_Outliers.tsv >> Population_Outlier.tsv

  Rscript --vanilla `which Get_Variant_Count_Threshold.R` Population_Outlier.tsv ${params.percent_outlier}
  """

}

sm_processed_outliers
  .spread(variant_count_threshold)
  .set{sm_processed_outliers_threshold}


process generate_masks {

  tag {"${SM}"}

  publishDir "${params.out}/Strain_Masks/Plots", mode: 'copy', pattern: "*_masked_plot.pdf"
  publishDir "${params.out}/Strain_Masks/Data/Loose_Masks", mode: 'copy', pattern: "*_Mask_Regions.tsv"
  publishDir "${params.out}/Strain_Masks/Data/Processed_Masks", mode: 'copy', pattern: "*_Mask_DF.tsv"

  input:
    set val(SM), file(windows), file(covbed), file(covind), file(outliers), file(sv_bed), file(pr_outliers),file(threshold) from sm_processed_outliers_threshold

  output:
    set val(SM), file("${SM}_Mask_Regions.tsv") into masked_region
    file("${SM}_Mask_DF.tsv") into processed_masked_df
    file("${SM}_masked_plot.pdf") into processed_masked_plots

  """
    thresh=`head -1 Count_Threshold.tsv`

    Rscript --vanilla `which Make_Windows_DL.R` ${pr_outliers} ${params.sv_fraction} \$thresh ${params.nsv_svcov} ${params.cov_svcov}
  """

}

sample_vcf
  .join(masked_region)
  .spread(smallvcf_masking)
  .set{vcf_to_mask}


process mask_vcf {

  tag {"${SM}"}

  publishDir "${params.out}/VCF/STRAIN/", mode: 'copy', pattern: "*_Mask_Annotated.vcf.gz"

  input:
    set val(SM), file(vcf), file(masks), file(popvcf), file(popindex) from vcf_to_mask

  output:
    set val(SM), file("${SM}_Mask_Annotated.vcf.gz"), file("${SM}_Mask_Annotated.vcf.gz.tbi") into sm_mask_annotated_vcfs
    set file("${SM}_Mask_Hard_filter.vcf.gz") into sm_mask_hf_vcfs
    set file("${SM}_Mask_Hard_filter.vcf.gz.tbi") into sm_mask_hf_index

  """
    awk -F":" '\$1=\$1' OFS="\\t" ${masks} |\\
    awk -F"-" '\$1=\$1' OFS="\\t" |\\
    awk '{
    if (\$2 == 0 )
      print \$1, "1", \$3, "MASKED";
    else
      print \$0, "MASKED"
    }' OFS="\\t" |\\
    bgzip > Mask.bed.gz

    tabix Mask.bed.gz

    bcftools view -s ${SM} ${popvcf} -Oz -o temp_sm.vcf.gz
    tabix -p vcf temp_sm.vcf.gz

    vcfanno ${params.vcfanno_config} temp_sm.vcf.gz |\\
    bcftools view -Oz -o ${SM}_Mask_Annotated.vcf.gz

    tabix -p vcf ${SM}_Mask_Annotated.vcf.gz

    bcftools view ${SM}_Mask_Annotated.vcf.gz |\\
    bcftools filter -e 'INFO/Divergent_Mask="MASKED"' --set-GTs . -Oz -o ${SM}_Mask_Hard_filter.vcf.gz

    tabix -p vcf ${SM}_Mask_Hard_filter.vcf.gz
  """

}

sm_mask_hf_vcfs
  .toSortedList()
  .set{ merged_mask_vcf }

sm_mask_hf_index
  .toSortedList()
  .set{ merged_mask_index }

process merge_conservative {
    
    publishDir "${params.out}/VCF", mode: 'copy'

    memory '64 GB'
    cpus 8

    input:
      file merged_mask_vcf
      file merged_mask_index

    output:
      set file("WI.Masked.vcf.gz"), file("WI.Masked.vcf.gz.tbi") into joint_conservative_masked

    """
      bcftools merge -m none -Oz -o WI.Masked.vcf.gz ${merged_mask_vcf} 

      tabix -p vcf WI.Masked.vcf.gz
    """

}

process masked_only_vcf {

    tag {"${SM}"}

    input:
      set val(SM), file(vcf), file(index) from sm_mask_annotated_vcfs

    output:
      set val(SM), file("${SM}_Masked_Regions.vcf.gz"), file("${SM}_Masked_Regions.vcf.gz.tbi") into sm_masked_vcfs
      file("${SM}_Masked_Regions.vcf.gz") into sm_masked_vcfs_to_merge
      file("${SM}_Masked_Regions.vcf.gz.tbi") into sm_masked_vcf_index

    """
      bcftools view ${vcf} |\\
      awk '\$0 ~ "#" || \$0 ~ "Divergent_Mask=MASKED" {print}' |\\
      bcftools view -Oz -o ${SM}_Masked_Regions.vcf.gz

      tabix -p vcf ${SM}_Masked_Regions.vcf.gz
    """

}

sm_masked_vcfs_to_merge
  .toSortedList()
  .set{ merged_masked_vcf }

sm_masked_vcf_index
  .toSortedList()
  .set{ merged_masked_index }


process merge_masked {
  
  publishDir "${params.out}/VCF", mode: 'copy'

  memory '64 GB'
  cpus 8

  input:
    file merged_masked_vcf
    file merged_masked_index

  output:
    set file("WI.Masked.Region.vcf.gz"), file("WI.Masked.Region.vcf.gz.tbi") into joint_masked_vcf

  """
    bcftools merge -m none -Oz -o WI.Masked.Region.vcf.gz ${merged_masked_vcf} 

    tabix -p vcf WI.Masked.Region.vcf.gz
  """

}

LD_CUTOFFS = [0.1, 0.2, 0.4, 0.6, 0.8]
Channel.from(LD_CUTOFFS)
  .spread(smallvcf_ldprune)
  .set{smallvcf_to_prune}


process ld_prune_vcf {

  tag {"${ld}"}

  publishDir "${params.out}/VCF", mode: 'copy'

  input:
    set val(ld), file(vcf), file(index) from smallvcf_to_prune

  output:
    set file("WI.LD_${ld}.vcf.gz"), file("WI.LD_${ld}.vcf.gz.tbi") into ldpruned_vcf

  when:
    params.ld_vcf

  """
      plink --vcf ${vcf} --snps-only --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 ${ld} --allow-extra-chr 
      
      awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | sort -k1,1d -k2,2n > markers.txt

      bcftools query -l ${vcf} | sort > sorted_samples.txt 

      bcftools view -v snps -S sorted_samples.txt -R markers.txt -Oz -o WI.LD_${ld}.vcf.gz ${vcf} 

      tabix -p vcf WI.LD_${ld}.vcf.gz
  """

}