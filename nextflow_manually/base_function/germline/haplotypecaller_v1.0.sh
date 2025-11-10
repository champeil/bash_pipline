#!/bin/bash
# this script is for call germline with haplotypecaller_GVCF mode
# author: laojp
# time: 2023.05.31
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: 
#  process_germline_haplotypecaller_gvcf [id] [input_dir] [output_dir] [input_suffix] [output_suffix]
#  process_germline_combinevcf [input_dir] [output_dir] [input_suffix]
#  process_germline_genotypegvcf [input_dir] [output_dir]
#  process_germline_vqsr [input_dir] [output_dir] [output_snp] [output_indel]
#  process_germline_vqsr [snp] [indel] [input_dir] [output_dir]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_germline_haplotypecaller_gvcf {
    echo -e "\t\t[base function] process_germline_haplotypecaller_gvcf"
    if [ ! -f ${3}/${1}${5}.g.vcf.gz ]; then
        ${gatk}/gatk HaplotypeCaller \
            -R ${reference_fasta} \
            -I ${2}/${1}${4}.bam \
            -O ${3}/${1}${5}.g.vcf.gz \
            -ERC GVCF \
            1>${3}/${1}${5}_log 2>&1
        process_exist_status_background ${1}
    fi
}

function process_germline_combinevcf {
    echo -e "\t\t[base function] process_gemline_combinevcf"
    if [ ! -f ${2}/cohort.g.vcf.gz ]; then
        ${gatk}/gatk CombineGVCFs \
            -R ${reference_fasta} \
            $(ls ${1}/*${3}*.g.vcf.gz | xargs -I {} echo "--variant {}") \
            -O ${2}/cohort.g.vcf.gz \
            1>${2}/combine_gvcf_log 2>&1 
        process_exist_status_background "combinegvcf"
    fi
}

function process_germline_genotypegvcf {
    echo -e "\t\t[base function] process_germline_genotypegvcf"
    if [ ! -f ${2}/genotype.vcf.gz ]; then 
        ${gatk}/gatk GenotypeGVCFs \
            -R ${reference_fasta} \
            -V ${1}/cohort.g.vcf.gz \
            -O ${2}/genotype.vcf.gz \
        1>${2}/genotype_log 2>&1 
        process_exist_status_background "genotype"
    fi
}

function process_germline_vqsr {
    echo -e "\t\t[base function] process_germline_vqsr"
    if [ ! -f ${2}/${3} ]; then
        ${gatk}/gatk VariantRecalibrator \
            -R ${reference_fasta} \
            -V ${1}/genotype.vcf.gz \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${snp_1000G} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            --an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            -O ${2}/germline_snp_recal \
            --tranches-file ${2}/germline_snp_tranches \
            --rscript-file ${2}/germline_snp_plot.R \
            1>${2}/vqsr_snp_log 2>&1
        process_exist_status_background "vqsr_snp"
        ${gatk}/gatk ApplyVQSR \
            -R ${reference_fasta} \
            -V ${1}/genotype.vcf.gz \
            -O ${2}/${3} \
            -ts-filter-level 99.0 \
            --tranches-file ${2}/germline_snp_tranches \
            --recal-file ${2}/germline_snp_recal \
            -mode SNP \
            1>>${2}/vqsr_snp_log 2>&1
        process_exist_status_background "vqsr_snp_apply"
    fi
    if [ ! -f ${2}/${4} ]; then
        ${gatk}/gatk VariantRecalibrator \
            -R ${reference_fasta} \
            -V ${1}/genotype.vcf.gz \
            --resource:mills,known=true,training=true,truth=true,prior=12.0 ${mill} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            --an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode INDEL \
            -O ${2}/germline_indel_recal \
            --tranches-file ${2}/germline_indel_tranches \
            --rscript-file ${2}/germline_indel_plot.R \
            1>${2}/vqsr_indel_log 2>&1
        process_exist_status_background "vqsr_indel"
        ${gatk}/gatk ApplyVQSR \
            -R ${reference_fasta} \
            -V ${1}/genotype.vcf.gz \
            -O ${2}/${4} \
            -ts-filter-level 99.0 \
            --tranches-file ${2}/germline_indel_tranches \
            --recal-file ${2}/germline_indel_recal \
            -mode INDEL \
            1>>${2}/vqsr_indel_log 2>&1
        process_exist_status_background "vqsr_indel_apply"
    fi
}

function process_combine_germline {
    if [ ! -f ${4}/combine_germline.vcf ]; then
        ${gatk}/gatk MergeVcfs \
            -I ${3}/${1} \
            -I ${3}/${2} \
            -O ${4}/combine_germline.vcf \
            1>${4}/combine_germline_log 2>&1
        process_exist_status_background "combine_germline"
    fi
}
