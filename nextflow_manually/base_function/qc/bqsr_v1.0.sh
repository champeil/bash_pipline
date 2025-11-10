#!/bin/bash
# this script is for wes seq: bqsr qc
# time: 2023.06.16
# position: SYSUCC bioinformatic platform
# author: laojp
# usage: process_bqsr [id] [input_dir] [output_dir] [suffix_input] [suffix_output]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_bqsr {
    echo -e "\t\t[base function] process_bqsr"
    if [ ! -f ${3}/${1}${5}.bam ]; then
        ${gatk}/gatk BaseRecalibrator \
            -R ${reference_fasta} \
            -I ${2}/${1}${4}.bam \
            --known-sites ${dbsnp} \
            --known-sites ${gold_indel} \
            -O ${3}/${1}_recal.table \
            1>${3}/${1}_log 2>&1
        process_exist_status_background ${1}
        ${gatk}/gatk ApplyBQSR \
            -R ${reference_fasta} \
            -I ${2}/${1}${4}.bam \
            -bqsr ${3}/${1}_recal.table \
            --create-output-bam-index true \
            -O ${3}/${1}${5}.bam \
            1>>${3}/${1}_log 2>&1
        process_exist_status_background ${1}
    fi
}