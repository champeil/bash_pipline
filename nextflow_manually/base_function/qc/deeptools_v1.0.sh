#!/bin/bash
# this script is for using the deeptools for qcing the file
# time: 2023.06.19
# position: SYSUCC bioinformatic platform
# author: laojp
# usage
#  process_alignmentsieve [id] [input_dir] [output_dir] [suffix_input] [suffix_output]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_alignmentsieve {
    echo -e "\t\t[base function] process_alignmentsieve"
    if [ ! -f ${3}/${1}${5}.bam ]; then
        ${deeptools}/alignmentSieve \
            --numberOfProcessors 8 \
            --ATACshift \
            --bam ${2}/${1}${4}.bam \
            -o ${3}/${1}${4}_temp.bam \
            1>${3}/${1}${4}_alignmentsieve_log 2>&1 
        process_exist_status_background ${1}
        ${samtools}/samtools sort -@ 8 -O bam -o ${3}/${1}${5}.bam ${3}/${1}${4}_temp.bam 1>>${3}/${1}${4}_alignmentsieve_log 2>&1
        process_exist_status_background ${1}
        ${samtools}/samtools index -@ 8 ${3}/${1}${5}.bam 1>>${3}/${1}${4}_alignmentsieve_log 2>&1
        process_exist_status_background ${1}
        rm ${3}/${1}${4}_temp.bam
        process_exist_status_background ${1}
    fi
}