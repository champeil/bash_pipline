#!/bin/bash
# this script is for markduplicated to remove duplicate
# time: 2023.06.13
# position: SYSUCC bioinformatic platform
# author: laojp
# usage: process_markduplicate [id] [input_dir] [output_dir] [suffix_input] [suffix_output]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_markduplicate {
    echo -e "\t\t[base function] process_markduplicate"
    if [ ! -f ${3}/${1}${5}.bam ]; then
        ${gatk}/gatk --java-options "-Xmx200G -Djava.io.tmpdir=${3}/${1}_temp" MarkDuplicates \
            -I ${2}/${1}${4}.bam \
            --REMOVE_DUPLICATES true \
            --CREATE_INDEX true \
            -O ${3}/${1}${5}.bam \
            -M ${3}/${1}.metrics \
            1>${3}/${1}${5}_log 2>&1
        process_exist_status_background ${1}
    fi
}