#!/bin/bash
# this script is for fastp to check quality
# author: laojp
# time: 2023.06.12
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_fastp [id] [input_dir] [output_dir] [input_suffix] [output_suffix]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_fastp {
    echo -e "\t\t[base function] process_fastp"
    if [ ! -f ${3}/${1}${5}.html ]; then
        case ${single_end} in
            "TRUE") 
                ${fastp}/fastp \
                    -i ${2}/${1}${4}.fq.gz \
                    -o ${3}/${1}${5}.fq.gz \
                    -j ${3}/${1}${5}.json \
                    -h ${3}/${1}${5}.html \
                    -w 8 \
                    1>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
            ;;
            "FALSE") 
                ${fastp}/fastp \
                    -i ${2}/${1}${4}_1.fq.gz -I ${2}/${1}${4}_2.fq.gz \
                    -o ${3}/${1}${5}_1.fq.gz -O ${3}/${1}${5}_2.fq.gz \
                    -j ${3}/${1}${5}.json \
                    -h ${3}/${1}${5}.html \
                    -w 8 \
                    1>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
            ;;
        esac
    fi
}
