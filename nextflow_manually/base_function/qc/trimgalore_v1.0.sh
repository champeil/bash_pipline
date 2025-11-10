#!/bin/bash
# this script is for trimgalore the raw file
# author: laojp
# time: 2023.03.31
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_trimgalore [sample] [input_dir] [output_dir] [suffix_input] [sufffix_output]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_trimgalore {
    echo -e "\t\t[base function] process_trimgalore"
    case ${single_end} in 
        "TRUE")
            if [ ! -f ${3}/${1}${5}.fq.gz ]; then
                ${trimgalore}/trim_galore \
                    -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 6 \
                    -o ${3}/ \
                    ${2}/${1}${4}.fq.gz \
                    1>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
                mv ${3}/${1}${4}_trimmed.fq.gz ${3}/${1}${5}.fq.gz
            fi
        ;;
        "FALSE")
            if [ -z "$(find ${3} -name "${1}${5}*fq.gz" -print -quit)" ]; then
                ${trimgalore}/trim_galore \
                    --paired \
                    -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 6 \
                    -o ${3}/ \
                    ${2}/${1}${4}_1.fq.gz \
                    ${2}/${1}${4}_2.fq.gz \
                    1>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
                mv ${3}/${1}${4}_1_val_1.fq.gz ${3}/${1}${5}_1.fq.gz
                mv ${3}/${1}${4}_2_val_2.fq.gz ${3}/${1}${5}_2.fq.gz
            fi
        ;;
    esac
}