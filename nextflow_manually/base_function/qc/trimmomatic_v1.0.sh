#!/bin/bash
# this script is for trimmomatic the raw_file
# author: laojp
# time: 2023.05.24
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_trimmomatic [sample] [input_dir] [output_dir] [suffix_input] [suffix_output]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_trimmomatic {
    echo -e "\t\t[base function] process_trimmomatic"
    case ${single_end} in 
        "TRUE")
            if [ ! -f ${3}/${1}${5}.fq.gz ]; then
                java -jar ${trimmomatic}/trimmomatic-0.39.jar \
                    SE -phred33 \
                    ${2}/${1}${4}.fq.gz ${3}/${1}${5}.fq.gz \
                    ILLUMINACLIP:${adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
                    1>${3}/${1}${5}_log 2>&1 
                process_exist_status_background ${1}
            fi
            
        ;;
        "FALSE")
            if [ -z "$(find ${3} -name "${1}${5}*fq.gz" -print -quit)" ]; then
                java -jar ${trimmomatic}/trimmomatic-0.39.jar \
                    PE -phred33 \
                    ${2}/${1}${4}_1.fq.gz ${2}/${1}${4}_2.fq.gz \
                    ${3}/${1}${5}_1.fq.gz ${3}/${1}${5}_unpaired_1.fq.gz \
                    ${3}/${1}${5}_2.fq.gz ${3}/${1}${5}_unpaired_2.fq.gz \
                    ILLUMINACLIP:${adapter_file}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                    1>${3}/${1}${5}_log 2>&1
                process_exist_status_background ${1}
            fi
            
        ;;
    esac
}

