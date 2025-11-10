#!/bin/bash
# this script is for create samplelist
# author: laojp
# time: 2023.05.31
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_create_sample_list [input_dir] [output_samplelistfile] [single | pair end]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_create_sample_list {
    echo -e "\t\t[base function] process_create_sample_list"
    if [ ! -f ${samplelist_file} ]; then
        case ${3} in
        "TRUE") 
            ls ${1}/*.fq.gz | while read id ; do 
                echo $(basename ${id} .fq.gz) >> ${samplelist_file}
                process_exist_status_background $(basename ${id} .fq.gz)
            done
        ;;
        "FALSE") 
            ls ${1}/*_1.fq.gz | while read id ; do 
                echo $(basename ${id} _1.fq.gz) >> ${samplelist_file}
                process_exist_status_background $(basename ${id} _1.fq.gz)
            done
        ;;
        esac
    fi
}