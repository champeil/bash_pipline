#!/bin/bash
# this script is for create the dir needed
# author: laojp
# time: 2023.05.24
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_create_dir [out_folder] [create dir array]
# ps: when using [create dir array] variant, need to add quote "${create_dir_name}" or else will seen as multi-parameter but not a string

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_create_dir {
        echo -e "\t\t[base function] process_create_dir"
        IFS=' ' read -r -a folder <<< "${2}"
        for content in "${folder[@]}"; do
                if [ ! -d ${1}/${content} ]; then
                        mkdir -p ${1}/${content}
                        process_exist_status_background ${content}
                fi
        done
}