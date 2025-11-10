#!/bin/bash
# this script is for reading the config.txt and read parameter
# author: laojp
# time: 2023.05.25
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_read_parameter [config] 
# ps
#   declare -g means create a global variant shared in the shell
#   declare -x = export means create a global variant in the environment variant 

# note
# belong to the base function output 

# need base function 
#   exit_state_status.sh

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_read_parameter {
    echo -e "\t\t[base function] process_read_parameter"
    while IFS= read -r line; do
        var_name=$(echo "$line" | awk -F ' = ' '{print $1}')
        var_value=$(echo "$line" | awk -F ' = ' '{print $2}')
        declare -g "${var_name}=${var_value}" # -g为创建全局变量
        process_exist_status_background ${var_name}
        #export ${var_name}="${var_value}" 
    done < <(cat ${1} | grep -vE "^$|^#")
} 