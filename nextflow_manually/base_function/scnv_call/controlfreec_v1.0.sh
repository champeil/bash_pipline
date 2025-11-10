#!/bin/bash
# this script is for the controlfreec from bqsr
# author: laojp
# time: 2023.12.07
# position: SYSUCC bioinformatic platform
# version: 1.0
# only support the reference of human and with the chr prefix
# usage: 
#   process_controlfreec [tumor_id] [normal_id] [input_dir] [output_dir] [suffix_input] [suffix_output]
#   process_create_controlfreec_ref [reference_fasta] [output_dir]


source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_controlfreec {

}

function process_create_controlfreec_ref {
    echo -e "\t\t[base function] process_create_controlfreec_ref"
    # first is fai of the reference_fasta
    if [ ! -f ${1}.fai ]; then
        ${samtools}/samtools faidx ${1}
        process_exist_status_background "reference_fai"
    fi
    # make the gc-content with gccount
    if [ ! -f ${2}/GC_profile.cnp ]; then
        echo 





}