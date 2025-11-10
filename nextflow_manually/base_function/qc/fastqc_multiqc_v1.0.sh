#!/bin/bash
# this script is for fastqc to multiqc to check the quality
# author: laojp
# time: 2023.05.24
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_fastqc_multiqc [input_dir] [output_dir]
# ps: need to point out the fastqc and multiqc dir 

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_fastqc_multiqc {
        echo -e "\t\t[base function] process_fastqc_multiqc"
        if [ ! -f ${2}/multiqc_report.html ]; then
                ${fastqc}/fastqc -o ${2} --threads 20 ${1}/*.fq.gz 1>${2}/qc_log 2>&1
                process_exist_status_background "fastqc"
                ${multiqc}/multiqc -o ${2} ${2}/*.zip 1>>${2}/qc_log 2>&1
                process_exist_status_background "multiqc"
        fi
}