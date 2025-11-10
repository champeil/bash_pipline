#!/bin/bash
# this script is for homer to detect the motif
# time: 2023.08.08
# position: SYSUCC bioinformatic platform
# author: laojp
# usage
#   process_homer [id] [input_dir] [output_dir]
#   process_homer [tredt_id] [control_id] [input_dir] [output_dir]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_homer {
    echo -e "\t\t[base function] process_homer"
    if [ ! -f ${3}/${1}/knownResults.txt ]; then
        mkdir -p ${3}/${1}
        perl ${homer}/findMotifsGenome.pl ${2}/${1}/${1}_summits.bed \
            ${reference_fasta} \
            ${3}/${1} \
            1>${3}/${1}/${1}_log 2>&1 
        process_exist_status_background ${1}
    fi
}

function process_homer_chip {
    echo -e "\t\t[base function] process_homer_chip"
    if [ ! -f ${4}/${1}-${2}/knownResults.txt ]; then
        mkdir -p ${4}/${1}-${2}
        perl ${homer}/findMotifsGenome.pl ${3}/${1}-${2}/${1}_summits.bed \
            ${reference_fasta} \
            ${4}/${1}-${2}/ \
            1>${4}/${1}-${2}/${1}_log 2>&1 
        process_exist_status_background ${1}-${2}
    fi
}