#!/bin/bash
# this script is for germline annotation
# author: laojp
# time: 2023.12.26
# position: SYSUCC bioinformatic platform
# uasge: process_vep_annotation [input vcf] [output vcf] [input dir] [output dir]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_vep_anno {
    echo -e "\t\t[base function] process_vep_anno"
    if [ ! -f ${4}/${2} ]; then
        if [[ ${vep_env} != "" ]]; then
            source ${conda_pos}/etc/profile.d/conda.sh
            conda activate ${vep_env}
        fi
        vep \
            --everything \
            --species ${species} \
            --input_file ${3}/${1} --format vcf \
            --output_file ${4}/${2} --vcf \
            --fasta ${reference_fasta} \
            --assembly ${ncbi_build} \
            --dir_cache ${vep_data} \
            --cache \
            --offline \
            --force_overwrite \
            --stats_file ${4}/$(basename ${3}/${1} .vcf)_stat.html \
            --warning_file ${4}/$(basename ${3}/${1} .vcf)_warning.html \
            --fork 10 \
            1>${4}/$(basename ${3}/${1} .vcf)_vep_log 2>&1
            process_exist_status_background $(basename ${3}/${1} .vcf)
    fi
}