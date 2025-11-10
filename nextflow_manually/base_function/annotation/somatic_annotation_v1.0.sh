#!/bin/bash
# this script is for wes_process: gatk[getpileupsummary-mutect2]
# time: 2023.06.20
# position: SYSUCC bioinformatic platform
# author: laojp
# usage: 
#  process_somatic_annotation [tumor_id] [normal_id] [input_dir] [output_dir] [suffix_input] [suffix_output]
#  process_merge_maf [input_dir] [output_dir] [suffix_input]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_somatic_annotation {
    echo -e "\t\t[base function] process_somatic_annotation"
    if [[ ${2} == "false" ]]; then # tumor only mode
        if [ ! -f ${4}/${1}${6}.maf ]; then 
            if [[ ${vep_env} != "" ]]; then
                source ${conda_pos}/etc/profile.d/conda.sh
                conda activate ${vep_env}
            fi
            ${conda_pos}/envs/${vep_env}/bin/perl ${vcf2maf}/vcf2maf.pl \
                --input-vcf ${3}/${1}${5}.vcf \
                --output-maf ${4}/${1}${6}.maf \
                --ref-fasta ${reference_fasta} \
                --tumor-id ${1} \
                --vep-path ${vep_path} \
                --vep-data ${vep_data} \ # must be {vep_data}/{species}/{vep_version}_{ncbi_build} format
                --species ${species} \
                --ncbi-build ${ncbi_build} \
                --vep-overwrite true \
                1>${4}/${1}${6}_log 2>&1
            process_exist_status_background ${1}
        fi
    else 
        if [ ! -f ${4}/${1}${6}.maf ]; then 
            if [[ ${vep_env} != "" ]]; then
                source ${conda_pos}/etc/profile.d/conda.sh
                conda activate ${vep_env}
            fi
            ${conda_pos}/envs/${vep_env}/bin/perl ${vcf2maf}/vcf2maf.pl \
                --input-vcf ${3}/${1}${5}.vcf \
                --output-maf ${4}/${1}${6}.maf \
                --ref-fasta ${reference_fasta} \
                --tumor-id ${1} \
                --normal-id ${2} \
                --vep-path ${vep_path} \
                --vep-data ${vep_data} \
                --species ${species} \
                --ncbi-build ${ncbi_build} \
                --vep-overwrite true \
                1>${4}/${1}${6}_log 2>&1 
            process_exist_status_background ${1}
        fi
    fi
}

function process_merge_maf {
    if [ ! -f ${2}/vep_merge.maf ]; then
        cat ${1}/*${3}.maf | grep -v "^#" | grep -v "^Hugo_Symbol" > ${2}/vep_merge_temp.maf
        ls ${2}/*${3}.maf | grep -v "vep_merge_temp" | head -n 1 | while read id ; do 
            grep "^Hugo_Symbol" ${id} > ${2}/vep_merge.maf 
        done
        cat ${2}/vep_merge_temp.maf >> ${2}/vep_merge.maf
        rm ${2}/vep_merge_temp.maf
    fi
}