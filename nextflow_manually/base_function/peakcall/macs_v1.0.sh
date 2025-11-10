#!/bin/bash
# this script is for peak calling for some species
# time: 2023.06.19
# position: SYSUCC bioinformatic platform
# author: laojp
# usage
#  process_macs [id] [input_dir] [output_dir] [suffix_input]
#  process_macs_chip [treat_id] [control_id] [input_dir] [output_dir] [suffix_input]
# update:
#  1. will support single-end file [-f]
#  2. will support species [-g]
#  3. will support the replicate chipseq (create tre-con dir)

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_macs_atac {
    echo -e "\t\t[base function] process_macs_atac"
    if [ ! -f ${3}/${1}/${1}_cutoff_analysis.txt ]; then 
        case ${species} in
            "human") 
                local species_macs="hs"
            ;;
            "mouse") 
                local species_macs="mm"
            ;;
        esac
        if [ ! -d ${3}/${1} ]; then 
            mkdir -p ${3}/${1}
        fi
        case ${single_end} in
            "TRUE")
                ${macs}/macs2 callpeak \
                    -f BAM \
                    -g ${species_macs} \
                    --keep-dup all \
                    --cutoff-analysis \
                    -t ${2}/${1}${4}.bam \
                    -n ${1} \
                    --outdir ${3}/${1} \
                    1>${3}/${1}/macs_log 2>&1 
                process_exist_status_background ${1}
            ;;
            "FALSE")
                ${macs}/macs2 callpeak \
                    -f BAMPE \
                    -g ${species_macs} \
                    --keep-dup all \
                    --cutoff-analysis \
                    -t ${2}/${1}${4}.bam \
                    -n ${1} \
                    --outdir ${3}/${1} \
                    1>${3}/${1}/macs_log 2>&1 
                process_exist_status_background ${1}
            ;;
        esac
    fi
}

function process_macs_chip {
    echo -e "\t\t[base function] process_macs_chip"
    if [ ! -f ${4}/${1}-${2}/${1}_cutoff_analysis.txt ]; then 
        case ${species} in
            "human") 
                local species_macs="hs"
            ;;
            "mouse") 
                local species_macs="mm"
            ;;
        esac
        if [ ! -d ${4}/${1}-${2} ]; then 
            mkdir -p ${4}/${1}-${2}
        fi
        case ${single_end} in
            "TRUE")
                ${macs}/macs2 callpeak \
                    -f BAM \
                    -g ${species_macs} \
                    --keep-dup all \
                    --cutoff-analysis \
                    -t ${3}/${1}*${5}.bam \
                    -c ${3}/${2}*${5}.bam \
                    -n ${1} \
                    --outdir ${4}/${1}-${2} \
                    1>${4}/${1}-${2}/macs_log 2>&1
                process_exist_status_background "${1} compare ${2}"
            ;;
            "FALSE")
                ${macs}/macs2 callpeak \
                    -f BAMPE \
                    -g ${species_macs} \
                    --keep-dup all \
                    --cutoff-analysis \
                    -t ${3}/${1}*${5}.bam \
                    -c ${3}/${2}*${5}.bam \
                    -n ${1} \
                    --outdir ${4}/${1}-${2} \
                    1>${4}/${1}-${2}/macs_log 2>&1 
                process_exist_status_background "${1} compare ${2}"
            ;;
        esac
    fi
}
