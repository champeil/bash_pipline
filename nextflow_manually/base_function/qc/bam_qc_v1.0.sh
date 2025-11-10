#!/bin/bash
# this script is for samtools and bamqc of bam
# time: 2023.06.13
# position: SYSUCC bioinformatic platform
# author: laojp
# usage
#   process_bam_qc_check [id] [input_dir] [output_dir] [suffix_input]
#   process_bwa_bam_filter [id] [input_dir] [output_dir] [suffix_input] [suffix_output]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_bam_qc_check {
    echo -e "\t\t[base function] process_bam_qc_check"
    if [ ! -f ${3}/${1}${4}_bamqc.html ]; then
        ${bamqc}/bamqc -o ${3} ${2}/${1}${4}.bam 1>${3}/${1}${4}_bamqc_log 2>&1 
        process_exist_status_background ${1}
    fi
    if [ ! -f ${3}/${1}${4}_flagstat.txt ]; then 
        ${samtools}/samtools flagstat ${2}/${1}${4}.bam > ${3}/${1}${4}_flagstat.txt
        process_exist_status_background ${1}
    fi
}
function process_bwa_bam_filter {
    # filter 
    #  mitochodrial genome: reads mapped to MT chromosome
    #  improperly paired: -f 2: reads all mapped to chromosome properly (only for PE and remove in SE)
    #  low mapping quality: -q 20: minimum mapping quality
    #  ENCODE blacklisted regions: ${blacklist_file}

    # # Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
    # # Retain properly paired reads -f 2

    echo -e "\t\t[base function] process_bwa_bam_filter"
    if [ ! -f ${3}/${1}${5}.bam ]; then
        if [ ! -f ${3}/${1}${4}_in_blacklist.bam ]; then
            ${samtools}/samtools view -h ${2}/${1}${4}.bam | grep -v "MT" | ${samtools}/samtools view -bS -L ${blacklist_file} -U ${3}/${1}${4}_out_blacklist.bam > ${3}/${1}${4}_in_blacklist.bam \
            2>${3}/${1}${5}_filter_log 
            process_exist_status_background ${1}
        fi
        case ${single_end} in
            "FALSE")
                ${samtools}/samtools view -h -q 20 -f 2 -F 1804 ${3}/${1}${4}_out_blacklist.bam | ${samtools}/samtools sort -@ 10 -m 2G -o ${3}/${1}${5}.bam \
                2>>${3}/${1}${5}_filter_log 
                process_exist_status_background ${1}
            ;;
            "TRUE")
                ${samtools}/samtools view -h -q 20 -F 1804 ${3}/${1}${4}_out_blacklist.bam | ${samtools}/samtools sort -@ 10 -m 2G -o ${3}/${1}${5}.bam \
                2>>${3}/${1}${5}_filter_log 
                process_exist_status_background ${1}
            ;;
        esac
    fi
}
function process_bwa_bam_filter_chip {
    # filter 
    #  mitochodrial genome: reads mapped to MT chromosome
    #  ENCODE blacklisted regions: ${blacklist_file}

    echo -e "\t\t[base function] process_bwa_bam_filter_chip"
    if [ ! -f ${3}/${1}${5}.bam ]; then
        if [ ! -f ${3}/${1}${4}_in_blacklist.bam ]; then
            ${samtools}/samtools view -h ${2}/${1}${4}.bam | grep -v "MT" | ${samtools}/samtools view -bS -L ${blacklist_file} -U ${3}/${1}${4}_out_blacklist.bam > ${3}/${1}${4}_in_blacklist.bam \
            2>>${3}/${1}${5}_filter_log
            process_exist_status_background ${1}
        fi
        mv ${3}/${1}${4}_out_blacklist.bam ${3}/${1}${5}.bam
    fi
}
