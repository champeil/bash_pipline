#!/bin/bash
# this script is for mapping to the raw file by bwa
# author: laojp
# time: 2023.05.25
# position: SYSUCC bioinformatic platform
# version: 1.1
# usage: process_bwa [sample] [input_dir] [output_dir] [suffix_input] [suffix_output] 


source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_bwa {
    echo -e "\t\t[base function] process_bwa"
    case ${single_end} in
        "FALSE") 
            if [ ! -f ${3}/${1}${5}.bam ]; then
                ${bwa}/bwa mem -M -t 16 \
                    -R "@RG\tID:${1}\tSM:${1}\tLB:WXS\tPL:Illumina" \
                    ${reference_bwa} \
                    ${2}/${1}${4}_1.fq.gz \
                    ${2}/${1}${4}_2.fq.gz 2>${3}/${1}_log | ${samtools}/samtools sort -@ 10 -m 1G -o ${3}/${1}${5}.bam - \
                    1>>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
            fi
            
        ;;
        "TRUE") 
            if [ ! -f ${3}/${1}${5}.bam ]; then
                ${bwa}/bwa mem -M -t 16 \
                    -R "@RG\tID:${1}\tSM:${1}\tLB:WXS\tPL:Illumina" \
                    ${reference_bwa} \
                    ${2}/${1}${4}.fq.gz 2>${3}/${1}_log | ${samtools}/samtools sort -@ 10 -m 1G -o ${3}/${1}${5}.bam - \
                    1>>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
            fi
            
        ;;
        esac
}

function process_bwa_atac {
    echo -e "\t\t[base function] process_bwa_atac"
    case ${single_end} in
        "FALSE") 
            if [ ! -f ${3}/${1}${5}.bam ]; then
                ${bwa}/bwa mem -M -t 16 \
                    $reference_bwa \
                    ${2}/${1}${4}_1.fq.gz \
                    ${2}/${1}${4}_2.fq.gz 2>${3}/${1}_log | ${samtools}/samtools sort -@ 10 -m 1G -o ${3}/${1}${5}.bam - \
                    1>>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
            fi
            
        ;;
        "TRUE") 
            if [ ! -f ${3}/${1}${5}.bam ]; then
                ${bwa}/bwa mem -M -t 16 \
                    $reference_bwa \
                    ${2}/${1}${4}.fq.gz 2>${3}/${1}_log | ${samtools}/samtools sort -@ 10 -m 1G -o ${3}/${1}${5}.bam - \
                    1>>${3}/${1}_log 2>&1
                process_exist_status_background ${1}
            fi
            
        ;;
        esac
}

