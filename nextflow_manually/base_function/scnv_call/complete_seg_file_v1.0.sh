#!/bin/bash 
# this script is for complete segment means with genome
# author: laojp
# time: 2024.07.02
# position: SYSUCC bioinformatic platform
# usage: process_complete_seg [input_id] [input_dir] [output_dir]
#   we use fasta.fai to get the length of the genome, and complete with bedtools, and segment_mean=log2(ploidy_normal(2)/ploidy_normal(2))=0

function process_complete_seg {
    echo -e "\t\t[base function] process_complete_seg"
    if [ ! -f ${reference_fasta}.fai ]; then
        ${samtools}/samtools faidx ${reference_fasta}
    fi
    if [ ! -f ${3}/$(basename ${reference_fasta}).fai ]; then
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\b' ${reference_fasta}.fai > ${3}/$(basename ${reference_fasta}).fai
    fi
    if [ ! -f ${3}/${1}.seg ]; then
        sed '1d' ${2}/${1}.seg | cut -f 1 | sort -u | \
            while read id ; do 
                ${bedtools}/bedtools complement -i <(grep ${id} ${2}/${1}.seg | \
                awk -v FS="\t" -v OFS="\t" '{print $2,$3,$4,$1,$5,$6}') -g ${3}/$(basename ${reference_fasta}).fai | \
                awk -v FS="\t" -v OFS="\t" -v name=${id} '{print name,$1,$2,$3,"0\t0"}' >> ${3}/temp.seg  
            done
        sed '1d' ${2}/${1}.seg >> ${3}/temp.seg
        head -n 1 ${2}/${1}.seg > ${3}/${1}.seg
        sort -k1 -k2 -k3,3n ${3}/temp.seg >> ${3}/${1}.seg
        rm ${3}/temp.seg
        process_exist_status_background "segment file complement"
    fi
}

