#!/bin/bash
# this script is for cnv detecting in somatic
# time: 2023.06.16
# position: SYSUCC bioinformatic platform
# author: laojp
# description
#   we first prepare all necessary files
#       access_ref.bed: file for limit the background region
#       target.bed: then autobin all input bam file and decide the bin size automatically, then output target and anti-target based on background region
#       access_ref.cnn: this is for the tumor-only file, which would infer the cn with the assumtion that the netural is 0.5
#   we then handle each tumor-paired or tumor-only file respectively
#       if got normal-paired, then use normal as reference, scaled based on background and normal, then output the segment result
#       if tumor-only, then use background as reference, only scale for background and then output the segment result
# update
#   will support amption and WGS fulture
#   will support refflat file creating manually
#   will support vcf input and allelic specific cn infering
#   will support tumor-only with many normal, which will use pool normal samples stategy, but not suggest

# usage: 
#   process_cnvkit_ref [input_dir] [output_dir] [suffix_input] [suffix_output]
#   process_cnvkit_sample [tumor_id] [normal_id] [input_dir] [output_dir] [suffix_input] [suffix_output]
#   process_cnvkit_merge [input_dir] [output_dir] [input_suffix]
#   process_complete_seg [input_id] [input_dir] [output_dir]

function process_cnvkit_ref {
    echo -e "\t\t[base function] process_cnvkit_ref"
    if [[ ${cnvkit_env} != "" ]]; then
        source ${conda_pos}/etc/profile.d/conda.sh
        conda activate ${cnvkit_env}
    fi
    if [ ! -f ${2}/access_ref.bed ]; then
        ${cnvkit}/cnvkit.py access ${reference_fasta} -o ${2}/access_ref.bed \
            1>>${2}/ref_building_log 2>&1 
        process_exist_status_background "cnvkit ref_access"
    fi
    if [[ ${interval_list} != "" ]]; then # WES mode
        if [ ! -f ${2}/target.bed ]; then
            ${cnvkit}/cnvkit.py autobin ${1}/*${3}.bam \
                -m hybrid \
                -t ${interval_list} \
                -g ${2}/access_ref.bed \
                --annotate ${refflat} \
                --target-output-bed ${2}/target.bed \
                --antitarget-output-bed ${2}/antitarget.bed \
                1>>${2}/ref_building_log 2>&1
            process_exist_status_background "cnvkit ref_autobin"
        fi
        if [ ! -f ${2}/access_ref.cnn ]; then
            ${cnvkit}/cnvkit.py reference \
                -o ${2}/access_ref.cnn \
                -f ${reference_fasta} \
                -t ${2}/target.bed \
                -a ${2}/antitarget.bed \
                1>>${2}/ref_building_log 2>&1
            process_exist_status_background "cnvkit ref_reference"
        fi
    else
        if [ ! -f ${2}/target.bed ]; then # WGS mode
            ${cnvkit}/cnvkit.py autobin ${1}/*${3}.bam \
                -m wgs \
                -b 50000 \
                -g ${2}/access_ref.bed \
                --annotate ${refflat} \
                --target-output-bed ${2}/target.bed \
                --antitarget-output-bed ${2}/antitarget.bed \
                1>>${2}/ref_building_log 2>&1
            process_exist_status_background "cnvkit ref_autobin"
        fi
        if [ ! -f ${2}/access_ref.cnn ]; then
            ${cnvkit}/cnvkit.py reference \
                -o ${2}/access_ref.cnn \
                -f ${reference_fasta} \
                -t ${2}/target.bed \
                -a ${2}/antitarget.bed \
                1>>${2}/ref_building_log 2>&1
            process_exist_status_background "cnvkit ref_reference"
        fi
    fi
}


function process_cnvkit_sample {
    echo -e "\t\t[base function] process_cnvkit_sample"
    if [[ ${cnvkit_env} != "" ]]; then
        source ${conda_pos}/etc/profile.d/conda.sh
        conda activate ${cnvkit_env}
    fi
    if [ ! -f ${4}/${1}${6}_coverage.cnn ]; then
        ${cnvkit}/cnvkit.py coverage \
            ${3}/${1}${5}.bam \
            ${4}/antitarget.bed \
            -o ${4}/${1}${6}_anti_coverage.cnn \
            1>>${4}/${1}${6}_log 2>&1
        ${cnvkit}/cnvkit.py coverage \
            ${3}/${1}${5}.bam \
            ${4}/target.bed \
            -o ${4}/${1}${6}_coverage.cnn \
            1>>${4}/${1}${6}_log 2>&1 
        process_exist_status_background "${1} coverage"
    fi
    if [[ ${2} != "false" ]]; then
        if [ ! -f ${4}/${2}${6}_coverage.cnn ]; then
            ${cnvkit}/cnvkit.py coverage \
                ${3}/${2}${5}.bam \
                ${4}/antitarget.bed \
                -o ${4}/${2}${6}_anti_coverage.cnn \
                1>>${4}/${1}${6}_log 2>&1
            ${cnvkit}/cnvkit.py coverage \
                ${3}/${2}${5}.bam \
                ${4}/target.bed \
                -o ${4}/${2}${6}_coverage.cnn \
                1>>${4}/${1}${6}_log 2>&1
            process_exist_status_background "${2} coverage"
        fi
        if [ ! -f ${4}/${2}${6}_ref.cnn ]; then
            ${cnvkit}/cnvkit.py reference \
                -t ${4}/${2}${6}_coverage.cnn \
                -a ${4}/${2}${6}_anti_coverage.cnn \
                --fasta ${reference_fasta} \
                -o ${4}/${2}${6}_ref.cnn \
                1>>${4}/${1}${6}_log 2>&1
            process_exist_status_background "${2} coverage"
        fi
        if [ ! -f ${4}/${1}${6}_fix.cnr ]; then
            ${cnvkit}/cnvkit.py fix \
                ${4}/${1}${6}_coverage.cnn \
                ${4}/${1}${6}_anti_coverage.cnn \
                ${4}/${2}${6}_ref.cnn \
                -i ${1} \
                -o ${4}/${1}${6}_fix.cnr \
                1>>${4}/${1}${6}_log 2>&1
            process_exist_status_background "${1} fix"
        fi
    else
        if [ ! -f ${4}/${1}${6}_fix.cnr ]; then
            ${cnvkit}/cnvkit.py fix \
                ${4}/${1}${6}_coverage.cnn \
                ${4}/${1}${6}_anti_coverage.cnn \
                ${4}/access_ref.cnn \
                -i ${1} \
                -o ${4}/${1}${6}_fix.cnr \
                1>>${4}/${1}${6}_log 2>&1
            process_exist_status_background "${1} fix"
        fi
    fi
    if [ ! -f ${4}/${1}${6}_seg.cns ]; then
        ${cnvkit}/cnvkit.py segment \
            ${4}/${1}${6}_fix.cnr \
            -m cbs -t 0.05 \
            --drop-low-coverage \
            --drop-outliers 10 \
            -o ${4}/${1}${6}_seg.cns \
            1>>${4}/${1}${6}_log 2>&1
        process_exist_status_background "${1} segment"
    fi
    if [ ! -f ${4}/${1}${6}_scatter.pdf ]; then
        ${cnvkit}/cnvkit.py scatter \
            ${4}/${1}${6}_fix.cnr \
            -s ${4}/${1}${6}_seg.cns \
            -o ${4}/${1}${6}_scatter.pdf
            1>>${4}/${1}${6}_log 2>&1
        process_exist_status_background "${1} scatter"
    fi
    if [ ! -f ${4}/${1}${6}_diagram.pdf ]; then
        ${cnvkit}/cnvkit.py diagram \
            ${4}/${1}${6}_fix.cnr \
            --no-shift-xy \
            -s ${4}/${1}${6}_seg.cns \
            -o ${4}/${1}${6}_diagram.pdf
            1>>${4}/${1}${6}_log 2>&1
        process_exist_status_background "${1} diagram"
    fi
}

function process_cnvkit_merge {
    echo -e "\t\t[base function] process_cnvkit_merge"
    if [ ! -f ${2}/heatmap.pdf ]; then
        ${cnvkit}/cnvkit.py heatmap \
            ${1}/*_seg.cns -d \
            --no-shift-xy \
            -o ${2}/heatmap.pdf
            1>${2}/heatmap_log 2>&1
        process_exist_status_background "cnvkit heatmap"
    fi
    if [ ! -f ${2}/gistic.seg ]; then
        ${cnvkit}/cnvkit.py export seg \
            ${1}/*_seg.cns \
            -o ${2}/gistic.seg
            1>${2}/seg_output_log 2>&1
        sed -i 's/_seg//g' ${2}/gistic.seg
        process_exist_status_background "cnvkit seg_output"
    fi
}

function process_complement_seg {
    echo -e "\t\t[base function] process_complete_seg"
    if [ ! -f ${reference_fasta}.fai ]; then
        ${samtools}/samtools faidx ${reference_fasta}
    fi
    if [ ! -f ${3}/$(basename ${reference_fasta}).fai ]; then
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\b' ${reference_fasta}.fai > ${3}/$(basename ${reference_fasta}).fai
    fi
    if [ ! -f ${3}/${1}_complement.seg ]; then
        sed '1d' ${2}/${1}.seg | cut -f 1 | sort -u | \
            while read id ; do 
                ${bedtools}/bedtools complement -i <(grep ${id} ${2}/${1}.seg | \
                awk -v FS="\t" -v OFS="\t" '{print $2,$3,$4,$1,$5,$6}') -g ${3}/$(basename ${reference_fasta}).fai | \
                awk -v FS="\t" -v OFS="\t" -v name=${id} '{print name,$1,$2,$3,"0\t0"}' >> ${3}/temp.seg  
            done
        sed '1d' ${2}/${1}.seg >> ${3}/temp.seg
        head -n 1 ${2}/${1}.seg > ${3}/${1}_complement.seg
        sort -k1 -k2 -k3,3n ${3}/temp.seg >> ${3}/${1}_complement.seg
        rm ${3}/temp.seg
        process_exist_status_background "segment file complement"
    fi
}

