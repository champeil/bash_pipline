#!/bin/bash
# this script is for rsem quantify
# author: laojp
# position: SYSUCC bioinformatic platform
# time: 2023.06.13
# usage
#   process_rsem_star [id] [input_dir] [output_dir] [input_suffix] [output_suffix]
#   process_rsem_generate [input_dir] [output_dir] [input_suffix]
#   process_rsem_reference [output_dir]

function process_rsem_star {
    echo -e "\t\t[base function] process_rsem_star"
    if [ ! -f ${3}/${1}.genes.results ]; then
        case ${single_end} in
            "TRUE") 
                ${rsem}/rsem-calculate-expression \
                    --no-bam-output \
                    --alignments \
                    -p 16 \
                    -q ${2}/${1}${4}.bam \
                    ${reference_rsem} \
                    ${3}/${1}${5} \
                    1>${3}/log 2>&1
                process_exist_status_background ${1}
            ;;
            "FALSE") 
                ${rsem}/rsem-calculate-expression \
                    --paired-end \
                    --no-bam-output \
                    --alignments \
                    -p 16 \
                    -q ${2}/${1}${4}.bam \
                    ${reference_rsem} \
                    ${3}/${1}${5} \
                    1>${3}/log 2>&1
                process_exist_status_background ${1}
            ;;
            esac
        
    fi
}

function process_rsem_reference {
    echo -e "\t\t[base function] process_rsem_reference"
    if [ ! -d ${1} ]; then
        mkdir -p ${1}
    fi
    if [ ! -f ${1}/.idx.fa ]; then
        ${rsem}/rsem-prepare-reference \
            --gtf ${reference_gtf} \
            ${reference_fasta} \
            ${1}/ \
            -p 8 \
            1>${1}/log 2>&1 
        process_exist_status_background "rsem_reference_building"
    fi
}

function process_rsem_generate {
    echo -e "\t\t[base function] process_rsem_generate"
    if [ ! -f ${2}/output_matrix.txt ]; then
        ${rsem}/rsem-generate-data-matrix ${1}/*${3}*genes.results > ${2}/output_matrix.txt
        process_exist_status_background "rsem_generate"
    fi
}