#!/bin/bash
# this script is for star mapping 
# author: laojp
# position: SYSUCC bioinformatic platform
# time: 2023.06.12
# usage:
#   process_star [id] [input_dir] [output_dir] [input_suffix] [output_suffix] 
#   process_star_reference [output_dir] 

# update: will considerate [version] [reference] [single-paired-end] [species]

function process_star {
    echo -e "\t\t[base function] process_star"
    if [ ! -f ${3}/${1}${5}Aligned.sortedByCoord.out.bam ]; then
        case ${single_end} in
            "TRUE") 
                ${star}/STAR \
                    --runThreadN 10 \
                    --genomeDir ${reference_star} \
                    --readFilesIn ${2}/${1}${4}.fq.gz \
                    --readFilesCommand zcat \
                    --outFileNamePrefix ${3}/${1}${5} \
                    --outSAMtype BAM SortedByCoordinate \
                    --quantMode TranscriptomeSAM GeneCounts \
                    --outReadsUnmapped Fastx \
                    1>${3}/${1}${5}_log 2>&1 
                process_exist_status_background ${1}
            ;;
            "FALSE") 
                ${star}/STAR \
                    --runThreadN 10 \
                    --genomeDir ${reference_star} \
                    --readFilesIn ${2}/${1}${4}_1.fq.gz ${2}/${1}${4}_2.fq.gz \
                    --readFilesCommand zcat \
                    --outFileNamePrefix ${3}/${1}${5} \
                    --outSAMtype BAM SortedByCoordinate \
                    --quantMode TranscriptomeSAM GeneCounts \
                    --outReadsUnmapped Fastx/ \
                    1>${3}/${1}${5}_log 2>&1
                process_exist_status_background ${1}
            ;;
            esac
    fi
}

function process_star_reference {
    echo -e "\t\t[base function] process_star_reference"
    if [ ! -f ${1}/Genome ]; then
        if [ ! -d ${1} ]; then
            mkdir -p ${1}
        fi
        ${star}/STAR \
            --runMode genomeGenerate \
            --runThreadN 10 \
            --genomeDir ${1}/ \
            --genomeFastaFiles ${reference_fasta} \
            --sjdbGTFfile ${reference_gtf} \
            1>${1}/log 2>&1 
        process_exist_status_background "star_reference_building"
    fi
}