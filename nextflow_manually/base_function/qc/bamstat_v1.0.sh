# this script is for bamstat for the bam file
# author: laojp
# time: 2023.11.21
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_bam_stat [input] [input_dir] [output_dir] [input_suffix]

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_bam_stat {
    echo -e "\t\t[base function] process_bam_stat"
    if [ ! -f ${3}/${1}${4}_plot.html ]; then
        {
            ${samtools}/samtools depth ${2}/${1}${4}.bam > ${3}/${1}${4}_depth.txt
            ${samtools}/samtools stats --threads 10 ${2}/${1}${4}.bam > ${3}/${1}${4}_bamstat.txt
            ${samtools}/plot-bamstats -s ${reference_fasta} > ${3}/${1}${4}_ref_gc
            ${samtools}/plot-bamstats -r ${3}/${1}${4}_ref_gc -p ${3}/${1}${4}_plot ${3}/${1}${4}_bamstat.txt
        } 1>${3}/${1}${4}_log 2>&1 
        process_exist_status_background ${1}
    fi
}