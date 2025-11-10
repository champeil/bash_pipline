#!/bin/bash
#SBATCH --job-name="merip_seq"
#SBATCH --partition=fat-4820-Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

# this script is for merip_seq analysis
# author: laojp
# time: 2024.05.09
# position: SYSUCC bioinformatic platform
# reference: Zhang C, Tunes L, Hsieh MH, Wang P, Kumar A, Khadgi BB, Yang Y, Doxtader KA, Herrell E, Koczy O, Setlem R, Zhang X, Evers B, Wang Y, Xing C, Zhu H, Nam Y. Cancer mutations rewire the RNA methylation specificity of METTL3-METTL14. bioRxiv [Preprint]. 2023 Mar 16:2023.03.16.532618. doi: 10.1101/2023.03.16.532618. PMID: 36993753; PMCID: PMC10055151.

# first source the script
script_dir=$(realpath $0 | xargs dirname)
nextflow_manually_dir=/home/laojp/data/nextflow_manually/
source ${nextflow_manually_dir}/base_function/main_function/read_config_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/parallel_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/create_dir_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/merge_replicate_sample_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/create_samplelist_v2.0.sh
source ${nextflow_manually_dir}/base_function/main_function/output_complete_script_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/fastqc_multiqc_v1.0.sh

source ${nextflow_manually_dir}/base_function/qc/bam_qc_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/trimgalore_v1.0.sh
source ${nextflow_manually_dir}/base_function/map_to_reference/star_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/markduplicate_v1.0.sh
source ${nextflow_manually_dir}/base_function/analysis/fragment_plot_v1.0.sh
source ${nextflow_manually_dir}/base_function/analysis/deeptools_v1.0.sh
source ${nextflow_manually_dir}/base_function/analysis/faCount_v1.0.sh
process_pid_new

echo "start in $(date)"

echo -e "\tread_parameter......"
process_read_parameter ${1}

echo -e "\tcreate_dir......"
process_create_dir ${result_dir} "${create_dir_name}"

echo -e "\tmerge replicate sample"
process_merge_multiple_sample ${result_dir}/0.raw_data ${result_dir}/0.raw_data

echo -e "\tcreate_sample_list......"
process_create_sample_list ${result_dir}/0.raw_data

echo -e "\tfastqc_multiqc......"
process_fastqc_multiqc ${result_dir}/0.raw_data ${result_dir}/1.qc

echo -e "\ttringalore......"
while read id ; do 
    process_trimgalore ${id} ${result_dir}/0.raw_data ${result_dir}/2.trimgalore "" "_trimgalore" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tfastqc_multiqc......"
process_fastqc_multiqc ${result_dir}/2.trimgalore ${result_dir}/3.qc

echo -e "\tSTAR......"
process_star_reference ${reference_star}
while read id ; do 
    process_star ${id} ${result_dir}/2.trimgalore ${result_dir}/4.star "_trimgalore" "_star" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que 

echo -e "\tqc the bam file......"
while read id; do 
    {
        if [ ! -f ${result_dir}/5.qc/${id}_star.bam ]; then
            ${samtools}/samtools sort -@ 10 -m 1G ${result_dir}/4.star/${id}_starAligned.toTranscriptome.out.bam -o ${result_dir}/5.qc/${id}_star.bam
        fi
        process_markduplicate ${id} ${result_dir}/5.qc ${result_dir}/5.qc "_star" "_marked"
        ${samtools}/samtools view -h -q 256 -F 1804 ${result_dir}/5.qc/${id}_marked.bam | ${samtools}/samtools sort -@ 10 -m 2G -o ${result_dir}/5.qc/${id}_marked_filter.bam 1>${result_dir}/5.qc/${id}_marked_filter_log 2>&1
        ${samtools}/samtools flagstat ${result_dir}/5.qc/${id}_marked.bam > ${result_dir}/5.qc/${id}_marked_flagstat.txt
        ${samtools}/samtools flagstat ${result_dir}/5.qc/${id}_marked_filter.bam > ${result_dir}/5.qc/${id}_marked_filter_flags    tat.txt
    } &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tmac to detect the peaks......"
while read id; do 
    if [ ! -d ${result_dir}/6.macs/${id} ]; then
        mkdir -p ${result_dir}/6.macs/${id}
    fi
    if [ ! -f ${result_dir}/6.macs/${id}/${id}_summits.bed ]; then
        ${macs}/macs2 callpeak \
            -t /home/laojp/data/liaok/merip/5.qc/${id}_marked_filter.bam \
            -c /home/laojp/data/liaok/RNA/5.qc/${id}_marked_filter.bam \
            -g hs \
            --outdir ${result_dir}/6.macs/${id}/ \
            -n ${id} \
            -B \
            --keep-dup \
            --verbose 3 \
            -p 0.00001 \
            --nomodel --extsize 66 \
            1>${result_dir}/6.macs/${id}/log 2>&1 &
        process_pid_add $!
        wait_Que $!
    fi
done < <(cat ${samplelist_file})
finish_Que

echo "end in $(date)"