#!/bin/bash
#SBATCH --job-name="atac_seq"
#SBATCH --partition=fat-4820-Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

# this script is for atac_seq process: trimgalore-bwa-macs
# time: 2023.06.13
# position: SYSUCC bioinformatic platform
# author: laojp
# usage: bash atac_seq_main.sh [config]
# reference
#  DOI: 10.1186/s13059-020-1929-3: step to handle atac-seq
#  DOI: https://doi.org/10.1186/s12864-018-4559-3: to make sure the filter standard

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
source ${nextflow_manually_dir}/base_function/qc/trimmomatic_v1.0.sh
source ${nextflow_manually_dir}/base_function/map_to_reference/bwa_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/bam_qc_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/markduplicate_v1.0.sh
source ${nextflow_manually_dir}/base_function/analysis/fragment_plot_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/deeptools_v1.0.sh
source ${nextflow_manually_dir}/base_function/analysis/deeptools_v1.0.sh
source ${nextflow_manually_dir}/base_function/analysis/faCount_v1.0.sh
source ${nextflow_manually_dir}/base_function/peakcall/macs_v1.0.sh
source ${nextflow_manually_dir}/base_function/motifcall/homer_v1.0.sh
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

echo -e "\ttrimmomatic......"
while read id; do 
    process_trimmomatic ${id} \
        ${result_dir}/0.raw_data ${result_dir}/2.trimmomatic "" "_trimmomatic" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tfastqc_multiqc......"
process_fastqc_multiqc ${result_dir}/2.trimmomatic ${result_dir}/3.qc

echo -e "\tbwa mapping......"
while read id; do 
    process_bwa_atac ${id}  ${result_dir}/2.trimmomatic ${result_dir}/4.bwa "_trimmomatic" "_bwa" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tqc the bam file......"
while read id; do 
    {
        process_bam_qc_check ${id} ${result_dir}/4.bwa ${result_dir}/4.bwa "_bwa"
        process_bwa_bam_filter ${id} ${result_dir}/4.bwa ${result_dir}/5.qc "_bwa" "_bwa_filter"
        process_bam_qc_check ${id} ${result_dir}/5.qc ${result_dir}/5.qc "_bwa_filter"
        process_markduplicate ${id} ${result_dir}/5.qc ${result_dir}/5.qc "_bwa_filter" "_marked"
        process_bam_qc_check ${id} ${result_dir}/5.qc ${result_dir}/5.qc "_marked"
    } &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tbam file fragment distribution plot......"
while read id; do 
    process_fragment_distribution ${id} ${result_dir}/5.qc ${result_dir}/5.qc "_marked" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tbam file then shifting......"
while read id; do 
    process_alignmentsieve ${id} ${result_dir}/5.qc ${result_dir}/5.qc "_marked" "_marked_shift" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tshifted bam file then drawing bam coverage plot......"
while read id; do 
    process_bamcoverage ${id} ${result_dir}/5.qc ${result_dir}/5.qc "_marked_shift" "_marked_shift" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tpeak call using macs2......"
while read id; do 
    process_macs_atac ${id} ${result_dir}/5.qc ${result_dir}/6.macs "_marked_shift" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tmotif finding with homer"
while read id; do 
    process_homer ${id} ${result_dir}/6.macs ${result_dir}/7.homer &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\toutput the script as complete_script.sh......"
process_output_complete_script ${nextflow_manually_dir} ${0} ${1} ${result_dir}

echo "end in $(date)"

