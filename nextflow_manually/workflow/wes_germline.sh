#!/bin/bash
#SBATCH --job-name="germline"
#SBATCH --partition=fat-4820-Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

# this sctript is for running germline variants for each samples
# author: laojp
# time: 2023.09.19
# position: SYSUCC bioinformatic platform
# usage: wes_germline.sh [config file]

# first source the script
# first source the script
script_dir=$(realpath $0 | xargs dirname)
nextflow_manually_dir=/home/chenjy/chenyx/nextflow_manually/
source ${nextflow_manually_dir}/base_function/main_function/read_config_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/parallel_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/create_dir_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/merge_replicate_sample_v1.0.sh
source ${nextflow_manually_dir}/base_function/main_function/create_samplelist_v2.0.sh
source ${nextflow_manually_dir}/base_function/main_function/output_complete_script_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/fastqc_multiqc_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/trimgalore_v1.0.sh
source ${nextflow_manually_dir}/base_function/map_to_reference/bwa_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/markduplicate_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/bqsr_v1.0.sh
source ${nextflow_manually_dir}/base_function/germline/haplotypecaller_v1.0.sh
source ${nextflow_manually_dir}/base_function/annotation/germline_annotation_v1.0.sh
process_pid_new

echo "start in $(date)"

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
process_fastqc_multiqc ${result_dir}/2.trimmomatic ${result_dir}/3.qc

echo -e "\tbwa mapping......"
while read id; do 
    process_bwa ${id} \
        ${result_dir}/2.trimmomatic ${result_dir}/4.bwa "_trimgalore" "_bwa" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tmarkduplicate to qc......"
while read id; do 
    process_markduplicate ${id} ${result_dir}/4.bwa ${result_dir}/5.markduplicate "_bwa" "_marked" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tbqsr to qc......"
while read id; do 
    process_bqsr ${id} ${result_dir}/5.markduplicate ${result_dir}/6.bqsr "_marked" "_bqsr" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\thaplotype to output the gvcf......"
while read id; do 
    process_germline_haplotypecaller_gvcf ${id} ${result_dir}/6.bqsr ${result_dir}/7.haplotypecaller "_bqsr" "_haplotype" &
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

echo -e "\tcombine gvcf......"
process_germline_combinevcf ${result_dir}/7.haplotypecaller ${result_dir}/7.haplotypecaller "_haplotype"

echo -e "\tgenotype gvcf......"
process_germline_genotypegvcf ${result_dir}/7.haplotypecaller ${result_dir}/7.haplotypecaller

echo -e "\tvqsr to qc the vcf file......"
process_germline_vqsr ${result_dir}/7.haplotypecaller ${result_dir}/7.haplotypecaller snp.germline.vcf.gz indel.germline.vcf.gz

echo -e "\tthen annotate the vcf file"
process_vep_anno snp.germline.vcf.gz snp.germline_vep.vcf ${result_dir}/7.haplotypecaller ${result_dir}/8.vep
process_vep_anno indel.germline.vcf.gz indel.germline_vep.vcf ${result_dir}/7.haplotypecaller ${result_dir}/8.vep
process_combine_germline snp.germline_vep.vcf indel.germline_vep.vcf ${result_dir}/8.vep ${result_dir}/8.vep

echo -e "\toutput the script as complete_script.sh......"
process_output_complete_script ${nextflow_manually_dir} ${0} ${1} ${result_dir}

echo "end in $(date)"
