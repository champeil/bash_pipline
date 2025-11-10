#!/bin/bash
#SBATCH --job-name="somatic"
#SBATCH --partition=fat-4820-Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --error=/home/laojp/data/liaok/merip/error
#SBATCH --output=/home/laojp/data/liaok/merip/log

# this script is for building the same workflow of wes_seq.sh but the storage save version
# author: laojp
# time: 2023.09.19
# position: SYSUCC bioinformatic platform
# usage: wes_somatic.sh [config]
#   the same dir_construction  of base_function
#   modify the nextflow_manually_dir value
#   modify the configuration file

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
source ${nextflow_manually_dir}/base_function/qc/trimgalore_v1.0.sh
source ${nextflow_manually_dir}/base_function/map_to_reference/bwa_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/bamstat_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/markduplicate_v1.0.sh
source ${nextflow_manually_dir}/base_function/qc/bqsr_v1.0.sh
source ${nextflow_manually_dir}/base_function/ssnv_call/mutect2_v1.0.sh
source ${nextflow_manually_dir}/base_function/annotation/somatic_annotation_v1.0.sh
source ${nextflow_manually_dir}/base_function/scnv_call/cnvkit_v1.0.sh
source ${nextflow_manually_dir}/base_function/scnv_call/gistic_v1.0.sh

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

case ${single_sample_mode} in 
    "TRUE") 
        echo -e "\t##########################################################################################################"
        echo -e "\tWarning: single sample mode start, and the result of median will be deleted after done with the last step"
        echo -e "\tWarning: the median results include: trimgalore, bwa, markduplicate and keep bqsr only. Please use it in storage saving situation"
        echo -e "\tWarning: please check the log to detect the fail-runing samples, and the error-detect-auto-quit function will update fulture"
        echo -e "\t##########################################################################################################"
        while read id; do 
            {
                if [ ! -f ${result_dir}/4.bwa/${id}_bwa.bam ]; then
                    process_trimgalore ${id} ${result_dir}/0.raw_data ${result_dir}/2.trimgalore "" "_trimgalore"
                    process_bwa ${id} ${result_dir}/2.trimgalore ${result_dir}/4.bwa "_trimgalore" "_bwa"
                    if [[ -f ${result_dir}/4.bwa/${id}_bwa.bam && $(ls -sh ${result_dir}/4.bwa/${id}_bwa.bam | cut -f 1 -d " ") != 512 && $(tail -n 1 ${result_dir}/4.bwa/${id}_log | grep -c bam_sort_core) == 1 ]]; then
                        rm -rf ${result_dir}/2.trimgalore/${id}*fq.gz
                    fi
                fi
                if [ ! -f ${result_dir}/5.markduplicate/${id}_marked.bam ]; then
                    process_markduplicate ${id} ${result_dir}/4.bwa ${result_dir}/5.markduplicate "_bwa" "_marked"
                    if [[ -f ${result_dir}/5.markduplicate/${id}_marked.bam && $(ls -sh ${result_dir}/5.markduplicate/${id}_marked.bam | cut -f 1 -d " ") != 512 && $(tail -n 1 ${result_dir}/5.markduplicate/${id}_marked_log) == 0 ]]; then
                        rm -rf ${result_dir}/4.bwa/${id}*bam
                    fi
                fi
                if [ ! -f ${result_dir}/6.bqsr/${id}_bqsr.bam ]; then
                    process_bqsr ${id} ${result_dir}/5.markduplicate ${result_dir}/6.bqsr "_marked" "_bqsr"
                    if [[ -f ${result_dir}/6.bqsr/${id}_bqsr.bam && $(ls -sh ${result_dir}/6.bqsr/${id}_bqsr.bam | cut -f 1 -d " ") != 512 && $(tail -n 1 ${result_dir}/6.bqsr/${id}_log | grep -c "Runtime.totalMemory") == 1 ]]; then
                        rm -rf ${result_dir}/5.markduplicate/${id}_marked.bam
                    fi
                fi
            } &
            process_pid_add $!
            wait_Que $!
        done < <(cat ${samplelist_file})
        finish_Que
    ;;
    "FALSE")
        echo -e "\t##########################################################################################################"
        echo -e "\tmulti-sample mode started"
        echo -e "\t##########################################################################################################"

        echo -e "\ttrimgalore......"
        while read id ; do 
            process_trimgalore ${id} ${result_dir}/0.raw_data ${result_dir}/2.trimgalore "" "_trimgalore" &
            process_pid_add $!
            wait_Que $!
        done < <(cat ${samplelist_file})
        finish_Que

        echo -e "\tfastqc_multiqc......"
        process_fastqc_multiqc ${result_dir}/2.trimgalore ${result_dir}/3.qc

        echo -e "\tbwa mapping......" 
        while read id; do 
            process_bwa ${id} \
                ${result_dir}/2.trimgalore ${result_dir}/4.bwa "_trimgalore" "_bwa" &
            process_pid_add $!
            wait_Que $!
        done < <(cat ${samplelist_file})
        finish_Que

        echo -e "\t bam stat for bwa......"
        while read id; do 
            process_bam_stat ${id} ${result_dir}/4.bwa ${result_dir}/4.bwa "_bwa" &
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
    ;;
    esac

echo -e "\tgetpileupsummaries......"
while read id; do 
    process_getpileupsummaries ${id} ${result_dir}/6.bqsr ${result_dir}/7.mutect2 "_bqsr" "_getpileupsummaries" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que

if [[ ${panel_of_normal} == "" ]]; then
    echo -e "\tcreate pon of somatic......"
    while read id; do 
        treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
        control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
        process_mutect2 ${treat_id} ${control_id} ${result_dir}/6.bqsr ${result_dir}/7.mutect2 "_bqsr" "_pon" "_f1r2" "TRUE" &
        process_pid_add $!
        wait_Que $!
    done < <(cat ${treat_control_pair})
    finish_Que
    process_somatic_pannel_of_normal ${result_dir}/7.mutect2 ${result_dir}/7.mutect2 "_pon"
fi

echo -e "\tmutect2......"
while read id; do 
    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
    process_mutect2 ${treat_id} ${control_id} ${result_dir}/6.bqsr ${result_dir}/7.mutect2 "_bqsr" "_mutect2" "_f1r2" "FALSE" &
    wait_Que $!
done < <(cat ${treat_control_pair})
finish_Que

echo -e "\tcalculatecontamination......"
while read id; do 
    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
    process_calculatecontamination ${treat_id} ${control_id} ${result_dir}/7.mutect2 ${result_dir}/7.mutect2 "_getpileupsummaries" "_contamination" "_segment" &
    wait_Que $!
done < <(cat ${treat_control_pair})
finish_Que

echo -e "\tlearnreadorientationmodel......"
while read id; do 
    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
    process_learnreadorientationmodel ${treat_id} ${result_dir}/7.mutect2 ${result_dir}/7.mutect2 "_f1r2" "_read-orientation-model" &
    wait_Que $!
done < <(cat ${treat_control_pair})
finish_Que

echo -e "\tfilter vcf......"
while read id; do 
    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
    process_filtervcf ${treat_id} ${result_dir}/7.mutect2 ${result_dir}/7.mutect2 "_mutect2" "_filter" "_contamination" "_segment" "_read-orientation-model" &
    wait_Que $!
done < <(cat ${treat_control_pair})
finish_Que

echo -e "\tmaf file creation......"
while read id; do 
    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
    process_somatic_annotation ${treat_id} ${control_id} ${result_dir}/7.mutect2 ${result_dir}/8.maf "_filter" "_anno" &
    wait_Que $!
done < <(cat ${treat_control_pair})
finish_Que
process_merge_maf ${result_dir}/8.maf ${result_dir}/8.maf "_anno"

echo -e "\tcnvkit to call......"
process_cnvkit_ref ${result_dir}/6.bqsr ${result_dir}/9.cnvkit "_bqsr" "_cnvkit"
while read id; do 
    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
    process_cnvkit_sample ${treat_id} ${control_id} ${result_dir}/6.bqsr ${result_dir}/9.cnvkit "_bqsr" "_cnvkit" &
    wait_Que $!
done < <(cat ${treat_control_pair})
finish_Que
process_cnvkit_merge ${result_dir}/9.cnvkit ${result_dir}/9.cnvkit "_cnvkit"
process_complete_seg "gistic" ${result_dir}/9.cnvkit ${result_dir}/9.cnvkit

process_gistic "gistic" ${result_dir}/9.cnvkit ${result_dir}/10.gistic

echo -e "\toutput the script as complete_script.sh......"
process_output_complete_script ${nextflow_manually_dir} ${0} ${1} ${result_dir}

echo "end in $(date)"