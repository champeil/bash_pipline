#!/bin/bash
# this script is for wes_process: gatk[getpileupsummary-mutect2]
# time: 2023.06.16
# position: SYSUCC bioinformatic platform
# author: laojp
# usage: 
#  process_getpileupsummaries [id] [input_dir] [output_dir] [suffix_input] [suffix_output]
#  process_mutect2 [tumor_id] [normal_id] [inputdir] [output_dir] [suffix_input] [suffix_output] [suffix_f1r2] [pon_mode]
#   new: to support the tumor only mode
#  process_calculatecontamination [tumor_id] [normal_id] [inputdir] [output_dir] [suffix_input] [suffix_output] [suffix_segment_table]
#  process_learnreadorientationmodel [tumor_id] [inputdir] [outputdir] [suffix_f1r2] [suffix_output]
#  process_filtervcf [tumor_id] [inputdir] [outputdir] [suffix_input] [suffix_output] [suffix_contamination] [suffix_tumor_segment] [suffix_orienmodel]
#  process_pon_create [input_dir] [output_dir] [suffix_input]
#  because gatk genomicinportdb need a lot of device space, so we use the lower version of gatk to call pon and the latest version is 4.1.0.0

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_getpileupsummaries {
    echo -e "\t\t[base function] process_getpileupsummaries"
    if [ ! -f ${3}/${1}${5}.table ]; then
        ${gatk}/gatk GetPileupSummaries \
            -I ${2}/${1}${4}.bam \
            -V ${common_variant_af} \
            -L ${interval_list} \
            -O ${3}/${1}${5}.table \
            1>${3}/${1}${5}_log 2>&1
        process_exist_status_background ${1}
    fi
}

function process_mutect2 {
    echo -e "\t\t[base function] process_mutect2"
    case ${8} in
        "TRUE") 
            if [ ! -f ${4}/${2}${6}.vcf ]; then  # pon mode
                ${gatk}/gatk Mutect2 \
                -R ${reference_fasta} \
                -I ${3}/${2}${5}.bam \
                -tumor ${2} \
                -L ${interval_list} \
                --germline-resource ${af_only_gnomad} \
                -O ${4}/${2}${6}.vcf \
                1>${4}/${2}${6}_log 2>&1 
                process_exist_status_background ${1}
            fi
        ;;
        "FALSE")
            if [ ! -f ${4}/${1}${6}.vcf ]; then
                if [[ ${interval_list} != "" ]]; then # WES mode
                    if [[ ${2} == "false" ]]; then # tumor_only mode
                        ${gatk}/gatk Mutect2 \
                            -R ${reference_fasta} \
                            -I ${3}/${1}${5}.bam \
                            -tumor ${1} \
                            -L ${interval_list} \
                            --panel-of-normals ${panel_of_normal} \
                            --germline-resource ${af_only_gnomad} \
                            --f1r2-tar-gz ${4}/${1}${7}.tar.gz \
                            --genotype-germline-sites TRUE \
                            --genotype-pon-sites TRUE \
                            -O ${4}/${1}${6}.vcf \
                            1>${4}/${1}${6}_log 2>&1 
                            process_exist_status_background ${1} 
                    else # tumor-normal paired mode
                        ${gatk}/gatk Mutect2 \
                            -R ${reference_fasta} \
                            -I ${3}/${1}${5}.bam \
                            -I ${3}/${2}${5}.bam \
                            -normal ${2} \
                            -tumor ${1} \
                            -L ${interval_list} \
                            --panel-of-normals ${panel_of_normal} \
                            --germline-resource ${af_only_gnomad} \
                            --f1r2-tar-gz ${4}/${1}${7}.tar.gz \
                            --genotype-germline-sites TRUE \
                            --genotype-pon-sites TRUE \
                            -O ${4}/${1}${6}.vcf \
                            1>${4}/${1}${6}_log 2>&1 
                            process_exist_status_background ${1}
                    fi
                else # WGS mode
                    if [[ ${2} == "false" ]]; then # tumor_only mode
                        ${gatk}/gatk Mutect2 \
                            -R ${reference_fasta} \
                            -I ${3}/${1}${5}.bam \
                            -tumor ${1} \
                            --panel-of-normals ${panel_of_normal} \
                            --germline-resource ${af_only_gnomad} \
                            --f1r2-tar-gz ${4}/${1}${7}.tar.gz \
                            --genotype-germline-sites TRUE \
                            --genotype-pon-sites TRUE \
                            -O ${4}/${1}${6}.vcf \
                            1>${4}/${1}${6}_log 2>&1 
                            process_exist_status_background ${1} 
                    else # tumor-normal paired mode
                        ${gatk}/gatk Mutect2 \
                            -R ${reference_fasta} \
                            -I ${3}/${1}${5}.bam \
                            -I ${3}/${2}${5}.bam \
                            -normal ${2} \
                            -tumor ${1} \
                            --panel-of-normals ${panel_of_normal} \
                            --germline-resource ${af_only_gnomad} \
                            --f1r2-tar-gz ${4}/${1}${7}.tar.gz \
                            --genotype-germline-sites TRUE \
                            --genotype-pon-sites TRUE \
                            -O ${4}/${1}${6}.vcf \
                            1>${4}/${1}${6}_log 2>&1 
                            process_exist_status_background ${1}
                    fi
                fi
            fi
        ;;
    esac
}

function process_calculatecontamination {
    echo -e "\t\t[base function] process_calculatecontamination"
    if [ ! -f ${4}/${1}${6}.table ]; then
        if [[ ${2} != "false" ]]; then
            ${gatk}/gatk CalculateContamination \
                -I ${4}/${1}${5}.table \
                -matched ${4}/${2}${5}.table \
                -O ${4}/${1}${6}.table \
                --tumor-segmentation ${4}/${1}${7}.table \
                1>${4}/${1}${6}_log 2>&1
        else
            ${gatk}/gatk CalculateContamination \
                -I ${4}/${1}${5}.table \
                -O ${4}/${1}${6}.table \
                --tumor-segmentation ${4}/${1}${7}.table \
                1>${4}/${1}${6}_log 2>&1
        fi
        process_exist_status_background ${1}
    fi
}

function process_learnreadorientationmodel {
    echo -e "\t\t[base function] process_learnreadorientationmodel"
    if [ ! -f ${3}/${1}${5}.tar.gz ]; then
        ${gatk}/gatk LearnReadOrientationModel \
            -I ${2}/${1}${4}.tar.gz \
            -O ${3}/${1}${5}.tar.gz \
            1>${3}/${1}${5}_log 2>&1
        process_exist_status_background ${1}
    fi
}

function process_filtervcf {
    echo -e "\t\t[base function] process_filtervcf"
    if [ ! -f ${3}/${1}${5}.vcf ]; then
        ${gatk}/gatk FilterMutectCalls \
            -R ${reference_fasta} \
            -V ${2}/${1}${4}.vcf \
            -O ${3}/${1}_somatic.vcf \
            --contamination-table ${2}/${1}${6}.table \
            --tumor-segmentation ${2}/${1}${7}.table \
            --ob-priors ${2}/${1}${8}.tar.gz \
            1>${3}/${1}${5}_log 2>&1 
        process_exist_status_background ${1}
        cat ${3}/${1}_somatic.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0]=~/_/; print}}' > ${3}/${1}${5}.vcf
    fi
}


function process_somatic_pannel_of_normal {
    echo -e "\t\t[base function] process_somatic_pannel_of_normal"
    if [ ! -f ${2}/pon_list.args ]; then
        awk -v FS="\t" '{print $2}' ${treat_control_pair} | while read id ; do
            find ${1} -name "*${id}*${3}*vcf" | xargs realpath >> ${2}/pon_list.args
            process_exist_status_background ${id}
        done
    fi
    if [ ! -f ${2}/panel_of_normal.vcf ]; then
        ${gatk_pon}/gatk CreateSomaticPanelOfNormals \
            -vcfs ${2}/pon_list.args \
            -O ${2}/panel_of_normal.vcf \
            1>${2}/CreateSomaticPanelOfNormals_log 2>&1
            process_exist_status_background "CreateSomaticPanelOfNormals"
    fi
    declare -g "panel_of_normal=$(realpath ${2}/panel_of_normal.vcf)"
    process_exist_status_background "declare global pon"
}
