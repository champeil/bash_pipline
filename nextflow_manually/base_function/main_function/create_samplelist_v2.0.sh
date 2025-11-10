#!/bin/bash
# this script is for create samplelist
# author: laojp
# time: 2023.11.20
# position: SYSUCC bioinformatic platform
# version: 2.0
# usage: process_create_sample_list [input_dir]

# updated: to support single mode and support the multiple_sample mode
#   first to judge the tumor-only 
#   then to judge the multiple id
# updated: grep sample list from treat_control_pair file
# updated: judge ifelse treat_control_pair exist, else then create according to raw_dir

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_create_sample_list {
    echo -e "\t\t[base function] process_create_sample_list"
    if [ ! -f ${samplelist_file} ]; then
        if [[ ${treat_control_pair} != "" ]]; then
            while read id; do 
                treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
                control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
                if [[ ${control_id} == "false" ]]; then
                    echo -e "\t\t\t[sample] ${treat_id} is treat_only"
                    echo ${treat_id} >> ${samplelist_file}
                else # tumor + normal paired
                    echo -e "\t\t\t[sample] ${treat_id} is treat-control paired"
                    echo ${treat_id} >> ${samplelist_file}
                    echo ${control_id} >> ${samplelist_file}
                fi
            done < <(cat ${treat_control_pair})
            sort -u ${samplelist_file} -o ${samplelist_file}
            process_exist_status_background "create sample list"
        else # create according to 0.raw_data dir
            case ${single_end} in
                "TRUE") 
                    ls ${1}/*.fq.gz | while read id ; do 
                        echo $(basename ${id} .fq.gz) >> ${samplelist_file}
                        process_exist_status_background $(basename ${id} .fq.gz)
                    done
                ;;
                "FALSE") 
                    ls ${1}/*_1.fq.gz | while read id ; do 
                        echo $(basename ${id} _1.fq.gz) >> ${samplelist_file}
                        process_exist_status_background $(basename ${id} _1.fq.gz)
                    done
                ;;
                esac
        fi
    fi
}

