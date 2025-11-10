#!/bin/bash
# this script is for judge the exist status
# author: laojp
# time: 2023.05.25
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage
#   process_exist_status [content]
#   process_exist_status_background [content]
#   process_pid_new is runned in the head of a module or a workflow
#   process_pid_all [pid] : is runned after each script in the function
#   process_pid_kill : kill all pid in pid_group

# ps: about the output of the code status
#   when code block {} : only return the state of the last script of the block
#   when loop : only return the state of the last loop script
#   $! means the pid of the last script that runned in the background
#   $? means the running state of the last script
#       if a scirpt is runned in the forground, just need to script + process_exit_status
#       if a script is runned in the background, just need to script & + process_exit_status + add_pid

function process_exist_status {
    local status_id=$?

    if [ $status_id -eq 0 ]; then
        echo -e "\t\t\t[sample] ${1:-""}\tsuccessed"  # 输出成功的标记
    else
        echo -e "\t\t\t[sample] ${1:-""}\tfailed with status $status_id"  # 输出失败的标记
        exit $status_id
    fi
}

function process_exist_status_background {
    local status_id=$?

    if [ $status_id -eq 0 ]; then
        echo -e "\t\t\t[sample] ${1:-""}\tsuccessed"  # 输出成功的标记
    else
        echo -e "\t\t\t[sample] ${1:-""}\tfailed with status $status_id"  # 输出失败的标记
        if [[ -v pid_group ]]; then
            process_pid_kill
        fi
        exit $status_id
    fi
}




function process_pid_new {
    if [[ ! -v pid_group ]]; then
        declare -ag pid_group=()
    fi
}

function process_pid_renew {
    declare -ag pid_group=()
}

function process_pid_add {
    if [[ ! -v pid_group ]]; then
        declare -ag pid_group=()
    fi

    pid_group+=($1)
}

function process_pid_kill {
    for pid_mem in "${pid_group[@]}"; do
        if kill -O $pid_mem 2>/dev/null; then
            kill $pid_mem
        fi
    done
}