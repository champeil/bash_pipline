#!/bin/bash
# this script is for parallel the tasks and quit sequencely
# author: laojp
# time: 2023.05.24
# position: SYSUCC bioinformatic platform
# version: 1.0
# ps
#       please do not use | and use < <(cat......) instead
#               when use |, means shell1 finish and the result will be the input of the other shell2 (shell1 and shell2 are unique between each other)
#               when use < <, means all operations are construct in the same shell, means the results replace some variants of the same shells

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function PushQue {    # 将PID压入队
        Que="$Que $1"
        Nrun=$(($Nrun+1))
}
function GenQue {     # 更新队列
        OldQue=$Que
        Que=""; Nrun=0
        for PID in $OldQue; do
                if [[ -d /proc/$PID ]]; then
                        PushQue $PID
                fi
        done
}
function ChkQue {     # 检查队列
        OldQue=$Que
        for PID in $OldQue; do
                if [[ ! -d /proc/$PID ]]; then
                        GenQue; break
                fi
        done
        if [ "${OldQue}" = "${Que}" ]; then 
                old_finish_Que=1
        else
                old_finish_Que=0
        fi
}
function wait_Que { # wait_Que $!
        PID=$1
        PushQue ${PID}
        while [[ $Nrun -ge $Nproc ]]; do
                sleep 1s
                ChkQue
        done
}
function finish_Que {
        while [[ $Nrun -gt 0 ]]; do
                if [[ $old_finish_Que -eq 0 ]]; then
                        old_finish_Que=1
                fi
                sleep 1s
                ChkQue
        done
}