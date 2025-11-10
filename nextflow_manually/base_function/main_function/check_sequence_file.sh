#!/bin/bash
# this script is for  checking the sequence file
# author: laojp
# time: 2024.05.09
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: 
#   process_detect_pred [fq file]

# reference: https://cloud.tencent.com/developer/article/1443946
function process_detect_pred {
    echo -e "\t\t[base function] process_detect_pred"
    Pred_value=$(less $1 | head -n 1000 | awk '{if(NR%4==0) printf("%s",$0);}' | od -A n -t u1 -v \
                | awk 'BEGIN{min=100;max=0;} {for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END \
                    {if(max<=126 && min<59) print "Phred33"; \
                        else if(max>73 && min>=64) print "Phred64"; \
                        else if(min>=59 && min<64 && max>73) print "Solexa64"; \
                        else print "Unknown score encoding";}')
    declare -g "${Pred}=${Pred_value}"
    process_exist_status_background "Pred"
} 