#!/bin/bash
# this script is for quantify the count with featurecount
# author: laojp
# time: 2023.12.29
# position: SYSUCC bioinformatic platform
# usage
#   process_featurecount [id] [input_dir] [output_dir] [input_suffix] [output_suffix]
#   process_featurecount_merge [input_dir] [output_dir] [input_suffix] [output_suffix]
# fulture: will support each samples separately and then merge the specific samples together

source ${nextflow_manually_dir}/base_function/main_function/exit_state_status_v1.0.sh

function process_featurecount {
    echo -e "\t\t[base function] process_featurecount"
    if [ ! -f ${3}/${1}${5}.txt ]; then
        case ${single_end} in
            "TRUE") 
                ${featurecount}/featureCounts \
                    -t exon \
                    -g gene_id \
                    -a ${reference_gtf} \
                    -o ${3}/${1}${5}.txt \
                    ${2}/${1}${4}.bam \
                    1>${3}/${1}${5}_log 2>&1 
                process_exist_status_background ${1}
            ;;
            "FALSE") 
                ${featurecount}/featureCounts \
                    -p \
                    --countReadPairs \
                    -t exon \
                    -g gene_id \
                    -a ${reference_gtf} \
                    -o ${3}/${1}${5}.txt \
                    ${2}/${1}${4}.bam \
                    1>${3}/${1}${5}_log 2>&1 
                process_exist_status_background ${1}
            ;;
            esac
    fi
}

function process_featurecount_merge {
    echo -e "\t\t[base function] process_featurecount_merge"
    if [ ! -f ${2}/count${4}.txt ]; then
        ls ${1}/*${3}.txt | head -n 1 | xargs cat | grep -v "#" | awk -v FS="\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${2}/count${4}.txt
        ls ${1}/*${3}.txt | while read id ; do 
            name=$(basename ${id} ${3}.txt)
            paste ${2}/count${4}.txt <(grep -v "#" ${id} | awk -v name=${name} -v FS="\t" '{if(NR==1){$7=name};print $7}') > temp ;
            mv temp ${2}/count${4}.txt
        done
        process_exist_status_background "feturecount merging"
    fi
}