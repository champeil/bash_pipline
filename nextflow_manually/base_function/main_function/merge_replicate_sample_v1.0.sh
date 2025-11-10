# this script is for merging the fq.gz sample with the multiple samples
# author: laojp
# time: 2023.11.20
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: process_merge_multiple_sample [input_dir] [output_dir] 
# sort the sample, use the first name_mer as new and combine the files, to avoid the replicate of the same merging samples

function process_merge_multiple_sample {
    echo -e "\t\t[base function] process_merge_multiple_sample"
    if [[ ${treat_control_pair} != "" ]]; then
        case ${single_end} in
            "FALSE")
                while read id; do
                    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
                    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
                    if [[ ${treat_id} == *";"* ]]; then 
                        treat_id_merge=$(echo $(echo ${treat_id} | awk -v FS=";" '{for(i=1;i<=NF;i++){print $i}}' | sort | head -n 1)_mer)
                        if [ $(ls ${2}/${treat_id_merge}*.fq.gz 2>/dev/null | wc -l) == "0" ]; then
                            echo -e "\t\t\t${treat_id} need to merge replicate in treat: ${treat_id}"
                            echo ${treat_id} | awk -v FS=";" -v OFS="\t" '{for(i=1;i<=NF;i++){print $i}}' | while read id; do
                                cat "${1}/${id}_1.fq.gz" >> ${2}/${treat_id_merge}_1.fq.gz
                                cat "${1}/${id}_2.fq.gz" >> ${2}/${treat_id_merge}_2.fq.gz
                            done
                            sed -i 's/'"${treat_id}"'/'"${treat_id_merge}"'/g' ${treat_control_pair}
                            process_exist_status_background "${treat_id_merge} merge"
                        fi
                    elif [[ ${control_id} == *";"* ]]; then 
                        control_id_merge=$(echo $(echo ${control_id} | awk -v FS=";" '{for(i=1;i<=NF;i++){print $i}}' | sort | head -n 1)_mer)
                        if [ $(ls ${3}/${control_id_merge}*.fq.gz 2>/dev/null | wc -l) == "0" ]; then
                            echo -e "\t\t\t${treat_id} need to merge replicate in control: ${control_id}"
                            echo ${control_id} | awk -v FS=";" -v OFS="\t" '{for(i=1;i<=NF;i++){print $i}}' | while read id; do
                                cat "${1}/${id}_1.fq.gz" >> ${2}/${control_id_merge}_1.fq.gz
                                cat "${1}/${id}_2.fq.gz" >> ${2}/${control_id_merge}_2.fq.gz
                            done
                            sed -i 's/'"${control_id}"'/'"${control_id_merge}"'/g' ${treat_control_pair}
                            process_exist_status_background "${control_id_merge} merge"
                        fi
                    else 
                        echo -e "\t\t\t${treat_id} need not to merge replicate"
                    fi
                done < <(cat ${treat_control_pair})
            ;;
            "TRUE")
                while read id; do 
                    treat_id=$(echo ${id} | awk -v FS=" " '{print $1}')
                    control_id=$(echo ${id} | awk -v FS=" " '{print $2}')
                    if [[ ${treat_id} == *";"* ]]; then 
                        treat_id_merge=$(echo $(echo ${treat_id} | awk -v FS=";" '{for(i=1;i<=NF;i++){print $i}}' | sort | head -n 1)_mer)
                        if [ $(ls ${2}/${treat_id_merge}*.fq.gz 2>/dev/null | wc -l) == "0" ]; then
                            echo -e "\t\t\t${treat_id} need to merge replicate in treat: ${treat_id}"
                            echo ${treat_id} | aw
                            k -v FS=";" -v OFS="\t" '{for(i=1;i<=NF;i++){print $i}}' | while read id; do
                                cat "${1}/${id}.fq.gz" >> ${2}/${treat_id_merge}.fq.gz
                            done
                            sed -i 's/'"${treat_id}"'/'"${treat_id_merge}"'/g' ${treat_control_pair}
                            process_exist_status_background "${treat_id_merge} merge"
                        fi
                    elif [[ ${control_id} == *";"* ]]; then 
                        control_id_merge=$(echo $(echo ${control_id} | awk -v FS=";" '{for(i=1;i<=NF;i++){print $i}}' | sort | head -n 1)_mer)
                        if [ $(ls ${2}/${control_id_merge}*.fq.gz 2>/dev/null | wc -l) == "0" ]; then
                            echo -e "\t\t\t${treat_id} need to merge replicate in control: ${control_id}"
                            echo ${control_id} | awk -v FS=";" -v OFS="\t" '{for(i=1;i<=NF;i++){print $i}}' | while read id; do
                                cat "${1}/${id}.fq.gz" >> ${2}/${control_id_merge}.fq.gz
                            done
                            sed -i 's/'"${control_id}"'/'"${control_id_merge}"'/g' ${treat_control_pair}
                            process_exist_status_background "${control_id_merge} merge"
                        fi
                    else 
                        echo -e "\t\t\t${treat_id} need not to merge replicate"
                    fi
                done < <(cat ${treat_control_pair})
            ;;
        esac
    fi
}