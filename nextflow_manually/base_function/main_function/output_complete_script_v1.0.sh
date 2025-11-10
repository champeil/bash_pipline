# this script is for output the complete script for workflow
# author: laojp
# time: 2023.10.27
# position: SYSUCC bioinformatic platform
# usage: output_complete_script [nextflow_dir] [workflow_script] [configure_script] [output_dir]

function process_output_complete_script {
    if [ ! -f ${4}/complete_script.sh ]; then
        echo -e "\t\toutput the annotation......"
        echo -e "# this script is complete script outputed by [process_output_complete_script] base function" >> ${4}/complete_script.sh
        echo -e "# author: laojp" >> ${4}/complete_script.sh
        echo -e "# time: $(date)" >> ${4}/complete_script.sh
        echo -e "# position: SYSUCC bioinformatic platform" >> ${4}/complete_script.sh
        echo -e "# usage: bash complete_script.sh 1>log 2>&1" >> ${4}/complete_script.sh

        echo -e "\t\toutput the configure part......"
        echo -e "\n\n# first is the configure part" >> ${4}/complete_script.sh
        while IFS= read -r line; do
            var_name=$(echo "$line" | awk -F ' = ' '{print $1}')
            var_value=$(echo "$line" | awk -F ' = ' '{print $2}')
            echo -e  "${var_name}=${var_value}" >> ${4}/complete_script.sh
        done < <(cat ${3} | grep -vE "^$|^#")

        echo -e "\t\toutput the base function......"
        echo -e "\n\n# second is base function to use" >> ${4}/complete_script.sh
        grep "^\s*process" ${2} | sed 's/^[[:space:]]*//' | awk -v FS=" " '{print $1}' | sort | uniq | while read process_function ; do
            process_pos=$(grep -rlw ${process_function} ${1}/base_function/* | xargs realpath)
            sed -n '/^function '"${process_function}"'/,/^function process_/p' ${process_pos} | head -n -1 >> ${4}/complete_script.sh
        done

        echo -e "\t\tcreate pid new......"
        echo -e "\n\n# third create pid new" >> ${4}/complete_script.sh
        echo -e "process_pid_new" >> ${4}/complete_script.sh

        echo -e "\t\toutput main script......"
        echo -e "\n\n# fourth output main part of script" >> ${4}/complete_script.sh
        sed -n '/^echo "start in $(date)"/,/^echo "end in $(date)"/p' ${2} >> ${4}/complete_script.sh

        echo -e "\t\tdone......"
    fi
}
