#!/bin/bash
#!/bin/bash
#SBATCH --job-name="rnaseq_rsem"
#SBATCH --partition=fat-4820-Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

# this script is for rnseq [star+rsem]
# author: laojp
# time: 2023.06.12
# position: SYSUCC bioinformatic platform
# version: 1.0
# usage: rnaseq_star_rsem.sh [config.txt]

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
source ${nextflow_manually_dir}/base_function/qc/fastp_v1.0.sh
source ${nextflow_manually_dir}/base_function/map_to_reference/star_v1.0.sh
source ${nextflow_manually_dir}/base_function/quantify/rsem_v1.0.sh

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

echo -e "\tfastp......"
while read id ; do 
    process_fastp ${id} ${result_dir}/0.raw_data ${result_dir}/2.fastp "" "_fastp" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file}) 
finish_Que 

echo -e "\tfastqc_multiqc......"
process_fastqc_multiqc ${result_dir}/2.fastp ${result_dir}/3.qc

echo -e "\tSTAR......"
process_star_reference ${reference_star}
while read id ; do 
    process_star ${id} ${result_dir}/2.fastp ${result_dir}/4.star "_fastp" "_star" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})
finish_Que 

echo -e "\tRSEM......"
process_rsem_reference ${reference_rsem}
while read id ; do 
    process_rsem_star ${id} ${result_dir}/4.star ${result_dir}/5.rsem "_starAligned.toTranscriptome.out" "_count" &
    process_pid_add $!
    wait_Que $!
done < <(cat ${samplelist_file})  
finish_Que   
process_rsem_generate ${result_dir}/5.rsem ${result_dir}/5.rsem

echo -e "\toutput the script as complete_script.sh......"
process_output_complete_script ${nextflow_manually_dir} ${0} ${1} ${result_dir}

echo "end in $(date)"


