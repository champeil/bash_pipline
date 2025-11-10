# this script is for chipseq workflow
# author: laojp
# time: 2023.09.07
# position: SYSUCC bioinformatic platform
# usage wangzn_chipseq.sh [0_raw_dir]

script_dir=$(realpath $0 | xargs dirname)
ls ${1}/*1.fq.gz | while read id ; do
    echo $(basename ${id} _1.fq.gz) >> ${script_dir}/samplelist
done

echo "fastqc-multiqc start"
if [ -f ${script_dir}/1.qc/multiqc_report.html ]; then
{
    /home/laojp/software/fastqc/FastQC/fastqc ${1}/*fq.gz --threads 10 -o ${script_dir}/1.qc/
    /home/laojp/software/multiqc/bin/multiqc -o ${script_dir}/1.qc/ ${script_dir}/1.qc/*.zip
} 1>${script_dir}/1.qc/qc_log 2>&1 &
if

echo "trimmomatic---bwa mem---marked---macs"
cat ${script_dir}/samplelist | while read id ; do 
    {
        echo "${id} trimmomatic"
        if [ -z "$(find ${script_dir}/2.trimmomatic/ -name "${id}*fq.gz" -print -quit)" ]; then
            java -jar ~/software/trimmomatic/Trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar \
                PE -phred33 \
                ${1}/${id}_1.fq.gz ${1}/${id}_2.fq.gz \
                ${script_dir}/2.trimmomatic/${id}_1.fq.gz ${script_dir}/2.trimmomatic/${id}_unpaired_1.fq.gz \
                ${script_dir}/2.trimmomatic/${id}_2.fq.gz ${script_dir}/2.trimmomatic/${id}_unpaired_2.fq.gz \
                ILLUMINACLIP:/home/laojp/software/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                1>${script_dir}/${id}_trimmomatic_log 2>&1
        fi
        echo "${id} bwa"
        if [ ! -f ${script_dir}/4.bwa/${id}.bam ]; then
            eval "${bwa}/bwa mem -M -t 16 \
                /home/laojp/database/mouse/ensembl/fasta/GRCm38/ensembl_102/bwa/ \
                ${script_dir}/2.trimmomatic/${id}_1.fq.gz \
                ${script_dir}/2.trimmomatic/${id}_2.fq.gz 2>${script_dir}/4.bwa/${id}_log | samtools sort -@ 10 -m 1G -o ${script_dir}/4.bwa/${id}.bam - \
                1>>${script_dir}/4.bwa/${id}_log 2>&1"
        fi
        echo "${id} markduplicate"
        if [ ! -f ${script_dir}/5.qc/${id}_marked.bam ]; then
            /home/laojp/software/gatk_4.2.6.1/gatk-4.2.6.1/gatk --java-options "-Xmx100G -Djava.io.tmpdir=${script_dir}/5.qc/${id}_temp" MarkDuplicates \
                -I ${script_dir}/4.bwa/${id}.bam \
                --REMOVE_DUPLICATES true \
                --CREATE_INDEX true \
                -O ${script_dir}/5.qc/${id}_marked.bam \
                -M ${script_dir}/5.qc/${id}_marked.metrics 
        fi
        echo "${id} macs"
        if [ ! -f ${3}/${1}_cutoff_analysis.txt ]; then 
            /home/laojp/software/macs/bin/macs2 callpeak \
                -f BAMPE \
                -g "mm" \
                --keep-dup all \
                --cutoff-analysis \
                -t ${2}/${1}${6}.bam \
                -n ${1} \
                --outdir ${3}/${1} \
                1>${3}/${1}/macs_log 2>&1 
        fi

    }
done
    


