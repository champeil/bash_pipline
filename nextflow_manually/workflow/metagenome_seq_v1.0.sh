#!/bin/bash
#SBATCH --job-name="metagenome_seq"
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

# this script is for metagenome with the seq-file, no matter with rnaseq
# author: laojp
# time: 2024.05.20
# position: SYSUCC bioinformatic platform
# reference: A pan-cancer analysis of the microbiome in metastatic cancer DOI: 10.1016/j.cell.2024.03.021
# workflow
#   RNASEQ data------qc with fastp------mapped with STAR and output the unmapped reads------choose the reads both did not mapped sufficiently(samtools flag 12)
#   ------qc with fastp [phred quality score > 15, minimum length > 50 nt, trim adapter]
#   ------kraken2+pathseq metagenomic classification tool------Bracken2 to recomputes read assignments

