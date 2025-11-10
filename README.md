# bash_pipeline
bash_pipeline使用bash构建DNA_seq、RNA_seq、cuttag_seq(chip-seq, atac-seq)流程，支持步骤样本并行与sbatch提交。
与nextflow pipeline相比，缺少断点续分析、按样本独立并行与多sbatch节点分流等……运行速度远慢与nextflow，但是支持代码直接移植、容易上手、支持定制化、并支持单双端测序、tumor-only样本、t-n pair样本的运行，所以具有一定的使用价值

# usage
## running
- 下载文件夹
- 修改workflow主流程代码中，nextflow_manually_dir的路径为该文件夹
- 在结果目录下构建0.raw_data，并将原始文件放里面，修改文件结尾为\[sample_name\]_\[1\|2\].fq.gz格式
- 如果是DNA-seq测序，则需要tumor_normal_pair信息，没有表头，制表符分割，第一列为tumor，第二列为normal的sample_name
- 修改config内容
- sbatch workflow.sh
  
## error
- 该代码因为缺少断点续分析以及重新分析的功能，并且在同一个节点中运行多个样本，所以会容易出现错误
- 内存溢出：常见于markduplicate
- 解决方法：出现错误以后，寻找该步骤对应的样本，删除该样本该步骤及之后步骤所有的结果，重新提交sbatch
- 内存溢出解决方法：调整config中Nproc，降低数目，减少每一步中样本并行数目

# update
该代码因为学习了nextflow以后，发现很好用，所以不更新了，直接全系移植到nextflow流程，能够将运行时间压缩到原本的十分之一左右
详情请见：[nextflow_project](https://github.com/champeil/nextflow_project)
