#!/bin/bash
# this script is for install miniconda3 and configure the necessary packages
# author: laojp
# time: 2023.12.07
# position: SYSUCC bioinformatic platform
# usage: bash conda_config.sh [download_dir] [install_pos]

down_dir=$(realpath ${1})
install_pos=$(realpath ${2})

wget -c -t 0  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh ${down_dir}
bash ${down_dir}/Miniconda3-latest-Linux-x86_64.sh -b -p ${install_pos}/miniconda
${install_pos}/miniconda/bin/conda init $(echo $SHELL | awk -F "/" '{print $NF}')
echo "Successfully install conda"

