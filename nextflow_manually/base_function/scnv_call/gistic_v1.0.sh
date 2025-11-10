#!/bin/bash
# this script is for significant cnv detection
# time: 2023.06.16
# position: SYSUCC bioinformatic platform
# author: laojp
# usage: process_gistic [seg_suffix] [input_dir] [output_dir]

function process_gistic {
    echo -e "\t\t[base function] process_gistic"
    echo -e "#!/bin/bash" >> ${3}/gistic.sh
    echo -e "# this script is for gistic to detect copy number" >> ${3}/gistic.sh
    echo -e "# author: laojp\n# time: $(date)\n# position: SYSUCC bioinformatic platform\n" >> ${3}/gistic.sh

    echo -e "basedir=${3}" >> ${3}/gistic.sh
    echo -e "segfile=$(realpath ${2}/${1}.seg)" >> ${3}/gistic.sh
    echo -e "refgenefile=${gistic}/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat" >> ${3}/gistic.sh
    echo -e '${gistic}/gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.99 -armpeel 1 -savegene 1 -gcm extreme -qvt 0.05 -saveseg 1 -v 30 -ta 0.1 -td 0.1' >> ${3}/gistic.sh

    cp ${3}/gistic.sh ${gistic}/
    cd ${gistic}
    bash gistic.sh 1>${3}/gistic_log 2>&1 
    process_exist_status_background "gistic finished"
    cd - 
}