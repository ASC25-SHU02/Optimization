#!/bin/bash

refFile="../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ncrnaRefFile="../ref/Homo_sapiens.GRCh38.ncrna.fa"

if [ ! -f "$refFile" ] || [! -f "$ncrnaRefFile"]; then
    echo "It seems that the script \`script/pre.sh\` hasn't been run. MAKE SURE YOU ARE NOT UNDER ULTIMATE TEST!"
    bash script/pre.sh
else
    echo "Files exist. Directly run the file. You may on test mode"
fi


echo "
*********************************************************
**********************STAGE1*****************************
*********************************************************

"

time {
    bash script/stage1_SRR23538290.sh
    bash script/stage1_SRR23538291.sh
    bash script/stage1_SRR23538292.sh
}

echo "
*********************************************************
**********************STAGE2*****************************
*********************************************************

"

time {
    bash script/stage2.sh
}

read -p "You wanna check your answer?(yes/no): " answer
if [[ "$answer" == "yes" ]]; then
    cd check
    bash ./check.sh
fi