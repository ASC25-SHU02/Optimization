#!/bin/bash
echo \
"============================================================================
YOU ARE RUNNING THE SCRIPTS FROM TEAM ASCeleration, developed by Shuhuai Li
============================================================================

"
# check the files before cutseq.
function check_files { 
    refFile="../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    ncrnaRefFile="../ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa"

    if [ ! -f "$refFile" ] || [ ! -f "$ncrnaRefFile" ]; then
        echo "It seems that the script \`script/pre.sh\` hasn't been run. MAKE SURE YOU ARE NOT UNDER ULTIMATE TEST!"
        bash script/pre.sh
    else
        echo "Files exist. Directly run the file. You may on test mode"
    fi
}

check_files

echo "
*********************************************************
**********************STAGE1*****************************
*********************************************************

"

check_tools() {
	echo "Checking the tools..."
    for tool in cutseq hisat-3n samtools java python; do
        if ! command -v $tool &> /dev/null; then
            echo "$tool could not be found"
            exit 1
        fi
    done
	echo "All tools available"
}

check_tools

read -p "Do you want to enable timing for each step with tiny overhead? (yes/no): " enable_timing
export enable_timing

if [ "$enable_timing" == "yes" ]; then
	echo "Timing enabled. See logs in ./log/time.log"
    log_file="./log/stage1.sh.log"
    exec > >(tee -a $log_file)
    exec 2>&1
else 
	echo "Timing disabled. See logs in ./log/plain.log"
    log_file="./log/plain.log"
    exec > >(tee -a $log_file)
    exec 2>&1
fi

time {
    bash script/stage1_SRR23538290.sh enable_timing
    bash script/stage1_SRR23538291.sh enable_timing
    bash script/stage1_SRR23538292.sh enable_timing

    echo "@@@@@@@@@@@@@@@TIME SPENDT IN STAGE1@@@@@@@@@@@@@@@@"
}

echo "
*********************************************************
**********************STAGE2*****************************
*********************************************************

"

time {
    bash script/stage2.sh
    echo "@@@@@@@@@@@@@@@TIME SPENDT IN STAGE2@@@@@@@@@@@@@@@@"
}

read -p "You wanna check your answer?(yes/no): " answer
if [[ "$answer" == "yes" ]]; then
    cd check
    bash ./check.sh
fi