#!/bin/bash
echo \
"============================================================================
YOU ARE RUNNING THE SCRIPTS FROM TEAM ASCeleration, developed by Shuhuai Li
============================================================================

"
# create the log folder
mkdir -p log
# check the files before cutseq.
function check_files { 
    refFile="../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    ncrnaRefFile="../ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa.fai"

    if [ ! -f "$refFile" ] || [ ! -f "$ncrnaRefFile" ]; then
        echo "It seems that the script \`script/pre.sh\` hasn't been run. MAKE SURE YOU ARE NOT UNDER ULTIMATE TEST!"
        bash script/pre.sh
    else
        echo "Files exist. Directly run the file. You may on test mode"
    fi
}

check_files


function clear_files {
    for sample in SRR23538290 SRR23538291 SRR23538292; do
        find ../process/${sample} -type f -not -name ${sample}.fastq -delete
    done
    rm -rf ../process/WT
    mkdir ../process/WT
}

clear_files
echo "
*********************************************************
**********************STAGE1*****************************
*********************************************************
!!!!!!!!!!!!!!!!!!!cutseq begin!!!!!!!!!!!!!!!!!!!!!!!!!!
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

enable_timing="yes"
export enable_timing

if [ "$enable_timing" == "yes" ]; then
	echo "Timing enabled. See logs in ./log/time.log"
    log_file="./log/time.log"
    exec > >(tee -a $log_file)
    exec 2>&1
else 
	echo "Timing disabled. See logs in ./log/plain.log"
    log_file="./log/plain.log"
    exec > >(tee -a $log_file)
    exec 2>&1
fi

time {
time {
    for sample in SRR23538290 SRR23538291 SRR23538292; do
        bash script/stage1_${sample}.sh
    done

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
echo "!!!!!!!!!!!!!!!!!!!Workflow End!!!!!!!!!!!!!!!!!!!!!!!!!!"
}
answer="yes"
if [[ "$answer" == "yes" ]]; then
    cd check
    bash ./check.sh
fi
