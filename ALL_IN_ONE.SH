#!/bin/bash
echo \
"============================================================================
YOU ARE RUNNING THE SCRIPTS FROM TEAM ASCeleration, developed by Shuhuai Li
============================================================================

"
# check run mode
if [[ -n "$RUN_MODE" ]]; then
    if [ "$RUN_MODE" = "MANUAL" ]; then
        cd /asc25/Optimization
	bash
    elif [ "$RUN_MODE" = "RERUN" ]; then
        echo "Deleting all the generated files..."
        rm -rf /asc25/process/WT/
        rm -rf /asc25/process/SRR2353829*/*.bam /asc25/process/SRR2353829*/*.tsv /asc25/process/SRR2353829*/*.log /asc25/process/SRR2353829*/*.arrow /asc25/process/SRR2353829*/*.csi /asc25/process/SRR2353829*/*.bai /asc25/process/SRR2353829*/*.summary /asc25/process/SRR2353829*/*.mRNA.fastq /asc25/process/SRR2353829*/*.fastq_*
    elif [ "$RUN_MODE" = "AUTORUN" ]; then
        echo "Running the prebuild process..."
    elif [ "$RUN_MODE" = "PREPROCESS_ONLY" ]; then
        echo "Running the prebuild process only..."
    else
        echo -e "\033[1;31mERROR:\033[0m Invaild RUN_MODE"
        echo
	echo -e \
        "Available options: 
        MANUAL           The DEFAULT option, will NOT do anything when starting the container or script.
        RERUN            \033[1;31m!USE WITH CAUTION!\033[0m Will DELETE ALL OF the generated files.
        PREPROCESS_ONLY  Will only run the pre-process script.
        AUTORUN          Will run the full steps including prebuild and the ultimate test.
        "
        exit 1
    fi
else
    echo -e "\033[1;31mERROR:\033[0m RUN_MODE not defined."
    echo
    echo -e \
    "Available options: 
    MANUAL           The DEFAULT option, will NOT do anything when starting the container or script.
    RERUN            \033[1;31m!USE WITH CAUTION!\033[0m Will DELETE ALL OF the generated files.
    PREPROCESS_ONLY  Will only run the pre-process script.
    AUTORUN          Will run the full steps including prebuild and the ultimate test.
    "
fi

# check dataset
files=(
    "/asc25/process/SRR23538290/SRR23538290.fastq"
    "/asc25/process/SRR23538291/SRR23538291.fastq"
    "/asc25/process/SRR23538292/SRR23538292.fastq"
    "/asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    "/asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa"
)

missing_files=()

for file in "${files[@]}"; do
    if ! [ -e "$file" ]; then
        missing_files+=("$file")
    fi
done

if [ ${#missing_files[@]} -eq 0 ]; then
    echo "All datasets exist"
else
    echo "The following dataset do NOT exist:"
    for file in "${missing_files[@]}"; do
        echo "$file"
    done
    echo "Make sure you have mounted all the datasets correctly."
fi

# create the log folder
mkdir -p log
# check the files before cutseq.
function check_files { 
    refFile="/asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    ncrnaRefFile="/asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa.fai"

    if [ ! -f "$refFile" ] || [ ! -f "$ncrnaRefFile" ]; then
        echo "It seems that the script \`/asc25/Optimization/script/pre.sh\` hasn't been run. MAKE SURE YOU ARE NOT UNDER ULTIMATE TEST!"
        bash /asc25/Optimization/script/pre.sh
	if [ "$RUN_MODE" = "PREPROCESS_ONLY" ]; then
	    echo "Pre-process finished"
	    exit 0
	fi
    elif [ "$RUN_MODE" = "PREBUILD_ONLY" ]; then
	    exit 0
    else
        echo "Files exist. Directly run the file. You may on test mode"
    fi
}

check_files


function clear_files {
    for sample in SRR23538290 SRR23538291 SRR23538292; do
        find /asc25/process/${sample} -type f -not -name ${sample}.fastq -delete
    done
    if [[ -d "/asc25/process/WT" ]]; then
        rm -rf /asc25/process/WT
    else
        mkdir -p /asc25/process/WT
    fi
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
	echo "Timing enabled. See logs in /asc25/Optimization/log/time.log"
    log_file="/asc25/Optimization/log/time.log"
    exec > >(tee -a $log_file)
    exec 2>&1
else 
	echo "Timing disabled. See logs in /asc25/Optimization/log/plain.log"
    log_file="/asc25/Optimization/log/plain.log"
    exec > >(tee -a $log_file)
    exec 2>&1
fi

time {
time {
    for sample in SRR23538290 SRR23538291 SRR23538292; do
        bash /asc25/Optimization/script/stage1_${sample}.sh
    done

    echo "@@@@@@@@@@@@@@@TIME SPENDT IN STAGE1@@@@@@@@@@@@@@@@"
}

echo "
*********************************************************
**********************STAGE2*****************************
*********************************************************

"
time {
    bash /asc25/Optimization/script/stage2.sh
    echo "@@@@@@@@@@@@@@@TIME SPENDT IN STAGE2@@@@@@@@@@@@@@@@"
}
echo "!!!!!!!!!!!!!!!!!!!Workflow End!!!!!!!!!!!!!!!!!!!!!!!!!!"
}
answer="yes"
if [[ "$answer" == "yes" ]]; then
    cd /asc25/Optimization/check
    bash /asc25/Optimization/check/check.sh
fi
