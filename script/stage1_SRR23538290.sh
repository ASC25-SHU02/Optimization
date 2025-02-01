#!/bin/bash


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

echo "============================================================================"
echo "YOU ARE RUNNING THE SCRIPTS FROM TEAM ASCeleration, developed by Shuhuai Li"
echo "============================================================================"

read -p "Do you want to enable timing for each step with tiny overhead? (yes/no): " enable_timing

if [ "$enable_timing" == "yes" ]; then
	echo "**Timing enabled. See logs in ./log/stage1_SRR23538290.sh.log**"
    log_file="./log/stage1_SRR23538290.sh.log"
    exec > >(tee -a $log_file)
    exec 2>&1
    start_time=$(date +%s)
else 
	echo "You didn't enable tming"
fi

run_with_timing() {
	local cmd="$@"
    echo "Starting: $cmd"
    if [ "$enable_timing" == "yes" ]; then
        time eval "$@"
    else
        eval "$@"
    fi
	echo "##################...Ending...##################"
	echo
}

# Comment this command when running vtune, run it separately before the workflow
run_with_timing 'cutseq ../process/SRR23538290/SRR23538290.fastq \
        -t 20 \
        -A INLINE \
        -m 20 \
        --trim-polyA \
        --ensure-inline-barcode\
        -o ../process/SRR23538290/SRR23538290.fastq_cut \
        -s ../process/SRR23538290/SRR23538290.fastq_tooshort \
        -u ../process/SRR23538290/SRR23538290.fastq_untrimmed'

run_with_timing 'hisat-3n --index ../ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa \
    --summary-file ../process/SRR23538290/map2ncrna.output.summary \
    --new-summary \
    -q \
    -U ../process/SRR23538290/SRR23538290.fastq_cut\
    -p 16 \
    --base-change C,T \
    --mp 8,2 \
    --no-spliced-alignment \
    --directional-mapping \
    | samtools \
    view -@ 16 \
    -e '!flag.unmap' \
    -O BAM \
    -U ../process/SRR23538290/SRR23538290.ncrna.unmapped.bam \
    -o ../process/SRR23538290/SRR23538290.ncrna.mapped.bam'


run_with_timing 'samtools fastq \
	-@ 16 \
	-O ../process/SRR23538290/SRR23538290.ncrna.unmapped.bam \
	> ../process/SRR23538290/SRR23538290.mRNA.fastq'

run_with_timing 'hisat-3n \
	--index ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	-p 16 \
	--summary-file ../process/SRR23538290/map2genome.output.summary \
	--new-summary \
	-q \
	-U ../process/SRR23538290/SRR23538290.mRNA.fastq \
	--directional-mapping \
	--base-change C,T \
	--pen-noncansplice 20 \
	--mp 4,1 \
	| samtools view \
	-@ 16 \
	-e '!flag.unmap'\
       	-O BAM \
	-U ../process/SRR23538290/SRR23538290.mRNA.genome.unmapped.bam \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.bam'


run_with_timing 'samtools sort -@ 16 \
	--write-index \
	-O BAM \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.bam'


run_with_timing 'samtools view \
	-@ 20 \
	-F 3980 \
	-c ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam \
	> ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam.tsv'


run_with_timing 'java -server -Xms8G -Xmx40G -Xss100M \
	-Djava.io.tmpdir=../process/SRR23538290 \
	-jar ../UMICollapse/umicollapse.jar bam \
	-t 2 -T 16 --data naive --merge avgqual --two-pass \
	-i ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	> ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.log'


run_with_timing 'samtools index \
	-@ 8 \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam.bai'

run_with_timing 'samtools view -e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-p 4 \
	-u --alignments - --ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_unfiltered_uniq.tsv \
	--base-change C,T'


run_with_timing 'samtools view -e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table -p 4 \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_unfiltered_multi.tsv \
	--base-change C,T'


run_with_timing 'samtools view -@ 8 \
	-e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	-O BAM \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam'


run_with_timing 'samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-p 4 \
	-u --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_filtered_uniq.tsv\
	--base-change C,T'


run_with_timing 'samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table -p 4\
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_filtered_multi.tsv \
	--base-change C,T'


# You should also comment this command when vtune, and run it after finishing analysing the work
run_with_timing 'python bin/join_pileup.py \
	-i ../process/SRR23538290/SRR23538290_unfiltered_uniq.tsv \
	../process/SRR23538290/SRR23538290_unfiltered_multi.tsv \
	../process/SRR23538290/SRR23538290_filtered_uniq.tsv \
	../process/SRR23538290/SRR23538290_filtered_multi.tsv \
	-o ../process/SRR23538290/SRR23538290_genome.arrow'


if [ "$enable_timing" == "yes" ]; then
    end_time=$(date +%s)
    total_time=$((end_time - start_time))
    echo "Total execution time: $total_time seconds" | tee -a $log_file
fi

echo "FINISH ALL STEPS."