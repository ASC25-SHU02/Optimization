#!/bin/bash

source ./script/func.sh

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
    -p 20 \
    --base-change C,T \
    --mp 8,2 \
    --no-spliced-alignment \
    --directional-mapping \
    | samtools \
    view -@ 20 \
    -e '!flag.unmap' \
    -O BAM \
    -U ../process/SRR23538290/SRR23538290.ncrna.unmapped.bam \
    -o ../process/SRR23538290/SRR23538290.ncrna.mapped.bam'


run_with_timing 'samtools fastq \
	-@ 20 \
	-O ../process/SRR23538290/SRR23538290.ncrna.unmapped.bam \
	> ../process/SRR23538290/SRR23538290.mRNA.fastq'

run_with_timing 'hisat-3n \
	--index ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	-p 20 \
	--summary-file ../process/SRR23538290/map2genome.output.summary \
	--new-summary \
	-q \
	-U ../process/SRR23538290/SRR23538290.mRNA.fastq \
	--directional-mapping \
	--base-change C,T \
	--pen-noncansplice 20 \
	--mp 4,1 \
	| samtools view \
	-@ 20 \
	-e '!flag.unmap'\
       	-O BAM \
	-U ../process/SRR23538290/SRR23538290.mRNA.genome.unmapped.bam \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.bam'


run_with_timing 'samtools sort -@ 20 \
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
	-t 2 -T 20 --data naive --merge avgqual --two-pass \
	-i ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	> ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.log'


run_with_timing 'samtools index \
	-@ 20 \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam.bai'

run_with_timing 'samtools view -@ 20 \
	-e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	-O BAM \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam'

run_with_timing 'samtools view -e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-u --alignments - --ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_unfiltered_uniq.tsv \
	--base-change C,T &

	samtools view -e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_unfiltered_multi.tsv \
	--base-change C,T &

	samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-u --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_filtered_uniq.tsv\
	--base-change C,T &

    samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538290/SRR23538290_filtered_multi.tsv \
	--base-change C,T &
	
	wait'


# You should also comment this command when vtune, and run it after finishing analysing the work
run_with_timing 'python bin/join_pileup.py \
	-i ../process/SRR23538290/SRR23538290_unfiltered_uniq.tsv \
	../process/SRR23538290/SRR23538290_unfiltered_multi.tsv \
	../process/SRR23538290/SRR23538290_filtered_uniq.tsv \
	../process/SRR23538290/SRR23538290_filtered_multi.tsv \
	-o ../process/SRR23538290/SRR23538290_genome.arrow'
