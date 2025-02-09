#!/bin/bash
source ./script/func.sh

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
	
	samtools view -e "rlen<100000" \
	-h ../process/SRR23538291/SRR23538291.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-u --alignments - --ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538291/SRR23538291_unfiltered_uniq.tsv \
	--base-change C,T &

	samtools view -e "rlen<100000" \
	-h ../process/SRR23538291/SRR23538291.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538291/SRR23538291_unfiltered_multi.tsv \
	--base-change C,T &

	samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538291/SRR23538291.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-u --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538291/SRR23538291_filtered_uniq.tsv\
	--base-change C,T &

    samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538291/SRR23538291.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538291/SRR23538291_filtered_multi.tsv \
	--base-change C,T &
	
	samtools view -e "rlen<100000" \
	-h ../process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-u --alignments - --ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538292/SRR23538292_unfiltered_uniq.tsv \
	--base-change C,T &

	samtools view -e "rlen<100000" \
	-h ../process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538292/SRR23538292_unfiltered_multi.tsv \
	--base-change C,T &

	samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-u --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538292/SRR23538292_filtered_uniq.tsv\
	--base-change C,T &

    samtools view \
	-e "rlen<100000" \
	-h ../process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name ../process/SRR23538292/SRR23538292_filtered_multi.tsv \
	--base-change C,T &
	
	wait'


# You should also comment this command when vtune, and run it after finishing analysing the work
run_with_timing 'python bin/join_pileup.py \
	-i ../process/SRR23538290/SRR23538290_unfiltered_uniq.tsv \
	../process/SRR23538290/SRR23538290_unfiltered_multi.tsv \
	../process/SRR23538290/SRR23538290_filtered_uniq.tsv \
	../process/SRR23538290/SRR23538290_filtered_multi.tsv \
	-o ../process/SRR23538290/SRR23538290_genome.arrow &
    
    python bin/join_pileup.py \
	-i ../process/SRR23538291/SRR23538291_unfiltered_uniq.tsv \
	../process/SRR23538291/SRR23538291_unfiltered_multi.tsv \
	../process/SRR23538291/SRR23538291_filtered_uniq.tsv \
	../process/SRR23538291/SRR23538291_filtered_multi.tsv \
	-o ../process/SRR23538291/SRR23538291_genome.arrow &

    python bin/join_pileup.py \
	-i ../process/SRR23538292/SRR23538292_unfiltered_uniq.tsv \
	../process/SRR23538292/SRR23538292_unfiltered_multi.tsv \
	../process/SRR23538292/SRR23538292_filtered_uniq.tsv \
	../process/SRR23538292/SRR23538292_filtered_multi.tsv \
	-o ../process/SRR23538292/SRR23538292_genome.arrow &

wait'
