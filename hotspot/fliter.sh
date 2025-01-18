#!/bin/bash
#


# module load samtools/debugVer
# module load hisat-3n/lshVer


time {
	samtools view -e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	| /home/debian/Downloads/hisat-3n/hisat-3n-table -p 16 \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /dev/stdout \
	--base-change C,T \
	| cut -f 1,2,3,5,7 \
	| gzip -c \
	> ../process/SRR23538290/SRR23538290_unfiltered_multi.tsv.gz
} 2> timeInfo/filter.txt
