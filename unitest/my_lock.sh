#!/bin/bash
#


# module load samtools/debugVer
# module load testhisat-3n/lshVer

echo "That is what I need"

time {
#stdbuf -oL \
/home/debian/Downloads/hisat-3n/hisat-3n-table \
	-p 4 \
	-m --alignments ./data/Shrinked.in \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name  unitest/check/lock.tsv \
	--base-change C,T
}


# /dev/stdout
# unitest/check/lock.tsv
