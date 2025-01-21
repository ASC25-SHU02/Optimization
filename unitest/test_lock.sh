#!/bin/bash
#


# module load samtools/debugVer
# module load testhisat-3n/lshVer


time {
/home/debian/Downloads/hisat-3n/hisat-3n-table \
	-p 16 \
	-m --alignments ./data/Shrinked.in \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /dev/stdout \
	--base-change C,T \
	| cut -f 1,2,3,5,7 \
	> 'unitest/check/true.tsv'
}
