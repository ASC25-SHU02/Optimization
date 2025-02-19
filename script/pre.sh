#!/bin/bash

hisat-3n-build \
	-p 32 \
	--base-change C,T \
	/asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	/asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

samtools faidx /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' \
	/asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
	> /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf

hisat-3n-build \
	-p 16 \
	--base-change C,T \
	/asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa \
	/asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa

samtools faidx /asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa
