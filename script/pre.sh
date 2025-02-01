#!/bin/bash

hisat-3n-build \
	-p 32 \
	--base-change C,T \
	../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

samtools faidx ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' \
	../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
	> ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf

hisat-3n-build \
	-p 16 \
	--base-change C,T \
	../ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa \
	../ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa

samtools faidx ../ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa
