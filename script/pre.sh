#!/bin/bash

hisat-3n-build \
    -p 12 \
    --base-change C,T\
    ../process/rep1/SRR23538290.fasta \
    ../process/rep1/rep1

samtools faidx ../process/rep1/SRR23538290.fasta

awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' ../process/SRR23538290.fasta.fai > ../process/rep1/rep1.fa.saf 
