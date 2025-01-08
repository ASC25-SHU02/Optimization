#!/bin/bash

cutseq ../process/rep1/rep1.fasta \
        -t 20 \
        -A INLINE \
        -m 20 \
        --trim-polyA \
        --ensure-inline-barcode\
        -o ../process/rep1/rep1.fastq_cut \
        -s ../process/rep1/rep1.fastq_tooshort \
        -u ../process/rep1/rep1.fastq_untrimmed

hisat-3n -x ../process/rep1/rep1 \
    --summary-file ../process/rep1/map2poly.output.summary \
    --new-summary \
    -q \
    -U ../process/rep1/rep1.fastq_cut\
    -p 16 \
    --base-change C,T \
    --mp 8,2 \
    --no-spliced-alignment \
    --directional-mapping \
    --large-index \
    | samtools \
    view -@ 16 \
    -e '!flag.unmap' \
    -O BAM \
    -U ../process/rep1/rep1.ncrna.unmapped.bam \
    -o ../process/rep1/rep1.ncrna.mapped.bam
