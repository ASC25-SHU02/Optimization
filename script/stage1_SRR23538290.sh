#!/bin/bash

cutseq ../process/SRR23538290/SRR23538290.fasta \
        -t 20 \
        -A INLINE \
        -m 20 \
        --trim-polyA \
        --ensure-inline-barcode\
        -o ../process/SRR23538290/SRR23538290.fastq_cut \
        -s ../process/SRR23538290/SRR23538290.fastq_tooshort \
        -u ../process/SRR23538290/SRR23538290.fastq_untrimmed

hisat-3n --index ../ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa \
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
    -o ../process/SRR23538290/SRR23538290.ncrna.mapped.bam


samtools fastq \
	-@ 16 \
	-O ../process/SRR23538290/SRR23538290.ncrna.unmapped.bam \
	> ../process/SRR23538290/SRR23538290.mRNA.fastq

hisat-3n \
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
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.bam

echo "FINISH ALL HISAT-3N & SAMTOOL VIEW. GET: mapped.bam"

samtools sort -@ 16 \
	--write-index \
	-O BAM \
	-o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.bam

echo "FINISH [SAMTOOL] SORT. GET: mapped.sorted.bam"

samtools view \
	-@ 20 \
	-F 3980 \
	-c ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam \
	> ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam.tsv

echo "FINISH SAMTOOLS VIEW. GET: mRNA.genome.mapped.sorted.bam.tsv"

java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=../process/SRR23538290 -jar ../UMICollapse/umicollapse.jar bam -t 2 -T 16 --data naive --merge avgqual --two-pass -i ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.bam -o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam > ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.log

echo "FINISH JAVA"

samtools index \
	-@ 8 \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam.bai

echo "FINISH SAMTOOL INDEX"

samtools view -e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-p 16 \
	-u --alignments - --ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /dev/stdout \
	--base-change C,T \
	| cut -f 1,2,3,5,7 \
	| gzip -c \
	> ../process/SRR23538290/SRR23538290_unfiltered_uniq.tsv.gz

echo "FINISH SAMTOOL. Get: unfiltered_uniq.tsv.gz"

samtools view -e "rlen<100000" \
	-h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table -p 16 \
	-m --alignments - \
	--ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /dev/stdout \
	--base-change C,T \
	| cut -f 1,2,3,5,7 \
	| gzip -c \
	> ../process/SRR23538290/SRR23538290_unfiltered_multi.tsv.gz

echo "FINISH SAMTOOL. Get: unfiltered_multi.tsv.gz"

samtools view -@ 8 -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.bam -O BAM -o ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam

echo "FINISH SAMTOOL [XM] * 20 <=. Get: mapped.sorted.dedup.filtered.bam"

samtools view -e "rlen<100000" -h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam | hisat-3n-table -p 16 -u --alignments - --ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > ../process/SRR23538290/SRR23538290_filtered_uniq.tsv.gz

echo "FINISH SAMTOOL rlen<100000. Get: filtered_uniq.tsv.gz"

samtools view -e "rlen<100000" -h ../process/SRR23538290/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam | hisat-3n-table -p 16 -m --alignments - --ref ../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | gzip -c > ../process/SRR23538290/SRR23538290_filtered_multi.tsv.gz

echo "FINISH SAMTOOL rlen<100000. Get: filtered_multi.tsv.gz"
python bin/join_pileup.py -i ../process/SRR23538290/SRR23538290_unfiltered_uniq.tsv.gz ../process/SRR23538290/SRR23538290_unfiltered_multi.tsv.gz ../process/SRR23538290/SRR23538290_filtered_uniq.tsv.gz ../process/SRR23538290/SRR23538290_filtered_multi.tsv.gz -o ../process/SRR23538290/SRR23538290_genome.arrow

echo "FINISH PYTHON. Get: genome.arrow"
