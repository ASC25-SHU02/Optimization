#!/bin/bash

source ./script/func.sh

# Comment this command when running vtune, run it separately before the workflow
run_with_timing 'cutseq /asc25/process/SRR23538292/SRR23538292.fastq \
        -t 20 \
        -A INLINE \
        -m 20 \
        --trim-polyA \
        --ensure-inline-barcode\
        -o /asc25/process/SRR23538292/SRR23538292.fastq_cut \
        -s /asc25/process/SRR23538292/SRR23538292.fastq_tooshort \
        -u /asc25/process/SRR23538292/SRR23538292.fastq_untrimmed'

run_with_timing 'hisat-3n --index /asc25/ncrna_ref/Homo_sapiens.GRCh38.ncrna.fa \
    --summary-file /asc25/process/SRR23538292/map2ncrna.output.summary \
    --new-summary \
    -q \
    -U /asc25/process/SRR23538292/SRR23538292.fastq_cut\
    -p 20 \
    --base-change C,T \
    --mp 8,2 \
    --no-spliced-alignment \
    --directional-mapping \
    | samtools \
    view -@ 20 \
    -e '!flag.unmap' \
    -O BAM \
    -U /asc25/process/SRR23538292/SRR23538292.ncrna.unmapped.bam \
    -o /asc25/process/SRR23538292/SRR23538292.ncrna.mapped.bam'


run_with_timing 'samtools fastq \
	-@ 20 \
	-O /asc25/process/SRR23538292/SRR23538292.ncrna.unmapped.bam \
	> /asc25/process/SRR23538292/SRR23538292.mRNA.fastq'

run_with_timing 'hisat-3n \
	--index /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	-p 20 \
	--summary-file /asc25/process/SRR23538292/map2genome.output.summary \
	--new-summary \
	-q \
	-U /asc25/process/SRR23538292/SRR23538292.mRNA.fastq \
	--directional-mapping \
	--base-change C,T \
	--pen-noncansplice 20 \
	--mp 4,1 \
	| samtools view \
	-@ 20 \
	-e '!flag.unmap'\
       	-O BAM \
	-U /asc25/process/SRR23538292/SRR23538292.mRNA.genome.unmapped.bam \
	-o /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.bam'


run_with_timing 'samtools sort -@ 20 \
	--write-index \
	-O BAM \
	-o /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.bam \
	/asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.bam'


run_with_timing 'samtools view \
	-@ 20 \
	-F 3980 \
	-c /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.bam \
	> /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.bam.tsv'


run_with_timing 'java -server -Xms8G -Xmx40G -Xss100M \
	-Djava.io.tmpdir=/asc25/process/SRR23538292 \
	-jar /Tools/UMICollapse/umicollapse.jar bam \
	-t 2 -T 20 --data naive --merge avgqual --two-pass \
	-i /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.bam \
	-o /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam \
	> /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.log'


run_with_timing 'samtools index \
	-@ 20 \
	/asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam \
	/asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam.bai'

run_with_timing 'samtools view -@ 20 \
	-e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" \
	/asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam \
	-O BAM \
	-o /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.filtered.bam'

run_with_timing 'samtools view -e "rlen<100000" \
	-h /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-u --alignments - --ref /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /asc25/process/SRR23538292/SRR23538292_unfiltered_uniq.tsv \
	--base-change C,T &

	samtools view -e "rlen<100000" \
	-h /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /asc25/process/SRR23538292/SRR23538292_unfiltered_multi.tsv \
	--base-change C,T &

	samtools view \
	-e "rlen<100000" \
	-h /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-u --alignments - \
	--ref /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /asc25/process/SRR23538292/SRR23538292_filtered_uniq.tsv\
	--base-change C,T &

    samtools view \
	-e "rlen<100000" \
	-h /asc25/process/SRR23538292/SRR23538292.mRNA.genome.mapped.sorted.dedup.filtered.bam \
	| hisat-3n-table \
	-m --alignments - \
	--ref /asc25/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--output-name /asc25/process/SRR23538292/SRR23538292_filtered_multi.tsv \
	--base-change C,T &
	
	wait'


# You should also comment this command when vtune, and run it after finishing analysing the work
run_with_timing 'python bin/join_pileup.py \
	-i /asc25/process/SRR23538292/SRR23538292_unfiltered_uniq.tsv \
	/asc25/process/SRR23538292/SRR23538292_unfiltered_multi.tsv \
	/asc25/process/SRR23538292/SRR23538292_filtered_uniq.tsv \
	/asc25/process/SRR23538292/SRR23538292_filtered_multi.tsv \
	-o /asc25/process/SRR23538292/SRR23538292_genome.arrow'
