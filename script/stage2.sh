#!/bin/bash


python bin/group_pileup.py \
    -i ../process/SRR23538290/SRR23538290_genome.arrow \
    ../process/SRR23538291/SRR23538291_genome.arrow \
    ../process/SRR23538292/SRR23538292_genome.arrow \
    -o ../process/WT/WT.arrow

python bin/select_sites.py -i output/WT.arrow -o ../process/WT/WT.prefilter.tsv

python bin/filter_sites.py \
    -i ../process/SRR23538290/SRR23538290_genome.arrow \
    -m ../process/WT/WT.prefilter.tsv \
    -b ../process/SRR23538290/SRR23538290.bg.tsv \
    -o ../process/SRR23538290/SRR23538290.filtered.tsv

python bin/filter_sites.py \
    -i ../process/SRR23538291/SRR23538291_genome.arrow \
    -m ../process/WT/WT.prefilter.tsv \
    -b ../process/SRR23538291/SRR23538291.bg.tsv \
    -o ../process/SRR23538291/SRR23538291.filtered.tsv

python bin/filter_sites.py \
    -i ../process/SRR23538292/SRR23538292_genome.arrow \
    -m ../process/WT/WT.prefilter.tsv \
    -b ../process/SRR23538292/SRR23538292.bg.tsv \
    -o ../process/SRR23538292/SRR23538292.filtered.tsv
