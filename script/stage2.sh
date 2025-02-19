#!/bin/bash


python bin/group_pileup.py \
     -i /asc25/process/SRR23538290/SRR23538290_genome.arrow \
     /asc25/process/SRR23538291/SRR23538291_genome.arrow \
     /asc25/process/SRR23538292/SRR23538292_genome.arrow \
     -o /asc25/process/WT/WT.arrow

python bin/select_sites.py -i /asc25/process/WT/WT.arrow -o /asc25/process/WT/WT.prefilter.tsv

python bin/filter_sites.py \
    -i /asc25/process/SRR23538290/SRR23538290_genome.arrow \
    -m /asc25/process/WT/WT.prefilter.tsv \
    -b /asc25/process/SRR23538290/SRR23538290.bg.tsv \
    -o /asc25/process/SRR23538290/SRR23538290.filtered.tsv

python bin/filter_sites.py \
    -i /asc25/process/SRR23538291/SRR23538291_genome.arrow \
    -m /asc25/process/WT/WT.prefilter.tsv \
    -b /asc25/process/SRR23538291/SRR23538291.bg.tsv \
    -o /asc25/process/SRR23538291/SRR23538291.filtered.tsv

python bin/filter_sites.py \
    -i /asc25/process/SRR23538292/SRR23538292_genome.arrow \
    -m /asc25/process/WT/WT.prefilter.tsv \
    -b /asc25/process/SRR23538292/SRR23538292.bg.tsv \
    -o /asc25/process/SRR23538292/SRR23538292.filtered.tsv
