hisat-3n-build -p 32 --base-change C,T ../processed/rep1/SRR23538290.fasta ../processed/rep1/rep1
samtools faidx ../processed/rep1/SRR23538290.fasta

awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' ../processed//SRR23538290.fasta.fai > ../processed/rep1/rep1.fa.saf 
