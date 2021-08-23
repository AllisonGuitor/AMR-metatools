#!/bin/bash
echo using prefix $1

## The line below might need to be adjusted depending on the output of the sequencing data. 
gunzip $1/*R*fastq.gz

cat $1/*R1*fastq >> $1/R1.fastq
cat $1/*R2*fastq >> $1/R2.fastq

skewer -m pe -q 25 -Q 25 $1/R1.fastq $1/R2.fastq -o $1/skewer --threads 40

dedupe.sh in=$1/skewer-trimmed-pair1.fastq,$1/skewer-trimmed-pair2.fastq out=$1/trimmed-dedupe.fastq ac=f mlp=100 mop=100 -Xmx20g

bbsplitpairs.sh in=$1/trimmed-dedupe.fastq out=$1/trimmed-dedupe-R1.fastq out2=$1/trimmed-dedupe-R2.fastq 
