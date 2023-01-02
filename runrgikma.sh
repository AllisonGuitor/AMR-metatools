#!/bin/bash
#Number is a unique value associated with the sample you are running
NUMBER=$1

echo "=================================== RUN RGI BWT WITH KMA ON TWO PAIRED READS =========================="
## Options to include at this step, which aligner to use, how many threads to run, whether the RGI database
## is located locally or globally on the system (--local) and whether the variants/wildcard database
## will be included in the analysis (--include_wildcard)
## Start with the trimmed reads optional for the reads to be deduplicated and/or subsampled
## Edit paths to reflect location of reads as well as where the output will be saved.

rgi bwt \
--read_one ./path/to/trimmedanddeduplicatedreads1.fastq \
--read_two ./path/to/trimmedanddeduplicatedreads2.fastq \
--aligner kma \
--threads 20 \
--output_file ./path/to/output/$1 \
--debug --clean --local
echo "=================================== RUN RGI KMER_QUERY ON RGI*BWT ===================================="
## To generate predictions of potential bacterial hosts of resistance genes
## K-mers need to be loaded into RGI
## change the --local flag if needed and # of threads as appropriate

rgi kmer_query \
--input ./path/to/output/$1.sorted.length_100.bam \
--kmer_size 61 \
--bwt \
--threads 20 \
--output ./path/to/output/$1.sorted.length_100.kmer.bam \
--debug --local

echo "================================== Run assembly and RGI*main ======================================="
## First the metaspades option of SPAdes is used to produce a de novo assembly with default options.
## Next, RGI*main is run to predict ARGs from the assembly.
## Options to include at this step, how many threads to run, whether the RGI database
## is located locally or globally on the system (--local), whether the contigs are of low quality and
## whether to exclude the nudge feature. See RGI documentation on github for more details about parameters.
## Start with the trimmed reads optional for the reads to be deduplicated and/or subsampled
## Edit paths to reflect location of reads as well as where the output will be saved.

metaspades.py -t 24 -m 64 \
-1 ./path/to/trimmedanddeduplicatedreads1.fastq \
-2 ./path/to/trimmedanddeduplicatedreads2.fastq \
-o ./path/to/output/$1

cp ./path/to/output/$1/scaffolds.fasta ./path/to/output/$1.scaffolds.fasta

rgi main -i ./path/to/output/$1.scaffolds.fasta \
-o ./path/to/output/$1.scaffolds.fasta.output --local --debug \
-a blast --clean --exclude_nudge --low_quality

echo "================================  Assembly and RGI stats ========================================="
## Simple stats on the assembly and RGI*main results including the number of ARGs, number of Perfect
## Hits, number of Strict hits, and number of Strict hits that were nudged from Loose to strict.
## The RGI*main results are then filtered for AMR gene families

additionalscripts/stats.py ./path/to/output/$1.scaffolds.fasta

echo "Count RGI hits (+1)"
wc ./path/to/output/$1.scaffolds.fasta.output.txt

echo "Count Perfect Hits"
grep "Perfect" ./path/to/output/$1.scaffolds.fasta.output.txt | wc
echo "Count Strict Hits"
grep "Strict" ./path/to/output/$1.scaffolds.fasta.output.txt | wc
echo "Count Nudged hits"
grep "loose hit" ./path/to/output/$1.scaffolds.fasta.output.txt | wc

echo "Filtering RGI*main results for gene families"
python additionalscripts/card_counts_genefam_rgi.py ./path/to/output/$1.scaffolds.fasta.output.txt

mv genefamilies.csv ./path/to/output/$1-rgi-genefams.csv

mv genefamilyreadcounts.csv ./path/to/output/$1-rgi-genefamcounts.csv

echo "========================  Count Reads Gene_mapping_data =========================="
## First, look at the overall number of reads mapping to CARD with the overall_mapping_stats
## Then count the number of reads with at least 1 read mapping. This will be further filtered.
## Current filtering is for genes with at least 85% coverage and 50 reads or 100%
## coverage and 10 reads.

echo "Percent mapping to CARD"
cat /path/to/output/$1.overall_mapping_stats.txt

echo "Count RGI*BWT hits (+1)"
wc /path/to/output/$1.gene_mapping_data.txt
additionalscripts/filterbyreadcountskma.sh /path/to/output/$1


echo "Done."
