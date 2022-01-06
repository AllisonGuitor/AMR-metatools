#Number is a unique value associated with the sample you are running
NUMBER=$1

echo "=================================== RUN RGI BWT ON TWO PAIRED READS ==================================="
## Options to include at this step, which aligner to use, how many threads to run, whether the RGI database
## is located locally or globally on the system (--local) and whether the variants/wildcard database
## will be included in the analysis (--include_wildcard)
## Start with the trimmed reads optional for the reads to be deduplicated and/or subsampled
## Edit paths to reflect location of reads as well as where the output will be saved.

rgi bwt \
--read_one /path/to/trimmedanddeduplicatedreads1.fastq \
--read_two /path/to/trimmedanddeduplicatedreads2.fastq \
--aligner bowtie2 \
--threads 45 \
--output_file /path/to/output/$1 \
--debug --clean --local --include_wildcard


echo "=================================== RUN RGI KMER_QUERY ON RGI*BWT ==================================="
## To generate the final report output KMER query must be run
## K-mers need to be loaded into RGI
## change the --local flag if needed and # of threads as appropriate

rgi kmer_query \
--input /path/to/output/$1.sorted.length_100.bam \
--kmer_size 61 \
--bwt \
--threads 45 \
--output /path/to/output/$1.sorted.length_100.bam \
--debug --local

echo "=================================== GENERATE CLINICAL REPORT ==================================="
## There are a series of scripts stored under the scripts folder
## These need to to be accessible in the folder where the rgibwtscript is being run

python3 scripts/clinical_report.py \
--bwt_report /path/to/output/$1.gene_mapping_data.txt \
--kmer /path/to/output/$1.sorted.length_100.bam_61mer_analysis.gene.txt \
--output /path/to/output/$1


echo "=================================== PULL READS, ASSEMBLE, AND RUN RGI*MAIN ==================================="

python3 scripts/filter_pull_reads.py \
--bam_file /path/to/output/$1.sorted.length_100.bam \
--json_file /path/to/output/$1.families.json \
--read_one /path/to/trimmedanddeduplicatedreads1.fastq \
--read_two /path/to/trimmedanddeduplicatedreads2.fastq \
--output /path/to/output/$1

echo "===================================  CONSOLIDATE CLINICAL REPORT ==================================="

python3 scripts/consolidate_clinical_report.py \
--raw_clinical_report /path/to/output/$1.output.tsv \
--assembly_report /path/to/output/$1.assembly.json \
--rgi_report /path/to/output/$1.scaffolds.fasta.output.txt \
--bwt_report /path/to/output/$1.gene_mapping_data.txt \
--output /path/to/output/$1


echo "===================================  REPORT 1 ==================================="
python3 scripts/report_one.py \
--final_report /path/to/output/$1.final_output.tsv \
--output /path/to/output/$1

echo "========================  Count Reads Gene_mapping_data =========================="
cat /path/to/output/$1.overall_mapping_stats.txt
wc /path/to/output/$1.gene_mapping_data.txt
additionalscripts/filterbyreadcounts.sh /path/to/output/$1
wc /path/to/output/$1.final_report1.tsv
echo "========================  Assembly and RGI stats =========================="
additionalscripts/stats.py /path/to/output/$1.scaffolds.fasta

echo Count RGI hits (+1)
wc /path/to/output/$1.scaffolds.fasta.output.txt

echo Count Perfect Hits
grep "Perfect" /path/to/output/$1.scaffolds.fasta.output.txt | wc
echo Count Strict Hits
grep "Strict" /path/to/output/$1.scaffolds.fasta.output.txt | wc
echo Count Nudged hits
grep "loose hit" /path/to/output/$1.scaffolds.fasta.output.txt | wc


echo "Done."
