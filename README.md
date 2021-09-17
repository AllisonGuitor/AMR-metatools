**Meta-AMR analysis tools**

This repository contains the necessary scripts for running a metagenomic analysis of sequencing data to detect antibiotic resistance genes (ARGs).

This suite of commands is ideal for use with MiSeq or HiSeq Illumina sequencing data (2 x 250 bp reads ideally)
This pipeline can be used with either shotgun sequencing data or targeted capture data. 

**Requirements**

- Skewer v0.2.2 tested - https://github.com/relipmoc/skewer 
- BBTools/BBMap v 38.57 tested - dedupe.sh and bbsplitpairs.sh https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/ 
- Optional subsampling using seqtk (v 1.3 tested) https://github.com/lh3/seqtk
- RGI v 5.1.1 tested (instructions for download [here](https://github.com/arpcard/rgi)) Installing via Conda is recommended. 
- CARD database downloaded and loaded into RGI. 
-   Tested using CARD v 3.1.0 https://card.mcmaster.ca/download/0/broadstreet-v3.1.0.tar.bz2 and the prevalence/variants database 3.0.7 https://card.mcmaster.ca/download/6/prevalence-v3.0.7.tar.bz2 
- Follow the instructions [here](https://github.com/arpcard/rgi#id42) to load the CARD databases. If you think you will be using RGI and various version of CARD, it might be best to load the data locally using the --local flag. 
- Ensure you follow the instructions to load `Additional Reference Data for Metagenomics Analyses` including the steps for `Additional Reference Data for K-mer Pathogen-of-Origin Analyses`
- This pipeline uses RGI*BWT or the Metagenomic feature of RGI. More information on it can be found [here](https://github.com/arpcard/rgi#id51). After installing RGI, test that RGI*BWT works using: 
```
rgi bwt --help
```
If databases are installed localled, move to the folder that contains the `/localDB/` folder and check
```
rgi databases -v --all --local
```

**Step 1 - Trim reads**

Use the processreads.sh script to trim reads using skewer and remove duplicate reads (optional) using dedupe.sh. 
This bash script can be run in a folder above the folder containing the raw reads for a given sample. 

```
./processreads.sh $1
``` 

Where $1 is the name of the folder containing the raw reads. 

Other read trimming and deduplication software can be used at this step. 

**Optional Step - Subsampling**

Reads can be subsampled in any way if the user wants. We used `seqtk sample`. The following can be generated as a bash script. Example with subsampling to 50,000 reads (paired)

```
#!/bin/bash
echo using prefix $1
 
seqtk sample -s 11 $1/*trimmed-dedupe-*1.fastq 50000 > $1/trimmed-dedupe-R1-50K.fastq

seqtk sample -s 11 $1/*trimmed-dedupe-*2.fastq 50000 > $1/trimmed-dedupe-R2-50K.fastq
```

**Step 2 - Generate a `all.sh` file**

For processing many samples, it is ideal to create a for loop or file to run all samples consecutively. 

```
./runrgibwt.sh SAMPLEID > SAMPLEID.log 2>&1
```

**Step 3 - Run all.sh file**

Ensure the folder in which you are running the all.sh command contains the [runrgibwt.sh](meta-tools/runrgibwt.sh) script with the correct paths to your reads. 

- Lines [13/14](meta-tools/runrgibwt.sh#L13) should be changed to the direct path of the reads along with lines [49/50](meta-tools/runrgibwt.sh#L49). 
- Line [17](meta-tools/runrgibwt.sh#L17) and subsequent lines should direct to the output location. 
- Line 18 can be modified in case the CARD reference is installed globally (remove the --local flag) and if the wildcard/variants db is not included (remove the --wildcard flag)

If the database is not installed locally there are a few lines that need to be changed in the [filter_pull_reads.py script](meta-tools/scripts/filter_pull_reads.py)

- Line 68, 161 and 166 run rgi main/kmer query - remove the --local flag if the db is loaded globally
Also check the following lines: 
- Line 61 and Line 150 - ensure the correct path to metaspades and that the memory and threads flags are appropriate




