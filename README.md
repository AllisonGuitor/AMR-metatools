**Meta-AMR analysis tools**

This repository contains the necessary scripts for running a metagenomic analysis of data to detect antibiotic resistance genes (ARGs)

This suite of commands are ideal for use with MiSeq or HiSeq Illumina sequencing data (2 x 250 bp reads ideally)
This pipeline can be used with either shotgun sequencing data or targeted capture data. 

**Step 1 - Trim reads**

Use the processreads.sh script to trim reads using [skewer](https://doi.org/10.1186/1471-2105-15-182) (Version 0.2.2 was used for this work).
This bash script can be run in a folder above the folder containing the raw reads for a given sample. 

`./processreads.sh $1` 

Where $1 is the name of the folder containing the raw reads. 

Other read trimming software can be used at this step. 

**Step 2 - Generate a all.sh file**
For processing many samples, it is ideal to create a for loop or file to run all samples consecutively. An example file is included in this repository. 

