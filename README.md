**Meta-AMR analysis tools**

This repository contains the necessary scripts for running a metagenomic analysis of data to detect antibiotic resistance genes (ARGs)

This suite of commands are ideal for use with MiSeq or HiSeq Illumina sequencing data (2 x 250 bp reads ideally)
This pipeline can be used with either shotgun sequencing data or targeted capture data. 

**Requirements**

- Skewer v0.2.2 tested - https://github.com/relipmoc/skewer 
- BBTools/BBMap v 38.57 tested - dedupe.sh and bbsplitpairs.sh https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/ 
- Optional trimming using seqtk (v 1.3 tested) https://github.com/lh3/seqtk
- RGI v 5.1.1 tested (instructions for download [here](https://github.com/arpcard/rgi)) Installing via Conda is recommended. 
- CARD database downloaded and loaded into RGI. 
-   Tested using CARD v 3.1.0 https://card.mcmaster.ca/download/0/broadstreet-v3.1.0.tar.bz2 and the prevalence/variants database 3.0.7 https://card.mcmaster.ca/download/6/prevalence-v3.0.7.tar.bz2 
- Follow the instructions [here](https://github.com/arpcard/rgi#id42) to load the CARD databases. If you think you will be using RGI and various version of CARD it might be best to load the data locally using the --local flag. 
- Ensure you follow the instruction to load `Additional Reference Data for Metagenomics Analyses` including the steps for `Additional Reference Data for K-mer Pathogen-of-Origin Analyses`
- This pipeline uses RGI bwt or the Metagenomic feature of RGI. More information on it can be found [here](https://github.com/arpcard/rgi#id51). Test that rgi bwt works using 
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

**Step 2 - Generate a all.sh file**
For processing many samples, it is ideal to create a for loop or file to run all samples consecutively. An example file is included in this repository. 

