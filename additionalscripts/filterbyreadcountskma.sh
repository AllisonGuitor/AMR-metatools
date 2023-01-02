#!/bin/bash

## Additional script to filter the RGI*bwt kma mapping results based on number of reads
## mapping and the percent length of coverage of the reference gene.
## Current filters are for genes with at least 85% coverage and 50 reads OR
## 100 percent coverage and 10 reads. These cut-offs can be modified in the parsereadcountskma.py
## script.

echo sample $1

python additionalscripts/parsereadcountskma.py $1.gene_mapping_data.txt

mv parsedrgibwtresults-85perand50or100perand10.csv $1-gene_mapping_85perand50or100perand10.csv

python additionalscripts/getreadcounts.py $1-gene_mapping_85perand50or100perand10.csv
mv rgibwtreadcounts.csv $1-gene_mapping_85perand50or100perand10_readcounts.csv
python additionalscripts/card_counts_genefam.py $1-gene_mapping_85perand50or100perand10.csv
mv genefamilies.csv $1-genefamilies-85perand50or100perand10.csv
mv genefamilyreadcounts.csv $1-genefamilyreadcounts-85perand50or100perand10.csv
