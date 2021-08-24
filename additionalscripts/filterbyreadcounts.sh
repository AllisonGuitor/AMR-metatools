#!/bin/bash

echo sample $1
echo naming $2

python ~/Desktop/Programs/parsergibwtresultsreadcounts.py $1.gene_mapping_data.txt

mv parsedrgibwtresults-10.csv $1-gene_mapping_10.csv

mv parsedrgibwtresults-100.csv $1-gene_mapping_100.csv

python ~/Desktop/Programs/getreadcounts.py $1-gene_mapping_10.csv

mv rgibwtreadcounts.csv $1-gene_mapping_10_readcounts.csv

python ~/Desktop/Programs/getreadcounts.py $1-gene_mapping_100.csv

mv rgibwtreadcounts.csv $1-gene_mapping_100_readcounts.csv

python ~/Desktop/Programs/getpercentcoverage.py $1.gene_mapping_data.txt

mv rgibwtpercentcoverage.csv $1-percentcoverage.csv
