import sys
import csv
results = open(sys.argv[1]).readlines()

## Additional script to filter the RGI*bwt kma mapping results based on number of reads
## mapping and the percent length of coverage of the reference gene.
## Current filters are for genes with at least 85% coverage and 50 reads OR
## 100 percent coverage and 10 reads. These cut-offs can be modified below.

with open('parsedrgibwtresults-85perand50or100perand10.csv', 'w') as csvfile:
    body = []
    header = list()
    count = 0
    for i in results:
        if i.startswith("ARO"):
            t = i.strip().split('\t')
            header = t
            body.append(header)
        if not i.startswith("ARO"):
            t = i.replace(", ","_")
            t = t.replace(",", "_")
            t = t.strip().split('\t')
            percent = float(t[12])
            reads = float(t[11])
            if (percent >= 85 and reads >= 50) or (percent >=100 and reads >=10):
                count += 1
                subbody = t
                body.append(subbody)
    print ("ARO percent greater than 85 and 50 reads or 100 per and 10 reads:" + str(count))
    writer = csv.writer(csvfile)
    for n in body:
        writer.writerow(n)
