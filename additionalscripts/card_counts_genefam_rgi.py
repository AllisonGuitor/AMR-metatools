import sys
import csv
input = open(sys.argv[1], "r").readlines()

##Takes the RGI main output.txt results and counts the number of ARGs in each
## AMR gene family

sequences = {}
for l in input:
    l = l.strip().split("\t")
    if not l[0].startswith("ORF"):
        fam = l[16]
        if (fam in sequences):
            sequences[fam]["counts"] = sequences[fam]["counts"] + 1
            sequences[fam]["names"] = sequences[fam]["names"] + "," + str(l[8])
        else:
            sequences.update({
              fam: {
                "counts": 1,
                "names": "",
              }})
            sequences[fam]["names"] = str(l[8])

with open('genefamilies.csv', 'w') as csvfile:
    body = []
    count = 0
    header = "AMR Gene Family, Number of Genes, Gene Names"
    csvfile.write(header + '\n')
    writer = csv.writer(csvfile)
    for item in sequences:
        counts = sequences[item]["counts"]
        names = sequences[item]["names"]
        body.append([item, counts, names])
        count += 1
    print ("AMR GENE FAMILY COUNTS:" + str(count))
    for item in body:
        writer.writerow(item)

with open('genefamilyreadcounts.csv', 'w') as csvfile:
    body = []
    header = "AMR Gene Family, Reads"
    for item in sequences:
        counts = sequences[item]["counts"]
        body.append([item, counts])
    csvfile.write(header + '\n')
    writer = csv.writer(csvfile)
    for n in body:
        writer.writerow(n)
