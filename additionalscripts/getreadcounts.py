import sys
import csv
results = open(sys.argv[1]).readlines()
        
with open('rgibwtreadcounts.csv', 'w') as csvfile:
    body = []
    header = "Name, Reads"
    for i in results: 
        if not i.startswith("ARO"): 
            t = i.strip().split(',')
            reads = float(t[11])
            aro = t[1]
            body.append([aro, str(reads)])   
    csvfile.write(header + '\n')
    writer = csv.writer(csvfile)
    for n in body: 
        writer.writerow(n)
