import sys
import csv
results = open(sys.argv[1]).readlines()

with open('rgibwtpercentcoverage.csv', 'w') as csvfile:
    body = []
    header = "Name, Reads"
    count = 0 
    total = 0 
    for i in results: 
        if not i.startswith("ARO"): 
            t = i.strip().split('\t')
            percent = float(t[12])
            aro = t[1]
            body.append([aro, str(percent)])
            count += 1
            total += percent
    average = total / count         
    print ("Average coverage :" + str(average))
    csvfile.write(header + '\n')
    writer = csv.writer(csvfile)
    for n in body: 
        writer.writerow(n)
