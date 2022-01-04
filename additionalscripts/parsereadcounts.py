import sys
import csv
results = open(sys.argv[1]).readlines()

with open('parsedrgibwtresults-10.csv', 'w') as csvfile:
    body = []
    header = list()
    count = 0
    for i in results:
        if i.startswith("ARO"):
            t = i.strip().split('\t')
            header = t
            body.append(header)
        if not i.startswith("ARO"):
            t = i.strip().split('\t')
            reads = float(t[11])
            if reads >= 10:
                count += 1
                subbody = t
                body.append(subbody)
    print ("ARO greater than 10:" + str(count))
    writer = csv.writer(csvfile)
    for n in body:
        writer.writerow(n)

with open('parsedrgibwtresults-100.csv', 'w') as csvfile:
    body = []
    header = list()
    count = 0
    for i in results:
        if i.startswith("ARO"):
            t = i.strip().split('\t')
            header = t
            body.append(header)
        if not i.startswith("ARO"):
            t = i.strip().split('\t')
            reads = float(t[11])
            if reads >= 100:
                count += 1
                subbody = t
                body.append(subbody)
    print ("ARO greater than 100:" + str(count))
    writer = csv.writer(csvfile)
    for n in body:
        writer.writerow(n)
