import os, sys, json, csv, argparse
import pandas as pd
'''
# RGI pathogen ID columns

0 ARO term
1 Mapped reads with kmer DB hits
2 CARD*kmer Prediction
3 Single species (chromosome) reads
4 Single species (chromosome or plasmid) reads
5 Single species (plasmid) reads
6 Single species (no genomic info) reads
7 Single genus (chromosome) reads
8 Single genus (chromosome or plasmid) reads
9 Single genus (plasmid) reads
10 Single genus (no genomic info) reads
11 Promiscuous plasmid reads
12 Unknown taxonomy (chromosome) reads
13 Unknown taxonomy (chromosome or plasmid) reads
14 Unknown taxonomy (no genomic info) reads

'''
def get_kmer_predictions(term, f):
	"""
	Get kmer predictions for a given term
	"""
	kmer = {
	"Mapped reads with kmer DB hits": 0,
	"CARD*kmer Prediction": ""
	}
	with open(os.path.join(f), 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			if "ARO Term" not in row[0]:
				if term == row[0]:
					kmer["Mapped reads with kmer DB hits"] = int(row[1])
					kmer["CARD*kmer Prediction"] = row[2]
					return kmer

'''
# RGI BWT columns

 0 - 'ARO Term',
 1 - 'ARO Accession',
 2 - 'Reference Model Type',
 3 - 'Reference DB',
 4 - 'Alleles with Mapped Reads',
 5 - 'Reference Allele(s) Identity to CARD Reference Protein (%)',
 6 - 'Resistomes & Variants: Observed in Genome(s)',
 7 - 'Resistomes & Variants: Observed in Plasmid(s)',
 8 - 'Resistomes & Variants: Observed Pathogen(s)',
 9 - 'Completely Mapped Reads',
 10 - 'Mapped Reads with Flanking Sequence',
 11 - 'All Mapped Reads',
 12 - 'Average Percent Coverage',
 13 - 'Average Length Coverage (bp)',
 14 - 'Average MAPQ (Completely Mapped Reads)',
 15 - 'Number of Mapped Baits',
 16 - 'Number of Mapped Baits with Reads',
 17 - 'Average Number of reads per Bait',
 18 - 'Number of reads per Bait Coefficient of Variation (%)',
 19 - 'Number of reads mapping to baits and mapping to complete gene',
 20 - 'Number of reads mapping to baits and mapping to complete gene (%)',
 21 - 'Mate Pair Linkage (# reads)',
 22 - 'Reference Length',
 23 - 'AMR Gene Family',
 24 - 'Drug Class',
 25 - 'Resistance Mechanism'

'''

def get_top_pathogen(list_of_pathogens):
	"""
	Filter out top pathogen
	"""
	top = 0
	top_term = ""
	l = list_of_pathogens.split("; ")
	for i in l:
		if i:
			a = i.split(": ")
			if int(a[1]) > top:
				top = int(a[1])
				top_term = a[0]
	return top_term

def sort_by_amr_gene_family_pick_top_term(out_file, f):
    """
    Filter out top pathogen
    """
    data = {}

    with open(os.path.join(f), 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if "ARO Term" not in row[0]:
                cat = row[23].replace("'","").replace("(","").replace(")","")
                if cat not in data.keys():
                    data[cat] = []
                    data[cat].append(row[0])
                else:
                    data[cat].append(row[0])

    with open(os.path.join("{}.families.json".format(out_file)), "w") as outfile2:
        json.dump(data, outfile2)

def _sort_by_amr_gene_family_pick_top_term(out_file, f):
	"""
	Filter out top pathogen
	"""
	data = {}

	with open(os.path.join(f), 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			if "ARO Term" not in row[0]:
				cat = row[7].replace("'","").replace("(","").replace(")","")
				if cat not in data.keys():
					data[cat] = []
					data[cat].append(row[0])
				else:
					data[cat].append(row[0])

	with open(os.path.join("{}.families.json".format(out_file)), "w") as outfile2:
		json.dump(data, outfile2)


def main(args):
	"""
	Filter RGI-BWT by coverage and number of reads to an ARO term
	"""
	bwt_report = pd.read_csv(os.path.join(args.bwt),keep_default_na=False, na_values=[""],sep='\t')
	kmer_report = pd.read_csv(os.path.join(args.kmer),keep_default_na=False, na_values=[""],sep='\t')
	merged = pd.merge(left=bwt_report,right=kmer_report, how='left', left_on='ARO Term', right_on='ARO term', sort=True)
	coverage = float(args.coverage)
	reads = int(args.reads)
	out_file = os.path.join("{}.output.tsv".format(args.output))
	# if reads > 0:
	# 	merged = merged[merged["All Mapped Reads"] > reads]
	# if coverage:
	# 	merged = merged[merged["Average Percent Coverage"] >= coverage]
	merged.to_csv(os.path.join(out_file),sep="\t", index=False)

	sort_by_amr_gene_family_pick_top_term(args.output,out_file)

def _main(args):
	"""
	Filter RGI-BWT by coverage and number of reads to an ARO term
	"""
	hits = {}
	out_file = os.path.join("{}.output.tsv".format(args.output))
	with open(out_file, "w") as out_tab:
		writer = csv.writer(out_tab, delimiter='\t', dialect='excel')
		writer.writerow([
			"ARO Term",
			"ARO Accession",
			"Average Percent Coverage",
			"All Mapped Reads",
			"Mapped reads with kmer DB hits",
			"CARD*kmer Prediction",
			"Resistomes & Variants: Observed Pathogen(s)",
			"AMR Gene Family",
			"Drug Class",
			"Resistance Mechanism",
			"Mapped Reads with Flanking Sequence",
			"Observed in RGI bwt",
			"Observed in RGI kmer"
			])
		with open(os.path.join(args.bwt), 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t')
			for row in reader:
				if "ARO Term" not in row[0]:
					# print(row[0])
					if float(row[12]) >= float(args.coverage) and float(row[11]) >  float(args.reads):
							k = get_kmer_predictions(row[0],args.kmer)
							# if row[0] == "arlR":
								# print(row[0], k)
							if k is not None:
								if int(k["Mapped reads with kmer DB hits"]) > 0 and k["CARD*kmer Prediction"] != "":
									writer.writerow([
										row[0],
										row[1],
										row[12],
										row[11],
										k["Mapped reads with kmer DB hits"],
										get_top_pathogen(k["CARD*kmer Prediction"]),
										row[8],
										row[23],
										row[24],
										row[25],
										row[10],
										"YES",
										"YES"
									])
								else:
									# if row[0] == "arlR":
										# print(row[0], k)
										# print("whhat", )
									# pass
									writer.writerow([
										row[0],
										row[1],
										row[12],
										row[11],
										0,
										"",
										row[8],
										row[23],
										row[24],
										row[25],
										row[10],
										"YES",
										"YES"
									])
							else:
								# pass
								writer.writerow([
									row[0],
									row[1],
									row[12],
									row[11],
									0,
									"",
									row[8],
									row[23],
									row[24],
									row[25],
									row[10],
									"YES",
									"NO"
								])

	sort_by_amr_gene_family_pick_top_term(args.output,out_file)


def create_parser():
    parser = argparse.ArgumentParser(prog="clinical_report",description='creates card clinical report')
    parser.add_argument('--bwt_report', dest="bwt", required=True, help="rgi bwt gene mapping file (TABULAR)")
    parser.add_argument('--kmer_report', dest="kmer", required=True, help="rgi kmer report file (TABULAR)")

    parser.add_argument('--coverage', dest="coverage", type=int, default=50, help="filter by coverage (INT)")
    parser.add_argument('--reads', dest="reads", type=int, default=10, help="filter by number of reads (INT)")
    parser.add_argument('--output', dest="output", required=True, help="output file (TABULAR)")

    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
