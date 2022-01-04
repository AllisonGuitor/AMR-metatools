import os, sys, json, csv, argparse, subprocess
import pandas as pd

def check_top_hits(data, term):
	"""
	Get top hit
	"""
	top_hit = ""
	for i in data:
		if term in data[i]["filter_reads_for_genes"]:
			top_hit = data[i]["top_hit"]
	return top_hit

def get_genes_from_node(f,node, term):
	cmd = "cat {path_to_file} | grep {node} | cut -f9 ".format(
		path_to_file=f,
		node=node,
		term=term)

	output = subprocess.check_output(cmd, shell=True)

	l = output.decode().rstrip().split("\n")

	if term in l:
		l.remove(term)

	return l


def get_gene_info(f, accession, term):
	hits = []
	all_terms = {}
	# hit = {
	# 	"accession":accession,
	# 	"term": term,
	# 	"ORF_ID": "",
	# 	"Cut_Off": "",
	# 	"Best_Hit_Bitscore": "", 
	# 	"Best_Hit_ARO": "",
	# 	"Best_Identities": "",
	# 	"Percentage Length of Reference Sequence": "",
	# 	"Other genes on same node": [],
	# 	"Coverage": 0.0,
	# 	"Length": 0,
	# }
	with open(os.path.join(f), 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			# print(row)
			if "ORF_ID" not in row[0]:
				arr = row[1].split("_")
				all_terms.update(
				{ row[8]: {
					"accession": row[10],
					"Cut_Off": row[5],
					"Percentage Length of Reference Sequence": row[20],
					"Coverage": float(arr[5]),
					"Length": int(arr[3]),
					"Drug Class": row[14],
					"Resistance Mechanism":  row[15],
					"AMR Gene Family": row[16]
					# "Other genes on same node": get_genes_from_node(f,"{}_{}_".format(arr[0],arr[1]), term)
				}
				}
				)
				if accession == row[10]:
					hit = {}
					
					hit["accession"] = accession
					hit["term"] = term
					hit["ORF_ID"] = row[1]
					hit["Cut_Off"] = row[5]
					hit["Best_Hit_Bitscore"] = row[7]
					hit["Best_Hit_ARO"] = row[8]
					hit["Best_Identities"] = row[9]
					hit["Percentage Length of Reference Sequence"] = row[20]
					hit["Other genes on same node"] = get_genes_from_node(f,"{}_{}_".format(arr[0],arr[1]), term)
					hit["Coverage"] = float(arr[5])
					hit["Length"] = int(arr[3])
					hits.append(hit)
					# print(hits)
	# print(hits)
	return [hits, all_terms]

def main(args):
	try:
		with open(os.path.join(args.assembly_report), 'r') as jfile:
			j = json.load(jfile)
	except Exception as e:
		logger.error(e)
		exit()
	# join raw clinical report to rgi main assembly report
	raw_clinical_report = pd.read_csv(os.path.join(args.raw_clinical_report),keep_default_na=False, na_values=[""],sep='\t')
	rgi_report = pd.read_csv(os.path.join(args.rgi_report),keep_default_na=False, na_values=[""],sep='\t')
	merged = pd.merge(left=raw_clinical_report,right=rgi_report, how='left', left_on='ARO Accession', right_on='ARO', sort=True)

	# out_file = os.path.join("{}.output.tsv".format(args.output))
	out_file = os.path.join("{}.final_output.tsv".format(args.output))
	# if reads > 0:
	# 	merged = merged[merged["All Mapped Reads"] > reads]
	# if coverage:
	# 	merged = merged[merged["Average Percent Coverage"] >= coverage]

	# remove nudged
	merged = merged[merged["Nudged"] != True]
	# remove hits with no rgi
	merged = merged[~pd.isnull(merged["ARO"])]

	contig_info = {}
	# print(merged["Contig"].values)
	# print(merged["ARO Term"].values)
	terms = []
	for term in merged["ARO Term"].values:
		terms.append(term)

	for i in merged["Contig"].values:
		arr = i.split("_")
		term = terms.pop()
		# print(term, "=>", get_genes_from_node(args.rgi_report,"{}_{}_".format(arr[0],arr[1]), term))
		contig_info[term] = {"same_node": get_genes_from_node(args.rgi_report,"{}_{}_".format(arr[0],arr[1]), term),
		"length": arr[3],
		"coverage": arr[5],
		"top_hit": check_top_hits(j, term)
		}
	# print(merged["Contig"].split("_"))
	# print(len(contig_info.keys()))
	same_node = []
	coverage = []
	length = []
	top_hit = []
	# merged = merged.iloc[::-1]
	# modified = merged.assign(What = [])
	for index, row in merged.iloc[::-1].iterrows():
		# print(row["ARO Term"], row["Contig"])
		# print(contig_info[row["ARO Term"]])
		same_node.append('; '.join(contig_info[row["ARO Term"]]["same_node"]))
		coverage.append(contig_info[row["ARO Term"]]["coverage"])
		length.append(contig_info[row["ARO Term"]]["length"])
		top_hit.append(contig_info[row["ARO Term"]]["top_hit"])

	top_hit.reverse()

	merged["Same Node"] = same_node
	merged["Coverage"] = coverage
	merged["Length"] = length
	merged["top_hit"] = top_hit

	merged.to_csv(os.path.join(out_file),sep="\t", index=False)

	# annotate report with top hits from selected reads for assembly (separate column?)

def _main(args):
	hits_in_bwt_and_kmers = []
	hits_in_rgi = []
	try:
		with open(os.path.join(args.assembly_report), 'r') as jfile:
			j = json.load(jfile)
	except Exception as e:
		logger.error(e)
		exit()
	other = {}
	final = []
	with open(os.path.join(args.raw_clinical_report), 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			if "ARO Term" not in row[0] and row[0]:
				# print(">>>>", row[0])
				
				res = check_top_hits(j, row[0])

				if res != row[0] and res != "":
					# skip this hit
					# pass
					# print("??", row[0]) 
					other[row[0]] = row

				else:
					# print(row[0],"??")
					values = get_bwt_values(args.bwt, row[1])
					hits_in_bwt_and_kmers.append(row[0])
					info, a = get_gene_info(args.rgi_report, row[1], row[0])
					# if res == "":
						# print("[name]=>", row[0], "\n[res]=>", res, "\n[row]=>", row, "\n[info]=>",info,"\n")
						# print(info)
					hits_in_rgi = a
					if len(info) > 0:
						# print("[name]=>", row[0], "\n[res]=>", res, "\n[row]=>", row, "\n[info]=>",info,"\n")
						for hit in info:
							# print("??>>>>>", hit)
							# print("")
							observed_in_rgi_main = "NO"
							if hit["Best_Hit_ARO"] != "":
								observed_in_rgi_main = "YES"								
							# print(row)
							# print("[if]", row[0])
							# print("[name]=>", row[0])
							final.append([
								row[0], 
								row[1],
								hit["Cut_Off"], # Cut_Off
								observed_in_rgi_main, # Observed in RGI main
								row[11],
								row[12],
								values[12],
								row[2],
								row[3],
								row[4],
								row[5],
								row[6],
								row[7], 
								row[8],
								row[9],
								row[10],
								'; '.join(hit["Other genes on same node"]),
								hit["Coverage"],
								hit["Length"]
							])
					else:
						# print("[else]", row[0])
						# print("[name]=>", row[0], "\n[res]=>", res, "\n[row]=>", row, "\n[info]=>",info,"\n")
						# print("[name]=>", row[0])
						final.append([
							row[0], 
							row[1],
							"N/A",
							"NO",
							row[11],
							row[12],
							"N/A",
							row[2],
							row[3],
							row[4],
							row[5],
							row[6],
							row[7], 
							row[8],
							row[9],
							row[10],
							"",
							"N/A",
							"N/A"
						])

	#"""
	# hits in bwt and kmers but no hits for rgi main assembly
	missing_hits_assembly = (diff(hits_in_bwt_and_kmers, hits_in_rgi.keys()))

	# hits missing from kmers or bwt and rgi main assembly
	missing_hits_kmers = (diff(hits_in_rgi.keys(), hits_in_bwt_and_kmers))

	for i in missing_hits_kmers:
		# if i == "PmrF":
		# 	print("[what]",i)
		# 	print(other.keys())
			# print(other[i])
		values = get_bwt_values(args.bwt, hits_in_rgi[i]["accession"])
		if i in other.keys():
			final.append([
				i, 
				hits_in_rgi[i]["accession"],
				hits_in_rgi[i]["Cut_Off"],
				"YES",
				"YES",
				"NO",
				hits_in_rgi[i]["Percentage Length of Reference Sequence"], # other[i][2]???
				values[12],
				other[i][3],
				other[i][4],
				other[i][5],
				"N/A",
				hits_in_rgi[i]["Drug Class"],
				hits_in_rgi[i]["Resistance Mechanism"],
				hits_in_rgi[i]["AMR Gene Family"],
				"N/A",
				"",
				hits_in_rgi[i]["Coverage"],
				hits_in_rgi[i]["Length"]
			])
		else:
			# in bwt, kmers and but not in partial assembly
			# values = get_bwt_values(args.bwt, hits_in_rgi[i]["accession"])
			# print(">>>", i)
			if values is not None:
				final.append([
					i, 
					hits_in_rgi[i]["accession"],
					hits_in_rgi[i]["Cut_Off"],
					"YES",
					"YES",
					"NO",
					hits_in_rgi[i]["Percentage Length of Reference Sequence"], # other[i][2]???
					values[12],
					values[11],
					0,
					"",
					values[8],
					hits_in_rgi[i]["Drug Class"],
					hits_in_rgi[i]["Resistance Mechanism"],
					hits_in_rgi[i]["AMR Gene Family"],
					"N/A",
					"",
					hits_in_rgi[i]["Coverage"],
					hits_in_rgi[i]["Length"]
				])
			else:
				# only in rgi main assembly
				# print(hits_in_rgi[i])
				final.append([
					i, 
					hits_in_rgi[i]["accession"],
					hits_in_rgi[i]["Cut_Off"],
					"YES",
					"YES",
					"NO",
					hits_in_rgi[i]["Percentage Length of Reference Sequence"], # other[i][2]???
					"N/A",
					"N/A",
					0,
					"",
					"N/A",
					hits_in_rgi[i]["Drug Class"],
					hits_in_rgi[i]["Resistance Mechanism"],
					hits_in_rgi[i]["AMR Gene Family"],
					"N/A",
					"",
					hits_in_rgi[i]["Coverage"],
					hits_in_rgi[i]["Length"]
				])
			# exit
	#"""


	with open(os.path.join("{}.final_output.tsv".format(args.output)), "w") as out_tab:
		writer = csv.writer(out_tab, delimiter='\t', dialect='excel')
		writer.writerow([
			"ARO Term", 
			"ARO Accession",
			"RGI Main: Cut_Off", # RGI Detection Paradigm (Perfect, Strict, Loose)
			"RGI Main: Observed in Denovo Assembly",
			"RGI BWT: Observed in Read Alignment",
			"RGI Pathogen ID: Observed in Kmer Prediction",
			"RGI Main: Percentage Length of Reference Sequence",
			"RGI BWT: Average Percent Coverage",
			"RGI BWT: All Mapped Reads",
			"RGI Pathogen ID: Mapped reads with kmer DB hits",
			"RGI Pathogen ID: CARD*kmer Prediction",
			"RGI BWT: Observed Pathogen(s)",
			# "NCBI BLAST top hit", # blast groot results to NCBI
			"AMR Gene Family", 
			"Drug Class", 
			"Resistance Mechanism",
			'RGI BWT: Mapped Reads with Flanking Sequence',
			"RGI Main: Other genes on same node", # from full assembly
			"RGI Main: Assembly Coverage", # from full assembly
			"RGI Main: Assembly Gene length" # from full assembly
		])
		for item in final:
			# print(item)
			writer.writerow(item)
	"""
	with open(os.path.join("{}.final_output.tsv".format(args.output)), "w") as out_tab:
		writer = csv.writer(out_tab, delimiter='\t', dialect='excel')
		writer.writerow([
			"ARO Term", 
			"ARO Accession",
			"Cut_Off", # RGI Detection Paradigm (Perfect, Strict, Loose)
			"Observed in RGI main",
			"Observed in RGI bwt",
			"Observed in RGI kmer",
			"Average Percent Coverage",
			"All Mapped Reads",
			"Mapped reads with kmer DB hits",
			"CARD*kmer Prediction",
			"Resistomes & Variants: Observed Pathogen(s)",
			# "NCBI BLAST top hit", # blast groot results to NCBI
			"AMR Gene Family", 
			"Drug Class", 
			"Resistance Mechanism",
			'Mapped Reads with Flanking Sequence',
			"Other genes on same node", # from full assembly
			"Coverage", # from full assembly
			"Gene length" # from full assembly
		])

		with open(os.path.join(args.raw_clinical_report), 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t')
			for row in reader:
				if "ARO Term" not in row[0] and row[0]:
					# print(">>>>", row[0])
					
					res = check_top_hits(j, row[0])
					if res != row[0] and res != "":
						# skip this hit
						pass
						# print(row[0],"??")  
					else:
						# print(row[0],"??")
						hits_in_bwt_and_kmers.append(row[0])
						info, a = get_gene_info(args.rgi_report, row[1], row[0])
						hits_in_rgi = a
						for hit in info:
							# print(">>>>>", hit)
							# print("")
							observed_in_rgi_main = "NO"
							if hit["Best_Hit_ARO"] != "":
								observed_in_rgi_main = "YES"								
								
							writer.writerow([
								row[0], 
								row[1],
								hit["Cut_Off"], # Cut_Off
								observed_in_rgi_main, # Observed in RGI main
								row[11],
								row[12],
								row[2],
								row[3],
								row[4],
								row[5],
								row[6],
								row[7], 
								row[8],
								row[9],
								row[10],
								'; '.join(hit["Other genes on same node"]),
								hit["Coverage"],
								hit["Length"]
							])
	# print(hits_in_bwt_and_kmers)
	# print()
	# print(hits_in_rgi.keys())
	# hits in bwt and kmers but no hits for rgi main assembly
	print(diff(hits_in_bwt_and_kmers, hits_in_rgi.keys()))
	print()
	# hits missing from kmers but in bwt and rgi main assembly
	print(diff(hits_in_rgi.keys(), hits_in_bwt_and_kmers))
	"""

def get_bwt_values(f, term):
	with open(os.path.join(f), 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			if "ARO Term" not in row[0]:
				# if row[0].lower() == term.lower():
				if row[1] == term:
					return row

def diff(first, second):
        second = set(second)
        return [item for item in first if item not in second]

def create_parser():
    parser = argparse.ArgumentParser(prog="consolidate_clinical_report",description='finalise card clinical report')
    parser.add_argument('--raw_clinical_report', dest="raw_clinical_report", required=True, help="raw clinical report file (TABULAR)")
    parser.add_argument('--assembly_report', dest="assembly_report", required=True, help="assembly json file (JSON)")
    parser.add_argument('--rgi_report', dest="rgi_report", required=True, help="rgi main tabular file (TABULAR)")
    parser.add_argument('--bwt_report', dest="bwt", required=True, help="rgi bwt gene mapping file (TABULAR)")
    parser.add_argument('--output', dest="output", required=True, help="output file (TABULAR)")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()