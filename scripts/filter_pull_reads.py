import os, sys, json, csv, argparse
from operator import itemgetter, attrgetter
import shutil, glob

def pull_reads_from_bam(item,genes, f, read_one, read_two, mapq):
	"""
	Filter reads based on AMR gene family, assemble and run rgi main
	"""
	s = "|".join(genes)
	cat = item.replace(" ","_").replace(";","_")
	#'''
	read_list = "read_list_{}.txt".format(cat)

	# filter by mapq
	mapq_filter = ""

	if mapq > 0:
		mapq_filter = " | awk '$3>={}'".format(mapq)

	directory_path = os.path.join(cat)

	if os.path.exists(directory_path) and os.path.isdir(directory_path):
		print(">>>exists, remove it {}".format(directory_path))
		remove_directory(directory_path)
	else:
		print(">>>make dir: {}".format(directory_path))
		# create directory
		os.system("mkdir {}".format(cat))

	# extract reads name using ARO names
	os.system("samtools view {input_sorted_bam} | cut -f1,3,5 | grep -E \"{genes}\" {mapq_filter} | cut -f1 > {output_dir}/{read_list}".format(
		input_sorted_bam=f,
		genes=s,
		read_list=read_list,
		output_dir=cat,
		mapq_filter=mapq_filter
	))

	# extract reads my reads names from fastqs
	os.system("seqtk subseq {read_one} {output_dir}/{read_list} > {output_dir}/{read_one_fastq}.fastq".format(
		read_one=read_one,
		output_dir=cat,
		read_list=read_list,
		read_one_fastq="R1"
	))
	os.system("seqtk subseq {read_two} {output_dir}/{read_list} > {output_dir}/{read_two_fastq}.fastq".format(
		read_two=read_two,
		output_dir=cat,
		read_list=read_list,
		read_two_fastq="R2"
	))

	# re-pair reads
	os.system("fastq_pair -v {output_dir}/{read_one} -2 {output_dir}/{read_two}".format(
		output_dir=cat,
		read_one="R1.fastq",
		read_two="R2.fastq"
	))

	# assemble properly paired reads
	os.system("metaspades.py -t 40 -m 64 -1 {output_dir}/{read_one}.paired.fq -2 {output_dir}/{read_two}.paired.fq -o {output_dir}/metaspades_assembly".format(
		output_dir=cat,
		read_one="R1.fastq",
		read_two="R2.fastq"
	))
	#  rgi run on scaffolds
	# rgi main -i scaffolds.fasta -o scaffolds.fasta.output --debug -a diamond --clean --exclude_nudge
	os.system("rgi main -i {output_dir}/metaspades_assembly/{input_file} -o {output_dir}/{input_file}.output --local --debug -a blast --clean --exclude_nudge".format(
		output_dir=cat,
		scaffolds=cat,
		input_file="scaffolds.fasta"
	))

	# pathogen id
	# TODO:: run pathogen id on RGI*MAIN results and compare to RGI*BWT
	#'''

	# get top hit
	tab_file = "{output_dir}/{input_file}.output.txt".format(output_dir=cat,input_file="scaffolds.fasta")
	cc = item.replace("'","").replace("(","").replace(")","").strip()
	data = {}
	ordered = []
	all_hits = []
	if os.path.isfile(tab_file) == True:
		with open(os.path.join(tab_file), 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t')
			for row in reader:
				if "Best_Hit_ARO" not in row[5]:
					c = row[16].replace("'","").replace("(","").replace(")","").strip()
					if cc == c and row[8] in genes:
						criteria = 50
						if row[5] == "Perfect":
							criteria = 100

						if cc not in data.keys():
							data[cc] = []
							data[cc].append((
								row[8], criteria, row[9], row[20]
							))
						else:
							data[cc].append((
								row[8], criteria, row[9], row[20]
							))
						# print("[",cc,"] ===> ",
						# 	# Cut_Off
						# 	row[5],
						# 	# Best_Hit_ARO
						# 	row[8],
						# 	# Best_Identities
						# 	row[9],
						# 	# Percentage Length of Reference Sequence
						# 	row[20])
		if data:
			ordered_all = sorted(data[cc], key=itemgetter(1,2,3), reverse=True)
			# all_hits.append(ordered_all)
			ordered = ordered_all.pop(0)
			remove_directory(directory_path)
			return {"top_hit": ordered[0], "message": "success", "filter_reads_for_genes": genes, "assembled_genes_from_filtered_reads": data[cc]}
		else:
			remove_directory(directory_path)
			return {"top_hit":"", "message": "no rgi hits", "filter_reads_for_genes": genes, "assembled_genes_from_filtered_reads": ordered}
	else:
		remove_directory(directory_path)
		return {"top_hit": "", "message": "no assembly", "filter_reads_for_genes": genes, "assembled_genes_from_filtered_reads": ordered}

def remove_directory(directory_path):
	if os.path.exists(directory_path) and os.path.isdir(directory_path):
		try:
			print("Remove dir: {}".format(directory_path))
			shutil.rmtree(directory_path)
		except OSError as e:
			print("{} - {}" % (e.filename, e.strerror))
			exit("STOPPED")

def main(args):
	report = {}
	try:
		with open(os.path.join(args.json_file), 'r') as jfile:
			j = json.load(jfile)
	except Exception as e:
		print(e)
		exit()

	# assembly output directoty
	output_assembly = "{}.assembly".format(args.output)
	# create a directory for full assembly
	os.system("mkdir {}".format(output_assembly))

	# # assemble raw reads for chromosome
	os.system("metaspades.py -t 40 -m 64 --debug -1 {read_one} -2 {read_two} -o {output_dir}/metaspades_assembly".format(
		output_dir=output_assembly,
		read_one=args.read_one,
		read_two=args.read_two
	))
	# copy scaffolds.fasta to current directory
	cmd = "cp {}/metaspades_assembly/scaffolds.fasta {}.scaffolds.fasta".format(output_assembly, args.output)
	print(cmd)
	os.system(cmd)

	# run rgi with partial orfs
	rgi_cmd = "rgi main -i {}.scaffolds.fasta -o {}.scaffolds.fasta.output --debug --local --clean --low_quality --alignment_tool BLAST".format(args.output,args.output)
	print(rgi_cmd)
	os.system(rgi_cmd)

	# kmer predictions for full assembly
	rgi_kmer_cmd="rgi kmer_query --input {rgi_json} --output {output} --rgi --debug --local --kmer_size {kmer_size}".format(
		rgi_json="{}.scaffolds.fasta.output.json".format(args.output),
		output=args.output,
		kmer_size=61
	)
	print(rgi_kmer_cmd)
	os.system(rgi_kmer_cmd)

	for item in j:
		if len(j[item]) > 1:
			results = pull_reads_from_bam(item,j[item],args.bam_file, args.read_one, args.read_two, args.mapq)
			report.update({item: results})

	# write json
	with open(os.path.join("{}.assembly.json".format(args.output)), "w") as af:
		af.write(json.dumps(report,sort_keys=True))

def create_parser():
    parser = argparse.ArgumentParser(prog="filter_pull_reads",description='creates card clinical report')
    parser.add_argument('--json_file', dest="json_file", required=True, help="a file containing genes of interest (JSON)")
    parser.add_argument('--bam_file', dest="bam_file", required=True, help="bam alignment file from rgi bwt (BAM)")
    parser.add_argument('--read_one', dest="read_one", required=True, help="fastq file with forward paired-end reads (FASTQ)")
    parser.add_argument('--read_two', dest="read_two", required=True, help="fastq file with reverse paired-end reads (FASTQ)")
    parser.add_argument('--mapq', dest="mapq", required=False, type=int, default=0, help="filter reads by mapq (INT)")
    parser.add_argument('--output', dest="output", required=True, help="output file (TABULAR)")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
