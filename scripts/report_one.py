import os, sys, json, csv, argparse

def main(args):
	cmd = "cat {filename} | cut -f{columns} > {output_file}".format(
		filename=os.path.join(args.final_report),
		columns="1,2,3,5,6,7,8,9,10,12,13,14,15,22,24,25,26,27,28,29,30,31,32,38,43,44,45,47,50,51,52,53,54,55,62,65,66,67,68,69,70",
		output_file=os.path.join("{}.final_report1.tsv".format(args.output))
	)
	os.system(cmd)

def create_parser():
    parser = argparse.ArgumentParser(prog="clinical_report",description='creates card clinical report')
    parser.add_argument('--final_report', dest="final_report", required=True, help="final report file (TABULAR)")
    parser.add_argument('--output', dest="output", required=True, help="output file (TABULAR)")

    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()
