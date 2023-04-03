#!/usr/bin/python3

import openpyxl
import pandas as pd
import os, argparse, sys
pd.set_option("display.max_columns", None)


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Parser for CRISPRCasFinder results to find spacers
Provide the GFF file for extracting the spacer sequences in fasta format''')

parser.add_argument('FILE', action="store", help="GFF file containing CRISPRs")
args = parser.parse_args()
nodes = args.FILE


## SAVING THE DATAFRAME OF THE REBASE DATABASE

db = []

with open(nodes, 'r') as file:	
	for line in file:
		if "CRISPRspacer" in line:
			line = line.rstrip('\n')
			line = line.split("\t")
			long = line[8]
			start = line[3]
			end = line[4]
			length = int(end) - int(start) +1
			long = long.split(";")
			header1 = long[-1]
			header1 = header1.replace("ID=", "")
			header2 = long[-2]
			header2 = header2.replace("Parent=", "")
			seq = long[0]
			seq = seq.replace("sequence=", "")
			print(">"+header1+"\t"+str(length)+"\t"+header2)
			print(seq)
