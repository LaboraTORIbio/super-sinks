#! /usr/bin/python3

from bs4 import BeautifulSoup
import openpyxl
import pandas as pd
import os, argparse, sys
from functools import reduce
#pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option("display.max_rows", None)


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Parser for blastn results of CRISPR spacers against a sequence''')

parser.add_argument('DIRECTORY', action="store", help="Master directory containing the files")

args = parser.parse_args()

DirName = args.DIRECTORY
DirValidation = os.path.isdir(DirName)
if DirValidation not in ['True', True]:
	print("Error. Please, enter a valid master directory to parse. Try --help to see the parser's description. Exiting!")
	exit(1)
DirectoryNames = os.listdir(DirName)

blastnres = DirName+'/blastn_spacers_pOXA-48_K8.tsv'
print("Filtering "+blastnres)

with open(blastnres, 'r') as blastn:	
	for line in blastn:
		line = line.rstrip("\n")
		if line.startswith('sp'):
			fields = line.split("\t")
			if int(fields[3])/int(fields[4]) > 0.60: # Filtering by % of alignment
				print(line)

