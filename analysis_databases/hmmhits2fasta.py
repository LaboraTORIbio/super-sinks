#!/usr/bin/python3

import os, argparse, sys
import glob
import subprocess

PFfiles = glob.glob("PF*txt")

# Iterate through PF files
for file in PFfiles:
	outfile = file[:-4]
	outfile = "refseq_hmm/"+outfile+".faa"
	print("Generating file "+outfile)
	count=0
	
	with open(file, "r") as start_file:
		start_patterns = start_file.readlines()
		# Iterate through the start patterns
		for start_pattern in start_patterns:
			count+=1
			start_pattern = start_pattern.rstrip()
			if start_pattern.startswith(("ESCO", "KLPN")):
				name = start_pattern[:-11]
				name = "refseq/"+name+".faa"
				subprocess.call("echo \">\"" + start_pattern + " >> " + outfile, shell=True)
				#print("sed -n '/"+start_pattern+"/,/>/{/>/!p}' " + name + " >> " + outfile)
				subprocess.call("sed -n '/"+start_pattern+"/,/>/{/>/!p}' " + name + " >> " + outfile, shell=True)
			else:
				subprocess.call("echo \">\"" + start_pattern + " >> " + outfile, shell=True)
				#print("sed -n '/"+start_pattern+"/,/>/{/>/!p}' refseq/GCF_003204095.1.faa >> " + outfile)
				subprocess.call("sed -n '/"+start_pattern+"/,/>/{/>/!p}' refseq/GCF_003204095.1.faa >> " + outfile, shell=True)
	print("  Sequences added: " + str(count))
