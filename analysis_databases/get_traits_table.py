#!/usr/bin/python3

import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''This script takes three files as input:
i) File exported from R containing the ordered tree nodes
ii) File of strains that have a dependent intrait (e.g. presence of pOXA-48)
iii) File of strains that have an independent trait (e.g. presence of Klebsiella capsule or T6SS)
And prints a table with the strain name and presence/absence of the traits ordered as the tree nodes''')

parser.add_argument('-n', action="store", help="File of ordered tree nodes", dest="NODES", required=True)
parser.add_argument('-d', action="store", help="File of strains that have a dependent trait", dest="DEP", required=True)
parser.add_argument('-i', action="store", help="File of strains that have an independent trait", dest="INDEP", required=True)

args = parser.parse_args()

Nodes = args.NODES
Dep = args.DEP
Indep = args.INDEP


print("Strain\tIndependent\tDependent")


## SAVING IN A LIST THE NEW NAMES OF THE STRAINS THAT HAVE A DEPENDENT TRAIT

listdep = []

with open(Dep, 'r') as dep:
	try:
		for line in dep:
			line = line.rstrip()
			listdep.append(line)
			
	except:
		print("Error: could not process dependent traits file")
#print("There are "+str(len(listdep))+" strains with the dependent trait")


## SAVING IN A LIST THE NEW NAMES OF THE STRAINS THAT HAVE AN INDEPENDENT TRAIT

listindep = []

with open(Indep, 'r') as indep:
	try:
		for line in indep:
			line = line.rstrip()
			listindep.append(line)
			
	except:
		print("Error: could not process independent traits file")
#print("There are "+str(len(listindep))+" strains with the independent trait")


## LOOPING OVER NODE NAMES AND ASSINGNING PRESENCE/ABSENCE OF TRAITS

with open(Nodes, 'r') as nodes:
	try:
		for line in nodes:
			line = line.rstrip()
			if line in listindep and line in listdep:
				print(line + "\t1\t1")
			elif line in listindep and line not in listdep:
				print(line + "\t1\t0")
			elif line not in listindep and line in listdep:
				print(line + "\t0\t1")
			else:
				print(line + "\t0\t0")
	except:
		print("Error: could not process tree nodes file")
