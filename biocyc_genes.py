#!/usr/local/python

# python biocyc_genes.py <de gene table with fold change> <table with gene ids>
# outputs table that can be imported into biocyc with DE gene names in the first column, and fold change in the second
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu 


import sys, os

fc=open(str(sys.argv[1]),"r")
gn=open(str(sys.argv[2]),"r")
outfile=open("biocyc_final.txt","a")

gdict={}
for g in gn:
	splits=g.split("\t")
	peg=splits[0]
	genename=splits[1]
	if peg.strip() not in gdict.keys():
		gdict[peg.strip()]=[]
		gdict[peg.strip()].append(genename.strip())
	else:
		gdict[peg.strip()].append(genename.strip())
gn.close()

fdict={}
for f in fc:
	if "geneID" not in f:
		splits=f.split("\t")
		peg=splits[1]
		peg=peg.replace('"',"")
		fold=splits[2]
		fdict[peg.strip()]=fold.strip()
fc.close()

for k in fdict.keys():
	try:
		realgene=gdict[k]
	except KeyError:
		realgene=["NA"]
	for rg in realgene:
		fc=fdict[k]
		print >> outfile, rg+"\t"+fc
outfile.close()
