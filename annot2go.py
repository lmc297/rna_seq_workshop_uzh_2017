#!/usr/bin/python

# python annot2go.py <*.annot file output from blast2go>
# converts blast2go annot file to one we can use with goseq in R
# RNA-Seq Workshop, Unviersity of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

import sys, os

infile=open(str(sys.argv[1]),"r")
outfile=open(str(sys.argv[1]).split(".annot")[0]+"_goseq.annot","a")

for line in infile:
	splits=line.split("\t")
	peg=splits[0]
	pref=peg.split(".")[2]
	suff=peg.split(".")[-1]
	pad=suff.strip().zfill(4)
	finalpeg=pref.strip()+"."+pad.strip()
	go=splits[1]
	print >> outfile, finalpeg+"\t"+go.strip()
infile.close()
outfile.close()
