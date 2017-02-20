#!/usr/bin/python

# python get_annotation.py <an annotation_table.txt> <corresponding gene list>
# match your up/down regulated genes with the correspoding annotation
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

import sys, os

annotation_table=open(str(sys.argv[1]),"r")
de_genes=open(str(sys.argv[2]),"r")
outfile=open(str(sys.argv[2]).split(".txt")[0]+"_annotated.txt","a")

adict={}
for at in annotation_table:
	gene=at.split("\t")[0]
	pref=gene.split(".")[0]
	suff=gene.split(".")[1]
	padded=suff.strip().zfill(4)
	a=pref.strip()+"."+padded.strip()
	ann=at.split("\t")[1]
	adict[a]=ann.strip()
annotation_table.close()

for de in de_genes:
	try:
		peg=de.split("\t")[1]
		peg=peg.replace('"',"")
		annot=adict[peg.strip()]
		print >> outfile, peg.strip()+"\t"+annot.strip()
	except IndexError:
		print "Skipping header..."
print "Annotations finished!"
outfile.close()
