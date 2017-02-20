#!/usr/bin/python

# python rast2gff.py <rast gff3 file>
# Converts rast gff3 formal to gff for RNA-Seq
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

import sys, os

infile=open(str(sys.argv[1]),"r")
outfile=open(str(sys.argv[1]).split(".gff")[0]+"_final.gff","a")
outpeg=open(str(sys.argv[1]).split(".gff")[0]+"_annotation_table.txt","a")

testpeg=[]
for line in infile:
	if "##" not in line:
		splits=line.split("\t")
		s=splits[3]
		e=splits[4]
		strand=splits[6]
		description=splits[8]
		if "peg." in description:
			gid="peg."
		elif "rna." in description:
			gid="rna."
		suff=description.split(gid)[1]
		g=suff.split(";Name=")[1]
		pegnum=suff.split(";Name=")[0]
		if gid+pegnum.strip() not in testpeg:
			print >> outpeg, gid+pegnum.strip()+"\t"+g.strip()
			print >> outfile, s.strip()+"\t"+e.strip()+"\t"+strand.strip()+"\t"+gid+pegnum.strip()
			testpeg.append(gid+pegnum.strip())
infile.close()
outfile.close()
outpeg.close()
