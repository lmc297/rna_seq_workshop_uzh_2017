#!/usr/bin/python

# python gbk2gff.py <.gbk file>
# converts genbank (.gbk) file to a 4-column table
# column1 = some_gene_id
# column2 = gene start
# column3 = gene end
# column 4 = strand
# Laura Carroll, February 2017

import sys, os
from Bio import SeqIO

infile=open(str(sys.argv[1]),"r")
outfile=open(str(sys.argv[1]).split(".gbk")[0]+"_gbk_final.gff","a")
outpeg=open(str(sys.argv[1]).split(".gbk")[0]+"_gbk_annotation_table.txt","a")

proteins={}
for record in SeqIO.parse(infile,"genbank"):
	for feature in record.features:
		if feature.type=="CDS" or "RNA" in feature.type:
			f=feature.location
			xref=feature.qualifiers["db_xref"]
			product=feature.qualifiers["product"]
			for x,p in zip(xref,product):	
				try:
					if "peg." in x:
						gid="peg."
					elif "rna." in x:
						gid="rna."
					pin=x.split(gid)[1]
					pi=gid+pin.strip()
					if pi not in proteins.keys():
						proteins[pi]=[]
						proteins[pi].append(f.start)
						proteins[pi].append(f.end)
						proteins[pi].append(f.strand)
						proteins[pi].append(p.strip())
					else:
						print "duplicate"
				except NameError:
					if x not in proteins.keys():
						proteins[x]=[]
						proteins[x].append(f.start)
						proteins[x].append(f.end)
						proteins[x].append(f.strand)
						proteins[x].append(p.strip())
					else:
						print "duplicate" 
infile.close()

for key in sorted(proteins.keys(), key=lambda k: proteins[k][0]):
	line=proteins[key]
	s=line[0]
	e=line[1]
	strand=line[2]
	product=line[3]
	if strand==1:
		fs="+"
	elif strand==-1:
		fs="-"
	s=str(s).replace(">","").replace("<","")
	e=str(e).replace(">","").replace("<","")
	if s.strip()=="0":
		s="1"
	print >> outfile, "\t".join([str(s),str(e),fs,key])
	print >> outpeg, key+"\t"+product
outfile.close()	
outpeg.close()
