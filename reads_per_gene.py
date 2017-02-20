#!/usr/bin/python

# python reads_per_gene.py <gff file> <combined.cov matrix>
# count number of reads within gene interval
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll lmc297@cornell.edu

import sys, os, collections

gff=open(str(sys.argv[1]),"r")

gr={}
gstrand={}
for g in gff:
	gsplits=g.split("\t")
	s=gsplits[0]
	e=gsplits[1]
	r=range(int(s.strip()),int(e.strip())+1,1)
	strand=gsplits[2]
	fid=gsplits[3]
	fidpre=fid.split(".")[0]
	fidsuff=fid.split(".")[1]
	newsuff=fidsuff.strip().zfill(4)
	uid=fidpre.strip()+"."+newsuff.strip()
	gr[uid.strip()]=r
	gstrand[uid.strip()]=strand.strip()
gff.close()

outfile=open("final_coverage_table.txt","a")

covheader=open(str(sys.argv[2]),"r")
line1=covheader.readline()
nsamp=line1.strip().count("\t")+1
sampdict={}
for n in range(0,nsamp,1):
	sampdict[n]=0
print >> outfile, "geneID"+"\t"+"Start"+"\t"+"End"+"\t"+"Strand"+"\t"+line1.strip()
covheader.close()

incov=open(str(sys.argv[2]),"r")
cov=incov.readlines()

gr=collections.OrderedDict(sorted(gr.items()))

for g in gr.keys():
	generange=gr[g]
	sampdict=dict.fromkeys(sampdict,0)
	for pos in generange:
		cline=cov[pos]
		dsplits=cline.split("\t")
		for n,d in zip(sampdict.keys(),dsplits):
			reads=sampdict[n]
			sampdict[n]=reads+int(d.strip())
	finalline=[g.strip(),gr[g.strip()][0],gr[g.strip()][-1],gstrand[g.strip()]]
	for v in sampdict.values():
		finalline.append(v)
	#finalline.append([v for v in sampdict.values()])
	print >> outfile, "\t".join([str(f).strip() for f in finalline])
outfile.close()		
				
