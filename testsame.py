#!/usr/bin/python

import sys, os

infile1=open(str(sys.argv[1]),"r")
infile2=open(str(sys.argv[2]),"r")

for o,t in zip(infile1,infile2):
	s1=o.split("\t")
	s2=t.split("\t")
	for x,y in zip(s1[1:],s2[1:]):
		if x.strip()!=y.strip():
			print x.strip()+"\t"+y.strip()
infile1.close()
infile2.close()
