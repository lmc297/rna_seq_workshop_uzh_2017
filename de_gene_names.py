#!/usr/bin/python

# python de_gene_names.py <de_gene_list_annotated.txt> <protein_sequences.faa> <your reference db to blast against>
# get the sequences of proteins so we can run blast and get gene names
# to get de_gene_list_annotated.txt run the following command:
# cat upregulated_gene_list_annotated.txt downregulated_gene_list_annotated.txt > de_gene_list_annotated.txt
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

import sys, os
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

degenes=open(str(sys.argv[1]),"r")
infasta=open(str(sys.argv[2]),"r")
outfile=open("de_gene_seqs.fasta","a")

genelist=[]
for de in degenes:
	peg=de.split("\t")[0]
	genelist.append(peg.strip())
degenes.close()

for record in SeqIO.parse(infasta,"fasta"):
	seqid=str(record.id).strip()
	seqseq=str(record.seq).strip()
	pref=seqid.split(".")[2]
	suff=seqid.split(".")[-1]
	pad=suff.strip().zfill(4)
	testid=pref.strip()+"."+pad.strip()
	if testid in genelist:
		print >> outfile, ">"+testid
		print >> outfile, seqseq
infasta.close()
outfile.close()

os.system("makeblastdb -in "+str(sys.argv[3])+" -dbtype prot")
cline=NcbiblastpCommandline(query="de_gene_seqs.fasta", db=str(sys.argv[3]), out="blastp_results.xml",outfmt='"5"')
stdout, stderr = cline()
result_handle=open("blastp_results.xml")
blast_records=NCBIXML.parse(result_handle)
filtres=open("blastp_results_table.txt","a")
headerlist=["hit_number","alignment_title","query_id","db_gene_query","alignment_length","evalue","blast_bitscore","query_coverage","percent_idenity","db_gene_start","db_gene_end","genome_sequence_start","genome_sequence_end","db_gene_sequence","genome_sequence","match_sequence"]
print >> filtres, "\t".join([h for h in headerlist])
# loop through each record in the blast output
counter=0
evalue_thresh=1e-5
qcov_thresh=0
pident_thresh=0
for record in blast_records:
	for alignment in record.alignments:
	# for each hsp in the alignment
		for hsp in alignment.hsps:
			qcov=float(hsp.align_length)/float(record.query_length)*100
			pident=float(hsp.identities)/float(hsp.align_length)*100
			# add information for each HSP if the e-value is lower than the threshold
			if hsp.expect<=evalue_thresh and qcov>=qcov_thresh and pident>=pident_thresh:
				counter=counter+1
                                maxseq=hsp.sbjct
                                maxgen=str(record.query).split("_")[0]
                                maxallele=record.query
                                genefacts=[counter,alignment.title,record.query_id,record.query,hsp.align_length,hsp.expect,hsp.bits,qcov,pident,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,hsp.query,hsp.sbjct,hsp.match]
				print >> filtres, "\t".join([str(g).strip() for g in genefacts])
filtres.close()
newin=open("blastp_results_table.txt","r")
biocyc=open("biocyc_table.txt","a")
for line in newin:
	if "hit_number" not in line:
		splits=line.split("\t")
		geneline=splits[1]
		gene=geneline.split("GN=")[1]
		gene=gene.split(" ")[0]
		peg=splits[3]
		print >> biocyc, peg.strip()+"\t"+gene.strip()
newin.close()
biocyc.close()
