#!/bin/bash

# sh bwa_mem_loop.sh <path to reference.fasta> <path to directory with .fastq.gz reads>
# a simple loop to create bwa indexes and run bwa mem for single-end reads
# example for RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

# first, let's index our reference
bwa index $1

# now, let's move to our directory with .fastq.gz files
cd $2

# let's loop through all of our files in this directory with a .fastq.gz extension
for f in *.fastq.gz
do
# run bwa mem and save our output as a sam file with a name that matches our fastq.gz file
bwa mem $1 $f > ${f%.fastq.gz}.sam
done
