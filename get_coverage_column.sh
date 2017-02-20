#!/bin/bash

# sh get_coverage_column.sh <path to directory with sense/antisense bam files> <path to reference.fasta file>
# a simple loop to create 3-column coverage files from our sense and anti-sense bam files, and then create a single-column coverage file
# example for RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

# go to directory with sense and anti-sense bam files
cd $1

# loop through anti-sense files
for f in *.as.sorted.bam
do
# create 3-column coverage file for anti-sense reads
genomeCoverageBed -d -ibam $f -g $2 > ${f%.as.sorted.bam}_antisense_cov.txt
# create a 3-column coverage file for sense reads
genomeCoverageBed -d -ibam ${f%.as.sorted.bam}.s.sorted.bam -g $2 > ${f%.as.sorted.bam}_sense_cov.txt
# create 1-column coverage file from anti-sense coverage file
awk '{ print $3 }' ${f%.as.sorted.bam}_antisense_cov.txt > ${f%.as.sorted.bam}_antisense.cov
# create 1-column coverage file from sense coverage file
awk '{ print $3 }' ${f%.as.sorted.bam}_sense_cov.txt > ${f%.as.sorted.bam}_sense.cov
done
