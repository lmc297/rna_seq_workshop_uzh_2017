# goseq Workflow for Performing Gene Ontology (GO) Analysis of RNA-Seq Data
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

# Step 0: Set your working directory
# Replace the charcters between the quotation marks with the path to the directory where DE_list.txt is stored
#setwd("/path/to/DE_list.txt/directory/")

# Step 1: Install goseq, if necessary
# source("https://bioconductor.org/biocLite.R")
# biocLite("goseq")

# Step 2: Load goseq
library(goseq)

# Step 3: Load DE_list.txt, if necessary, as well as out _goseq.annot file
# if you still have DE_list saved in your environment, you can skip this step
# DE_list<-read.delim("DE_list.txt",header=TRUE,sep="\t")
h7858_go_terms <- read.delim("l_mono_b2g_goseq.annot",header=FALSE,sep="\t")
names(h7858_go_terms) <- c("locus","Go_term")

# Step 4: Define some variables to use with goseq
h7858_lengths<- DE_list[,3]
h7858_locus<- DE_list[,1]
DE<- DE_list[,2]
names(DE)=h7858_locus

# Step 5: Calculate a probability weighting function for our genes
# Because some of our genes are longer than others, we need to obtain a weight for each gene based on its length
pwf_DE= nullp(DE, bias.data=h7858_lengths)
# don't be alarmed by the warning message!

# Step 6: Get enriched GO terms, and correct for multiple testing
statsDE= goseq(pwf_DE, gene2cat= h7858_go_terms)
Go_terms_enriched= statsDE$category[p.adjust(statsDE$over_represented_pvalue, method= "BH")<0.05]
