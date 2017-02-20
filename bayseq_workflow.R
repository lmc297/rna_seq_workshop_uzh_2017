# baySeq Workflow for Detecting Differentially-Transcripted Genes
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

# Step 0: Set your working directory
# Replace the charcters between the quotation marks with the path to the directory where final_coverage_table.txt is stored
#setwd("/path/to/final_coverage_table.txt/directory/")

# Step 1: Install baySeq (if necessary)
#source("https://bioconductor.org/biocLite.R")
#biocLite("baySeq")

# Step 2: Load baySeq
library(baySeq)

# Step 3: Load our coverage table

infile<-read.delim("final_coverage_table.txt",header=TRUE,sep="\t")
# Maybe you only want to test protein-encoding genes (no rna). To do that, run:
infile<-subset(infile,!grepl("rna",infile$geneID))

# Step 4: Reshape our coverage table
# install the dplyr and stringr packages, if necessary
#install.packages("dplyr")
#install.packages("stringr")
# load the dplyr and stringr packages
library(dplyr)
library(stringr)
# First, we'll add a column with the length of each gene
infile$Length <- infile$End-infile$Start+1
# Reorder the columns, putting the Length column after End
ordered.table <- cbind(infile[,1:3],infile$Length,infile[,4:20])
# Rename Length column header
names(ordered.table)[4] <- "Length"
# Now let's select the sense samples and positive-strand genes
# First, let's create a subtable with the positive-stranded genes
strand.pos <- subset(ordered.table,ordered.table$Strand=="+")
# Now, let's get rid of the antisense columns by selecting only the sense columns
strand.pos <- select(strand.pos, 1:5, contains("_sense.cov"))
# Finally, let's change the names of the columns to get rid of the suffix
names(strand.pos) <- gsub("_sense.cov","",names(strand.pos))
# Now let's do the same thing for the anti-sense and minus-strand genes
# First, let's create a subtable with the minus-strand genes
strand.neg <- subset(ordered.table,ordered.table$Strand=="-")
# Now, let's get rid of the sense columns by selecting only the antisense columns
strand.neg <- select(strand.neg, 1:5, contains("_antisense.cov"))
# Finally, let's change the names of the columns to get rid of the suffix
names(strand.neg) <- gsub("_antisense.cov","",names(strand.neg))
# Now, let's merge the pos/sense and neg/antisense tables together
coverage.table <- rbind(strand.pos,strand.neg)
# Finally, let's sort by our geneIDs, and get rid of any genes that have zero mapped reads
ctable <- coverage.table[order(coverage.table$geneID),]
ctable <- ctable[rowSums(ctable[,6:ncol(ctable)])>0, ]


# Step 5: Use baySeq to Identify Differentially-Transcribed Genes
# First, let's store some subsets of our table as variables that we can use later
locus <- data.frame(ctable[,1])
lengths <- data.frame(ctable[,4])
coverage <- as.matrix(ctable[,6:ncol(ctable)])
# let's assign "1" to our control replicates (bhi) and "2" to our cold-smoked salmon replicates
# our first 4 columns are bhi, so we'll give them "1"
# our last 4 columns are css, so we'll give them "2"
reps <- c("1","1","1","1","2","2","2","2")
# now we are going to define our model
# we expect that many of the genes in our control and case will not show differences in transcript levels
# this is our "Not Differentially Expressed" (NDE) group, which we will code as all 1s
# we also expect (and hope!) that some genes will be differentially expressed in our case relative to our control
# this is our "Differentially Expressed" (DE) group, which we'll code as 1s for our control, and 2s for our case
groups <- list(NDE=c(1,1,1,1,1,1,1,1), DE=c(1,1,1,1,2,2,2,2))
# now we will create a "countData" object with the variables we've described above
CD <- new("countData", data = coverage, seglens = lengths, replicates = reps, groups = groups)
# we'll have baySeq caluculate our library sizes for us
libsizes(CD)<-getLibsizes(CD)
# remember the MA plots? let's make one of our raw data
par(mfrow=c(1,2))
plotMA.CD(CD, samplesA = "1", samplesB = "2", normaliseData = FALSE)
# we'll add a red line where M=0 (the log2 ratio = 0), just for fun
abline(0,0,col="red")
# now let's see how this changes when we normalize our data
plotMA.CD(CD, samplesA = "1", samplesB = "2", normaliseData = TRUE)
abline(0,0,col="red")
# our data definitely seem to be clustering closer to 0 after normalizing!
# now, let's add our gene names to our CD object
CD@annotation <- locus
# now, we're going to get some priors
# we're going to assume that our counts are generated from an underlying negative binomial distribution
# we're going to use baySeq to estimate the parameters of this distribution using quasi-likelihood methods
# the creators of baySeq recommend setting samplesize=10000, but here we'll use a sample size of 1000 to estimate the priors for the sake of time
CD<-getPriors.NB(CD, samplesize=1000, estimation = "QL", cl=NULL)
# whew! now that we have our priors, let's calculate posterior probabilities that each of our genes belong to either of our models
CD<-getLikelihoods(CD, cl=NULL, pET="BIC")
# now let's look at the genes that have the highest likelihood of being differentially expressed!
top_counts_DE<-topCounts(CD, group="DE", number=nrow(ctable), normaliseData=TRUE)
# we can plot the posterior likelihoods against the log-ratios for the two groups
par(mfrow=c(1,1))
plotPosteriors(CD, group = "DE")
abline(0.95,0,col="red")
# let's add a column containing fold change to our table
# we'll caluculate this by taking the average normalized counts for cold smoked salmon divided by the average normalized counts for bhi
top_counts_DE$FC <- (rowSums(top_counts_DE[,6:9])/4)/(rowSums(top_counts_DE[,2:5])/4)
# we'll define a gene as being upregulated if its FDR < 0.05 and fold change > 2.5
upregulated <- which(top_counts_DE$FDR.DE<0.05 & top_counts_DE$FC >= 2.5)
# we'll define a gene as being downregulatd if its FDR < 0.05 and its fold change < 0.4
downregulated <- which(top_counts_DE$FDR.DE<0.05 & top_counts_DE$FC <= 0.4)
# let's get the names of upregulated genes
uplist <- top_counts_DE$annotation[upregulated]
write.table(uplist,file="upregulated_gene_list.txt",sep="\t")
# names of downregulated genes
downlist <- top_counts_DE$annotation[downregulated]
write.table(downlist,file="downregulated_gene_list.txt",sep="\t")

# Step 6: Output table for use with goseq (later!)
# first create a vector of 1s and 0s
# 0 if a gene was not differentially expressed, 1 if it was
de.vec <- as.numeric(ctable$geneID %in% uplist|ctable$geneID %in% downlist)
DE_list <- data.frame(ctable$geneID,de.vec,ctable$Length)
names(DE_list)<-c("locus","DE","Length")
# we can write this to a file, if desired 
# write.table(DE_list,file="DE_list.txt",sep="\t")

# Step 7: Output table for use with BioCyc (later!)
biocyc_de <- data.frame(top_counts_DE$annotation,top_counts_DE$FC)
biocyc_de <- biocyc_de[biocyc_de$top_counts_DE.annotation%in%uplist|biocyc_de$top_counts_DE.annotation%in%downlist,]
names(biocyc_de) <- c("geneID","FC")
write.table(biocyc_de,file="biocyc_de.txt",sep="\t")
