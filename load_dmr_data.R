library(GenomicRanges)
library(magrittr)
library(pheatmap)

## Load the mCG levels for DMRs from CSV files
memory_dmr_mCG <- read.csv(file = "memory_dmr_mCG.csv", row.names = 1)
ips_dmr_mCG <- read.csv(file = "ips_dmr_mCG.csv", row.names = 1)
 
# Plot DMRs ---------------------------------------------------------------

# Heatmap. Note that the default settings are not best here. Play with them to explore. 
pheatmap(memory_dmr_mCG)

### Now make the heatmap for the ips DMRs

pheatmap(ips_dmr_mCG)

# Attempted to understand clusering method
?pheatmap

# Made heatmap using previously covered information 
pheatmap(ips_dmr_mCG, 
         clustering_distance_rows = "correlation",
         na.rm = TRUE)

# Did the same for memory DMR heatmap
pheatmap(memory_dmr_mCG, 
         clustering_distance_rows = "correlation",
         na.rm = TRUE)

# Encountered an error (this bogie has got some moves) This was the code that caused the error
pheatmap(memory_dmr_mCG, 
         clustering_distance_rows = "correlation",
         na.rm = TRUE)

# Erased all the irratating little letter values in dataset
memory_dmr_mCG <- memory_dmr_mCG[is.finite(rowSums(memory_dmr_mCG)),]

# Tried again without a need for .rm function
pheatmap(memory_dmr_mCG, 
         clustering_distance_rows = "correlation")

# I was wrong, .rm included
pheatmap(memory_dmr_mCG,
         clustering_distance_rows = "correlation",
         na.rm = TRUE)       

# Set up everything to build final heat map for memory DMRs
pheatmap(memory_dmr_mCG[ ,c(1,5,7,10,8,15,14,12,13)],
         clustering_distance_rows = "correlation",
         cluster_rows=F,
         cluster_cols=F,
         na.rm = TRUE) 

# Set up everything to build final heat map for iDMRs
pheatmap(ips_dmr_mCG[ ,c(1,5,7,10,8,15,14,12,13)],
         clustering_distance_rows = "correlation",
         cluster_rows=F,
         cluster_cols=F,
         na.rm = TRUE)   

### Attempted to process data further and hopefully get setup to identify DMRs in NSC samples using DMRs seen in iPSC samples
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dmrseq")
a

library(dmrseq)
?dmrseq
library(bsseq)

# Attemped to open iPSC DMR file vis the new package
infile <- system.file("ips_dmr_mCG",
                      package = 'bsseq')
bismarkBSseq <- read.csv(file = "ips_dmr_mCG")

# Checked to see if data opened was correct
head(ips_dmr_mCG)

# Unsure if the output for this was appropriate as did not fully match what was seen in the guide
ips_dmr_mCG <- BSseq(chr = chr, pos = pos,
            M = M, Cov = Cov, 
            sampleNames = sampleNames)
show(ips_dmr_mCG)
?show         

# I do not understand what this code is trying to do and unable to adjust it to work on my dataset                
pData(ips_dmr_mCG)$CellType <- celltype
pData(ips_dmr_mCG)$Replicate <- substr(sampleNames, 
                              nchar(sampleNames), nchar(sampleNames))


#PRACTICE
data("BS.chr21")
print(celltype)
?dmrseq

# Added column to data (unsure why this was needed?)
pData(BS.chr21)$CellType <- celltype
pData(BS.chr21)$Replicate <- substr('Rep', 
                              nchar(1), nchar(2))

pData(BS.chr21)

# Set testcovariate to be cell type
testCovariate <- 'CellType'

# Followed code exactly to filter no coverage values and only select cell types of interest
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(BS.chr21, type="Cov")==0) == 0)
sample.idx <- which(pData(BS.chr21)$CellType %in% c("imr90", "h1"))

bs.filtered <- BS.chr21[loci.idx, sample.idx]

# Now moving on to attempting identification of DMRs
testCovariate <- "CellType"
regions <- dmrseq(bs=BS.chr21[240001:260000,],
                  cutoff = 0.05,
                  testCovariate=testCovariate)

show(regions)
library("BiocParallel")


?dmrseq

# ?dmrseq says that when blocksize is TRUE (minimum consideration = 5000bp)
# all other parameters form the DMR analysis on large scale methyaltion blocks is specified

#input this code to run comparison of methylation on larger blocks
testCovariate <- "CellType"
blocks <- dmrseq(bs=BS.chr21[120001:125000,],
                 cutoff = 0.05,
                 testCovariate=testCovariate,
                 block = TRUE,
                 minInSpan = 500,
                 bpSpan = 5e4,
                 maxGapSmooth = 1e6,
                 maxGap = 5e3)

head(blocks)

# Number of significant regions 
sum(regions$qval < 0.05)

# Proportion of significant regions hypermethylated
sigRegions <- regions[regions$qval < 0.05,]
sum(sigRegions$stat > 0) / length(sigRegions)

# A bit tricky with the plotting
?plotDMRs

annoTrack <- getAnnot("hg18")
yes

plotDMRs(BS.chr21, regions=regions[1,], testCovariate="CellType",
         annoTrack="hg18")


plotDMRs(BS.chr21, regions=blocks[1,], testCovariate="CellType",
         annoTrack="hg18")

# extraction of mean methylation differences
rawDiff <- meanDiff(BS.chr21, dmrs=sigRegions, testCovariate="CellType")
str(rawDiff)
# This was the end of practice stuff







# RETRY ON MY DATA
ips_dmr_mCG <- read.csv(file = "ips_dmr_mCG.csv", row.names = 1)

head(ips_dmr_mCG)

# My issue is here, I want to get 'ips_dmr_mCG' into a format that can be opened by pData to resemble the practice data
# How should I go about putting this data into a format that can be read this way?
pData(ips_dmr_mCG)



