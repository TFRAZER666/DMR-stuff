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
pheatmap(memory_dmr_mCG[ ,c(1,4,5,10,8,15,14,12,13)],
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

# Attemped to open memory DMR file vis the new package
infile <- system.file("memory_dmr_mCG.csv",
                      package = 'bsseq')
bismarkBSseq <- read.csv(file = "memory_dmr_mCG.csv")

# Checked to see if data opened was correct
head(memory_dmr_mCG)

# Unsure if the output for this was appropriate as did not fully match what was seen in the guide
memory_dmr_mCG <- BSseq(chr = chr, pos = pos,
            M = M, Cov = Cov, 
            sampleNames = sampleNames)
show(memory_dmr_mCG)
?show         

# I do not understand what this code is trying to do and unable to adjust it to work on my dataset                
pData(memory_dmr_mCG)$CellType <- celltype
pData(memory_dmr_mCG)$Replicate <- substr(sampleNames, 
                              nchar(sampleNames), nchar(sampleNames))

pData(memory_dmr_mCG)
