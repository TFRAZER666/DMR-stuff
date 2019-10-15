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

