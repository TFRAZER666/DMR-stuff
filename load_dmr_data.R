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

