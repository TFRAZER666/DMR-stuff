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




#some potentially useful commands 

ips_dmr_mCG <- read.csv(file = "ips_dmr_mCG.csv", row.names = 1)

head(ips_dmr_mCG)

colnames(ips_dmr_mCG)

show(ips_dmr_mCG)
myvars <- c("h9_WGBS", "RL417_P3_32F_Primed.rmdup")
newdata <- ips_dmr_mCG[ ,myvars]
show(newdata)

subset(newdata, (c("h9") - c("32p")) > 0.1)





# Selection for hypermethylated iDMRs
esc <- c(1)
ipsc <- c(5,7)
keep <- (memory_dmr_mCG[ ,esc]) > rowMeans(memory_dmr_mCG[ ,ipsc])
hyper_newdata <- memory_dmr_mCG[keep, ]
show(hyper_newdata)
Finaldata <- c("h9_WGBS", 
               "RL418_P10_32F_Primed.rmdup", 
               "RL703_all_merge", 
               "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged", 
               "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge",
               "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge", 
               "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge",
               "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge",
               "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge")
PLEASEdata <- hyper_newdata[ ,Finaldata]
show(PLEASEdata)
boxplot(PLEASEdata)
pheatmap(PLEASEdata, show_rownames = FALSE)

### GGplot the data
library(reshape2)
library(ggplot2)

hyper_newdata_melt <- melt(PLEASEdata)
#ERROR IS OK HERE ^^
ggplot(hyper_newdata_melt, aes(x=value)) +
  geom_histogram() + 
  facet_grid(variable~.) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(hyper_newdata_melt, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# USE THIS PLOT AS WELL ^^

head(hyper_newdata_melt)
table(hyper_newdata_melt$variable)




#SHOW SPECIFIC GRAPHS
#Here are the 38F commands 
keep_row <- c("h9_WGBS", "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge", "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")
keep_row <- c("h9_WGBS", "RL703_all_merge", "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")
#32F commands
keep_row <- c("h9_WGBS", "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge", "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")
keep_row <- c("h9_WGBS", "RL418_P10_32F_Primed.rmdup", "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")


sub_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_row, ]
ggplot(sub_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#32 NSCS
primed <- c("RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge")
NP <- c("RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primed]
NSCNP <- PLEASEdata[ ,NP]
t.test(NSCprimed, NSCNP, paired=TRUE)

#32 iPSCs
primedstem <- c("RL418_P10_32F_Primed.rmdup")
NPstem<- c("RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")
iPSCprimed <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
t.test(iPSCprimed, iPSCNP, paired=TRUE)

#38F NSCs
primed <- c("RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge")
NP <- c("RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primed]
NSCNP <- PLEASEdata[ ,NP]
t.test(NSCprimed, NSCNP, paired=TRUE)

#38F iPSCs
primedstem <- c("RL703_all_merge")
NPstem<- c("RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")
iPSCprimed <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
t.test(iPSCprimed, iPSCNP, paired=TRUE)







########## THE FINAL ANALYSIS
#HYPOMETHYLATED MEMORY DMRS
memory_dmr_mCG <- read.csv(file = "memory_dmr_mCG.csv", row.names = 1)
esc <- c(1)
ipsc <- c(5,7)
keep <- (memory_dmr_mCG[ ,esc]) > rowMeans(memory_dmr_mCG[ ,ipsc])
hyper_newdata <- memory_dmr_mCG[keep, ]
show(hyper_newdata)
Finaldata <- c("h9_WGBS", 
               "RL418_P10_32F_Primed.rmdup", 
               "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged", 
               "RL703_all_merge",
               "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge", 
               "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge",
               "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge", 
               "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge",
               "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")
               
PLEASEdata <- hyper_newdata[ ,Finaldata]
pheatmap(PLEASEdata, show_rownames = FALSE)


hyper_newdata_melt <- melt(PLEASEdata)
#ERROR IS OK HERE ^^
ggplot(hyper_newdata_melt, aes(x=value)) +
  geom_histogram() + 
  facet_grid(variable~.) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(hyper_newdata_melt, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# USE THIS PLOT AS WELL ^^

pheatmap(PLEASEdata[ ,c(1,2,3,4,5,6,7,8,9)],
         clustering_distance_rows = "correlation",
         cluster_rows=F,
         cluster_cols=F,
         na.rm = TRUE) 

### Comparisons drawn between 38ipsc
keep_38iPSCrow <- c("h9_WGBS", "RL703_all_merge", "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL703_all_merge")
NPstem<- c("RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")
iPSCprimed2 <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed2, iPSCNP, paired = TRUE)
?wilcox.test


### Comparisons drawn between 38NSCs
keep_38NSCrow <- c("h9_WGBS", "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge", "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge")
NPNSC<- c("RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP, paired = TRUE)


### Comparisons drawn between 32ipsc
keep_32iPSCrow <- c("h9_WGBS", "RL418_P10_32F_Primed.rmdup", "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL418_P10_32F_Primed.rmdup")
NPstem<- c("RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")
iPSCprimed <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed, iPSCNP, paired = TRUE)

### Comparisons drawn between 32NSCs
keep_32NSCrow <- c("h9_WGBS", "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge", "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge")
NPNSC<- c("RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP1 <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP1, paired = TRUE)








########## THE FINAL ANALYSIS
# HYPERMETHYLATION MEMORY DMRS
memory_dmr_mCG <- read.csv(file = "memory_dmr_mCG.csv", row.names = 1)
esc <- c(1)
ipsc <- c(5,7)
keep <- (memory_dmr_mCG[ ,esc]) < rowMeans(memory_dmr_mCG[ ,ipsc])
hyper_newdata <- memory_dmr_mCG[keep, ]
show(hyper_newdata)
Finaldata <- c("h9_WGBS", 
               "RL418_P10_32F_Primed.rmdup", 
               "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged", 
               "RL703_all_merge",
               "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge", 
               "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge",
               "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge", 
               "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge",
               "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")

PLEASEdata <- hyper_newdata[ ,Finaldata]
pheatmap(PLEASEdata, show_rownames = FALSE)


hyper_newdata_melt <- melt(PLEASEdata)
#ERROR IS OK HERE ^^
ggplot(hyper_newdata_melt, aes(x=value)) +
  geom_histogram() + 
  facet_grid(variable~.) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(hyper_newdata_melt, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# USE THIS PLOT AS WELL ^^

pheatmap(PLEASEdata[ ,c(1,2,3,4,5,6,7,8,9)],
         clustering_distance_rows = "correlation",
         cluster_rows=F,
         cluster_cols=F,
         na.rm = TRUE) 

### Comparisons drawn between 38ipsc
keep_38iPSCrow <- c("h9_WGBS", "RL703_all_merge", "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL703_all_merge")
NPstem<- c("RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")
iPSCprimed2 <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed2, iPSCNP, paired = TRUE)
?wilcox.test


### Comparisons drawn between 38NSCs
keep_38NSCrow <- c("h9_WGBS", "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge", "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge")
NPNSC<- c("RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP, paired = TRUE)


### Comparisons drawn between 32ipsc
keep_32iPSCrow <- c("h9_WGBS", "RL418_P10_32F_Primed.rmdup", "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL418_P10_32F_Primed.rmdup")
NPstem<- c("RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")
iPSCprimed <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed, iPSCNP, paired = TRUE)

### Comparisons drawn between 32NSCs
keep_32NSCrow <- c("h9_WGBS", "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge", "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge")
NPNSC<- c("RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP1 <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP1, paired = TRUE)






########## THE FINAL ANALYSIS
#HYPOMETHYLATED IDMRS
ips_dmr_mCG <- read.csv(file = "ips_dmr_mCG.csv", row.names = 1)
esc <- c(1)
ipsc <- c(5,7)
keep <- (ips_dmr_mCG[ ,esc]) > rowMeans(ips_dmr_mCG[ ,ipsc])
hyper_newdata <- ips_dmr_mCG[keep, ]
show(hyper_newdata)
Finaldata <- c("h9_WGBS", 
               "RL418_P10_32F_Primed.rmdup", 
               "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged", 
               "RL703_all_merge",
               "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge", 
               "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge",
               "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge", 
               "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge",
               "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")

PLEASEdata <- hyper_newdata[ ,Finaldata]
pheatmap(PLEASEdata, show_rownames = FALSE)


hyper_newdata_melt <- melt(PLEASEdata)
#ERROR IS OK HERE ^^
ggplot(hyper_newdata_melt, aes(x=value)) +
  geom_histogram() + 
  facet_grid(variable~.) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(hyper_newdata_melt, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# USE THIS PLOT AS WELL ^^

pheatmap(PLEASEdata[ ,c(1,2,3,4,5,6,7,8,9)],
         clustering_distance_rows = "correlation",
         cluster_rows=F,
         cluster_cols=F,
         na.rm = TRUE) 

### Comparisons drawn between 38ipsc
keep_38iPSCrow <- c("h9_WGBS", "RL703_all_merge", "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL703_all_merge")
NPstem<- c("RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")
iPSCprimed2 <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed2, iPSCNP, paired = TRUE)
?wilcox.test


### Comparisons drawn between 38NSCs
keep_38NSCrow <- c("h9_WGBS", "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge", "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge")
NPNSC<- c("RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP, paired = TRUE)


### Comparisons drawn between 32ipsc
keep_32iPSCrow <- c("h9_WGBS", "RL418_P10_32F_Primed.rmdup", "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL418_P10_32F_Primed.rmdup")
NPstem<- c("RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")
iPSCprimed <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed, iPSCNP, paired = TRUE)

### Comparisons drawn between 32NSCs
keep_32NSCrow <- c("h9_WGBS", "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge", "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge")
NPNSC<- c("RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP1 <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP1, paired = TRUE)









########## THE FINAL ANALYSIS
#HYPERMETHYLATED IDMRS
ips_dmr_mCG <- read.csv(file = "ips_dmr_mCG.csv", row.names = 1)
esc <- c(1)
ipsc <- c(5,7)
keep <- (ips_dmr_mCG[ ,esc]) < rowMeans(ips_dmr_mCG[ ,ipsc])
hyper_newdata <- ips_dmr_mCG[keep, ]
show(hyper_newdata)
Finaldata <- c("h9_WGBS", 
               "RL418_P10_32F_Primed.rmdup", 
               "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged", 
               "RL703_all_merge",
               "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge", 
               "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge",
               "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge", 
               "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge",
               "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")

PLEASEdata <- hyper_newdata[ ,Finaldata]
pheatmap(PLEASEdata, show_rownames = FALSE)


hyper_newdata_melt <- melt(PLEASEdata)
#ERROR IS OK HERE ^^
ggplot(hyper_newdata_melt, aes(x=value)) +
  geom_histogram() + 
  facet_grid(variable~.) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(hyper_newdata_melt, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# USE THIS PLOT AS WELL ^^

pheatmap(PLEASEdata[ ,c(1,2,3,4,5,6,7,8,9)],
         clustering_distance_rows = "correlation",
         cluster_rows=F,
         cluster_cols=F,
         na.rm = TRUE) 

### Comparisons drawn between 38ipsc
keep_38iPSCrow <- c("h9_WGBS", "RL703_all_merge", "RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL703_all_merge")
NPstem<- c("RL837_2017_09_28_P12_plus_10_38F_SmithR_to_E8_merge")
iPSCprimed2 <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed2, iPSCNP, paired = TRUE)
?wilcox.test


### Comparisons drawn between 38NSCs
keep_38NSCrow <- c("h9_WGBS", "RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge", "RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_38NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1812_2019_09_19_38F_NSCs_primed_p8_S3_all_lanes_merge")
NPNSC<- c("RL1811_2019_09_19_38F_NSCs_N2P_p8_S2_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP, paired = TRUE)


### Comparisons drawn between 32ipsc
keep_32iPSCrow <- c("h9_WGBS", "RL418_P10_32F_Primed.rmdup", "RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32iPSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedstem <- c("RL418_P10_32F_Primed.rmdup")
NPstem<- c("RL936_32F_P11_plus_10_32F_SmithR_to_E8_rep1_all_merged")
iPSCprimed <- PLEASEdata[ ,primedstem]
iPSCNP <- PLEASEdata[ ,NPstem]
wilcox.test(iPSCprimed, iPSCNP, paired = TRUE)

### Comparisons drawn between 32NSCs
keep_32NSCrow <- c("h9_WGBS", "RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge", "RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")

new_dat <- hyper_newdata_melt[hyper_newdata_melt$variable %in% keep_32NSCrow, ]
ggplot(new_dat, aes(x=variable, y = value)) +
  geom_boxplot(notch = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

?geom_boxplot
primedNSC <- c("RL1813_2019_09_19_32F_NSCs_primed_p7_S4_all_lanes_merge")
NPNSC<- c("RL1810_2019_09_19_32F_NSCs_N2P_p8_S1_all_lanes_merge")
NSCprimed <- PLEASEdata[ ,primedNSC]
NSCNP1 <- PLEASEdata[ ,NPNSC]
wilcox.test(NSCprimed, NSCNP1, paired = TRUE)















