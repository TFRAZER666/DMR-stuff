library(GenomicRanges)
library(bsseq)
library(magrittr)
library(pheatmap)

# Import the methylation data ---------------------------------------------

# Create a list of the BSseq objects that have the CG methylation data
bs_obj_list <- list.files(path = "~/Desktop/No", pattern = "_CG_bsseq_obj.Rds", full.names = TRUE)


# Read a Bs_seq boject from .Rds file
read_bs_obj <- function(rds_path){
    
    bs_obj <- readRDS(file = rds_path)
    
    return(bs_obj)
}


# Read these objects into the R session in the form of a list
all_bs_obj <- lapply(bs_obj_list, read_bs_obj)

# combine the list of objects into one BSseq object for downstream analyses
all_bs_obj <- combineList(all_bs_obj)

# Import the DMR coordinates ----------------------------------------------

##### A function that imports a bed file to a GRanges object
# Only takes the first 3 columns from a bed file
bed_to_gr <- function(bed_path){
    dat <- fread(bed_path, sep = "\t", header = FALSE)
    gr <- GRanges(seqnames = dat$V1,
                  ranges = IRanges(start = dat$V2,
                                   end = dat$V3))
    return(gr)
    
}

memory_gr <-  bed_to_gr("mcg_dmr_memory.bed")
ips_gr <-  bed_to_gr("mcg_dmr_ips.bed")

# Calculate methylation levels for DMRs -----------------------------------

# Get the methylation calls for each sample as a matrix
memory_m <- getCoverage(all_bs_obj, regions = memory_gr, type = "M", what = "perRegionTotal")

# Get the coverage of each sample 
memory_cov <- getCoverage(all_bs_obj, regions = memory_gr, type = "Cov", what = "perRegionTotal")

# Calculate methylation levels for DMRs as mCG/CG
memory_mCG <- memory_m / memory_cov


#####  Also try and make the mCG matrix by taking the per-base average instead of
#####  the per-region weighted total to see if it influences the resuts. 

# Plot DMRs ---------------------------------------------------------------

# Heatmap. Note that the default settings are not best here. Play with them to explore. 
pheatmap(memory_mCG)

### Now make the heatmap for the ips DMRs

