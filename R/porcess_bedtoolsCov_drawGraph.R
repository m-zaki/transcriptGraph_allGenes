library(binr)
library(dplyr)
library(ggplot2)

source("R/functions.R")

workDir <- "data/bedtoolsCov_output/example_data"

# List all the files 
files <- list.files(workDir, pattern="*.txt",full.names = TRUE) # List all .txt files
names(files) <- sapply(files, function(fn) { gsub('.txt', '', basename(fn)) })

perc_sample_list <- list()
mean_bin_list <- list()

for (r in seq_along(files)){
  # Read the an example data which contains example ouput from bedtools coverage (across 20 genes)
  dat <- read.delim(files[r], header=F)
  
  # Divide data into list based on individual gene 
  # This step takes the longest on a big data (anyway of doing this with dplyr??)
  dat_l <- split( dat , f = dat$V4 ) #column V4 is the ensembl gene_id
  
  # Create some temporary objects
  perc_cov_list <- list()
  gene_name <- vector() # Store the gene_name
  val_cov  <- vector() # Temporatly store the percentage of genes covered of individual gene
  
  # Use the function "per.covered"
  for (i in seq_along(dat_l)) {
    perc_cov_list[[i]] <- perc.covered(input = dat_l, min.depth = 1)
    gene_name[i] <- names(perc_cov_list[[i]])
    val_cov[i] <- as.vector(perc_cov_list[[i]])
  }
  
  perc_cov_df <- data.frame(gene_name = gene_name,
                            coverage_perc = val_cov)
  colnames(perc_cov_df)[2] <- paste0("perc_", names(files)[r]) # Rename the columns so each sample has a different name
  perc_sample_list[[r]] <- perc_cov_df
  
  
  # Use the function "gene.bin.mean"
  mean_bin_list[[i]] <- gene.bin.mean(input = dat_l)
  
  # Unlist the object
  bin_cov <- do.call("rbind", mean_bin_list)
  colnames(bin_cov)[2] <- paste0("meanCov_", names(files)[r]) # Rename the columns so each sample has a different name
  
  saveRDS(bin_cov, file = paste0(names(files[r]) ,"_cov.rds"))
}
  







