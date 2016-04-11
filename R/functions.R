# Function to extract the percentage of gene covered
perc.covered <- function(input, min.depth) {
  dat_l <- input[[i]]
  gene_id <- as.character(dat_l$V4) # column V4 has the gene_id
  cov_value <- dat_l$V9 # column V9 contains the coverage
  
  # Getting the percentage of gene covered more than 1 base
    perc_cov <- round(sum(cov_value >= min.depth) / length(cov_value) * 100, 2)
  names(perc_cov) <- unique(gene_id)
  
  return(perc_cov)
}

# Function to extract the mean value of each bin
gene.bin.mean <- function(input) { 
  dat_l <- input[[i]]
  gene_strand <- as.character(dat_l$V6) # column V6 has the strand infomation
  gene_name <- as.character(dat_l$V5) # column V5 has the name
  gene_id <- as.character(dat_l$V4) # column V4 has the gene_id
  cov_value <- dat_l$V9 # column V9 contains the coverage
  
  # Eg - want to bin the data into 100 bins
  binfunc <- bins(1:nrow(dat_l), target.bins = 100, max.breaks = 100)
  # Extract result from binfunction
  binlength <- as.numeric(binfunc$binct) # This contains the number of continous data that makes one bin
  # Create a category of the bin
  binlist <- list()
  for (k in seq_along(binlength)){
    binlist[[k]] <-  rep(k, each=binlength[k])
  }
  bincateg <- unlist(binlist) 
  # Data frame containing the coverage and the bin it belongs to
  df_cov <- data.frame(value = cov_value,
                       bin = bincateg,
                       gene_id = gene_id, 
                       gene_name = gene_name,
                       strand = gene_strand)
  # Get the mean of each bin
  cov_mean <- as.data.frame(group_by(df_cov, bin) %>% 
                              summarise(mean_cov = mean(value),
                                        gene_id = unique(gene_id),
                                        gene_name = unique(gene_name),
                                        strand = unique(gene_strand)))
  return(cov_mean)
}