# Stage 3 
# Drawing the actual graph 

library(binr)
library(dplyr)
library(ggplot2)

# ------------------------------------- #
# Processing in loops
# ------------------------------------- #

workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/20160329_transcriptGraph"

# Read the an example data which contains example ouput from bedtools coverage (across 20 genes)
dat <- read.delim(file.path(workDir, "data/bedtoolsCov_output/eg20.txt"), header=F)

# Divide data into list based on individual gene 
# This step takes the longest on a big data (anyway of doing this with dplyr??)
dat_l <- split( dat , f = dat$V4 ) #column V4 is the ensembl gene_id

# Create empty list for loop
cov_list <- list()
perc_cov <- vector()

# Make everything into a loop 
# Refer previous github commit for the script version which does not deal with loops
for (i in seq_along(dat_l)){
  gene_strand <- as.character(dat_l[[i]]$V6) # column V6 has the strand infomation
  gene_name <- as.character(dat_l[[i]]$V5) # column V5 has the name
  gene_id <- as.character(dat_l[[i]]$V4) # column V4 has the gene_id
  cov_value <- dat_l[[i]]$V9 # column V9 contains the coverage
  
  # Getting the percentage of gene covered more than 1 base
  perc_cov[i] <- round(sum(cov_value > 1) / length(cov_value) * 100, 2)
  names(perc_cov)[i] <- unique(gene_id)
  
  # Eg - want to bin the data into 100 bins
  binfunc <- bins(1:nrow(dat_l[[i]]), target.bins = 100, max.breaks = 100)
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
  cov_list[[i]] <- as.data.frame(group_by(df_cov, bin) %>% 
                                   summarise(mean_cov = mean(value),
                                             gene_id = unique(gene_id),
                                             gene_name = unique(gene_name),
                                             strand = unique(gene_strand)))
}

# Unlist the cov_list into a data frame
bin_cov <- do.call("rbind", cov_list)

# ----------- For one gene BEGIN
# If we are interested to plot just for one gene
pos_gene <- filter(bin_cov, gene_name == "Eif4g3") %>% group_by(bin) %>% 
  summarise(numOfgene = length(bin),
            totalValue = sum(mean_cov))

ggplot(pos_gene, aes(x=bin, y=totalValue)) +
  geom_line() 
# ----------- For one gene END

# Now plotting for all gene
# Need to seperate genes according to strand because ;
# for positive gene the transcirpt orientation is : 5' ----------- 3' 
# for negative gene the transcirpt orientation is : 3' ----------- 5' 
# Therefore need to reverese the direction of negative gene
pos_gene <- filter(bin_cov, strand == "+") %>% group_by(bin) %>% 
  summarise(numOfgene = length(bin),
            totalValue = sum(mean_cov))

neg_gene <- filter(bin_cov, strand == "-") %>% group_by(bin) %>% 
  summarise(numOfgene = length(bin),
            totalValue = sum(mean_cov))

# Combine the coverage from positive gene and negative gene according to thier bins
totalCov <- pos_gene$totalValue + rev(neg_gene$totalValue)
totalGene <- pos_gene$numOfgene + rev(neg_gene$numOfgene) # Total number of genes (needs the number to calculate mean)

meanCov <- totalCov / totalGene

df_bin <- data.frame(bin = 1:100,
                     mean_cov <- meanCov)

ggplot(df_bin, aes(x=bin, y=rev(mean_cov))) +
  geom_line() +
  ylab("Read count") +
  xlab("5' to 3' coverage")


# ----------------------------- #
# Check strand bias
# ----------------------------- #

# Combine the coverage from positive gene and negative gene according to thier bins
totalCov_strand <- c(pos_gene$totalValue , rev(neg_gene$totalValue))
totalGene_strand <- c(pos_gene$numOfgene , rev(neg_gene$numOfgene)) # Total number of genes (needs the number to calculate mean)

meanCov_strand <- totalCov_strand / totalGene_strand

df_bin_strand <- data.frame(bin = rep(1:100,2),
                            mean_cov = meanCov_strand,
                            strand = c(rep("pos", 100),
                                       rep("neg", 100)))


ggplot(df_bin_strand, aes(x=bin, y=rev(mean_cov), colour = strand)) +
  geom_line() +
  ylab("Read count") +
  xlab("5' to 3' coverage")
