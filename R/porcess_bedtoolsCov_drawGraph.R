# Stage 3 
# Drawing the actual graph 

workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/20160329_transcriptGraph"

library(dplyr)
library(ggplot2)

# Now use the bed file to ; 
# 1) Generate exonCoverage percentage
# 2) Generae 5' to 3' graph

# Now use bedtools coverage to calculate individual base coverage
# module load apps/bedtools/2.25.0/gcc-4.4.7 

# Convert bam to bed
## bedtools coverage -a toy.bed -b sample.bed -d  -s 

#  bedtools coverage -a /scratch/zaki/file/ref/zaki_customBED/EnsDb.Mmusculus.v75_proteinCoding_example.bed -b WT_r1.bed -d  -s > eg.txt

# -a = the bed originated from  the BAM
# -b = The EnsDb reference 
# -d = report per base coverage
# -s = only report things with the same strand 
# -S = report opposite strandess (for duTP protocol)

# Two genes coverage
dat <- read.delim(file.path(workDir, "bedtoolsCov_output/eg.txt"), header=F)

# Divide data into list based on individual gene
dat_l <- split( dat , f = ee$V5 ) #column V5 is the ensbl gene_id

# For one gene
gene_strand <- dat_l[[1]]$V6 # column V6 has the strand infomation
gene_name <- dat_l[[1]]$V5 # column V9 has the gene name
cov_value <- dat_l[[1]]$V9 # column V9 contains the coverage
# Eg - want to bin the data into 100 bins
binfunc <- bins(1:nrow(dat_l[[1]]), target.bins = 100, max.breaks = 100)
# Extract result from binfunction
binlength <- as.numeric(binfunc$binct) # This contains the number of continous data that makes one bin
# Create a category of the bin
binlist <- list()
for (i in seq_along(binlength)){
  binlist[[i]] <-  rep(i, each=binlength[i])
}
bincateg <- unlist(binlist) 
# Data frame containing the coverage and the bin it belongs to
df_cov <- data.frame(value = cov_value,
                     bin = bincateg,
                     gene_name = gene_name,
                     strand = gene_strand)
# Get the mean of each bin
df_bin <- as.data.frame(group_by(df_cov, bin) %>% 
                          summarise(mean_cov = mean(value),
                                    gene_name = unique(gene_name),
                                    strand = unique(strand)))

plot(x=df_bin$bin, y=rev(df_bin$mean))

ggplot(df_bin, aes(x=bin, y=rev(mean_cov))) +
  geom_line() 

# ------------------------------------- #
# Processing in loops
# ------------------------------------- #

workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/20160329_transcriptGraph"

# Two genes coverage
dat <- read.delim(file.path(workDir, "bedtoolsCov_output/eg20.txt"), header=F)


# Divide data into list based on individual gene
dat_l <- split( dat , f = dat$V4 ) #column V4 is the ensmbl gene_id

# Splid data with dplyr
dat_plyr <- group_by(dat, V4)

cov_list <- list()
perc_cov <- vector()

# Make everything into a loop
for (i in 1:length(group_size(dat_plyr))){
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

# For one gene BEGIN
pos_gene <- filter(bin_cov, gene_name == "Gm12169") %>% group_by(bin) %>% 
  summarise(numOfgene = length(bin),
            totalValue = sum(mean_cov))

ggplot(pos_gene, aes(x=bin, y=totalValue)) +
  geom_line() 
# For one gene END

# Now for all gene
pos_gene <- filter(bin_cov, strand == "+") %>% group_by(bin) %>% 
  summarise(numOfgene = length(bin),
            totalValue = sum(mean_cov))

neg_gene <- filter(bin_cov, strand == "-") %>% group_by(bin) %>% 
  summarise(numOfgene = length(bin),
            totalValue = sum(mean_cov))

totalCov <- pos_gene$totalValue + rev(neg_gene$totalValue)

totalGene <- pos_gene$numOfgene + rev(neg_gene$numOfgene)

meanCov <- totalCov / totalGene

plot(x=1:100, y=meanCov)

df_bin <- data.frame(bin = 1:100,
                     mean_cov <- meanCov)

ggplot(df_bin, aes(x=bin, y=rev(mean_cov))) +
  geom_line() 


