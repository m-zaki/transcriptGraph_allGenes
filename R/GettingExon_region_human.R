# Generating a reference bed file which contains the exonic regions of all protein coding genes
# 30/03/16


library(GenomicAlignments)
library(dplyr)

library(EnsDb.Hsapiens.v75)

workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/20160329_transcriptGraph"
workDir_output <- file.path(workDir, "output")

# For all previous QC and methods refer to the GettingExon_region_mouse.R script

# --------------------------------------------# 
# Try a different way to extract the data
# Select the transcript with the longest exonic length (longest exonic base pair)
# --------------------------------------------#
# 1st) Extract the exonic region for all transcript
el2 <- exons(EnsDb.Hsapiens.v75, 
             columns=c("tx_biotype","gene_id", "tx_id", "gene_name", "gene_biotype"),
             filter=list(SeqnameFilter(c(1:22, "X", "Y")),
                         TxbiotypeFilter("protein_coding")))

# Convert to df
el2_df <- as.data.frame(el2, row.names=1:length(el2))

# Group all the exonic region into transcript
by_transcript <- group_by(el2_df, tx_id, gene_id) %>% 
  summarise(width = sum(width))

# Select the transcript with the longest exonic lenght
trans_df2 <- as.data.frame(by_transcript)
trans_df2 <- arrange(trans_df2, gene_id, -width)

# Select the longest unique one (should be at the top)
longest_exon_df <- distinct(trans_df2, gene_id)
longest_exon <- longest_exon_df$tx_id

# 2nd) Select the transcripts which has the longest exonic length
el3 <- exons(EnsDb.Hsapiens.v75, 
             columns=c("tx_biotype","gene_id", "tx_id", "gene_name", "gene_biotype", "exon_idx"),
             filter=list(SeqnameFilter(c(1:22, "X", "Y")),
                         TxbiotypeFilter("protein_coding"),
                         TxidFilter(longest_exon),
                         GenebiotypeFilter("protein_coding")))

# Sort the Granges object
el3 <- sortSeqlevels(el3)
el3 <- sort(el3)

el3_df <- as.data.frame(el3)

# Now make the data frame more tidy
el3_bed <- data.frame(seqnames = el3_df$seqnames,
                      start = el3_df$start,
                      end = el3_df$end,
                      gene_id = as.character(el3_df$gene_id),
                      gene_name = el3_df$gene_name,
                      strand = el3_df$strand,  # Make sure strand is in the 6th column
                      exon = el3_df$exon_idx)
# Sanity check again
# 1) Check which gene has the shortest length
el3_bed_edit <- mutate(el3_bed, width = (end - start) + 1)
sanity4 <- group_by(el3_bed_edit, gene_id) %>% summarise(length = sum(width))
head(arrange(sanity4, -desc(length)))
#gene_id length
#(fctr)  (dbl)
#1 ENSG00000268141     21
#2 ENSG00000268556     30
#3 ENSG00000269531     30
#4 ENSG00000269766     30
#5 ENSG00000268208     33
#6 ENSG00000269182     33

# 2) What is the shortest exon length
head(arrange(el3_bed_edit, -desc(width)), n =20)

# Save the bed file to disk on cluster
# /data/stemcell/sczaki/file/ref/allExon/zaki_created
write.table(el3_bed, file.path(workDir_output, "EnsDb.Hsapiens.v75_proteinCoding_transcript.bed"),
            quote=F, sep= "\t", row.names=F)

# Select only genes more than 100bp in length
df_100bp <- filter(sanity4, length > 100)
genes_100bp <- as.character(df_100bp$gene_id)
# dim(sanity4) - 20138 genes
# dim(df_100bp)
# 20056 / 20138 - genes are more than 100 bp

el_bed_100bp <- filter(el3_bed, gene_id %in% genes_100bp)
write.table(el_bed_100bp, file.path(workDir_output, "EnsDb.Hsapiens.v75_proteinCoding_transcript_100bp.bed"),
            quote=F, sep= "\t", row.names=F)

# Example bed of 20 genes ; 10 pos 10 neg
neg_random <- filter(el_bed_100bp, strand == "-") %>% sample_n(10, replace = TRUE)
pos_random <- filter(el_bed_100bp, strand == "+") %>% sample_n(10, replace = TRUE)

random_20 <- c(as.character(neg_random$gene_id), as.character(pos_random$gene_id))


bed_20 <- filter(el3_bed, gene_id %in% random_20 )
write.table(bed_20, file.path(workDir_output, "EnsDb.Hsapiens.v75_proteinCoding_example_20.bed"),
            quote=F, sep= "\t", row.names=F)

# Transfer the bed file to cluster
