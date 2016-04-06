# Generating a refernece bed file which contains the exonic regions of all protein coding genes
# 30/03/16


library(GenomicAlignments)
library(dplyr)

workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/20160329_transcriptGraph"

# --------------------------------------------# 
# Generate a reference of all non-overlapping exonic regions
# --------------------------------------------# 
# If using mm10 reference
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# OR 
#library(OrganismDbi)
#library(Mus.musculus)
# txdb <- makeOrganismDbFromTxDb(Mus.musculus))

# If using mm10 Ens genes
#library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
# mmens <- makeOrganismDbFromTxDb(txdb)
#txs <- exons(mmens, columns=c("TXNAME","GENEID"))

# If using Ensemble package directly
# Getting gene name https://support.bioconductor.org/p/74074/
# http://bioconductor.org/packages/release/bioc/html/ensembldb.html
library(EnsDb.Mmusculus.v75)
# If interested in exon region ; 
#txs <- exons(EnsDb.Mmusculus.v75, 
#             columns=c("tx_biotype","gene_id", "tx_id", "gene_name"),
#             filter=list(TxbiotypeFilter("protein_coding")))
#df <- as.data.frame(txs,
#                    row.names = paste0("g_", 1:length(txs)))

# ------------------------------------------------------------# 
# OLD METHOD BEGINS - DO NOT USE THIS METHOD
# ------------------------------------------------------------# 

# Since we are interested in non-overlaping exon region & only protein coding transcripts
DJE <- disjointExons(EnsDb.Mmusculus.v75, 
                     columns=c("tx_biotype","gene_id", "tx_id", "gene_name"),
                     filter=list(SeqnameFilter(c(1:20, "X", "Y")),
                                 TxbiotypeFilter("protein_coding")))
# SeqnameFilter - select chr 1 to 20 and sex chr
# TxbiotypeFilter - select only specfified transcript type

# Sanity check
# 1) Check to see if Hoxa9 (only known transcript was extracted)
df <- as.data.frame(DJE)
h <- df[grep("ENSMUSG00000038227", df$gene_id),]
# shows that only portein coding transcript was included

# 2) Shortest exon region?

# Now make the data frame more tidy
bed <- data.frame(seqnames = df$seqnames,
                  start = df$start,
                  end = df$end,
                  gene_id = as.character(df$gene_id),
                  exonic_part = df$exonic_part,
                  strand = df$strand) # Make sure strand is in the 6th column
bed_edit <- mutate(bed, width = (end - start) + 1)
sanity2 <- group_by(bed_edit, gene_id) %>% summarise(length = sum(width))
head(arrange(sanity2, -desc(length)))
#gene_id length
#(fctr)  (dbl)
#1 ENSMUSG00000036948      2
#2 ENSMUSG00000048910      2
#3 ENSMUSG00000090852      2
#4 ENSMUSG00000098650      3
#5 ENSMUSG00000050347     11
#6 ENSMUSG00000093575     11

# Does it make sense that the gene ENSMUSG00000036948 has a length of 2 base pair?
filter(bed, gene_id == "ENSMUSG00000036948")
# Checking ensembl website - http://feb2014.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000036948;r=5:138259656-138264046
# The gene should have more genes - therefore the script above gives in accurate results

# --------------------------------------------# 
# Try a different way to extract the data
# Select the longest transcript
# --------------------------------------------# 
# 1st stage) Get the longest transcript for each protein coding transcripts
listColumns(EnsDb.Mmusculus.v75)
txs <- transcripts(EnsDb.Mmusculus.v75, 
                   columns=c("tx_biotype","gene_id", "tx_id", "gene_name"),
                   filter=list(SeqnameFilter(c(1:20, "X", "Y")),
                               TxbiotypeFilter("protein_coding")))
# Convert to data frame
trans_df <- as.data.frame(txs)

# Sort 1st by gene_id then by transcript width (note the longest transcript should be 1st)
trans_df <- arrange(trans_df, gene_id, -width)

# Get the longest transcipt
longest_trans_df <- distinct(trans_df, gene_id)
longest_trans <- longest_trans_df$tx_id

# 2nd stage) Extract the exon region only for the longest transcipts
el <- exons(EnsDb.Mmusculus.v75, 
            columns=c("tx_biotype","gene_id", "tx_id", "gene_name", "gene_biotype"),
            filter=list(SeqnameFilter(c(1:20, "X", "Y")),
                        TxbiotypeFilter("protein_coding"),
                        TxidFilter(longest_trans)))
el_df <- as.data.frame(el)
# Now make the data frame more tidy
el_bed <- data.frame(seqnames = el_df$seqnames,
                     start = el_df$start,
                     end = el_df$end,
                     gene_id = as.character(el_df$gene_id),
                     gene_name = el_df$gene_name,
                     strand = el_df$strand) # Make sure strand is in the 6th column

# Sanity check again
el_bed_edit <- mutate(el_bed, width = (end - start) + 1)

sanity3 <- group_by(el_bed_edit, gene_id) %>% summarise(length = sum(width))
head(arrange(sanity3, -desc(length)))
#gene_id length
#(fctr)  (dbl)
#1 ENSMUSG00000095386     24
#2 ENSMUSG00000091041     27
#3 ENSMUSG00000090441     36
#4 ENSMUSG00000091768     39
#5 ENSMUSG00000079333     45
#6 ENSMUSG00000091109     45

# Checking ensembl website - http://feb2014.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000095386;r=4:104857329-104859137;t=ENSMUST00000179793
# ENSMUSG00000095386 has 24 base 

# The methods looks okay but one slight problem ;
# We are selecting the longest transcript, but the longest transcript does not necessarly means the transcipt which has the longest exon sequence
# Example - Kat6a ; http://feb2014.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000031540;r=8:22859535-22943259
# Longest transcript -  ENSMUST00000110696 - transcirpt length: 83725 || exonic length: 9041
# 2nd longest transcript -  ENSMUST00000044331 - transcirpt length: 83711 || exonic length: 9113

# ------------------------------------------------------------# 
# OLD METHOD ENDS 
# ------------------------------------------------------------# 

# --------------------------------------------# 
# Try a different way to extract the data
# Select the transcript with the longest exonic length (longest exonic base pair)
# --------------------------------------------#
# 1st) Extract the exonic region for all transcript
el2 <- exons(EnsDb.Mmusculus.v75, 
             columns=c("tx_biotype","gene_id", "tx_id", "gene_name", "gene_biotype"),
             filter=list(SeqnameFilter(c(1:20, "X", "Y")),
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
el3 <- exons(EnsDb.Mmusculus.v75, 
             columns=c("tx_biotype","gene_id", "tx_id", "gene_name", "gene_biotype", "exon_idx"),
             filter=list(SeqnameFilter(c(1:20, "X", "Y")),
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
#head(arrange(sanity4, -desc(length)))
# All seems okay

# 2) What is the shortest exon length
head(arrange(el3_bed_edit, -desc(width)), n =20)

# Start and end seem to be oppoite in ensembl website ;
# Eg - 
# filter(el3_bed_edit, gene_id == "ENSMUSG00000042564")

#seqnames    start      end            gene_id gene_name strand width
#1        15 79658681 79658956 ENSMUSG00000042564   Fam227a      -   276

# On website http://feb2014.archive.ensembl.org/Mus_musculus/Transcript/Exons?db=core;g=ENSMUSG00000042564;r=15:79609576-79658956;t=ENSMUST00000109648
#No.	Exon / Intron	    Start	     End	   Length	
#1	ENSMUSE00000895010	79658956	79658681	276

# This will be important when plotting graph 

# Save the bed file to disk on cluster
# /data/stemcell/sczaki/file/ref/allExon/zaki_created
write.table(el3_bed, file.path(workDir, "EnsDb.Mmusculus.v75_proteinCoding_transcript.bed"),
            quote=F, sep= "\t", row.names=F)

# Select only genes more than 100bp in length
df_100bp <- filter(sanity4, length > 100)
genes_100bp <- as.character(df_100bp$gene_id)
# dim(sanity4) - 22527 genes
# dim(df_100bp)
# 22472 / 22527 - genes are more than 100 bp

el_bed_100bp <- filter(el3_bed, gene_id %in% genes_100bp)
write.table(el_bed_100bp, file.path(workDir, "EnsDb.Mmusculus.v75_proteinCoding_transcript_100bp.bed"),
            quote=F, sep= "\t", row.names=F)

# Example bed to run initial coverage matrix 
# Hoxa9 [ENSMUSG00000038227] (- strand) & Lpar5 [ENSMUSG00000067714] (+ strand)

eg_bed <- filter(el3_bed, gene_name %in% c("Hoxa9", "Lpar5"))
write.table(eg_bed, file.path(workDir, "EnsDb.Mmusculus.v75_proteinCoding_example.bed"),
            quote=F, sep= "\t", row.names=F)

# Example bed of 20 genes ; 10 pos 10 neg

neg_random <- filter(el_bed_100bp, strand == "-") %>% sample_n(10, replace = TRUE)
pos_random <- filter(el_bed_100bp, strand == "+") %>% sample_n(10, replace = TRUE)

random_20 <- c(as.character(neg_random$gene_id), as.character(pos_random$gene_id))


bed_20 <- filter(el3_bed, gene_id %in% random_20 )
write.table(bed_20, file.path(workDir, "EnsDb.Mmusculus.v75_proteinCoding_example_20.bed"),
            quote=F, sep= "\t", row.names=F)

# Transfer the bed file to cluster
