# 2nd stage of plotting the 5' to 3' graph involves using bedtools

# Load bedtools module
module load apps/bedtools/2.25.0/gcc-4.4.7 

# Convert bam to bed
bedtools bamtobed -i input.bam > input.bed

# Use bedtools coverage with the following parameters

# -a = The EnsDb bed file generated from stage 1
# -b = the bed file which orginated from a BAM file
# -d = report per base coverage
##### Choose one of below depending on the experiment type
# -s = only report sequence with the same strand 
# -S = report opposite strandess (for duTP protocol)


bedtools coverage -a EnsDb.Mmusculus.v75_proteinCoding_example_20.bed -b input.bed -d -s > eg_20.txt
# If bed file contains header, delete the header 
