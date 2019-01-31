## PART 2: Load .txt file with mapped read count and remove genes + samples with no valuable data
cntmiR21 <- read.table(file="input/star_genes_erc.counts.txt")


# Remove features (certain rownames) we no longer need
colnames(cntmiR21) <- strsplit2(colnames(cntmiR21),"_")[,1]
htseq_drop_rows <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique","ERCC-")
drop <- rowSums(sapply(htseq_drop_rows, grepl, rownames(cntmiR21)))>0
cntmiR21all <- cntmiR21[!drop,]

# Remove stuff we no longer need
rm(cntmiR21, drop, htseq_drop_rows)


# Print dimensions of 
dim(cntmiR21all)

# Keep genes with more than 1 read in at least 2 samples
keep <- rowSums(cntmiR21all>=1)>=2
table(keep)
cntmiR21all <- cntmiR21all[keep,]
dim(cntmiR21all)

# Drop low read count (less than 5 reads) samples. First show these
hist(colSums(cntmiR21all>5), breaks=20)

# Print colSums output to spot the samples with less counts
colSums(cntmiR21all>5)

# Drop samples where less than 5000 genes were detected with more than 5 reads
drop <- colSums(cntmiR21all>5)<5000
table(drop)
cntmiR21all <- cntmiR21all[,!drop]

rm(drop, keep)