##
##    mir21_norm_reads.R
##
##    Description here
##

##    DESCRIPTION OF OPERATION FLOW:
##      Part 1: 
##

##      Part 1: Install required packages

## PART 1: Load required packages

# Load required packages
library(DESeq2)
library(gplots)
library(limma)

## PART 2: source the file generating the dataframe. 
source("R_scripts/shared/Datasets_microglia.R")


## PART 3: Setup the DESeq2 model and calculate the rlog values

### DESeq2 Setup
colDataDMM <- data.frame(condition=rep(NA,length(colnames(cntmiR21all))), row.names=colnames(cntmiR21all))
colDataDMM$condition[grep("m....",rownames(colDataDMM))] <- "CP4"
colDataDMM$condition[grep("pGHR8.1",rownames(colDataDMM))] <- "TP4"
colDataDMM$condition[grep("pGHR9.1",rownames(colDataDMM))] <- "TP4"
colDataDMM$condition[grep("pGHR10.1",rownames(colDataDMM))] <- "TP4"
colDataDMM$condition[grep("pGHR8.2",rownames(colDataDMM))] <- "TP5"
colDataDMM$condition[grep("pGHR9.2",rownames(colDataDMM))] <- "TP5"
colDataDMM$condition[grep("pGHR10.2",rownames(colDataDMM))] <- "TP5"

dds <- DESeqDataSetFromMatrix(countData = cntmiR21all,colData = colDataDMM,design = ~condition)

dds <- DESeq(dds)

rm(colDataDMM, cntmiR21all)

## For the supplementary tables we work with the normalized gene counts
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- normalized_counts[,c("mock6", "mock7", "mock8",
                                          "pGHR8.1", "pGHR9.1", "pGHR10.1",
                                          "pGHR8.2", "pGHR9.2", "pGHR10.2")]


CP4_normalized_counts <- as.data.frame(normalized_counts[, c("mock6", "mock7", "mock8")])
TP4_normalized_counts <- as.data.frame(normalized_counts[, c("pGHR8.1", "pGHR9.1", "pGHR10.1")])
TP5_normalized_counts <- as.data.frame(normalized_counts[, c("pGHR8.2", "pGHR9.2", "pGHR10.2")])

CP4_normalized_counts$mean_ncounts = rowMeans(CP4_normalized_counts[,colnames(CP4_normalized_counts)], na.rm=FALSE)
TP4_normalized_counts$mean_ncounts = rowMeans(TP4_normalized_counts[,colnames(TP4_normalized_counts)], na.rm=FALSE)
TP5_normalized_counts$mean_ncounts = rowMeans(TP5_normalized_counts[,colnames(TP5_normalized_counts)], na.rm=FALSE)

CP4_normalized_counts = CP4_normalized_counts[c("mean_ncounts")]
TP4_normalized_counts = TP4_normalized_counts[c("mean_ncounts")]
TP5_normalized_counts = TP5_normalized_counts[c("mean_ncounts")]

Avg_norm_reads = data.frame(CP4_normalized_counts, TP4_normalized_counts, TP5_normalized_counts)
colnames(Avg_norm_reads) <- c("CP4_Average", "TP4_Average", "TP5_Average")

# Extract subset and write files to disk
source("R_scripts/shared/Groups_for_analysis.R")

write.table(normalized_counts,file="output/tables/norm_reads/All_samples_individual_180202_ncounts.txt",col.names=T,quote=F,sep="\t")
write.table(Avg_norm_reads,file="output/tables/norm_reads/Group_averages_180202_ncounts.txt",col.names=T,quote=F,sep="\t")

# Subset for miRTarBase_miR21_strong set
individual_samples_mir21_targets = subset(normalized_counts, rownames(normalized_counts) %in% miRTarBase_miR21_strong_order_CP4_TP5)
individual_samples_mir21_targets = individual_samples_mir21_targets[match(miRTarBase_miR21_strong_order_CP4_TP5, rownames(individual_samples_mir21_targets)),]
write.table(individual_samples_mir21_targets, file="output/tables/norm_reads/subsets/miRTarBase_miR21_strong/All_samples_individual_miRTarBase_miR21_strong_180202_ncounts.txt",col.names=T,quote=F,sep="\t")

Avg_norm_reads_mir21_targets = subset(Avg_norm_reads, rownames(Avg_norm_reads) %in% miRTarBase_miR21_strong_order_CP4_TP5)
Avg_norm_reads_mir21_targets = Avg_norm_reads_mir21_targets[match(miRTarBase_miR21_strong_order_CP4_TP5, rownames(Avg_norm_reads_mir21_targets)),]
write.table(Avg_norm_reads_mir21_targets, file="output/tables/norm_reads/subsets/miRTarBase_miR21_strong/Group_averages_miRTarBase_miR21_strong_180202_ncounts.txt",col.names=T,quote=F,sep="\t")


