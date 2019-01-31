##
##    mir21_DE.R
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

# For comparative and descriptive analysis we use the rlog (see manual for descriptions). blind=FALSE to keep sample group info for calculating normalizations.  
rldMM <- rlog(dds, blind = FALSE)


## PART 5: DIFFERENTIAL EXPRESSION

# CP4 vs. TP4
dds_resCP4TP4 <- results(dds, c("condition","CP4","TP4"))
dds_resCP4TP4 <- dds_resCP4TP4[order(dds_resCP4TP4$padj),]
DESeq2::plotMA(dds_resCP4TP4,  alpha = 0.05, main="MA plot CP4 vs TP4. Alpha level 0.05", ylim=c(-10,10))
write.table(dds_resCP4TP4,file="output/tables/dds_resCP4_TP4_ordered_170821.txt",col.names=T,quote=F,sep="\t")

# CP4 vs. TP5
dds_resCP4TP5 <- results(dds, c("condition","CP4","TP5"))
dds_resCP4TP5 <- dds_resCP4TP5[order(dds_resCP4TP5$padj),]
DESeq2::plotMA(dds_resCP4TP5,  alpha = 0.05, main="MA plot CP4 vs TP5. Alpha level 0.05", ylim=c(-10,10))
write.table(dds_resCP4TP5,file="output/tables/dds_resCP4_TP5_ordered_170821.txt",col.names=T,quote=F,sep="\t")

#TP4 vs. TP5
dds_resTP4TP5 <- results(dds, c("condition","TP4","TP5"))
dds_resTP4TP5 <- dds_resTP4TP5[order(dds_resTP4TP5$padj),]
DESeq2::plotMA(dds_resTP4TP5,  alpha = 0.05, main="MA plot TP4 vs TP5. Alpha level 0.05", ylim=c(-10,10))
write.table(dds_resTP4TP5,file="output/tables/dds_resTP4_TP5_ordered_170821.txt",col.names=T,quote=F,sep="\t")


# TP4 vs. CP4
dds_resTP4CP4 <- results(dds, c("condition","TP4","CP4"))
dds_resTP4CP4 <- dds_resTP4CP4[order(dds_resTP4CP4$padj),]
DESeq2::plotMA(dds_resTP4CP4,  alpha = 0.05, main="MA plot TP4 vs CP4. Alpha level 0.05", ylim=c(-10,10))
write.table(dds_resTP4CP4,file="output/tables/dds_resTP4_CP4_ordered_170821.txt",col.names=T,quote=F,sep="\t")

# TP5 vs. CP4
dds_resTP5CP4 <- results(dds, c("condition","TP5","CP4"))
dds_resTP5CP4 <- dds_resTP5CP4[order(dds_resTP5CP4$padj),]
DESeq2::plotMA(dds_resTP5CP4,  alpha = 0.05, main="MA plot TP5 vs CP4. Alpha level 0.05", ylim=c(-10,10))
write.table(dds_resTP5CP4,file="output/tables/dds_resTP5_CP4_ordered_170821.txt",col.names=T,quote=F,sep="\t")

#TP5 vs. TP4
dds_resTP5TP4 <- results(dds, c("condition","TP5","TP4"))
dds_resTP5TP4 <- dds_resTP5TP4[order(dds_resTP5TP4$padj),]
DESeq2::plotMA(dds_resTP5TP4,  alpha = 0.05, main="MA plot TP5 vs TP4. Alpha level 0.05", ylim=c(-10,10))
write.table(dds_resTP5TP4,file="output/tables/dds_resTP5_TP4_ordered_170821.txt",col.names=T,quote=F,sep="\t")
