##
##    Mir21_targets_heatmaps.R
##
##    This script
##

##    DESCRIPTION OF OPERATION FLOW:
##      Part 1: 
##

## PART 1: Load required packages
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

dds$condition <- factor(dds$condition,levels=c("CP4",
                                               "TP5",
                                               "TP4"))


dds <- DESeq(dds)

# For comparative and descriptive analysis we use the rlog (see manual for descriptions). blind=FALSE to keep sample group info for calculating normalizations.  
rldMM <- rlog(dds, blind = FALSE)

rld56_values <- assay(rldMM)

rm(cntmiR21all, colDataDMM)

# Reorder columns for graphing
colnames(rld56_values)
rld56_values <- rld56_values[,c("mock6", "mock7", "mock8", "pGHR10.1", "pGHR8.1", "pGHR9.1", "pGHR10.2", "pGHR8.2", "pGHR9.2")]


# load groups of genes as manually currated in the Groups_for_analysis.R file
source("R_scripts/shared/Groups_for_analysis.R")

# Now obtain the right DE order for EV-GFPpos vs WT microglia
dds_resT5C4 <- results(dds, c("condition","TP5","CP4"))
DESeq2::plotMA(dds_resT5C4,  alpha = 0.05, main="MA plot TP5 vs CP4", ylim=c(-10,10))
dds_resT5C4Ord_log2fold <- dds_resT5C4[order(dds_resT5C4$log2FoldChange),]


## PART 6: Setup plotting variables for the heatmaps
my_palette  = colorRampPalette(c("blue4","white", "red"))(n = 32)
col_breaks  = c(seq(2,14, length=33)) 
colsep      = c(0,3,6)
sepwidth    = c(0.05,0.05)
sepcolor    = "white"
cexCol      = 1.0
cexRow      = 0.8
lmat        = rbind(c(0,3),c(2,1),c(0,4))
lwid        = c(1.5,1) # second position for column width
lhei        = c(1.5,0.7,0.3) # second position for row height



# Also work on graph with all 9 samples
rld56_values_mir21_targets = subset(rld56_values, rownames(rld56_values) %in% Targets_miR21_strong_validated_181011)
rld56_values_mir21_targets = as.data.frame(rld56_values_mir21_targets)

dds_resT5C4Ord_log2fold_mir21_targets = subset(dds_resT5C4Ord_log2fold, rownames(dds_resT5C4Ord_log2fold) %in% Targets_miR21_strong_validated_181011)

write.table(dds_resT5C4Ord_log2fold_mir21_targets,file="output/tables/subsets/mir21_targets_strong/dds_resTP5_CP4_mir21_targets_validated_ordered_181011.txt",col.names=T,quote=F,sep="\t")

mir21_targets_order_T5C4 = rownames(dds_resT5C4Ord_log2fold_mir21_targets)
rld56_values_mir21_targets = rld56_values_mir21_targets[match(mir21_targets_order_T5C4, rownames(rld56_values_mir21_targets)),]


## PART 7: Perform the actual plotting and saving of the files
pdf(file="output/pdf/rlog_heatmaps/miR21_targets_strong/181011_mouse_miR21_targets_validated_sort_ALL_TP5CP4.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rld56_values_mir21_targets),
           main="miR21 Targets Strong T5 vs C4", # heat map title
           col=my_palette,
           breaks = c(seq(-2,2, length=33)),
           colsep=colsep,
           rowsep=c(1:nrow(rld56_values_mir21_targets)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", # turns off trace lines inside the heat map
           density.info="none", # turns off density plot inside color legend
           margins =c(25,12),     # widens margins around plot
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           #keysize=0.05,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=-7.4,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()






## Order T4 vs. C4

# Now obtain the right DE order for EV-GFPpos vs WT microglia
dds_resT4C4 <- results(dds, c("condition","TP4","CP4"))
DESeq2::plotMA(dds_resT4C4,  alpha = 0.05, main="MA plot TP4 vs CP4", ylim=c(-10,10))
dds_resT4C4Ord_log2fold <- dds_resT4C4[order(dds_resT4C4$log2FoldChange),]

# Also work on graph with all 9 samples
rld56_values_mir21_targets = subset(rld56_values, rownames(rld56_values) %in% Targets_miR21_strong_validated_181011)
rld56_values_mir21_targets = as.data.frame(rld56_values_mir21_targets)

dds_resT4C4Ord_log2fold_mir21_targets = subset(dds_resT4C4Ord_log2fold, rownames(dds_resT4C4Ord_log2fold) %in% Targets_miR21_strong_validated_181011)
write.table(dds_resT4C4Ord_log2fold_mir21_targets,file="output/tables/subsets/mir21_targets_strong/dds_resTP4_CP4_mir21_targets_validated_ordered_181011.txt",col.names=T,quote=F,sep="\t")

mir21_targets_order_T4C4 = rownames(dds_resT4C4Ord_log2fold_mir21_targets)
rld56_values_mir21_targets = rld56_values_mir21_targets[match(mir21_targets_order_T4C4, rownames(rld56_values_mir21_targets)),]


## PART 6: Setup plotting variables for the heatmaps
my_palette  = colorRampPalette(c("blue4","white", "red"))(n = 32)
col_breaks  = c(seq(2,14, length=33)) 
colsep      = c(0,3,6)
sepwidth    = c(0.05,0.05)
sepcolor    = "white"
cexCol      = 1.0
cexRow      = 0.8
lmat        = rbind(c(0,3),c(2,1),c(0,4))
lwid        = c(1.5,1) # second position for column width
lhei        = c(1.5,0.7,0.3) # second position for row height


## PART 7: Perform the actual plotting and saving of the files
pdf(file="output/pdf/rlog_heatmaps/miR21_targets_strong/181011_mouse_miR21_targets_validated_sort_ALL_TP4CP4.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rld56_values_mir21_targets),
           main="miR21 Targets Strong T4 vs C4", # heat map title
           col=my_palette,
           breaks = c(seq(-2,2, length=33)),
           colsep=colsep,
           rowsep=c(1:nrow(rld56_values_mir21_targets)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", # turns off trace lines inside the heat map
           density.info="none", # turns off density plot inside color legend
           margins =c(25,12),     # widens margins around plot
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           #keysize=0.05,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=-7.4,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()




## Order T5 vs. T4

# Now obtain the right DE order for EV-GFPpos vs WT microglia
dds_resT5T4 <- results(dds, c("condition","TP5","TP4"))
DESeq2::plotMA(dds_resT5T4,  alpha = 0.05, main="MA plot TP5 vs TP4", ylim=c(-10,10))
dds_resT5T4Ord_log2fold <- dds_resT5T4[order(dds_resT5T4$log2FoldChange),]

# Also work on graph with all 9 samples
rld56_values_mir21_targets = subset(rld56_values, rownames(rld56_values) %in% Targets_miR21_strong_validated_181011)
rld56_values_mir21_targets = as.data.frame(rld56_values_mir21_targets)

dds_resT5T4Ord_log2fold_mir21_targets = subset(dds_resT5T4Ord_log2fold, rownames(dds_resT5T4Ord_log2fold) %in% Targets_miR21_strong_validated_181011)
write.table(dds_resT5T4Ord_log2fold_mir21_targets,file="output/tables/subsets/mir21_targets_strong/dds_resTP5_TP4_mir21_targets_validated_ordered_181011.txt",col.names=T,quote=F,sep="\t")

mir21_targets_order_T5T4 = rownames(dds_resT5T4Ord_log2fold_mir21_targets)
rld56_values_mir21_targets = rld56_values_mir21_targets[match(mir21_targets_order_T5T4, rownames(rld56_values_mir21_targets)),]


## PART 6: Setup plotting variables for the heatmaps
my_palette  = colorRampPalette(c("blue4","white", "red"))(n = 32)
col_breaks  = c(seq(2,14, length=33)) 
colsep      = c(0,3,6)
sepwidth    = c(0.05,0.05)
sepcolor    = "white"
cexCol      = 1.0
cexRow      = 0.8
lmat        = rbind(c(0,3),c(2,1),c(0,4))
lwid        = c(1.5,1) # second position for column width
lhei        = c(1.5,0.7,0.3) # second position for row height


## PART 7: Perform the actual plotting and saving of the files
pdf(file="output/pdf/rlog_heatmaps/miR21_targets_strong/181011_mouse_miR21_targets_validated_sort_ALL_TP5TP4.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rld56_values_mir21_targets),
           main="miR21 Targets Strong T5 vs. T4", # heat map title
           col=my_palette,
           breaks = c(seq(-2,2, length=33)),
           colsep=colsep,
           rowsep=c(1:nrow(rld56_values_mir21_targets)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", # turns off trace lines inside the heat map
           density.info="none", # turns off density plot inside color legend
           margins =c(25,12),     # widens margins around plot
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           #keysize=0.05,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=-7.4,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()
