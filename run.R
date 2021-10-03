library(ggplot2)
library(gridExtra)
library(DESeq2)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(RColorBrewer)
library(factoextra)
library(Seurat)
library(readxl)
library(ggpubr)
library(factoextra)
library(readxl)
library( DESeq2 )
library( genefilter )
library( statmod )
library(ggplot2)
library(Rtsne) 
                  
setwd("/home/hiri/Desktop/PhD_Projects/2020-09-21-Christian")
dir.create("output")
                          
## Importing the count table
temp <- list.files("./temp", pattern = "out.txt")
col <- gsub("_S.*", "", gsub("Aligned.out.txt", "", temp))
time <- gsub("_.*","",col)
nCell <- gsub("_.*","",gsub(".*d_","",gsub(".*h_","",col)))
                      
filenames <- list.files(path = "./temp", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)[,6]})
unique_table <- sapply(datalist, cbind)
colnames(unique_table) <- col
rownames(unique_table) <- rownames(read.table("./temp/0h_1_A3_1Aligned.out.txt", row.names = 1, header = TRUE))
              
write.csv(unique_table, file = "raw_counts.csv")
              
colData <- data.frame(col, time = time, nCell = nCell)
              
              
#############################################
#############################################
############################################
## Adding the number of detected genes and library size
colData$n.gene <- colSums(unique_table > 5)
colData$lib.size <- colSums(unique_table)
write.csv(colData, file="colData.csv")
              
              
colData50 <- colData[colData$nCell == 50,]
colData1 <- colData[colData$nCell == 1,]
unique_table1 <- unique_table[,colData1$col]
unique_table50 <- unique_table[,colData50$col]          
              
              
colData1$time <- factor(colData1$time,levels = c("0h", "4h", "12h", "24h", "7d"))
              
p1 <- ggplot(data = colData1, aes(x = time, y = n.gene, color = time)) +
   geom_jitter() + 
   geom_violin(alpha = 0) +
   theme_classic() +
   ylab("Number of detected genes") +
   xlab("") +
   theme(legend.position = "none")
              
p2 <- ggplot(data = colData1, aes(x = time, y = lib.size, color = time)) +
   geom_jitter() + 
   geom_violin(alpha = 0) +
   theme_classic() +
   ylab("library size") +
   xlab("") +
   theme(legend.position = "none")  +
   scale_y_continuous(trans='log10')
              
p <- grid.arrange(p1, p2, ncol=1)
ggsave("./output/gene_lib_1.eps", Fig1c, width = 6, height = 5)
                          
ggplot(data = colData1, aes(x = lib.size, y=n.gene, color = time)) +
   geom_jitter() +
   theme_classic() +
   theme() +
   ylab("Number of detected genes") +
   xlab("Library Size") +
   geom_hline(yintercept = 200) +
   geom_hline(yintercept = 1000) +
   geom_vline(xintercept = 200000) +
   geom_vline(xintercept = 2000000)
              
ggsave("./output/gene_lib_cor_1.eps", width = 5, height = 5)
                          
p1 <- ggplot(data = colData50, aes(x = time, y = n.gene, color = time)) +
   geom_jitter() + 
   geom_violin(alpha = 0) +
   theme_classic() +
   ylab("Number of detected genes") +
   xlab("") +
   theme(legend.position = "none")
              
p2 <- ggplot(data = colData50, aes(x = time, y = lib.size, color = time)) +
   geom_jitter() + 
   geom_violin(alpha = 0) +
   theme_classic() +
   ylab("library size") +
   xlab("") +
   theme(legend.position = "none")  +
   scale_y_continuous(trans='log10')
              
Fig1c <- grid.arrange(p1, p2, ncol=1)
ggsave("./output/gene_lib_50.png", Fig1c, width = 6, height = 5)
              
ggplot(data = colData50, aes(x = lib.size, y=n.gene, color = time)) +
    geom_jitter() +
    theme_classic() +
    theme() +
    ylab("Number of detected genes") +
    xlab("Library Size") 
              
ggsave("./output/gene_lib_cor_50.png", width = 5, height = 5) 
              


            
########################################
################# Single ###############
########################################

colData_filt <- colData1[colData1$n.gene > 200 , ]
colData_filt <- colData_filt[colData_filt$n.gene < 1000 , ]
colData_filt <- colData_filt[colData_filt$lib.size > 200000 , ]
colData_filt <- colData_filt[colData_filt$lib.size < 2000000 , ]
colData_filt <- colData_filt[-c(7),]
            
colData_filt$col <- as.character(colData_filt$col)
colData_filt$time <- as.character(colData_filt$time)
colData_filt$nCell <- as.character(colData_filt$nCell)
rownames(colData_filt) <- colData_filt$col
outlier <- c("0h_1_B5_1", "0h_1_G2", "0h_1_C5_1")
colData_filt <- colData_filt[!rownames(colData_filt) %in% outlier, ]
            
colData_write <- colData_filt
rownames(colData_write) <- colData_write$col
outlier <- c("0h_1_B5_1", "0h_1_G2", "0h_1_C5_1")
colData_write <- colData_write[!rownames(colData_write) %in% outlier, ]
write.csv(colData_write, file = "colData_w.csv")
            
            
            
single_table <- unique_table1[c(1:57,157:11807), as.character(colData_filt$col)]
pos <- grep("mVAT",rownames(single_table))  
pos <- c(pos,11190)
single_table <- single_table[-pos, ]
single_table <- single_table[rowSums(single_table) > 20, ]
            
                        
### Normalizing data
sfSingle <- estimateSizeFactorsForMatrix(single_table)
nSingle <- t(t(single_table) / sfSingle)
            
##########################################################################
## Selecting the top variable genes and log transformation
topVarGenesSingle <- head(order(rowVars(nSingle), decreasing = TRUE),1000)
nSingle_fin <- log(nSingle[topVarGenesSingle,] + 1)
### Performing PCA (single cell)
res.pca <- prcomp(t(nSingle_fin), scale = FALSE)
#colData1$time <- factor(colData1$time,levels = c("0h", "4h", "12h", "24h", "7d"))
ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], time = factor(colData_filt$time,levels = c("0h", "4h", "12h", "24h", "7d"))), aes(x = PC1, y = PC2, color = time)) +
      geom_jitter(size = 2) +
      theme_classic() +
      theme() +
      geom_hline(yintercept=0, linetype="dashed") +
      geom_vline(xintercept=0, linetype="dashed")
              
ggsave("./output/pca.eps", width = 6.5, height = 5)
              
              
write.csv(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], time = colData_filt$time), file = "coordinate.csv")
              
              
              
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
ggsave("./output/scree_plot.png", width = 5, height = 5)

              
df <- data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], time = colData_filt$time)  
df <- df[!df$time == "7d",]
              
ggplot(df, aes(x = PC1, y = PC2, color = time)) +
    geom_jitter(size = 2) +
    theme_classic() +
    theme() +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed")
              
ggsave("./output/pca_7d_removed.png", width = 5, height = 5)
              
write.csv(df, file = "coordinate_d7_removed.csv")
              
              
### outliers === 0h_1_B5_1 0h_1_G2 0h_1_C5_1
#s <- rownames(df) 
outlier <- c("0h_1_B5_1", "0h_1_G2", "0h_1_C5_1")
df <- df[!rownames(df) %in% outlier, ]
              
              
ggplot(df, aes(x = PC1, y = PC2, color = time)) +
    geom_jitter(size = 2) +
    theme_classic() +
    theme() +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed")
              
ggsave("./output/pca_7d_outlier_removed.png", width = 5, height = 5)
write.csv(df, file = "coordinate_d7_outlier_removed.csv")
              
#############################################################
#################### DE 7d vs 24h ###########################
#############################################################
colData_pair <- colData_filt[colData_filt$time == "7d" | colData_filt$time == "24h",]
rownames(colData_pair) <- colData_pair$col
colData_pair <- colData_pair[!rownames(colData_pair) %in% outlier, ]
tab <- single_table[,rownames(colData_pair)]
dds <- DESeqDataSetFromMatrix(countData = tab,
       colData = colData_pair,
       design= ~time)
dds <- DESeq(dds)
resultsNames(dds)
res_7d_vs_24h <- results(dds, name = "time_7d_vs_24h")
res_7d_vs_24h <- data.frame(res_7d_vs_24h)
res_7d_vs_24h <- res_7d_vs_24h[complete.cases(res_7d_vs_24h),]
write.csv(res_7d_vs_24h, file = "res_7d_vs_24h.csv")
              
res_7d_vs_24h_hm <- res_7d_vs_24h[res_7d_vs_24h$padj < 0.01,]
positive <- res_7d_vs_24h_hm[res_7d_vs_24h_hm$log2FoldChange > 2,]
negative <- res_7d_vs_24h_hm[res_7d_vs_24h_hm$log2FoldChange < -2,]
              
hm_gene <- c(rownames(positive), rownames(negative))
norm <- nSingle[hm_gene, rownames(colData_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
group <- factor(colData_DE$time)
ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "24h" = "#F8766D", "4h" = "#6a3d9a")))
pdf("./output/heatmap_res_7d_vs_24h.pdf", width = 14, height = 10)
p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 2), show_column_names = TRUE)
print(p)
dev.off()
              
              
#############################################################
#################### DE 7d vs 0h ###########################
#############################################################
colData_pair <- colData_filt[colData_filt$time == "7d" | colData_filt$time == "0h",]
rownames(colData_pair) <- colData_pair$col
colData_pair <- colData_pair[!rownames(colData_pair) %in% outlier, ]
tab <- single_table[,rownames(colData_pair)]
dds <- DESeqDataSetFromMatrix(countData = tab,
       colData = colData_pair,
       design= ~time)
dds <- DESeq(dds)
resultsNames(dds)
res_7d_vs_0h <- results(dds, name = "time_7d_vs_0h")
res_7d_vs_0h <- data.frame(res_7d_vs_0h)
res_7d_vs_0h <- res_7d_vs_0h[complete.cases(res_7d_vs_0h),]
write.csv(res_7d_vs_0h, file = "res_7d_vs_0h.csv")
              
res_7d_vs_0h_hm <- res_7d_vs_0h[res_7d_vs_0h$padj < 0.01,]
positive <- res_7d_vs_0h_hm[res_7d_vs_0h_hm$log2FoldChange > 2,]
negative <- res_7d_vs_0h_hm[res_7d_vs_0h_hm$log2FoldChange < -2,]
              
hm_gene <- c(rownames(positive), rownames(negative))
norm <- nSingle[hm_gene, rownames(colData_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
group <- factor(colData_DE$time)
ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "0h" = "#F8766D", "4h" = "#6a3d9a")))
pdf("./output/heatmap_res_7d_vs_0h.pdf", width = 14, height = 10)
p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 2), show_column_names = TRUE)
print(p)
dev.off()
              
              
              
#############################################################
#############################################################
############################################################
             
## performing DE analysis
colData_DE <- colData_filt#[!colData_filt$time == "7d",]
rownames(colData_DE) <- colData_DE$col
colData_DE <- colData_DE[!rownames(colData_DE) %in% outlier, ]
              
h0 <- factor(gsub("7d","rest",gsub("4h","rest",gsub("12h","rest",gsub("24h", "rest", gsub("0h", "aim", colData_DE$time))))))
h4 <- factor(gsub("7d","rest",gsub("4h","aim",gsub("12h","rest",gsub("24h", "rest", gsub("0h", "rest", colData_DE$time))))))
h12 <- factor(gsub("7d","rest",gsub("4h","rest",gsub("12h","aim",gsub("24h", "rest", gsub("0h", "rest", colData_DE$time))))))
h24 <- factor(gsub("7d","rest",gsub("4h","rest",gsub("12h","rest",gsub("24h", "aim", gsub("0h", "rest", colData_DE$time))))))
d7 <- factor(gsub("7d","aim",gsub("4h","rest",gsub("12h","rest",gsub("24h", "rest", gsub("0h", "rest", colData_DE$time))))))
colData_DE <- data.frame(colData_DE, h0, h4, h12, h24, d7)
#  cts <- single_table[,as.character(colData_DE$col) ]
#cts <- cts[rowSums(cts) > 0, ]
#cts <- cts[rowSums(cts > 0) > 5, ]
              
#  col_data <- data.frame(colData_DE, condition = factor(condition))
#  rownames(col_data) <- col_data$col
#  dds <- DESeqDataSetFromMatrix(countData = cts,
#          colData = col_data,
#          design= ~condition)
#  resultsNames(dds) 
#h4 <- as.character(colData_DE[colData_DE$time == "4h", ]$col)
#h0 <- as.character(colData_DE[colData_DE$time == "0h", ]$col)
#h12 <- as.character(colData_DE[colData_DE$time == "12h", ]$col)
#h24 <- as.character(colData_DE[colData_DE$time == "24h", ]$col)
#ht_order <- c(h0, h4, h12, h24)  
single_table <- single_table[,rownames(colData_DE)]
#condition <- factor(c(rep("early", 71), rep("late", 49)))
condition <- factor(colData_DE$time)
            
              
library(scde) 
names(condition) <- as.character(colnames(single_table)) 
sg <- condition
cd <-   single_table
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 8, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
            
            
### 0h
#  group <- factor(colData_DE$h0)
#  names(group) <- rownames(o.ifm)
#  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
#  de_0h <- ediff[order(ediff$Z, decreasing  =  TRUE), ]
#  write.csv(de_0h, file = "de_0h.csv")
            
            
### 4h
#  group <- factor(colData_DE$h4)
#  names(group) <- rownames(o.ifm)
#  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
#  de_4h <- ediff[order(ediff$Z, decreasing  =  TRUE), ]
#  write.csv(de_4h, file = "de_4h.csv")
            
            
### 12h
#  group <- factor(colData_DE$h12)
#  names(group) <- rownames(o.ifm)
#  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
#  de_12h <- ediff[order(ediff$Z, decreasing  =  TRUE), ]
#  write.csv(de_12h, file = "de_12h.csv")
            
            
### 24h
#  group <- factor(colData_DE$h24)
#  names(group) <- rownames(o.ifm)
#  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
#  de_24h <- ediff[order(ediff$Z, decreasing  =  TRUE), ]
#  write.csv(de_24h, file = "de_24h.csv")
            
            
 ### 7d
 #group <- factor(colData_DE$d7)
 #names(group) <- rownames(o.ifm)
 #ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
 #de_d7 <- ediff[order(ediff$Z, decreasing  =  TRUE), ]
 #write.csv(de_d7, file = "de_d7.csv")
            
            
 #  de_ht_0h <- de_0h[de_0h$cZ > 1.65, ]
 #  de_ht_4h <- de_4h[de_4h$cZ > 1.65, ]
 #  de_ht_12h <- de_12h[de_12h$cZ > 1.65, ]
 #  de_ht_24h <- de_24h[de_24h$cZ > 1.65, ]
 #de_ht_d7 <- de_d7[de_d7$cZ > 1.96, ]
            
 #  de_sig <- rbind(de_ht_0h, de_ht_4h, de_ht_12h, de_ht_24h)
 #  de_sig <- rownames(de_sig)
 #  de_sig <- unique(de_sig)
            
            
            
            
 #de_a <- de[de$cZ > 1.96, ]
 #de_b <- de[de$cZ < -1.96, ]
 #de_sig <- rbind(de_a, de_b)
 #de_sig <- rownames(de_sig)
            
 ### creating the heatmap
 h4 <- as.character(colData_DE[colData_DE$time == "4h", ]$col)
 h0 <- as.character(colData_DE[colData_DE$time == "0h", ]$col)
 h12 <- as.character(colData_DE[colData_DE$time == "12h", ]$col)
 h24 <- as.character(colData_DE[colData_DE$time == "24h", ]$col)
 d7 <- as.character(colData_DE[colData_DE$time == "7d", ]$col)
 ht_order <- c(h0, h4, h12, h24, d7)
            
            
 #  norm <- nSingle[de_sig, ht_order]
 #  scaled <- t(apply(norm, 1, scale))
 #  colnames(scaled) <- colnames(norm)
 #group <- factor(colData_DE$time)
 #ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "24h" = "#F8766D", "4h" = "#6a3d9a")))
 #  pdf("./output/heatmap_comparision.pdf", width = 14, height = 7)
 #  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 4), show_column_names = TRUE)
 #  print(p)
 #  dev.off()
              
            
            
            
            
 #norm <- nSingle[hm_gene, ht_order]
 #  scaled <- t(apply(norm, 1, scale))
 #  colnames(scaled) <- colnames(norm)
 #group <- factor(colData_DE$time)
 #ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "24h" = "#F8766D", "4h" = "#6a3d9a")))
 #  pdf("./output/heatmap_comparision.pdf", width = 14, height = 7)
 #  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 4), show_column_names = TRUE)
 #  print(p)
 #  dev.off()
            
 ##################################################################
 ######################## Performing DEseq2 #######################
 ##################################################################
 condition <- colData_DE
          
 ##h0
 dds <- DESeqDataSetFromMatrix(countData = single_table,
      colData = condition,
      design= ~h0)
 dds <- DESeq(dds)
 resultsNames(dds)
 res_h0 <- results(dds, name = "h0_rest_vs_aim")
 res_h0 <- data.frame(res_h0)
 res_h0 <- res_h0[complete.cases(res_h0),]
 write.csv(res_h0, file = "res_h0.csv")
            
            
 #### finding genes that are low at the beginning but low at the end
 ### creating the violinplot
 vln_norm <- nSingle[,rownames(colData_DE)]
 de_gene <- res_h0[res_h0$padj < 0.01,]
 pos <- de_gene[de_gene$log2FoldChange > 2,]
            
            
 ## violinplot
 colData_DE$gene <- log10(vln_norm["Tb927.1.290",] + 1)
 ggplot(colData_DE, aes(x = time, y = gene, color = time)) +
       geom_violin() +
       geom_jitter() +
       theme_classic()
            
            
            
##h4
dds <- DESeqDataSetFromMatrix(countData = single_table,
       colData = condition,
       design= ~h4)
dds <- DESeq(dds)
resultsNames(dds)
res_h4 <- results(dds, name = "h4_rest_vs_aim")
res_h4 <- data.frame(res_h4)
res_h4 <- res_h4[complete.cases(res_h4),]
            
write.csv(res_h4, file = "res_h4.csv")
##h12
dds <- DESeqDataSetFromMatrix(countData = single_table,
       colData = condition,
       design= ~h12)
dds <- DESeq(dds)
resultsNames(dds)
res_h12 <- results(dds, name = "h12_rest_vs_aim")
res_h12 <- data.frame(res_h12)
res_h12 <- res_h12[complete.cases(res_h12),]
            
          
write.csv(res_h12, file = "res_h12.csv")
##h24
dds <- DESeqDataSetFromMatrix(countData = single_table,
       colData = condition,
       design= ~h24)
dds <- DESeq(dds)
resultsNames(dds)
res_h24 <- results(dds, name = "h24_rest_vs_aim")
res_h24 <- data.frame(res_h24)
res_h24 <- res_h24[complete.cases(res_h24),]
              
write.csv(res_h24, file = "res_h24.csv")
            
            
##d7
dds <- DESeqDataSetFromMatrix(countData = single_table,
       colData = condition,
       design= ~d7)
dds <- DESeq(dds)
resultsNames(dds)
res_d7 <- results(dds, name = "d7_rest_vs_aim")
res_d7 <- data.frame(res_d7)
res_d7 <- res_d7[complete.cases(res_d7),]
           
write.csv(res_d7, file = "res_d7.csv") 
            
            
            
res_h0_hm <- res_h0[res_h0$padj < 0.01,]
res_h0_hm <- res_h0_hm[res_h0_hm$log2FoldChange < 2 | res_h0_hm$log2FoldChange > 2,]
      
res_h4_hm <- res_h4[res_h4$padj < 0.01,]
res_h4_hm <- res_h4_hm[res_h4_hm$log2FoldChange < 2 | res_h4_hm$log2FoldChange > 2,] 
            
res_h12_hm <- res_h12[res_h12$padj < 0.01,]
res_h12_hm <- res_h12_hm[res_h12_hm$log2FoldChange < 2 | res_h12_hm$log2FoldChange > 2,] 
            
res_h24_hm <- res_h24[res_h24$padj < 0.01,]
res_h24_hm <- res_h24_hm[res_h24_hm$log2FoldChange < 2 | res_h24_hm$log2FoldChange > 2,]
            
res_d7_hm <- res_d7[res_d7$padj < 0.01,]
res_d7_hm <- res_d7_hm[res_d7_hm$log2FoldChange < 2 | res_d7_hm$log2FoldChange > 2,]
        
gene_selected <- read.csv(file = "All RBPs for Heatmap.txt", header = FALSE)
hm_gene <- as.character(gene_selected$V1)
        
#hm_gene <- unique(c(rownames(res_h0_hm), rownames(res_h4_hm), rownames(res_h12_hm), rownames(res_h24_hm), rownames(res_d7_hm)))
norm <- nSingle[hm_gene, ht_order]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
group <- factor(colData_DE$time)
ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "24h" = "#F8766D", "4h" = "#6a3d9a")))
pdf("./output/heatmap_comparision5.pdf", width = 14, height = 10)
p <- Heatmap(scaled, cluster_columns = TRUE, cluster_rows = TRUE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 2), show_column_names = TRUE)
print(p)
dev.off()
  
  
  
  ########################################################
  ########################################################
  res_h0_hm <- res_h0[res_h0$padj < 0.01,]
  res_h0_hm <- res_h0_hm[res_h0_hm$log2FoldChange < 2,]
  
  res_h4_hm <- res_h4[res_h4$padj < 0.01,]
  res_h4_hm <- res_h4_hm[res_h4_hm$log2FoldChange < 2,] 
  
  res_h12_hm <- res_h12[res_h12$padj < 0.01,]
  res_h12_hm <- res_h12_hm[res_h12_hm$log2FoldChange < 2,] 
  
  res_h24_hm <- res_h24[res_h24$padj < 0.01,]
  res_h24_hm <- res_h24_hm[res_h24_hm$log2FoldChange < 2,] 
  
  res_d7_hm <- res_d7[res_d7$padj < 0.01,]
  res_d7_hm <- res_d7_hm[res_d7_hm$log2FoldChange < 2,]
  
  
  
  hm_gene <- unique(c(rownames(res_h0_hm), rownames(res_h4_hm), rownames(res_h12_hm), rownames(res_h24_hm), rownames(res_d7_hm)))
  
  
  norm <- nSingle[hm_gene, ht_order]
  scaled <- t(apply(norm, 1, scale))
  colnames(scaled) <- colnames(norm)
  #group <- factor(colData_DE$time)
  #ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "24h" = "#F8766D", "4h" = "#6a3d9a")))
  pdf("./output/heatmap_comparision2.pdf", width = 14, height = 10)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 2), show_column_names = TRUE)
  print(p)
    
  
  ################ Performing DE analysis ... early vs late ###################
  #colData_DE <- colData_DE[!colData_DE$time == "7d",]
  #cd <- single_table[,as.character(colData_DE$col) ]
  pair <- factor(gsub("4h","early",gsub("12h","late",gsub("24h", "late", gsub("0h", "early", colData_DE$time)))))
  colData_DE <- cbind(colData_DE,pair)
  
  
  group <- factor(colData_DE$pair)
  names(group) <- rownames(o.ifm)
  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
  de_early <- ediff[order(ediff$Z, decreasing  =  TRUE), ]
  write.csv(de_early, file = "de_early_vs_late.csv")
  
  de_ht_pos <- de_early[de_early$cZ > 1.96, ]
  de_ht_neg <- de_early[de_early$cZ < -1.96, ]
  de_sig <- rbind(de_ht_pos, de_ht_neg)
  de_sig <- rownames(de_sig)
  de_sig <- unique(de_sig)
  
  
  h4 <- as.character(colData_DE[colData_DE$time == "4h", ]$col)
  h0 <- as.character(colData_DE[colData_DE$time == "0h", ]$col)
  h12 <- as.character(colData_DE[colData_DE$time == "12h", ]$col)
  h24 <- as.character(colData_DE[colData_DE$time == "24h", ]$col)
  d7 <- as.character(colData_DE[colData_DE$time == "7d", ]$col)
  ht_order <- c(h0, h4, h12, h24)
  
  
  norm <- nSingle[de_sig, ht_order]
  scaled <- t(apply(norm, 1, scale))
  colnames(scaled) <- colnames(norm)
  #group <- factor(colData_DE$time)
  #ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "24h" = "#F8766D", "4h" = "#6a3d9a")))
  pdf("./output/heatmap_early_vs_late.pdf", width = 14, height = 7)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 4), show_column_names = TRUE)
  print(p)
  dev.off()

############################################################################
#################### DE early vs late DESeq ################################
  ######################################################################
  pair <- factor(gsub("4h","early",gsub("12h","late",gsub("24h", "late", gsub("0h", "early", colData_DE$time))))) 
  colData_DE <- cbind(colData_DE,pair)
  
  condition <- colData_DE
  dds <- DESeqDataSetFromMatrix(countData = single_table,
                                colData = condition,
                                design= ~pair)
  dds <- DESeq(dds)
  resultsNames(dds)
  res_pair <- results(dds, name = "pair_late_vs_early")
  res_pair <- data.frame(res_pair)
  res_pair <- res_pair[complete.cases(res_pair),]
  write.csv(res_pair, file = "res_early_vs_late.csv")
  
  res_pair_hm <- res_pair[res_pair$padj < 0.01,]
  pos <- res_pair_hm[res_pair_hm$log2FoldChange > 5,] 
  neg <- res_pair_hm[res_pair_hm$log2FoldChange < 5,]
  
  hm_gene <- c(rownames(pos), rownames(neg))
  
  
  norm <- nSingle[hm_gene, ht_order]
  scaled <- t(apply(norm, 1, scale))
  colnames(scaled) <- colnames(norm)
  #group <- factor(colData_DE$time)
  #ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("0h" = "#619cff", "12h" = "#00BA38", "24h" = "#F8766D", "4h" = "#6a3d9a")))
  pdf("./output/heatmap_early_vs_late_2.pdf", width = 14, height = 7)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), row_names_gp = gpar(fontsize = 2), show_column_names = TRUE)
  print(p)
  dev.off()
  
##########################################################################
############################# Pairwise comparision ########################
##############################################################################
  ### 7d vs 0h $ 7d vs 24h
  
  
  
  
  
#############################################################################  
    ### PCA on 50 cells
  colData_filt <- colData[colData$nCell == "50", ]
  
  single_table <- unique_table50[c(1:57,157:11807), as.character(colData_filt$col)]
  pos <- grep("mVAT",rownames(single_table))  
  pos <- c(pos,11190)
  single_table <- single_table[-pos, ]
  single_table <- single_table[rowSums(single_table) > 20, ]
  
  
  ### Normalizing data
  sfSingle <- estimateSizeFactorsForMatrix(single_table)
  nSingle <- t(t(single_table) / sfSingle)
  
  ##########################################################################
  ## Selecting the top variable genes and log transformation
  topVarGenesSingle <- head(order(rowVars(nSingle), decreasing = TRUE),1000)
  nSingle_fin <- log(nSingle[topVarGenesSingle,] + 1)
  ### Performing PCA (single cell)
  res.pca <- prcomp(t(nSingle_fin), scale = FALSE)
  ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], time = colData_filt$time), aes(x = PC1, y = PC2, color = time)) +
    geom_jitter(size = 2) +
    theme_classic() +
    theme() +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed")
  
  ggsave("./output/pca50.png", width = 5, height = 5)
  
  ### running DE analysis
  colData_DE <- colData_filt
  condition <- factor(c("one", "one", "one", "one", "one", "one", "two", "two", "two", "two", "two", "two", "two", "two", "two", "two", "two"))
  
  cts <- as.matrix(single_table[,as.character(colData_DE$col) ])
  #cts <- cts[rowSums(cts) > 0, ]
  #cts <- cts[rowSums(cts > 0) > 5, ]
  
  col_data <- data.frame(colData_DE, condition = factor(condition))
  rownames(col_data) <- col_data$col
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = col_data,
                                design= ~condition)
  resultsNames(dds)  
  
  
  
  
  
  
