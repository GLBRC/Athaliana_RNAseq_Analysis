library(pheatmap)
library(dendextend)
library(dichromat)
library(factoextra)
library(cluster)
library(dplyr)
library(UpSetR)
library(grid)
setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/XCC_Exp/")

xcc_data <- read.table(file = "XCC_FPKM_values_FC_pAdj_SIG.txt", header = T, sep = "\t", row.names = 1)
head(xcc_data)
summary(xcc_data)
xcc_data_noNA <- na.omit(xcc_data)

xcc_data_log2FC <- xcc_data_noNA[ , c("Col_Mock_vs_CSLC_Mock_Log2_FC", "Col_XCC_vs_CSLC_XCC_Log2_FC", "Col_Mock_vs_Col_XCC_Log2_FC", "CSLC_Mock_vs_CSLC_XCC_Log2_FC")]
xcc_sig_data <- xcc_data_noNA[ , c("CSLC_Mock_vs_CSLC_XCC_SIG", "Col_Mock_vs_Col_XCC_SIG", "Col_XCC_vs_CSLC_XCC_SIG", "Col_Mock_vs_CSLC_Mock_SIG")]
xcc_ann_colors <- list(CSLC_Mock_vs_CSLC_XCC_SIG=c("Yes" = "darkgreen",
                                                 "no" = "grey"),
                       Col_Mock_vs_Col_XCC_SIG=c("Yes" = "darkgreen",
                                                    "no" = "grey"),
                       Col_XCC_vs_CSLC_XCC_SIG=c("Yes" = "darkgreen",
                                                   "no" = "grey"),
                       Col_Mock_vs_CSLC_Mock_SIG=c("Yes" = "darkgreen",
                                                 "no" = "grey")
)
breaks = seq(-10,10,1)

pheatmap(xcc_data_log2FC, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "average",
         main = "WT (col) and Mutant (CSLC)\nplus Bacteria (XCC) minus Bacteria (Mock)\nLog2 DEseq2 Fold Change",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("gold", "green4", "purple4"))(20),
         annotation_row = xcc_sig_data,
         annotation_colors = xcc_ann_colors,
         breaks = breaks
)

only_sig_xcc_data_noNA <- filter(xcc_data_noNA, Col_Mock_vs_Col_XCC_SIG == "Yes" 
                                 | Col_Mock_vs_CSLC_Mock_SIG == "Yes" 
                                 | CSLC_Mock_vs_CSLC_XCC_SIG == "Yes" 
                                 | Col_XCC_vs_CSLC_XCC_SIG == "Yes")
only_sig_xcc_data_noNA_logFC <- only_sig_xcc_data_noNA[ , c("Col_XCC_vs_CSLC_XCC_Log2_FC", "Col_Mock_vs_CSLC_Mock_Log2_FC", "CSLC_Mock_vs_CSLC_XCC_Log2_FC", "Col_Mock_vs_Col_XCC_Log2_FC")]
breaks = seq(-10,10,1)

pheatmap(only_sig_xcc_data_noNA_logFC, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "average",
         main = "WT (col) and Mutant (CSLC)\nplus Bacteria (XCC) minus Bacteria (Mock)\nLog2 DEseq2 Fold Change SIG in at least one Comparison",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = xcc_sig_data,
         annotation_colors = xcc_ann_colors,
         gaps_col=c(2), breaks = breaks
)

pheatmap(only_sig_xcc_data_noNA_logFC, kmeans_k = 6)

only_sig_xcc_col_mockVSxcc <- filter(xcc_data_noNA, Col_Mock_vs_Col_XCC_SIG == "Yes", (Col_Mock_vs_Col_XCC_Log2_FC >= 1 | Col_Mock_vs_Col_XCC_Log2_FC <= -1))
only_sig_xcc_cslc_mockVSxcc <- filter(xcc_data_noNA, CSLC_Mock_vs_CSLC_XCC_SIG == "Yes", (CSLC_Mock_vs_CSLC_XCC_Log2_FC >= 1 | CSLC_Mock_vs_CSLC_XCC_Log2_FC <= -1))
only_sig_xcc_mock_colVScslc <- filter(xcc_data_noNA, Col_Mock_vs_CSLC_Mock_SIG == "Yes", (Col_Mock_vs_CSLC_Mock_Log2_FC >= 1 | Col_Mock_vs_CSLC_Mock_Log2_FC <= -1))
only_sig_xcc_xcc_colVScslc <- filter(xcc_data_noNA, Col_XCC_vs_CSLC_XCC_SIG == "Yes", (Col_XCC_vs_CSLC_XCC_Log2_FC >=1 | Col_XCC_vs_CSLC_XCC_Log2_FC <= -1))

only_sig_xcc_col_mockVSxcc_genes <- row.names(only_sig_xcc_col_mockVSxcc)
only_sig_xcc_cslc_mockVSxcc_genes <- row.names(only_sig_xcc_cslc_mockVSxcc)
only_sig_xcc_mock_colVScslc_genes <- row.names(only_sig_xcc_mock_colVScslc)
only_sig_xcc_xcc_colVScslc_genes <- row.names(only_sig_xcc_xcc_colVScslc)

geneList <- list(col_mock_vs_xcc = row.names(only_sig_xcc_col_mockVSxcc),
                 cslc_mock_vs_xcc = row.names(only_sig_xcc_cslc_mockVSxcc),
                 mock_col_vs_cslc = row.names(only_sig_xcc_mock_colVScslc),
                 xcc_col_vs_cslc = row.names(only_sig_xcc_xcc_colVScslc))
head(fromList(geneList))

#Get the lists from the upset plot and intersections

interactionList <- fromList(geneList)

uniqueGeneNames <- data.frame(gene=unique(unlist(geneList)))
head(uniqueGeneNames)

GeneComparison_List <- lapply(geneList, function(x){
  data.frame(gene=x)
}) %>%  bind_rows(.id = "path")
head(GeneComparison_List)

geneList_GeneName <- lapply(uniqueGeneNames$gene,function(x){
  intersection <- GeneComparison_List %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  data.frame(gene = x, int = intersection)
}) %>% bind_rows()
head(geneList_GeneName)
dim(geneList_GeneName)

geneList_GeneName %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

write.table(geneList_GeneName, file = "Intersection_of_FC_Genes_UpSet_Graph_20250929.txt", quote = F, sep = "\t")

upset(fromList(geneList),
      sets = c("col_mock_vs_xcc", 
               "cslc_mock_vs_xcc", 
               "mock_col_vs_cslc", 
               "xcc_col_vs_cslc"),
      order.by = "freq",
      sets.bar.color = "grey10", line.size = 1, point.size = 3,
      text.scale = 1.2, matrix.color = "darkgreen", main.bar.color = "grey10",
      sets.x.label = "DEGs per Comparison",
      mainbar.y.label = "Number of Total DEGs",
      mainbar.y.max = 8500,
      set_size.scale_max = 10000
      )
grid.text("ALL FC DEGs", x = 0.15, y = 0.95)

only_sig_xcc_col_mockVSxcc_increase <- filter(xcc_data_noNA, Col_Mock_vs_Col_XCC_SIG == "Yes", (Col_Mock_vs_Col_XCC_Log2_FC >= 1))
only_sig_xcc_cslc_mockVSxcc_increase <- filter(xcc_data_noNA, CSLC_Mock_vs_CSLC_XCC_SIG == "Yes", (CSLC_Mock_vs_CSLC_XCC_Log2_FC >= 1))
only_sig_xcc_mock_colVScslc_increase <- filter(xcc_data_noNA, Col_Mock_vs_CSLC_Mock_SIG == "Yes", (Col_Mock_vs_CSLC_Mock_Log2_FC >= 1))
only_sig_xcc_xcc_colVScslc_increase <- filter(xcc_data_noNA, Col_XCC_vs_CSLC_XCC_SIG == "Yes", (Col_XCC_vs_CSLC_XCC_Log2_FC >=1))

only_sig_xcc_col_mockVSxcc_increase_genes <- row.names(only_sig_xcc_col_mockVSxcc_increase)
only_sig_xcc_cslc_mockVSxcC_increasec_genes <- row.names(only_sig_xcc_cslc_mockVSxcc_increase)
only_sig_xcc_mock_colVScslc_increase_genes <- row.names(only_sig_xcc_mock_colVScslc_increase)
only_sig_xcc_xcc_colVScslc_increase_genes <- row.names(only_sig_xcc_xcc_colVScslc_increase)

geneList <- list(col_mock_vs_xcc = row.names(only_sig_xcc_col_mockVSxcc_increase),
                 cslc_mock_vs_xcc = row.names(only_sig_xcc_cslc_mockVSxcc_increase),
                 mock_col_vs_cslc = row.names(only_sig_xcc_mock_colVScslc_increase),
                 xcc_col_vs_cslc = row.names(only_sig_xcc_xcc_colVScslc_increase))
head(fromList(geneList))

#Get the lists from the upset plot and intersections

interactionList <- fromList(geneList)

uniqueGeneNames <- data.frame(gene=unique(unlist(geneList)))
head(uniqueGeneNames)

GeneComparison_List <- lapply(geneList, function(x){
  data.frame(gene=x)
}) %>%  bind_rows(.id = "path")
head(GeneComparison_List)

geneList_GeneName <- lapply(uniqueGeneNames$gene,function(x){
  intersection <- GeneComparison_List %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  data.frame(gene = x, int = intersection)
}) %>% bind_rows()
head(geneList_GeneName)
dim(geneList_GeneName)

geneList_GeneName %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

write.table(geneList_GeneName, file = "Intersection_of_INCREASE_FC_Genes_UpSet_Graph_20250929.txt", quote = F, sep = "\t")


upset(fromList(geneList),
      sets = c("col_mock_vs_xcc", 
               "cslc_mock_vs_xcc", 
               "mock_col_vs_cslc", 
               "xcc_col_vs_cslc"),
      order.by = "freq",
      sets.bar.color = "grey10", line.size = 1, point.size = 3,
      text.scale = 1.2, matrix.color = "darkgreen", main.bar.color = "grey10",
      sets.x.label = "DEGs per Comparison",
      mainbar.y.label = "Number of Total DEGs",
      mainbar.y.max = 4000,
      set_size.scale_max = 4500
)
grid.text("Positive FC DEGs", x = 0.15, y = 0.95)

only_sig_xcc_col_mockVSxcc_decrease <- filter(xcc_data_noNA, Col_Mock_vs_Col_XCC_SIG == "Yes", (Col_Mock_vs_Col_XCC_Log2_FC <= 1))
only_sig_xcc_cslc_mockVSxcc_decrease <- filter(xcc_data_noNA, CSLC_Mock_vs_CSLC_XCC_SIG == "Yes", (CSLC_Mock_vs_CSLC_XCC_Log2_FC <= 1))
only_sig_xcc_mock_colVScslc_decrease <- filter(xcc_data_noNA, Col_Mock_vs_CSLC_Mock_SIG == "Yes", (Col_Mock_vs_CSLC_Mock_Log2_FC <= 1))
only_sig_xcc_xcc_colVScslc_decrease <- filter(xcc_data_noNA, Col_XCC_vs_CSLC_XCC_SIG == "Yes", (Col_XCC_vs_CSLC_XCC_Log2_FC <=1))

only_sig_xcc_col_mockVSxcc_decrease_genes <- row.names(only_sig_xcc_col_mockVSxcc_decrease)
only_sig_xcc_cslc_mockVSxcc_decrease_genes <- row.names(only_sig_xcc_cslc_mockVSxcc_decrease)
only_sig_xcc_mock_colVScslc_decrease_genes <- row.names(only_sig_xcc_mock_colVScslc_decrease)
only_sig_xcc_xcc_colVScslc_decrease_genes <- row.names(only_sig_xcc_xcc_colVScslc_decrease)

geneList <- list(col_mock_vs_xcc = row.names(only_sig_xcc_col_mockVSxcc_decrease),
                 cslc_mock_vs_xcc = row.names(only_sig_xcc_cslc_mockVSxcc_decrease),
                 mock_col_vs_cslc = row.names(only_sig_xcc_mock_colVScslc_decrease),
                 xcc_col_vs_cslc = row.names(only_sig_xcc_xcc_colVScslc_decrease))
head(fromList(geneList))

#Get the lists from the upset plot and intersections

interactionList <- fromList(geneList)

uniqueGeneNames <- data.frame(gene=unique(unlist(geneList)))
head(uniqueGeneNames)

GeneComparison_List <- lapply(geneList, function(x){
  data.frame(gene=x)
}) %>%  bind_rows(.id = "path")
head(GeneComparison_List)

geneList_GeneName <- lapply(uniqueGeneNames$gene,function(x){
  intersection <- GeneComparison_List %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  data.frame(gene = x, int = intersection)
}) %>% bind_rows()
head(geneList_GeneName)
dim(geneList_GeneName)

geneList_GeneName %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

write.table(geneList_GeneName, file = "Intersection_of_DECREASE_FC_Genes_UpSet_Graph_20250929.txt", quote = F, sep = "\t")

upset(fromList(geneList),
      sets = c("col_mock_vs_xcc", 
               "cslc_mock_vs_xcc", 
               "mock_col_vs_cslc", 
               "xcc_col_vs_cslc"),
      order.by = "freq",
      sets.bar.color = "grey10", line.size = 1, point.size = 3,
      text.scale = 1.2, matrix.color = "darkgreen", main.bar.color = "grey10",
      sets.x.label = "DEGs per Comparison",
      mainbar.y.label = "Number of Total DEGs",
      mainbar.y.max = 5500,
      set_size.scale_max = 5500
)
grid.text("Negative FC DEGs", x = 0.15, y = 0.95)

colnames(xcc_data_noNA)

xcc_data_noNA_FPKM_means <- xcc_data_noNA[ , c("Col_Mock_AVG", "Mut_Mock_AVG", "Col_XCC_AVG", "Mut_XCC_AVG")]
xcc_data_noNA_FPKM_means_log2 <- log(xcc_data_noNA_FPKM_means,2)

Col_Mock_AVG_values <- xcc_data_noNA_FPKM_means_log2$Col_Mock_AVG
Col_Mock_AVG_values_noInf <- ifelse(is.finite(Col_Mock_AVG_values), Col_Mock_AVG_values, 0)
Mut_Mock_AVG_values <- xcc_data_noNA_FPKM_means_log2$Mut_Mock_AVG
Mut_Mock_AVG_values_noInf <- ifelse(is.finite(Mut_Mock_AVG_values), Mut_Mock_AVG_values, 0)
Col_XCC_AVG_values_noInf <- ifelse(is.finite(xcc_data_noNA_FPKM_means_log2$Col_XCC_AVG), xcc_data_noNA_FPKM_means_log2$Col_XCC_AVG, 0)
Mut_XCC_AVG_values_noInf <- ifelse(is.finite(xcc_data_noNA_FPKM_means_log2$Mut_XCC_AVG), xcc_data_noNA_FPKM_means_log2$Mut_XCC_AVG, 0)

xcc_data_noNA_FPKM_means_log2_noInf <- xcc_data_noNA_FPKM_means_log2[!is.infinite(rowSums(xcc_data_noNA_FPKM_means_log2)),]

#### Z-score ####
colNames <- colnames(xcc_data_noNA_FPKM_means_log2_noInf)
samples <- setdiff(names(xcc_data_noNA_FPKM_means_log2_noInf), "geneID")
xcc_FPKM_matrix <- as.matrix(xcc_data_noNA_FPKM_means_log2_noInf)
rowMeans <- rowMeans(xcc_FPKM_matrix)
rowStDev <- apply(xcc_FPKM_matrix, 1, sd)
rowStDev2 <- apply(xcc_FPKM_matrix, 1, sd); rowStDev2[is.na(rowStDev2) | rowStDev2 == 0] <- 1
zscore <- sweep(sweep(xcc_FPKM_matrix, 1, rowMeans, "-"), 1, rowStDev2, "/")
xcc_FPKM_log2_noInf_Zscore <- zscore

write.table(xcc_FPKM_log2_noInf_Zscore, file = "XCC_FPKM_log2_noInf_zscore.txt")

fviz_nbclust(xcc_FPKM_log2_noInf_Zscore, kmeans, method = "wss")
fviz_nbclust(xcc_FPKM_log2_noInf_Zscore, kmeans, method = "silhouette")

colNames_notNorm <- colnames(xcc_data_noNA_FPKM_means)
samples_nonNorm <- setdiff(names(xcc_data_noNA_FPKM_means), "geneID")
xcc_FPKM_matrix_nonNorm <- as.matrix(xcc_data_noNA_FPKM_means)
rowMeans_nonNorm <- rowMeans(xcc_FPKM_matrix_nonNorm)
rowStDev_nonNorm <- apply(xcc_FPKM_matrix_nonNorm, 1, sd)
rowStDev2_nonNorm <- apply(xcc_FPKM_matrix_nonNorm, 1, sd); rowStDev2_nonNorm[is.na(rowStDev2_nonNorm) | rowStDev2_nonNorm == 0] <- 1
zscore_nonNorm <- sweep(sweep(xcc_FPKM_matrix_nonNorm, 1, rowMeans_nonNorm, "-"), 1, rowStDev2_nonNorm, "/")
xcc_FPKM_nonNorm_Zscore <- zscore_nonNorm

fviz_nbclust(xcc_FPKM_nonNorm_Zscore, kmeans, method = "wss")
fviz_nbclust(xcc_FPKM_nonNorm_Zscore, kmeans, method = "silhouette")

##### Zscore Kmeans 2 clusters #####

pheatmap(xcc_FPKM_log2_noInf_Zscore, kmeans_k = 2, cluster_cols = F)
xcc_FPKM_log2_noInf_Zscore_out_2_clusters <- pheatmap(xcc_FPKM_log2_noInf_Zscore, kmeans_k = 2, cluster_cols = F, display_numbers = T)
clusters6_xcc_FPKM_log2_noInf_Zscore_out_2clusters <- xcc_FPKM_log2_noInf_Zscore_out_2_clusters$kmeans$cluster
write.table(clusters6_xcc_FPKM_log2_noInf_Zscore_out_2clusters, file = "Kmeans_2_FPKM_log2_noInf_Zscore_Clusters_XCC.txt")

xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters <- read.table(file = "XCC_FPKM_Log2_noInf_ZScore_cluster_2Clusters_sorted.txt", header = T, sep = "\t", row.names = 1)
xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters_values <- xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters[ , c("Col_Mock_AVG", "Mut_Mock_AVG", "Col_XCC_AVG", "Mut_XCC_AVG")]
xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters_kMeanClusters <- read.table(file = "Kmeans_2_FPKM_log2_noInf_Zscore_Clusters_XCC.txt", header = T, sep = "\t", row.names = 1)
xcc_annotation_col <- data.frame(Cluster = factor(xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters_kMeanClusters$Zscore_Kmeans_2_Cluster))
rownames(xcc_annotation_col) <- rownames(xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters_kMeanClusters)
colnames(xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters_values) <- c("Col_Mock_Avg", "CSLC_Mock_Avg", "Col_XCC_Avg", "CSLC_XCC_Avg")

kmeans_colors <- list(Kmeans_2_Cluster=c("1" = "blue",
                                         "2" = "red"
                                         ))

pheatmap(xcc_FPKM_log2_Zscore_kmeans_sorted_2clusters_values, cluster_rows = F, cluster_cols = F,
         main = "XCC Experiment Z-score Log2 KMeans 2 Clusters",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = xcc_annotation_col,
         annotation_colors = kmeans_colors,
         gaps_col=c(1,2),
         breaks=seq(-2,2, length.out=21), legend = T
)

library(ggplot2)
library(GGally)
library(dplyr)
library(patchwork)
mean_per_cluster_zscore <- read.table(file = "for_plotting_FPKM_average_zScore_per_cluster_2Clusters.txt", header = T, sep = "\t")
cluster1_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_1, group = 1)) +
  geom_line( color = "blue", lwd = 1) +
  geom_point( color = "blue", shape = 21, fill = "blue", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 1 (n=9609)")
cluster1_plot

cluster2_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_2, group = 1)) +
  geom_line( color = "red", lwd = 1) +
  geom_point( color = "red", shape = 21, fill = "red", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 2 (n=8214)")
cluster2_plot

cluster1_plot + cluster2_plot + plot_annotation(title = "Kmeans Z-Score 2 Clusters Graphs - XCC 2025-09-15")
ggsave("XCC_Kmeans_Zscore_expressionClusters_2Clusters_20250915.pdf")
ggsave("XCC_Kmeans_Zscore_expressionClusters_2Clusters_20250915.eps")

setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/XCC_Exp")
goBPdata <- read.table(file = "XCC_topTenGO_BP_2_Clusters_Zscore_negLogPval.txt", header = T, row.names = 1, sep = "\t")
colnames(goBPdata) <- c("Cluster1", "Cluster2")
cluster_list <- c("Cluster1", "Cluster2")
annotation_col <- data.frame(Cluster = factor(colnames(goBPdata)))
rownames(annotation_col) <- colnames(goBPdata)
ann_colors_col <- list(Cluster=c("Cluster1" = "blue",
                                 "Cluster2" = "red"
                                 ))
goBPplot <- pheatmap(goBPdata, cluster_rows = T, cluster_cols = F,
                     clustering_distance_rows = 'euclidean',
                     clustering_method = "complete",
                     main = "        GO:BP Top 10 per 2 Zscore Clusters XCC\n(neg log10 adj P-value)",
                     show_rownames = T, show_colnames = T,
                     border_color = "grey50", fontsize_col = 12,
                     color = colorRampPalette(c("white", "lightblue1", "darkblue"))(60),
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors_col,
                     breaks=seq(0,60, length.out=61),
                     labels_col=c("1","2"),
                     cutree_rows = 6, cellwidth = 20, cellheight = 10
)

pdf("XCC_GO_BP_top10_per_Cluster_2clusters_Zscore_Heatmap_20290918.pdf", height = 8, width = 10)
goBPplot
dev.off()



###### Zscore Kmeans 6 clusters #####
pheatmap(xcc_FPKM_log2_noInf_Zscore, kmeans_k = 6, cluster_cols = F)
xcc_FPKM_log2_noInf_Zscore_out <- pheatmap(xcc_FPKM_log2_noInf_Zscore, kmeans_k = 6, cluster_cols = F, display_numbers = T)
clusters6_xcc_FPKM_log2_noInf_Zscore_out <- xcc_FPKM_log2_noInf_Zscore_out$kmeans$cluster
write.table(clusters6_xcc_FPKM_log2_noInf_Zscore_out, file = "Kmeans_6_FPKM_log2_noInf_Zscore_Clusters_XCC.txt")

xcc_FPKM_log2_Zscore_kmeans_sorted <- read.table(file = "XCC_FPKM_Log2_noInf_ZScore_cluster_sorted.txt", header = T, sep = "\t", row.names = 1)
xcc_FPKM_log2_Zscore_kmeans_sorted_values <- xcc_FPKM_log2_Zscore_kmeans_sorted[ , c("Col_Mock_AVG", "Mut_Mock_AVG", "Col_XCC_AVG", "Mut_XCC_AVG")]
xcc_FPKM_log2_Zscore_kmeans_sorted_kMeanClusters <- read.table(file = "Kmeans_6_FPKM_log2_noInf_Zscore_Clusters_XCC.txt", header = T, sep = "\t", row.names = 1)
xcc_annotation_col <- data.frame(Cluster = factor(xcc_FPKM_log2_Zscore_kmeans_sorted_kMeanClusters$Zscore_Kmeans_6_Cluster))
rownames(xcc_annotation_col) <- rownames(xcc_FPKM_log2_Zscore_kmeans_sorted_kMeanClusters)
colnames(xcc_FPKM_log2_Zscore_kmeans_sorted_values) <- c("Col_Mock_Avg", "CSLC_Mock_Avg", "Col_XCC_Avg", "CSLC_XCC_Avg")

kmeans_colors <- list(Kmeans_6_Cluster=c("1" = "blue",
                                         "2" = "red",
                                         "3" = "yellow",
                                         "4" = "orange",
                                         "5" = "green",
                                         "6" = "purple"))

pheatmap(xcc_FPKM_log2_Zscore_kmeans_sorted_values, cluster_rows = F, cluster_cols = F,
         main = "XCC Experiment Z-score Log2 KMean Clusters",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = xcc_annotation_col,
         annotation_colors = kmeans_colors,
         gaps_col=c(1,2),
         breaks=seq(-2,2, length.out=21), legend = T
)

library(ggplot2)
library(GGally)
library(dplyr)
library(patchwork)
mean_per_cluster_zscore <- read.table(file = "for_plotting_FPKM_average_zScore_per_cluster.txt", header = T, sep = "\t")
cluster1_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_1, group = 1)) +
  geom_line( color = "blue", lwd = 1) +
  geom_point( color = "blue", shape = 21, fill = "blue", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 1 (n=866)")
cluster1_plot

cluster2_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_2, group = 1)) +
  geom_line( color = "red", lwd = 1) +
  geom_point( color = "red", shape = 21, fill = "red", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 2 (n=461)")
cluster2_plot

cluster3_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_3, group = 1)) +
  geom_line( color = "yellow", lwd = 1) +
  geom_point( color = "yellow", shape = 21, fill = "yellow", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 3 (n=991)")
cluster3_plot

cluster4_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_4, group = 1)) +
  geom_line( color = "orange", lwd = 1) +
  geom_point( color = "orange", shape = 21, fill = "orange", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 4 (n=759)")
cluster4_plot

cluster5_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_5, group = 1)) +
  geom_line( color = "green", lwd = 1) +
  geom_point( color = "green", shape = 21, fill = "green", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 5 (n=6797)")
cluster5_plot

cluster6_plot <- ggplot(mean_per_cluster_zscore, aes(x = Sample, y = Cluster_6, group = 1)) +
  geom_line( color = "purple", lwd = 1) +
  geom_point( color = "purple", shape = 21, fill = "purple", size = 5) +
  ylim(-1.1, 1) +
  labs(x = "Sample Condition", y = "Z-score Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 6 (n=7949)")
cluster6_plot

cluster1_plot + cluster2_plot + cluster3_plot + cluster4_plot + cluster5_plot + cluster6_plot + plot_annotation(title = "Kmeans Z-Score Cluster Graphs - XCC 2025-09-15")
ggsave("XCC_Kmeans_Zscore_expressionClusters_20250915.pdf")
ggsave("XCC_Kmeans_Zscore_expressionClusters_20250915.eps")

##### Mean Center #######

rowmean_log2 <- apply(xcc_data_noNA_FPKM_means_log2,1,mean)
xcc_data_noNA_FPKM_MeanCenter_log2 <- xcc_data_noNA_FPKM_means_log2 - rowmean_log2

write.table(xcc_data_noNA_FPKM_MeanCenter_log2, file = "XCC_MeanCenter_FPKM_Log2.txt")

rowmean <- apply(xcc_data_noNA_FPKM_means,1,mean)
xcc_data_noNA_FPKM_MeanCenter <- xcc_data_noNA_FPKM_means - rowmean
xcc_data_noNA_FPKM_MeanCenter_thenLog2 <- log(xcc_data_noNA_FPKM_MeanCenter,2)

write.table(xcc_data_noNA_FPKM_MeanCenter_thenLog2, file = "XCC_Linear_MeanCenter_FPKM_Log2.txt")

xcc_data_noNA_FPKM_MeanCenter_Log2_data <- read.table(file = "XCC_MeanCenter_FPKM_Log2.txt", header = T, row.names = 1, sep = "\t")
xcc_data_noNA_FPKM_MeanCenter_Log2_data_noNA <- na.omit(xcc_data_noNA_FPKM_MeanCenter_Log2_data)

fviz_nbclust(xcc_data_noNA_FPKM_MeanCenter_Log2_data_noNA, kmeans, method = "wss")
fviz_nbclust(xcc_data_noNA_FPKM_MeanCenter_Log2_data_noNA, kmeans, method = "silhouette")

pheatmap(xcc_data_noNA_FPKM_MeanCenter_Log2_data_noNA, kmeans_k = 6, cluster_cols = F)
xcc_data_noNA_FPKM_MeanCenter_Log2_data_noNA_out <- pheatmap(xcc_data_noNA_FPKM_MeanCenter_Log2_data_noNA, kmeans_k = 6, cluster_cols = F, display_numbers = T)
clusters6_meanCenter_XCC_log2 <- xcc_data_noNA_FPKM_MeanCenter_Log2_data_noNA_out$kmeans$cluster
write.table(clusters6_meanCenter_XCC_log2, file = "Kmeans_6_meanCenter_Clusters_XCC.txt")

xcc_meanCenter_kmeans_sorted <- read.table(file = "XCC_MeanCenter_FPKM_Log2_cluster_sorted.txt", header = T, sep = "\t", row.names = 1)
xcc_meanCenter_kmeans_sorted_values <- xcc_meanCenter_kmeans_sorted[ , c("Col_Mock_AVG", "Mut_Mock_AVG", "Col_XCC_AVG", "Mut_XCC_AVG")]
xcc_meanCenter_kmeans_sorted_kMeanClusters <- read.table(file = "Kmeans_6_meanCenter_Clusters_XCC.txt", header = T, sep = "\t", row.names = 1)
xcc_annotation_col <- data.frame(Cluster = factor(xcc_meanCenter_kmeans_sorted$Kmeans_6_cluster))
rownames(xcc_annotation_col) <- rownames(xcc_meanCenter_kmeans_sorted)
colnames(xcc_meanCenter_kmeans_sorted_values) <- c("Col_Mock_Avg", "CSLC_Mock_Avg", "Col_XCC_Avg", "CSLC_XCC_Avg")

kmeans_colors <- list(Kmeans_6_Cluster=c("1" = "blue",
                                         "2" = "red",
                                         "3" = "yellow",
                                         "4" = "orange",
                                         "5" = "green",
                                         "6" = "purple"))

pheatmap(xcc_meanCenter_kmeans_sorted_values, cluster_rows = F, cluster_cols = F,
         main = "XCC Experiment Mean Center KMean Clusters",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = xcc_annotation_col,
         annotation_colors = kmeans_colors,
         gaps_col=c(1,2),
         breaks=seq(-2,2, length.out=21), legend = T
)

library(ggplot2)
library(GGally)
library(dplyr)
library(patchwork)
mean_per_cluster <- read.table(file = "for_plotting_FPKM_average_mean_center_per_cluster.txt", header = T, sep = "\t")
cluster1_plot <- ggplot(mean_per_cluster, aes(x = Sample, y = Cluster_1, group = 1)) +
  geom_line( color = "blue", lwd = 1) +
  geom_point( color = "blue", shape = 21, fill = "blue", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Mean Centered FPKM Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 1")
cluster1_plot

cluster2_plot <- ggplot(mean_per_cluster, aes(x = Sample, y = Cluster_2, group = 1)) +
  geom_line( color = "red", lwd = 1) +
  geom_point( color = "red", shape = 21, fill = "red", size = 5) +
  ylim(-3, 3) +
  labs(x = "Sample Condition", y = "Mean Centered FPKM Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 2")
cluster2_plot

cluster3_plot <- ggplot(mean_per_cluster, aes(x = Sample, y = Cluster_3, group = 1)) +
  geom_line( color = "yellow", lwd = 1) +
  geom_point( color = "yellow", shape = 21, fill = "yellow", size = 5) +
  ylim(-1, 1) +
  labs(x = "Sample Condition", y = "Mean Centered FPKM Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 3")
cluster3_plot

cluster4_plot <- ggplot(mean_per_cluster, aes(x = Sample, y = Cluster_4, group = 1)) +
  geom_line( color = "orange", lwd = 1) +
  geom_point( color = "orange", shape = 21, fill = "orange", size = 5) +
  ylim(-2, 2) +
  labs(x = "Sample Condition", y = "Mean Centered FPKM Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 4")
cluster4_plot

cluster5_plot <- ggplot(mean_per_cluster, aes(x = Sample, y = Cluster_5, group = 1)) +
  geom_line( color = "green", lwd = 1) +
  geom_point( color = "green", shape = 21, fill = "green", size = 5) +
  ylim(-2, 2) +
  labs(x = "Sample Condition", y = "Mean Centered FPKM Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 5")
cluster5_plot

cluster6_plot <- ggplot(mean_per_cluster, aes(x = Sample, y = Cluster_6, group = 1)) +
  geom_line( color = "purple", lwd = 1) +
  geom_point( color = "purple", shape = 21, fill = "purple", size = 5) +
  ylim(-3, 3) +
  labs(x = "Sample Condition", y = "Mean Centered FPKM Average") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Cluster 6")
cluster6_plot

cluster1_plot + cluster2_plot + cluster3_plot + cluster4_plot + cluster5_plot + cluster6_plot + plot_annotation(title = "Kmeans Mean Center Cluster Graphs - XCC 2025-09-10")
ggsave("XCC_Kmeans_expressionClusters_20250910.pdf")
ggsave("XCC_Kmeans_expressionClusters_20250910.eps")

setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/XCC_Exp/gProfiler_results/")
goBPdata <- read.table(file = "XCC_topTenGO_BP_Cluster_negLogPval.txt", header = T, row.names = 1, sep = "\t")
colnames(goBPdata) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")
cluster_list <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")
annotation_col <- data.frame(Cluster = factor(colnames(goBPdata)))
rownames(annotation_col) <- colnames(goBPdata)
ann_colors_col <- list(Cluster=c("Cluster1" = "blue",
                                 "Cluster2" = "red",
                                 "Cluster3" = "yellow",
                                 "Cluster4" = "orange",
                                 "Cluster5" = "green",
                                 "Cluster6" = "purple"))
goBPplot <- pheatmap(goBPdata, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "complete",
         main = "        GO:BP Top 10 per Cluster XCC\n(neg log10 adj P-value)",
         show_rownames = T, show_colnames = T,
         border_color = "grey50", fontsize_col = 12,
         color = colorRampPalette(c("white", "lightblue1", "darkblue"))(60),
         annotation_col = annotation_col,
         annotation_colors = ann_colors_col,
         breaks=seq(0,60, length.out=61),
         labels_col=c("1","2","3","4","5","6"),
         cutree_rows = 6, cellwidth = 20, cellheight = 10
)

pdf("XCC_GO_BP_top10_per_Cluster_Heatmap_20290910.pdf", height = 10, width = 10)
goBPplot
dev.off()
 