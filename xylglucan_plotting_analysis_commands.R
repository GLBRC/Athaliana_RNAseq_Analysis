library(gplots)
library(viridis)

setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/Xylglucan_Exp/")

xylG_data <- read.table(file = "Xylglucan_FPKM_values.txt", header = T, sep = "\t", row.names = 1)

xylG_data_log2 <- log(xylG_data,2)

rowmean <- apply(xylG_data,1,mean)
rowmean_log2 <- apply(xylG_data_log2,1,mean)

xylG_meanCenter <- xylG_data - rowmean
head(xylG_meanCenter)

xylG_meanCenter_log2 <- xylG_data_log2 - rowmean_log2
head(xylG_meanCenter_log2)

xylG_meanCenter_matrix <- as.matrix(xylG_meanCenter)
xylG_meanCenter_log2_matrix <- as.matrix(xylG_meanCenter_log2)

#meanCenter time course cluster and heatmap

col_gradient3 = colorRampPalette(c("purple", "thistle1", "white", "white", "white", "orange", "darkred"))
breaks=c(-3,-2,-1,0,1,2,3)
gradient1 = colorpanel( sum(breaks[-1] <= -1 ), "purple", "thistle1")
gradient2 = colorpanel( sum(breaks[-1] > -1 & breaks[-1] < 1), "white", "white")
gradient3 = colorpanel( sum(breaks[-1] >= 1), "orange", "red")
hm.colors = c(gradient1, gradient2, gradient3)

breaks=c(-4,-2,-1,0,1,2,4)
pdf(file = "XylG_MeanCenter_Heatmap_Cluster_Log2.pdf")
heatmap.2(na.omit(xylG_meanCenter_log2_matrix), density.info = "none",
          trace = "none", symbreaks = F, breaks = breaks, dendrogram = "row",
          col = col_gradient3, Colv = FALSE, Rowv = TRUE, margins=c(12,9), scale = "none",
          cexRow=0.5, labRow = F, labCol = c("0H_1", "0H_2", "0H_3", "3H_1", "3H_2", "3H_3", "9H_1", "9H_2", "9H_3"),
          main = "Mean Center Log2 Xylglucan\nTime Course\nHeatmap and Clustering",
          key.xlab = "Log2 Mean Center", key.title = NA,
          )
dev.off()


library(pheatmap)
library(dendextend)
library(dichromat)
library(factoextra)
library(cluster)

setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/Xylglucan_Exp/")

xylG_FPKM_data <- read.table(file = "Xylglucan_FPKM_values_sig_1_comparison_ONLY.txt", header = T, sep = "\t", row.names = 1)

xylG_FPKM_data_noNA <- na.omit(xylG_FPKM_data)

avg_xylG_FPKM_data <- xylG_FPKM_data_noNA[ , c("Avg_0hr", "Avg_3hr", "Avg_9hr")]

avg_xylG_FPKM_data_log2 <- log(avg_xylG_FPKM_data,2)

rowmean <- apply(avg_xylG_FPKM_data,1,mean)
rowmean_log2 <- apply(avg_xylG_FPKM_data_log2,1,mean)

avg_xylG_FPKM_data_meanCenter <- avg_xylG_FPKM_data - rowmean
head(avg_xylG_FPKM_data_meanCenter)

avg_xylG_FPKM_data_log2_meanCenter <- avg_xylG_FPKM_data_log2 - rowmean_log2
head(avg_xylG_FPKM_data_log2_meanCenter)

avg_xylG_FPKM_data_log2_meanCenter_noNA <- na.omit(avg_xylG_FPKM_data_log2_meanCenter)

sig_data <- xylG_FPKM_data_noNA[ , c('FC_Sig_3v9', 'FC_Sig_0v9', 'FC_Sig_0v3')]
mean_center_colors <- list(FC_Sig_0v3=c("Yes" = "darkgreen",
                                     "no" = "grey"),
                           FC_Sig_0v9=c("Yes" = "darkgreen",
                                     "no" = "grey"),
                           FC_Sig_3v9=c("Yes" = "darkgreen",
                                     "no" = "grey"))

pheatmap(avg_xylG_FPKM_data_log2_meanCenter_noNA, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "complete",
         main = "Xylglucan Mean Center Comparisons\n(Sig in at least one comparison)",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = sig_data,
         annotation_colors = mean_center_colors,
         gaps_col=c(1,2), breaks=seq(-4,4, length.out=21)
)

data <- read.table("fc_fdr_sig_xylG.txt", header = T, row.names = "gene")
pheatmap(data)
data_noNA <- na.omit(data)
nrow(data_noNA)
nrow(data)

fc_data <- data_noNA[ , c('Log2FC_0v3', 'Log2FC_0v9', 'Log2FC_3v9')]
sig_0v3 <- data_noNA[ , c('Sig_0v3')]
sig_0v9 <- data_noNA[ , c('Sig_0v9')]
sig_3v9 <- data_noNA[ , c('Sig_3v9')]
sig_data <- data_noNA[ , c('Sig_0v3', 'Sig_0v9', 'Sig_3v9')]

ann_colors <- list(Sig_0v3=c("Sig" = "darkgreen",
                             "NS" = "grey"),
                   Sig_0v9=c("Sig" = "darkgreen",
                             "NS" = "grey"),
                   Sig_3v9=c("Sig" = "darkgreen",
                             "NS" = "grey"))

pheatmap(fc_data, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "complete",
         main = "Xylglucan Time Course Fold Change",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = sig_data,
         annotation_colors = ann_colors,
         gaps_col=c(1,2), cutree_rows = 6
         )

fc_data_vs0 <- data_noNA[ , c('Log2FC_0v3', 'Log2FC_0v9')]
sig_data_v0 <- data_noNA[ , c('Sig_0v9', 'Sig_0v3')]
ann_colors_v0 <- list(Sig_0v3=c("Sig" = "darkgreen",
                             "NS" = "grey"),
                   Sig_0v9=c("Sig" = "darkgreen",
                             "NS" = "grey")
                   )


#pheatmap(fc_data_vs0, kmeans_k = 3, culster_cols=F)
#kmeans_out <- pheatmap(fc_data_vs0, kmeans_k = 3, cluster_cols=F, display_numbers=T)
#clusters <- kmeans_out$kmeans$cluster
#write.table(clusters, file = "Kmeans_3_Clusters.txt")

#pheatmap(fc_data_vs0, kmeans_k = 6, culster_cols=F)
#kmeans6_out <- pheatmap(fc_data_vs0, kmeans_k = 6, cluster_cols=F, display_numbers=T)
#clusters6 <- kmeans6_out$kmeans$cluster
#write.table(clusters6, file = "Kmeans_6_Clusters.txt")

pdf(file = "XylG_FoldChange_Heatmap_Cluster_Log2_pheatmap.pdf")
pheatmap(fc_data_vs0, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "complete",
         main = "Xylglucan Time Course Fold Change",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = sig_data_v0,
         annotation_colors = ann_colors_v0,
         gaps_col=c(1,2),
         breaks=seq(-4,4, length.out=21)
)
dev.off()


fviz_nbclust(fc_data_vs0, kmeans, method = "wss")
fviz_nbclust(fc_data_vs0, kmeans, method = "silhouette")
fviz_nbclust(fc_data_vs0, kmeans, method = "gap_stat")
gap_stat <- clusGap(fc_data_vs0, FUN=kmeans, nstart=25, K.max=10, B=50)
fviz_gap_stat(gap_stat)

webshot::install_phantomjs(force=T)
library(heatmaply)

heatmaply(fc_data_vs0, Colv = FALSE, limits = c(-10,10), na.value = "gray50",
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
            low = "purple4", high = "darkgoldenrod3", midpoint = 0, limits = c(-10,10)),
          xlab = "Fold Change", ylab = "Gene", 
          main = "Xylglucan Time Course Fold Change", 
          file = "xylglucan_vs_0_cluster_heatmap.html", width = 800, height = 2500)

xylG_data <- read.table(file = "Xylglucan_FPKM_values.txt", header = T, sep = "\t", row.names = 1)
xylG_data_noNA <- na.omit(xylG_data)
xylG_data_noNA_means <- xylG_data_noNA[ , c('Avg_0hr', 'Avg_3hr', 'Avg_9hr')]
xylG_data_noNA_means_log2 <- log(xylG_data_noNA_means,2)
xylG_data_noNA_means_log2_zscore <- sapply(xylG_data_noNA_means_log2, function(xylG_data_noNA_means_log2) (xylG_data_noNA_means_log2 - mean(xylG_data_noNA_means_log2)) / sd(xylG_data_noNA_means_log2))

fviz_nbclust(xylG_data_noNA_means_log2, kmeans, method = "wss")
fviz_nbclust(xylG_data_noNA_means_log2, kmeans, method = "silhouette")
fviz_nbclust(xylG_data_noNA_means_log2, kmeans, method = "gap_stat")

fviz_nbclust(xylG_data_noNA_means_log2_zscore, kmeans, method = "wss")
fviz_nbclust(xylG_data_noNA_means_log2_zscore, kmeans, method = "silhouette")
fviz_nbclust(xylG_data_noNA_means_log2_zscore, kmeans, method = "gap_stat")

rowmean <- apply(xylG_data_noNA_means,1,mean)
rowmean_log2 <- apply(xylG_data_noNA_means_log2,1,mean)

xylG_meanCenter <- xylG_data_noNA_means - rowmean
head(xylG_meanCenter)
xylG_meanCenter_log2 <- xylG_data_noNA_means_log2 - rowmean_log2

write.table(xylG_meanCenter_log2, file="XylG_MeanCenter_TimeCourse_TimePoints.txt")

fviz_nbclust(xylG_meanCenter_log2, kmeans, method = "wss")
fviz_nbclust(xylG_meanCenter_log2, kmeans, method = "silhouette")


pheatmap(xylG_meanCenter_log2, kmeans_k = 5, culster_cols=F)
kmeans5_meanCenter_out <- pheatmap(xylG_meanCenter_log2, kmeans_k = 5, cluster_cols=F, display_numbers=T)
clusters5_meanCenter <- kmeans5_meanCenter_out$kmeans$cluster
write.table(clusters5_meanCenter, file = "Kmeans_5_meanCenter_Clusters.txt")

kmeans_5_meanCenter_gene_list <- read.table(file = "Kmeans_5_meanCenter_Clusters.txt", header = T, sep = "\t", row.names = 1)

sig_data_v0 <- kmeans_5_meanCenter_gene_list[ , c("Kmeans_5_Cluster", 'Sig_3v9', 'Sig_0v9', 'Sig_0v3')]
ann_colors_v0 <- list(Sig_0v3=c("Sig" = "darkgreen",
                                "NS" = "grey"),
                      Sig_0v9=c("Sig" = "darkgreen",
                                "NS" = "grey"),
                      Sig_3v9=c("Sig" = "darkgreen",
                                "NS" = "grey"),
                      Kmeans_5_Cluster=c("1" = "blue",
                                         "2" = "red",
                                         "3" = "yellow",
                                         "4" = "orange",
                                         "5" = "green")
                      )

pheatmap(xylG_meanCenter_log2, cluster_cols=F,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         breaks=seq(-4,4, length.out=21))

pdf(file = "XylG_MeanCenter_timeCourse_Heatmap_Cluster_Log2_pheatmap.pdf")
pheatmap(xylG_meanCenter_log2, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "complete",
         main = "Xylglucan Time Course Mean Center",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = sig_data_v0,
         annotation_colors = ann_colors_v0,
         gaps_col=c(1,2),
         breaks=seq(-4,4, length.out=21)
)
dev.off()

meanCenter_kMeans_sorted <- read.table(file = "XylG_MeanCenter_TimeCourse_TimePoints.txt", header = T, sep = "\t", row.names = 1)
meanCenter_kMeans_sorted_values <- meanCenter_kMeans_sorted[ , c("Avg_0hr", "Avg_3hr", "Avg_9hr")]
meanCenter_kMeans_sorted_kMeanClusters <- read.table(file = "XylG_MeanCenter_TimeCourse_TimePoints_ONLY_CLUSTERS.txt", header = T, sep = "\t", row.names = 1)

ann_colors_v0 <- list(Kmeans_5_Cluster=c("1" = "blue",
                                         "2" = "red",
                                         "3" = "yellow",
                                         "4" = "orange",
                                         "5" = "green"))

pheatmap(meanCenter_kMeans_sorted_values, cluster_rows = F, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "complete",
         main = "Xylglucan Time Course Mean Center KMean Clusters",
         show_rownames = F,
         border_color = F, fontsize_col = 12,
         color = colorRampPalette(c("purple", "white", "goldenrod"))(20),
         annotation_row = meanCenter_kMeans_sorted_kMeanClusters,
         annotation_colors = ann_colors_v0,
         gaps_col=c(1,2),
         breaks=seq(-2,2, length.out=21)
)

library(ggplot2)
library(GGally)
library(dplyr)
library(patchwork)
setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/Xylglucan_Exp/")
mean_per_cluster <- read.table(file = "for_plotting_average_mean_center_per_cluster.txt", header = T, sep = "\t")
plot(mean_per_cluster$Cluster_1)
cluster1 <- mean_per_cluster[ , c("Hour", "Cluster_1")]
cluster1_plot <- ggplot(mean_per_cluster, aes(x = Hour, y = Cluster_1)) +
  geom_line( color = "blue", lwd = 1) +
  geom_point( color = "blue", shape = 21, fill = "blue", size = 5) +
  xlim(0, 10) +
  ylim(-2, 2) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 1")
cluster1_plot

cluster2_plot <- ggplot(mean_per_cluster, aes(x = Hour, y = Cluster_2)) +
  geom_line( color = "red", lwd = 1) +
  geom_point( color = "red", shape = 21, fill = "red", size = 5) +
  xlim(0, 10) +
  ylim(-1, 1) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 2")
cluster2_plot

cluster3_plot <- ggplot(mean_per_cluster, aes(x = Hour, y = Cluster_3)) +
  geom_line( color = "yellow", lwd = 1) +
  geom_point( color = "yellow", shape = 21, fill = "yellow", size = 5) +
  xlim(0, 10) +
  ylim(-1, 1) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 3")
cluster3_plot

cluster4_plot <- ggplot(mean_per_cluster, aes(x = Hour, y = Cluster_4)) +
  geom_line( color = "orange", lwd = 1) +
  geom_point( color = "orange", shape = 21, fill = "orange", size = 5) +
  xlim(0, 10) +
  ylim(-1, 1) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 4")
cluster4_plot

cluster5_plot <- ggplot(mean_per_cluster, aes(x = Hour, y = Cluster_5)) +
  geom_line( color = "green", lwd = 1) +
  geom_point( color = "green", shape = 21, fill = "green", size = 5) +
  xlim(0, 10) +
  ylim(-2, 2) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 5")
cluster5_plot

cluster1_plot + cluster2_plot + cluster3_plot + cluster4_plot + cluster5_plot + plot_annotation(title = "Kmeans Mean Center Cluster Graphs - Xylglucan 2025-09-09")
ggsave("Xylglucan_Kmeans_expressionClusters_20250909.pdf")
ggsave("Xylglucan_Kmeans_expressionClusters_20250909.eps")

setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/Xylglucan_Exp/gProfiler_results/")
goBPdata <- read.table(file = "topTenGO_BP_Cluster_negLogPval.txt", header = T, sep = "\t", row.names = 1)
colnames(goBPdata) = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
cluster_list <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
annotation_col <- data.frame(Cluster = factor(colnames(goBPdata)))
rownames(annotation_col) <- colnames(goBPdata)
ann_colors_col <- list(Cluster=c("Cluster1" = "blue",
                                         "Cluster2" = "red",
                                         "Cluster3" = "yellow",
                                         "Cluster4" = "orange",
                                         "Cluster5" = "green"))

pheatmap(goBPdata, cluster_rows = T, cluster_cols = F,
         clustering_distance_rows = 'euclidean',
         clustering_method = "complete",
         main = "                   GO:BP Top 10 per Cluster (neg log10 adj P-value)",
         show_rownames = T, show_colnames = T,
         border_color = "grey90", fontsize_col = 12,
         color = colorRampPalette(c("white", "lightblue1", "darkblue"))(60),
         annotation_col = annotation_col,
         annotation_colors = ann_colors_col,
         breaks=seq(0,60, length.out=61),
         labels_col=c("1","2","3","4","5"),
         cutree_rows = 5
         )


# TPM VALUES ######################

setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/Xylglucan_Exp/")

xylG_tpm_data <- read.table(file = "xylglucan_tpm_counts_avg_sigValues.txt", header = T, sep = "\t", row.names = 1)
xylG_tpm_data_noNA <- na.omit(xylG_tpm_data)
xylG_tpm_data_noNA_means <- xylG_tpm_data_noNA[ , c('Avg_0hr', 'Avg_3hr', 'Avg_9hr')]
xylG_tpm_data_noNA_means_log2 <- log(xylG_tpm_data_noNA_means,2)

rowmean_tpm <- apply(xylG_tpm_data_noNA_means,1,mean)
rowmean_tpm_log2 <- apply(xylG_tpm_data_noNA_means_log2,1,mean)

xylG_tpm_meanCenter <- xylG_tpm_data_noNA_means - rowmean_tpm
head(xylG_tpm_meanCenter)
xylG_tpm_meanCenter_log2 <- xylG_tpm_data_noNA_means_log2 - rowmean_tpm_log2
head(xylG_tpm_meanCenter_log2)

write.table(xylG_tpm_meanCenter_log2, file="XylG_TPM_LOG2_MeanCenter_TimeCourse_TimePoints.txt")

pheatmap(xylG_tpm_meanCenter_log2, kmeans_k = 5, culster_cols=F)
kmeans5_tpm_meanCenter_out <- pheatmap(xylG_tpm_meanCenter_log2, kmeans_k = 5, cluster_cols=F, display_numbers=T)
clusters5_tpm_meanCenter <- kmeans5_tpm_meanCenter_out$kmeans$cluster
write.table(clusters5_tpm_meanCenter, file = "Kmeans_5_tpm_log2_meanCenter_Clusters.txt")

library(ggplot2)
library(GGally)
library(dplyr)
library(patchwork)
setwd("~/OneDrive - UW-Madison/Documents/Sang-Jin_Project/Xylglucan_Exp/")
mean_per_cluster_tpm <- read.table(file = "for_plotting_TPM_average_mean_center_per_cluster.txt", header = T, sep = "\t")
plot(mean_per_cluster$Cluster_1)
cluster1 <- mean_per_cluster[ , c("Hour", "Cluster_1")]
cluster1_plot <- ggplot(mean_per_cluster_tpm, aes(x = Hour, y = Cluster_1)) +
  geom_line( color = "blue", lwd = 1) +
  geom_point( color = "blue", shape = 21, fill = "blue", size = 5) +
  xlim(0, 10) +
  ylim(3, 5) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 1")
cluster1_plot

cluster2_plot <- ggplot(mean_per_cluster_tpm, aes(x = Hour, y = Cluster_2)) +
  geom_line( color = "red", lwd = 1) +
  geom_point( color = "red", shape = 21, fill = "red", size = 5) +
  xlim(0, 10) +
  ylim(2, 4) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 2")
cluster2_plot

cluster3_plot <- ggplot(mean_per_cluster_tpm, aes(x = Hour, y = Cluster_3)) +
  geom_line( color = "yellow", lwd = 1) +
  geom_point( color = "yellow", shape = 21, fill = "yellow", size = 5) +
  xlim(0, 10) +
  ylim(1, 6) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 3")
cluster3_plot

cluster4_plot <- ggplot(mean_per_cluster_tpm, aes(x = Hour, y = Cluster_4)) +
  geom_line( color = "orange", lwd = 1) +
  geom_point( color = "orange", shape = 21, fill = "orange", size = 5) +
  xlim(0, 10) +
  ylim(2, 4) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 4")
cluster4_plot

cluster5_plot <- ggplot(mean_per_cluster_tpm, aes(x = Hour, y = Cluster_5)) +
  geom_line( color = "green", lwd = 1) +
  geom_point( color = "green", shape = 21, fill = "green", size = 5) +
  xlim(0, 10) +
  ylim(1, 5) +
  labs(x = "Sample time (Hour)", y = "Mean Centered FPKM Average") +
  theme_classic() +
  ggtitle("Cluster 5")
cluster5_plot

cluster1_plot + cluster2_plot + cluster3_plot + cluster4_plot + cluster5_plot + plot_annotation(title = "Kmeans TPM Mean Center Cluster Graphs - Xylglucan 2025-09-10")
ggsave("Xylglucan_TPM_Kmeans_expressionClusters_20250910.pdf")
ggsave("Xylglucan_TPM_Kmeans_expressionClusters_20250910.eps")
