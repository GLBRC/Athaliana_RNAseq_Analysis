#install BiocManager
#install RColorBrewer
#install pheatmap
#install gplots
#install ggplot2
#install glmpca
#install plotly
#BiocManager::install("DESeq2")
#BiocManager::install("genefilter")


# 1. Col – Mock condition
# 2. Cslc – mock condition
# 3. Col – Xcc inoculation
# 4. Cslc – Xcc inoculation
# 5. Col – 3hr xyloglucan treatment
# 6. Col – 9 hr xyloglucan treatment

# We have two experiments.
# How Col and cslc plant respond to Xcc inoculation (sample 1,2,3,4)
# - pairwise condition (col-mock v col-xcc,  cslc-mock v cslc-xcc)
# - pairwise geno (col-mock v cslc-mock, col-xcc v cslc-xcc)
# - interaction  (geno-specific condition changes) - different in cslc-mock/xcc vs col-mock/xcc)

# How Xyloglucan trigger immune response (sample 1,5,6) – 1 will serve as 0hr control.
# - pairwise (0 v 3,  0 v 9, 3 v 9)


library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("pheatmap")
library("glmpca")
sessionInfo()

total_counts <- read.delim(file.path(getwd(),"counts.txt"), comment.char="#",header=TRUE,sep="\t")
total_counts$Chr = NULL
total_counts$Start = NULL
total_counts$End = NULL
total_counts$Strand = NULL
total_counts$Length = NULL
rownames(total_counts) <- total_counts$Geneid
total_counts$Geneid = NULL

# Full samples
colData <- data.frame(
  row.names=names(total_counts),
  geno = c(rep("col",9),rep("cslc",3),rep("col",3),rep("cslc",3)),
  cond = c(rep("xylo",6),rep("mock",6),rep("xcc",6)),
  time = c(rep("3hr",3),rep("9hr",3),rep("0hr",6),rep("24hr",6))
)
colData$group <- paste0(colData$cond,":",colData$geno,":",colData$time)
dds <- DESeqDataSetFromMatrix(
  countData = total_counts,
  colData = colData,
  design = ~group)

# pre-filtering - 10 reads in 3 or more samples
nrow(dds)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)

# distances - VST normalization
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Basic PCA
plotPCA(vsd, intgroup=c("cond","geno","time"))

# GLM-PCA - 2 dimensions
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$group <- paste0(dds$cond,"_",dds$geno,"_",dds$time)
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group)) +
  geom_point(size =3) + coord_fixed() + ggtitle("GLM-PCA - Generalized PCA")

# GLM-PCA - 3 dimensions
gpca <- glmpca(counts(dds), L=3)
gpca.dat <- gpca$factors
gpca.dat$group <- paste0(dds$cond,"_",dds$geno,"_",dds$time)
library(plotly)
library(htmlwidgets)
glm_plot = plot_ly(gpca.dat, x=~dim1, y=~dim2, z=~dim3, type="scatter3d", color=~group, mode="markers")
saveWidget(glm_plot, "glm_pca_3d.html", selfcontained = TRUE)


# Experiment one
# How Xyloglucan trigger immune response (sample 1,5,6) – 1 will serve as 0hr control.
# - pairwise (0 v 3,  0 v 9, 3 v 9)

total_counts_1 = total_counts[,(1:9)]
colData1 <- data.frame(
  row.names=names(total_counts_1),
  time = c(rep("3hr",3),rep("9hr",3),rep("0hr",3))
)
dds1 <- DESeqDataSetFromMatrix(
  countData = total_counts_1,
  colData = colData1,
  design = ~time)

#rld1 <- rlog(dds1)
#plotPCA(rld1, intgroup=c("time"))

dds1 <- DESeq(dds1)
res <- results(dds1)
res1_0_3 = results(dds1, contrast=c("time", "3hr", "0hr"))
res1_0_9 = results(dds1, contrast=c("time", "9hr", "0hr"))
res1_3_9 = results(dds1, contrast=c("time", "9hr", "3hr"))

summary(res1_0_3)
summary(res1_0_9)
summary(res1_3_9)

write.csv(as.data.frame(res1_0_3), file="0hr_v_3hr.csv")
write.csv(as.data.frame(res1_0_9), file="0hr_v_9hr.csv")
write.csv(as.data.frame(res1_3_9), file="3hr_v_9hr.csv")


# Experiment two
# How Col and cslc plant respond to Xcc inoculation (sample 1,2,3,4)
# - pairwise condition (col-mock v col-xcc,  cslc-mock v cslc-xcc)
# - pairwise geno (col-mock v cslc-mock, col-xcc v cslc-xcc)
# - interaction  (geno-specific condition changes) - different in cslc-mock/xcc vs col-mock/xcc)

total_counts_2 = total_counts[,(7:18)]

colData2 <- data.frame(
  row.names=names(total_counts_2),
  geno = c(rep("col",3),rep("cslc",3),rep("col",3),rep("cslc",3)),
  cond = c(rep("mock",6),rep("xcc",6))
)
colData2$group <- paste0(colData2$geno,":",colData2$cond)
dds2 <- DESeqDataSetFromMatrix(
  countData = total_counts_2,
  colData = colData2,
  design = ~group)

#rld2 <- rlog(dds2)
#plotPCA(rld2, intgroup=c("group"))


dds2 <- DESeq(dds2)

res2_col_mock_v_xcc = results(dds2, contrast=c("group", "col:xcc", "col:mock"))
res2_cslc_mock_v_xcc = results(dds2, contrast=c("group", "cslc:xcc", "cslc:mock"))
res2_mock_col_v_cslc = results(dds2, contrast=c("group", "cslc:mock", "col:mock"))
res2_xcc_col_v_cslc = results(dds2, contrast=c("group", "cslc:xcc", "col:xcc"))

write.csv(as.data.frame(res2_col_mock_v_xcc), file="col_xcc_v_mock.csv")
write.csv(as.data.frame(res2_cslc_mock_v_xcc), file="cslc_xcc_v_mock.csv")
write.csv(as.data.frame(res2_mock_col_v_cslc), file="mock_col_v_cslc.csv")
write.csv(as.data.frame(res2_xcc_col_v_cslc), file="xcc_col_v_cslc.csv")


# interaction term - extra condition effect of xcc in cslc, compared to col
dds2i <- DESeqDataSetFromMatrix(
  countData = total_counts_2,
  colData = colData2,
  design = ~geno + cond + geno:cond)
dds2i <- DESeq(dds2i)

res2_interaction = results(dds2i, name="genocslc.condxcc")

summary(res2_interaction)
write.csv(as.data.frame(res2_interaction), file="geno_cond_interaction.csv")

      
