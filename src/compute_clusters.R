#!/bin/env Rscript
#
# Load gene sets into a single list and store as an rda file
#
suppressMessages(library(coop))
suppressMessages(library(dynamicTreeCut))
suppressMessages(library(NMF))
suppressMessages(library(viridis))
suppressMessages(library(readr))

set.seed(1)

# load matrix of gene set intersection sizes
intersections <- readRDS(snakemake@input[[1]])

# compute correlation matrix
cor_mat <- pcor(intersections)

# ignore negative correlations
cor_mat[cor_mat < 0] <- 0

# store correlation matrix
saveRDS(cor_mat, snakemake@output$correlations)

# perform hierarchical clustering
hc <- hclust(as.dist(1 - cor_mat), method = 'average')
clusters <- cutreeDynamicTree(hc)

# generate heatmaps of the original intersection data as well as the correlation matrix
ind <- sample(1:nrow(cor_mat), 1000)

png(snakemake@output$intersection_heatmap, width = 1080, height = 1080)
aheatmap(log2(intersections[ind, ind] + 1), color = viridis(100),
         annRow = factor(clusters[ind]),
         main = "Gene Set Intersection Size Heatmap (Log2)")
dev.off()

png(snakemake@output$correlation_heatmap, width = 1080, height = 1080)
aheatmap(cor_mat[ind, ind], color = magma(100),
         annRow = factor(clusters[ind]),
         main = "Gene Set Correlation Heatmap")
dev.off()

# save cluster assignments
res <- data.frame(gene_set = rownames(cor_mat), cluster = clusters)

write_tsv(res, path = snakemake@output[['clusters']])
