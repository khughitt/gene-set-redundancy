#!/bin/env Rscript
#
# Load gene sets into a single list and store as an rda file
#
suppressMessages(readr)

# load gene sets
gene_sets <- readRDS(snakemake@input[[1]])

# build gene set metadata file with <collection, gene_set, size> info for each gene set
mdat <- NULL

for (collection in names(gene_sets)) {
  gene_set_names <- names(gene_sets[[collection]])
  sizes <- unlist(lapply(gene_sets[[collection]], length))

  mdat <- rbind(mdat, data.frame(collection, gene_set = gene_set_names, size = sizes))
}

# remove collection name prior to combining gene sets into a single list
names(gene_sets) <- NULL

# unpack nested list
gene_sets <- unlist(gene_sets, recursive = FALSE)

if (any(duplicated(names(gene_sets)))) {
  stop("Duplicate gene set names encountered!")
}

unions <- matrix(NA, nrow = length(gene_sets), ncol = length(gene_sets))

for (i in 1:length(gene_sets)) {
  for (j in 1:length(gene_sets)) {
    message(sprintf("Iteration: %d, %d...", i, j))

    if (is.na(unions[i, j])) {
      union_len <- length(unique(c(gene_sets[[i]], gene_sets[[j]])))
      unions[i, j] <- union_len
      unions[j, i] <- union_len
    }
  }
}

rownames(unions) <- names(gene_sets)
colnames(unions) <- names(gene_sets)

# measure # shared gene_sets in top gene sets
intersections <- crossprod(table(stack(gene_sets)))

# convert to frequencies
jaccard <- intersections / unions

# save results
write_tsv(mdat, path = snakemake@output[['metadata']])

saveRDS(intersections, file = snakemake@output[['counts']])
saveRDS(jaccard, file = snakemake@output[['jaccard']])
