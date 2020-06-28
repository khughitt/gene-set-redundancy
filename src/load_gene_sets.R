#!/bin/env Rscript
#
# Load gene sets into a single list and store as an rda file
#
suppressMessages(library(GSEABase))

# load gene sets
gmt_paths <- Sys.glob(snakemake@config$gmts)

gene_sets <- lapply(gmt_paths, function(x) {
  res <- geneIds(getGmt(x))
  lapply(res, function(gset) { gset[gset != ''] })
})
names(gene_sets) <- tools::file_path_sans_ext(basename(gmt_paths))

# unpack nested list
#gene_sets <- unlist(gene_sets, recursive = FALSE)

saveRDS(gene_sets, file = snakemake@output[[1]])
