"""
Gene Set Redundancy Analysis Pipeline

KH 2020.06.16
"""
from os.path import join
import pandas as pd

# base input and output data dirs;
# output comes from the separate feature association ("fassoc") pipeline
out_dir = join(config["output_dir"], config["name"], config["version"])

rule compute_clusters:
    input:
        join(out_dir, "overlap", "count_matrix.rda")
    output:
        correlations=join(out_dir, "overlap", "correlation_matrix.rda"),
        clusters=join(out_dir, "clusters", "clusters.tsv"),
        intersection_heatmap=join(out_dir, "plots", "intersection_heatmap.png"),
        correlation_heatmap=join(out_dir, "plots", "correlation_heatmap.png")
    script:
        "src/compute_clusters.R"

rule compute_overlap:
    input:
        join(out_dir, "gene_sets", "gene_sets.rda")
    output:
        counts=join(out_dir, "overlap", "count_matrix.rda")
        jaccard=join(out_dir, "overlap", "jaccard_matrix.rda"),
        metadata=join(out_dir, "metadata.tsv")
    script:
        "src/compute_overlap.R"

rule load_gmts:
    output:
        join(out_dir, "gene_sets", "gene_sets.rda")
    script:
        "src/load_gene_sets.R"


