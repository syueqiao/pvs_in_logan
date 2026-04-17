Scripts for phylogenetic tree construction, visualization, and host metadata analysis.

- `L1_based_trees_and_annotation.R` — Main tree visualization and annotation: reads tree files, builds host classification heatmaps, generates circular tree layouts, and defines custom `gheatmap()` and `tree_subset()` functions used downstream.
- `test_21k_host.R` — Library-level host metadata resolution (SRA metadata, BioSample API, manual curation), three-way novelty classification (ncbi_virus/blastn/novel), NCBI host concordance analysis, Fisher's exact tests, summary figures, and master table generation.

**Run order**: `L1_based_trees_and_annotation.R` must be run first; `test_21k_host.R` depends on objects it creates (e.g., `all_tree_metadata`, `final_tree`).
