Case study scripts for detailed visualization of individual PV genomes.

- `gene_maps_for_pub.R` — Gene map visualizations for pangolin, rhinoceros, lizard, human, and salmon PV case studies. ORFs were manually formatted to gggenes-compliant style; jellyroll regions were manually annotated. Uses AlphaFold JSON data for structural confidence scoring.
- `blastn_hits_vis.R` — BLASTn hit visualization and gene map rendering for novel contigs. Includes strand-aware coordinate flipping and the `plot_genome_map()` function.
