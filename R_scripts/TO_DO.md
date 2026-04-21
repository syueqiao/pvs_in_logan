# R Scripts - Manual Fixes Checklist

Use this checklist to track manual fixes needed before publication.

---

## Script Reorganization (2026-02-17)

The following scripts were **deleted** and replaced by consolidated versions:

| Deleted Script | Replaced By |
|---|---|
| `02_sequence_characterization/bugging2.R` | `02_sequence_characterization/investigation_of_characteristics.R` |
| `02_sequence_characterization/hmmer_results.R` | (content folded into other scripts) |
| `02_sequence_characterization/reclass_viz.R` | (content folded into other scripts) |
| `02_sequence_characterization/sum_stats_longcontig.R` | (content folded into other scripts) |
| `03_phylogenetics/bugging.R` | `03_phylogenetics/L1_based_trees_and_annotation.R` |
| `03_phylogenetics/errmm_twees.R` | `03_phylogenetics/L1_based_trees_and_annotation.R` |
| `03_phylogenetics/more_trees.R` | `03_phylogenetics/L1_based_trees_and_annotation.R` |
| `03_phylogenetics/twees.R` | `03_phylogenetics/L1_based_trees_and_annotation.R` |
| `04_geographic_analysis/poolygons.R` | `04_geographic_analysis/polygon_analyses.R` |

---

## Resolved Bugs (Fixed)

- [x] **04_geographic_analysis/vert_counts.R:2** - Fix typo: `gbnstall.packages("ggpmisc")` should be `install.packages("ggpmisc")`
- [x] **04_geographic_analysis/vert_counts.R:114** - Syntax error: two statements on same line
- [x] **02_sequence_characterization/L1_pathracer_mining.R:59** - Fixed typo: `FALSEALSE` -> `FALSE`
- [x] **05_case_studies/gene_maps_for_pub.R:128-129** - Fixed syntax error: double comma in `tree_subset()` call
- [x] **03_phylogenetics/L1_based_trees_and_annotation.R:225** - Fixed typo: `"ouputs/"` -> `"outputs/"`
- [x] **01_database_construction/logan_3_novelty.R:56** - Fixed wrong variable: `all_seq_geo_fil_low_conf_hits` -> `blastn_results_fil_low_conf_hits_all`
- [x] **02_sequence_characterization/investigation_of_characteristics.R:1,175** - Removed `%%R` Jupyter magic commands

---

## Critical Bugs (Will Crash)

### Undefined Variables / Missing Dependencies

- [x] **03_phylogenetics/L1_based_trees_and_annotation.R:42** - `all_genus_curated_thin` is used but never defined or read in. Must be loaded from another script or a file.
- [x] **03_phylogenetics/L1_based_trees_and_annotation.R:893** - `all_pr_addt_orf_e_5_mat_pa_acc` is used but never defined. Must come from a prior script (likely sequence characterization).
- [x] **03_phylogenetics/L1_based_trees_and_annotation.R:1007-1012** - `output_df` and `file_name` are undefined. This `dates` plot code block appears to be leftover fragment code - remove or define the variables.
- [x] **05_case_studies/gene_maps_for_pub.R** - Multiple variables depend on `L1_based_trees_and_annotation.R` being run first: `final_tree`, `addtl_ano_tree_broad`, `addtl_ano_tree_broad_vals`, `biosamp_data_anno`, `all_genus_curated_thin`, `gheatmap()`, `tree_subset()`. Add a comment documenting this dependency.

### File Path Issues

- [x] **05_case_studies/gene_maps_for_pub.R** - `SRR13789839_2669_human_genome.txt` was missing `files/` prefix
- [x] **05_case_studies/gene_maps_for_pub.R** - `SRR20078264_4021_salmon_genome.txt` was missing `files/` prefix
- [x] **05_case_studies/gene_maps_for_pub.R** - `generic_pv.txt` was missing `files/` prefix
- [x] **05_case_studies/gene_maps_for_pub.R** - AlphaFold JSON files were missing `files/` prefix
- [x] **04_geographic_analysis/polygon_analyses.R** - `hits_library_biosample.list` (formerly `who_puts_vlookup_man.list`) was missing `files/` prefix
- [x] **04_geographic_analysis/polygon_analyses.R** - `geo_data_annotation_for_all_biosamps.txt` was missing `files/` prefix
- [x] **04_geographic_analysis/polygon_analyses.R** - `cluster_number_key_for_geo.tsv` — now generated in-script, read removed
- [x] **03_phylogenetics/L1_based_trees_and_annotation.R** - `wot` and `all_known_geo.tsv` — code block removed in prior cleanup

---

## Non-Critical Bugs / Logic Issues

- [x] **02_sequence_characterization/investigation_of_characteristics.R:9** - `fw_rev_contigs <- unique(all_ncbi_pv_sto_sto_fr_split$contig)` - accessing `$contig` on a list of data frames (result of `split()`). Removed — dead code, variable was never used downstream.
- [x] **04_geographic_analysis/polygon_analyses.R** - Added `library(RColorBrewer)` for `brewer.pal()`
- [x] **03_phylogenetics/L1_based_trees_and_annotation.R:988** - `geom_point(aes(x = lat_lon.X[1], y = lat_lon.X[2]...))` - code block was removed in prior cleanup (file is now 776 lines)

---

## Missing R Package Dependencies

These packages are used but not declared with `library()` in the scripts that need them:

- [x] **05_case_studies/gene_maps_for_pub.R** - Added `library(gggibbous)` for `geom_moon()`
- [x] **04_geographic_analysis/polygon_analyses.R** - `library(ggpmisc)` now present (line 203)
- [x] **03_phylogenetics/L1_based_trees_and_annotation.R** - `sf` and `ggmosaic` no longer used in current version

### Verify: `geom_aline()` Function

`geom_aline()` is used in `L1_based_trees_and_annotation.R` and `gene_maps_for_pub.R` but is not a standard ggtree/ggplot2 function. Verify it works in your environment. If it's from a custom source or dev version of ggtree, document where it comes from.

---

## Code Cleanup

- [x] **05_case_studies/gene_maps_for_pub.R** - Removed localhost URLs (leftover from Jupyter)
- [x] **03_phylogenetics/L1_based_trees_and_annotation.R** - `"wee"` label and loose code block removed in prior cleanup
- [x] **04_geographic_analysis/polygon_analyses.R** - Commented-out code was cleaned up in prior sessions (file is now 252 lines)

---

## Script Dependencies

### Phylogenetics + Case Studies Pipeline

Scripts have specific dependencies. Run in this order:

1. **L1_based_trees_and_annotation.R** (run first)
   - Reads: `files/2026.02.17.fasttree.nhx`, `files/2026.02.12all_genus_addtl_anno_fixed.csv`, `files/all_hits_info.list`, `files/tree_seqs.list`, `files/all_runs.csv`
   - Requires pre-loaded: `all_genus_curated_thin` (source TBD)
   - Creates: `final_tree`, `addtl_ano_tree_broad`, `addtl_ano_tree_broad_vals`, `sra_metadata`, `all_pvs_mapping` (from map_gen_for_pub.R), `biosamp_data_anno`, custom `gheatmap()` and `tree_subset()` functions
   - Also depends on: `map_gen_for_pub.R` (for `all_pvs_mapping` at line 232)

2. **gene_maps_for_pub.R** (depends on L1_based_trees_and_annotation.R)
   - Requires: `final_tree`, `addtl_ano_tree_broad`, `addtl_ano_tree_broad_vals`, `biosamp_data_anno`, `all_genus_curated_thin`, custom `gheatmap()` and `tree_subset()` functions

### Geographic Analysis Pipeline

1. **map_gen_for_pub.R** (run first)
   - Creates: `all_pvs_mapping`, `all_novels`, `known_pvs`, `geo_data_priority`, `world`

2. **polygon_analyses.R** (depends on map_gen_for_pub.R)
   - Requires: `all_pvs_mapping`, `all_novels` from map_gen_for_pub.R
   - Creates: `grid_sf`, `world`, `data_wide`, `giant_geo_table_grid_id_geometry`
   - Reads: `my_sf_data.gpkg`

3. **geo_analysis.R** (depends on polygon_analyses.R)
   - Requires: `giant_geo_table_grid_id_geometry`, `grid_sf`, `world`, `data_wide`
   - Creates: `sf_object_joined_prop_rarefy_5k_rowmeans_joined_s`

4. **biomes.R** (depends on polygon_analyses.R + geo_analysis.R)
   - Reads: `my_sf_data.gpkg`
   - Requires: `data_wide`, `grid_sf`, `world`, `all_novels`, `final_tree`

5. **vert_counts.R** (standalone - no dependencies)

### Sequence Characterization Pipeline

1. **investigation_of_characteristics.R** (depends on pre-loaded environment)
   - Requires: `all_ncbi_pv_sto_sto`, `NCBI_BI_contigs`, `all_ncbi_pv_sto_sto_BI_orfs`, `all_ncbi_pv_sto_sto_split`, `ncbi_tags_manual` - all from a prior session or logan_3.R
   - Sources: `00_utilities/hmmsearch_utils.R` (for `hmmsearch_clean()` at line 178)

---

## Update Source Statements

Add this line to scripts that use `hmmsearch_clean()`:

```r
source("00_utilities/hmmsearch_utils.R")
```

Scripts that need this update:
- [x] `01_database_construction/logan_3.R`
- [x] `02_sequence_characterization/ncbi_L1_character.R`
- [x] `02_sequence_characterization/pathracer_fullness_hmm.R`
- [x] `02_sequence_characterization/investigation_of_characteristics.R` (line 178 uses `hmmsearch_clean()`)

---

## Optional Improvements

- [x] Add README.md to each subdirectory explaining the scripts
- [ ] Add file headers with author/date/description
- [ ] Create a setup/initialization script that loads shared variables (`world`, `sra_metadata`, `final_tree`, `all_novels`, etc.) used across multiple scripts
- [x] Standardize all file paths to use `files/` or `outputs/` prefix consistently — fixed 8 bare ggsave paths in gene_maps_for_pub.R, 8 bare write.table paths in pathracer_fullness_hmm.R
- [x] `all_genus_curated_thin` — only referenced in commented-out lines in gene_maps_for_pub.R; no active code depends on it
