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
- [ ] **05_case_studies/gene_maps_for_pub.R** - Multiple variables depend on `L1_based_trees_and_annotation.R` being run first: `final_tree`, `addtl_ano_tree_broad`, `addtl_ano_tree_broad_vals`, `biosamp_data_anno`, `all_genus_curated_thin`, `gheatmap()`, `tree_subset()`. Add a comment documenting this dependency.

### File Path Issues

- [ ] **05_case_studies/gene_maps_for_pub.R:6** - `SRR25256522_663139_pangolin_genome.txt` is missing `files/` prefix (inconsistent with other genome file reads)
- [ ] **05_case_studies/gene_maps_for_pub.R:347** - `SRR13789839_2669_human_genome.txt` is missing `files/` prefix
- [ ] **05_case_studies/gene_maps_for_pub.R:470** - `SRR20078264_4021_salmon_genome.txt` is missing `files/` prefix
- [ ] **05_case_studies/gene_maps_for_pub.R:629** - `generic_pv.txt` is missing `files/` prefix
- [ ] **05_case_studies/gene_maps_for_pub.R:274,443,565,598** - AlphaFold JSON files missing `files/` prefix (`fold_legless_zard_*.json`, `fold_human_*.json`, `fold_salmon_*.json`, `fold_pango_*.json`)
- [ ] **04_geographic_analysis/polygon_analyses.R:19** - `cluster_number_key_for_geo.tsv` missing `files/` prefix
- [ ] **04_geographic_analysis/polygon_analyses.R:26** - `who_puts_vlookup_man.list` missing `files/` prefix
- [ ] **04_geographic_analysis/polygon_analyses.R:30** - `geo_data_annotation_for_all_biosamps.txt` missing `files/` prefix
- [ ] **03_phylogenetics/L1_based_trees_and_annotation.R:923-924** - `wot` and `all_known_geo.tsv` missing `files/` prefix

---

## Non-Critical Bugs / Logic Issues

- [ ] **02_sequence_characterization/investigation_of_characteristics.R:9** - `fw_rev_contigs <- unique(all_ncbi_pv_sto_sto_fr_split$contig)` - accessing `$contig` on a list of data frames (result of `split()`). Should probably be accessing a specific element e.g. `$fw$contig` or using `lapply`.
- [ ] **04_geographic_analysis/polygon_analyses.R:380** - Uses `brewer.pal()` without loading `RColorBrewer` package
- [ ] **03_phylogenetics/L1_based_trees_and_annotation.R:988** - `geom_point(aes(x = lat_lon.X[1], y = lat_lon.X[2]...))` - likely should be `lat_lon.X` and `lat_lon.Y` columns, not indexing into X twice

---

## Missing R Package Dependencies

These packages are used but not declared with `library()` in the scripts that need them:

- [ ] **05_case_studies/gene_maps_for_pub.R** - Missing `library(gggibbous)` for `geom_moon()` (line 527)
- [ ] **04_geographic_analysis/polygon_analyses.R** - Missing `library(ggpmisc)` for `stat_poly_eq()` (line 827)
- [ ] **03_phylogenetics/L1_based_trees_and_annotation.R** - Missing `library(sf)` at top (used on line 956)
- [ ] **03_phylogenetics/L1_based_trees_and_annotation.R** - Missing `library(ggmosaic)` (used on line 979)

### Verify: `geom_aline()` Function

`geom_aline()` is used in `L1_based_trees_and_annotation.R` and `gene_maps_for_pub.R` but is not a standard ggtree/ggplot2 function. Verify it works in your environment. If it's from a custom source or dev version of ggtree, document where it comes from.

---

## Code Cleanup

- [ ] **05_case_studies/gene_maps_for_pub.R:236** - Localhost URL in comment (leftover from Jupyter): `http://127.0.0.1:40769/...`
- [ ] **05_case_studies/gene_maps_for_pub.R:635** - Same localhost URL in comment
- [ ] **03_phylogenetics/L1_based_trees_and_annotation.R:994** - y-axis label is `"wee"` - placeholder?
- [ ] **03_phylogenetics/L1_based_trees_and_annotation.R:1007-1066** - Loose code block at end of file (references `output_df`, `file_name` which don't exist). Appears to be leftover exploratory code - consider removing.
- [ ] **04_geographic_analysis/polygon_analyses.R** - Massive amount of commented-out code (lines ~420-828). Consider removing if no longer needed.

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
- [ ] Create a setup/initialization script that loads shared variables (`world`, `sra_metadata`, `final_tree`, `all_novels`, `all_genus_curated_thin`, etc.) used across multiple scripts
- [ ] Standardize all file paths to use `files/` prefix consistently
- [ ] Define `all_genus_curated_thin` - currently a mystery dependency used in multiple scripts but never explicitly created/read
