# R Scripts Analysis Report

This report analyzes 37 R scripts for the papillomavirus (PV) discovery pipeline to help prepare them as supplementary material for publication.

---

## 1. Script Categorization

### PV-Related Core Scripts (Recommended for Publication)

| Script | Purpose | Priority |
|--------|---------|----------|
| `logan_2.R` | HMMER output analysis for PV database iteration (E1, E2, L1, L2 domain identification) | High |
| `logan_2_v2.R` | Extended version with L1 characterization via B/I domains | High |
| `logan_3.R` | Feb 7 PV analysis with novelty search, full L1 identification | High |
| `pv_2_stats.R` | Diamond output processing for PV hits, metadata joining | High |
| `pro_file_parsing.R` | Diamond .pro file parsing with `clean_inputs()` function | High |
| `ncbi_L1_character.R` | NCBI L1 sequence characterization, domain order verification | High |
| `errmm_twees.R` | Main tree visualization (1017 lines), gheatmap override, annotations | High |
| `gene_maps_for_pub.R` | Case study visualizations (pangolin, rhino, lizard, human, salmon) | High |
| `biomes.R` | Biome/ecology analysis with WWF ecoregion shapefiles | Medium |
| `geo_analysis.R` | Geographic resampling analysis | Medium |
| `poolygons.R` | Geographic analysis with sf/terra for PV sampling density maps | Medium |
| `e6_e7.R` | E6/E7 gene evolution analysis, all-vs-all comparison | Medium |
| `E6_andE7_synthesis.R` | E6/E7 novel sequence finding via HMMER | Medium |
| `reclass_viz.R` | Accessory gene distribution and late/early gene ratios | Medium |
| `sum_stats_longcontig.R` | Summary statistics for clustered ORFs, presence/absence heatmaps | Medium |
| `bugging.R` | Tree annotation with host metadata | Low |
| `bugging2.R` | L1 domain inversion analysis (detailed) | Low |
| `more_trees.R` | Tree annotation with SRA and NCBI metadata | Low |
| `genome_graphs.R` | Gene visualization using gggenes for various organisms | Low |
| `vert_counts.R` | Vertebrate PV counts vs sequencing data analysis | Low |
| `unsticker_better_vis.R` | Diamond .pro file visualization with gggenes | Low |

### Scripts That Require Manual Review

| Script | Issue | Action Needed |
|--------|-------|---------------|
| `pv_cluster_investigation.R` | References undefined variable `askansgag_2` | Check if this is generated elsewhere |
| `hmmer_results.R` | 511 lines, complex analysis | Review for completeness |
| `L1_pathracer_mining.R` | Pathracer pipeline setup | May need pathracer tool availability |
| `pathracer_fullness_hmm.R` | Pathracer output analysis | Depends on pathracer outputs |
| `twees.R` | Mixed PV/petase content | Remove petase sections if publishing |

### Non-PV Scripts (Consider Excluding)

| Script | Content | Recommendation |
|--------|---------|----------------|
| `dge_test.R` | DESeq2 differential gene expression | Exclude unless relevant to HPV E7 |
| `dge_test2.R` | DESeq2 analysis (similar to above) | Exclude |
| `misc_coding.R` | MMTV and anellovirus analysis | Exclude (not PV) |
| `recount_test.R` | GEO RNA-seq data retrieval functions | Helper functions only, exclude |
| `matrix_generation.R` | Uses Windows paths, incomplete | Exclude or fix paths |
| `idk_what_this_is.R` | Petase project references | Exclude |

---

## 2. Duplicate/Overlapping Scripts

### Highly Overlapping Scripts (Consolidate)

**Group A: HMMER Parsing and Logan Analysis**
- `logan_2.R`, `logan_2_v2.R`, `logan_3.R`
- All contain the `hmmsearch_clean()` function
- **Recommendation**: Keep `logan_3.R` as most complete, or create a single unified script

**Group B: Pathracer Analysis**
- `pathracer_fullness_hmm.R`, `ermmmm_fullness.R`, `idk_what_this_is.R`
- Very similar content
- **Recommendation**: Keep `pathracer_fullness_hmm.R` only

**Group C: Tree Visualization**
- `errmm_twees.R`, `twees.R`, `more_trees.R`, `bugging.R`
- All use ggtree for phylogenetic visualization
- **Recommendation**: Keep `errmm_twees.R` as primary, reference others for specific analyses

**Group D: E6/E7 Analysis**
- `e6_e7.R`, `E6_andE7_synthesis.R`, `pv_cluster_investigation.R`
- Related E6/E7 analyses
- **Recommendation**: Consolidate into single script if publishing together

---

## 3. Common Function: `hmmsearch_clean()`

This function appears in **8+ scripts** with minor variations. Consider extracting to a shared utilities file.

```r
# Standard version (from logan_2.R)
hmmsearch_clean <- function(input_domtbl){
  input_domtbl_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>%
    na.omit() %>% .[,1:22]
  input_domtbl_scan <- input_domtbl_scan[!(is.na(input_domtbl_scan$V21) |
    input_domtbl_scan$V21=="" | !(input_domtbl_scan$V2 == "-")), ]

  colnames(input_domtbl_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc",
    "tlen", "eval_full", "score_full", "bias_full", "#", "of", "c_eval", "i_eval",
    "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from",
    "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}
```

---

## 4. Bugs and Issues Found

### Critical Issues

| Script | Line | Issue |
|--------|------|-------|
| `pv_cluster_investigation.R` | 36, 45, 77, 88 | References undefined `askansgag_2$aaa` |
| `pv_cluster_investigation.R` | 214 | Empty `write.table()` call |
| `vert_counts.R` | 2 | Typo: `gbnstall.packages("ggpmisc")` should be `install.packages()` |
| `vert_counts.R` | 114 | Syntax error: `n_occur <- data.frame(table(vocabulary$id))` on same line as another statement |
| `reclass_viz.R` | 98 | References undefined `el_ratio` (should be `le_ratio`) |
| `genome_graphs.R` | 62-68 | Uses `manidae` data frame instead of `rat` for rat plot |
| `genome_graphs.R` | 140-148 | Uses `manidae` instead of proper data |

### Minor Issues

| Script | Line | Issue |
|--------|------|-------|
| `matrix_generation.R` | Multiple | Windows-style paths (`C:\\...`) won't work on Linux/Mac |
| `pro_file_parsing.R` | Multiple | Uses `browser()` debugging calls left in code |
| `more_trees.R` | 16, 40 | Uses `sep = "/t"` instead of `sep = "\t"` for tab separator |
| `poolygons.R` | 820 | Typo: `ge  om_sf` split across lines |

---

## 5. Required Input Files

### High-Priority Files (Core Pipeline)

| File | Used By | Format |
|------|---------|--------|
| `pv2_fasta_star_sto.domtbl` | logan_2.R, logan_2_v2.R | HMMER domtbl |
| `feb_7_pv_orf.domtbl` | logan_3.R | HMMER domtbl |
| `all_all_e6_fil.tsv` | e6_e7.R | TSV (Diamond output) |
| `E1_lengths.txt`, `E2_lengths.txt`, etc. | logan_2.R | CSV (protein lengths) |
| `concat_pv2.csv.gz` | pv_2_stats.R | Compressed CSV |
| `all_SRA_md.csv` | pv_2_stats.R | SRA metadata |
| `L1_B_super5.hmm`, `L1_I_super5.hmm` | ncbi_L1_character.R | HMM profiles |

### Tree Files

| File | Used By |
|------|---------|
| `*.treefile` | errmm_twees.R, twees.R, more_trees.R |
| `ncbi_and_novel_tree_final_*.nhx` | Various tree scripts |
| `FastTree_output_tree*.nhx` | more_trees.R |

### Metadata Files

| File | Used By |
|------|---------|
| `tree_meta.csv` | errmm_twees.R |
| `ncbi_info.txt` | more_trees.R |
| `sra_metadata_all_sra.txt` | more_trees.R |
| `geo_data_annotation_for_all_biosamps.txt` | poolygons.R |

### Geographic/Biome Files

| File | Used By |
|------|---------|
| `wwf_terr_ecos.shp` | biomes.R |
| `table_1_js.csv` | poolygons.R |
| `all_geo_grids_table_50km_by_50km_new_proj.tsv` | poolygons.R |

---

## 6. Recommended File Structure for Publication

```
R_scripts/
├── 00_utilities/
│   └── hmmsearch_clean.R              # Shared parsing function
│
├── 01_database_construction/
│   ├── logan_analysis.R               # Consolidated from logan_2/2v2/3.R
│   └── pv_diamond_processing.R        # From pv_2_stats.R
│
├── 02_sequence_characterization/
│   ├── L1_characterization.R          # From ncbi_L1_character.R
│   ├── e6_e7_analysis.R               # Consolidated from e6_e7.R + E6_andE7_synthesis.R
│   └── gene_presence_absence.R        # From sum_stats_longcontig.R + reclass_viz.R
│
├── 03_phylogenetics/
│   ├── tree_visualization.R           # Main from errmm_twees.R
│   └── tree_annotation.R              # From bugging.R, more_trees.R
│
├── 04_geographic_analysis/
│   ├── biome_analysis.R               # From biomes.R
│   ├── geographic_resampling.R        # From geo_analysis.R
│   └── spatial_mapping.R              # From poolygons.R
│
├── 05_case_studies/
│   └── gene_maps_for_publication.R    # From gene_maps_for_pub.R
│
└── data/
    ├── README.md                      # Document all required input files
    └── [input files listed above]
```

---

## 7. Action Items Before Publication

### Must Fix
1. [ ] Fix typo in `vert_counts.R:2` - `gbnstall.packages` -> `install.packages`
2. [ ] Fix undefined variable `askansgag_2` in `pv_cluster_investigation.R`
3. [ ] Fix `el_ratio` typo in `reclass_viz.R:98`
4. [ ] Fix tab separator typos (`/t` -> `\t`) in `more_trees.R`
5. [ ] Fix data frame name errors in `genome_graphs.R` (manidae used instead of rat)

### Should Do
1. [ ] Extract `hmmsearch_clean()` to shared utilities file
2. [ ] Remove petase-related code from `twees.R` if publishing
3. [ ] Convert Windows paths in `matrix_generation.R` to relative paths
4. [ ] Remove `browser()` calls from `pro_file_parsing.R`
5. [ ] Add file existence checks before `read.table()` calls
6. [ ] Document required input files with example formats

### Consolidation
1. [ ] Merge `logan_2.R`, `logan_2_v2.R`, `logan_3.R` into single script
2. [ ] Merge pathracer analysis scripts
3. [ ] Merge E6/E7 analysis scripts

---

## 8. Required R Packages

```r
# Core packages
library(tidyverse)      # dplyr, tidyr, ggplot2, stringr, etc.
library(reshape2)       # melt/dcast functions

# Phylogenetics
library(ggtree)         # Tree visualization
library(ape)            # Phylogenetic analysis
library(Polychrome)     # Color palette generation

# Visualization
library(viridis)        # Color palettes
library(ggpubr)         # Publication-ready plots
library(pheatmap)       # Heatmaps
library(gggenes)        # Gene arrow diagrams
library(ggbeeswarm)     # Beeswarm plots
library(ggupset)        # Upset plots
library(ggrepel)        # Label repelling

# Geographic/Spatial
library(sf)             # Simple features
library(terra)          # Raster data
library(rnaturalearth)  # World maps
library(rnaturalearthdata)

# Statistics
library(iNEXT)          # Rarefaction
library(matrixStats)    # Row statistics

# Bioinformatics (optional)
library(rentrez)        # NCBI Entrez
library(DESeq2)         # Differential expression (dge_test only)
```

---

## Summary

- **37 total scripts** analyzed
- **~25 PV-related** scripts recommended for publication
- **~6 non-PV scripts** can be excluded
- **3-4 script groups** have significant overlap and should be consolidated
- **5 critical bugs** identified that need fixing
- **40+ input files** referenced that need to be documented/provided

The scripts represent a comprehensive PV discovery pipeline with good analytical depth. The main improvements needed are consolidation of duplicate code, fixing a few bugs, and providing clear documentation of required input files.
