# R Scripts - Manual Fixes Checklist

Use this checklist to track manual fixes needed before publication.

---

## Critical Bugs to Fix

- [x] **04_geographic_analysis/vert_counts.R:2** - Fix typo: `gbnstall.packages("ggpmisc")` should be `install.packages("ggpmisc")`

- [x] **04_geographic_analysis/vert_counts.R:114** - Syntax error: two statements on same line. Split into separate lines.

- [x] **02_sequence_characterization/reclass_viz.R:98** - Fix typo: `el_ratio` should be `le_ratio`

- [x] **03_phylogenetics/more_trees.R:16,40** - Fix tab separator: `sep = "/t"` should be `sep = "\t"`

- [x] **05_case_studies/genome_graphs.R:62-68** - Wrong data frame used: `manidae` should be `rat` for the rat genome plot

---

## Code Cleanup

- [x] **01_database_construction/pro_file_parsing.R** - Remove `browser()` debugging calls

- [ ] **03_phylogenetics/twees.R** - Remove or comment out petase-related code sections if not relevant to PV publication

---

## Script Dependencies (Geographic Analysis)

The geographic analysis scripts have specific dependencies. Run in this order:

1. **poolygons.R** (run first)
   - Creates: `grid_sf`, `world`, `data_wide`, `giant_geo_table_grid_id_geometry`
   - Outputs: `my_sf_data.gpkg`
   - **Note**: Line 17 and 22 reference `all_pvs_mapping_binned` and `all_pvs_mapping` - these need to be loaded first (check if these come from another script or external source)

   #comes from map_gen_for_pub.R

2. **geo_analysis.R** (depends on poolygons.R)
   - Requires: `giant_geo_table_grid_id_geometry`, `grid_sf`, `world`, `data_wide`

3. **biomes.R** (depends on poolygons.R)
   - Reads: `my_sf_data.gpkg` (output from poolygons.R)
   - Requires: `data_wide`, `grid_sf`, `world`

4. **vert_counts.R** (standalone - no dependencies on other scripts)

---

## Update Source Statements

Add this line to scripts that use `hmmsearch_clean()`:

```r
source("00_utilities/hmmsearch_utils.R")
```

Scripts that need this update:
- [ ] `01_database_construction/logan_3.R`
- [ ] `02_sequence_characterization/ncbi_L1_character.R`
- [ ] `02_sequence_characterization/pathracer_fullness_hmm.R`
- [ ] `02_sequence_characterization/hmmer_results.R`

---

## Input Files to Document/Provide

### HMMER Output Files
- [x] `feb_7_pv_fil_form.aa.domtbl` (logan_3.R)
- [x] `feb_7_pv_fil_form_L1_BI.domtbl` (logan_3.R)
- [ ] Various `.domtbl` files referenced in scripts

### Tree Files
- [ ] `*.treefile` (IQtree outputs)
- [ ] `*.nhx` files

### Metadata Files
- [ ] `all_SRA_md.csv`
- [ ] `tree_meta.csv`
- [ ] `ncbi_info.txt`
- [x] `lib_source_big.txt` (logan_3.R)

### Geographic Data
- [x] `Ecoregions2017.shp` (biomes.R)
- [x] `geo_data_annotation_for_all_biosamps.txt` (poolygons.R)
- [x] `cluster_number_key_for_geo.tsv` (poolygons.R)
- [x] `who_puts_vlookup_man.list` (poolygons.R)

---

## Optional Improvements

- [ ] Add README.md to each subdirectory explaining the scripts
- [ ] Add file headers with author/date/description
