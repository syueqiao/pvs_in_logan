# R Scripts - Manual Fixes Checklist

Use this checklist to track manual fixes needed before publication.

---

## Critical Bugs to Fix

- [ ] **04_geographic_analysis/vert_counts.R:2** - Fix typo: `gbnstall.packages("ggpmisc")` → `install.packages("ggpmisc")`

- [ ] **02_sequence_characterization/pv_cluster_investigation.R** - Fix undefined variable `askansgag_2` (appears on lines 36, 45, 77, 88). Either define this variable or update references to the correct data frame name.

- [ ] **02_sequence_characterization/reclass_viz.R:98** - Fix typo: `el_ratio` → `le_ratio`

- [ ] **03_phylogenetics/more_trees.R:16,40** - Fix tab separator: `sep = "/t"` → `sep = "\t"`

- [ ] **05_case_studies/genome_graphs.R:62-68** - Wrong data frame used: `manidae` should be `rat` for the rat genome plot

---

## Code Cleanup

- [ ] **01_database_construction/pro_file_parsing.R** - Remove `browser()` debugging calls

- [ ] **03_phylogenetics/twees.R** - Remove or comment out petase-related code sections (lines referencing petase project) if not relevant to PV publication

---

## Scripts to Consider Consolidating

These script groups have significant overlap. Consider merging if you want a cleaner repository:

### Logan Analysis Scripts
- `01_database_construction/logan_2.R`
- `01_database_construction/logan_2_v2.R`
- `01_database_construction/logan_3.R`

**Recommendation**: Keep `logan_3.R` as most complete, or merge into single `logan_analysis.R`

### Pathracer Scripts
- `02_sequence_characterization/pathracer_fullness_hmm.R`
- `02_sequence_characterization/ermmmm_fullness.R`

**Recommendation**: Keep `pathracer_fullness_hmm.R`, remove `ermmmm_fullness.R`

### E6/E7 Scripts
- `02_sequence_characterization/e6_e7.R`
- `02_sequence_characterization/E6_andE7_synthesis.R`
- `02_sequence_characterization/pv_cluster_investigation.R`

**Recommendation**: Review for complementary vs duplicate analyses

---

## Update Source Statements

After fixing bugs, add this line to scripts that use `hmmsearch_clean()`:

```r
source("00_utilities/hmmsearch_utils.R")
```

Scripts that need this update:
- [ ] `01_database_construction/logan_2.R`
- [ ] `01_database_construction/logan_2_v2.R`
- [ ] `01_database_construction/logan_3.R`
- [ ] `02_sequence_characterization/ncbi_L1_character.R`
- [ ] `02_sequence_characterization/pathracer_fullness_hmm.R`
- [ ] `02_sequence_characterization/ermmmm_fullness.R`
- [ ] `02_sequence_characterization/hmmer_results.R`
- [ ] `02_sequence_characterization/E6_andE7_synthesis.R`

---

## Input Files to Document/Provide

Before publication, ensure these key input files are available or documented:

### HMMER Output Files
- [ ] `pv2_fasta_star_sto.domtbl`
- [ ] `feb_7_pv_orf.domtbl`
- [ ] Various `.domtbl` files referenced in scripts

### Tree Files
- [ ] `*.treefile` (IQtree outputs)
- [ ] `*.nhx` files

### Metadata Files
- [ ] `all_SRA_md.csv`
- [ ] `tree_meta.csv`
- [ ] `ncbi_info.txt`

### Geographic Data
- [ ] `wwf_terr_ecos.shp` (WWF ecoregions shapefile)
- [ ] Geolocation metadata files

---

## Optional Improvements

- [ ] Add README.md files to each subdirectory explaining the scripts
- [ ] Add consistent file headers with author/date/description
- [ ] Standardize variable naming conventions across scripts
- [ ] Add input/output documentation to each script header
