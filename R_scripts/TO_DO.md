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

- [x] **03_phylogenetics/twees.R** - Remove or comment out petase-related code sections if not relevant to PV publication

---

## Script Dependencies (Geographic Analysis)

The geographic analysis scripts have specific dependencies. Run in this order:

1. **poolygons.R** (run first)
   - Creates: `grid_sf`, `world`, `data_wide`, `giant_geo_table_grid_id_geometry`
   - Outputs: `my_sf_data.gpkg`
   - **Note**: Line 17 and 22 reference `all_pvs_mapping_binned` and `all_pvs_mapping` - these need to be loaded first (check if these come from another script or external source)

   #comes from map_gen_for_pub.R!

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
- [x] `01_database_construction/logan_3.R`
- [x] `02_sequence_characterization/ncbi_L1_character.R`
- [x] `02_sequence_characterization/pathracer_fullness_hmm.R`
- [x] `02_sequence_characterization/hmmer_results.R`

---

## Input Files to Document/Provide

### HMMER Output Files
- [x] `feb_7_pv_fil_form.aa.domtbl` (logan_3.R)
- [x] `feb_7_pv_fil_form_L1_BI.domtbl` (logan_3.R)
- [x] `PVfam2_scan_tbl.out` (hmmer_results.R)
- [x] `feb_7_L1_nuc_hmmer_env.tsv` (logan_3.R)
- [x] `feb7_all_pr_pilot_outputs_sort_centroids_huh.tsv.csv` (pilot_pr_test.R)
- [x] `vert3_sum.txt` (pathracer_fullness_hmm.R)

### Tree Files (IQtree and FastTree outputs)
- [x] `w_bs.treefile` (twees.R)
- [x] `FastTree_output_tree (3).nhx` (twees.R, more_trees.R)
- [x] `fasttree_ncbi_novel.nhx` (twees.R, more_trees.R)
- [x] `ncbi_and_novel_for_tree_40.aln.treefile` (twees.R)
- [x] `tree3.aln.treefile` (errmm_twees.R)
- [x] `final_pr_ncbi_and_novel_sort_trim.aln.treefile` (errmm_twees.R)

### SRA/Library Metadata Files
**Note**: `all_SRA_md.csv` and `tree_meta.csv` mentioned in older docs have been replaced by:
- [x] `3400_sra_metadata.txt` (pv_2_stats.R) - Main SRA metadata file
- [x] `lib_source_big.txt` (logan_3.R, errmm_twees.R) - Library source information
- [x] `sra_metadata_all_sra.txt` (more_trees.R)
- [x] `SRA_info.csv` (twees.R)
- [x] `ncbi_and_novel_sra_info.txt` (twees.R, bugging.R)
- [x] `all_runs.csv` (errmm_twees.R) - BioSample data

### NCBI/Taxonomy Metadata
- [x] `ncbi_info.txt` (multiple scripts)
- [x] `ncbi_info_fill_redo.txt` (twees.R, more_trees.R)
- [x] `ncbi_tags_manual.txt` (twees.R, more_trees.R)
- [x] `animal_reference_clones.txt` (twees.R)
- [x] `pap_genus_ncbi_virus.tsv` (errmm_twees.R) - BLAST taxonomy assignments
- [x] `tax_lineages.txt` (errmm_twees.R) - Taxonomy lineage data
- [x] `all_genus_curated.tsv.txt` (errmm_twees.R)
- [x] `2025.10.29all_genus_addtl_anno.tsv.txt` (errmm_twees.R)
- [x] `2025.11.28.type_missing_manual.tsv.txt` (errmm_twees.R)
- [x] `final_novel_pvs.headers` (errmm_twees.R)
- [x] `library_annotations_for_missing.csv` (errmm_twees.R)
- [x] `host_table_ncbi.txt` (bugging2.R, ncbi_L1_character.R)

### Geographic Data
- [x] `Ecoregions2017.shp` (biomes.R)
- [x] `geo_data_annotation_for_all_biosamps.txt` (poolygons.R, map_gen_for_pub.R)
- [x] `cluster_number_key_for_geo.tsv` (poolygons.R)
- [x] `who_puts_vlookup_man.list` (poolygons.R, map_gen_for_pub.R)
- [x] `all_full_l1s_logan_and_pr_blastn.tsv` (map_gen_for_pub.R)
- [x] `run_in_blastn.list` (map_gen_for_pub.R)
- [x] `table_1_js.csv` (poolygons.R)
- [x] `wot` (errmm_twees.R) - Geographic data for novel sequences
- [x] `all_known_geo.tsv` (errmm_twees.R)
- [x] `2025.10.10biosamp_novel_headers_geo_all.tsv.txt` (errmm_twees.R)
- [x] `hdi_2023.txt` (biomes.R) - Human Development Index data

### Sequence Accessions & Clustering
- [x] `pv_accession.acc` (hmmer_results.R)
- [x] `herp_acc.acc` (hmmer_results.R)
- [x] `poly_acc.acc` (hmmer_results.R)
- [x] `SRR_clean_acc` (hmmer_results.R)
- [x] `complete_nucleotide_acc_herps.acc` (hmmer_results.R)
- [x] `animal_pave_genomes.txt` (hmmer_results.R)
- [x] `pv_accessions_sort.txt` (hmmer_results.R)
- [x] `acc_pfam_pacc_sort.txt` (hmmer_results.R)
- [x] `ncbi_JR_centroids_acc.nt` (bugging2.R, ncbi_L1_character.R)
- [x] `L1_cluster_nt_acc_annot_NAs_key.txt` (bugging2.R, ncbi_L1_character.R)
- [x] `blacklist.txt` (bugging2.R, ncbi_L1_character.R) - QC-failed sequences

### Logan/Diamond Run Outputs
- [x] `july1_3400_run.pro` (pv_2_stats.R) - First Logan run
- [x] `july5_3400_run.pro` (pv_2_stats.R) - Second Logan run
- [x] `pv_only_evail_fil.pro` (pv_2_stats.R)
- [x] `pv_only_sto_sto_lengths.txt` (pv_2_stats.R)

### Classification/Characterization Data
- [x] `files/all_pv.txt` through `files/e10_pv.txt` (reclass_viz.R) - E4-E10 classification
- [x] `files/risks.txt` (reclass_viz.R)
- [x] `sequences_species.csv` (hmmer_results.R)
- [x] `anello_clus.csv` (hmmer_results.R)
- [x] `L1_clus.csv` (hmmer_results.R)

### Case Study Genome Files
- [x] `SRR25256522_663139_pangolin_genome.txt` (gene_maps_for_pub.R)
- [x] `SRR10902309_46767_rhino_genome.txt` (gene_maps_for_pub.R)
- [x] `files/SRR22028468_199653_zard_genome.txt` (gene_maps_for_pub.R)
- [x] `SRR13789839_2669_human_genome.txt` (gene_maps_for_pub.R)
- [x] `SRR20078264_4021_salmon_genome.txt` (gene_maps_for_pub.R)
- [x] `generic_pv.txt` (gene_maps_for_pub.R)

### Vertebrate/Host Analysis
- [x] `files/log_mb_counts.csv` (vert_counts.R)
- [x] `files/uniq.counts` (vert_counts.R)
- [x] `files/vert_counts_fr.txt` (vert_counts.R)
- [x] `verts_info_clean_final_0306.txt` (logan_3.R)
- [x] `chordata_nomushomo.csv` (hmmer_results.R)

### Other Data Files
- [x] `do_sparql_results.csv` (logan_search_eorfs.R) - Disease Ontology query results
- [x] `L1_library_metadata.txt` (L1_pathracer_mining.R)
- [x] `acc_fail.txt` (L1_pathracer_mining.R)

---

## Optional Improvements

- [ ] Add README.md to each subdirectory explaining the scripts
- [ ] Add file headers with author/date/description
