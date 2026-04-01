# Input Files Catalog

This directory contains all input data files for the pvs_in_logan pipeline. Files are organized by type below.

> **Note**: Some R scripts reference files in this directory without the `files/` prefix. See `R_scripts/TO_DO.md` for tracked path issues.

## HMMER Domain Tables (.domtbl)

| File | Size | Description | Used by |
|------|------|-------------|---------|
| `feb_7_pv_fil_form.aa.domtbl` | ~484 MB | HMMER domain hits for all filtered PV amino acid ORFs | `logan_3.R` |
| `feb_7_pv_fil_form_L1_BI.domtbl` | ~161 MB | HMMER domain hits for L1 B/I region models | `logan_3.R` |
| `all_logan.domtbl` | - | HMMER domain hits for all Logan contigs (accessory genes) | `logan_3.R` |
| `all_pr_need_hmmer.pep.domtbl` | ~394 MB | HMMER domain hits for pathracer protein sequences | `logan_3.R` |
| `im_stupid.domtbl` | ~30 MB | HMMER hits for first clustering iteration | - |
| `im_stupid_2.domtbl` | ~30 MB | HMMER hits for second clustering iteration | `logan_3_novelty.R` |
| `pv_only_sto_sto.domtbl` | ~382 MB | HMMER hits for PV-only Stockholm alignment | - |

## USEARCH Cluster Files (.uc)

| File | Size | Description | Used by |
|------|------|-------------|---------|
| `pvdb2_v1_sort_clusters.uc` | ~190 MB | First-pass clustering of PV database v2 | - |
| `pvdb2_v1_sort_clusters_2.uc` | ~187 MB | Second-pass clustering of PV database v2 (NCBI + Logan ORFs) | `logan_3.R` |
| `all_full_l1s_logan_and_pr.nuc.uc` | - | Nucleotide clustering of full L1 sequences | - |
| `all_full_l1s_logan_and_pr_90.uc` | - | 90% identity clustering of full L1 sequences (main clustering for analysis) | `test_21k_host.R` |
| `all_pv_seqs.uc` | - | Clustering of all PV sequences | - |

## Phylogenetic Trees (.treefile, .nhx)

| File | Description | Used by |
|------|-------------|---------|
| `2026.03.24.tree_seqs.aln.treefile` | IQ-TREE output, main tree for figures (March 2026) | `L1_based_trees_and_annotation.R` |
| `2026.02.17.fasttree.nhx` | FastTree output with bootstrap values (Feb 2026) | - |
| `fasttree_ncbi_novel.nhx` | FastTree combining NCBI refs and novel sequences | - |
| `all_sequences_for_tree.aln.treefile` | Full alignment tree | - |
| `final_pr_ncbi_and_novel_sort_trim.aln.treefile` | Trimmed/sorted NCBI + novel tree | - |
| `ncbi_and_novel_for_tree_40.aln.treefile` | Small 40-sequence subset tree | - |
| `tree3.aln.treefile` | Third iteration tree | - |
| `w_bs.treefile` | Tree with bootstrap support | - |

## Alignment / BLAST Tables

| File | Description | Used by |
|------|-------------|---------|
| `alignment_envs.table` | Alignment of NCBI virus refs to Logan centroids (defines ncbi_virus clusters) | `logan_3_novelty.R`, `test_21k_host.R` |
| `alignment_ncbi_centroids_and_pr.table` | Alignment for NCBI centroids + novel PVs | `logan_3_novelty.R` |
| `alignment_ncbi_centroids.table` | Alignment filtered for NCBI centroid sequences | - |
| `alignment.table` | General query-subject alignment table | - |
| `blastn_hits_centroids.tsv` | BLASTN of Logan centroids vs NCBI nt (defines blastn clusters) | `logan_3_novelty.R`, `test_21k_host.R` |
| `all_full_l1s_logan_and_pr_blastn.tsv` | BLASTN of full L1 sequences vs reference database | `map_gen_for_pub.R` |
| `feb_7_L1_nuc_hmmer_env.tsv` | HMMER nucleotide envelope coordinates for L1 | `logan_3.R` |

## Metadata and Annotation Files

| File | Description | Used by |
|------|-------------|---------|
| `all_hits_info.list` | Comprehensive per-sequence metadata (library, coordinates, etc.) | `L1_based_trees_and_annotation.R`, `test_21k_host.R` |
| `L1_library_metadata.txt` | SRA/BioProject/taxonomy metadata per L1 library | `L1_pathracer_mining.R`, `test_21k_host.R` |
| `all_21323_libraries_metadata.list` | Metadata for all 21,323 SRA libraries screened | `test_21k_host.R` |
| `all_runs.csv` | Complete SRA run metadata | `L1_based_trees_and_annotation.R` |
| `lib_source_big.txt` | SRA library source type (GENOMIC, METAGENOMIC, etc.) | `logan_3.R` |
| `SRA_info.csv` | SRA run information including platform and taxonomy | - |
| `sra_metadata_all_sra.txt` | Complete SRA metadata dump | - |
| `2026.02.12all_genus_addtl_anno_fixed.csv` | Corrected genus-level annotations | `L1_based_trees_and_annotation.R` |
| `2025.10.29all_genus_addtl_anno.tsv.txt` | Additional genus annotations (Oct 2025) | - |
| `all_genus_curated.tsv.txt` | Hand-curated genus annotations | - |
| `2025.11.28.type_missing_manual.tsv.txt` | Manual annotations for missing types | `L1_based_trees_and_annotation.R` (commented) |
| `2025.10.10biosamp_novel_headers_geo_all.tsv.txt` | Biosample metadata with geographic info | - |
| `library_annotations_for_missing.csv` | Annotations for sequences missing metadata | - |

## Taxonomy and Host Files

| File | Description | Used by |
|------|-------------|---------|
| `host_table_ncbi.txt` | NCBI host taxonomy for papillomaviruses | `investigation_of_characteristics.R`, `ncbi_L1_character.R` |
| `tax_lineages.txt` | Full NCBI taxonomic lineages | - |
| `sequences_species.csv` | Sequence-to-host-species mapping with genus | - |
| `pap_genus_ncbi_virus.tsv` | PV genus info from NCBI Virus | - |
| `ncbi_tags_manual.txt` | Manual taxonomic tags | - |
| `ncbi_info.txt` | NCBI sequence info with taxon names | - |
| `ncbi_info_fill_redo.txt` | Corrected NCBI sequence info | - |
| `chordata_nomushomo.csv` | Chordata list excluding rodents/primates | - |
| `risks.txt` | HPV risk stratification data | - |
| `do_sparql_results.csv` | Disease ontology SPARQL query results | - |

## Geographic Data

| File | Description | Used by |
|------|-------------|---------|
| `geo_data_annotation_for_all_biosamps.txt` | Geographic + biosample annotations | `map_gen_for_pub.R`, `polygon_analyses.R` |
| `who_puts_vlookup_man.list` | Manual VLOOKUP for biosample-to-geo mapping | `map_gen_for_pub.R`, `polygon_analyses.R` |
| `grids_sf.gpkg` | GeoPackage with grid geometries for spatial analysis | `polygon_analyses.R` |
| `all_known_geo.tsv` | Known geographic locations for validated sequences | - |
| `cluster_number_key_for_geo.tsv` | Cluster number to geographic/taxonomic mapping | `polygon_analyses.R` |
| `hdi_2023.txt` | Human Development Index 2023 by country | `biomes.R` |

## Case Study Genome Files

| File | Description | Used by |
|------|-------------|---------|
| `SRR25256564_207500_pangolin_genome.txt` | Gene annotations for white-bellied pangolin PV | `gene_maps_for_pub.R` |
| `SRR10902309_46767_rhino_genome.txt` | Gene annotations for white rhinoceros PV | `gene_maps_for_pub.R` |
| `SRR22028468_199653_zard_genome.txt` | Gene annotations for legless lizard PV | `gene_maps_for_pub.R` |
| `SRR13789839_2669_human_genome.txt` | Gene annotations for human PV | `gene_maps_for_pub.R` |
| `SRR20078264_4021_salmon_genome.txt` | Gene annotations for salmon PV | `gene_maps_for_pub.R` |
| `fold_pango_srr25256564_207500_full_data_0.json` | AlphaFold structure prediction for pangolin PV | `gene_maps_for_pub.R` |
| `fold_rhino_srr10902309_46767_full_data_0.json` | AlphaFold structure prediction for rhino PV | `gene_maps_for_pub.R` |

## Sequence Lists and Accessions

| File | Description | Used by |
|------|-------------|---------|
| `L1_cluster_nt_acc_annot_NAs_key.txt` | Key mapping for L1 cluster nucleotide accessions | `investigation_of_characteristics.R`, `ncbi_L1_character.R` |
| `ncbi_JR_centroids_acc.nt` | NCBI centroid nucleotide accessions (FASTA headers) | `investigation_of_characteristics.R`, `ncbi_L1_character.R` |
| `new_seq.txt` | Newly identified NCBI sequences (475 accessions from second download) | `logan_3.R` |
| `run_in_blastn.list` | Sequences used in BLASTN analysis | `map_gen_for_pub.R` |
| `344_list.txt` | List of 344 accession numbers | - |
| `tree_seqs.list` | Sequences used in tree construction | - |
| `final_novel_pvs.headers` | FASTA headers for final novel PVs | - |
| `SRR_clean_acc` | Cleaned SRR accession list | - |
| `ncbi_accesions.list` | Complete NCBI PV accession list | - |
| `pv_accession.acc` | PV accession numbers | - |
| `pv_accessions_sort.txt` | Sorted PV accessions with FASTA headers | - |
| `complete_nucleotide_acc_herps.acc` | Herpesvirus nucleotide accessions (for comparison) | - |
| `acc_fail.txt` | Failed accessions that couldn't be processed | `L1_pathracer_mining.R` |

## Reference Databases

| File | Description | Used by |
|------|-------------|---------|
| `all_pv.txt` | All PV reference sequences | - |
| `animal_pave_genomes.txt` | Animal PV genomes from PAVE database | - |
| `animal_reference_clones.txt` | Reference clones with full taxonomy | - |
| `e5_pv.txt` through `e10_pv.txt` | HPV protein reference sequences (E5-E10) | - |

## Other Files

| File | Description | Used by |
|------|-------------|---------|
| `L1_clus.csv` | L1 clustering statistics | - |
| `anello_clus.csv` | Anellovirus clustering data (for comparison) | - |
| `log_mb_counts.csv` | Log-scale MB counts vs PV counts | `vert_counts.R` |
| `verts_info_clean_final_0306.txt` | Vertebrate sequence info (March 6) | `logan_3.R` |
| `uniq.counts` | Unique sequence counts across samples | `vert_counts.R` |
| `vert_counts_fr.txt` | Vertebrate counts data | `vert_counts.R` |
| `july1_3400_run.pro`, `july5_3400_run.pro` | ProML phylogenetic reconstruction outputs | - |
| `pvs.jpg` | Image file | - |
