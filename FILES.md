# Large Data Files

All large data files required for this analysis are tracked via **Git LFS**.
No separate Zenodo download is needed — cloning the repository with `git lfs pull` will retrieve them.

## Git LFS-tracked files

### HMMER domain tables (`files/`)

| File | Size | Description |
|------|------|-------------|
| `feb_7_pv_fil_form.aa.domtbl` | 484 MB | HMMER search of all translated Logan contigs against PV HMM profiles |
| `all_pr_need_hmmer.pep.domtbl` | 394 MB | HMMER search of Pathracer-recovered peptides |
| `feb_7_pv_fil_form_L1_BI.domtbl` | 161 MB | HMMER search for L1 B/I domains in Logan contigs |
| `pv_all_L1.domtbl` | 138 MB | HMMER search of all PV L1 sequences |
| `mega_results.domtbl` | 88 MB | HMMER search of Pathracer mega-assembly results |
| `all_logan.domtbl` | 23 MB | HMMER search of all Logan contigs |
| `L1_BI_hmmer_iter2.domtbl` | 30 MB | Iteration 2 HMMER search for L1 B/I domains |

### Clustering and alignment files (`files/`)

| File | Size | Description |
|------|------|-------------|
| `pvdb2_v1_sort_clusters.uc` | 190 MB | USEARCH 90% identity clustering of all PV sequences |
| `all_full_l1s_logan_and_pr_blastn.tsv` | 47 MB | BLASTn results: full-length L1s against NCBI nt |
| `alignment_ncbi_centroids_and_pr.table` | 21 MB | Alignment table: NCBI centroids + Pathracer sequences |

### Metadata (`files/`)

| File | Size | Description |
|------|------|-------------|
| `L1_library_metadata.txt` | 14 MB | SRA library metadata for all L1-containing libraries |

### Outputs (`outputs/`)

| File | Size | Description |
|------|------|-------------|
| `2026.01.13.biomes_joined_sampling.gpkg` | 51 MB | GeoPackage with biome-grid intersection data |

## Reproduction instructions

1. Clone the repository: `git clone <repo-url>`
2. Pull LFS files: `git lfs pull`
3. Follow the run order described in `CLAUDE.md` and subdirectory READMEs
