# CLAUDE.md - AI Assistant Guide for pvs_in_logan

## Project Overview

This repository documents a **papillomavirus (PV) discovery and characterization pipeline** that uses the Logan unsticker tool to screen petabases of sequencing data. The project identifies novel and known papillomaviruses from large metagenomic datasets (3400+ SRA libraries) and performs phylogenetic analysis.

**Primary Goal**: Generate and iteratively improve papillomavirus databases by:
1. Creating initial PV databases from NCBI and PAVE sequences
2. Running Logan searches against large metagenomic datasets
3. Identifying novel papillomaviruses with < 90% nucleotide identity to known sequences
4. Building phylogenetic trees to classify and visualize PV diversity

## Repository Structure

```
pvs_in_logan/
├── CLAUDE.md                          # This file - AI assistant guide
├── README.md                          # Project overview
├── R_SCRIPTS_ANALYSIS.md              # Detailed R scripts analysis report
├── R_scripts/                         # All R analysis scripts
│   ├── TO_DO.md                       # Bug tracker and checklist
│   ├── 00_utilities/
│   │   └── hmmsearch_utils.R          # Shared HMMER parsing functions
│   ├── 01_database_construction/
│   │   ├── logan_3.R                  # Feb 7 PV analysis, novelty search, full L1 ID
│   │   └── logan_3_novelty.R          # Updated novelty search with blastn centroids
│   ├── 02_sequence_characterization/
│   │   ├── investigation_of_characteristics.R  # L1 domain inversion/coverage analysis
│   │   ├── L1_pathracer_mining.R      # Pathracer pipeline setup
│   │   ├── ncbi_L1_character.R        # NCBI L1 characterization
│   │   ├── pathracer_fullness_hmm.R   # Pathracer output analysis
│   │   ├── pilot_pr_test.R            # Pathracer pilot testing
│   │   └── vert_counts.R             # Vertebrate PV counts analysis
│   ├── 03_phylogenetics/
│   │   └── L1_based_trees_and_annotation.R  # Main tree viz, annotation, host classification
│   ├── 04_geographic_analysis/
│   │   ├── biomes.R                   # Biome/ecology analysis (WWF ecoregions)
│   │   ├── geo_analysis.R             # Geographic resampling analysis
│   │   ├── map_gen_for_pub.R          # Publication-ready maps and novelty search
│   │   └── polygon_analyses.R         # Spatial mapping (sf/terra), grid analysis
│   └── 05_case_studies/
│       └── gene_maps_for_pub.R        # Case study gene maps (pangolin, rhino, etc.)
├── files/                             # Input data files
│   ├── *.txt, *.tsv, *.csv           # Metadata, annotations, accession lists
│   ├── *.treefile, *.nhx             # Phylogenetic tree files
│   ├── *.domtbl                      # HMMER domain table outputs
│   ├── *.acc                          # Accession lists
│   └── *.jpg, *.json                 # Images and AlphaFold data
└── outputs/                           # Script-generated output files
    ├── *.list, *.txt                  # Filtered sequence/library lists
    ├── *.bed                          # Coordinate files
    └── *.png                          # Generated visualizations
```

## Technology Stack

### Languages
- **R** (primary): Data manipulation and visualization via tidyverse, ggtree
- **Bash** (secondary): Bioinformatic tool orchestration and file processing

### External Bioinformatic Tools (Required)
```bash
# Sequence processing
seqkit          # Sequence manipulation
usearch         # Clustering and alignment (90% identity threshold)
getorf          # ORF finding (EMBOSS suite)
transeq         # Sequence translation (EMBOSS suite)

# HMM analysis
hmmsearch       # Profile HMM searching (HMMER v3.4)
hmmbuild        # Build HMM profiles
esl-reformat    # Reformat alignment files (Easel)

# Alignment and phylogenetics
muscle5.1       # Multiple sequence alignment
iqtree2         # Phylogenetic tree building
FastTree        # Fast approximate phylogenetic trees

# Sequence comparison
blastn          # NCBI nucleotide BLAST
diamond         # Fast protein alignment (via Logan)

# Utilities
bedtools        # Genome arithmetic
zstd            # Compression
aws s3          # Data transfer from cloud storage
```

### R Packages
```r
# Core
library(tidyverse)           # Data manipulation (dplyr, tidyr, ggplot2, stringr)
library(reshape2)            # melt/dcast functions

# Phylogenetics
library(ggtree)              # Tree visualization
library(treeio)              # Tree I/O
library(ape)                 # Phylogenetic analysis
library(Polychrome)          # Color palette generation
library(ggnewscale)          # Multiple scales in ggplot

# Visualization
library(viridis)             # Color palettes
library(ggpubr)              # Publication-ready plots
library(pheatmap)            # Heatmaps
library(gggenes)             # Gene arrow diagrams
library(ggbeeswarm)          # Beeswarm plots
library(ggupset)             # Upset plots
library(ggrepel)             # Label repelling
library(ggpmisc)             # Polynomial equation annotations
library(gggibbous)           # Moon/pie geoms (geom_moon)
library(ggmosaic)            # Mosaic plots
library(extrafont)           # Font management

# Geographic/Spatial
library(sf)                  # Simple features
library(terra)               # Raster data
library(rnaturalearth)       # World maps
library(rnaturalearthdata)   # Map data
library(tmap)                # Thematic maps

# Statistics
library(iNEXT)               # Rarefaction
library(matrixStats)         # Row statistics

# Bioinformatics
library(rentrez)             # NCBI Entrez API
library(jsonlite)            # JSON parsing (AlphaFold data)
```

## Pipeline Phases (R Scripts)

### Phase 1: Shared Utilities (`00_utilities/`)
- `hmmsearch_utils.R` - Defines `hmmsearch_clean()`, `hmmsearch_clean_info()`, and `clean_inputs()` for parsing HMMER and Diamond outputs
- **Must be sourced** by scripts in `01_*`, `02_*` that use these functions

### Phase 2: Database Construction (`01_database_construction/`)
- `logan_3.R` - Core pipeline: HMMER filtering, L1 B/I domain identification, novelty search, contig characterization
- `logan_3_novelty.R` - Updated novelty search using blastn centroids
- Sources `00_utilities/hmmsearch_utils.R`

### Phase 3: Sequence Characterization (`02_sequence_characterization/`)
- `investigation_of_characteristics.R` - L1 domain inversion analysis, domain coverage, missing B/I domain detection
- `ncbi_L1_character.R` - Characterize NCBI L1 sequences, verify domain order
- `vert_counts.R` - Vertebrate PV counts vs sequencing depth analysis
- `L1_pathracer_mining.R`, `pathracer_fullness_hmm.R`, `pilot_pr_test.R` - Pathracer-related analysis

### Phase 4: Phylogenetics (`03_phylogenetics/`)
- `L1_based_trees_and_annotation.R` - **Main tree visualization and annotation**: gheatmap override, host classification, geographic mapping, novelty status, library source analysis, ENA/NCBI metadata lookup, circular tree layout

### Phase 5: Geographic Analysis (`04_geographic_analysis/`)
Scripts have **run-order dependencies** (see TO_DO.md for details):
1. `map_gen_for_pub.R` -> creates `all_pvs_mapping`, `all_novels`, `known_pvs`
2. `polygon_analyses.R` -> creates `grid_sf`, `world`, `data_wide`, reads `my_sf_data.gpkg`
3. `geo_analysis.R` -> requires objects from polygon_analyses.R
4. `biomes.R` -> requires objects from polygon_analyses.R + geo_analysis.R
5. `vert_counts.R` - standalone

### Phase 6: Case Studies (`05_case_studies/`)
- `gene_maps_for_pub.R` - Gene map visualizations for pangolin, rhinoceros, lizard, human, and salmon PV case studies. Uses AlphaFold JSON data for structural confidence.

## Key Analysis Parameters

| Parameter | Value | Context |
|-----------|-------|---------|
| HMMER F1 | 0.01 | Filtering threshold |
| HMMER T | 20 | Score threshold |
| Clustering | 90% | AA and nucleotide identity |
| Diamond e-value | 10^-10 | High-confidence Logan hits |
| BLAST e-value | 10^-4 | Novelty search threshold |
| Novelty cutoff | < 90% | Nucleotide identity to NCBI nt |
| Bootstrap | 1000 | IQtree replicates |

## Conventions for AI Assistants

### R Code Style
```r
# Use tidyverse piping
data %>%
  filter(evalue < 1e-10) %>%
  group_by(sequence_id) %>%
  summarize(count = n())

# Custom functions follow snake_case
hmmsearch_clean <- function(file_path) { ... }

# ggtree annotation with %<+% operator
tree %<+% annotation_df
```

### File Path Conventions
- Input data files: `files/` directory
- Script outputs: `outputs/` directory
- Scripts reference files with relative paths from the project root
- **Note**: Some scripts have inconsistent path prefixes - see TO_DO.md

### Script Execution Notes
- Many scripts are **interactive/exploratory** - designed to be run cell-by-cell in RStudio
- Several scripts share variables via the global R environment (not standalone)
- Geographic analysis scripts have strict run-order dependencies
- Most scripts assume the working directory is the project root (`pvs_in_logan/`)

### When Modifying Code
1. **Read context**: Scripts build on each other - understand data flow and variable dependencies
2. **Preserve thresholds**: The e-value and identity cutoffs are biologically motivated
3. **Check TO_DO.md**: Known bugs are tracked there - avoid introducing duplicates
4. **Mind file paths**: Scripts reference `files/` for input, `outputs/` for generated files
5. **Source utilities**: Scripts using `hmmsearch_clean()` must source `00_utilities/hmmsearch_utils.R`

### Common Tasks

**Adding new sequences to analysis**:
1. Add to appropriate input file in `files/`
2. Re-run clustering step with usearch
3. Update annotation/metadata files
4. Re-run phylogenetic analysis if needed

**Modifying filtering thresholds**:
1. Document biological rationale
2. Track how changes affect downstream counts
3. Consider impact on novelty classification

**Updating visualizations**:
1. Use ggtree for phylogenetic trees
2. Generate PDF/PNG outputs for publication quality
3. Use Polychrome for consistent color schemes with many categories

## Data Scale Reference
- 3400+ SRA libraries screened
- ~2.7 million unique contigs retrieved from Logan
- 307K+ L1-annotated ORFs
- 14,048 full-length L1 candidates
- 386 putatively novel clusters (as of 2026-04, per `test_21k_host.R` three-way classification)

## Caveats and Limitations
- **External dependencies**: Requires installation of bioinformatic tools not tracked in repo
- **Manual curation**: Some steps require manually generated/curated files
- **Large files**: Some data files stored with Git LFS
- **No automated tests**: Pipeline is exploratory - validate outputs manually
- **Path dependencies**: Scripts assume project root as working directory
- **Interactive scripts**: Most scripts share state via global environment, not standalone
- **Known bugs**: See `R_scripts/TO_DO.md` for tracked issues

## Quick Commands Reference

```bash
# Common seqkit operations
seqkit stats input.fa                    # Get sequence statistics
seqkit seq -m 100 input.fa > filtered.fa # Filter by length
seqkit grep -f ids.txt input.fa          # Extract by ID list

# HMMER annotation
hmmsearch --domtblout out.domtbl -T 20 profile.hmm seqs.fa

# Clustering
usearch -cluster_fast input.fa -id 0.9 -centroids out.fa

# Tree building
iqtree2 -s alignment.fa -B 1000 --keep-ident
```

## Contact and Resources
- Repository issues: Document problems in `R_scripts/TO_DO.md`
- Missing files: Check `files/` directory and README for data sourcing instructions
- PAVE database: https://pave.niaid.nih.gov/
- NCBI Virus: https://www.ncbi.nlm.nih.gov/labs/virus/
