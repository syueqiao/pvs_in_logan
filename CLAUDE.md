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
└── notebooks/                         # Main analysis directory
    ├── 2024.backlog.p1.first_pv_db.ipynb      # Phase 1: Initial PV database
    ├── 2024.backlog.p2.hmms.ipynb             # Phase 2: HMM profile development
    ├── 2025.01.24.second_pv_db.ipynb          # Phase 3: PVDB2 generation
    ├── 2025.02.07.process_logan_output.ipynb  # Phase 4: Logan output processing
    ├── 2025.02.07.statistics_and_interpretation.ipynb  # Phase 5: Data analysis
    ├── 2025.02.14.novelty_search.ipynb        # Phase 6: Novel sequence ID
    ├── 2025.02.20.updated_tree.ipynb          # Phase 7: Phylogenetic trees
    ├── data/                                  # Processed data files
    │   ├── *_envcoords_*.txt                  # Jellyroll domain coordinates
    │   ├── *_metadata.txt                     # SRA and NCBI metadata
    │   ├── blacklist.txt                      # QC-failed sequences
    │   └── *.txt                              # Various intermediate files
    └── *.pdf                                  # Visualization outputs
```

## Technology Stack

### Languages
- **R** (primary): Data manipulation and visualization via tidyverse, ggtree
- **Bash** (secondary): Bioinformatic tool orchestration and file processing
- **Jupyter Notebooks**: Documentation and code execution environment

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
tidyverse       # Data manipulation (dplyr, tidyr, ggplot2)
ggtree          # Phylogenetic tree visualization
Polychrome      # Color palette generation for large trees
```

## Pipeline Phases

### Phase 1: Initial Database Creation (`2024.backlog.p1.first_pv_db.ipynb`)
- Downloads PV sequences from NCBI Virus portal
- Performs ORF finding using stop-stop method
- Runs HMMER against Pfam models (e-value < 10^-9)
- QC against host genomes (human, yeast, mouse, E. coli)
- Uses logan_unsticker to mask contaminants

### Phase 2: HMM Profile Development (`2024.backlog.p2.hmms.ipynb`)
- Downloads PAVE L1 protein sequences
- Generates MSAs for L1 jellyroll domains (B, I, CD, E, F, GH)
- Creates HMM profiles for each domain
- Characterizes sequences for domain completeness

### Phase 3: Database Iteration (`2025.01.24.second_pv_db.ipynb`)
- Processes first Logan run output
- Filters Diamond results (e < 10^-10)
- Combines PVDB1 with novel SRA-derived sequences
- Clusters at 90% AA identity
- Output: `pvdb2_90_final.fa`

### Phase 4-5: Logan Output Processing (`2025.02.07.*.ipynb`)
- Downloads results from AWS S3 (~2.7M unique contigs)
- Filters high-confidence hits (e < 10^-10)
- ORF calling and HMM annotation (307K L1 hits)
- Identifies 14,048 full-length L1 candidates

### Phase 6-7: Novelty and Phylogenetics (`2025.02.14-20.*.ipynb`)
- BLAST filtering (e < 10^-4)
- Identifies < 90% identity sequences as novel
- Clusters 240 novel sequences to 234 representatives
- IQtree phylogenetics (1000 bootstrap replicates)
- Final tree: 1039 sequences (805 NCBI + 234 novel)

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

### Notebook Organization
- **Date-based naming**: `YYYY.MM.DD.phase_description.ipynb`
- **Section markers**: Use "## P1", "## P2", etc. for pipeline phases within notebooks
- **Markdown cells**: Provide biological context and methodology rationale

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

### Bash Code Style
```bash
# Complex piping with comments explaining each step
cat input.fa | \
  seqkit seq -m 100 | \     # Filter minimum length
  seqkit rmdup | \           # Remove duplicates
  awk '{print $1}' > out.fa  # Extract headers
```

### File Naming Conventions
- **Intermediate files**: `descriptive_name_vX.ext` (e.g., `L1_B_super5.hmm`)
- **Data outputs**: `[date]_[content]_[status].ext`
- **Final outputs**: Clear descriptive names with format suffix

### When Modifying Code
1. **Read context**: Each notebook builds on previous ones - understand the data flow
2. **Preserve thresholds**: The e-value and identity cutoffs are biologically motivated
3. **Document changes**: Add markdown cells explaining modifications
4. **Test incrementally**: Notebooks are exploratory - run cells individually
5. **Mind file paths**: Many paths reference `notebooks/data/` - maintain consistency

### Common Tasks

**Adding new sequences to analysis**:
1. Add to appropriate input file in `notebooks/data/`
2. Re-run clustering step with usearch
3. Update annotation/metadata files
4. Re-run phylogenetic analysis if needed

**Modifying filtering thresholds**:
1. Document biological rationale in notebook
2. Track how changes affect downstream counts
3. Consider impact on novelty classification

**Updating visualizations**:
1. Use ggtree for phylogenetic trees
2. Generate PDF outputs for publication quality
3. Use Polychrome for consistent color schemes with many categories

## Data Scale Reference
- 3400+ SRA libraries screened
- ~2.7 million unique contigs retrieved from Logan
- 307K+ L1-annotated ORFs
- 14,048 full-length L1 candidates
- 240 putatively novel sequences (234 after clustering)

## Caveats and Limitations
- **External dependencies**: Requires installation of bioinformatic tools not tracked in repo
- **Manual curation**: Some steps require manually generated/curated files
- **Large files**: Some data files need to be downloaded from external sources (NCBI, AWS S3)
- **No automated tests**: Pipeline is exploratory - validate outputs manually
- **Path dependencies**: Some hardcoded paths may need adjustment for different environments

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
- Repository issues: Document problems in notebooks with clear descriptions
- Missing files: Check `notebooks/data/` folder and README for data sourcing instructions
- PAVE database: https://pave.niaid.nih.gov/
- NCBI Virus: https://www.ncbi.nlm.nih.gov/labs/virus/
