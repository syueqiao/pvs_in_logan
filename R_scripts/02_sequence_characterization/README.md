Scripts for characterizing PV sequences recovered from Logan and NCBI.

- `investigation_of_characteristics.R` — L1 domain inversion analysis, domain coverage, and missing B/I domain detection across Logan contigs.
- `ncbi_L1_character.R` — Characterization of NCBI L1 sequences; verifies domain order and completeness.
- `L1_pathracer_mining.R` — Pathracer pipeline setup for recovering full-length L1 genes from short-read data.
- `pathracer_fullness_hmm.R` — Analysis of Pathracer output: HMM domain completeness assessment.
- `pv_gene_proportions_by_library.R` — PV gene (E1–E7, L1, L2) proportions weighted by ka.f across SRA library types. Collapses E1/E2 subdomains via averaging. Includes chi-squared test for DNA vs RNA gene composition.
- `fetch_missing_library_source.R` — Batch queries the ENA Portal API to retrieve library source metadata for accessions missing from local files. Merges with existing metadata and saves combined output.
- `vert_counts.R` — Vertebrate PV counts vs sequencing depth analysis (standalone, no dependencies).
