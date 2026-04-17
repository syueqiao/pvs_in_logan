Scripts for constructing and iterating on the papillomavirus database from Logan outputs.

- `logan_3.R` — Core pipeline: HMMER filtering, L1 B/I domain identification, novelty search, contig characterization, and coverage analysis.
- `logan_3_novelty.R` — Updated novelty search using blastn against clustered centroids; identifies novel PV clusters at <90% nucleotide identity.

Both scripts source `00_utilities/hmmsearch_utils.R` for HMMER output parsing.
