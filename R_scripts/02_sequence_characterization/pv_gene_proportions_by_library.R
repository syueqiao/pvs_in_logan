library(tidyverse)

# =====================================================================
# Proportion of PV genes (HMM domains) weighted by ka.f, by library type
# Recreates the stacked bar chart: PV gene proportions x SRA library source
# =====================================================================

# --- 1. Read domtbl and extract ka.f + HMM gene from each hit ---
domtbl_raw <- read.table("files/feb_7_pv_fil_form.aa.domtbl",
                         comment.char = "#", fill = TRUE,
                         stringsAsFactors = FALSE)

# Columns: V1 = target name, V4 = HMM query name, V7 = E-value, V8 = score
domtbl <- domtbl_raw %>%
  select(orf_id = V1, hmm = V4, evalue = V7, score = V8) %>%
  # Extract SRR accession from orf_id (everything before first _)
  mutate(srr = sub("^([SED]RR[0-9]+)_.*", "\\1", orf_id)) %>%
  # Extract ka.f value from the header: pattern is _ka_f_<number>
  mutate(ka_f = as.numeric(sub(".*_ka_f_([0-9.]+).*", "\\1", orf_id))) %>%
  # Clean up HMM names to match figure labels
  mutate(hmm_clean = case_when(
    hmm == "Late_protein_L1"      ~ "L1",
    hmm == "Late_protein_L2"      ~ "L2",
    hmm == "PPV_E1_DBD"           ~ "E1_DBD",
    hmm == "PPV_E1_C"             ~ "E1_C",
    hmm == "PPV_E1_N"             ~ "E1_N",
    hmm == "PPV_E2_C"             ~ "E2_C",
    hmm == "PPV_E2_N"             ~ "E2_N",
    hmm == "Pap_E4"               ~ "E4",
    hmm == "Papilloma_E5"         ~ "E5",
    hmm == "E6"                   ~ "E6",
    hmm == "E7"                   ~ "E7",
    hmm == "ncbi_and_novel_L1_JR" ~ "L1",
    TRUE                          ~ hmm
  ))

cat("Parsed", nrow(domtbl), "domain hits from domtbl\n")
cat("Unique ORFs:", n_distinct(domtbl$orf_id), "\n")
cat("Unique SRRs:", n_distinct(domtbl$srr), "\n")

# --- 2. For each ORF, keep only the best-scoring HMM hit ---
# (an ORF may hit multiple HMMs; keep the top one by score)
best_hits <- domtbl %>%
  group_by(orf_id) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Best hits per ORF:", nrow(best_hits), "\n")

# --- 3. Read library source metadata ---
# Use the combined file if available (run fetch_missing_library_source.R first),
# otherwise fall back to lib_source_big.txt
lib_source_file <- ifelse(file.exists("files/lib_source_combined.txt"),
                          "files/lib_source_combined.txt",
                          "files/lib_source_big.txt")
cat("Using library source file:", lib_source_file, "\n")

lib_source <- read.table(lib_source_file, sep = "\t",
                         header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("srr", "library_source"),
                         fill = TRUE, quote = "")

# Clean up and assign display categories matching the figure
lib_source <- lib_source %>%
  filter(library_source != "", !is.na(library_source)) %>%
  mutate(
    lib_type = case_when(
      library_source == "GENOMIC"                  ~ "Genomic",
      library_source == "GENOMIC SINGLE CELL"      ~ "Genomic Single Cell",
      library_source == "METAGENOMIC"              ~ "Metagenomic",
      library_source == "TRANSCRIPTOMIC"           ~ "Transcriptomic",
      library_source == "TRANSCRIPTOMIC SINGLE CELL" ~ "Transcriptomic Single Cell",
      library_source == "VIRAL RNA"                ~ "Viral RNA",
      library_source == "METATRANSCRIPTOMIC"       ~ "Meta-transcriptomic",
      library_source == "SYNTHETIC"                ~ "Synthetic",
      TRUE                                         ~ "Other"
    ),
    lib_group = case_when(
      lib_type %in% c("Genomic", "Genomic Single Cell", "Metagenomic") ~ "DNA",
      lib_type %in% c("Transcriptomic", "Transcriptomic Single Cell",
                       "Viral RNA", "Meta-transcriptomic")             ~ "RNA",
      TRUE                                                              ~ "MISC."
    )
  )

cat("Library source entries:", nrow(lib_source), "\n")

# --- 4. Join hits to library metadata ---
hits_lib <- best_hits %>%
  inner_join(lib_source, by = "srr")

cat("Hits with library info:", nrow(hits_lib), "\n")
cat("Hits without library info (dropped):",
    nrow(best_hits) - nrow(hits_lib), "\n")

# --- 5. Collapse E1/E2 subdomains, then sum ka.f per gene per library type ---
# E1 has 3 HMM profiles (E1_N, E1_DBD, E1_C) and E2 has 2 (E2_N, E2_C),
# while L1/L2/E6/E7 etc. each have 1. To avoid inflating E1/E2's share,
# average ka.f across subdomains within each (lib_type x orf's parent contig).
#
# Approach: for each ORF, assign a parent gene, then within each
# (lib_type, parent_gene, contig), average ka.f across subdomains.

hits_lib <- hits_lib %>%
  mutate(gene_collapsed = case_when(
    hmm_clean %in% c("E1_N", "E1_DBD", "E1_C") ~ "E1",
    hmm_clean %in% c("E2_N", "E2_C")           ~ "E2",
    TRUE                                         ~ hmm_clean
  ))

# Average ka.f across E1/E2 subdomains per contig per library type,
# so a contig with hits to E1_N + E1_DBD + E1_C contributes one E1 value
# (the mean of the three ka.f values, which are often identical since ka.f
# is a contig-level metric). Single-domain genes pass through unchanged.
gene_by_lib <- hits_lib %>%
  # Extract contig ID (SRR + contig number, before _ka_f_)
  mutate(contig_id = sub("_ka_f_.*", "", orf_id)) %>%
  group_by(lib_type, lib_group, gene_collapsed, contig_id) %>%
  summarise(mean_kaf = mean(ka_f, na.rm = TRUE), .groups = "drop") %>%
  # Now sum across contigs
  group_by(lib_type, lib_group, gene_collapsed) %>%
  summarise(total_kaf = sum(mean_kaf, na.rm = TRUE), .groups = "drop") %>%
  # Compute proportions within each lib_type
  group_by(lib_type) %>%
  mutate(proportion = total_kaf / sum(total_kaf)) %>%
  ungroup()

# --- 6. Set factor levels for consistent ordering ---
hmm_levels <- c("L1", "L2", "E7", "E6", "E5", "E4", "E1", "E2")

# Library type order (left to right, grouped by DNA/RNA/MISC)
lib_type_levels <- c("Genomic", "Genomic Single Cell", "Metagenomic",
                     "Transcriptomic", "Transcriptomic Single Cell",
                     "Viral RNA", "Meta-\ntranscriptomic",
                     "Synthetic", "Unknown", "Other")

# Rename Meta-transcriptomic with line break for plotting
gene_by_lib <- gene_by_lib %>%
  mutate(lib_type = ifelse(lib_type == "Meta-transcriptomic",
                           "Meta-\ntranscriptomic", lib_type))

gene_by_lib <- gene_by_lib %>%
  mutate(gene_collapsed = factor(gene_collapsed, levels = rev(hmm_levels)),
         lib_type = factor(lib_type, levels = lib_type_levels))

# --- 7. Color palette: greens for L genes, greyscale for E genes ---
# Consistent with gene_maps_for_pub.R color scheme
hmm_colors <- c(
  "L1" = "#46D452",
  "L2" = "#7cd586",
  "E7" = "#6B7D6E",
  "E6" = "#839886",
  "E5" = "#99A599",
  "E4" = "#ADB8AF",
  "E1" = "#C2CCC4",
  "E2" = "#D6DDD7"
)

# Group annotation bar colors
group_colors <- c("DNA" = "#D98030", "RNA" = "#30B050", "MISC." = "#A0C8E0")

# --- 8. Build the plot ---

# Create group annotation data for colored bars beneath x-axis
group_annot <- gene_by_lib %>%
  distinct(lib_type, lib_group) %>%
  filter(!is.na(lib_type)) %>%
  mutate(lib_group = factor(lib_group, levels = c("DNA", "RNA", "MISC.")))

p <- ggplot(gene_by_lib, aes(x = lib_type, y = proportion, fill = gene_collapsed)) +
  geom_col(color = NA, width = 0.8) +
  # Facet by group with free x scales so DNA/RNA/MISC are visually separated
  facet_grid(~ lib_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = hmm_colors, name = "HMM",
                    breaks = hmm_levels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion",
       title = "Proportion of PV genes from SRA library types") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0.8, "cm"),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

p

ggsave("outputs/pv_gene_proportions_by_library_type.png", p,
       width = 28, height = 18, units = "cm", dpi = 300, bg = "white")

cat("\nSaved: outputs/pv_gene_proportions_by_library_type.png\n")

# --- 9. Chi-squared test: DNA vs RNA gene composition ---
# Aggregate ka.f totals by broad group (DNA vs RNA) and gene
dna_rna <- hits_lib %>%
  filter(lib_group %in% c("DNA", "RNA")) %>%
  mutate(contig_id = sub("_ka_f_.*", "", orf_id)) %>%
  group_by(lib_group, gene_collapsed, contig_id) %>%
  summarise(mean_kaf = mean(ka_f, na.rm = TRUE), .groups = "drop") %>%
  group_by(lib_group, gene_collapsed) %>%
  summarise(total_kaf = sum(mean_kaf, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = lib_group, values_from = total_kaf, values_fill = 0)

# Build contingency table and test
mat <- as.matrix(dna_rna[, c("DNA", "RNA")])
rownames(mat) <- dna_rna$gene_collapsed

chi_test <- chisq.test(mat)
cat("\n--- Chi-squared test: DNA vs RNA gene proportions (all 8 genes) ---\n")
print(chi_test)
cat("p-value:", format.pval(chi_test$p.value, digits = 3), "\n")

# --- 10. E vs L proportion test: DNA vs RNA ---
# Collapse genes into E (E1-E7) vs L (L1, L2), then test whether the
# E:L ratio differs between DNA and RNA libraries
el_table <- hits_lib %>%
  filter(lib_group %in% c("DNA", "RNA")) %>%
  mutate(gene_class = ifelse(grepl("^L", gene_collapsed), "L", "E"),
         contig_id = sub("_ka_f_.*", "", orf_id)) %>%
  group_by(lib_group, gene_class, contig_id) %>%
  summarise(mean_kaf = mean(ka_f, na.rm = TRUE), .groups = "drop") %>%
  group_by(lib_group, gene_class) %>%
  summarise(total_kaf = sum(mean_kaf, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = lib_group, values_from = total_kaf, values_fill = 0)

el_mat <- as.matrix(el_table[, c("DNA", "RNA")])
rownames(el_mat) <- el_table$gene_class

cat("\n--- E vs L proportions by library type ---\n")
cat("Contingency table (ka.f sums):\n")
print(el_mat)

# Proportions
dna_prop <- el_mat[, "DNA"] / sum(el_mat[, "DNA"])
rna_prop <- el_mat[, "RNA"] / sum(el_mat[, "RNA"])
cat(sprintf("\nDNA: %.1f%% E, %.1f%% L\n", dna_prop["E"] * 100, dna_prop["L"] * 100))
cat(sprintf("RNA: %.1f%% E, %.1f%% L\n", rna_prop["E"] * 100, rna_prop["L"] * 100))

# Chi-squared test (appropriate here — large cell counts make Fisher's
# exact test computationally intractable, and chi-squared is well-powered)
el_chi <- chisq.test(el_mat)
cat(sprintf("\nChi-squared test: X² = %.1f, df = %d, p = %s\n",
            el_chi$statistic, el_chi$parameter,
            format.pval(el_chi$p.value, digits = 3)))
