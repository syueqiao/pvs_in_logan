library(tidyverse)
library(gggenes)

# =====================================================================
# Visualize ORF gene maps with top BLAST (nr) hit annotations beneath
# =====================================================================

plot_genome_map <- function(gene_df, contig_id, title_label, genome_length) {
was_flipped <- sum(gene_df$start > gene_df$end) > sum(gene_df$start < gene_df$end)

 if (was_flipped) {
    gene_df <- gene_df %>%
      mutate(new_start = genome_length - start + 1,
             new_end   = genome_length - end   + 1,
             start = new_start, end = new_end) %>%
      select(-new_start, -new_end)
  }
  # Orient the contig so the PV genome reads 5'->3' left-to-right.
  # If most ORFs are encoded on the reverse strand of the contig
  # (start > end in native contig coords), reverse-complement the
  # contig so all ORFs become forward-strand and arrows point right.
  # After this step all ORFs have start < end and the canonical PV
  # reading direction (..., L2, L1, URR, E6, E7, E1, E2, ...) runs L->R.
  if (sum(gene_df$start > gene_df$end) > sum(gene_df$start < gene_df$end)) {
    gene_df <- gene_df %>%
      mutate(new_start = genome_length - start + 1,
             new_end   = genome_length - end   + 1,
             start = new_start, end = new_end) %>%
      select(-new_start, -new_end)
  }

  gene_layer <- gene_df %>%
    transmute(molecule = contig_id, gene, start, end, forward = TRUE)

  # Nudge gene labels left/right when midpoints are too close
  # Uses iterative passes so cascading overlaps get resolved
  gene_label_df <- gene_layer %>%
    mutate(mid = (start + end) / 2) %>%
    arrange(mid)
  gene_label_df$label_x <- gene_label_df$mid
  gene_label_df$label_y <- 1.24
  min_gap <- genome_length * 0.07
  for (pass in 1:5) {
    any_nudge <- FALSE
    for (i in 2:nrow(gene_label_df)) {
      gap <- abs(gene_label_df$label_x[i] - gene_label_df$label_x[i - 1])
      if (gap < min_gap) {
        nudge <- (min_gap - gap) / 2 + min_gap * 0.15
        gene_label_df$label_x[i - 1] <- gene_label_df$label_x[i - 1] - nudge
        gene_label_df$label_x[i]     <- gene_label_df$label_x[i] + nudge
        any_nudge <- TRUE
      }
    }
    if (!any_nudge) break
  }

  blast_df <- gene_df %>%
    filter(!is.na(top_hit), top_hit != "") %>%
    mutate(mid = (start + end) / 2,
           hit_label = paste0(top_hit, "\n", pct_id, "% | ", evalue))

  # Stagger labels: greedily assign y-levels avoiding x-overlap at the same level
  # Approximate label width as fraction of genome span (scales with text length)
  y_levels <- c(0.70, 0.35, 0.00, -0.35)
  blast_df <- blast_df %>% arrange(desc(mid))
  label_width <- nchar(blast_df$top_hit) * genome_length * 0.012  # rough per-char width
  assigned_y <- numeric(nrow(blast_df))
  # Track x-ranges used at each y-level
  used <- lapply(y_levels, function(x) data.frame(lo = numeric(), hi = numeric()))
  for (i in seq_len(nrow(blast_df))) {
    lo_i <- blast_df$mid[i] - label_width[i] / 2
    hi_i <- blast_df$mid[i] + label_width[i] / 2
    placed <- FALSE
    for (j in seq_along(y_levels)) {
      if (nrow(used[[j]]) == 0 || all(lo_i > used[[j]]$hi | hi_i < used[[j]]$lo)) {
        assigned_y[i] <- y_levels[j]
        used[[j]] <- rbind(used[[j]], data.frame(lo = lo_i, hi = hi_i))
        placed <- TRUE
        break
      }
    }
    if (!placed) assigned_y[i] <- y_levels[length(y_levels)] - 0.24 * i
  }
  blast_df$label_y <- assigned_y

  color_match_genes <- c(
    "L1" = "#46D452", "L2" = "#99A599",
    "E1" = "#99A599", "E1_a" = "#99A599", "E1_b" = "#99A599", "E2" = "#99A599",
    "E4" = "#CBD7CB", "E5" = "#CBD7CB", "E6" = "#CBD7CB", "E7" = "#CBD7CB"
  )

  p <- ggplot(gene_layer, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = TRUE)) +
  geom_segment(aes(x = 1, xend = genome_length, y = molecule, yend = molecule),
               inherit.aes = FALSE,
               data = data.frame(molecule = contig_id, genome_length = genome_length),
               color = "grey60", linewidth = 0.5) +
               geom_segment(data = data.frame(molecule = contig_id, genome_length = genome_length),
               aes(x = 1, xend = genome_length, y = -Inf, yend = -Inf),
               inherit.aes = FALSE,
               color = "grey70", linewidth = 0.3) +
    geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),
                    arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
    geom_text(data = gene_label_df,
              aes(x = label_x, y = label_y, label = gene),
              inherit.aes = FALSE,
              fontface = "bold", family = "Noto Sans", size = 5) +
    geom_segment(data = blast_df,
                 aes(x = mid, xend = mid, y = 0.93, yend = label_y + 0.05),
                 inherit.aes = FALSE, linetype = "dotted", color = "grey60", linewidth = 0.3) +
    geom_label(data = blast_df,
               aes(x = mid, y = label_y, label = hit_label),
               inherit.aes = FALSE,
               size = 3.5, family = "Noto Sans", lineheight = 0.85,
               fill = "white", label.size = 0.2, label.padding = unit(2, "pt"),
               color = "grey30") +
    scale_fill_manual(values = color_match_genes) +
    scale_x_continuous(
    expand = expansion(mult = 0.08),
    breaks = if (was_flipped) {
      original_breaks <- seq(0, genome_length, by = 1000)
      genome_length - original_breaks + 1
    } else {
      seq(0, genome_length, by = 1000)
    },
    labels = if (was_flipped) {
      function(x) format(round(genome_length - x + 1), big.mark = ",")
    } else {
      function(x) format(round(x), big.mark = ",")
    }
  ) +
    # Title anchored at the left edge of the genome (x = 1 on normal axis)
    annotate("text", x = 1, y = 1.50,
             label = paste0(contig_id, ": ", title_label, "  (",
                            format(genome_length, big.mark = ","), " nt)"),
             hjust = 0, fontface = "bold", family = "Noto Sans", size = 5) +
    coord_cartesian(ylim = c(min(blast_df$label_y) - 0.1, 1.42), clip = "off") +
    labs(x = NULL, y = NULL) +
    guides(fill = "none") +
    theme_bw() +
    theme(text = element_text(family = "Noto Sans", size = 12),
          plot.margin = margin(2, 16, 2, 16, "pt"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_blank())

  return(p)
}

# =====================================================================
# --- Pangolin PV: SRR25256564_207500 ---
# =====================================================================
pangolin_blast <- tribble(
  ~gene,  ~start, ~end,  ~top_hit,                                       ~accession,       ~evalue,   ~pct_id, ~coverage,
  "E6",   7196,   6612,  "Phocoena phocoena papillomavirus 4",           "YP_006470634.1", "9e-35",   47.01,   69,
  "E1",   6496,   4619,  "Tursiops truncatus papillomavirus 8",          "AVH76832.1",     "0.0",     52.11,   97,
  "E2",   4698,   3310,  "Papillomavirus meriones7946",                  "DBA47531.1",     "1e-54",   49.47,   57,
  "L2",   3145,   1730,  "Phocoena phocoena papillomavirus 4",           "YP_006470638.1", "2e-49",   41.63,   94,
  "L1",   1950,   265,   "Manis pentadactyla papillomavirus 1",          "DAZ92280.1",     "0",       77.37,   86
)

pangolin_map <- plot_genome_map(
  pangolin_blast,
  contig_id    = "SRR25256564_207500",
  title_label  = "White-bellied pangolin PV",
  genome_length = 7200
)

# pangolin_map
# # ggsave("outputs/pangolin_genome_map.png", pangolin_map,
#        width = 28, height = 8, units = "cm", dpi = 300, bg = "white")

# =====================================================================
# --- Rhinoceros PV: SRR10902309_46767 ---
# =====================================================================
rhino_blast <- tribble(
  ~gene,  ~start, ~end,  ~top_hit,                                       ~accession,       ~evalue,   ~pct_id, ~coverage,
  "E6",   2,      556,   "Bos taurus papillomavirus",                    "XZE12012.1",     "5e-12",   38.46,   63,
  "E1",   766,    2631,  "Equus asinus papillomavirus 1",                "YP_009021883.1", "9e-177",  50.59,   93,
  "E2",   2537,   3790,  "Equus asinus papillomavirus 2",                "QBN12682.1",     "9e-71",   37.14,   96,
  "E4",   3123,   3518,  "Deltapapillomavirus 5",                        "ALO50064.1",     "3e-04",   40.00,   73,
  "L2",   4033,   5646,  "Equus caballus papillomavirus 2",              "ADZ73318.1",     "7e-34",   34.45,   98,
  "L1",   5543,   7126,  "Bat papillomavirus",                           "WXG28265.1",     "1e-158",  47.05,   95
)

rhino_map <- plot_genome_map(
  rhino_blast,
  contig_id    = "SRR10902309_46767",
  title_label  = "Rhinoceros PV",
  genome_length = 7126
)

# rhino_map
# ggsave("outputs/rhino_genome_map.png", rhino_map,
      #  width = 28, height = 10, units = "cm", dpi = 300, bg = "white")

# =====================================================================
# --- Lizard PV: SRR22028468_199653 ---
# =====================================================================
lizard_blast <- tribble(
  ~gene,  ~start, ~end,  ~top_hit,                                       ~accession,       ~evalue,   ~pct_id, ~coverage,
  "E1",   1133,   3097,  "Papillomaviridae sp.",                         "XQE69105.1",     "3e-144",  40.41,   98,
  "E2",   3069,   4085,  "Papillomaviridae sp.",                         "XQE69106.1",     "2e-27",   32.37,   71,
  "E4",   3181,   3408,  "Serinus canaria papillomavirus 1",             "YP_009551919.1", "2e-17",   50.98,   29,
  "L2",   4078,   5481,  "Rhinolophus bat papillomavirus 1",             "WXG28054.1",     "8e-04",   50.00,   12,
  "L1",   5441,   6958,  "Macaca mulatta papillomavirus 3",              "YP_009553631.1", "4e-137",  43.88,   93
)

lizard_map <- plot_genome_map(
  lizard_blast,
  contig_id    = "SRR22028468_199653",
  title_label  = "Lizard PV",
  genome_length = 6958
)

# lizard_map
# ggsave("outputs/lizard_genome_map.png", lizard_map,
#        width = 28, height = 8, units = "cm", dpi = 300, bg = "white")

# =====================================================================
# --- Human PV: SRR13789839_2669 ---
# =====================================================================
human_blast <- tribble(
  ~gene,   ~start, ~end,  ~top_hit,                                      ~accession,       ~evalue,   ~pct_id, ~coverage,
  "E1_a",  7333,   7148,  "Human papillomavirus",                        "AYA94407.1",     "1e-29",   87.10,   100,
  "E1_b",  1673,   3,     "Gammapapillomavirus 7",                       "ATQ38256.1",     "0.0",     87.31,   96,
  "E7",    1887,   1591,  "Gammapapillomavirus 7",                       "ATQ38255.1",     "1e-34",   83.16,   96,
  "E6",    2342,   1878,  "Gammapapillomavirus 7",                       "ATQ38254.1",     "3e-75",   78.99,   89,
  "L1",    4422,   2782,  "Human papillomavirus",                        "AYA94410.1",     "0.0",     91.11,   99,
  "L2",    5989,   4349,  "Human papillomavirus type 203",               "YP_009552753.1", "0.0",     93.97,   97,
  "E4",    6734,   6204,  "Human papillomavirus",                        "AYA94104.1",     "1e-31",   57.42,   87,
  "E2",    7215,   5962,  "Human papillomavirus",                        "AYA94408.1",     "0.0",     76.16,   98
)

human_map <- plot_genome_map(
  human_blast,
  contig_id    = "SRR13789839_2669",
  title_label  = "Human PV",
  genome_length = 7333
)

# human_map
# ggsave("outputs/human_genome_map.png", human_map,
      #  width = 28, height = 8, units = "cm", dpi = 300, bg = "white")

# =====================================================================
# --- Salmon PV: SRR20078264_4021 ---
# =====================================================================
salmon_blast <- tribble(
  ~gene,   ~start, ~end,  ~top_hit,                                      ~accession,       ~evalue,   ~pct_id, ~coverage,
  "E1_a",  2,      418,   "Papillomavirus astyanax5827",                 "DBA48550.1",     "3e-55",   65.69,   99,
  "E2",    411,    1565,  "Papillomaviridae sp. Haddock_c45",            "QYW06004.1",     "5e-65",   37.82,   97,
  "L2",    1504,   2898,  "Papillomavirus gadus6361",                    "DBA49308.1",     "2e-25",   29.75,   92,
  "L1",    2870,   4516,  "Papillomaviridae sp. Haddock_c45",            "QYW06007.1",     "0.0",     55.89,   100,
  "E1_b",  4545,   6110,  "Papillomaviridae sp. Haddock_c2655",         "QYW06014.1",     "1e-83",   40.37,   80
)

salmon_map <- plot_genome_map(
  salmon_blast,
  contig_id    = "SRR20078264_4021",
  title_label  = "Salmon PV",
  genome_length = 6110
)

# salmon_map
# ggsave("outputs/salmon_genome_map.png", salmon_map,
      #  width = 28, height = 8, units = "cm", dpi = 300, bg = "white")

# =====================================================================
# --- Combined figure: all genome maps A-E ---
# =====================================================================
library(ggpubr)

# Add bottom margin to each panel for spacing between rows
add_spacing <- function(p) p + theme(plot.margin = margin(2, 5, 15, 5, "pt"))

combined_genome_maps <- ggarrange(
  add_spacing(pangolin_map), add_spacing(lizard_map), add_spacing(rhino_map),
  add_spacing(human_map), add_spacing(salmon_map),
  ncol = 1, nrow = 5,
  labels = c("A", "B", "C", "D", "E"),
  font.label = list(size = 16, family = "Noto Sans")
)

ggsave("outputs/genome_maps_combined.png", combined_genome_maps,
       width = 32, height = 42, units = "cm", dpi = 300, bg = "white")
