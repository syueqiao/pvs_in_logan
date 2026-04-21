library(tidyverse)
library(sf)
library(extrafont)

# Requires objects from geo_analysis.R (which sources geo_setup.R and loads polygon_analyses.R outputs):
#   - sf_object_joined_prop_rarefy_5k_rowmeans_joined: data_wide + resampling z_score/p_val
#   - giant_geo_table_grid_id_geometry: grid cells with Freq (SRA sampling effort)
#   - data_wide: grid_group, known, novel, novel_species_proportion

sf_object_joined_prop_rarefy_5k_rowmeans_joined <- readRDS(
  "outputs/geo_intermediates/sf_object_joined_prop_rarefy_5k_rowmeans_joined_s.rds"
)
sf_object_joined_prop_rarefy_5k_rowmeans_joined_moran <- dplyr::select(sf_object_joined_prop_rarefy_5k_rowmeans_joined, geometry, novel_species_proportion, z_score, mean, p_val, grid_group, known, novel)

# Build unified grid dataframe by joining sampling effort (Freq) to the resampling results
grid_df <- sf_object_joined_prop_rarefy_5k_rowmeans_joined_moran %>%
  st_drop_geometry() %>%
  as.data.frame()

# cat("Columns in grid_df:", paste(names(grid_df), collapse = ", "), "\n")

# Add SRA sampling effort from giant_geo_table_grid_id_geometry
# Handle possible .x suffix from prior joins in session
freq_cols <- names(giant_geo_table_grid_id_geometry)
freq_col <- if ("Freq" %in% freq_cols) "Freq" else if ("Freq.x" %in% freq_cols) "Freq.x" else stop("No Freq column found")
var_col <- if ("Var1" %in% freq_cols) "Var1" else stop("No Var1 column found")

freq_lookup <- giant_geo_table_grid_id_geometry %>%
  st_drop_geometry() %>%
  as.data.frame()
freq_lookup <- data.frame(Var1 = freq_lookup[[var_col]], Freq = freq_lookup[[freq_col]])

grid_df$grid_group <- as.numeric(paste(grid_df$grid_group))
grid_df <- left_join(grid_df, freq_lookup, by = c("grid_group" = "Var1"))

cat("Rows in grid_df:", nrow(grid_df), "\n")
cat("Freq NAs:", sum(is.na(grid_df$Freq)), "/", nrow(grid_df), "\n")
cat("p_val NAs:", sum(is.na(grid_df$p_val)), "/", nrow(grid_df), "\n")

# Filter to grid squares that actually have PVs
grid_pv <- grid_df %>%
  dplyr::filter((known + novel) > 0)

# Flag enriched squares using uncorrected p_val from permutation test
grid_pv <- grid_pv %>%
  mutate(enriched = p_val < 0.05)

cat("Grid squares with PVs:", nrow(grid_pv), "\n")
cat("Enriched (p_val < 0.05):", sum(grid_pv$enriched, na.rm = TRUE), "\n")


# Get grid cell centroids with latitude
grid_with_geom <- st_sf(grid_pv, geometry = st_geometry(grid_sf)[grid_pv$grid_group])


# ---- Option 8: Moran's I and Getis-Ord Gi* on original 50km grid ----
# Tests spatial autocorrelation of novel enrichment across all PV-positive grid cells
# This formalizes "enrichment clusters in East Africa and South America" into a real statistic

library(spdep)

# Use the grid_with_geom sf object (50km grid cells with PVs)
# kNN (k=8) — appropriate for unevenly distributed data; every cell gets neighbors
grid_coords <- st_coordinates(st_centroid(grid_with_geom))
knn_grid    <- knearneigh(grid_coords, k = 16)
nb_grid     <- knn2nb(knn_grid)
wts_grid    <- nb2listw(nb_grid, style = "W")
wts_grid_B  <- nb2listw(nb_grid, style = "B")

# 8a: Moran's I on novel species proportion (continuous)
library(spdep)
moran_prop <- moran.test(grid_pv$novel_species_proportion, wts_grid, na.action = na.omit)
cat("\n==== Moran's I on 50km grid (novel species proportion) ====\n")
cat(paste0("  I = ", round(moran_prop$estimate["Moran I statistic"], 3),
           ", p = ", format.pval(moran_prop$p.value, digits = 3), "\n"))

# 8b: Moran's I on z-scores from permutation test
moran_z <- moran.test(grid_pv$z_score, wts_grid, na.action = na.omit)
cat(paste0("\nMoran's I on z-scores: I = ", round(moran_z$estimate["Moran I statistic"], 3),
           ", p = ", format.pval(moran_z$p.value, digits = 3), "\n"))

# 8c: Getis-Ord Gi* hotspot analysis
# Identifies individual grid cells that are significant hot/cold spots
# Requires binary weights (not row-standardized)
gi_star <- localG(grid_pv$novel_species_proportion, wts_grid_B)
grid_pv$gi_star <- as.numeric(gi_star)
grid_pv$gi_pval <- 2 * pnorm(abs(grid_pv$gi_star), lower.tail = FALSE)  # two-tailed
grid_pv$hotspot <- NULL
grid_pv$hotspot <- case_when(
  grid_pv$gi_star > 1.96 & grid_pv$novel_species_proportion > 0 ~ "Hotspot (high novelty)",
  grid_pv$gi_star < -1.96                                        ~ "Coldspot (low novelty)",
  TRUE                                                           ~ "Not significant"
)

cat("\nGetis-Ord Gi* hotspot summary:\n")
print(table(grid_pv$hotspot))

# Map hotspots
hotspot_sf <- st_sf(grid_pv, geometry = st_geometry(grid_sf)[grid_pv$grid_group])

p8 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey90") +
  geom_sf(data = st_buffer(filter(hotspot_sf, hotspot == "Not significant"), 50000),
          fill = "grey70", color = NA, alpha = 0.3) +
  geom_sf(data = st_buffer(filter(hotspot_sf, hotspot == "Hotspot (high novelty)"), 50000),
          fill = "#3333E7", color = NA, alpha = 0.8) +
  geom_sf(data = st_buffer(filter(hotspot_sf, hotspot == "Coldspot (low novelty)"), 50000),
          fill = "#3333E7", color = NA, alpha = 0.8) +
  coord_sf(crs = "ESRI:54009") +
  labs(title = paste0("Getis-Ord Gi* hotspot analysis of novel PV enrichment\n",
                      "(Moran's I = ", round(moran_prop$estimate["Moran I statistic"], 3),
                      ", p = ", format.pval(moran_prop$p.value, 3), ")")) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Noto Sans"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

print(p8)



ggsave("outputs/hotspot_map.png", p8, width = 30, height = 18, units = "cm", bg = "transparent")


# ---- Option 9: Hotspot hull map ----
# Draws convex hulls around spatial clusters of Gi* hotspot cells
# instead of plotting individual grid squares

library(dbscan)

# Extract hotspot cell centroids
hotspot_cells <- hotspot_sf %>% filter(hotspot == "Hotspot (high novelty)")
hotspot_pts <- st_centroid(hotspot_cells)

# Cluster nearby hotspot cells using DBSCAN (eps = max gap between cluster members)
coords_hs <- st_coordinates(hotspot_pts)
clust <- dbscan(coords_hs, eps = 500000, minPts = 2)  # 1000km gap tolerance
hotspot_pts$cluster <- clust$cluster

# Drop noise points (cluster == 0 means not part of any cluster)
hotspot_pts_clustered <- hotspot_pts %>% filter(cluster > 0)

# Draw rounded hulls: union points per cluster, buffer to create smooth blob shapes
hotspot_hulls <- hotspot_pts_clustered %>%
  group_by(cluster) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_buffer(500000)  # 1000km buffer creates rounded blob shapes

p9 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey90") +
  geom_sf(data = hotspot_hulls, fill = "#3333E7", color = NA,
          alpha = 0.2, linewidth = 0) +
  # geom_sf(data = coldspot_hulls, fill = "#3333E7", color = "#3333E7",
  #         alpha = 0.3, linewidth = 0.8) +
  # Overlay individual hotspot/coldspot points for reference
  geom_sf(data = st_centroid(filter(hotspot_sf, hotspot == "Hotspot (high novelty)")),
          color = "#3333E7", size = 0.5, alpha = 0.8) +
  coord_sf(crs = "ESRI:54009") +
  labs(title = paste0("Gi* hotspot/coldspot clusters (convex hulls)\n",
                      "(Moran's I = ", round(moran_prop$estimate["Moran I statistic"], 3),
                      ", p = ", format.pval(moran_prop$p.value, 3), ")")) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Noto Sans"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

print(p9)
ggsave("outputs/hotspot_hull_map.png", p9, width = 30, height = 18, units = "cm", bg = "transparent")

#check k appropriateness

for (k in c(4, 8, 12, 16, 20, 24, 28, 32, 48)) {
  knn <- knearneigh(grid_coords, k = k)
  nb <- knn2nb(knn)
  wts <- nb2listw(nb, style = "W")
  mi <- moran.test(grid_pv$z_score, wts, na.action = na.omit)
  cat("k =", k, "  I =", round(mi$estimate[1], 4), 
      "  p =", formatC(mi$p.value, format = "e", digits = 2), "\n")
}

#final overlay


sf_novel_sorted_b  <- st_as_sf(arrange(sf_object_joined_prop_rarefy_5k_rowmeans_joined, novel_species_proportion))
sf_novel_zero_b    <- filter(sf_novel_sorted_b, novel_species_proportion == 0)
sf_novel_nonzero_b <- filter(sf_novel_sorted_b, novel_species_proportion > 0)

vals_prop_all_overlay <- ggplot() +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  # Bottom layer: sampling density (red)
  geom_sf(data = giant_geo_table_grid_id_geometry,
          aes(fill = log(Freq)), lwd = 0, color = "transparent") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "OrRd"), na.value = "transparent") +
  ggnewscale::new_scale_fill() +
  # Zero novelty: sampled but no novel PVs found
  geom_sf(data = st_buffer(sf_novel_zero_b, 75000),
          fill = "grey60", color = "white", lwd = 0.2, alpha = 0.6) +
  # Non-zero novelty: blue gradient, higher values drawn on top
  geom_sf(data = st_buffer(sf_novel_nonzero_b, 100000),
          aes(fill = novel_species_proportion), lwd = 0.2, color = "white", alpha = 0.8) +
  scale_fill_gradientn(colors = c("#86b3c6", "#3333E7"), na.value = "transparent") +
  coord_sf(crs = "ESRI:54009") +
  theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none')

plot_sample_n_log10_hulls <- vals_prop_all_overlay +
  geom_sf(data = hotspot_hulls,
          fill = "#3333E7", color = NA, alpha = 0.2, linewidth = 0) +
  geom_sf(data = st_centroid(filter(hotspot_sf, hotspot == "Hotspot (high novelty)")),
          color = "#3333E7", size = 0, alpha = 0)

plot_sample_n_log10_hulls

ggsave("outputs/2026.02.24.plot_sample_n_log10_hulls.png", plot_sample_n_log10_hulls,
       width = 20, height = 20, units = "cm", limitsize = FALSE, bg = "transparent", dpi = 600)


# Reference 100km circle around the first hotspot point for scale
ref_pt     <- st_centroid(filter(hotspot_sf, hotspot == "Hotspot (high novelty)"))[1, ]
ref_circle <- st_buffer(ref_pt, 2500000)

p8 +
  geom_sf(data = ref_circle, fill = NA, color = "red", linewidth = 0.8, linetype = "dashed")
