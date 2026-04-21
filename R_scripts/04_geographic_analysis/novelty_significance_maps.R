# novelty_significance_maps.R — Novelty proportion significance map
# Extracted from polygon_analyses.R (was lines 269-303)
# Run AFTER: polygon_analyses.R AND geo_analysis.R
#
# Dependencies:
#   - data_wide (from polygon_analyses.R)
#   - sf_object_joined_prop_rarefy_5k_rowmeans_joined_s (from geo_analysis.R)
#   - world (from geo_setup.R)

source("R_scripts/04_geographic_analysis/geo_setup.R")
library(ggnewscale)

if (!exists("data_wide")) {
  data_wide <- readRDS("outputs/geo_intermediates/data_wide.rds")
}
if (!exists("sf_object_joined_prop_rarefy_5k_rowmeans_joined_s")) {
  sf_object_joined_prop_rarefy_5k_rowmeans_joined_s <- readRDS(
    "outputs/geo_intermediates/sf_object_joined_prop_rarefy_5k_rowmeans_joined_s.rds")
}
if (!exists("giant_geo_table_grid_id_geometry")) {
  giant_geo_table_grid_id_geometry <- readRDS(
    "outputs/geo_intermediates/giant_geo_table_grid_id_geometry.rds")
}

# Significance borders: grid cells with p < 0.05
borders <- st_buffer(
  st_as_sf(filter(sf_object_joined_prop_rarefy_5k_rowmeans_joined_s, p_val < 0.05)),
  100000
)

plot_novel_proportion <- ggplot() +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  # Bottom layer: SRA sampling density (orange)
  geom_sf(data = giant_geo_table_grid_id_geometry,
          aes(fill = log(Freq)), lwd = 0, color = "transparent") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "OrRd"),
                       na.value = "transparent", name = "log(SRA libraries)") +
  new_scale_fill() +
  # Middle layer: novelty proportion (blue) with dark borders and transparency
  geom_sf(data = st_buffer(st_as_sf(sf_object_joined_prop_rarefy_5k_rowmeans_joined_s), 100000),
          aes(fill = novel_proportion), color = "#ffffff", linewidth = 0.3, alpha = 0.7) +
  scale_fill_gradientn(colors = c("#343434", "#8EB0C8", "#7573DA", "#3333E7"),
                       name = "Proportion of novel types",
                       breaks = c(0, 0.25, 0.5, 1),
                       limits = c(0, 1)) +
  # Top layer: significance borders (p < 0.05)
  geom_sf(data = borders,
          aes(geometry = geom), color = "#00ff99", fill = NA, linewidth = 0.6) +
  coord_sf(crs = "ESRI:54009") +
  theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

plot_novel_proportion

ggsave("outputs/plot_novel_proportion_2.png", plot_novel_proportion,
       width = 30, height = 20, units = "cm", limitsize = F, bg = 'transparent')
