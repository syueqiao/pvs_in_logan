# novelty_significance_maps.R — Novelty proportion significance map
# Extracted from polygon_analyses.R (was lines 269-303)
# Run AFTER: polygon_analyses.R AND geo_analysis.R
#
# Dependencies:
#   - data_wide (from polygon_analyses.R)
#   - sf_object_joined_prop_rarefy_5k_rowmeans_joined_s (from geo_analysis.R)
#   - world (from geo_setup.R)

source("R_scripts/04_geographic_analysis/geo_setup.R")

if (!exists("data_wide")) {
  data_wide <- readRDS("outputs/geo_intermediates/data_wide.rds")
}
if (!exists("sf_object_joined_prop_rarefy_5k_rowmeans_joined_s")) {
  sf_object_joined_prop_rarefy_5k_rowmeans_joined_s <- readRDS(
    "outputs/geo_intermediates/sf_object_joined_prop_rarefy_5k_rowmeans_joined_s.rds")
}


sf_object_joined_prop_rarefy_5k_rowmeans_joined_s

plot_novel_proportion <- ggplot(data = st_buffer(st_as_sf(sf_object_joined_prop_rarefy_5k_rowmeans_joined_s), 50000), aes(fill = novel_species_proportion)) +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(shape=16, color = NA, alpha = 0.8) +
  scale_fill_gradientn(colors = c("grey70", "#3333E7")) +
  # geom_sf (data =borders,
  #         aes(geometry = geometry), color = "#E7338D", fill = NA) +
  coord_sf(crs = "ESRI:54009") +
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

plot_novel_proportion

ggsave("outputs/plot_novel_proportion_2.png", plot_novel_proportion, width = 30, height = 20, units = "cm", limitsize = F,bg='transparent')
