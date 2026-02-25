#cleaned geo analysis
source("R_scripts/04_geographic_analysis/geo_setup.R")
library(iNEXT)
library(ggplot2)
library(raster)
library(reshape2)
library(viridis)
library(matrixStats)
library(extrafont)
library(paletteer)

# Load from environment (if polygon_analyses.R was run) or from saved RDS
if (!exists("giant_geo_table_grid_id_geometry")) {
  giant_geo_table_grid_id_geometry <- readRDS("outputs/geo_intermediates/giant_geo_table_grid_id_geometry.rds")
}
if (!exists("data_wide")) {
  data_wide <- readRDS("outputs/geo_intermediates/data_wide.rds")
}

#by resampling the geo data, we can generate a "mean" of whatever parameter we are interested in
#requires the sf_object information generated in polygon_analyses.R

##resampling
#10k iterations, or else takes too long
 

giant_geo_table_grid_id_geometry$Var1 <- as.numeric(paste0(giant_geo_table_grid_id_geometry$Var1))

aa <- st_sf(Var1 = seq_len(nrow(grid_sf)), geometry = st_geometry(grid_sf))

giant_geo_table_grid_id_geometry <- left_join(aa, st_drop_geometry(giant_geo_table_grid_id_geometry), by = "Var1")

giant_geo_table_grid_id_geometry[is.na(giant_geo_table_grid_id_geometry)] <- 0


#set arbitrary seed
set.seed(102938501)

sf_object_joined_sampling_rarefy_5k <- as.data.frame(giant_geo_table_grid_id_geometry$Freq)
sf_object_joined_sampling_rarefy_5k[is.na(sf_object_joined_sampling_rarefy_5k)] <- 0
colnames(sf_object_joined_sampling_rarefy_5k) <- c("V1")
sf_object_joined_sampling_rarefy_5k$V1 <- as.numeric(sf_object_joined_sampling_rarefy_5k$V1)

for (i in 1:2000) {
  # Randomly sample one value from x and assign it to y
  sf_object_joined_sampling_rarefy_5k[[i]] <- sample(sf_object_joined_sampling_rarefy_5k$V1, replace = F)
}


sf_object_joined_sampling_rarefy_5k_rowmeans <- sf_object_joined_sampling_rarefy_5k %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_longer(-row) %>%
  group_by(row) %>%   
  summarize(mean = mean(value),
            sd = sd(value)) %>%
  bind_cols(sf_object_joined_sampling_rarefy_5k, .)


sf_object_joined_sampling_rarefy_5k_rowmeans <- cbind(dplyr::select(giant_geo_table_grid_id_geometry, geometry, Var1, Freq) , sf_object_joined_sampling_rarefy_5k_rowmeans)

#calculate stat values#calculate stat valuesNULL
sf_object_joined_sampling_rarefy_5k_rowmeans$z_score <- (as.numeric(sf_object_joined_sampling_rarefy_5k_rowmeans$Freq) - sf_object_joined_sampling_rarefy_5k_rowmeans$mean) / sf_object_joined_sampling_rarefy_5k_rowmeans$sd
sf_object_joined_sampling_rarefy_5k_rowmeans$p_val <- pnorm(sf_object_joined_sampling_rarefy_5k_rowmeans$z_score, mean = 0, sd = 1, lower.tail = F)


#smaller df to look at manually
sf_object_joined_sampling_rarefy_5k_rowmeans_joined_zs <- dplyr::select(sf_object_joined_sampling_rarefy_5k_rowmeans, geometry, Freq, z_score, p_val, mean)
# sf_object_joined_summed_by_long_samples_joined_zs$p_val <- pnorm(sf_object_joined_summed_by_long_samples_joined_zs$z_score, mean = 0, sd = 1, lower.tail = F)

#graph the ones below a certain cutoff value set by cut()

sf_object_joined_sampling_rarefy_5k_rowmeans_joined_zs_fil <- filter(sf_object_joined_sampling_rarefy_5k_rowmeans_joined_zs, p_val < 0.05)

p_vals <- ggplot(data = sf_object_joined_sampling_rarefy_5k_rowmeans_joined_zs_fil) +
  geom_sf(data = world, fill = 'grey80', color = 'grey70') +
  geom_sf(fill = "red", color = "transparent", aes(fill = 'red', color = "red", geometry = geometry)) + 
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme_bw() + theme(legend.position = 'none')

ggsave("outputs/p_vals_test.png", p_vals, width = 20, height = 15, units = "cm", limitsize = F)


mean_sampling_dense <- ggplot(data = sf_object_joined_sampling_rarefy_5k_rowmeans, aes(fill = 'transparent', color = "transparent")) +
  geom_sf(data = world, fill = 'grey80', color = 'grey70') +
  geom_sf(aes(fill = mean, color = "transparent", geometry = geometry), lwd = 0,  na.value = "transparent", color = NA) +
  scale_fill_gradientn(colors = c(viridis(option = 'plasma', n = 10, direction = -1)), na.value = "transparent") +
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme_bw() + theme(legend.position = 'none')

ggsave("outputs/mean_sampling_dense_test.png", mean_sampling_dense, width = 20, height = 15, units = "cm", limitsize = F)

#can also use the "static" means and sds, but will not for this analysis
# mean <- sum(sf_object_joined_summed_by_long_samples_joined$sample_n)/length(sf_object_joined_summed_by_long_samples_joined$polygon_long)
#write to a df the
# long_values <- sf_object_joined_summed_by_long_samples_joined_zs
###############################################################################################do the same for the latitude################################################################################################

#set arbitrary seed
set.seed(102938501)

sf_object_joined_novel_prop_rarefy_5k <- as.data.frame(data_wide$novel_species_proportion)
colnames(sf_object_joined_novel_prop_rarefy_5k) <- c("novel_proportion")
for (i in 1:2000) {
  # Randomly sample one value from x and assign it to y
  sf_object_joined_novel_prop_rarefy_5k[[i]] <- sample(sf_object_joined_novel_prop_rarefy_5k$novel_proportion, replace = F)
}

sf_object_joined_novel_prop_rarefy_5k[is.na(sf_object_joined_novel_prop_rarefy_5k)] <- 0

#add original information, since the first two columns get overwritten
# sf_object_joined_summed_by_long_samples_joined <- cbind(select(sf_object_joined_sampling_rarefy,sample_n, geometry), sf_object_joined_sampling_rarefy_5k)

#calc the row means, so for every grid point, an average of 10k random assignments of a y value
sf_object_joined_novel_prop_rarefy_5k_rowmeans <- sf_object_joined_novel_prop_rarefy_5k %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_longer(-row) %>%
  group_by(row) %>%
  summarize(mean = mean(value),
            sd = sd(value)) %>%
  bind_cols(sf_object_joined_novel_prop_rarefy_5k, .)

sf_object_joined_novel_prop_rarefy_5k_rowmeans[is.na(sf_object_joined_novel_prop_rarefy_5k_rowmeans)] <- 0
# sf_object_joined_prop_rarefy$novel_proportion <- sf_object_joined_prop_rarefy$species_obvs_novel/sf_object_joined_prop_rarefy$species_obvs_all

sf_object_joined_prop_rarefy_5k_rowmeans_joined <- cbind(data_wide, sf_object_joined_novel_prop_rarefy_5k_rowmeans)



# Calculate row standard deviations
#ignore the row means column that was just generated, and the original information columns
#plot the distribution over 10k samples
#generate some nice colors
# sf_object_joined_sampling_rarefy_5k_rowmeans_joined$id_col <- c(1:length(sf_object_joined_sampling_rarefy_5k_rowmeans_joined$V2))

# mean_graphs <- ggplot(sf_object_joined_sampling_rarefy_5k_rowmeans_joined, aes(y=mean, x = id_col)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.3) +
#   theme_classic()# Apply the classic theme
# ggsave("mean_graphs.png", mean_graphs, width = 15, height = 15, units = "cm", limitsize = F)


#calculate stat values
sf_object_joined_prop_rarefy_5k_rowmeans_joined$z_score <- (sf_object_joined_prop_rarefy_5k_rowmeans_joined$novel_species_proportion - sf_object_joined_prop_rarefy_5k_rowmeans_joined$mean) / sf_object_joined_prop_rarefy_5k_rowmeans_joined$sd
sf_object_joined_prop_rarefy_5k_rowmeans_joined$p_val <- pnorm(sf_object_joined_prop_rarefy_5k_rowmeans_joined$z_score, mean = 0, sd = 1, lower.tail = F)

#smaller df to look at manually
sf_object_joined_prop_rarefy_5k_rowmeans_joined_s <- dplyr::select(sf_object_joined_prop_rarefy_5k_rowmeans_joined, geometry, novel_species_proportion, z_score, mean, p_val, grid_group)
# sf_object_joined_summed_by_long_samples_joined_zs$p_val <- pnorm(sf_object_joined_summed_by_long_samples_joined_zs$z_score, mean = 0, sd = 1, lower.tail = F)

#graph the ones below a certain cutoff value set by cut()
# sf_object_joined_prop_rarefy_5k_rowmeans_joined_na <- sf_object_joined_prop_rarefy_5k_rowmeans_joined[!is.na(sf_object_joined_prop_rarefy_5k_rowmeans_joined$p_val),]
# 
# vals_prop_p_0.05_means <- ggplot(data =  sf_object_joined_prop_rarefy_5k_rowmeans_joined, aes(fill = 'transparent', color = "transparent")) +
#   geom_sf(data = world, fill = 'grey80', color = 'grey70') +
#   geom_sf(aes(fill = mean, geometry = geometry), color = "transparent", lwd = 0) +
#   scale_fill_gradientn(limits = c(0.0, 0.2), 
#     colors = rev(c("#3333E7", "#737389")),
#     values = scales::rescale(c(0, 0.2))) + # Optional: set the overall limits of the scale
#   coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
#   theme_bw() 
# ggsave("p_vals_prop_p_0.05_2_means.png", vals_prop_p_0.05_means, width = 20, height = 15, units = "cm", limitsize = F)
# 
# 
# vals_prop_p_0.05 <- ggplot(data =  st_buffer(st_as_sf(sf_object_joined_prop_rarefy_5k_rowmeans_joined), 25000), aes(fill = 'transparent', color = "transparent")) +
#   geom_sf(data = world, fill = 'grey95', color = 'grey80') +
#   geom_sf(aes(fill = cut(p_val, c(-Inf, 0.05, Inf)), geometry = geometry), lwd = 0, color = "grey20") +
#   scale_fill_manual(values = c("#3333E7", "#A3A3C8"), na.value = "transparent") +
#   coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
#   theme(legend.position = 'none') +
#   theme_bw()
# 
# ggsave("outputs/p_vals_prop_p_0.05_test.png", vals_prop_p_0.05, width = 20, height = 15, units = "cm", limitsize = F)

# Split into zero and non-zero novelty for distinct coloring
sf_novel_sorted <- st_as_sf(arrange(sf_object_joined_prop_rarefy_5k_rowmeans_joined, novel_species_proportion))
sf_novel_zero    <- filter(sf_novel_sorted, novel_species_proportion == 0)
sf_novel_nonzero <- filter(sf_novel_sorted, novel_species_proportion > 0)

vals_prop_all <- ggplot() +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  # Zero novelty: distinct warm-grey so it reads as "sampled but none found"
  geom_sf(data = st_buffer(sf_novel_zero, 75000),
          fill = "grey60", color = "white", lwd = 0.2, alpha = 0.6) +
  # Non-zero novelty: blue gradient, higher values drawn on top
  geom_sf(data = st_buffer(sf_novel_nonzero, 100000),
          aes(fill = novel_species_proportion), color = "white", lwd = 0.2, alpha = 0.8) +
  scale_fill_gradientn(colors = c("#86b3c6", "#3333E7"), na.value = "transparent") +
  coord_sf(crs = "ESRI:54009") +
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggsave("outputs/vals_prop_all_split.png", vals_prop_all, width = 20, height = 15, units = "cm", limitsize = F)


vals_prop_all_overlay <- ggplot() +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  # Bottom layer: sampling density (red)
  geom_sf(data = giant_geo_table_grid_id_geometry,
          aes(fill = log(Freq)), lwd = 0, color = "transparent") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "OrRd"), na.value = "transparent") +
  ggnewscale::new_scale_fill() +
  # Zero novelty: distinct warm-grey so it reads as "sampled but none found"
  geom_sf(data = st_buffer(sf_novel_zero, 75000),
          fill = "grey20", color = "white", lwd = 0.2, alpha = 0.6) +
  # Non-zero novelty: blue gradient, higher values drawn on top
  geom_sf(data = st_buffer(sf_novel_nonzero, 75000),
          aes(fill = novel_species_proportion), color = "white", lwd = 0.2, alpha = 0.8) +
  scale_fill_gradientn(colors = c("#a8c4d4", "#5f89b0", "#3333E7"), na.value = "transparent") +
  coord_sf(crs = "ESRI:54009") +
  theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggsave("outputs/vals_prop_all_overlay_split.png", vals_prop_all_overlay, width = 20, height = 15, units = "cm", limitsize = F)


vals_prop_p_0.05 <- ggplot(data =  st_buffer(st_as_sf(sf_object_joined_prop_rarefy_5k_rowmeans_joined), 50000), aes(fill = 'transparent', color = "transparent")) +
  geom_sf(data = world, fill = 'grey95', color = 'grey80') +
  geom_sf(aes(geometry = geometry), lwd = 0, fill = "grey70", color = c("grey20")) +
  geom_sf(fill = '#3333E7', data = st_buffer(st_as_sf(filter(sf_object_joined_prop_rarefy_5k_rowmeans_joined, p_val < 0.05)), 50000), aes(geometry = geometry), lwd = 0, color = "grey20", alpha = 1) +
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme(legend.position = 'none') +
  theme_bw()

ggsave("outputs/p_vals_prop_p_0.05_test.png", vals_prop_p_0.05, width = 20, height = 15, units = "cm", limitsize = F)

# Save key object for novelty_significance_maps.R
dir.create("outputs/geo_intermediates", showWarnings = FALSE, recursive = TRUE)
saveRDS(sf_object_joined_prop_rarefy_5k_rowmeans_joined_s,
        "outputs/geo_intermediates/sf_object_joined_prop_rarefy_5k_rowmeans_joined_s.rds")
