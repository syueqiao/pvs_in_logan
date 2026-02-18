# install.packages("iNEXT")
source("R_scripts/04_geographic_analysis/geo_setup.R")
library(iNEXT)
library(ggplot2)
library(raster)

# Load from environment (if map_gen_for_pub.R was run) or from saved RDS
if (!exists("all_pvs_mapping")) {
  all_pvs_mapping <- readRDS("outputs/geo_intermediates/all_pvs_mapping.rds")
}
if (!exists("all_novels")) {
  all_novels <- readRDS("outputs/geo_intermediates/all_novels.rds")
}

# all_pvs_mapping
# 
# 
# 
# all_pvs_mapping_binned <- all_pvs_mapping_binned[complete.cases(all_pvs_mapping_binned[ , 24]),]
# 
# all_pvs_mapping_binned_groups <- as.data.frame(table(all_pvs_mapping_binned$group_id))
# hist(all_pvs_mapping_binned_groups$Freq, breaks = 50)

# write.table(all_pvs_mapping_binned$V1.x, "all_pvs_mapping_binned_id.tsv", sep = "\t", quote = F, row.names = F, col.names = F)

# cluster_number_key_for_geo is now generated from st_intersects below (after grid intersection)


all_hits_library_biosample <- read.table("files/who_puts_vlookup_man.list", sep = "\t", header = T, fill = T)
all_pvs_mapping <- left_join(all_pvs_mapping, all_hits_library_biosample, by = c("Run"))


geo_data_annotation_for_all_biosamps <- read.csv("files/geo_data_annotation_for_all_biosamps.txt", header = F)
all_pvs_mapping <- left_join(all_pvs_mapping, geo_data_annotation_for_all_biosamps, by = c("BioSample" = "V1"))



#select clean ver

all_pvs_mapping_clean <- dplyr::select(all_pvs_mapping, V1.x, Run, status, LibrarySource.y, BioSample, V2, V3.y, V4.y, V2.y)
all_pvs_mapping_clean$V2 <- factor(all_pvs_mapping_clean$V2, levels = c("lat_lon", "geographic location (latitude),geographic location (longitude)", 
                                              "geo_loc_name_sam", "geo_loc_name", 
                                              "geo_loc_name_country_calc", "geographic location (region and locality)", "region", "birth_location", "INSDC center name"))

all_pvs_mapping_clean <- all_pvs_mapping_clean[order(all_pvs_mapping_clean$V2),]

all_pvs_mapping_clean <- all_pvs_mapping_clean[!duplicated(all_pvs_mapping_clean[,c('V1.x')]),]

# world and grid_sf loaded by geo_setup.R
library(tmap)
library(terra)

# --- Visualization ---
# ggplot() +
#   geom_sf(data = grid_sf, fill = NA, color = "steelblue", size = 0.05) +
#   theme_minimal() +
#   labs(title = "50km Equal Area Grid (Terra Method)",
#        subtitle = "Correctly densified to avoid vertical line collapse")
# 
# # grid <- st_make_grid(cellsize = c(50000, 50000), crs = "ESRI:54009")
# 
# ggplot() +
#   geom_sf(data = grid_sf) +
#   coord_sf(crs = "ESRI:54009")
# # Dataframe with latlong coordinates:

# all_pvs_mapping
wkb_all = structure((all_pvs_mapping_clean$V4.y), class = "WKB")
all_pvs_mapping_clean$geometry <- st_as_sfc(wkb_all, EWKB = TRUE) %>% st_transform("ESRI:54009")

dsf <- sf::st_as_sf(all_pvs_mapping_clean, crs="ESRI:54009")

# grid_sf loaded by geo_setup.R

int <- sf::st_intersects(dsf, grid_sf)
pv_geo_grids_table <- as.data.frame(table(unlist(int)))
pv_geo_grids_table$Var1 <- as.numeric(paste0(pv_geo_grids_table$Var1))

all_pvs_mapping_clean$grid_group <- int

# Generate cluster_number_updated from grid intersection
# Each PV gets the grid cell it falls in; priority selection already applied via dedup above
all_pvs_mapping_clean$cluster_number_updated <- sapply(int, function(x) if (length(x) > 0) x[1] else NA)

cluster_number_key_for_geo <- data.frame(
  cluster_number_updated = all_pvs_mapping_clean$cluster_number_updated,
  V1.x = all_pvs_mapping_clean$V1.x
)
cluster_number_key_for_geo <- cluster_number_key_for_geo[!is.na(cluster_number_key_for_geo$cluster_number_updated), ]



write.table(cluster_number_key_for_geo, "files/cluster_number_key_for_geo.tsv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# giant_geo_table_grid_id_geometry: re-index old gpkg sample counts to new terra grid
# The old grids_sf.gpkg was indexed to a different grid; re-map via centroid intersection
old_gpkg <- st_read("files/grids_sf.gpkg")
old_gpkg <- st_as_sf(old_gpkg)
st_crs(old_gpkg) <- "ESRI:54009"

# Find which new grid cell each old cell's centroid falls in
old_centroids <- st_centroid(old_gpkg)
old_to_new <- st_intersects(old_centroids, grid_sf)
old_gpkg$new_grid_id <- sapply(old_to_new, function(x) if (length(x) > 0) x[1] else NA)

# Aggregate sample counts (Freq) by new grid cell ID
giant_geo_table_grid_id_geometry <- old_gpkg %>%
  filter(!is.na(new_grid_id)) %>%
  as.data.frame() %>%
  group_by(new_grid_id) %>%
  summarise(Freq = sum(Freq, na.rm = TRUE), .groups = "drop") %>%
  rename(Var1 = new_grid_id)

# Attach new grid geometry
giant_geo_table_grid_id_geometry <- left_join(
  data.frame(Var1 = grid_sf$id),
  giant_geo_table_grid_id_geometry,
  by = "Var1"
)
giant_geo_table_grid_id_geometry$Freq[is.na(giant_geo_table_grid_id_geometry$Freq)] <- 0
giant_geo_table_grid_id_geometry <- st_sf(giant_geo_table_grid_id_geometry, geometry = st_geometry(grid_sf))


all_pvs_mapping_clean$grid_group <- as.numeric(paste(all_pvs_mapping_clean$cluster_number_updated))

all_geo_grids_table_pvs_joined <- left_join(giant_geo_table_grid_id_geometry, st_drop_geometry(dplyr::select(all_pvs_mapping_clean, grid_group, V1.x)), by = c("Var1" = "grid_group"))
# giant_geo_table_grouping <- cbind(giant_geo_table, unlist(int_all))
# all_geo_grids_table_pvs_joined[is.na(all_geo_grids_table_pvs_joined)] <- 0
# all_geo_grids_table_pvs_joined$observation_rate <- all_geo_grids_table_pvs_joined$Freq.y/all_geo_grids_table_pvs_joined$Freq.x
# colnames(all_geo_grids_table_pvs_joined) <- c("grid_no", "sample_n", "pv_obvs_n", "pv_obs_per_sample")


# all_pvs_mapping_binned_clusters_collapse <-  all_pvs_mapping[!duplicated(all_pvs_mapping$cluster_number), ]


# all_pvs_mapping_binned_clusters_collapse <- all_pvs_mapping_binned_clusters_collapse[!is.na(all_pvs_mapping_binned_clusters_collapse$V4.y),]


dsf_novel <- sf::st_as_sf(filter(all_pvs_mapping_clean, status == "novel"), crs="ESRI:54009")

dsf_all_species <- sf::st_as_sf(all_pvs_mapping_clean, crs="ESRI:54009")

int_novels <- sf::st_intersects(dsf_novel, grid_sf)
pv_geo_grids_table_novel <- as.data.frame(table(unlist(int_novels)))
colnames(pv_geo_grids_table_novel) <- c("grid_no", "species_obvs_novel")


int_all_species <- sf::st_intersects(dsf_all_species, grid_sf)
pv_geo_grids_table_all_species <- as.data.frame(table(unlist(int_all_species)))
colnames(pv_geo_grids_table_all_species) <- c("grid_no", "species_obvs_all")

pv_geo_grids_table_novel$grid_no <- as.numeric(paste0(pv_geo_grids_table_novel$grid_no))
pv_geo_grids_table_all_species$grid_no <- as.numeric(paste0(pv_geo_grids_table_all_species$grid_no))
all_geo_grids_table_pvs_joined_add_novel <- left_join(all_geo_grids_table_pvs_joined, pv_geo_grids_table_novel, by = c("Var1" = "grid_no"))

all_geo_grids_table_pvs_joined_add_all <- left_join(all_geo_grids_table_pvs_joined_add_novel, pv_geo_grids_table_all_species, by = c("Var1" = "grid_no"))


all_geo_grids_table_pvs_joined_add_all[is.na(all_geo_grids_table_pvs_joined_add_all)] <- 0  

all_geo_grids_table_pvs_joined_add_all$novel_proportion <- all_geo_grids_table_pvs_joined_add_all$species_obvs_novel/all_geo_grids_table_pvs_joined_add_all$species_obvs_all
all_geo_grids_table_pvs_joined_add_all$species_per_sample <- all_geo_grids_table_pvs_joined_add_all$species_obvs_all/as.numeric(all_geo_grids_table_pvs_joined_add_all$Freq)


# all_pvs_mapping_binned_clusters %>%
#   ggplot() +
#   geom_polygon(aes(x = long, y = lat, group = group), data = world_map, fill = "lightgrey", color = "lightgrey") +
#   geom_point(aes(x = all_pvs_mapping_binned$lat_lon[,1], y = all_pvs_mapping_binned$lat_lon[,2], color = cluster_number), size = 2, alpha = 0.5, stroke = NA) +
#   coord_fixed() +
#   scale_color_viridis_b() + 
#   theme_classic() +
#   theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95"), legend.position = 'none') + 
#   geom_vline(xintercept = vertical_line_positions, 
#              linetype = "solid",  # Optional: change line type (e.g., "dotted", "solid")
#              color = "blue",      # Optional: change line color
#              size = 0.1)          + 
#   geom_hline(yintercept = vertical_line_positions, 
#              linetype = "solid",  # Optional: change line type (e.g., "dotted", "solid")
#              color = "red",      # Optional: change line color
#              size = 0.1)          

# length(grid)

# grid_sf already has st_as_sf + id from geo_setup.R

all_geo_grids_table_pvs_joined_add_all$grid_numeric <- as.numeric(paste0(all_geo_grids_table_pvs_joined_add_all$Var1))

sf_object_joined <- left_join(grid_sf, as.data.frame(all_geo_grids_table_pvs_joined_add_all), by = c("id" = "grid_numeric"))
sf_object_joined$sample_n_log10 <- as.numeric(log10( as.numeric(sf_object_joined$Freq)))
# class(sf_object_joined)
library(viridis)

# df_fixed$Freq <- as.numeric(df_fixed$Freq)

plot_sample_n_log10 <- ggplot(data = giant_geo_table_grid_id_geometry) +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(aes(fill = log(Freq), color = "transparent"), lwd = 0,  color = "transparent") +
  scale_fill_gradientn(colors = c( viridis(option = 'mako', n = 10, direction = -1)), na.value = "transparent") +
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

ggsave("outputs/2026.02.18.plot_sample_n_log10_new_test.png", plot_sample_n_log10, width = 20, height = 15, units = "cm", limitsize = F)

giant_geo_table_grid_id_geometry$Var1 <- as.numeric(giant_geo_table_grid_id_geometry$Var1)


all_pvs_mapping_clean_geometry <- left_join(all_pvs_mapping_clean, as.data.frame(giant_geo_table_grid_id_geometry)[, c("Var1", "Freq")], by = c("grid_group" = "Var1"))
all_pvs_mapping_clean_geometry$status_2 <- ifelse(all_pvs_mapping_clean_geometry$V1.x %in% all_novels$V1, "novel", "known")
all_pvs_mapping_clean_geometry_linears <- dplyr::select(all_pvs_mapping_clean_geometry, geometry, Freq, V1.x, status_2)

library(ggpmisc)
result <- all_pvs_mapping_clean_geometry_linears %>%
  group_by(geometry, status_2) %>%
  summarise(count = n())

result_joined <- unique(st_join(st_as_sf(result), st_as_sf(dplyr::select(all_pvs_mapping_clean_geometry_linears, geometry, Freq)), largest = T))

result_joined_plot <- dplyr::select(result_joined, -status_2) %>%
  group_by(across(-count)) %>%
  summarise(count_plot = sum(count), .groups = "drop") 



result_joined_plot_gg <- ggplot(data = result_joined_plot, aes(x=Freq, y = count_plot)) + 
  geom_point(color = "#46D452", alpha = 0.8) + 
  geom_smooth(method = 'lm', color = "grey20") + 
  stat_poly_eq(mapping = use_label("R2", "p")) +
  theme_classic() + scale_x_log10() + scale_y_log10() + theme(text = element_text(family = "Noto Sans"))

result_joined_plot_gg

ggsave("outputs/result_joined_plot.png", result_joined_plot_gg, width = 10, height = 10, units = "cm", limitsize = F,bg='transparent')

# Count novel/known PVs per grid cell (drop geometry to avoid unique() treating each point as distinct)
all_pvs_mapping_clean_geometry_novel_props <- st_drop_geometry(
  dplyr::select(all_pvs_mapping_clean_geometry, status_2, cluster_number_updated, grid_group)
)
all_pvs_mapping_clean_geometry_novel_props_table <- as.data.frame(with(all_pvs_mapping_clean_geometry_novel_props, table(grid_group, status_2)))

data_wide <- all_pvs_mapping_clean_geometry_novel_props_table %>%
  pivot_wider(
    id_cols = grid_group,          # Columns to keep as identifiers in the wide format
    names_from = status_2, # Column whose values will become new column names
    values_from = Freq,   # Column whose values will populate the new columns
    values_fill = 0     # Fill NA values (where a category is absent) with 0
  )

data_wide$novel_species_proportion <- data_wide$novel/(data_wide$known + data_wide$novel)
data_wide$grid_group <- as.numeric(paste(data_wide$grid_group))

# Rejoin grid cell geometry (one geometry per grid_group)
grid_geom_lookup <- st_sf(grid_group = grid_sf$id, geometry = st_geometry(grid_sf))
data_wide <- left_join(data_wide, grid_geom_lookup, by = "grid_group")

# Save key objects for downstream scripts
dir.create("outputs/geo_intermediates", showWarnings = FALSE, recursive = TRUE)
saveRDS(data_wide, "outputs/geo_intermediates/data_wide.rds")
saveRDS(all_pvs_mapping_clean, "outputs/geo_intermediates/all_pvs_mapping_clean.rds")
saveRDS(sf_object_joined, "outputs/geo_intermediates/sf_object_joined.rds")
saveRDS(giant_geo_table_grid_id_geometry, "outputs/geo_intermediates/giant_geo_table_grid_id_geometry.rds")
