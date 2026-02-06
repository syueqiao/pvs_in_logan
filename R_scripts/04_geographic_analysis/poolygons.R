# install.packages("iNEXT")
library(iNEXT)
library(ggplot2)
library(raster)
library(tidyverse)
library(sf)

# all_pvs_mapping
# 
# 
# 
# all_pvs_mapping_binned <- all_pvs_mapping_binned[complete.cases(all_pvs_mapping_binned[ , 24]),]
# 
# all_pvs_mapping_binned_groups <- as.data.frame(table(all_pvs_mapping_binned$group_id))
# hist(all_pvs_mapping_binned_groups$Freq, breaks = 50)

write.table(all_pvs_mapping_binned$V1.x, "all_pvs_mapping_binned_id.tsv", sep = "\t", quote = F, row.names = F, col.names = F)

cluster_number_key_for_geo <- read.table("cluster_number_key_for_geo.tsv", sep = " ", header = F)
colnames(cluster_number_key_for_geo) <- c("cluster_number_updated", "V1.x")

all_pvs_mapping <- left_join(all_pvs_mapping, cluster_number_key_for_geo, by = c("V1.x"))
all_pvs_mapping$cluster_number <- NULL


all_hits_library_biosample <- read.table("who_puts_vlookup_man.list", sep = "\t", header = T, fill = T)
all_pvs_mapping <- left_join(all_pvs_mapping, all_hits_library_biosample, by = c("Run"))


geo_data_annotation_for_all_biosamps <- read.csv("geo_data_annotation_for_all_biosamps.txt", header = F)
all_pvs_mapping <- left_join(all_pvs_mapping, geo_data_annotation_for_all_biosamps, by = c("BioSample" = "V1"))



#select clean ver

all_pvs_mapping_clean <- dplyr::select(all_pvs_mapping, V1.x, Run, status, cluster_number_updated, LibrarySource.y, BioSample, V2, V3.y, V4.y, V2.y)
all_pvs_mapping_clean$V2 <- factor(all_pvs_mapping_clean$V2, levels = c("lat_lon", "geographic location (latitude),geographic location (longitude)", 
                                              "geo_loc_name_sam", "geo_loc_name", 
                                              "geo_loc_name_country_calc", "geographic location (region and locality)", "region", "birth_location", "INSDC center name"))

all_pvs_mapping_clean <- all_pvs_mapping_clean[order(all_pvs_mapping_clean$V2),]

all_pvs_mapping_clean <- all_pvs_mapping_clean[!duplicated(all_pvs_mapping_clean[,c('V1.x')]),]

# generate grids?
library(rnaturalearthdata)
library(rnaturalearth)
library(tmap)

world <- ne_countries(scale = "medium", returnclass = "sf")
world <-  st_transform(world, "ESRI:54009")
st_is_longlat(world)

bounding_world <- st_bbox(world, crs = "ESRI:54009")

ggplot(data = world) +
  geom_sf() + coord_sf(datum=st_crs("ESRI:54009"))
library(terra)

proj_crs <- "ESRI:54009"
cell_size <- 50000  # 50km

# 2. Create the World Vector
# We use a slightly smaller extent (-179.9, 89.9) to act as a safety buffer
# against the pole singularities, similar to the sf method.
v <- vect(ext(-179.9, 179.9, -89.9, 89.9), crs="+proj=longlat")

# --- THE FIX: DENSIFY ---
# We add a point every 100km (approx 1 degree) along the lines.
# This forces the projection to draw the curve instead of a straight line between poles.
v_dense <- densify(v, interval = 1) 

# 3. Project to Mollweide
world_ellipse <- project(v_dense, proj_crs)

# 4. Create Raster Grid
# Define the raster based on the extent of the projected ellipse
r <- rast(world_ellipse, res = cell_size)

# 5. Masking
# Set values to create the grid cells
values(r) <- 1:ncell(r)

# Mask out the empty corners (ocean outside the ellipse)
r_masked <- mask(r, world_ellipse)

# 6. Convert to Polygons
# This creates the vector grid from the valid pixels
grid_terra <- as.polygons(r_masked, dissolve = FALSE)

# Convert to sf for plotting
grid_sf <- st_as_sf(grid_terra)

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



int <- sf::st_intersects(dsf, grid_sf)
pv_geo_grids_table <- as.data.frame(table(unlist(int)))
pv_geo_grids_table$Var1 <- as.numeric(paste0(pv_geo_grids_table$Var1))

all_pvs_mapping_clean$grid_group <- int
# all_pvs_mapping_binned$grid_group <- unlist(unname(int))
# all_pvs_mapping_binned$grid_group <- as.numeric(all_pvs_mapping_binned$grid_group)

# ggplot() +
#   geom_sf(data = grid) + geom_sf(data = world, fill = "NA") +
#   stat_sf_coordinates(data = all_pvs_mapping_clean, aes(geometry = geometry)) +
#   theme_classic()
  
# try to find density
#oh loooooord
  
#THIS SECTION DOES NOT NEED TO BE RUN. SIMPLY READ IN AT L191
# giant_geo_table <- head(read.table("table_1_js.csv", sep = "\t", header = F))
# head(giant_geo_table)
# readlines <- read.table("table_1_js.csv", nrows = 1000, header = F, sep = "\t")

# giant_geo_table <- readlines[!is.na(readlines$V3),]
# giant_wkb = structure((giant_geo_table$V3), class = "WKB")
# giant_wkb$geometry <- st_as_sfc(giant_wkb, EWKB = TRUE) %>% st_transform("ESRI:54009")

# giant_dsf <- sf::st_as_sf(giant_wkb$geometry, crs="ESRI:54009")


# int_all <- sf::st_intersects(giant_dsf, grid_sf)
# asdasfasf <- int_all[lapply(int_all,length)>0]
# asdasfasf <- as.data.frame(unlist(int_all))
# asdasfasf$row.id <- as.numeric(row.names(asdasfasf))


# giant_geo_table$row.id <- as.numeric(row.names(giant_geo_table))
# giant_geo_table_grid_id <- left_join(giant_geo_table, asdasfasf, by = "row.id")

# giant_geo_table_grid_id_table <- as.data.frame(table(giant_geo_table_grid_id$`unlist(int_all)`))



# giant_geo_table <- giant_geo_table[!is.na(giant_geo_table$V3),]
# giant_wkb = structure((giant_geo_table$V3), class = "WKB")
# giant_wkb$geometry <- st_as_sfc(giant_wkb, EWKB = TRUE) %>% st_transform("ESRI:54009")

# giant_dsf <- sf::st_as_sf(giant_wkb$geometry, crs="ESRI:54009")


# int_all <- sf::st_intersects(giant_dsf, grid_sf)
# asdasfasf <- int_all[lapply(int_all,length)>0]
# asdasfasf <- as.data.frame(unlist(int_all))
# asdasfasf$row.id <- as.numeric(row.names(asdasfasf))


# giant_geo_table$row.id <- as.numeric(row.names(giant_geo_table))
# giant_geo_table_grid_id <- left_join(giant_geo_table, asdasfasf, by = "row.id")

# giant_geo_table_grid_id_table <- as.data.frame(table(giant_geo_table_grid_id$`unlist(int_all)`))


# aa <- as.data.frame(grid_sf)
# aa$Var1 <- rownames(aa)

# giant_geo_table_grid_id_geometry <- left_join(giant_geo_table_grid_id_table, aa, by = c("Var1"))
# write.table(giant_geo_table_grid_id_geometry, "all_geo_grids_table_50km_by_50km_new_proj.tsv", sep = "\t", row.names = F, col.names = T)


# all_geo_grids_table <- as.data.frame(table(unlist(int_all)))
# all_geo_grids_table$Var1 <- as.numeric(paste(all_geo_grids_table$Var1))

# write.table(all_geo_grids_table, "all_geo_grids_table_50km_by_50km_addtl.tsv", sep = "\t", row.names = F, col.names = T)
# giant_geo_table_grid_id_geometry <- st_read("all_geo_grids_table_50km_by_50km_new_proj.tsv")
# st_write(obj = giant_geo_table_grid_id_geometry, dsn = "my_sf_data.gpkg", driver = "GPKG", append=FALSE)

giant_geo_table_grid_id_geometry <- st_read("my_sf_data.gpkg")


giant_geo_table_grid_id_geometry <- st_as_sf(giant_geo_table_grid_id_geometry)
st_crs(giant_geo_table_grid_id_geometry) <- "ESRI:54009"
giant_geo_table_grid_id_geometry$Var1 <- as.numeric(giant_geo_table_grid_id_geometry$Var1)


all_pvs_mapping_clean$grid_group <- as.numeric(paste(all_pvs_mapping_clean$grid_group))

all_geo_grids_table_pvs_joined <- left_join(giant_geo_table_grid_id_geometry, dplyr::select(all_pvs_mapping_clean, grid_group, geometry), by = c("Var1" = "grid_group"))
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

grid_sf <- st_as_sf(grid_sf)
grid_sf$id <- as.numeric(row.names(grid_sf))

# sf_object <- st_sf(id = 1:length(grid_sf), geometry = grid_sf$geometry)

all_geo_grids_table_pvs_joined_add_all$grid_numeric <- as.numeric(paste0(all_geo_grids_table_pvs_joined_add_all$Var1))

sf_object_joined <- left_join(grid_sf, as.data.frame(all_geo_grids_table_pvs_joined_add_all), by = c("id" = "grid_numeric"))
sf_object_joined$sample_n_log10 <- as.numeric(log10( as.numeric(sf_object_joined$Freq)))

# class(sf_object_joined)
library(viridis)

# df_fixed$Freq <- as.numeric(df_fixed$Freq)

plot_sample_n_log10 <- ggplot(data = giant_geo_table_grid_id_geometry) +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(aes(fill = log(Freq), color = "transparent", geometry = geom), lwd = 0,  color = "transparent") +
  scale_fill_gradientn(colors = c( viridis(option = 'mako', n = 10, direction = -1)), na.value = "transparent") +
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())
ggsave("2026.01.08.plot_sample_n_log10_new_test.png", plot_sample_n_log10, width = 20, height = 15, units = "cm", limitsize = F)

giant_geo_table_grid_id_geometry$Var1 <- as.numeric(giant_geo_table_grid_id_geometry$Var1)

all_pvs_mapping_clean_geometry <- left_join(all_pvs_mapping_clean, giant_geo_table_grid_id_geometry, by = c("grid_group" = "Var1"))
all_pvs_mapping_clean_geometry$status_2 <- ifelse(all_pvs_mapping_clean_geometry$V1.x %in% all_novels$V1, "novel", "known")
all_pvs_mapping_clean_geometry_linears <- dplyr::select(all_pvs_mapping_clean_geometry, geom, Freq, V1.x, status_2)

library(ggpmisc)
result <- all_pvs_mapping_clean_geometry_linears %>%
  group_by(geom, status_2) %>%
  summarise(count = n())

result_joined <- unique(st_join(st_as_sf(result), st_as_sf(dplyr::select(all_pvs_mapping_clean_geometry_linears, geom, Freq)), largest = T))

result_joined_plot <- dplyr::select(result_joined, -status_2) %>%
  group_by(across(-count)) %>%
  summarise(count_plot = sum(count), .groups = "drop") 



result_joined_plot_gg <- ggplot(data = result_joined_plot, aes(x=Freq, y = count_plot)) + 
  geom_point(color = "#46D452", alpha = 0.8) + 
  geom_smooth(method = 'lm', color = "grey20") + 
  stat_poly_eq(mapping = use_label("R2", "p")) +
  theme_classic() + scale_x_log10() + scale_y_log10() + theme(text = element_text(family = "Noto Sans"))

result_joined_plot_gg

ggsave("result_joined_plot.png", result_joined_plot_gg, width = 10, height = 10, units = "cm", limitsize = F,bg='transparent')
#5710 observations not assigned to a polygon

#something is wrong
all_pvs_mapping_clean_geometry_novel_props <-unique(dplyr::select(all_pvs_mapping_clean_geometry, geom, status_2, cluster_number_updated, grid_group))
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

data_wide <- unique(left_join(data_wide, dplyr::select(all_pvs_mapping_clean_geometry_novel_props, grid_group, geom), by = c("grid_group")))

data_wide_sig <- st_join(st_as_sf(data_wide), st_as_sf(sf_object_joined_prop_rarefy_5k_rowmeans_joined_s), largest = T)

library(sf)

borders <- st_buffer(st_as_sf(filter(data_wide_sig, p_val < 0.05)), 75000)




plot_novel_proportion <- ggplot(data = st_buffer(st_as_sf(data_wide_sig), 100000), aes(fill = novel_species_proportion)) +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(color = NA) + 
  scale_fill_gradientn(colors = c("grey70", "#3333E7"))

plot_novel_proportion <- ggplot(data = st_buffer(st_as_sf(data_wide_sig), 50000), aes(fill = novel_species_proportion)) +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(shape=16, color = NA, alpha = 0.8) + 
  scale_fill_gradientn(colors = c("grey70", "#3333E7")) +
  geom_sf (data =borders, 
          aes(geometry = geom), color = "#E7338D", fill = NA) +
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

plot_novel_proportion

ggsave("plot_novel_proportion_2.png", plot_novel_proportion, width = 30, height = 20, units = "cm", limitsize = F,bg='transparent')

# 
# plot_species_per_sample <- ggplot(data = filter(sf_object_joined[!is.na(sf_object_joined$species_per_sample),], species_obvs_all > 0), aes(color = NA)) +
#   geom_sf(data = world_map,  fill = 'grey10', color = 'grey20') +
#   geom_sf(aes(fill = log(species_per_sample), color = NA), alpha = 0.9, lwd = 0, color = NA) +
#   scale_fill_gradientn(colors = rev(terrain.colors(10)), na.value = "transparent") +
#   coord_sf(crs = 4326) +
#   theme_dark()# Explicitly set coord_sf CRS
# 
# ggsave("plot_species_per_sample.png", plot_species_per_sample, width = 20, height = 15, units = "cm", limitsize = F,bg='transparent')
# 
# sf_object_joined$novel_species_per_obvs <- sf_object_joined$species_obvs_novel/sf_object_joined$pv_obvs_n


novel_species_per_obvs <- ggplot(data = filter(sf_object_joined[!is.na(sf_object_joined$species_per_sample),], species_obvs_all > 0), aes(color = NA)) +
  geom_sf(data = world_map, fill = 'grey10', color = 'grey20') +
  geom_sf(aes(fill = novel_species_per_obvs, color = NA), alpha = 0.9, lwd = 0, color = NA) +
  scale_fill_gradientn(colors = rev(brewer.pal(9, "YlGn")), na.value = "transparent") +
  coord_sf(crs = 4326) +
  theme_dark()# Explicitly set coord_sf CRS

ggsave("novel_species_per_obvs.png", novel_species_per_obvs, width = 20, height = 15, units = "cm", limitsize = F,bg='transparent')

aa <- filter(sf_object_joined, species_obvs_all > 0)

library(viridis)

##different color scheme
ggplot(data = st_buffer(st_as_sf(data_wide_sig), 50000), aes(fill = novel_species_proportion)) +
  geom_sf(data = world, fill = 'grey40', color = 'grey60') +
  geom_sf(shape=16, color = NA, alpha = 0.9) + 
  scale_fill_gradientn(colors = plasma(n = 5)) +
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),panel.background = element_rect(fill = "grey25"),
        plot.background = element_rect(fill = "grey25"))


ggplot(data = st_buffer(st_as_sf(filter(data_wide_sig, p_val < 0.05)), 50000), aes(fill = novel_species_proportion)) +
  geom_sf(data = world, fill = 'grey50', color = 'grey60') +
  geom_sf(shape=16, color = NA, alpha = 0.9) + 
  scale_fill_gradientn(limits = c(0, 1), colors = plasma(n = 5)) +
  coord_sf(crs = "ESRI:54009") + # Explicitly set coord_sf CRSs
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),panel.background = element_rect(fill = "grey25"),
        plot.background = element_rect(fill = "grey25"))


# ##finally do for alpha diversity?
# all_pvs_mapping_binned_clusters$lon <- all_pvs_mapping_binned_clusters$lat_lon[,1]
# all_pvs_mapping_binned_clusters$lat <- all_pvs_mapping_binned_clusters$lat_lon[,2]
# all_pvs_mapping_binned_clusters_no_na <- all_pvs_mapping_binned_clusters[!is.na(all_pvs_mapping_binned_clusters$lon),]
# all_pvs_mapping_binned_clusters_geom <- sf::st_as_sf(all_pvs_mapping_binned_clusters_no_na, 
#                         coords=c("lon", "lat"), crs=4326)
# 
# 
# 
# 
# all_pvs_mapping_binned_clusters_ints <- sf::st_intersects(all_pvs_mapping_binned_clusters_geom, grid)
# all_pvs_mapping_binned_clusters_geom$grid_no <- unlist(all_pvs_mapping_binned_clusters_ints)
# 
# all_pvs_mapping_binned_clusters_geom_counts <- dcast(all_pvs_mapping_binned_clusters_geom, cluster_number ~ grid_no,
#                             value.var = "grid_no", fun.aggregate = length)
# 
# all_pvs_mapping_binned_clusters_geom_counts[,-1]
# all_pvs_mapping_binned_clusters_geom_counts_inext$AsyEst$`11572` <- NULL
# 
# all_pvs_mapping_binned_clusters_geom_counts_inext <- iNEXT(all_pvs_mapping_binned_clusters_geom_counts[,-1], q=c(0,1), datatype="abundance")
# all_pvs_mapping_binned_clusters_geom_counts_inext_rm <- filter(all_pvs_mapping_binned_clusters_geom_counts_inext$AsyEst, !Assemblage == '11572')
#                                                                
# sf_object_joined_inext <- left_join(sf_object_joined, filter(all_pvs_mapping_binned_clusters_geom_counts_inext$AsyEst, Diversity == "Simpson diversity"), by = c("grid_no" = "Assemblage"))
# 
# 
# shannon_obs <- ggplot(data = filter(sf_object_joined_inext[!is.na(sf_object_joined_inext$species_per_sample),], species_obvs_all > 0), aes(color = NA)) +
#   geom_sf(data = world_map, fill = 'grey10', color = 'grey20') +
#   geom_sf(aes(fill = log(Observed), color = NA), alpha = 0.9, lwd = 0, color = NA) +
#   scale_fill_gradientn(colors = terrain.colors(10), na.value = "transparent") +
#   coord_sf(crs = 4326) +
#   theme_dark()# Explicitly set coord_sf CRS
# 
# ggsave("simpson_obs.png", shannon_obs, width = 20, height = 15, units = "cm", limitsize = F,bg='transparent')
#try rarefact
# n <- all_pvs_mapping_binned_clusters_geom_counts$cluster_number
# 
# all_pvs_mapping_binned_clusters_geom_counts_t <- as.data.frame(t(all_pvs_mapping_binned_clusters_geom_counts[,-1]))
# colnames(all_pvs_mapping_binned_clusters_geom_counts_t) <- n
# all_pvs_mapping_binned_clusters_geom_counts_t_10 <- all_pvs_mapping_binned_clusters_geom_counts_t[-which(rowSums(all_pvs_mapping_binned_clusters_geom_counts_t[sapply(all_pvs_mapping_binned_clusters_geom_counts_t, is.numeric)]) < 100),]
# 
# 
# S <- specnumber(all_pvs_mapping_binned_clusters_geom_counts_t_10) # observed number of species
# 
# raremax <- min(rowSums(all_pvs_mapping_binned_clusters_geom_counts_t_10))
# Srare <- rarefy(all_pvs_mapping_binned_clusters_geom_counts_t_10, raremax)
# 
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# 
# rarecurve(all_pvs_mapping_binned_clusters_geom_counts_t_10, step = 20, sample = raremax, col = "blue", cex = 0.6)
# 
# # Apply the function to each dataframe in the list
# one_to_one <- as.data.frame(sf_object_joined_inext[[19]][[13004]][[1]])
# 
# polygon_sf <- one_to_one %>%
#   st_as_sf(coords = c("V1", "V2"), crs = 4326) %>%
#   summarise(do_union = TRUE) %>%
#   st_cast("POLYGON")
# 
# ggplot(world_map) +
#   geom_sf(color="grey20",fill="grey95",size=.3) +
#   theme_bw() + geom_sf(data = polygon_sf)
# 
# 
# ##spatial analysis
# sf_object_joined$polygon_midpoints <- st_centroid(sf_object_joined$geometry)
# plot(polygon_midpoints)
# sf_object_joined$polygon_long <- sapply(sf_object_joined$polygon_midpoints, "[[", 1)
# sf_object_joined$polygon_lat <- sapply(sf_object_joined$polygon_midpoints, "[[", 2)
# 
# ggplot(sf_object_joined, aes(x=polygon_long, y=log10(sample_n))) +
#   geom_col()
# 
# sf_object_joined_long <- select(sf_object_joined, polygon_long, sample_n)
# sf_object_joined_long$geometry <- NULL
# 
# sf_object_joined_long <- sf_object_joined_long %>% replace_na(list(sample_n = 0))
# sf_object_joined_long$polygon_long <- as.numeric(sf_object_joined_long$polygon_long)
# 
# sf_object_joined_summed_by_long <- aggregate(sample_n ~ as.factor(signif(polygon_long, digits = 5)), sf_object_joined_long, sum)
# colnames(sf_object_joined_summed_by_long) <- c("polygon_long", "sample_n")
# 
# sf_object_joined_summed_by_long_nz <- filter(sf_object_joined_summed_by_long, sample_n >0)
# sf_object_joined_summed_by_long_nz$sample_n_log <- log(sf_object_joined_summed_by_long_nz$sample_n)
# 
# 
# ggplot(sf_object_joined_summed_by_long, aes(x=polygon_long, y=sample_n)) +
#   geom_col()
# ##resampling
# 
# sf_object_joined_summed_by_long_samples <- sf_object_joined_summed_by_long
# colnames(sf_object_joined_summed_by_long_samples) <- c("one", "two")
# 
# set.seed(102938501)
# 
# for (i in 1:10000) {
#   # Randomly sample one value from x and assign it to y
#   sf_object_joined_summed_by_long_samples[[i]] <- sample(sf_object_joined_summed_by_long_samples$two, replace = F)
# }
# 
# population <- sample(sf_object_joined_summed_by_long_samples$two, 1000000)
# 
# 
# sf_object_joined_summed_by_long_samples_joined <- cbind(sf_object_joined_summed_by_long, sf_object_joined_summed_by_long_samples)
# 
# ggplot(sf_object_joined_summed_by_long_samples_joined, aes(x=polygon_long, y=V866)) +
#   geom_col()
# 
# data_long <- sf_object_joined_summed_by_long_samples_joined %>%
#   pivot_longer(cols = -c(polygon_long, sample_n), names_to = "variable", values_to = "value")
# 
# subset <- slice_sample(data_long, n = 10000)
# 
# ggplot(subset, aes(x = polygon_long, y = log(value), color = variable)) +
#   geom_point() +
#   labs(title = "Overlaying Multiple Y Variables",
#        y = "Value") + theme(legend.position = 'none')
# 
# 
# 
# # Calculate row means
# sf_object_joined_summed_by_long_samples_joined$row_means <- rowMeans(sf_object_joined_summed_by_long_samples_joined[, -(1:2)])
# sf_object_joined_summed_by_long_samples_joined$row_means_log <- rowMeans(log(sf_object_joined_summed_by_long_samples_joined[, -(1:2)]))
# 
# 
# # Install and load matrixStats if not already installed
# # install.packages("matrixStats")
# library(matrixStats)
# 
# # Calculate row standard deviations
# sf_object_joined_summed_by_long_samples_joined$row_sds <- rowSds(as.matrix(sf_object_joined_summed_by_long_samples_joined[, -c(1, 2, 1003)])) # rowSds works on matrices
# 
# 
# # Calculate row standard errors
# sf_object_joined_summed_by_long_samples_joined$row_ses <- sf_object_joined_summed_by_long_samples_joined$row_sds / sqrt(ncol(sf_object_joined_summed_by_long_samples_joined[, -(1:2)])) # ncol(df) gives number of columns (data points per row)
# 
# ggplot(sf_object_joined_summed_by_long_samples_joined, aes(x = polygon_long, y = log(row_means))) +
#   geom_point() + theme(legend.position = 'none')
# 
# 
# ggplot(sf_object_joined_summed_by_long_samples_joined, aes(x = polygon_long,  y = as.numeric(row_means))) +
#   geom_point() +
#   geom_errorbar(aes(ymin = row_means-row_ses, ymax = row_means+row_ses), alpha = 0.3)+
#   geom_point(aes(x = as.numeric(polygon_long),  y = as.numeric(sample_n), colour = 'hotpink')) + scale_y_log10() +
#   ylim(0, NA)
# 
# 
# 
# sf_object_joined_summed_by_long_samples_joined$z_score <- (sf_object_joined_summed_by_long_samples_joined$sample_n - sf_object_joined_summed_by_long_samples_joined$row_means) / sf_object_joined_summed_by_long_samples_joined$row_sds
# sf_object_joined_summed_by_long_samples_joined_zs <- select(sf_object_joined_summed_by_long_samples_joined, polygon_long, sample_n, z_score)
# sf_object_joined_summed_by_long_samples_joined$p_val <- pnorm(sf_object_joined_summed_by_long_samples_joined$z_score, mean = 0, sd = 1, lower.tail = F)
# sf_object_joined_summed_by_long_samples_joined_zs$p_val <- pnorm(sf_object_joined_summed_by_long_samples_joined_zs$z_score, mean = 0, sd = 1, lower.tail = F)
# 
# ggplot(sf_object_joined_summed_by_long_samples_joined, aes(x = polygon_long,  y = as.numeric(row_means))) +
#   geom_point() +
#   geom_errorbar(aes(ymin = row_means-row_sds, ymax = row_means+row_sds), alpha = 0.3)+
#   geom_point(aes(x = as.numeric(polygon_long),  y = as.numeric(sample_n), colour = cut(p_val, c(-Inf, 0.1, Inf)))) +
#   geom_hline(yintercept = c())
#   ylim(0, NA)
# 
# # ppwoot + scale_y_continuous()
# # r(aes(ymax = row_means + row_sds, ymin = row_means - row_sds),
#                 # position = "dodge")
#   
# mean <- sum(sf_object_joined_summed_by_long_samples_joined$sample_n)/length(sf_object_joined_summed_by_long_samples_joined$polygon_long)
# sd <- sd(sf_object_joined_summed_by_long_samples_joined$sample_n)
# 
# ggplot(sf_object_joined_summed_by_long_samples_joined, aes(x = polygon_long,  y = as.numeric(row_means))) +
#   geom_point() +
#   geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.3)+
#   geom_point(aes(x = as.numeric(polygon_long),  y = as.numeric(sample_n), colour = cut(p_val, c(-Inf, 0.1, Inf)))) +
#   geom_hline(yintercept = c())
# ylim(0, NA)
# 
# long_values <- select(sf_object_joined_summed_by_long_samples_joined, z_score,p_val, polygon_long)
# 
# #####AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# 
# sf_object_joined_lat <- select(sf_object_joined, polygon_lat, sample_n)
# sf_object_joined_lat$geometry <- NULL
# 
# sf_object_joined_lat <- sf_object_joined_lat %>% replace_na(list(sample_n = 0))
# sf_object_joined_lat$polygon_lat <- as.numeric(sf_object_joined_lat$polygon_lat)
# 
# sf_object_joined_summed_by_lat <- aggregate(sample_n ~ as.factor(signif(polygon_lat, digits = 2)), sf_object_joined_lat, sum)
# colnames(sf_object_joined_summed_by_lat) <- c("polygon_lat", "sample_n")
# 
# sf_object_joined_summed_by_lat_nz <- filter(sf_object_joined_summed_by_lat, sample_n >0)
# # sf_object_joined_summed_by_lat_nz$sample_n_log <- log(sf_object_joined_summed_by_lat_nz$sample_n)
# 
# 
# # ggplot(sf_object_joined_summed_by_lat, aes(x=polygon_lat, y=sample_n)) +
#   geomC_col()
# ##resampling
# 
# sf_object_joined_summed_by_lat_samples <- sf_object_joined_summed_by_lat
# # colnames(sf_object_joined_summed_by_lat_samples) <- c("one", "two")
# 
# set.seed(102938501)
# 
# # for (i in 1:10000) {
# #   # Randomly sample one value from x and assign it to y
# #   sf_object_joined_summed_by_lat_samples[[i]] <- sample(sf_object_joined_summed_by_lat_samples$two, replace = F)
# # }
# 
# # population <- sample(sf_object_joined_summed_by_lat_samples$two, 1000000)
# 
# 
# # sf_object_joined_summed_by_lat_samples_joined <- cbind(sf_object_joined_summed_by_lat, sf_object_joined_summed_by_lat_samples)
# 
# # ggplot(sf_object_joined_summed_by_lat_samples_joined, aes(x=polygon_lat, y=V866)) +
#   # geom_col()
# 
# # data_lat <- sf_object_joined_summed_by_lat_samples_joined %>%
# #   pivot_later(cols = -c(polygon_lat, sample_n), names_to = "variable", values_to = "value")
# # 
# # subset <- slice_sample(data_lat, n = 10000)
# # 
# # ggplot(subset, aes(x = polygon_lat, y = log(value), color = variable)) +
# #   geom_point() +
# #   labs(title = "Overlaying Multiple Y Variables",
# #        y = "Value") + theme(legend.position = 'none')
# 
# 
# 
# # # Calculate row means
# # sf_object_joined_summed_by_lat_samples_joined$row_means <- rowMeans(sf_object_joined_summed_by_lat_samples_joined[, -(1:2)])
# # sf_object_joined_summed_by_lat_samples_joined$row_means_log <- rowMeans(log(sf_object_joined_summed_by_lat_samples_joined[, -(1:2)]))
# # 
# # 
# # # Install and load matrixStats if not already installed
# # # install.packages("matrixStats")
# # library(matrixStats)
# # 
# # # Calculate row standard deviations
# # sf_object_joined_summed_by_lat_samples_joined$row_sds <- rowSds(as.matrix(sf_object_joined_summed_by_lat_samples_joined[, -c(1, 2, 1003)])) # rowSds works on matrices
# # 
# # 
# # # Calculate row standard errors
# # sf_object_joined_summed_by_lat_samples_joined$row_ses <- sf_object_joined_summed_by_lat_samples_joined$row_sds / sqrt(ncol(sf_object_joined_summed_by_lat_samples_joined[, -(1:2)])) # ncol(df) gives number of columns (data points per row)
# # 
# # ggplot(sf_object_joined_summed_by_lat_samples_joined, aes(x = polygon_lat, y = log(row_means))) +
# #   geom_point() + theme(legend.position = 'none')
# 
# mean_lat <- sum(sf_object_joined_summed_by_lat_samples$sample_n)/length(sf_object_joined_summed_by_lat_samples$polygon_lat)
# sd_lat <- sd(sf_object_joined_summed_by_lat_samples$sample_n)
# 
# 
# 
# sf_object_joined_summed_by_lat_samples$z_score <- (sf_object_joined_summed_by_lat_samples$sample_n - mean_lat) / sd_lat
# sf_object_joined_summed_by_lat_samples_zs <- select(sf_object_joined_summed_by_lat_samples, polygon_lat, sample_n, z_score)
# sf_object_joined_summed_by_lat_samples$p_val <- pnorm(sf_object_joined_summed_by_lat_samples$z_score, mean = 0, sd = 1, lower.tail = F)
# sf_object_joined_summed_by_lat_samples_zs$p_val <- pnorm(sf_object_joined_summed_by_lat_samples_zs$z_score, mean = 0, sd = 1, lower.tail = F)
# 
# ggplot(sf_object_joined_summed_by_lat_samples, aes(x = polygon_lat,  y = mean_lat)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = mean_lat-sd_lat, ymax = mean_lat+sd_lat), alpha = 0.3)+
#   geom_point(aes(x = as.numeric(polygon_lat),  y = as.numeric(sample_n), colour = cut(p_val, c(-Inf, 0.1, Inf)))) +
#   geom_hline(yintercept = c())
# ylim(0, NA)
# 
# # ppwoot + scale_y_continuous()
# # r(aes(ymax = row_means + row_sds, ymin = row_means - row_sds),
# # position = "dodge")
# 
# 
# ggplot(sf_object_joined_summed_by_lat_samples, aes(x = polygon_lat,  y = mean_lat)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = mean_lat-sd_lat, ymax = mean_lat+sd_lat), alpha = 0.3)+
#   geom_point(aes(x = as.numeric(polygon_lat),  y = as.numeric(sample_n), colour = cut(p_val, c(-Inf, 0.1, Inf)))) +
#   geom_hline(yintercept = c())
# 
# lat_sig <- filter(sf_object_joined_summed_by_lat_samples_zs, p_val < 0.1)
# long_sig <- filter(sf_object_joined_summed_by_long_samples_joined_zs, p_val < 0.1)
# 
# ggplot(world_map) +
#   geom_sf(color="grey20",fill="grey95",size=.3) +
#   geom_hline(yintercept = as.numeric(paste0(lat_sig$polygon_lat))) +
#   geom_vline(xintercept = as.numeric(paste0(long_sig$polygon_long))) +
#   theme_bw()
# 
# lat_values <- select(sf_object_joined_summed_by_lat_samples, z_score,p_val, polygon_lat)
# 
# ####apply to the PVs?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########?########
# 
# sf_object_joined_novel_prop <- select(sf_object_joined, polygon_long, species_obvs_novel, species_obvs_all)
# sf_object_joined_novel_prop$geometry <- NULL
# 
# # sf_object_joined_long <- sf_object_joined_novel_prop %>% replace_na(list(sample_n = 0))
# sf_object_joined_novel_prop$polygon_long <- as.numeric(sf_object_joined_novel_prop$polygon_long)
# 
# sf_object_joined_novel_prop <- aggregate(cbind(species_obvs_novel,species_obvs_all) ~ as.factor(signif(polygon_long, digits = 5)), sf_object_joined_novel_prop, sum)
# colnames(sf_object_joined_novel_prop) <- c("polygon_long", "species_obvs_novel", "species_obvs_all")
# 
# ggplot(sf_object_joined_novel_prop, aes(x=polygon_long, y=species_obvs_novel/species_obvs_all)) +
#   geom_col()
# ##resampling
# 
# sf_object_joined_novel_prop_samples <- sf_object_joined_novel_prop
# # colnames(sf_object_joined_summed_by_long_samples) <- c("one", "two")
# 
# set.seed(102938501)
# 
# sf_object_joined_novel_prop_samples_nz <- filter(sf_object_joined_novel_prop_samples, species_obvs_all > 0)
# 
# mean_prop <- sum(sf_object_joined_novel_prop_samples_nz$species_obvs_novel/sf_object_joined_novel_prop_samples_nz$species_obvs_all)/length(sf_object_joined_novel_prop_samples$polygon_long)
# sd_prop <- sd(sf_object_joined_novel_prop_samples_nz$species_obvs_novel/sf_object_joined_novel_prop_samples_nz$species_obvs_all)
# 
# 
# 
# sf_object_joined_novel_prop_samples_nz$z_score <- ((sf_object_joined_novel_prop_samples_nz$species_obvs_novel/sf_object_joined_novel_prop_samples_nz$species_obvs_all) - mean_prop) / sd_prop
# # sf_object_joined_summed_by_long_samples_zs <- select(sf_object_joined_summed_by_long_samples, polygon_long, sample_n, z_score)
# sf_object_joined_novel_prop_samples_nz$p_val <- pnorm(sf_object_joined_novel_prop_samples_nz$z_score, mean = 0, sd = 1, lower.tail = F)
# # sf_object_joined_summed_by_long_samples_zs$p_val <- pnorm(sf_object_joined_summed_by_long_samples_zs$z_score, mean = 0, sd = 1, lower.tail = F)
# 
# ggplot(sf_object_joined_novel_prop_samples, aes(x = as.numeric(paste0(polygon_long)),  y = mean_prop)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = mean_prop-sd_prop, ymax = mean_prop+sd_prop), alpha = 0.3)+
#   geom_point(data = sf_object_joined_novel_prop_samples_nz, aes(x =  as.numeric(paste0(polygon_long)),  y = as.numeric(species_obvs_novel/species_obvs_all), colour = cut(p_val, c(-Inf, 0.05, Inf)))) +
#   geom_hline(yintercept = c())
# # ppwoot + scale_y_continuous()
# # r(aes(ymax = row_means + row_sds, ymin = row_means - row_sds),
# # position = "dodge")
# 
# 
# ggplot(sf_object_joined_summed_by_long_samples, aes(x = polygon_long,  y = mean_long)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = mean_long-sd_long, ymax = mean_long+sd_long), alpha = 0.3)+
#   geom_point(aes(x = as.numeric(polygon_long),  y = as.numeric(sample_n), colour = cut(p_val, c(-Inf, 0.05, Inf)))) +
#   geom_hline(yintercept = c())
# 
# long_sig <- filter(sf_object_joined_summed_by_long_samples_zs, p_val < 0.1)
# long_sig <- filter(sf_object_joined_summed_by_long_samples_joined_zs, p_val < 0.1)
# 
# ggplot(world_map) +
#   geom_sf(color="grey20",fill="grey95",size=.3) +
#   geom_hline(yintercept = as.numeric(paste0(long_sig$polygon_long))) +
#   geom_vline(xintercept = as.numeric(paste0(long_sig$polygon_long))) +
#   theme_bw()
# 
# ##############################################################################################################################################################################################################################################
# 
# sf_object_joined_novel_prop <- select(sf_object_joined, polygon_lat, species_obvs_novel, species_obvs_all)
# sf_object_joined_novel_prop$geometry <- NULL
# 
# # sf_object_joined_lat <- sf_object_joined_novel_prop %>% replace_na(list(sample_n = 0))
# sf_object_joined_novel_prop$polygon_lat <- as.numeric(sf_object_joined_novel_prop$polygon_lat)
# 
# sf_object_joined_novel_prop <- aggregate(cbind(species_obvs_novel,species_obvs_all) ~ as.factor(signif(polygon_lat, digits = 6)), sf_object_joined_novel_prop, sum)
# colnames(sf_object_joined_novel_prop) <- c("polygon_lat", "species_obvs_novel", "species_obvs_all")
# 
# # ggplot(sf_object_joined_novel_prop, aes(x=polygon_lat, y=species_obvs_novel/species_obvs_all)) +
#   # geom_col()
# ##resampling
# 
# sf_object_joined_novel_prop_samples <- sf_object_joined_novel_prop
# # colnames(sf_object_joined_summed_by_lat_samples) <- c("one", "two")
# 
# set.seed(102938501)
# 
# sf_object_joined_novel_prop_samples_nz <- filter(sf_object_joined_novel_prop_samples, species_obvs_all > 0)
# 
# mean_prop <- sum(sf_object_joined_novel_prop_samples_nz$species_obvs_novel/sf_object_joined_novel_prop_samples_nz$species_obvs_all)/length(sf_object_joined_novel_prop_samples$polygon_lat)
# sd_prop <- sd(sf_object_joined_novel_prop_samples_nz$species_obvs_novel/sf_object_joined_novel_prop_samples_nz$species_obvs_all)
# 
# 
# 
# sf_object_joined_novel_prop_samples_nz$z_score <- ((sf_object_joined_novel_prop_samples_nz$species_obvs_novel/sf_object_joined_novel_prop_samples_nz$species_obvs_all) - mean_prop) / sd_prop
# # sf_object_joined_summed_by_lat_samples_zs <- select(sf_object_joined_summed_by_lat_samples, polygon_lat, sample_n, z_score)
# sf_object_joined_novel_prop_samples_nz$p_val <- pnorm(sf_object_joined_novel_prop_samples_nz$z_score, mean = 0, sd = 1, lower.tail = F)
# # sf_object_joined_summed_by_lat_samples_zs$p_val <- pnorm(sf_object_joined_summed_by_lat_samples_zs$z_score, mean = 0, sd = 1, lower.tail = F)
# 
# ggplot(sf_object_joined_novel_prop_samples, aes(x = as.numeric(paste0(polygon_lat)),  y = mean_prop)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = mean_prop-sd_prop, ymax = mean_prop+sd_prop), alpha = 0.3)+
#   geom_point(data = sf_object_joined_novel_prop_samples_nz, aes(x =  as.numeric(paste0(polygon_lat)),  y = as.numeric(species_obvs_novel/species_obvs_all), colour = cut(p_val, c(-Inf, 0.05, Inf)))) +
#   geom_hline(yintercept = c())
# 
# long_value_prop <- filter(sf_object_joined_novel_prop_samples_nz, p_val < 0.05)
# lat_value_prop <- filter(sf_object_joined_novel_prop_samples_nz, p_val < 0.05)
# 
# # target_latitudes <- lat_values$polygon_lat
# 
# # Create a line geometry spanning the longitude range of your map at the target latitude
# # You might need to adjust the xmin and xmax based on your specific map's extent
# latitude_line <- st_linestring(matrix(c(-180, target_latitude, 180, target_latitude), ncol = 2, byrow = TRUE)) %>%
#   st_sfc(crs = 4326) # Assign WGS84 CRS
# 
# custom_graticule <- st_graticule(
#   lat = as.numeric(paste0(filter(lat_values, p_val < 0.05)$polygon_lat)),
#   lon = as.numeric(paste0(filter(long_values, p_val < 0.05)$polygon_long)))
# 
# ggplot(data = filter(sf_object_joined, sample_n_log10 > 2)) +
#   geom_sf(data = world_map, fill = 'grey30', color = 'grey20') +
#   geom_sf(aes(fill = sample_n_log10, color = "transparent",  alpha = 0.8), lwd = 0, color = "transparent") +
#   scale_fill_gradientn(colors = (viridis(option = 'plasma', n = 10)), na.value = "transparent") +
#   geom_sf(data = custom_graticule, color = "darkblue", size = 0.5) + # Custom graticule
#   coord_sf(crs = 4326)+ theme_dark()
# 
# custom_graticule_prop <- st_graticule(
#   lat = as.numeric(paste0(filter(lat_value_prop, p_val < 0.05)$polygon_lat)),
#   lon = as.numeric(paste0(filter(long_value_prop, p_val < 0.05)$polygon_long)))
# 
# ggplot(data = sf_object_joined[!is.na(sf_object_joined$novel_proportion),]) +
#   ge                                                                                                                                                                                                                                                                                                                             om_sf(data = world_map, fill = 'grey30', color = 'grey20') +
#   geom_sf(aes(fill = novel_proportion, color = "transparent"), lwd = 0, color = "transparent") +
#   scale_fill_gradientn(colors = heat.colors(10), na.value = "transparent") +
#   geom_sf(data = custom_graticule_prop, color = "darkblue", size = 0.5) + # Custom graticule
#   coord_sf(crs = 4326)+ theme_dark()




ggplot(data = sf_object_joined, aes(x = sample_n, y = pv_obvs_n)) +
  geom_point() + geom_smooth(method = 'lm') + stat_poly_eq(mapping = use_label("R2", "p")) + scale_y_log10()+scale_x_log10() + theme_classic()
