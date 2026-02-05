
library(rnaturalearthdata)
library(rnaturalearth)

###ecology?
giant_geo_table_grid_id_geometry <- st_read("my_sf_data.gpkg")

data_wide_joining <- st_as_sf(data_wide)  %>% st_transform("ESRI:54009")
ocean_polygons <- st_as_sf(ne_download(scale = 50, type = 'ocean', category = 'physical', returnclass = "sf")) %>% st_transform("ESRI:54009") %>%
  rename(featurecla_ocean = featurecla)
fw_polygons <- st_as_sf(ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")) %>% st_transform("ESRI:54009")%>%
  rename(featurecla_fw = featurecla)
fw_polygons_2 <- st_as_sf(ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = "sf")) %>% st_transform("ESRI:54009")%>%
  rename(featurecla_fw_2 = featurecla)

water <- bind_rows(dplyr::select(ocean_polygons, geometry, featurecla_ocean), dplyr::select(fw_polygons, geometry, featurecla_fw), dplyr::select(fw_polygons_2, geometry, featurecla_fw_2))

# my_shapefile <- read_sf("Global_200_Terrestrial.shp")
my_shapefile_2 <- read_sf("Ecoregions2017.shp")
# my_shapefile_fw <- st_as_sf(read_sf("Global_200_Marine.shp")) %>% st_transform("ESRI:54009")
# my_shapefile_fw_2 <- st_as_sf(read_sf("Global_200_Marine.shp")) %>% st_transform("ESRI:54009")

# my_shapefile_aqua_2 <- read_sf("WWF_Global_200_Ecoregions.shp")

# my_shapefile <- st_as_sf(my_shapefile) %>% st_transform("ESRI:54009")
my_shapefile_2 <- st_as_sf(my_shapefile_2) %>% st_transform("ESRI:54009")
my_shapefile_2$BIOME_NAME <- gsub("N/A", NA, my_shapefile_2$BIOME_NAME)
my_shapefile_2_clean <- dplyr::select(my_shapefile_2, geometry, BIOME_NAME)

all_cats <- bind_rows(water, my_shapefile_2_clean)
all_cats$area <- st_area(all_cats)

all_cats$area_biome <- coalesce(all_cats$BIOME_NAME, all_cats$featurecla_ocean, all_cats$featurecla_fw_2, all_cats$featurecla_fw)
all_cats$area <- parse_number(as.character(all_cats$area))

all_cats_area <- all_cats %>%
  group_by(area_biome) %>%
  summarise(total_area = sum(area), .groups = 'drop')

all_cats_area$area_biome <- replace_na(all_cats_area$area_biome, "Terrestrial Other")
# my_shapefile_aqua <- st_as_sf(my_shapefile_aqua)  %>% st_transform("ESRI:54009")
# my_shapefile_aqua_2 <- st_as_sf(my_shapefile_aqua_2)  %>% st_transform("ESRI:54009")

# biomes_joined <- st_join(data_wide_joining, all_cats, join = st_intersects, largest = TRUE)

sf_use_s2(FALSE)

# Create centroids of your grid cells
grid_centroids <- st_centroid(grid_sf)

# Find which biome each centroid falls into
grid_sf$biome_id <- st_intersects(grid_centroids, all_cats) %>%
  sapply(function(x) if(length(x) > 0) x[1] else NA)

# Add biome name (adjust column name as needed)
grid_sf$biome_name <- all_cats$area_biome[grid_sf$biome_id]

grid_sf$biome_name[is.na(grid_sf$biome_name)] <- "Terrestrial Other"


# Calculate total area per biome
biome_summary <- grid_sf %>%
  st_drop_geometry() %>%
  group_by(biome_name) %>%
  summarise(
    n_cells = n(),
    area_km2 = n() * 2500  # 50x50 = 2500 km² per cell
  ) %>%
  arrange(desc(area_km2))

print(biome_summary)

biome_summary$biome_name[is.na(biome_summary$biome_name)] <- "Terrestrial Other"

grid_sampling <- grid_sf

data_wide_joining_centroids <- st_centroid(data_wide_joining)


data_wide_joining$biome_id <- st_intersects(data_wide_joining_centroids, all_cats) %>%
  sapply(function(x) if(length(x) > 0) x[1] else NA)

data_wide_joining$biome_name <- all_cats$area_biome[data_wide_joining$biome_id]

data_wide_joining$biome_name[is.na(data_wide_joining$biome_name)] <- "Terrestrial Other"

biomes_joined <- data_wide_joining

# biomes_joined$guess_1 <- coalesce(biomes_joined$BIOME_NAME, biomes_joined$featurecla_ocean, biomes_joined$featurecla_fw_2, biomes_joined$featurecla_fw)
# biomes_joined_aqua_guess <- st_join(biomes_joined_aqua, my_shapefile, join = st_nearest_feature)
# biomes_joined_aqua_guess$guess_2 <- coalesce(biomes_joined_aqua$guess_1,as.character(biomes_joined_aqua_guess$G200_BIOME.y))

# 
# biomes <- c(1:14)
# biomes <- as.data.frame(c(1:14))
# biomes$desc_2 <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
# colnames(biomes) <- c("num", "desc")
# biomes$num <- as.character(biomes$num)
# 
# # library(data.table)
# 
# biomes_joined_aqua_guess <- left_join(biomes_joined_aqua, biomes, by = c("guess_1" = "num"))

# biomes_joined_aqua_guess$plot_cat <- coalesce(as.character(biomes_joined_aqua_guess$desc),biomes_joined_aqua_guess$guess_1)

summed_biomes <- dplyr::select(biomes_joined, biome_name, known, novel) %>%
  group_by(biome_name) %>% summarise(across(1:2, sum))

summed_biomes$proportions <- summed_biomes$novel/summed_biomes$known

plot(summed_biomes$proportions)

library(data.table)

long <- melt(setDT(summed_biomes), id.vars = c("biome_name")) %>% filter(!variable == "geom") %>% filter(!variable == "proportions")
long$variable <- relevel(long$variable, 'novel')


aa <- ggplot(data=long, aes(x=fct_reorder(biome_name, as.numeric(value)), y=as.numeric(value), fill=variable)) +
  geom_bar(stat="identity") + coord_flip()

aa


# marine_list <- my_shapefile_aqua$G200_REGIO
# fw_list <- my_shapefile_aqua_2$G200_REGIO
# 
# biomes_joined_aqua_guess_nz <- biomes_joined_aqua_guess
# biomes_joined_aqua_guess_nz <- biomes_joined_aqua_guess_nz %>%
#   mutate(
#     color = ifelse(plot_cat %in% marine_list, "#489690", ifelse(plot_cat %in% fw_list, "#6E96D0", ifelse(plot_cat == "Ocean", "#09BEE9", NA)))
#   )

biomes_list <- c("Tropical & Subtropical Moist Broadleaf Forests", 
                 "Tropical & Subtropical Dry Broadleaf Forests",
                 "Tropical & Subtropical Coniferous Forests",
                 "Temperate Broadleaf & Mixed Forests",
                 "Temperate Conifer Forests",
                 "Boreal Forests/Taiga",
                 "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                 "Temperate Grasslands, Savannas & Shrublands",
                 "Flooded Grasslands & Savannas",
                 "Montane Grasslands & Shrublands",
                 "Tundra",
                 "Mediterranean Forests, Woodlands & Scrub",
                 "Deserts & Xeric Shrublands",
                 "Mangroves", "Terrestrial Other", "Ocean", "Lake", "Alkaline Lake", "Reservoir")
biomes_color <- c('#008346','#9DCC00', '#C4B72E',
                  '#015C31','#20BC10','#FFA8BB'
                  ,'#FAD505','#8F7C00','#5AD488',
                  '#993E01','#C20088', '#81A527'
                  ,'#FFA405','#FFCC99', "#624E38", "#09BEE9","#6E96D0","#284A7C", "#8CD4D4")
biomes_color_key <- as.data.frame(cbind(biomes_list, biomes_color))

summed_biomes_colored <- left_join(summed_biomes, biomes_color_key, by = c("biome_name" = "biomes_list"))
summed_biomes_colored$proportions <- summed_biomes_colored$novel/summed_biomes_colored$known

summed_biomes_colored_props <- ggplot(data=summed_biomes_colored, 
                                          aes(x=fct_reorder(biome_name, as.numeric(proportions)), 
                                              y=as.numeric(proportions), fill=biomes_color)) +
  geom_text(aes(label = signif(proportions, 2)), hjust = -0.1, size = 2, fontface = 'bold') +
  geom_bar(stat="identity") + 
  scale_fill_identity(guide='legend', labels=summed_biomes_colored$biome_name, breaks=summed_biomes_colored$biomes_color) +
  coord_flip() + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

ggsave("2026.01.13summed_biomes_colored_props.png", summed_biomes_colored_props, width = 30, height = 20, units = "cm", limitsize = F ,bg='transparent')


# biomes_joined_aqua_guess_nz_color_2 <- biomes_joined_aqua_guess_nz_color_2 %>%
#   mutate(
#     plot_cat = ifelse(plot_cat %in% marine_list, "Marine", ifelse(plot_cat %in% fw_list, "Freshwater", plot_cat))
#   )

summed_biomes_colored_long <- summed_biomes_colored %>%
  pivot_longer(
    cols = c(known, novel),        # Columns to pivot
    names_to = "variable", # Name of the new key column
    values_to = "value"    # Name of the new value column
  )

summed_biomes_colored_long <- left_join(summed_biomes_colored_long, dplyr::select(summed_biomes_colored, biome_name), by = c("biome_name"))

summed_biomes_colored_long_sum <- summed_biomes_colored_long %>%
  group_by(biome_name, proportions) %>%
  summarize(sumz = sum(value), .groups = 'drop')

summed_biomes_colored_long <- unique(left_join(summed_biomes_colored_long, summed_biomes_colored_long_sum, by = c("biome_name")))

summed_biomes_colored_long$proportions.x[seq_len(nrow(summed_biomes_colored_long)) %% 2 == 0] <- NA

summed_biomes_colored_long_plot <- ggplot(data=summed_biomes_colored_long, 
       aes(x=fct_reorder(biome_name, as.numeric(value)), y=as.numeric(value), fill=biomes_color)) +
  geom_bar(stat="identity") + 
  geom_text(aes(y = sumz, x = fct_reorder(biome_name, as.numeric(sumz)), label = signif(proportions.x, 2)), size = 2, hjust = -0.3, fontface = 'bold') +
  geom_bar(data=filter(summed_biomes_colored_long, variable == 'known'), stat="identity", fill = "grey90", alpha = 0.8) + 
  scale_fill_identity(guide='legend', labels=summed_biomes_colored_long$biome_name, breaks=summed_biomes_colored_long$biomes_color) +
theme_bw() + coord_flip() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

summed_biomes_colored_long_plot

ggsave("2026.01.13.summed_biomes_colored_long_plot.png", summed_biomes_colored_long_plot, width = 40, height = 20, units = "cm", limitsize = F,bg='transparent')


# summed_biomes_colored_long_plot_zoom <- summed_biomes_colored_long_plot + 
#   coord_flip(ylim = c(0, 10)) +
#   labs(title = "Zoomed View (0-10)")  +
#   theme_bw() +theme(text = element_text(family = "Noto Sans"))  +
#   theme(axis.title = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.border = element_blank()) 
# 
# ggsave("summed_biomes_colored_long_plot_zoom.png", summed_biomes_colored_long_plot_zoom, width = 40, height = 20, units = "cm", limitsize = F,bg='transparent')
# 


psh__pv <- ggplot(data = st_buffer(st_as_sf(summed_biomes_colored_long), 0), aes(color = NA, fill = biomes_color)) + 
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(color = "transparent", alpha = 0.8, aes(geometry = geom)) +
  scale_fill_identity(guide='legend', labels=summed_biomes_colored_long$biome_name, breaks=summed_biomes_colored_long$biomes_color) +
  coord_sf(crs = "ESRI:54009") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right')

ggsave("2026.01.13.psh__pv.png", psh__pv, width = 40, height = 30, units = "cm", limitsize = F,bg='transparent')


aaa <- ggplot(data=summed_biomes_colored_long, aes(x=fct_reorder(biome_name, as.numeric(sumz)), y=as.numeric(sumz))) +
  geom_bar(stat="identity") + coord_flip()

aaa

summed_biomes_sampling_color_pv <- left_join(summed_biomes_colored_long, biomes_color_key, by = c("biome_name" = "biomes_list"))

summed_biomes_colored_long_plot_sampling_pv <- ggplot(data=filter(summed_biomes_sampling_color_pv, !is.na(proportions.x)) , aes(x=fct_reorder(biome_name, as.numeric(sumz)), y=as.numeric(sumz), fill =biomes_color.x)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label = sumz), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide='legend', labels=summed_biomes_sampling_color_pv$biome_name, breaks=summed_biomes_sampling_color_pv$biomes_color.x) +
  coord_flip() + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

summed_biomes_colored_long_plot_sampling_pv

ggsave("2026.01.13.summed_biomes_colored_long_plot_sampling_pv.png", summed_biomes_colored_long_plot_sampling_pv, width = 40, height = 20, units = "cm", limitsize = F,bg='transparent')


# biomes_joined_aqua_guess_nz_color_2$plot_cat <- gsub("Salween River|Gulf of Alaska Coastal Rivers", "Terrestrial Other", biomes_joined_aqua_guess_nz_color_2$plot_cat)

#normalize by sampling?##########?##########?##########?##########?##########?##########?##########?##########?##########?##########?##########

df_fixed_joining <- st_as_sf(giant_geo_table_grid_id_geometry)  %>% st_transform("ESRI:54009")

df_fixed_joining <- filter(df_fixed_joining, Freq > 0)

df_fixed_joining_centroids <- st_centroid(df_fixed_joining)

df_fixed_joining$biome_id <- st_intersects(df_fixed_joining_centroids, all_cats) %>%
  sapply(function(x) if(length(x) > 0) x[1] else NA)

# Add biome name (adjust column name as needed)
df_fixed_joining$biome_name <- all_cats$area_biome[df_fixed_joining$biome_id]

df_fixed_joining$biome_name[is.na(df_fixed_joining$biome_name)] <- "Terrestrial Other"

biomes_joined_sampling <- df_fixed_joining

st_write(obj = biomes_joined_sampling, dsn = "2026.01.13.biomes_joined_sampling.gpkg", driver = "GPKG", append=FALSE)

# biomes_joined_sampling <- st_read("biomes_joined_sampling.gpkg")

# biomes_joined_sampling$guess_1 <- coalesce(biomes_joined_sampling$BIOME_NAME, biomes_joined_sampling$featurecla_ocean, biomes_joined_sampling$featurecla_fw_2, biomes_joined_sampling$featurecla_fw)

# biomes_joined_sampling$guess_1[is.na(biomes_joined_sampling$guess_1)] <- "Terrestrial Other"

biomes_joined_sampling_nz <- filter(biomes_joined_sampling, Freq >0)

biomes_list <- c("Tropical & Subtropical Moist Broadleaf Forests", 
                 "Tropical & Subtropical Dry Broadleaf Forests",
                 "Tropical & Subtropical Coniferous Forests",
                 "Temperate Broadleaf & Mixed Forests",
                 "Temperate Conifer Forests",
                 "Boreal Forests/Taiga",
                 "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                 "Temperate Grasslands, Savannas & Shrublands",
                 "Flooded Grasslands & Savannas",
                 "Montane Grasslands & Shrublands",
                 "Tundra",
                 "Mediterranean Forests, Woodlands & Scrub",
                 "Deserts & Xeric Shrublands",
                 "Mangroves", "Terrestrial Other", "Ocean", "Lake", "Alkaline Lake", "Reservoir")
biomes_color <- c('#008346','#9DCC00', '#C4B72E',
                  '#015C31','#20BC10','#FFA8BB'
                  ,'#FAD505','#8F7C00','#5AD488',
                  '#993E01','#C20088', '#81A527'
                  ,'#FFA405','#FFCC99', "#624E38", "#09BEE9","#6E96D0","#284A7C", "#8CD4D4")

biomes_color_key <- as.data.frame(cbind(biomes_list, biomes_color))

# biomes_joined_sampling_nz <- biomes_joined_sampling_nz %>% rename(plot_cat = guess_1)

biomes_joined_sampling_nz_color <- left_join(biomes_joined_sampling_nz, biomes_color_key, by = c("biome_name" = "biomes_list"))

# biomes_joined_aqua_sampling_guess_nz_color_2 <- biomes_joined_aqua_sampling_guess_nz_color_2 %>%
#   mutate(
#     plot_cat = ifelse(plot_cat %in% marine_list, "Marine", ifelse(plot_cat %in% fw_list, "Freshwater", plot_cat))
#   )

psh <- ggplot(data = biomes_joined_sampling_nz_color, aes(color = NA, fill = biomes_color)) + 
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(color = "transparent", aes(geometry = geom)) +
  scale_fill_identity(guide='legend', labels=biomes_joined_sampling_nz_color$biome_name, breaks=biomes_joined_sampling_nz_color$biomes_color) +
  coord_sf(crs = "ESRI:54009") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right')

psh

ggsave("2026.01.13.psh.png", psh, width = 40, height = 30, units = "cm", limitsize = F,bg='transparent')


summed_biomes_sampling <- dplyr::select(biomes_joined_sampling_nz_color, biome_name, Freq, biomes_color) %>%
  group_by(biome_name) %>%
  summarize(freq_count = sum(Freq, na.rm = TRUE))

aa <- ggplot(data=summed_biomes_sampling, aes(x=fct_reorder(biome_name, as.numeric(freq_count)), y=as.numeric(freq_count))) +
  geom_bar(stat="identity") + coord_flip()

aa

summed_biomes_sampling_color <- left_join(summed_biomes_sampling, biomes_color_key, by = c("biome_name" = "biomes_list"))

summed_biomes_colored_long_plot_sampling <- ggplot(data=summed_biomes_sampling_color, aes(x=fct_reorder(biome_name, as.numeric(freq_count)), y=as.numeric(freq_count), fill =biomes_color)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label = freq_count), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide='legend', labels=summed_biomes_sampling_color$biome_name, breaks=summed_biomes_sampling_color$biomes_color) +
  coord_flip(clip = "off") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(), legend.background = element_rect(fill = "transparent"),             # Transparent legend bg
        legend.box.background = element_rect(fill = "transparent"),         # Transparent legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA), # Transparent panel bg (the area where data is plotted)
        plot.background = element_rect(fill = "transparent", colour = NA) ) 

summed_biomes_colored_long_plot_sampling

ggsave("2026.01.13.summed_biomes_colored_long_plot_sampling.png", summed_biomes_colored_long_plot_sampling, width = 30, height = 20, units = "cm", limitsize = F, bg='transparent')

# summed_biomes_colored_long_plot_zoom_sampling <- summed_biomes_colored_long_plot_sampling + 
#   coord_flip(ylim = c(0, 5000)) +
#   labs(title = "Zoomed View (0-10)")  +
#   theme_bw() +theme(text = element_text(family = "Noto Sans"))  +
#   theme(axis.title = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.border = element_blank())
# 
# summed_biomes_colored_long_plot_zoom_sampling
# 
# ggsave("summed_biomes_colored_long_plot_zoom.png", summed_biomes_colored_long_plot_zoom, width = 40, height = 20, units = "cm", limitsize = F,bg='transparent')


##normalize
# long$plot_cat <- gsub("Salween River|Gulf of Alaska Coastal Rivers", "Terrestrial Other", long$plot_cat)


long_normalized <- left_join(long, summed_biomes_sampling, by = "biome_name")
long_normalized$species_per_10k <- (as.numeric(long_normalized$value)/(as.numeric(long_normalized$freq_count)/10000))

ggplot(data=long_normalized, aes(x=fct_reorder(biome_name,as.numeric(species_per_10k)) , y=as.numeric(species_per_10k), fill=variable)) +
  geom_bar(stat="identity") + coord_flip()

library(jsonlite)
palma <- read.csv("https://ourworldindata.org/grapher/palma-ratio-s90s40-ratio.csv?v=1&csvType=full&useColumnShortNames=true")

# Fetch the metadata
metadata <- fromJSON("https://ourworldindata.org/grapher/palma-ratio-s90s40-ratio.metadata.json?v=1&csvType=full&useColumnShortNames=true")

palma_recent <- palma %>%
  group_by(Entity) %>%
  arrange(desc(Year)) %>%
  slice(1)

world_countries <- ne_countries(scale = "medium", returnclass = "sf")
world_countries <- select(world_countries, admin, geometry, adm0_a3)

palma_recent_geo <- left_join(palma_recent, world_countries, by = c("Entity" = "admin"))
palma_recent_geo <- left_join(palma_recent_geo, world_countries, by = c("Code" = "adm0_a3"))

palma_recent_geo_nas <- palma_recent_geo[!is.na(palma_recent_geo$admin),]
palma_recent_geo_nas$geometry <- coalesce(palma_recent_geo_nas$geometry.x, palma_recent_geo_nas$geometry.y)
palma_recent_geo_nas$geometry.x <- NULL
palma_recent_geo_nas$geometry.y <- NULL

df_palma <- st_as_sf(giant_geo_table_grid_id_geometry)  %>% st_transform("ESRI:54009")
palma_recent_geo_nas <- st_as_sf(palma_recent_geo_nas) %>% st_transform("ESRI:54009")

df_palma <- st_join(df_palma, palma_recent_geo_nas, join = st_intersects)
df_palma_info <- select(df_palma, Entity, Code, palma_ratio_pretax, Freq) %>%
  group_by(Entity, Code, palma_ratio_pretax) %>%
  summarize(total_freq = sum(Freq), .groups = 'drop')

palma_plot <- ggplot(data = df_palma_info, aes(x = log(total_freq), y = log(palma_ratio_pretax))) + geom_point() + geom_smooth(method = 'lm') + stat_poly_eq(mapping = use_label("R2", "p")) +
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  + geom_text_repel(aes(label = Entity)) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(), legend.background = element_rect(fill = "transparent"),             # Transparent legend bg
        legend.box.background = element_rect(fill = "transparent"),         # Transparent legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA), # Transparent panel bg (the area where data is plotted)
        plot.background = element_rect(fill = "transparent", colour = NA) ) 

ggsave("palma_plot.png", palma_plot, width = 20, height = 20, units = "cm", limitsize = F,bg='transparent')


hdi_2023 <- read.table("hdi_2023.txt", sep = "\t")
hdi_2023$V1 <- country_name(x= hdi_2023$V1, to="ISO3", fuzzy_match = T)

df_hdi <- left_join(as.data.frame(df_palma), hdi_2023, by = c("Code" = "V1"))

library(countries)
df_hdi$code3 <- country_name(x= df_hdi$admin, to="ISO3", fuzzy_match = T)

df_hdi_info <- select(df_hdi, Entity, Code, V2, Freq) %>%
  group_by(Entity, Code, V2) %>%
  summarize(total_freq = sum(Freq), .groups = 'drop')

hdi <- ggplot(data = df_hdi_info, aes(y = total_freq+0.01, x = V2)) + 
  geom_point() + 
  geom_smooth(method = 'lm') +
  geom_text_repel(aes(label = Code)) +
  stat_poly_eq(mapping = use_label("R2", "p"))  + theme_bw() +
  theme(panel.border = element_blank(), legend.background = element_rect(fill = "transparent"),             # Transparent legend bg
        legend.box.background = element_rect(fill = "transparent"),         # Transparent legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA), # Transparent panel bg (the area where data is plotted)
        plot.background = element_rect(fill = "transparent", colour = NA) ) 

ggsave("hdi_plot.png", hdi, width = 20, height = 20, units = "cm", limitsize = F,bg='transparent')

hdi_cut <- ggplot(data = df_hdi_info, aes(y = total_freq+0.01, x = V2)) + 
  ylim(c(0, 300000)) + 
  geom_point() + 
  geom_smooth(method = 'lm') +
  geom_text_repel(aes(label = Code)) +
  stat_poly_eq(mapping = use_label("R2", "p"))  + theme_bw() +
  theme(panel.border = element_blank(), legend.background = element_rect(fill = "transparent"),             # Transparent legend bg
        legend.box.background = element_rect(fill = "transparent"),         # Transparent legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA), # Transparent panel bg (the area where data is plotted)
        plot.background = element_rect(fill = "transparent", colour = NA) ) 

hdi_cut

ggsave("hdi_plot_cut.png", hdi_cut, width = 20, height = 20, units = "cm", limitsize = F,bg='transparent')

#ordered biome
biome_count_order <- levels(fct_reorder(summed_biomes_sampling_color$biome_name, as.numeric(summed_biomes_sampling_color$freq_count)))

summed_biomes_sampling_color_pv_ordered <- filter(summed_biomes_sampling_color_pv, !is.na(proportions.y))

resevoir <- c("Reservoir", NA, "white", NA, NA, NA,0, NA)
alkaline <- c("Alkaline Lake", NA, "white", NA, NA, NA, 0, NA)
summed_biomes_sampling_color_pv_ordered <- rbind(summed_biomes_sampling_color_pv_ordered[,-2], alkaline, resevoir)

summed_biomes_sampling_color_pv_ordered <- summed_biomes_sampling_color_pv_ordered %>%
  mutate(biome_name = factor(biome_name, levels = biome_count_order))


summed_biomes_colored_long_plot_sampling_pv_ordered <- ggplot(data=filter(summed_biomes_sampling_color_pv_ordered, !is.na(proportions.x)), aes(x=biome_name, y=as.numeric(sumz), fill =biomes_color.x)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label = as.numeric(sumz)), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide='legend', labels=summed_biomes_sampling_color$biome_name, breaks=summed_biomes_sampling_color$biomes_color) +
  coord_flip(clip = "off") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

summed_biomes_colored_long_plot_sampling_pv_ordered

ggsave("2026.01.13.summed_biomes_colored_long_plot_sampling_pv_ordered.png", summed_biomes_colored_long_plot_sampling_pv_ordered, width = 30, height = 20, units = "cm", limitsize = F,bg='transparent')


summed_biomes_colored_long_plot_sampling_pv_ordered_grey <- ggplot(data=filter(summed_biomes_sampling_color_pv_ordered, !is.na(proportions.x)), aes(x=biome_name, y=as.numeric(sumz), fill =biomes_color.x)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label = as.numeric(signif(as.numeric(proportions.y), 2))), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide='legend', labels=summed_biomes_sampling_color$biome_name, breaks=summed_biomes_sampling_color$biomes_color) +
  coord_flip(clip = "off") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

knwon <- filter(summed_biomes_colored_long, variable == 'known')

wuhhh <- summed_biomes_colored_long_plot_sampling_pv_ordered_grey + geom_bar(data=filter(summed_biomes_colored_long, variable == 'known'), stat="identity", fill = "grey90", alpha = 0.8, aes(y = value))


ggsave("2026.01.13.summed_biomes_colored_long_plot_sampling_pv_ordered_grey.png", wuhhh, width = 30, height = 20, units = "cm", limitsize = F,bg='transparent')


#normalize sampling based on biome area

summed_biomes_sampling_area <- left_join(summed_biomes_sampling_color, st_drop_geometry(biome_summary), c("biome_name" = "biome_name"))
summed_biomes_sampling_area$sampling_over_area <- summed_biomes_sampling_area$freq_count/summed_biomes_sampling_area$area_km2
# 
# ggplot(data = summed_biomes_sampling_area, aes(x = area_km2, y = freq_count)) + 
#   scale_x_log10() +
#   scale_y_log10() +
#   geom_smooth(method = 'lm') +
#   geom_point() + geom_text(aes(label = guess_1))
# 
# 
# biomes_list <- c("Tropical & Subtropical Moist Broadleaf Forests", 
#                  "Tropical & Subtropical Dry Broadleaf Forests",
#                  "Tropical & Subtropical Coniferous Forests",
#                  "Temperate Broadleaf & Mixed Forests",
#                  "Temperate Conifer Forests",
#                  "Boreal Forests/Taiga",
#                  "Tropical & Subtropical Grasslands, Savannas & Shrublands",
#                  "Temperate Grasslands, Savannas & Shrublands",
#                  "Flooded Grasslands & Savannas",
#                  "Montane Grasslands & Shrublands",
#                  "Tundra",
#                  "Mediterranean Forests, Woodlands & Scrub",
#                  "Deserts & Xeric Shrublands",
#                  "Mangroves", "Terrestrial Other", "Ocean", "Lake", "Alkaline Lake", "Reservoir")
# biomes_color <- c('#008346','#9DCC00', '#C4B72E',
#                   '#015C31','#20BC10','#FFA8BB'
#                   ,'#FAD505','#8F7C00','#5AD488',
#                   '#993E01','#C20088', '#81A527'
#                   ,'#FFA405','#FFCC99', "#624E38", "#09BEE9","#6E96D0","#284A7C", "#8CD4D4")
# biomes_color_key <- as.data.frame(cbind(biomes_list, biomes_color))
# 

land <- ne_download(scale = 50, type = 'land', category = 'physical', returnclass = "sf")
land <- st_transform(land, st_crs(grid_sf))

# Check which grid cells overlap with land
intersects_list <- st_intersects(grid_sf, land)

grid_sf$overlaps_land <- lengths(intersects_list) > 0

grid_sf <- grid_sf %>%
  mutate(biome_name = case_when(
    biome_name == "Terrestrial Other" & !overlaps_land ~ "Ocean",
    TRUE ~ biome_name
  ))

summed_biomes_colored_world <- left_join(grid_sf, biomes_color_key, by = c("biome_name" = "biomes_list"))

wurld_biomes <- ggplot(data =summed_biomes_colored_world, aes(geometry = geometry, fill=biomes_color)) + 
  geom_sf(color = "transparent") + 
  scale_fill_identity(guide='legend', labels=summed_biomes_colored_world$biome_name, breaks=summed_biomes_colored_world$biomes_color) +
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank())


library(ggrepel)


lin <- ggplot(data = summed_biomes_sampling_area, aes(x = as.numeric(area_km2), y = as.numeric(freq_count))) + 
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = 'lm', color = 'black', fill = 'grey85') +
    geom_point(aes(color = biomes_color)) + geom_text_repel(aes(label = biome_name, color = biomes_color)) + 
  scale_color_identity(
    name = "Biomes",               # Legend title
    breaks = summed_biomes_sampling_area$biome_name, # Which values to include in the legend
    labels = summed_biomes_sampling_area$biome_name               # Request a legend
  ) + 
  stat_poly_eq(mapping = use_label("R2", "p")) +
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank())

lin 

ggsave("2026.01.12.change_da_world.png", wurld_biomes, width = 40, height = 30, units = "cm", limitsize = F,bg='transparent')
ggsave("2026.01.13.lin.png", lin, width = 20, height = 20, units = "cm", limitsize = F,bg='transparent')

