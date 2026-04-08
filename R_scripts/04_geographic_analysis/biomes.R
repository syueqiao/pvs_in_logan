source("R_scripts/04_geographic_analysis/geo_setup.R")

# Load from environment (if polygon_analyses.R was run) or from saved RDS
if (!exists("data_wide")) {
  data_wide <- readRDS("outputs/geo_intermediates/data_wide.rds")
}
if (!exists("giant_geo_table_grid_id_geometry")) {
  giant_geo_table_grid_id_geometry <- readRDS("outputs/geo_intermediates/giant_geo_table_grid_id_geometry.rds")
}

# Biome color key (used throughout this script)
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

###ecology?

data_wide_joining <- st_as_sf(data_wide)  %>% st_transform("ESRI:54009")
ocean_polygons <- st_as_sf(ne_download(scale = 50, type = 'ocean', category = 'physical', returnclass = "sf")) %>% st_transform("ESRI:54009") %>%
  rename(featurecla_ocean = featurecla)
fw_polygons <- st_as_sf(ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")) %>% st_transform("ESRI:54009")%>%
  rename(featurecla_fw = featurecla)
fw_polygons_2 <- st_as_sf(ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = "sf")) %>% st_transform("ESRI:54009")%>%
  rename(featurecla_fw_2 = featurecla)

water <- bind_rows(dplyr::select(ocean_polygons, geometry, featurecla_ocean), dplyr::select(fw_polygons, geometry, featurecla_fw), dplyr::select(fw_polygons_2, geometry, featurecla_fw_2))

# my_shapefile <- read_sf("Global_200_Terrestrial.shp")
#the 2017 EcoRegions can be downloaded at https://ecoregions.appspot.com/
# The 2017 EcoRegions shapefile can be downloaded at https://ecoregions.appspot.com/
# Place Ecoregions2017.shp (and associated .dbf, .prj, .shx) in the files/ directory
my_shapefile_2 <- read_sf("files/Ecoregions2017.shp")
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

sf_use_s2(FALSE)
#calculate the area using simplified geometries for performance
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

# plot(summed_biomes$proportions)

library(data.table)

long <- melt(setDT(summed_biomes), id.vars = c("biome_name")) %>% filter(!variable == "geometry") %>% filter(!variable == "proportions")
long$variable <- relevel(long$variable, 'novel')


test_graph <- ggplot(data=long, aes(x=fct_reorder(biome_name, as.numeric(value)), y=as.numeric(value), fill=variable)) +
  geom_bar(stat="identity") + coord_flip()

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

ggsave("outputs/2026.01.13summed_biomes_colored_props.png", summed_biomes_colored_props, width = 30, height = 20, units = "cm", limitsize = F ,bg='transparent')


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

ggsave("outputs/2026.01.13.summed_biomes_colored_long_plot.png", summed_biomes_colored_long_plot, width = 40, height = 20, units = "cm", limitsize = F,bg='transparent')


psh__pv <- ggplot(data = st_buffer(st_as_sf(summed_biomes_colored_long), 50000), aes(color = NA, fill = biomes_color)) + 
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(color = "transparent", alpha = 0.8, aes(geometry = geometry)) +
  scale_fill_identity(guide='legend', labels=summed_biomes_colored_long$biome_name, breaks=summed_biomes_colored_long$biomes_color) +
  coord_sf(crs = "ESRI:54009") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right')

ggsave("outputs/2026.02.18.psh__pv_2.png", psh__pv, width = 40, height = 30, units = "cm", limitsize = F,bg='transparent')

biomes_joined$known_novel <- biomes_joined$known + biomes_joined$novel
biomes_joined <- left_join(biomes_joined, biomes_color_key, by = c("biome_name" = "biomes_list"))

psh__pv_2 <- ggplot(data = st_buffer(st_as_sf(biomes_joined), 50000), aes(color = NA, fill = biomes_color)) + 
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(color = "transparent", aes(geometry = geometry, alpha = 0.1+log10(known_novel))) +
  scale_fill_identity(guide='legend', labels=biomes_joined$biome_name, breaks=biomes_joined$biomes_color) +
  coord_sf(crs = "ESRI:54009") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right')

psh__pv_2
ggsave("outputs/2026.02.18.psh__pv_alpha.png", psh__pv_2, width = 40, height = 30, units = "cm", limitsize = F,bg='transparent')

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

sum(summed_biomes_sampling_color_pv$value)

ggsave("outputs/2026.01.13.summed_biomes_colored_long_plot_sampling_pv.png", summed_biomes_colored_long_plot_sampling_pv, width = 40, height = 20, units = "cm", limitsize = F,bg='transparent')


# biomes_joined_aqua_guess_nz_color_2$plot_cat <- gsub("Salween River|Gulf of Alaska Coastal Rivers", "Terrestrial Other", biomes_joined_aqua_guess_nz_color_2$plot_cat)

#normalize by sampling?##########?##########?##########?##########?##########?##########?##########?##########?##########?##########?##########

df_fixed_joining <- st_as_sf(giant_geo_table_grid_id_geometry)  %>% st_transform("ESRI:54009")

# df_fixed_joining <- filter(df_fixed_joining, Freq.x > 0)

df_fixed_joining_centroids <- st_centroid(df_fixed_joining)

df_fixed_joining$biome_id <- st_intersects(df_fixed_joining_centroids, all_cats) %>%
  sapply(function(x) if(length(x) > 0) x[1] else NA)

# Add biome name (adjust column name as needed)
df_fixed_joining$biome_name <- all_cats$area_biome[df_fixed_joining$biome_id]

df_fixed_joining$biome_name[is.na(df_fixed_joining$biome_name)] <- "Terrestrial Other"

biomes_joined_sampling <- df_fixed_joining

st_write(obj = biomes_joined_sampling, dsn = "outputs/2026.01.13.biomes_joined_sampling.gpkg", driver = "GPKG", append=FALSE)

# biomes_joined_sampling <- st_read("biomes_joined_sampling.gpkg")

# biomes_joined_sampling$guess_1 <- coalesce(biomes_joined_sampling$BIOME_NAME, biomes_joined_sampling$featurecla_ocean, biomes_joined_sampling$featurecla_fw_2, biomes_joined_sampling$featurecla_fw)

# biomes_joined_sampling$guess_1[is.na(biomes_joined_sampling$guess_1)] <- "Terrestrial Other"

biomes_joined_sampling_nz <- filter(biomes_joined_sampling,  Freq > 0) 

head(biomes_joined_sampling)
# biomes_color_key defined at top of script

# biomes_joined_sampling_nz <- biomes_joined_sampling_nz %>% rename(plot_cat = guess_1)

biomes_joined_sampling_nz_color <- left_join(biomes_joined_sampling_nz, biomes_color_key, by = c("biome_name" = "biomes_list"))

psh <- ggplot(data = biomes_joined_sampling_nz_color, aes(color = NA, fill = biomes_color)) + 
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(color = "transparent", aes(geometry = geometry)) +
  scale_fill_identity(guide='legend', labels=biomes_joined_sampling_nz_color$biome_name, breaks=biomes_joined_sampling_nz_color$biomes_color) +
  coord_sf(crs = "ESRI:54009") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right')

psh

ggsave("outputs/2026.02.18.psh.png", psh, width = 40, height = 30, units = "cm", limitsize = F,bg='transparent')


summed_biomes_sampling <- dplyr::select(biomes_joined_sampling_nz_color, biome_name, Freq, biomes_color) %>%
  group_by(biome_name) %>%
  summarize(freq_count = sum(Freq, na.rm = TRUE))

biome_freq_plot <- ggplot(data=summed_biomes_sampling, aes(x=fct_reorder(biome_name, as.numeric(freq_count)), y=as.numeric(freq_count))) +
  geom_bar(stat="identity") + coord_flip()

biome_freq_plot

summed_biomes_sampling_color <- left_join(summed_biomes_sampling, biomes_color_key, by = c("biome_name" = "biomes_list"))

summed_biomes_colored_long_plot_sampling <- ggplot(data=summed_biomes_sampling_color, aes(x=fct_reorder(biome_name, as.numeric(freq_count)), y=as.numeric(freq_count), fill =biomes_color)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label = scales::comma(freq_count)), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide='legend', labels=summed_biomes_sampling_color$biome_name, breaks=summed_biomes_sampling_color$biomes_color) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  coord_flip(clip = "off") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(), legend.background = element_rect(fill = "transparent"),             # Transparent legend bg
        legend.box.background = element_rect(fill = "transparent"),         # Transparent legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA), # Transparent panel bg (the area where data is plotted)
        plot.background = element_rect(fill = "transparent", colour = NA), legend.position = 'none' ) 

summed_biomes_colored_long_plot_sampling

ggsave("outputs/2026.02.18.summed_biomes_colored_long_plot_sampling.png", summed_biomes_colored_long_plot_sampling, width = 20, height = 10, units = "cm", limitsize = F, bg='transparent')
ggsave("outputs/2026.03.24.summed_biomes_colored_long_plot_sampling.png", summed_biomes_colored_long_plot_sampling, width = 15, height = 8, units = "cm", limitsize = F, bg='transparent')

##normalize


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

ggsave("outputs/palma_plot.png", palma_plot, width = 20, height = 20, units = "cm", limitsize = F,bg='transparent')


hdi_2023 <- read.table("files/hdi_2023.txt", sep = "\t")
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

ggsave("outputs/hdi_plot.png", hdi, width = 20, height = 20, units = "cm", limitsize = F,bg='transparent')

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

ggsave("outputs/hdi_plot_cut.png", hdi_cut, width = 20, height = 20, units = "cm", limitsize = F,bg='transparent')

#ordered biome
biome_count_order <- levels(fct_reorder(summed_biomes_sampling_color$biome_name, as.numeric(summed_biomes_sampling_color$freq_count)))

summed_biomes_sampling_color_pv_ordered <- filter(summed_biomes_sampling_color_pv, !is.na(proportions.y))

resevoir <- c("Reservoir", "bwip", "white", NA, NA, NA,0, NA)
alkaline <- c("Alkaline Lake", "bwip", "white", NA, NA, NA, 0, NA)
summed_biomes_sampling_color_pv_ordered <- rbind(summed_biomes_sampling_color_pv_ordered[,-2], alkaline, resevoir)

summed_biomes_sampling_color_pv_ordered <- summed_biomes_sampling_color_pv_ordered %>%
  mutate(biome_name = factor(biome_name, levels = biome_count_order))


summed_biomes_colored_long_plot_sampling_pv_ordered <- ggplot(data=filter(summed_biomes_sampling_color_pv_ordered, !is.na(proportions.x)), aes(x=biome_name, y=as.numeric(sumz), fill =biomes_color.x)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label = scales::comma(as.numeric(sumz))), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide='legend', labels=summed_biomes_sampling_color$biome_name, breaks=summed_biomes_sampling_color$biomes_color) +
  coord_flip(clip = "off") + theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(), legend.position = 'none') 

summed_biomes_colored_long_plot_sampling_pv_ordered
sum(as.numeric(summed_biomes_sampling_color_pv_ordered$value), na.rm = T)

ggsave("outputs/2026.02.18.summed_biomes_colored_long_plot_sampling_pv_ordered.png", summed_biomes_colored_long_plot_sampling_pv_ordered, width = 30, height = 10, units = "cm", limitsize = F,bg='transparent')
ggsave("outputs/2026.02.18.summed_biomes_colored_long_plot_sampling_pv_ordered_short.png", summed_biomes_colored_long_plot_sampling_pv_ordered, width = 15, height = 8, units = "cm", limitsize = F,bg='transparent')


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


ggsave("outputs/2026.01.13.summed_biomes_colored_long_plot_sampling_pv_ordered_grey.png", wuhhh, width = 30, height = 20, units = "cm", limitsize = F,bg='transparent')


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
library(ggExtra)
  library(ggpmisc)

lin <- ggplot(data = summed_biomes_sampling_area, aes(x = as.numeric(area_km2), y = as.numeric(freq_count))) + 
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = 'lm', color = 'black', fill = 'grey85') +
    geom_point(aes(color = biomes_color)) + geom_text_repel(aes(label = biome_name, color = biomes_color), size = 3) + 
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

ggsave("outputs/2026.01.12.change_da_world.png", wurld_biomes, width = 40, height = 30, units = "cm", limitsize = F,bg='transparent')
ggsave("outputs/2026.01.13.lin_medium.png", lin, width = 15, height = 12, units = "cm", limitsize = F,bg='transparent')

####overlay###

###overlay: sampling biomes (faded) + PV locations (vivid)####

psh_sampling_pv_overlay <- ggplot() +
  # World basemap
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  # Underlay: all sampled grid cells colored by biome, low alpha (grey overtone)
  geom_sf(data = biomes_joined_sampling_nz_color,
          aes(geometry = geometry, fill = biomes_color),
          color = "transparent", alpha = 0.1) +
  geom_sf(data = biomes_joined_sampling_nz_color,
          aes(geometry = geometry), fill = "grey40",
          color = "transparent", alpha = 0.2) +
  # Overlay: grid cells where PVs were detected, full biome color
  geom_sf(data = st_buffer(biomes_joined, 100000),
          aes(geometry = geometry, fill = biomes_color),
          color = "white", alpha = 0.9, linewidth = 0.25) +
  scale_fill_identity(
    guide = 'legend',
    labels = biomes_color_key$biomes_list,
    breaks = biomes_color_key$biomes_color
  ) +
  coord_sf(crs = "ESRI:54009") +
  theme_bw() +
  theme(text = element_text(family = "Noto Sans"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right')

psh_sampling_pv_overlay

ggsave("outputs/2026.02.24.psh_sampling_pv_overlay.png", psh_sampling_pv_overlay,
       width = 40, height = 30, units = "cm", limitsize = FALSE, bg = 'transparent')

###############################################################################################
# Same analysis but with "INSDC center name" geolocations filtered out
###############################################################################################

# Load all_pvs_mapping_clean to get geo source info (V2 column)
if (!exists("all_pvs_mapping_clean")) {
  all_pvs_mapping_clean <- readRDS("outputs/geo_intermediates/all_pvs_mapping_clean.rds")
}

# Remove individual PVs whose geo source is "INSDC center name", then rebuild counts
# (Previously only removed grid cells where ALL PVs were INSDC — this removes each INSDC-sourced PV)
all_pvs_noinsdc <- all_pvs_mapping_clean %>%
  filter(V2 != "INSDC center name")

cat("PVs removed (INSDC center name source):", nrow(all_pvs_mapping_clean) - nrow(all_pvs_noinsdc), "\n")

# Rebuild data_wide without INSDC PVs
noinsdc_props <- st_drop_geometry(
  dplyr::select(st_as_sf(all_pvs_noinsdc), status, cluster_number_updated, grid_group)
)
noinsdc_table <- as.data.frame(with(noinsdc_props, table(grid_group, status)))

data_wide_noinsdc <- noinsdc_table %>%
  pivot_wider(id_cols = grid_group, names_from = status, values_from = Freq, values_fill = 0)
data_wide_noinsdc$novel_species_proportion <- data_wide_noinsdc$novel / (data_wide_noinsdc$known + data_wide_noinsdc$novel)
data_wide_noinsdc$grid_group <- as.numeric(paste(data_wide_noinsdc$grid_group))

# Rejoin grid geometry
grid_geom_lookup <- st_sf(grid_group = grid_sf$id, geometry = st_geometry(grid_sf))
data_wide_noinsdc <- left_join(data_wide_noinsdc, grid_geom_lookup, by = "grid_group")

# Drop grid cells with 0 PVs after filtering
data_wide_noinsdc <- data_wide_noinsdc %>% filter((known + novel) > 0)

data_wide_noinsdc_joining <- st_as_sf(data_wide_noinsdc) %>% st_transform("ESRI:54009")

data_wide_noinsdc_joining_centroids <- st_centroid(data_wide_noinsdc_joining)
data_wide_noinsdc_joining$biome_id <- st_intersects(data_wide_noinsdc_joining_centroids, all_cats) %>%
  sapply(function(x) if(length(x) > 0) x[1] else NA)
data_wide_noinsdc_joining$biome_name <- all_cats$area_biome[data_wide_noinsdc_joining$biome_id]
data_wide_noinsdc_joining$biome_name[is.na(data_wide_noinsdc_joining$biome_name)] <- "Terrestrial Other"

biomes_joined_noinsdc <- data_wide_noinsdc_joining

# Summarise by biome
summed_biomes_noinsdc <- dplyr::select(biomes_joined_noinsdc, biome_name, known, novel) %>%
  group_by(biome_name) %>% summarise(across(1:2, sum))
summed_biomes_noinsdc$proportions <- summed_biomes_noinsdc$novel / summed_biomes_noinsdc$known

# Bar plot
summed_biomes_noinsdc_colored <- left_join(summed_biomes_noinsdc, biomes_color_key, by = c("biome_name" = "biomes_list"))
summed_biomes_noinsdc_colored$proportions <- summed_biomes_noinsdc_colored$novel / summed_biomes_noinsdc_colored$known

summed_biomes_noinsdc_colored_long <- summed_biomes_noinsdc_colored %>%
  pivot_longer(cols = c(known, novel), names_to = "variable", values_to = "value")
summed_biomes_noinsdc_colored_long <- left_join(summed_biomes_noinsdc_colored_long, dplyr::select(summed_biomes_noinsdc_colored, biome_name), by = "biome_name")
summed_biomes_noinsdc_colored_long_sum <- summed_biomes_noinsdc_colored_long %>%
  group_by(biome_name, proportions) %>%
  summarize(sumz = sum(value), .groups = 'drop')
summed_biomes_noinsdc_colored_long <- unique(left_join(summed_biomes_noinsdc_colored_long, summed_biomes_noinsdc_colored_long_sum, by = "biome_name"))
summed_biomes_noinsdc_colored_long$proportions.x[seq_len(nrow(summed_biomes_noinsdc_colored_long)) %% 2 == 0] <- NA

summed_biomes_noinsdc_long_plot <- ggplot(data = summed_biomes_noinsdc_colored_long,
       aes(x = fct_reorder(biome_name, as.numeric(value)), y = as.numeric(value), fill = biomes_color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = sumz, x = fct_reorder(biome_name, as.numeric(sumz)), label = signif(proportions.x, 2)), size = 2, hjust = -0.3, fontface = 'bold') +
  geom_bar(data = filter(summed_biomes_noinsdc_colored_long, variable == 'known'), stat = "identity", fill = "grey90", alpha = 0.8) +
  scale_fill_identity(guide = 'legend', labels = summed_biomes_noinsdc_colored_long$biome_name, breaks = summed_biomes_noinsdc_colored_long$biomes_color) +
  theme_bw() + coord_flip() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())

ggsave("outputs/2026.02.18.summed_biomes_noinsdc_long_plot.png", summed_biomes_noinsdc_long_plot, width = 40, height = 20, units = "cm", limitsize = F, bg = 'transparent')

# Proportion bar plot (no INSDC)
summed_biomes_noinsdc_colored_props <- ggplot(data = summed_biomes_noinsdc_colored,
       aes(x = fct_reorder(biome_name, as.numeric(proportions)), y = as.numeric(proportions), fill = biomes_color)) +
  geom_text(aes(label = signif(proportions, 2)), hjust = -0.1, size = 2, fontface = 'bold') +
  geom_bar(stat = "identity") +
  scale_fill_identity(guide = 'legend', labels = summed_biomes_noinsdc_colored$biome_name, breaks = summed_biomes_noinsdc_colored$biomes_color) +
  coord_flip() + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())

ggsave("outputs/2026.02.18.summed_biomes_noinsdc_props.png", summed_biomes_noinsdc_colored_props, width = 30, height = 20, units = "cm", limitsize = F, bg = 'transparent')

# Map plot (no INSDC)
biomes_joined_noinsdc$known_novel <- biomes_joined_noinsdc$known + biomes_joined_noinsdc$novel
biomes_joined_noinsdc <- left_join(biomes_joined_noinsdc, biomes_color_key, by = c("biome_name" = "biomes_list"))

psh__pv_noinsdc <- ggplot(data = st_buffer(st_as_sf(biomes_joined_noinsdc), 50000), aes(color = NA, fill = biomes_color)) +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(color = "transparent", aes(geometry = geometry, alpha = known_novel)) +
  scale_fill_identity(guide = 'legend', labels = biomes_joined_noinsdc$biome_name, breaks = biomes_joined_noinsdc$biomes_color) +
  coord_sf(crs = "ESRI:54009") + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), legend.position = 'right')

psh__pv_noinsdc

ggsave("outputs/2026.02.18.psh__pv_noinsdc.png", psh__pv_noinsdc, width = 40, height = 30, units = "cm", limitsize = F, bg = 'transparent')

# PV sampling bar plot — ordered by total PVs (no INSDC)
summed_biomes_noinsdc_sampling_color_pv <- left_join(summed_biomes_noinsdc_colored_long, biomes_color_key, by = c("biome_name" = "biomes_list"))

summed_biomes_noinsdc_long_plot_sampling_pv <- ggplot(
  data = filter(summed_biomes_noinsdc_sampling_color_pv, !is.na(proportions.x)),
  aes(x = fct_reorder(biome_name, as.numeric(sumz)), y = as.numeric(sumz), fill = biomes_color.x)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sumz), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide = 'legend', labels = summed_biomes_noinsdc_sampling_color_pv$biome_name, breaks = summed_biomes_noinsdc_sampling_color_pv$biomes_color.x) +
  coord_flip() + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())

cat("Total PVs (no INSDC):", sum(summed_biomes_noinsdc_sampling_color_pv$value, na.rm = TRUE), "\n")

ggsave("outputs/2026.02.18.summed_biomes_noinsdc_long_plot_sampling_pv.png", summed_biomes_noinsdc_long_plot_sampling_pv, width = 40, height = 20, units = "cm", limitsize = F, bg = 'transparent')

# Ordered PV bar plot by sampling count order (no INSDC)
# Use the same biome_count_order from the unfiltered sampling data
summed_biomes_noinsdc_sampling_color_pv_ordered <- as.data.frame(st_drop_geometry(filter(summed_biomes_noinsdc_sampling_color_pv, !is.na(proportions.y))))
summed_biomes_noinsdc_sampling_color_pv_ordered <- summed_biomes_noinsdc_sampling_color_pv_ordered[, -2]
summed_biomes_noinsdc_sampling_color_pv_ordered <- dplyr::select(summed_biomes_noinsdc_sampling_color_pv_ordered, 
                                                            biome_name, variable, biomes_color.x, sumz)

resevoir_ni <- c("Reservoir", "bwip", "white", 0)
alkaline_ni <- c("Alkaline Lake", "bwip", "white", 0)
summed_biomes_noinsdc_sampling_color_pv_ordered <- rbind(summed_biomes_noinsdc_sampling_color_pv_ordered, alkaline_ni, resevoir_ni)
                      
          summed_biomes_noinsdc_sampling_color_pv_ordered <- summed_biomes_noinsdc_sampling_color_pv_ordered %>%
            mutate(biome_name = factor(biome_name, levels = biome_count_order))

summed_biomes_noinsdc_pv_ordered_dedup <- unique(summed_biomes_noinsdc_sampling_color_pv_ordered[, c("biome_name", "biomes_color.x", "sumz")])
summed_biomes_noinsdc_long_plot_sampling_pv_ordered <- ggplot(
  data = summed_biomes_noinsdc_pv_ordered_dedup,
  aes(x = biome_name, y = as.numeric(sumz), fill = biomes_color.x)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = as.numeric(sumz)), hjust = 0, nudge_y = 50, size = 3, fontface = 'bold') +
  scale_fill_identity(guide = 'legend', labels = summed_biomes_sampling_color$biome_name, breaks = summed_biomes_sampling_color$biomes_color) +
  coord_flip(clip = "off") + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())

ggsave("outputs/2026.02.18.summed_biomes_noinsdc_long_plot_sampling_pv_ordered.png", summed_biomes_noinsdc_long_plot_sampling_pv_ordered, width = 30, height = 20, units = "cm", limitsize = F, bg = 'transparent')

# Grey overlay version — ordered, with known greyed out (no INSDC)
summed_biomes_noinsdc_long_plot_sampling_pv_ordered_grey <- ggplot(
  data = filter(summed_biomes_noinsdc_sampling_color_pv_ordered, !is.na(proportions.x)),
  aes(x = biome_name, y = as.numeric(sumz), fill = biomes_color.x)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = as.numeric(signif(as.numeric(proportions.y), 2))), hjust = -0.1, size = 3, fontface = 'bold') +
  scale_fill_identity(guide = 'legend', labels = summed_biomes_sampling_color$biome_name, breaks = summed_biomes_sampling_color$biomes_color) +
  coord_flip(clip = "off") + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())

wuhhh_noinsdc <- summed_biomes_noinsdc_long_plot_sampling_pv_ordered_grey +
  geom_bar(data = filter(summed_biomes_noinsdc_colored_long, variable == 'known'), stat = "identity", fill = "grey90", alpha = 0.8, aes(y = value))

ggsave("outputs/2026.02.18.summed_biomes_noinsdc_long_plot_sampling_pv_ordered_grey.png", wuhhh_noinsdc, width = 30, height = 20, units = "cm", limitsize = F, bg = 'transparent')

# Normalized by sampling — species per 10k samples (no INSDC)
long_noinsdc <- melt(setDT(summed_biomes_noinsdc), id.vars = c("biome_name")) %>%
  filter(!variable == "geometry") %>% filter(!variable == "proportions")
long_noinsdc$variable <- relevel(long_noinsdc$variable, 'novel')
long_noinsdc_normalized <- left_join(long_noinsdc, summed_biomes_sampling, by = "biome_name")
long_noinsdc_normalized$species_per_10k <- (as.numeric(long_noinsdc_normalized$value) / (as.numeric(long_noinsdc_normalized$freq_count) / 10000))

long_noinsdc_normalized_plot <- ggplot(data = long_noinsdc_normalized,
  aes(x = fct_reorder(biome_name, as.numeric(species_per_10k)), y = as.numeric(species_per_10k), fill = variable)) +
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())

ggsave("outputs/2026.02.18.long_noinsdc_normalized_plot.png", long_noinsdc_normalized_plot, width = 30, height = 20, units = "cm", limitsize = F, bg = 'transparent')

##final additional panel

proportions_df <- left_join(summed_biomes_sampling_color, summed_biomes_sampling_color_pv_ordered, by = c("biome_name"))
proportions_df$discovery_ratio <- as.numeric(proportions_df$sumz)/proportions_df$freq_count

proportions_df <- proportions_df %>%
  mutate(biome_name = factor(biome_name, levels = biome_count_order))

proportions_df <- unique(dplyr::select(proportions_df, biome_name, discovery_ratio, biomes_color))

summed_biomes_colored_long_plot_sampling <- ggplot(data=proportions_df, aes(x=biome_name, y=as.numeric(discovery_ratio), fill =biomes_color)) +
geom_bar(stat="identity") +
geom_text(aes(label = ifelse(discovery_ratio == 0, "0", formatC(discovery_ratio, format = "e", digits = 1)),
          hjust = 0, size = 3, fontface = 'bold'), hjust = 0, size = 3, fontface = 'bold') +
  scale_fill_identity(guide = 'legend', labels = summed_biomes_colored_long_plot_sampling$biome_name, breaks = summed_biomes_colored_long_plot_sampling$biomes_color) +
  coord_flip() + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)))+
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())

summed_biomes_colored_long_plot_sampling

ggsave("outputs/2026.03.06.ordered_ratios.png", summed_biomes_colored_long_plot_sampling, width = 20, height = 10, units = "cm", limitsize = F, bg = 'transparent')
