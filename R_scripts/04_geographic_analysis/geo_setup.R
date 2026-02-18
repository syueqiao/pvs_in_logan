# geo_setup.R — Shared data loading for geographic analysis scripts
# Source this at the top of any geo analysis script.
# Provides: world, grid_sf
#
# Run order for the full geographic analysis pipeline:
#   1. map_gen_for_pub.R        (novelty classification + distribution map)
#   2. polygon_analyses.R       (grid intersection, species proportions)
#   3. geo_analysis.R           (resampling statistics)
#   4. novelty_significance_maps.R (significance overlay — needs geo_analysis.R output)
#   5. biomes.R                 (biome/ecology analysis — needs polygon_analyses.R output)

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

# --- World map (Mollweide projection) ---
# Used by: all geographic analysis scripts
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, "ESRI:54009")

# --- 50km equal-area grid (terra method) ---
# Mollweide projection, densified to avoid pole singularity artifacts
proj_crs <- "ESRI:54009"
cell_size <- 50000  # 50km

v <- vect(ext(-179.9, 179.9, -89.9, 89.9), crs = "+proj=longlat")
v_dense <- densify(v, interval = 1)
world_ellipse <- project(v_dense, proj_crs)

r <- rast(world_ellipse, res = cell_size)
values(r) <- 1:ncell(r)
r_masked <- mask(r, world_ellipse)

grid_terra <- as.polygons(r_masked, dissolve = FALSE)
grid_sf <- st_as_sf(grid_terra)
grid_sf$id <- seq_len(nrow(grid_sf))

message("geo_setup.R loaded: world, grid_sf (", nrow(grid_sf), " cells)")
