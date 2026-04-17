Scripts for geographic and biome-based distribution analysis of PVs. These scripts have run-order dependencies.

**Run order:**

1. `geo_setup.R` — Shared data loading (world map, grid). Sourced by other scripts.
2. `map_gen_for_pub.R` — Novelty classification with geographic coordinates; generates publication-ready distribution maps. Creates `all_pvs_mapping`, `all_novels`, `known_pvs`.
3. `polygon_analyses.R` — Spatial grid intersection and species proportions. Creates `grid_sf`, `data_wide`, `giant_geo_table_grid_id_geometry`. Depends on `map_gen_for_pub.R`.
4. `geo_analysis.R` — Geographic resampling statistics. Depends on `polygon_analyses.R`.
5. `novelty_significance_maps.R` — Significance overlay maps. Depends on `geo_analysis.R`.
6. `biomes.R` — Biome/ecology analysis using WWF ecoregions. Depends on `polygon_analyses.R` and `geo_analysis.R`. Requires `Ecoregions2017.shp` in `files/` (download from https://ecoregions.appspot.com/).
7. `biomes_stats_test.R` — Statistical testing for biome associations.
