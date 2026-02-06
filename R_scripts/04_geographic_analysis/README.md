Scripts here are used to understand the geographic and biome-based distribution of PVs. This folder has run order dependencies: 

Claude:
```
1. **poolygons.R** (run first)
   - Creates: `grid_sf`, `world`, `data_wide`, `giant_geo_table_grid_id_geometry`
   - Outputs: `grids_sf.gpkg`
   - **Note**: Line 17 and 22 reference `all_pvs_mapping_binned` and `all_pvs_mapping`, which is generated in `map_gen_for_pub.R`

2. **geo_analysis.R** (depends on poolygons.R)
   - Requires: `giant_geo_table_grid_id_geometry`, `grid_sf`, `world`, `data_wide`

3. **biomes.R** (depends on poolygons.R)
   - Reads: `grids_sf.gpkg` (output from poolygons.R)
   - Requires: `data_wide`, `grid_sf`, `world`
```
I apologize for how interconnected things are!