#geo data
source("R_scripts/04_geographic_analysis/geo_setup.R")

###for novelty search
all_seq_geo <- read.table("files/all_full_l1s_logan_and_pr_blastn.tsv", sep = "\t", header = F)
all_seq_geo$library <- gsub("Score.*", "", all_seq_geo$V1)
#calculate qcov
all_seq_geo$qcov <- abs(all_seq_geo$V2 - all_seq_geo$V3)/all_seq_geo$V4

#evals fil 
all_seq_geo_fil <- filter(all_seq_geo, V10 < 0.00001)
all_seq_geo_fil_low_conf <- filter(all_seq_geo, V10 > 0.00001)

what_hits <- as.data.frame(table(all_seq_geo$V1))


#look for those with highest % iden first in confident hits
all_seq_geo_sliced <- all_seq_geo_fil %>% group_by(V1) %>% slice_max(n = 1, V9)
#under 90% iden
all_seq_geo_sliced_90 <- filter(all_seq_geo_sliced, V9 < 90)
all_seq_geo_sliced_90 = all_seq_geo_sliced_90[!duplicated(all_seq_geo_sliced_90$V1),]

all_seq_geo_less_90_nt <- unique(all_seq_geo_sliced_90$V1)
# all_seq_geo_qcov_less50 <- filter(all_seq_geo_qcov_less40, !V1 %in% all_seq_geo_less_90_nt)
length(unique(all_seq_geo_less_90_nt))


# write.table(all_seq_geo_less_90_nt, "all_seq_geo_less_90_nt_ugh.txt", quote = F, col.names = F, row.names = F)

hist(all_seq_geo_sliced_90$V9)


##some are not here!
length(unique(all_seq_geo$V1))
# 859 of 942, so 83 are unaccounted for

#create list of the ones that were hit, in general
hit_list <- unique(all_seq_geo_fil$V1)
# write.table(hit_list, "all_seq_geo_hit_list.txt", quote = F, col.names = F, row.names = F)

#look for low conf hits that were not represented in either novel already, or other filter set
all_seq_geo_fil_low_conf_hits <- filter(all_seq_geo_fil_low_conf, !V1 %in% all_seq_geo_less_90_nt)
all_seq_geo_fil_low_conf_hits <- filter(all_seq_geo_fil_low_conf_hits, !V1 %in% all_seq_geo_fil$V1)
length(unique(all_seq_geo_fil_low_conf_hits$library))

all_seq_geo_fil_low_conf_hits_list <- unique(all_seq_geo_fil_low_conf_hits$V1)

length(unique(all_seq_geo_fil_low_conf_hits$V1))

# write.table(all_seq_geo_fil_low_conf_hits_list, "all_seq_geo_fil_low_conf_hits_list_ugh.txt", quote = F, col.names = F, row.names = F)
#group all the ones with hits that are novel
all_hit_novel <- bind_rows(all_seq_geo_fil_low_conf_hits, all_seq_geo_sliced_90)
# write.table(all_hit_novel, "all_hit_novel.txt", quote = F, col.names = F, row.names = F)
length(unique(all_hit_novel$V1))

####read in geo data and full list of hits
all_input_into_blastn <- read.table("files/run_in_blastn.list", sep = "\t", header = F)
all_input_into_blastn$V1 <- gsub(">", "", all_input_into_blastn$V1)
all_input_into_blastn_joined <- left_join(all_input_into_blastn, all_seq_geo, by = c("V1"))
all_input_into_blastn_joined_min_eval <- all_input_into_blastn_joined %>% group_by(V1) %>% slice_min(n = 1, V10) %>% slice_max(n = 1, V9)
all_input_into_blastn_joined_min_eval <- all_input_into_blastn_joined_min_eval[!duplicated(all_input_into_blastn_joined_min_eval[c(1,9,10)]),]
length(unique(all_input_into_blastn_joined_min_eval$V1))

no_hits_in_blastn <- all_input_into_blastn_joined_min_eval %>% filter(if_any(everything(), is.na))
known_pvs <- filter(all_input_into_blastn_joined_min_eval, V9 > 90)
known_pvs$Run <- gsub("_.*", "", known_pvs$V1)
known_pvs$Run <- gsub("Score.*", "", known_pvs$Run)

novel_pvs <- filter(all_input_into_blastn_joined_min_eval, V9 < 90)

all_novels <- bind_rows(no_hits_in_blastn, novel_pvs)
all_novels$Run <- gsub("_.*", "", all_novels$V1)
all_novels$Run <- gsub("Score.*", "", all_novels$Run)

length(unique(all_novels$V1))

#add geo data to known and unknown
geo_data_annotation_for_all_biosamps <- read.csv("files/geo_data_annotation_for_all_biosamps.txt", header = F)


all_hits_library_biosample <- read.table("files/who_puts_vlookup_man.list", sep = "\t", header = T, fill = T)


geo_data <- left_join(geo_data_annotation_for_all_biosamps, all_hits_library_biosample, by = c("V1" = "BioSample"))
geo_data <- dplyr::select(geo_data, V1, V2, V4, Run, LibrarySource)


geo_data$V2 <- factor(geo_data$V2, levels = c("lat_lon", "geographic location (latitude),geographic location (longitude)", 
                                              "geo_loc_name_sam", "geo_loc_name", 
                                              "geo_loc_name_country_calc", "geographic location (region and locality)", "region", "birth_location", "INSDC center name"))

geo_data <- geo_data[order(geo_data$V2),]

geo_data_priority <- geo_data[!duplicated(geo_data[,c('Run')]),]

library(sf)

# Example WKB string (hexadecimal representation of a POINT)
# wkb = structure(list("0101000020E61000005859DB140F776140DC4944F817D84140"), class = "WKB")
# st_as_sfc(wkb,EWKB = TRUE) %>% st_coordinates()

wkb_obj <- structure(as.list(geo_data_priority$V4), class = "WKB")

geo_data_priority$geometry <- st_as_sfc(wkb_obj, EWKB = TRUE)

geo_data_priority <- st_as_sf(geo_data_priority, geometry = geo_data_priority$geometry, crs = 4326)
geo_data_priority <- st_transform(geo_data_priority, crs = "ESRI:54009")

all_novels_w_geo <- left_join(all_novels, geo_data_priority, by = c("Run"))
all_novels_w_geo$status <- c("novel")

all_knowns_w_geo <- left_join(known_pvs, geo_data_priority, by = c("Run"))
all_knowns_w_geo$status <- c("known")

all_pvs_mapping <- bind_rows(all_novels_w_geo, all_knowns_w_geo)

sum(is.na(all_pvs_mapping$V1.y))

# world is loaded by geo_setup.R

world_plot <- all_pvs_mapping %>%
  ggplot(aes(color= NA)) +
  geom_sf(data = world, fill = 'grey95', color = 'grey90') +
  geom_sf(aes(geometry = geometry,  color = status, lwd = 0, fill = NA), lwd = 0, alpha = 0.5, stroke = 0, size =1, shape=16) +
  scale_color_manual(values = c("grey40", "#3333E7")) +
  coord_sf(crs = "ESRI:54009") +
  theme_bw() + theme(text = element_text(family = "Noto Sans"))  +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) 

world_plot

ggsave("outputs/2026.01.09geo_data_moll.png", world_plot, width = 25, height = 8, units = "cm", limitsize = F,bg='transparent')

all_pvs_mapping$LibrarySource <- gsub("SYNTHETIC", "OTHER", all_pvs_mapping$LibrarySource)
type_table <- as.data.frame(table(filter(all_pvs_mapping, status == "novel")$LibrarySource))
type_table = type_table[-1,]

# Save key objects for downstream scripts
dir.create("outputs/geo_intermediates", showWarnings = FALSE, recursive = TRUE)
saveRDS(all_pvs_mapping, "outputs/geo_intermediates/all_pvs_mapping.rds")
saveRDS(all_novels, "outputs/geo_intermediates/all_novels.rds")
