library(tidyverse)
library(gggenes)
library(jsonlite)

addtl_ano <- read.table("2025.10.10all_genus_addtl_anno.tsv.txt", sep = "\t", header = T)
addtl_ano$lib <- gsub("_.*", "", addtl_ano$label)
addtl_ano$lib <- gsub("Score.*", "", addtl_ano$lib)

addtl_ano <- left_join(addtl_ano, sra_metadata, by = c("lib" = "V1"))
addtl_tip <- filter(addtl_ano, isTip == "TRUE")
addtl_tip$broader_gen <- gsub("Cattle", "Bovine", addtl_tip$broader_gen)


addtl_tip$broader_gen <- factor(addtl_tip$broader_gen, levels = c("Human", "Non-human Primate", "Rodent", 
                                                                  "Cetacean", "Bovine", 
                                                                  "Porcine", "Cervine", "Pinnipeds", "Canine", "Feline", "Mustelid", "Equine", "Pangolin", "Bat", "Avian", "Lizard", "Frog", "Ray-finned Fish", "Other"))

addtl_tip$broader_gen <- factor(addtl_tip$broader_gen, levels = c("Human", "Bat", "Bovine", 
                                                                  "Canine", "Ray-finned Fish", 
                                                                  "Rodent", "Avian", "Cetacean", "Cervine", "Frog", "Pangolin", "Feline", "Non-human Primate", "Equine", "Pinnipeds", "Lizard", "Porcine", "Mustelid", "Other"))

addtl_tip_novel <- subset(addtl_tip, grepl("[SED][R]{2}", addtl_tip$label))

addtl_tip_novel_tb <- as.data.frame(table(addtl_tip_novel$broader_gen))
levels(factor(addtl_tip_novel_tb$broader_gen)) 

addtl_tip_novel_tb_ord <- addtl_tip_novel_tb %>%
  mutate(Var1 = forcats::fct_reorder(factor(Var1), Freq))

addtl_ano <- left_join(addtl_ano, sra_metadata, by = c("lib" = "V1"))
addtl_ano$generalization <- iconv(addtl_ano$generalization, from = "UTF-8", to = "ASCII", sub = "")
addtl_ano <- filter(addtl_ano, isTip == "TRUE")


addtl_ano_tree_broad <- select(addtl_tip, label, broader_gen)

rownames(addtl_ano_tree_broad) <- addtl_ano_tree_broad$label
addtl_ano_tree_broad <- select(addtl_ano_tree_broad, broader_gen)


c20 <- c(
  "#D44746", "#F0B2B2", # red
  "#D4B046",
  "#175F67", # purple
  "#8C46D4", # orange
  "#E146D4",
  "#882020", "#D28036", # lt pink
  "#44C8D6", # lt purple
  "#A3E4EB", # lt orange
  "#4F1E80",
  "#941888", "#F0E0B2", "#C099E7", "#3092D4",
  "#53741D", "#C9E349", "#9BBED1", "#CCCCCC"
)

names(c20) <- NULL
color_key_broad <- unique(as.data.frame(addtl_ano_tree_broad$broader_gen))


ordered_gen_broad <- c("Human", "Non-human Primate", "Rodent", 
                       "Cetacean", "Bovine", 
                       "Porcine", "Cervine", "Pinnipeds", "Canine", "Feline", "Mustelid", "Equine", "Pangolin", "Bat", "Avian", "Lizard", "Frog", "Ray-finned Fish", "Other")

color_key_broad <- as.data.frame(color_key_broad[order(sapply(color_key_broad$`addtl_ano_tree_broad$broader_gen`, function(x) which(x == ordered_gen_broad))), ])

color_key_broad$color_values <- c20
colnames(color_key_broad) <- c("gen_label", "color_values")

addtl_ano_tree_broad_vals<- left_join(addtl_ano_tree_broad, color_key_broad, by = c("broader_gen" = "gen_label"))


##################################################PANGO################################
#case study visualizations
pango_contig <- read.table("SRR25256522_663139_pangolin_genome.txt", sep ="\t", header = T)

#set consistent colors for each gene
color_match_genes = setNames(c("#46D452","#99A599","#99A599", "#99A599", "#CBD7CB", "#CBD7CB", "#CBD7CB",  "#CBD7CB", "#B2F0B7"), c("L1", "L2", "E1", "E2", "E4", "E5", "E6", "E7", "JR"))

#SRR25256522_663139	JR	3617	2675	-

ggplot(pango_contig, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = T)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
  geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_text(aes(x = end - ((end-start)/2), y = 2, label = gene, fontface = 'bold', family = "Noto Sans"))  +  scale_fill_manual(values = color_match_genes) +
  scale_fill_manual(values = color_match_genes) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes() + theme(text = element_text(family="Noto Sans", size = 10)) + scale_x_reverse()

ggsave("pango_map.png", width = 20, height = 10, units = "cm", limitsize = F)

ggplot(pango_contig, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = T)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
  geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  # geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_text(aes(x = end - ((end-start)/2), y = 1.2, label = gene, fontface = 'bold', family = "Noto Sans"))  +  scale_fill_manual(values = color_match_genes) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes() + theme(text = element_text(family="Noto Sans", size = 10)) + scale_x_reverse()

ggsave("pango_map_2.png", width = 20, height = 10, units = "cm", limitsize = F)




world_map <- map_data("world")
china_map <- subset(world_map, region=="China")
china_map_sub <- subset(china_map, region=="Xiamen")

library(raster)
library(geodata)
library(tidyterra)

china_outer <- gadm("CHN", level = 0, path = tempdir())

china_bounds <- gadm("CHN", level = 1, path = tempdir())
china_bounds_f <- fortify(china_bounds)

ggplot() +
 geom_map(data = world_map, map = world_map,
             aes(long, lat, map_id = region), color = "grey70", 
             fill = "grey80") +
  geom_sf(data = china_outer, fill = "#EAE5D2", color = "grey70", linewidth = 0.5) +
  geom_sf(data = subset(china_bounds_f, NAME_1 == "Yunnan"), fill = "#F0E0B2", color = "grey60", linewidth = 0.5) +
  coord_sf(xlim = c(70, 160), ylim = c(15, 65), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70")) 

ggplot() +
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80") +
  geom_sf(data = china_outer, fill = "grey70", color = "grey70", linewidth = 1) +
  geom_sf(data = subset(china_bounds_f, NAME_1 == "Yunnan"), fill = "#F0E0B2", color = "grey70", linewidth = 1) +
  # coord_sf(xlim = c(70, 160), ylim = c(15, 65), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70")) 
  
ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = china_map, fill = "lightgrey", color = "lightgrey") +
  # geom_map(map=china_map,aes(map_id=china_map_sub, x=lon, y=lat),fill = "cornflowerblue", colour = "gray") +
  # geom_point(aes(x = geo_prior_all$lat_lons$X, y = geo_prior_all$lat_lons$Y, color = geo_prior_all$status), size = 1.5, alpha = 0.5, stroke = NA) +
  # scale_color_manual(values = c("grey40", "#3333E7")) + 
  # coord_fixed() +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  xlab("") +
  ylab("wee") 


biosamp_data_anno %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = world_map, fill = "lightgrey", color = "lightgrey") +
  geom_point(aes(x = lat_lon.X, y = lat_lon.Y), color = "#4646d4", size = 1, alpha = 0.5) +
  scale_color_identity() +
  coord_fixed() +
  xlab("") +
  ylab("wee") + theme_bw()


##################################################RHINO################################
#case study visualizations
rhino_contig <- read.table("SRR10902309_46767_rhino_genome.txt", sep ="\t", header = T)

#set consistent colors for each gene
color_match_genes = setNames(c("#46D452","#99A599","#99A599", "#99A599", "#CBD7CB", "#CBD7CB", "#CBD7CB",  "#CBD7CB", "#B2F0B7"), c("L1", "L2", "E1", "E2", "E4", "E5", "E6", "E7", "JR"))

#SRR25256522_663139	JR	3617	2675	-

ggplot(rhino_contig, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = T)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
  geom_gene_arrow(aes(xmin = 5543+126, xmax = 5543+1310, y = "SRR10902309_46767"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  # geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_text(aes(x = end - ((end-start)/2), y = 1.2, label = gene, fontface = 'bold', family = "Noto Sans"))  +  scale_fill_manual(values = color_match_genes) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes() + theme(text = element_text(family="Noto Sans", size = 10))

ggsave("rhino_map.png", width = 20, height = 10, units = "cm", limitsize = F)

rhino_subset <- tree_subset(final_tree, "SRR10902309_46767_ka_f_48.102_L_46767_", levels_back = 2)
rhino_tree <- ggtree(rhino_subset)
rhino_tree$data$label <- gsub("Edges.*", "", rhino_tree$data$label)
rhino_tree$data$label <- gsub("ka_f.*", "", rhino_tree$data$label)

rhino_tree$data$status <- ifelse(grepl("[E|S|D][R][R]", rhino_tree$data$label), "novel", "ncbi")
# rhino_tree$data <- left_join(rhino_tree$data, all_genus_curated_thin, by = c("label" = "V1_2"))
# rhino_tree$data$x10 <- iconv(rhino_tree$data$x10, from = "UTF-8", to = "ASCII", sub = "")

rhino1 <- rhino_tree + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.003)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(aes(subset = (node != "10")),color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans") +
  geom_tiplab(aes(subset = (node == "10")), color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans", fontface = 'bold', nudge_x = 0.0001) +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  theme(legend.position= "none", text = element_text(family = "Noto Sans"))

rhino1

rhino2 <- gheatmap(rhino1, addtl_ano_tree_broad, width=0.05, font.size=0,  offset = 1,  color = NA) +
  scale_fill_manual(breaks=addtl_ano_tree_broad_vals$broader_gen, 
                    values=addtl_ano_tree_broad_vals$color_values, name="genotype") + 
  theme(legend.position = 'right', text = element_text(family = "Noto Sans"))

ggsave("rhino2_subtree.jpg", rhino2, width = 25, height =10, units = "cm", limitsize = F)

pae_data <- fromJSON("fold_rhino_srr10902309_46767_full_data_0.json")
pae_matrix <- pae_data$pae


pae_df <- as.data.frame(pae_matrix)
pae_df$Residue_i <- seq_len(nrow(pae_df))
pae_df <- pivot_longer(pae_df, cols = -Residue_i, names_to = "Residue_j", values_to = "PAE_Error")
pae_df$Residue_j <- as.numeric(gsub("V", "", pae_df$Residue_j))
colnames(pae_df) <- c("Residue_i", "Residue_j", "PAE_Error")

pae <- ggplot(data = pae_df, aes(x = Residue_i, y = Residue_j, fill = PAE_Error)) +
  geom_tile() +
  scale_fill_gradient(low = "darkgreen", high = "white", limit = c(0, 31.7)) + # Emulate AlphaFold's color scheme
  labs(title = "Predicted Aligned Error (PAE) Plot",
       x = "Residue",
       y = "Residue",
       fill = "PAE Error (Å)") +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 10)) +
  coord_fixed() + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0)) # Ensures the tiles are square

ggsave("rhino_pae.png", pae, width = 10, height = 10, units = "cm", limitsize = F)

world_map <- map_data("world")


library(raster)
library(geodata)
library(tidyterra)

bwa_outer <- gadm("BWA", level = 0, path = tempdir())

bwa_outer <- fortify(bwa_outer)

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80") +
  geom_sf(data = bwa_outer, fill = "grey60", color = "grey70", linewidth = 1) +
  coord_sf(xlim = c(-10, 60), ylim = c(20, -40), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70")) 

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80") +
  geom_sf(data = bwa_outer, fill = "grey60", color = "grey70", linewidth = 1) +
  # coord_sf(xlim = c(-10, 60), ylim = c(20, -40), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70"))
  
  


##################################################ZARD################################

zard_contig <- read.table("SRR22028468_199653_zard_genome.txt", sep ="\t", header = T)

#set consistent colors for each gene
color_match_genes = setNames(c("#46D452","#99A599","#99A599", "#99A599", "#CBD7CB", "#CBD7CB", "#CBD7CB",  "#CBD7CB", "#B2F0B7"), c("L1", "L2", "E1", "E2", "E4", "E5", "E6", "E7", "JR"))

#SRR25256522_663139	JR	3617	2675	-http://127.0.0.1:40769/graphics/plot_zoom_png?width=1810&height=920

ggplot(zard_contig, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = T)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
  geom_gene_arrow(aes(xmin = 5441+69, xmax = 5441+1301, y = "SRR22028468_199653"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  # geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_text(aes(x = end - ((end-start)/2), y = 1.2, label = gene, fontface = 'bold', family = "Noto Sans"))  +  scale_fill_manual(values = color_match_genes) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes() + theme(text = element_text(family="Noto Sans", size = 10))

ggsave("zard_map.png", width = 20, height = 10, units = "cm", limitsize = F)

zard_subset <- tree_subset(final_tree, "SRR22028468_199653_ka_f_18.451_L_19965", levels_back = 1)
zard_tree <- ggtree(zard_subset)
zard_tree$data$label <- gsub("Edges.*", "", zard_tree$data$label)
zard_tree$data$label <- gsub("ka_f.*", "", zard_tree$data$label)

zard_tree$data$status <- ifelse(grepl("[E|S|D][R][R]", zard_tree$data$label), "novel", "ncbi")
# zard_tree$data <- left_join(zard_tree$data, all_genus_curated_thin, by = c("label" = "V1_2"))
# zard_tree$data$x10 <- iconv(zard_tree$data$x10, from = "UTF-8", to = "ASCII", sub = "")

zard1 <- zard_tree + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.003)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(aes(subset = (node != "11")),color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans") +
  geom_tiplab(aes(subset = (node == "11")), color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans", fontface = 'bold', nudge_x = 0.0001) +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  theme(legend.position= "none", text = element_text(family = "Noto Sans"))

zard1

zard2 <- gheatmap(zard1, addtl_ano_tree_broad, width=0.05, font.size=0,  offset = 1,  color = NA) +
  scale_fill_manual(breaks=addtl_ano_tree_broad_vals$broader_gen, 
                    values=addtl_ano_tree_broad_vals$color_values, name="genotype") + theme(legend.position = 'right', text = element_text(family = "Noto Sans"))

ggsave("zard2_subtree.jpg", zard2, width = 25, height =10, units = "cm", limitsize = F)

pae_data <- fromJSON("fold_legless_zard_srr22028468_199653_full_data_0.json")
pae_matrix <- pae_data$pae


pae_df <- as.data.frame(pae_matrix)
pae_df$Residue_i <- seq_len(nrow(pae_df))
pae_df <- pivot_longer(pae_df, cols = -Residue_i, names_to = "Residue_j", values_to = "PAE_Error")
pae_df$Residue_j <- as.numeric(gsub("V", "", pae_df$Residue_j))
colnames(pae_df) <- c("Residue_i", "Residue_j", "PAE_Error")

pae <- ggplot(data = pae_df, aes(x = Residue_i, y = Residue_j, fill = PAE_Error)) +
  geom_tile() +
  scale_fill_gradient(low = "darkgreen", high = "white", limit = c(0, 31.7)) + # Emulate AlphaFold's color scheme
  labs(title = "Predicted Aligned Error (PAE) Plot",
       x = "Residue",
       y = "Residue",
       fill = "PAE Error (Å)") +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 10)) +
  coord_fixed() + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0)) # Ensures the tiles are square

ggsave("zard_pae.png", pae, width = 10, height = 10, units = "cm", limitsize = F)

vnm_outer <- gadm("VNM", level = 0, path = tempdir())
vnm_inner <- gadm("VNM", level = 1, path = tempdir())


vnm_outer <- fortify(vnm_outer)
vnm_inner <- fortify(vnm_inner)


ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80") +
  geom_sf(data = vnm_outer, fill = "#B5C096", color = "grey70", linewidth = 0.5) +
  geom_sf(data = vnm_inner, fill = NA, color = "grey60", linewidth = 0.5) +
  geom_point(aes(x=107.49, y=10.50), shape = 23, size = 3, fill = "#53741D", color = "grey60", stroke = 1) +
  coord_sf(xlim = c(90, 120), ylim = c(5, 30), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70")) 

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80", linewidth = 0.5) +
  geom_sf(data = vnm_outer, fill = "#53741D", color = "grey70", linewidth = 0)+
  # coord_sf(xlim = c(-10, 60), ylim = c(20, -40), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70"))


##################################################HUMAN################################

human_contig <- read.table("SRR13789839_2669_human_genome.txt", sep ="\t", header = T)

#set consistent colors for each gene
color_match_genes = setNames(c("#46D452","#99A599","#99A599", "#99A599", "#CBD7CB", "#CBD7CB", "#CBD7CB",  "#CBD7CB", "#B2F0B7"), c("L1", "L2", "E1", "E2", "E4", "E5", "E6", "E7", "JR"))

#SRR25256522_663139	JR	3617	2675	-

ggplot(human_contig, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = T)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
  geom_gene_arrow(aes(xmin = 2779+138, xmax = 2779+1343, y = "SRR13789839_2669"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  # geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_text(aes(x = end - ((end-start)/2), y = 1.2, label = gene, fontface = 'bold', family = "Noto Sans"))  +  scale_fill_manual(values = color_match_genes) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes() + theme(text = element_text(family="Noto Sans", size = 10)) + scale_x_reverse()

ggsave("human_map.png", width = 20, height = 10, units = "cm", limitsize = F)

usa_outer <- gadm("USA", level = 0, path = tempdir())
usa_conn <- gadm("USA", level = 1, path = tempdir())


usa_outer <- fortify(usa_outer)
usa_conn <- fortify(usa_conn)


ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80") +
  geom_sf(data = usa_outer, fill = "#E2C4C4", color = "grey70") +
  # geom_point(aes(x=107.49, y=10.50), shape = 23, size = 5, fill = "#53741D", color = "grey60", stroke = 1) +
  # coord_sf(xlim = c(90, 130), ylim = c(-10, 30), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70")) 

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80", linewidth = 0.5) +
  geom_sf(data = usa_outer, fill = "#E2C4C4", color = "grey70") +
  geom_sf(data = usa_conn, fill = NA, color = "grey70") +
  geom_point(aes(y=41.7, x=-72.79), shape = 21, size = 5, fill = "#D44746", color = "grey60", stroke = 1) +
  geom_point(aes(y=41.7, x=-72.79), shape = 21, size = 2, fill = "#E2C4C4", color = "grey60", stroke = 1) +
   coord_sf(xlim = c(-130, -60), ylim = c(15, 60), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70"))



human_subset <- tree_subset(final_tree, "SRR13789839_2669_ka_f_16.185___138-134", levels_back = 4)
human_tree <- ggtree(human_subset)
human_tree$data$label <- gsub("Edges.*", "", human_tree$data$label)
human_tree$data$label <- gsub("ka_f.*", "", human_tree$data$label)

human_tree$data$status <- ifelse(grepl("[E|S|D][R][R]", human_tree$data$label), "novel", "ncbi")
# human_tree$data <- left_join(human_tree$data, all_genus_curated_thin, by = c("label" = "V1_2"))
# human_tree$data$x10 <- iconv(human_tree$data$x10, from = "UTF-8", to = "ASCII", sub = "")

human1 <- human_tree + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.003)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(aes(subset = (node != "8")),color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans") +
  geom_tiplab(aes(subset = (node == "8")), color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans", fontface = 'bold', nudge_x = 0.0001) +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  theme(legend.position= "none", text = element_text(family = "Noto Sans"))

human1 

human2 <- gheatmap(human1, addtl_ano_tree_broad, width=0.05, font.size=0,  offset = 0.5,  color = NA) +
  scale_fill_manual(breaks=addtl_ano_tree_broad$broader_gen, 
                     values=addtl_ano_tree_broad_vals$color_values, name="genotype") + theme(legend.position = 'right', text = element_text(family = "Noto Sans"))

ggsave("human2_subtree.jpg", human2, width = 25, height =10, units = "cm", limitsize = F)


pae_data <- fromJSON("fold_human_srr13789839_2669_full_data_0.json")
pae_matrix <- pae_data$pae


pae_df <- as.data.frame(pae_matrix)
pae_df$Residue_i <- seq_len(nrow(pae_df))
pae_df <- pivot_longer(pae_df, cols = -Residue_i, names_to = "Residue_j", values_to = "PAE_Error")
pae_df$Residue_j <- as.numeric(gsub("V", "", pae_df$Residue_j))
colnames(pae_df) <- c("Residue_i", "Residue_j", "PAE_Error")

pae <- ggplot(data = pae_df, aes(x = Residue_i, y = Residue_j, fill = PAE_Error)) +
  geom_tile() +
  scale_fill_gradient(low = "darkgreen", high = "white", limit = c(0, 31.7)) + # Emulate AlphaFold's color scheme
  labs(title = "Predicted Aligned Error (PAE) Plot",
       x = "Residue",
       y = "Residue",
       fill = "PAE Error (Å)") +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 10)) +
  coord_fixed() + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0)) # Ensures the tiles are square

ggsave("human_pae.png", pae, width = 10, height = 10, units = "cm", limitsize = F)


##################################################SALMON###########################################


salmon_contig <- read.table("SRR20078264_4021_salmon_genome.txt", sep ="\t", header = T)

#set consistent colors for each gene
color_match_genes = setNames(c("#46D452","#99A599","#99A599", "#99A599", "#CBD7CB", "#CBD7CB", "#CBD7CB",  "#CBD7CB", "#B2F0B7"), c("L1", "L2", "E1", "E2", "E4", "E5", "E6", "E7", "JR"))

#SRR25256522_663139	JR	3617	2675	-

ggplot(salmon_contig, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = T)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
  geom_gene_arrow(aes(xmin = 2870+66, xmax = 2870+1346, y = "SRR20078264_4021"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  geom_gene_arrow(aes(xmin = 1504, xmax = 2901, y = "SRR20078264_4021"), arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white",fill = "grey30", alpha = 0.3) +
  # geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_text(aes(x = end - ((end-start)/2), y = 1.2, label = gene, fontface = 'bold', family = "Noto Sans"))  +  scale_fill_manual(values = color_match_genes) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes() + theme(text = element_text(family="Noto Sans", size = 10))

ggsave("salmon_map.png", width = 20, height = 10, units = "cm", limitsize = F)

usa_outer <- gadm("USA", level = 0, path = tempdir())
usa_conn <- gadm("USA", level = 1, path = tempdir())


usa_outer <- fortify(usa_outer)
usa_conn <- fortify(usa_conn)


ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80") +
  geom_sf(data = usa_outer, fill = "#BBD1D7", color = "grey70") +
  # geom_point(aes(x=107.49, y=10.50), shape = 23, size = 5, fill = "#53741D", color = "grey60", stroke = 1) +
  # coord_sf(xlim = c(90, 130), ylim = c(-10, 30), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70")) 

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 23.5, linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -23.5, linetype = "dashed", color = "grey90") +
  
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region), color = "grey70", 
           fill = "grey80", linewidth = 0.5) +
  geom_sf(data = usa_outer, fill = "#E2EBEE", color = "grey70") +
  geom_sf(data = usa_conn, fill = "#E2EBEE", color = "grey70") +
  geom_moon(aes(y=45.3535293, x=-116.4035803, ratio = 0.2, right = F), size = 7, fill = "#A3C3D5", color = "grey60", stroke = 0.4) +
  coord_sf(xlim = c(-130, -60), ylim = c(15, 60), expand = T) +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 20), 
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 2),  
        axis.line.x = element_line(color = "grey70", linewidth = 1),
        axis.line.y = element_line(color = "grey70", linewidth = 1),
        axis.text.y=element_text(colour="grey70"),
        axis.text.x=element_text(colour="grey70"),
        axis.ticks = element_line(colour="grey70"))



salmon_subset <- tree_subset(final_tree, "SRR20078264_4021_ka_f_10.293_L_4021_L_", levels_back = 2)
salmon_tree <- ggtree(salmon_subset)
salmon_tree$data$label <- gsub("Edges.*", "", salmon_tree$data$label)
salmon_tree$data$label <- gsub("ka_f.*", "", salmon_tree$data$label)

salmon_tree$data$status <- ifelse(grepl("[E|S|D][R][R]", salmon_tree$data$label), "novel", "ncbi")
# salmon_tree$data <- left_join(salmon_tree$data, all_genus_curated_thin, by = c("label" = "V1_2"))
# salmon_tree$data$x10 <- iconv(salmon_tree$data$x10, from = "UTF-8", to = "ASCII", sub = "")

salmon1 <- salmon_tree + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.003)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(aes(subset = (node != "8")),color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans") +
  geom_tiplab(aes(subset = (node == "8")), color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans", fontface = 'bold', nudge_x = 0.0001) +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  theme(legend.position= "none", text = element_text(family = "Noto Sans"))

salmon1 

salmon2 <- gheatmap(salmon1, addtl_ano_tree_broad, width=0.05, font.size=0,  offset = 1,  color = NA) +
  scale_fill_manual(breaks=addtl_ano_tree_broad$broader_gen, 
                    values=addtl_ano_tree_broad_vals$color_values, name="genotype") + theme(legend.position = 'right', text = element_text(family = "Noto Sans"))

ggsave("salmon2_subtree.jpg", salmon2, width = 25, height =10, units = "cm", limitsize = F)

pae_data <- fromJSON("fold_salmon_srr20078264_4021_full_data_0.json")
pae_matrix <- pae_data$pae


pae_df <- as.data.frame(pae_matrix)
pae_df$Residue_i <- seq_len(nrow(pae_df))
pae_df <- pivot_longer(pae_df, cols = -Residue_i, names_to = "Residue_j", values_to = "PAE_Error")
pae_df$Residue_j <- as.numeric(gsub("V", "", pae_df$Residue_j))
colnames(pae_df) <- c("Residue_i", "Residue_j", "PAE_Error")

pae <- ggplot(data = pae_df, aes(x = Residue_i, y = Residue_j, fill = PAE_Error)) +
  geom_tile() +
  scale_fill_gradient(low = "darkgreen", high = "white", limit = c(0, 31.7)) + # Emulate AlphaFold's color scheme
  labs(title = "Predicted Aligned Error (PAE) Plot",
       x = "Residue",
       y = "Residue",
       fill = "PAE Error (Å)") +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 10)) +
  coord_fixed() + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0)) # Ensures the tiles are square

ggsave("salmon_pae.png", pae, width = 10, height = 10, units = "cm", limitsize = F)


##################################################END CASE STUDIES###########################################




#alphafold matrix
# install.packages("jsonlite")
library(jsonlite)

pae_data <- fromJSON("fold_pango_srr25256522_663139_full_data_0.json")
pae_matrix <- pae_data$pae


pae_df <- as.data.frame(pae_matrix)
pae_df$Residue_i <- seq_len(nrow(pae_df))
pae_df <- pivot_longer(pae_df, cols = -Residue_i, names_to = "Residue_j", values_to = "PAE_Error")
pae_df$Residue_j <- as.numeric(gsub("V", "", pae_df$Residue_j))
colnames(pae_df) <- c("Residue_i", "Residue_j", "PAE_Error")

pae <- ggplot(data = pae_df, aes(x = Residue_i, y = Residue_j, fill = PAE_Error)) +
  geom_tile() +
  scale_fill_gradient(low = "darkgreen", high = "white", limit = c(0, 31.7)) + # Emulate AlphaFold's color scheme
  labs(title = "Predicted Aligned Error (PAE) Plot",
       x = "Residue",
       y = "Residue",
       fill = "PAE Error (Å)") +
  theme_classic() +
  theme(text = element_text(family="Noto Sans", size = 10)) +
  coord_fixed() + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0)) # Ensures the tiles are square

ggsave("pango_pae.png", pae, width = 10, height = 10, units = "cm", limitsize = F)







  #generate generic map for linearized PV genome

generic_pv <- read.table("generic_pv.txt", sep ="\t", header = T)


#set consistent colors for each gene
color_match_genes = setNames(c("#46D452","#99A599","#99A599", "#99A599", "#CBD7CB", "#CBD7CB", "#CBD7CB",  "#CBD7CB", "#B2F0B7"), c("L1", "L2", "E1", "E2", "E4", "E5", "E6", "E7", "JR"))

#SRR25256522_663139	JR	3617	2675	-http://127.0.0.1:40769/graphics/plot_zoom_png?width=1810&height=920

ggplot(generic_pv, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = T)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm"), colour = "white", alpha = 0.8) +
  geom_gene_arrow(aes(xmin = 5441+69, xmax = 5441+1301, y = "PV"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  # geom_gene_arrow(aes(xmin = 3617, xmax = 2675, y = "SRR25256522_663139"), color="white", fill = "#B2F0B7", arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(5, "mm"),  alpha=0.4) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_text(aes(x = end - ((end-start)/2), y = 1.2, label = gene, fontface = 'bold', family = "Noto Sans"))  +  scale_fill_manual(values = color_match_genes) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes() + theme(text = element_text(family="Noto Sans", size = 10))

ggsave("generic_map.png", width = 20, height = 10, units = "cm", limitsize = F)


