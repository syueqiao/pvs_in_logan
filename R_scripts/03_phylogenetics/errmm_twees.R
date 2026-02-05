library(tidyverse)
library(Polychrome)
library(ggtree)
library(ggplot2)
library(magrittr)
library(ggnewscale)
library(treeio)
library(ape)


is.waive <- function(x) {
  return(inherits(x, "waiver"))
}



#' Encapsulate path
#'
#' This function is used to encapsulate paths to allow for behaviour switches.
#' One use is for example when plotting. The plot_all method will encapsulate a
#' path so that plots may be saved to a directory structure. Other plot methods,
#' e.g. plot_model_performance do not encapsulate a path, and if the user calls
#' these functions directly, the plot may be written to the provided path
#' instead of a directory structure.
#'
#' @return encapsulated_path object
#' @md
#' @keywords internal
encapsulate_path <- function(path) {
  structure(path, class = "encapsulated_path")
}



#to visualize the tree properly
updated_tree3 <- read.tree("tree3.aln.treefile")
#bootstrap value is saved as "label" in the tree 
tree <- ggtree(updated_tree3) + geom_tiplab() + geom_text(aes(label=label), hjust=-.3)

#to get the species labels for the tips, split into pave and SRA samples
tip_labs_ncbi_novel <- as.data.frame(updated_tree3[["tip.label"]])
#1039 entries,
colnames(tip_labs_ncbi_novel) <- c("ncbi_and_novel_tip")
SRA_only <- subset(ncbi_tags_manual_and_sra, grepl("[SED][R]{2}", ncbi_tags_manual_and_sra$Newick_label))

SRA_only <- as.data.frame(grepl("[SED][R]{2}", ncbi_tags_manual_and_sra$Newick_label))
length(SRA_only$`grep("[SED][R]{2}", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T)`)

ncbi_tags_manual_and_sra$Newick_label <- sub("-", "\\.", ncbi_tags_manual_and_sra$Newick_label)

p <- ggtree(updated_tree3, layout = "rectangular") %<+% SRA_only
p$data$status <- grepl("[SED][R]{2}", p$data$label)

# pp <- ggtree(p, aes(color = status), layout = "circular") %<+% SRA_only

# P50 = createPalette(251,  c("#ff0000", "#00ff00", "#0000ff"))
p + theme(legend.position= "none") 

P50 = createPalette(305,  c("#ff0000", "#00ff00", "#0000ff"))

names(P50) <- NULL

#v2 of graph with colored branches and various other values
pp <- p + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.035)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(color = "black",size = 2, hjust = 0, align= T, linetype = "blank") +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  theme(legend.position= "none")

pp
# known_ncbi_and_novel <- as.data.frame(grep("ncbi_and_novelase.|MGY", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
# known_ncbi_and_novel$species <- c("known")
# colnames(known_ncbi_and_novel) <- c('Newick_label', 'status') 

ppp <- pp %<+% SRA_merge_info
pppp <- ppp + geom_tippoint() + geom_text(aes(label = host), size = 4, hjust = -0.05, color = "black", fontface = "bold") +  geom_tiplab(size = 2, hjust = -1.5, color = "white", alpha = 0) 

ppppp <- pppp %<+% ncbi_tags_manual
pppppp <- ppppp + geom_tippoint() + geom_text(aes(label = uhhh), size = 4, hjust = -0.05, color = "black") + geom_tiplab(size = 2, hjust = -2, color = "black")
# 
ggsave("square+color.pdf", width = 80, height = 200, units = "cm", limitsize = F)

#
# #deal with SRA stuff, separate out for a clean input list
# SRA_ncbi_and_novel <- as.data.frame(grep("[SED][R]{2}", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
# colnames(SRA_ncbi_and_novel) <- c('ncbi_and_novel_SRA')
# SRA_ncbi_and_novel$ncbi_and_novel_SRA_info <- SRA_ncbi_and_novel$ncbi_and_novel_SRA
# SRA_ncbi_and_novel <- SRA_ncbi_and_novel %>% separate(ncbi_and_novel_SRA, into = c("library", "info"), sep="_")
# #grep with this
# write.table(SRA_ncbi_and_novel$library, "ncbi_and_novel_SRA.txt", sep = "/t", col.names = F, row.names = F, quote = F)

head(ncbi_tags_manual_and_sra)
tail(tip_labs_ncbi_novel)

#generate a dataframe that can be joined back into the data set in the tree to add annotation information
#using the data df as a guide
#blast outputn to assign predicted genus
blast_pap_tax <- read.table("pap_genus_ncbi_virus.tsv", sep = "\t", header = F)
#tax id to lineage file
tax_id_line <- read.table("tax_lineages.txt", header = F, fill = T, sep = "|")

blast_pap_tax_line <- left_join(blast_pap_tax, tax_id_line, by = c("V13" = "V1"))
#filter for high confidence hits
blast_pap_tax_line_fil <- filter(blast_pap_tax_line, V11 < 0.00001)
blast_pap_tax_line_fil_sep <- separate(data = blast_pap_tax_line_fil, col = V3.y, into = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14"), sep = ";")

#group and count assignments
table(blast_pap_tax_line_fil_sep$x10)
blast_pap_tax_line_fil_sep$V1_2 <- gsub("OriginScaffoldPath=_1.*", "", blast_pap_tax_line_fil_sep$V1)
wee <- blast_pap_tax_line_fil_sep %>% group_by(V1_2) %>% count(x10)
wee2 <- wee %>% group_by(V1_2) %>% slice_max(n)
#get "confident" ones
genus_conf_more5 <- filter(wee2, n > 5)
genus_conf_more5$x10 <- sub("\t", "unk", genus_conf_more5$x10)
write.table(wee2, "all_genus.tsv", sep = "\t", col.names = T, row.names = F)

all_genus_curated <- read.table("all_genus_curated.tsv.txt", sep = "\t", header = T)
all_genus_curated$x10 <- iconv(all_genus_curated$x10, from = "UTF-8", to = "ASCII", sub = "")

all_genus_curated$V1_2 <- gsub("=", "_", all_genus_curated$V1_2)
all_genus_curated$V1_2 <- gsub("Edges.*", "", all_genus_curated$V1_2)
all_genus_curated$V1_2 <- gsub("ka_f.*", "", all_genus_curated$V1_2)

all_genus_curated_thin <- select(all_genus_curated, V1_2, x10)


final_tree <- read.tree("final_pr_ncbi_and_novel_sort_trim.aln.treefile")
tree <- ggtree(final_tree)
tree$data$label <- gsub("Edges.*", "", tree$data$label)
tree$data$label <- gsub("ka_f.*", "", tree$data$label)


tree$data <- left_join(tree$data, all_genus_curated_thin, by = c("label" = "V1_2"))
tree$data$x10 <- iconv(tree$data$x10, from = "UTF-8", to = "ASCII", sub = "")
tree$data[1092, 10] <- "Zhangixalus_dugritei_associated"

tree$data$status <- ifelse(grepl("[E|S|D][R][R]", tree$data$label), "novel", "ncbi")

#for some reason need to set the is.waive function manually now???????
pp <- tree + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.035)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(color = "black",size = 2, hjust = 0, align= T, linetype = "blank") +
  # geom_text(aes(label = node), hjust = -0.5) +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  theme(legend.position= "none")

genus <- as.data.frame(tree$data[grepl("TRUE",tree$data$isTip),])
genus <- select(genus, label, x10)
rownames(genus) <- genus$label 
genus <- select(genus, x10)

g_color = createPalette(length(genus$x10),  c("#ff0000", "#00ff00", "#0000ff"))
names(g_color) <- NULL


gheatmap(pp, weeeeee, width=0.05, font.size=3, 
         colnames_angle=-45,  offset = 0.4,  color = NA) +
  scale_fill_manual(breaks=weeeeee$x10, 
                    values=g_color, name="genotype") + theme(legend.position = 'none' ) +
  geom_tiplab(aes(label=x10), size = 2, color = "grey66",  align= T, linetype = "blank", nudge_x = 1)





ggsave("ughhhhhhhh2.png", width = 80, height = 200, units = "cm", limitsize = F)
write.table(pp$data, "all_genus_addtl_anno.tsv", sep = "\t", col.names = T, row.names = F)

sra_metadata <- read.table("lib_source_big.txt", sep = "\t", fill = T, header = F)
sra_metadata["144301", "V2"] <- "GENOMIC"
sra_metadata <- unique(sra_metadata)

addtl_ano <- read.table("2025.10.29all_genus_addtl_anno.tsv.txt", sep = "\t", header = T)
addtl_ano$lib <- gsub("_.*", "", addtl_ano$label)
addtl_ano$lib <- gsub("Score.*", "", addtl_ano$lib)

addtl_ano <- left_join(addtl_ano, sra_metadata, by = c("lib" = "V1"))
addtl_tip <- filter(addtl_ano, isTip == "TRUE")
addtl_tip$broader_gen <- gsub("Cattle", "Bovine", addtl_tip$broader_gen)

addtl_tip[is.na(addtl_tip$broader_gen)] <- "Other"

addtl_tip$broader_gen <- factor(addtl_tip$broader_gen, levels = c("Human", "Bat", "Bovine",
"Canine", "Ray-finned Fish",
"Rodent", "Avian", 
"Cetacean", "Amphibian", 
"Cervine", "Pangolin", 
"Feline", "Non-human Primate", 
"Equine", "Pinnipeds", 
"Lizard", "Porcine", "Other"))

addtl_tip[is.na(addtl_tip)] <- "Other"


summary_table <- addtl_tip %>%
  group_by(generalization, status) %>%
  summarize(
    Count = n() # n() counts the number of rows in each group
  )

df_wide <- pivot_wider(
  summary_table,
  names_from = status,
  values_from = Count
)

df_wide[is.na(df_wide)] <- 0

                         
write.table(df_wide, "all_organism_group_assignments.tsv", sep = '\t', quote = F, row.names = F, col.names = F)

ggplot(addtl_tip, aes(x=broader_gen, fill = fct_rev(status))) + 
  geom_bar(stat="count") + scale_fill_manual(values = c("#4646d4", "grey")) + theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  labs(title = "Known and Novel PVs\nAssociated with Generalized Host Categories",
       x = "Category",
       y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(text = element_text(family="Noto Sans", size = 30))
 
ggsave("2026.01.27.host_group_chart.png", width = 30, height = 20, units = "cm", limitsize = F)


ggplot(filter(addtl_tip, broader_gen == "Other"), aes(x=fct_infreq(generalization) , fill =  fct_rev(status))) + 
  geom_bar(stat="count") + scale_fill_manual(values = c("#4646d4", "grey")) + theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  labs(title = "Known and Novel PVs Associated with Generalized Host Categories",
       x = "Category",
       y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(text = element_text(family="Noto Sans", size = 12))
       

ggsave("2026.01.21.host_group_chart_ZOOM.png", width = 15, height = 10, units = "cm", limitsize = F)

addtl_tip_novel <- subset(addtl_tip, grepl("[SED][R]{2}", addtl_tip$label))

addtl_tip_novel_tb <- as.data.frame(table(addtl_tip_novel$broader_gen))

addtl_tip_novel_tb_ord <- addtl_tip_novel_tb %>%
  mutate(Var1 = forcats::fct_reorder(factor(Var1), Freq))

levels(factor(addtl_tip_novel_tb_ord$Var1)) 


ggplot(addtl_tip_novel_tb_ord, aes(x=Var1, y=Freq)) +
  geom_segment( aes(x=Var1, xend=Var1, y=0, yend=Freq), color="grey20") +
  geom_point( color="#4646d4", size=2, alpha=0.9) +
  theme_bw() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+
  labs(title = "Novel Papillomaviruses by Associated Species",
       x = "Associated Species",
       y = "Count")

ggsave("uasdasfasfasasfasdasfasaf.png", width = 15, height = 25, units = "cm", limitsize = F)

ggplot(addtl_tip, aes(x=specific_species, fill = status)) + 
  geom_bar(stat="count") + scale_fill_manual(values = c("grey", "#4646d4")) + theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  labs(title = "Known and Novel PVs\nAssociated with Specific Species",
       x = "Category",
       y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(text = element_text(family="Noto Sans", size = 20))

novel_in_novel_hosts <- addtl_tip %>% 
  group_by(specific_species) %>% 
  mutate(status_count = !("ncbi" %in% status)) %>% filter(status_count == TRUE) %>% filter(grepl("[SED][R]{2}",label))

write.table(novel_in_novel_hosts, "novel_in_novel_hosts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)


novel_headers <- read.table("final_novel_pvs.headers", fill = T, header = F)
novel_headers$V1 <- gsub(">", "", novel_headers$V1)
novel_headers$V1 <- gsub("_.*", "", novel_headers$V1)
novel_headers$V1 <- gsub("Score.*", "", novel_headers$V1)

novel_headers_type <- left_join(novel_headers, sra_metadata, by = c("V1" = "V1"))
write.table(novel_headers_type$V1, "all_novel_libraries.tsv", sep = "\t", col.names = F, row.names = F, quote = F)


type_table <- as.data.frame(table(novel_headers_type$LibrarySource))
type_table = type_table[-1,]

#genomic = 146, METAGENOMIC = 121, METATRANSCRIPTOMIC = 2, OTHER = 1, VIRAL RNA = 36, TRANSCRIPTOMIC = 22

pie(c(146, 121, 36, 22, 2, 1), labels = c("GENOMIC", "METAGENOMIC", "VIRAL RNA", "TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "OTHER"), col = c("#C55B25", "#F9CCB5", "#668C6E", "#59D86A", "#C8F2CD", "grey40"))

ggplot(type_table, aes(x=reorder(Var1, Freq), y=Freq, fill = Var1)) + 
  geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = c("#101048", "#1C1C7E", "#E3E3F9", "gray", "#D2D2F6", "#B9B9F1")) +
  coord_flip()+ 
  theme(text = element_text(family = "Noto Sans", size = 28))+ 
  theme(axis.text.y = element_text(hjust=1))

ggsave("breakdown.jpg", width = 35, height = 20, units = "cm", limitsize = F)

all_pvs_mapping$library <- gsub("_.*", "", all_pvs_mapping$V1.x)
all_pvs_mapping$library <- gsub("Score.*", "", all_pvs_mapping$library)


novel_headers_type <- left_join(all_pvs_mapping, sra_metadata, by = c("library" = "V1"))
novel_headers_type$collapsed_type <-coalesce(novel_headers_type$LibrarySource, novel_headers_type$V2)

novel_headers_type_missing_type <- 
  novel_headers_type %>% 
  filter(is.na(collapsed_type))

novel_headers_type_missing_type$library <- gsub("^", "\"", novel_headers_type_missing_type$library)
novel_headers_type_missing_type$library <- gsub("$", "\"", novel_headers_type_missing_type$library)



paste0(novel_headers_type_missing_type$library, collapse=",")

write.table(paste0(novel_headers_type_missing_type$library, collapse=","), "novel_headers_type_missing_type.txt", sep = ",", col.names = F, row.names = F, quote = T)

addtl_annotations_for_missing <- read.table("library_annotations_for_missing.csv", sep = ",")
novel_headers_type <- left_join(novel_headers_type, addtl_annotations_for_missing, by = c("library" = "V1"))

novel_headers_type$collapsed_type_2 <- coalesce(novel_headers_type$collapsed_type, novel_headers_type$V2.y.y)

write.table(novel_headers_type, "type_missing_manual.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
novel_headers_type_MANUAL <- read.table("2025.11.28.type_missing_manual.tsv.txt", sep = "\t", header = T)
novel_headers_type_MANUAL_table <- as.data.frame(table(novel_headers_type_MANUAL$b.8, novel_headers_type_MANUAL$a.7))

library(extrafont)

ggplot(novel_headers_type_MANUAL_table, aes(x=reorder(Var1, Freq), y=Freq, fill = fct_rev(Var2))) + 
  geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = c("#4646d4", "grey")) +
  coord_flip()+ 
  theme(text = element_text(family = "Noto Sans", size = 25))+ 
  theme(axis.text.y = element_text(hjust=1))

ggsave("breakdown_new_small.jpg", width = 30, height = 15, units = "cm", limitsize = F)


  biosamp_data <- read.table("all_runs.csv", sep = ",", header = T, fill = T)
biosamp_novel_headers <- left_join(novel_headers_type_MANUAL, biosamp_data, by = c("a" = "Run"))
write.table(biosamp_novel_headers$BioSample, "all_novel_biosamples.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

addtl_ano <- left_join(addtl_ano, sra_metadata, by = c("lib" = "V1"))
addtl_ano$generalization <- iconv(addtl_ano$generalization, from = "UTF-8", to = "ASCII", sub = "")
addtl_ano <- filter(addtl_ano, isTip == "TRUE")


addtl_ano_tree_broad <- select(addtl_tip, label, broader_gen)

rownames(addtl_ano_tree_broad) <- addtl_ano_tree_broad$label
addtl_ano_tree_broad <- select(addtl_ano_tree_broad, broader_gen)


c19 <- c(
  "#D44746", "#F0B2B2", # red
  "#D4B046",
  "#175F67", # purple
  "#8C46D4", # orange
  "#E146D4",
  "#882020", "#D28036", # lt pink
  "#44C8D6", # lt purple
  "#A3E4EB", # lt orange
  "#941888", "#F0E0B2", "#C099E7", "#3092D4",
  "#53741D", "#C9E349", "#9BBED1", "#CCCCCC"
)

names(c20) <- NULL
color_key_broad <- unique(as.data.frame(addtl_ano_tree_broad$broader_gen))


ordered_gen_broad <- c("Human", "Non-human Primate", "Rodent", 
                                         "Cetacean", "Bovine", 
                                         "Porcine", "Cervine", "Pinnipeds", "Canine", "Feline", "Equine", "Pangolin", "Bat", "Avian", "Lizard", "Amphibian", "Ray-finned Fish", "Other")

color_key_broad <- as.data.frame(color_key_broad[order(sapply(color_key_broad$`addtl_ano_tree_broad$broader_gen`, function(x) which(x == ordered_gen_broad))), ])

color_key_broad$color_values <- c19
colnames(color_key_broad) <- c("gen_label", "color_values")

addtl_ano_tree_broad_vals<- left_join(addtl_ano_tree_broad, color_key_broad, by = c("broader_gen" = "gen_label"))

gheatmap <- function(p, data, offset=0, width=1, low="green", high="red", color="white",
                     colnames=TRUE, colnames_position="bottom", colnames_angle=0, colnames_level=NULL,
                     colnames_offset_x = 0, colnames_offset_y = 0, font.size=4, family="",
                     hjust=0.5, legend_title = "value", custom_column_labels = NULL) {
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## if (is.null(width)) {
  ##     width <- (p$data$x %>% range %>% diff)/30
  ## }
  
  ## convert width to width of each cell
  width <- width * (p$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)
  
  isTip <- x <- from <- to <- custom_labels <- NULL
  
  ## handle the display of heatmap on collapsed nodes
  ## https://github.com/GuangchuangYu/ggtree/issues/242
  ## extract data on leaves (& on collapsed internal nodes) 
  ## (the latter is extracted only when the input data has data on collapsed
  ## internal nodes)
  df <- p$data
  nodeCo <- intersect(df %>% dplyr::filter(is.na(x)) %>% 
                        select(.data$parent, .data$node) %>% unlist(), 
                      df %>% dplyr::filter(!is.na(x)) %>% 
                        select(.data$parent, .data$node) %>% unlist())
  labCo <- df %>% dplyr::filter(.data$node %in% nodeCo) %>% 
    select(.data$label) %>% unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- df$label %in% selCo
  
  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm=TRUE) + offset
  
  dd <- as.data.frame(data)
  ## dd$lab <- rownames(dd)
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## dd <- dd[lab, , drop=FALSE]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  
  
  dd$y <- sort(df$y)
  dd$lab <- lab
  ## dd <- melt(dd, id=c("lab", "y"))
  dd <- gather(dd, variable, value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels=colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels=colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  data_axis <- data.frame(from=dd$variable, to=V2)
  data_axis <- unique(data_axis)
  
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  if (is.null(color)) {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value), width=width, inherit.aes=FALSE)
  } else {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value), width=width, color=color, inherit.aes=FALSE)
  }
  if (is(dd$value,"numeric")) {
    p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
  } else {
    p2 <- p2 + scale_fill_discrete(na.value=NA, name = legend_title) #"white")
  }
  
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    data_axis$y <- y
    data_axis[[".panel"]] <- factor("Tree")
    # if custom column annotations are provided
    if (!is.null(custom_column_labels)) {
      # assess the type of input for the custom column annotation
      # either a vector or a named vector with positions for specific names
      if (is.null(names(custom_column_labels))) {
        if (length(custom_column_labels) > nrow(data_axis)) {
          warning(paste("Input column label vector has more elements than there are columns.",
                        "\n", "Using the first ", nrow(data_axis)," elements as labels", sep=""))
          data_axis[["custom_labels"]] <- as.character(custom_column_labels[1:nrow(data_axis)])
        } else if (length(custom_column_labels) < nrow(data_axis)) {
          warning(paste("Input column label vector has fewer elements than there are columns.",
                        "\n", "Using all available labels, n = ",
                        length(custom_column_labels), sep=""))
          data_axis[["custom_labels"]] <- as.character(c(custom_column_labels,
                                                         rep("", nrow(data_axis) - length(custom_column_labels))))
        } else {
          data_axis[["custom_labels"]] <- custom_column_labels
        }
      } else {
        if (!is.null(colnames_level)) {
          # use the colnames levels if available
          # otherwise use the default order provided by the data frame
          vector_order <- colnames_level
          
        } else {
          vector_order <- as.character(data_axis$from)
        }
        for (elem in custom_column_labels) {
          vector_order[which(vector_order == elem)] = names(which(custom_column_labels == elem))
        }
        data_axis[["custom_labels"]] <- vector_order
      }
      p2 <- p2 + geom_text(data=data_axis, aes(x=to, y = y, label=custom_labels),
                           size=font.size, family=family, inherit.aes = FALSE, angle=colnames_angle,
                           nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
    } else {
      p2 <- p2 + geom_text(data=data_axis, aes(x=to, y = y, label=from), size=font.size, family=family,
                           inherit.aes = FALSE, angle=colnames_angle,
                           nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
    }
  }
  p2 <- p2 + theme(legend.position="right")
  ## p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
  
  if (!colnames) {
    ## https://github.com/GuangchuangYu/ggtree/issues/204
    p2 <- p2 + scale_y_continuous(expand = c(0,0))
  }
  
  attr(p2, "data_axis") <- data_axis
  return(p2)
}

gheatmap(pp, addtl_ano_tree_broad, width=0.05, font.size=0,  offset = 0.4,  color = NA) +
  scale_fill_manual(breaks=addtl_ano_tree_broad_vals$broader_gen, 
                    values=addtl_ano_tree_broad_vals$color_values, name="genotype") + theme(legend.position = 'right')


ggsave("2025.11.18_fixed_errors_new_broad_colors_LARGE.pdf", width = 50, height =200, units = "cm", limitsize = F)
ggsave("2025.11.18_fixed_errors_new_broad_colors_LARGE.png", width = 50, height =200, units = "cm", limitsize = F)

###subset for pango#######

########Subset tree objects by related nodes########
#'
#' This function allows for a tree object to be subset by specifying a
#' node and returns all related nodes within a selected number of
#' levels
#'
#' @param tree a tree object of class phylo
#' @param node either a tip label or a node number for the given
#' tree that will be the focus of the subsetted tree
#' @param levels_back a number specifying how many nodes back from
#' the selected node the subsetted tree should include
#' @param group_node whether add grouping information of selected node
#' @param group_name group name (default 'group') for storing grouping information if group_node = TRUE
#' @param root_edge If TRUE (by default), set root.edge to path length of orginal root to the root of subset tree
#'
#' @details This function will take a tree and a specified node from
#' that tree and subset the tree showing all relatives back to a specified
#' number of nodes. This function allows for a combination of
#' \code{ancestor} and \code{offspring} to return a subsetted
#' tree that is of class phylo. This allows for easy graphing of the tree
#' with \code{ggtree}
#'
#' @examples
#' \dontrun{
#'   nwk <- system.file("extdata", "sample.nwk", package="treeio")
#'   tree <- read.tree(nwk)
#'
#'   sub_tree <- tree_subset(tree, node = "A", levels_back = 3)
#'   ggtree(sub_tree) + geom_tiplab() + geom_nodelab()
#' }
#'
#' @rdname tree_subset
#' @export
tree_subset <- function(tree, node, levels_back = 5, group_node = TRUE,
                        group_name = "group", root_edge = TRUE){
  UseMethod("tree_subset")
}


#' @method tree_subset phylo
#' @rdname tree_subset
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#'   nwk <- system.file("extdata", "sample.nwk", package="treeio")
#'   tree <- read.tree(nwk)
#'
#'   sub_tree <- tree_subset(tree, node = "A", levels_back = 3)
#'   ggtree(sub_tree) + geom_tiplab() + geom_nodelab()
#' }
#'
#' @importFrom utils tail
#' @importFrom utils head
#' @importFrom rlang quo .data
#' @export
tree_subset.phylo <- function(tree, node, levels_back = 5, group_node = TRUE,
                              group_name = "group", root_edge = TRUE){
  
  x <- tree_subset_internal(tree = tree, node = node, levels_back = levels_back, root_edge = root_edge)
  
  
  ## This drops all of the tips that are not included in group_nodes
  subtree <- drop.tip(tree, tree$tip.label[-x$subset_nodes], rooted = TRUE)
  
  if (group_node) subtree <- groupOTU.phylo(subtree, .node = x$group_labels, group_name = group_name)
  
  subtree$root.edge <- x$root.edge
  
  return(subtree)
}

##' @importFrom dplyr filter
##' @importFrom dplyr pull
##' @importFrom dplyr bind_rows
tree_subset_internal <- function(tree, node, levels_back = 5, root_edge = TRUE) {
  
  ## error catching to ensure the tree input is of class phylo
  ## if (class(tree) %in% c("phylo", "treedata")) {
  ##   tree_df <- tidytree::as_tibble(tree)
  ## } else {
  ##   stop("tree must be of class 'phylo'")
  ## }
  
  ## error catching to ensure the levels_back input is numeric
  ## or can be converted to numeric
  if (!is.numeric(levels_back)) {
    levels_back <- as.numeric(levels_back)
    if (is.na(levels_back)) stop("'levels_back' must be of class numeric")
  }
  
  tree_df <- tidytree::as_tibble(tree)
  
  selected_node <- node
  
  is_tip <- tree_df %>%
    dplyr::mutate(isTip = !.data$node %in% .data$parent) %>%
    dplyr::filter(.data$node == selected_node | .data$label == selected_node) %>%
    dplyr::pull(.data$isTip)
  
  if (is_tip & levels_back == 0){
    stop("The selected node (", selected_node, ") is a tip. 'levels_back' must be > 0",
         call. = FALSE)
  }
  
  if (is_tip) {
    group_labels <- tree_df %>%
      dplyr::filter(.data$node == selected_node | .data$label == selected_node) %>%
      dplyr::pull(.data$label)
  } else {
    group_labels <- tree_df %>%
      tidytree::offspring(selected_node) %>%
      dplyr::filter(!.data$node %in% .data$parent) %>%
      dplyr::pull(.data$label)
  }
  
  ## This pipeline returns the tip labels of all nodes related to
  ## the specified node
  ##
  ## The tail/head combo isolates the base node of the subsetted tree
  ## as the output from ancestor lists the closest parent nodes of a
  ## given node from the bototm up.
  ##
  ## It then finds all of the offspring of that parent node. From there
  ## it filters to include only tip and then pulls the labels.
  
  if (levels_back == 0) {
    new_root_node <- selected_node
  } else {
    new_root_node <- tidytree::ancestor(tree_df, selected_node) %>%
      tail(levels_back) %>%
      head(1) %>%
      dplyr::pull(.data$node)
  }
  
  subset_labels <- tidytree::offspring(tree_df, new_root_node) %>%
    dplyr::filter(!.data$node %in% .data$parent) %>%
    dplyr::pull(.data$label)
  
  ## This finds the nodes associated with the labels pulled
  subset_nodes <- which(tree$tip.label %in% subset_labels)
  
  root.edge <- NULL
  if (is.null(tree$edge.length)) {
    root_edge <- FALSE
    ## if not branch length info, no need to determine root.edge
  }
  
  if (root_edge) {
    root.edge <- ancestor(tree_df, new_root_node) %>%
      bind_rows(dplyr::filter(tree_df, node == new_root_node)) %>%
      pull(.data$branch.length) %>%
      sum(na.rm = TRUE)
    if (root.edge == 0)
      root.edge <- NULL
  }
  
  return(list(
    subset_nodes = subset_nodes,
    new_root_node = new_root_node,
    group_labels = group_labels,
    root.edge = root.edge
  ))
}

#' @method tree_subset treedata
#' @rdname tree_subset
#' @importFrom magrittr %>%
#'
#' @export
tree_subset.treedata <- function(tree, node, levels_back = 5, group_node = TRUE,
                                 group_name = "group", root_edge = TRUE){
  
  x <- tree_subset_internal(tree = tree@phylo, node = node, levels_back = levels_back)
  
  subtree <- drop.tip(tree, tree@phylo$tip.label[-x$subset_nodes], rooted = TRUE)
  
  if (group_node) subtree <- groupOTU(subtree, .node = x$group_labels, group_name = group_name)
  
  subtree@phylo$root.edge <- x$root.edge
  
  return(subtree)
}
###function defined
##anotha one
##' @importFrom ape which.edge
gfocus <- function(phy, focus, group_name, focus_label=NULL,
                   overlap="overwrite", connect = FALSE) {
  
  ## see https://goo.gl/VMMVhi for connect parameter
  
  overlap <- match.arg(overlap, c("origin", "overwrite", "abandon"))
  
  if (is.factor(focus)) {
    focus <- as.character(focus)
  }
  
  if (is.character(focus)) {
    focus <- which(phy$tip.label %in% focus)
  }
  
  n <- getNodeNum(phy)
  if (is.null(attr(phy, group_name))) {
    foc <- rep(0, n)
  } else {
    foc <- attr(phy, group_name)
  }
  i <- max(suppressWarnings(as.numeric(foc)), na.rm=TRUE) + 1
  if (is.null(focus_label)) {
    focus_label <- i
  }
  
  induced_edge <- phy$edge[which.edge(phy, focus),]
  
  hit <- unique(as.vector(induced_edge))
  if (overlap == "origin") {
    sn <- hit[is.na(foc[hit]) | foc[hit] == 0]
  } else if (overlap == "abandon") {
    idx <- !is.na(foc[hit]) & foc[hit] != 0
    foc[hit[idx]] <- NA
    sn <- hit[!idx]
  } else {
    sn <- hit
  }
  
  if (length(sn) > 0 && connect) {
    if (sum(table(induced_edge[,1]) > 1) == 1) {
      sn <- focus
    }
  }
  
  if (length(sn) > 0) {
    foc[sn] <- focus_label
  }
  
  attr(phy, group_name) <- foc
  phy
}

##' @method groupOTU phylo
##' @export
groupOTU.phylo <- function(.data, .node, group_name="group", ...) {
  phy <- .data
  focus <- .node
  attr(phy, group_name) <- NULL
  if ( is(focus, "list") ) {
    for (i in seq_along(focus)) {
      phy <- gfocus(phy, focus[[i]], group_name, names(focus)[i], ...)
    }
  } else {
    phy <- gfocus(phy, focus, group_name, ...)
  }
  res <- attr(phy, group_name)
  res[is.na(res)] <- 0
  attr(phy, group_name) <- factor(res)
  return(phy)
}

##' @method groupOTU treedata
##' @export
groupOTU.treedata <- function(.data, .node, group_name = "group", ...) {
  .data@phylo <- groupOTU(as.phylo(.data), .node, group_name, ...)
  return(.data)
}


##' calculate total number of nodes
##'
##'
##' @title getNodeNum
##' @param tree tree object
##' @return number
##' @export
##' @examples
##' getNodeNum(rtree(30))
##' @author Guangchuang Yu
getNodeNum <- function(tree) {
  Nnode(tree, internal.only=FALSE)
}
########
library(treeio)

pango_subset <- tree_subset(final_tree, "SRR25256522_663139_ka_f_8.130___246-14", levels_back = 4)
pango_tree <- ggtree(pango_subset)
pango_tree$data$label <- gsub("Edges.*", "", pango_tree$data$label)
pango_tree$data$label <- gsub("ka_f.*", "", pango_tree$data$label)

pango_tree$data$status <- ifelse(grepl("[E|S|D][R][R]", pango_tree$data$label), "novel", "ncbi")
pango_tree$data <- left_join(pango_tree$data, all_genus_curated_thin, by = c("label" = "V1_2"))
pango_tree$data$x10 <- iconv(pango_tree$data$x10, from = "UTF-8", to = "ASCII", sub = "")

pango1 <- pango_tree + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.003)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(aes(subset = (node != "9")),color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans") +
  geom_tiplab(aes(subset = (node == "9")), color = "black",size = 5, hjust = 0, align= T, linetype = "blank", family = "Noto Sans", fontface = 'bold', nudge_x = 0.0001) +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  theme(legend.position= "none", text = element_text(family = "Noto Sans"))

pango1

pango2 <- gheatmap(pango1, addtl_ano_tree_broad, width=0.05, font.size=0,  offset = 0.3,  color = NA) +
  scale_fill_manual(breaks=addtl_ano_tree_broad_vals$broader_gen, 
                    values=addtl_ano_tree_broad_vals$color_values, name="genotype") + theme(legend.position = 'right', text = element_text(family = "Noto Sans"))

ggsave("pango2_subtree.jpg", pango2, width = 30, height =7, units = "cm", limitsize = F)


#matrix for p/a of accessory genes

all_pr_addt_orf_e_5_mat_pa_acc$contig <- gsub("$", "_", all_pr_addt_orf_e_5_mat_pa_acc$contig)
rownames(all_pr_addt_orf_e_5_mat_pa_acc) <- all_pr_addt_orf_e_5_mat_pa_acc$contig
all_pr_addt_orf_e_5_mat_pa_acc$contig <- NULL
all_pr_addt_orf_e_5_mat_pa_acc <- all_pr_addt_orf_e_5_mat_pa_acc[,c(3,4,1,6,2,7,5)]
all_pr_addt_orf_e_5_mat_pa_acc$E6 <- gsub("1", "E6", all_pr_addt_orf_e_5_mat_pa_acc$E6)
all_pr_addt_orf_e_5_mat_pa_acc$E7 <- gsub("1", "E7", all_pr_addt_orf_e_5_mat_pa_acc$E7)
all_pr_addt_orf_e_5_mat_pa_acc$E1 <- gsub("1", "E1", all_pr_addt_orf_e_5_mat_pa_acc$E1)http://127.0.0.1:10571/graphics/plot_zoom_png?width=727&height=4425
all_pr_addt_orf_e_5_mat_pa_acc$Pap_E4 <- gsub("1", "E4", all_pr_addt_orf_e_5_mat_pa_acc$Pap_E4)
all_pr_addt_orf_e_5_mat_pa_acc$E2 <- gsub("1", "E2", all_pr_addt_orf_e_5_mat_pa_acc$E2)
all_pr_addt_orf_e_5_mat_pa_acc$Papilloma_E5 <- gsub("1", "E5", all_pr_addt_orf_e_5_mat_pa_acc$Papilloma_E5)
all_pr_addt_orf_e_5_mat_pa_acc$Late_protein_L2 <- gsub("1", "L2", all_pr_addt_orf_e_5_mat_pa_acc$Late_protein_L2)


ppp <- pp + new_scale_fill()

pppp <- gheatmap(ppp, all_pr_addt_orf_e_5_mat_pa_acc, width=0.1, font.size=2, 
               colnames_angle=90, color = "grey90", hjust = 2, offset = 3 ) +
  scale_fill_manual(na.value="white", breaks=c("E6", "E7", "E1", "E4", "E2", "E5", "L2", "0"), 
                    values=c("#CBD7CB", "#CBD7CB", "#839183", "#CBD7CB", "#839183", "#CBD7CB","#839183", "white"), name="genes") +
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA) # color = NA removes the border
  )

pppp

ggsave("acc_genes.jpg", width = 50, height =200, units = "cm", limitsize = F)

geo_data <- read.table("wot", sep = ",", header = F)
geo_data_known <- read.table("all_known_geo.tsv", sep = ",", header = F, fill = T)

#using the hierarchy that Teo considered:
#coordinates > geo_loc_name_sam > geo_loc_name_country_calc > geo_loc_name_country_continent_calc
unique(geo_data$V2)
geo_data$V2 <- factor(geo_data$V2, levels = c("lat_lon", "geographic location (latitude),geographic location (longitude)", 
                                              "geo_loc_name_sam", "geo_loc_name", 
                                              "geo_loc_name_country_calc", "geographic location (region and locality)", "region", "birth_location", "INSDC center name"))

geo_data_known$V2 <- factor(geo_data_known$V2, levels = c("lat_lon", "geographic location (latitude),geographic location (longitude)", 
                                              "geo_loc_name_sam", "geo_loc_name", 
                                              "geo_loc_name_country_calc", "geographic location (region and locality)", "region", "birth_location", "INSDC center name"))


geo_data <- geo_data[order(geo_data$V2),]
geo_data_known <- geo_data_known[order(geo_data_known$V2),]


geo_prior <- geo_data[!duplicated(geo_data[,c('V1')]),]
geo_data_known <- geo_data_known[!duplicated(geo_data_known[,c('V1')]),]

geo_data_known <- geo_data_known[order(geo_data_known$V2),]
geo_data_known$status <- c("known")
geo_prior$status <- c("novel")


wkb = structure((geo_data_known$V4), class = "WKB")
geo_data_known$lat_lon <- st_as_sfc(wkb,EWKB = TRUE) %>% st_coordinates()


#need to convert the string to lat long
# install.packages("sf")
library(sf)

# Example WKB string (hexadecimal representation of a POINT)
# wkb = structure(list("0101000020E61000005859DB140F776140DC4944F817D84140"), class = "WKB")
# st_as_sfc(wkb,EWKB = TRUE) %>% st_coordinates()

wkb = structure((geo_prior$V4), class = "WKB")
geo_prior$lat_lon <- st_as_sfc(wkb,EWKB = TRUE) %>% st_coordinates()

geo_prior_all <- as.data.frame(rbind(geo_data_known, geo_prior))
geo_prior_all <- geo_prior_all %>% unnest_wider(lat_lon, names_sep = "_")
colnames(geo_prior_all) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "status", "lat_lon.X", "lat_lon.Y")

biosamp_novel_headers_geo <- left_join(biosamp_novel_headers, geo_prior, by = c("BioSample" = "V1"))
biosamp_novel_headers_geo_all <- left_join(biosamp_novel_headers_geo, addtl_ano, by = c("a" = "lib"))
#manually add information for the generalization
write.table(biosamp_novel_headers_geo_all, "biosamp_novel_headers_geo_all.tsv", sep = "\t", col.names  = T, row.names = F)
biosamp_data_anno <- read.table("2025.10.10biosamp_novel_headers_geo_all.tsv.txt", sep = "\t", header = T, fill = T)

biosamp_data_anno <- unique(biosamp_data_anno)
#map stuff
library(ggmosaic)

world_map <- map_data("world")



geo_prior_all %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), data = world_map, fill = "lightgrey", color = "lightgrey") +
    geom_point(aes(x = lat_lon.X[1], y = lat_lon.X[2], color = status), size = 1.5, alpha = 0.5, stroke = NA) +
    scale_color_manual(values = c("grey40", "#3333E7")) + 
    coord_fixed() +
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
  
  
  dates <- output_df %>%
    ggplot() +
    geom_histogram(aes(x=collection_date_sam_parsed), binwidth = 10) +
    scale_x_date(breaks="2 years") +
    ggtitle(paste0(file_name)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  tree_unroot <- ggtree(final_tree, layout = "circular")
  
  tree_unroot$data$label <- gsub("Edges.*", "", tree_unroot$data$label)
  tree_unroot$data$label <- gsub("ka_f.*", "", tree_unroot$data$label)
  
  
  tree_unroot$data <- left_join(tree_unroot$data, all_genus_curated_thin, by = c("label" = "V1_2"))
  tree_unroot$data$x10 <- iconv(tree_unroot$data$x10, from = "UTF-8", to = "ASCII", sub = "")
  tree_unroot$data[1092, 10] <- "Zhangixalus_dugritei_associated"
  
  tree_unroot$data$status <- ifelse(grepl("[E|S|D][R][R]", tree_unroot$data$label), "novel", "ncbi")
  
  tree_unroot$data$color_values <- ifelse(grepl("ncbi", tree_unroot$data$status), "grey", "#3333E7")
  
  pp_unroot <- tree_unroot %<+% select(tree_unroot$data, label, status) + aes(colour=I(tree_unroot$data$color_values))+
    # geom_tippoint() +
      
    # geom_tiplab(color = "black",size = 2, hjust = 0, linetype = "blank") +
    # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
    # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
    theme(legend.position= "none")
  
  addtl_ano_tree_broad$label <- rownames(addtl_ano_tree_broad)

  novelty_status <- select(tree_unroot$data, label, status)
  addtl_ano_tree_broad <- left_join(addtl_ano_tree_broad, novelty_status, by = c("label" = "label"))
  rownames(addtl_ano_tree_broad) <- addtl_ano_tree_broad$label 
  addtl_ano_tree_broad$label <- NULL
  addtl_ano_tree_broad$stats_color_vals <- ifelse(grepl("novel", addtl_ano_tree_broad$status), "#4646d4", "grey")
  addtl_ano_tree_broad_stats_color_vals <- select(addtl_ano_tree_broad, status, stats_color_vals)
  

  addtl_ano_tree_broad_host <- select(addtl_ano_tree_broad, broader_gen)
  
  p1 <- gheatmap(pp_unroot, addtl_ano_tree_broad_host, offset=1, width=0.2,
           colnames_angle=-45, color = NA, font.size=0) + 
    scale_fill_manual(breaks=addtl_ano_tree_broad$broader_gen, 
                      values=addtl_ano_tree_broad_vals$color_values, 
                      name="status") + 
    theme(legend.position = 'none' )
  
  addtl_ano_tree_broad_stats <- select(addtl_ano_tree_broad, status)
  
  p2 <- p1 + new_scale_fill()
  
  p3 <- gheatmap(p2, addtl_ano_tree_broad_stats, width=0.2, font.size=0, 
           colnames_angle=-45, color = NA) +
    scale_fill_manual(breaks=addtl_ano_tree_broad_stats$status, 
                      values=addtl_ano_tree_broad_stats_color_vals$stats_color_vals, name="genotype") + theme(legend.position = 'none' )
  
  
  ggsave("transp_test_p3.png", width = 20, height =20, units = "cm", limitsize = F)
  
  
  
    
  
