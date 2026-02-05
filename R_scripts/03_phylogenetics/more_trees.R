ncbi_and_novel_tree_final <- read.tree("FastTree_output_tree (3).nhx")
#bootstrap value is saved as "label" in the tree 
ggtree(ncbi_and_novel_tree_final) + geom_tiplab() + geom_text(aes(label=label), hjust=-.3)

#to get the species labels for the tips, split into pave and SRA samples
tip_labs_ncbi_novel <- as.data.frame(ncbi_and_novel_tree_final[["tip.label"]])
colnames(tip_labs_ncbi_novel) <- c("ncbi_and_novel_tip")

#deal with SRA stuff
SRA_ncbi_and_novel <- as.data.frame(grep("[SED][R]{2}", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
colnames(SRA_ncbi_and_novel) <- c('ncbi_and_novel_SRA')
SRA_ncbi_and_novel$ncbi_and_novel_SRA_info <- SRA_ncbi_and_novel$ncbi_and_novel_SRA
SRA_ncbi_and_novel <- SRA_ncbi_and_novel %>% separate(ncbi_and_novel_SRA, into = c("library", "info"), sep="_")
SRA_ncbi_and_novel$library <- gsub("Score.*", "", SRA_ncbi_and_novel$library)

write.table(SRA_ncbi_and_novel$library, "2025.03.27.ncbi_and_novel_SRA.txt", sep = "/t", col.names = F, row.names = F, quote = T)
#do some stuff with entrez
SRA_entrez_nn <- read.table("sra_metadata_all_sra.txt", sep = "\t", header = F)
colnames(SRA_entrez_nn) <- c("library", "species")
SRA_entrez_nn <- SRA_entrez_nn[!duplicated(SRA_entrez_nn), ]

SRA_merge <- left_join(SRA_ncbi_and_novel, SRA_entrez_nn, by = "library")
SRA_merge$pr_lib <- gsub("_.*", "", SRA_merge$ncbi_and_novel_SRA_info)
SRA_merge <- select(SRA_merge, ncbi_and_novel_SRA_info, species)
#for ncbi sequences
ncbi_tips <- filter(tip_labs_ncbi_novel, !ncbi_and_novel_tip %in% SRA_ncbi_and_novel$ncbi_and_novel_SRA_info)
ncbi_tips$ncbi_acc <- sub("-.*", "", ncbi_tips$ncbi_and_novel_tip)
ncbi_tips$ncbi_acc <- as.character(ncbi_tips$ncbi_acc)

ncbi_nn <- read.table("ncbi_info.txt", sep = "\t", header = F)
colnames(ncbi_nn) <- c('ncbi_acc', 'host') 
ncbi_nn$host <- as.character(ncbi_nn$host)
ncbi_nn$ncbi_acc <- as.character(ncbi_nn$ncbi_acc)
ncbi_nn$ncbi_acc <- stringr::str_trim(ncbi_nn$ncbi_acc)

ncbi_merge <- left_join(ncbi_tips, ncbi_nn, by = "ncbi_acc")

ncbi_host <- subset(ncbi_merge, trimws(host) !="")
ncbi_nohost <- subset(ncbi_merge, trimws(host) == "")
write.table(ncbi_nohost$ncbi_acc, "ncbi_nh.txt", sep = "/t", col.names = F, row.names = F, quote = F)

#fill with next best thing, species
ncbi_nh <- read.table("ncbi_info_fill_redo.txt", sep = "\t", header = F)
colnames(ncbi_nh) <- c('ncbi_acc', 'host') 
ncbi_nh$ncbi_acc <- stringr::str_trim(ncbi_nh$ncbi_acc)

ncbi_nohost$host <- NULL
ncbi_merge_nohost <- left_join(ncbi_nohost, ncbi_nh, by = "ncbi_acc")

ncbi_novel_rbind <- rbind(ncbi_merge_nohost, ncbi_host)

#combine the two labels
SRA_merge_info <- select(SRA_merge, ncbi_and_novel_SRA_info, species)
colnames(SRA_merge_info) <- c('Newick_label', 'host') 
SRA_merge_info_merge2 <- SRA_merge_info
colnames(SRA_merge_info_merge2) <- c('Newick_label', 'status')
SRA_merge_info_merge3 <- SRA_merge_info_merge2
colnames(SRA_merge_info_merge3) <- c('Newick_label', 'status') 


ncbi_novel_rbind_info <- select(ncbi_novel_rbind, ncbi_and_novel_tip, host)
colnames(ncbi_novel_rbind_info) <- c('Newick_label', 'species') 
ncbi_novel_rbind_info_merge <- ncbi_novel_rbind_info
colnames(ncbi_novel_rbind_info_merge) <- c('Newick_label', 'status') 


ncbi_tags_manual <- read.table("ncbi_tags_manual.txt", sep = "\t", header = T)
colnames(ncbi_tags_manual) <- c('Newick_label', 'status') 
ncbi_tags_manual_and_sra <- rbind(SRA_merge_info_merge3, ncbi_tags_manual)
colnames(ncbi_tags_manual) <- c('Newick_label', 'uhhh') 


known_ncbi_and_novel_all <- rbind(ncbi_novel_rbind_info_merge, SRA_merge_info_merge2)
colnames(known_ncbi_and_novel_all) <- c('Newick_label', 'status')

length(unique(ncbi_tags_manual_and_sra$status))
#add this information back to the newick tree
p <- ggtree(ncbi_and_novel_tree_final, aes(color = status)) %<+% ncbi_tags_manual_and_sra

P50 = createPalette(1200,  c("#ff0000", "#00ff00", "#0000ff"))

names(P50) <- NULL

#v2 of graph with colored branches and various other values
pp <-  p + geom_tippoint() + 
  theme(legend.position= c(0,0.8), legend.justification = c(0, 1)) +
  guides(color = guide_legend(ncol = 1)) +
  # geom_text(aes(label = species), size = 4, hjust = -0.05, color = "black") +
  # geom_tiplab(size = 2, hjust = -2.5, color = "black") +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  scale_color_manual(values = P50)


pp

# known_ncbi_and_novel <- as.data.frame(grep("ncbi_and_novelase.|MGY", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
# known_ncbi_and_novel$species <- c("known")
# colnames(known_ncbi_and_novel) <- c('Newick_label', 'status') 

ppp <- pp %<+% SRA_merge_info
pppp <- ppp + geom_tippoint() + geom_text(aes(label = host), size = 4, hjust = -0.05, color = "black", fontface = "bold") +  geom_tiplab(size = 2, hjust = -1.5, color = "white", alpha = 0) 

ppppp <- pppp %<+% ncbi_tags_manual
pppppp <- ppppp + geom_tippoint() + geom_text(aes(label = uhhh), size = 4, hjust = -0.05, color = "black") + geom_tiplab(size = 2, hjust = -2, color = "black")
# 

ggsave("2025.03.28.ncbi_and_novel_new.pdf", width = 250, height = 350, units = "cm", limitsize = F)

ncbi_and_novel_tree_final_pyuu <- read.tree("fasttree_ncbi_novel.nhx")


tr = rtree(15)
d=fortify(ncbi_and_novel_tree_final_pyuu)
dd = subset(d, isTip)
beep <- as.data.frame(d$label[order(dd$y, decreasing=TRUE)])
write.table(beep, "order.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

ggsave("ncbi_and_novel_uhhh.pdf", width = 150, height = 350, units = "cm", limitsize = F)