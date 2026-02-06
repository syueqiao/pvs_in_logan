%%R
%%R
#read in NCBI host file
ncbi_nn <- read.table("files/ncbi_info.txt", sep = "\t", header = F)
colnames(ncbi_nn) <- c('ncbi_acc', 'host') 
ncbi_nn$host <- as.character(ncbi_nn$host)
ncbi_nn$ncbi_acc <- as.character(ncbi_nn$ncbi_acc)
ncbi_nn$ncbi_acc <- stringr::str_trim(ncbi_nn$ncbi_acc)

ncbi_merge <- left_join(ncbi_tips, ncbi_nn, by = "ncbi_acc")

ncbi_host <- subset(ncbi_merge, trimws(host) !="")
ncbi_nohost <- subset(ncbi_merge, trimws(host) == "")
#get ones that didnt have host annotation
write.table(ncbi_nohost$ncbi_acc, "outputs/ncbi_nh.txt", sep = "/t", col.names = F, row.names = F, quote = F)
#using sra metadata file, grep two columns that are the sra library, and associated host id
#what this means is, get the list of novel libraries from above, input into https://www.ncbi.nlm.nih.gov/sites/batchentrez and download metadata from webportal
#save as tab delim file
SRA_entrez_nn <- read.table("files/ncbi_and_novel_sra_info.txt", sep = "\t", header = F)
#clean up a little
colnames(SRA_entrez_nn) <- c("library", "species")
SRA_entrez_nn <- SRA_entrez_nn[!duplicated(SRA_entrez_nn), ]

#annotate information on to the sra-only dataframe
SRA_merge <- left_join(SRA_ncbi_and_novel, SRA_entrez_nn, by = "library")
SRA_merge <- select(SRA_merge, ncbi_and_novel_SRA_info, species)

#repeat process for ncbi sequences
ncbi_tips <- filter(tip_labs_ncbi_novel, !ncbi_and_novel_tip %in% SRA_ncbi_and_novel$ncbi_and_novel_SRA_info)
ncbi_tips$ncbi_acc <- sub("-.*", "", ncbi_tips$ncbi_and_novel_tip)
ncbi_tips$ncbi_acc <- as.character(ncbi_tips$ncbi_acc)
write.table(ncbi_tips$ncbi_acc, "outputs/TESTncbi_tips.txt", sep = "/t", col.names = F, row.names = F, quote = F)

ncbi_nohost <- subset(ncbi_merge, trimws(host) == "")

#fill with next best thing, species
ncbi_nh <- read.table("files/ncbi_info_fill.txt", sep = "\t", header = F)
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

write.table(ncbi_novel_rbind_info_merge, "outputs/TESTncbi_tags_manual.txt", sep = "/t", col.names = F, row.names = F, quote = F)
#for those that are species, go back and manually annotate by searching up the PV species (if the name is non informative, eg. alphapapillomavirus 49)

#read back in for complete ncbi stuff
colnames(ncbi_tags_manual) <- c('Newick_label', 'status') 
ncbi_tags_manual_and_sra <- rbind(SRA_merge_info_merge3, ncbi_tags_manual)
colnames(ncbi_tags_manual) <- c('Newick_label', 'uhhh') 

#finally, bind with the SRA stuff
known_ncbi_and_novel_all <- rbind(ncbi_novel_rbind_info_merge, SRA_merge_info_merge2)
colnames(known_ncbi_and_novel_all) <- c('Newick_label', 'status')

#add this information back to the newick tree
p <- ggtree(ncbi_and_novel_tree, aes(color = status)) %<+% ncbi_tags_manual_and_sra

P50 = createPalette(251,  c("#ff0000", "#00ff00", "#0000ff"))

names(P50) <- NULL

#v2 of graph with colored branches and various other values
pp <- p + geom_tippoint() + 
  # geom_text(aes(label = species), size = 4, hjust = -0.05, color = "black") +
  # geom_tiplab(size = 2, hjust = -2.5, color = "black") +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  # scale_color_manual(values = P50) +
  theme(legend.position= "none")+ 
  scale_color_manual(values = P50)

# known_ncbi_and_novel <- as.data.frame(grep("ncbi_and_novelase.|MGY", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
# known_ncbi_and_novel$species <- c("known")
# colnames(known_ncbi_and_novel) <- c('Newick_label', 'status') 

ppp <- pp %<+% SRA_merge_info
pppp <- ppp + geom_tippoint() + geom_text(aes(label = host), size = 4, hjust = -0.05, color = "black", fontface = "bold") +  geom_tiplab(size = 2, hjust = -1.5, color = "white", alpha = 0) 

ppppp <- pppp %<+% ncbi_tags_manual
pppppp <- ppppp + geom_tippoint() + geom_text(aes(label = uhhh), size = 4, hjust = -0.05, color = "black") + geom_tiplab(size = 2, hjust = -2, color = "black")
# 

ggsave("outputs/2025.02.18.ncbi_and_novel_new.pdf", width = 150, height = 350, units = "cm", limitsize = F)