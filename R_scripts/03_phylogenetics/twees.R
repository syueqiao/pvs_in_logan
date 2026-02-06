#twees
library(tidyverse)
library(ggtree)
library(Polychrome)
library(colorspace)

L1_tree <- read.tree("w_bs.treefile")
#bootstrap value is saved as "label" in the tree 
ggtree(L1_tree) + geom_tiplab() + geom_text(aes(label=label), hjust=-.3)

#to get the species labels for the tips, split into pave and SRA samples
L1_tip_labs <- as.data.frame(L1_tree[["tip.label"]])
colnames(L1_tip_labs) <- c("L1_tip")

#deal with pave samples first
#format for a join
pave_L1 <- as.data.frame(grep("[SED][R]{2}", L1_tip_labs$L1_tip, value = T, invert = T ))
colnames(pave_L1) <- c("L1_pave_name")

pave_L1$L1_pave <- pave_L1$L1_pave_name
pave_L1 <- pave_L1 %>% separate(L1_pave, into = c("L1_pave_key", "L1_pave_end"), sep="(?=-)")


L1_info <- read.table("animal_reference_clones.txt", sep = "\t", header = F, na.strings = "")
names(L1_info)[names(L1_info) == 'V1'] <- 'L1_pave_key'
L1_merge <- left_join(pave_L1, L1_info, by = "L1_pave_key")
L1_merge <- L1_merge %>%
  mutate(V6=replace_na(V6, "Homo sapiens"))
#done, should be the base of the datasheet with SRA stuff added on
#deal with SRA stuff
SRA_L1 <- as.data.frame(grep("[SED][R]{2}", L1_tip_labs$L1_tip, value = T ))
colnames(SRA_L1) <- c('L1_SRA')
SRA_L1$L1_SRA_info <- SRA_L1$L1_SRA
SRA_L1 <- SRA_L1 %>% separate(L1_SRA, into = c("library", "info"), sep="_")
#do some stuff with entrez
write.table(SRA_L1$library, "SRA_L1_library.txt", quote = F, row.names = F, col.names = F)
SRA_entrez <- read.table("SRA_info.csv", sep = ",", header = T)

SRA_merge <- left_join(SRA_L1, SRA_entrez, by = "library")

#merge the the two dataframes with just the tip names and the species
SRA_rbind <- select(SRA_merge, L1_SRA_info, V6)
colnames(SRA_rbind) <- c("Newick_label", "species")
pave_rbind <- select(L1_merge, L1_pave_name, V6)
colnames(pave_rbind) <- c("Newick_label", "species")

all_rbind <- rbind(SRA_rbind, pave_rbind)
all_rbind$species_sep <- all_rbind$species
all_rbind <- all_rbind %>% separate(species_sep, into = c("grouped", "end"), sep=" ")
  
#add this information back to the newick tree

p <- ggtree(L1_tree, aes(color = grouped)) %<+% all_rbind

P50 = createPalette(126,  c("#ff0000", "#00ff00", "#0000ff"))

names(P50) <- NULL

#v2 of graph with colored branches and various other values
p + geom_tippoint() + 
  theme(legend.position= c(0.2,0.8), legend.justification = c(0.5, 1)) +
  guides(color = guide_legend(ncol = 1)) +
  geom_text(aes(label = grouped), size = 2, hjust = -0.2, color = "black") +
  geom_tiplab(size = 2, hjust = -0.5, color = "black") +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  geom_nodelab(label = L1_tree$node.label, geom = 'text', size = 1.5) +
  scale_color_manual(values = P50)

ggsave("test.pdf", width = 40, height = 135, units = "cm", limitsize = F)

ncbi_and_novel_tree <- read.tree("FastTree_output_tree (3).nhx")
#bootstrap value is saved as "label" in the tree 
ggtree(ncbi_and_novel_tree) + geom_tiplab() + geom_text(aes(label=label), hjust=-.3)

#to get the species labels for the tips, split into pave and SRA samples
tip_labs_ncbi_novel <- as.data.frame(ncbi_and_novel_tree[["tip.label"]])
colnames(tip_labs_ncbi_novel) <- c("ncbi_and_novel_tip")

#deal with SRA stuff
SRA_ncbi_and_novel <- as.data.frame(grep("[SED][R]{2}", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
colnames(SRA_ncbi_and_novel) <- c('ncbi_and_novel_SRA')
SRA_ncbi_and_novel$ncbi_and_novel_SRA_info <- SRA_ncbi_and_novel$ncbi_and_novel_SRA
SRA_ncbi_and_novel <- SRA_ncbi_and_novel %>% separate(ncbi_and_novel_SRA, into = c("library", "info"), sep="_")
# write.table(SRA_ncbi_and_novel$library, "2025.02.18.ncbi_and_novel_SRA.txt", sep = "/t", col.names = F, row.names = F, quote = T)
#do some stuff with entrez
SRA_entrez_nn <- read.table("ncbi_and_novel_sra_info.txt", sep = "\t", header = F)
colnames(SRA_entrez_nn) <- c("library", "species")
SRA_entrez_nn <- SRA_entrez_nn[!duplicated(SRA_entrez_nn), ]

SRA_merge <- left_join(SRA_ncbi_and_novel, SRA_entrez_nn, by = "library")
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
# write.table(ncbi_nohost$ncbi_acc, "ncbi_nh.txt", sep = "/t", col.names = F, row.names = F, quote = F)

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

ggsave("2025.02.18.ncbi_and_novel_new.pdf", width = 150, height = 350, units = "cm", limitsize = F)

ncbi_and_novel_tree_pyuu <- read.tree("fasttree_ncbi_novel.nhx")


tr = rtree(15)
d=fortify(ncbi_and_novel_tree_pyuu)
dd = subset(d, isTip)
beep <- as.data.frame(d$label[order(dd$y, decreasing=TRUE)])
write.table(beep, "order.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

ggsave("ncbi_and_novel_uhhh.pdf", width = 150, height = 350, units = "cm", limitsize = F)

#to make host comparisons
write.table(ncbi_novel_rbind_info, "ncbi_tags.txt", col.names = F, row.names = F, quote = F, sep = "\t")
#goes up but used this cleaned up version with correct scientific names to create the annotated species tree
ggplot(ncbi_tags_manual, aes(x=species)) + geom_histogram(stat="count")
ncbi_tags_manual_human <- filter(ncbi_tags_manual, species == "Homo sapiens")
ncbi_tags_manual_else <- filter(ncbi_tags_manual, !species == "Homo sapiens")
length(unique(ncbi_tags_manual_else$species))

#make radial
p_round <- ggtree(ncbi_and_novel_tree, aes(color = status), layout='daylight') %<+% known_ncbi_and_novel_all
pp_round <- p_round + geom_tippoint() + 
  theme(legend.position= c(0,0.8), legend.justification = c(0, 1)) +
  guides(color = guide_legend(ncol = 5)) +
  # geom_text(aes(label = species), size = 4, hjust = -0.05, color = "black") +
  # geom_tiplab(size = 2, hjust = -2.5, color = "black") +f
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  scale_color_manual(values = P50)

# known_ncbi_and_novel <- as.data.frame(grep("ncbi_and_novelase.|MGY", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
# known_ncbi_and_novel$species <- c("known")
# colnames(known_ncbi_and_novel) <- c('Newick_label', 'status') 

ppp_round <- pp_round %<+% SRA_merge_info
pppp_round <- ppp_round + geom_tippoint() + geom_text(aes(label = host), size = 4, hjust = -0.05, color = "black", fontface = "bold") +  geom_tiplab(size = 2, hjust = -1.5, color = "white", alpha = 0) 

ppppp_round <- pppp_round %<+% ncbi_tags_manual
# 

ppppp_round + geom_tippoint() + geom_text(aes(label = species), size = 4, hjust = -0.05, color = "black") + geom_tiplab(size = 2, hjust = -2, color = "black") 



ggsave("ncbi_and_novel_daylight.pdf", width = 150, height = 350, units = "cm", limitsize = F)
pp_round + guides(fill="none")

p_round_manual <- ggtree(ncbi_and_novel_tree, aes(color = species), layout='circular') %<+% ncbi_tags_manual_and_sra
pp_round_manual <- p_round_manual + geom_tippoint() + 
  theme(legend.position= c(0,0.8), legend.justification = c(0, 1)) +
  guides(color = guide_legend(ncol = 5)) +
  # geom_text(aes(label = species), size = 4, hjust = -0.05, color = "black") +
  # geom_tiplab(size = 2, hjust = -2.5, color = "black") +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  scale_color_manual(values = P50)
ppp_round_manual <- pp_round_manual %<+% SRA_merge_info


pp_round_manual
ggsave("ncbi_and_novel_round_manual_legend.png", width = 150, height = 200, units = "cm", limitsize = F)

#new tree for 40% aa identity

aa_40_tree <- read.tree("ncbi_and_novel_for_tree_40.aln.treefile")
g <- ggtree(aa_40_tree, aes(color = status))
#bootstrap value is saved as "label" in the tree 
ncbi_tags_manual_and_sra_cust <- ncbi_tags_manual_and_sra
ncbi_tags_manual_and_sra_cust$Newick_label <- gsub("-", "\\.", ncbi_tags_manual_and_sra$Newick_label)
gg <- g %<+% ncbi_tags_manual_and_sra_cust

pp <- gg + geom_tiplab() + geom_text(aes(label=status, size = 20), vjust = 1.5) + geom_tippoint() + 
  theme(legend.position= "none") +
  guides(color = guide_legend(ncol = 1)) +
  # geom_text(aes(label = species), size = 4, hjust = -0.05, color = "black") +
  # geom_tiplab(size = 2, hjust = -2.5, color = "black") +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  scale_color_manual(values = P50)
pp

ggsave("aa_40_tree.pdf", width = 80, height = 25, units = "cm", limitsize = F)

aa_40_tree_pyuu <- read.tree("ncbi_and_novel_for_tree_40.aln.treefile")

d=fortify(aa_40_tree_pyuu)
dd = subset(d, isTip)
beep <- as.data.frame(d$label[order(dd$y, decreasing=TRUE)])
write.table(beep, "order.tsv", sep = "\t", quote = F, col.names = F, row.names = F)



ppp <- pp %<+% SRA_merge_info
pppp <- ppp + geom_tippoint() + geom_text(aes(label = host), size = 4, hjust = -0.05, color = "black", fontface = "bold") +  geom_tiplab(size = 2, hjust = -1.5, color = "white", alpha = 0) 

ppppp <- pppp %<+% ncbi_tags_manual
pyup <- ppppp + geom_tippoint() + geom_text(aes(label = uhhh), size = 4, hjust = -0.05, color = "black") + geom_tiplab(size = 2, hjust = -2, color = "black") 



P50 = createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))

names(P50) <- NULL

#v2 of graph with colored branches and various other values
pp <- p + geom_tippoint() + 
  theme(legend.position= c(0,0.8), legend.justification = c(0, 1)) +
  guides(color = guide_legend(ncol = 1)) +
  # geom_text(aes(label = species), size = 4, hjust = -0.05, color = "black") +
  # geom_tiplab(size = 2, hjust = -2.5, color = "black") +
  # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
  geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
  scale_color_manual(values = P50)

# known_ncbi_and_novel <- as.data.frame(grep("ncbi_and_novelase.|MGY", tip_labs_ncbi_novel$ncbi_and_novel_tip, value = T ))
# known_ncbi_and_novel$species <- c("known")
# colnames(known_ncbi_and_novel) <- c('Newick_label', 'status') 

ppp <- pp %<+% SRA_merge_info
pppp <- ppp + geom_tippoint() + geom_text(aes(label = host), size = 4, hjust = -0.05, color = "black", fontface = "bold") +  geom_tiplab(size = 2, hjust = -1.5, color = "white", alpha = 0) 

ppppp <- pppp %<+% ncbi_tags_manual
# 
ncbi_and_novel_tree_pyuu <- read.tree("fasttree_ncbi_novel.nhx")

pyup <- ppppp + geom_tippoint() + geom_text(aes(label = uhhh), size = 4, hjust = -0.05, color = "black") + geom_tiplab(size = 2, hjust = -2, color = "black") 


tree <- rcoal(5)
tree$tip.label <- c("REF-1","REF-2","S1","S2","S3")
full_tree    <- ggtree(tree) + xlim(0,2) + geom_tiplab()
full_tree
#Since we know the first two are reference, use the order
limited_tree1 <- ggtree(tree) + xlim(0,2) + geom_tiplab(aes(subset=(node %in% c(1,2))))
limited_tree1
#But we also know the names are different, so we can use that as well
#grepl() searches for 'RE' in the label, and returns true if found
#we use aes(subset=()) to only print those with the subset condition == TRUE
limited_tree2 <- ggtree(tree) + xlim(0,2) + geom_tiplab(aes(
  subset=(grepl('RE',label,fixed=TRUE))==TRUE))
limited_tree2
