#e6 e7 l1 trees
library(ggtree)
library("devtools")
install.packages("ggnewscale")
install_github('acarafat/tangler')
library(tangler)

l1_tree <- read.tree("l1_for_tree.aln.treefile")
e6_tree <- read.tree("e6_for_tree.aln.treefile")

#bootstrap value is saved as "label" in the tree 
l1 <- ggtree(l1_tree) + geom_tiplab()
e6 <- ggtree(e6_tree) + geom_tiplab()


#to get the species labels for the tips, split into pave and SRA samples
l1_tree_cols <- as.data.frame(l1_tree[["tip.label"]])
colnames(l1_tree_cols) <- c("label")

##

ncolors <- 1078
# coolors = createPalette(ncolors,  c("#ff0000", "#00ff00", "#0000ff"))
names(coolors) <- NULL
l1_tree_cols <- as.data.frame(cbind(l1_tree_cols$label, coolors))

# Load meta
# meta=read.csv('tree_meta.csv', header=T) 

# Annotate trees
ggtree(l1_tree) %<+% l1_tree_cols + geom_tippoint(pch=16, aes(col=coolors))+ theme(legend.position = "none") +
  scale_color_manual(values = l1_tree_cols$coolors)

#cluster_lowe6_sep from earlier e6/e7 processing

l1_tree_cols$coolors_2 <- with(l1_tree_cols, ifelse(V1 %in% askansgag_2$aaa, "#ff0000", l1_tree_cols$coolors)) 

l1_tree_cols_new <- as.data.frame(cbind(l1_tree_cols$V1, l1_tree_cols$coolors_2))

p <- ggtree(l1_tree) %<+% l1_tree_cols_new
 
pp <- p + geom_tippoint(pch=16, aes(col=p$data$label)) +
  scale_color_manual(values = p$data$V2) + theme(legend.position = "none") + geom_tiplab(size = 2, hjust = -2.5, color = "black")

l1_tree_cols_new_new <- filter(l1_tree_cols_new, V1 %in% askansgag_2$aaa)
l1_tree_cols_new_new$V2 <- c("IN GROUP")

ppp <- pp %<+% l1_tree_cols_new_new
pppp <- ppp + geom_text(aes(label = ppp$data$V2.y), size = 4, hjust = -0.05, color = "black", fontface = "bold") 

pppp



##################e6_tree#################
e6_tree_cols <- as.data.frame(e6_tree[["tip.label"]])
colnames(e6_tree_cols) <- c("label")

##

ncolors <- length(e6_tree_cols$label)
coolors = createPalette(ncolors,  c("#ff0000", "#00ff00", "#0000ff"))
names(coolors) <- NULL
e6_tree_cols <- as.data.frame(cbind(e6_tree_cols$label, coolors))

# Load meta
# meta=read.csv('tree_meta.csv', header=T) 

# Annotate trees
ggtree(e6_tree) %<+% e6_tree_cols + geom_tippoint(pch=16, aes(col=coolors))+ theme(legend.position = "none") +
  scale_color_manual(values = e6_tree_cols$coolors)

#cluster_lowe6_sep from earlier e6/e7 processing

common_contigs <- intersect(str_trim(cluster_lowe6_sep$id2), str_trim(cluster_lowe6_sep$id3))

e6_tree_cols$coolors_2 <- with(e6_tree_cols, ifelse(V1 %in% askansgag_2$aaa, "#ff0000", e6_tree_cols$coolors)) 

e6_tree_cols_new <- as.data.frame(cbind(e6_tree_cols$V1, e6_tree_cols$coolors_2))

p <- ggtree(e6_tree) %<+% e6_tree_cols_new

pp <- p + geom_tippoint(pch=5, aes(col=p$data$label)) +
  scale_color_manual(values = p$data$V2) + theme(legend.position = "none") + geom_tiplab(size = 2, hjust = -2.5, color = "black")

pp

e6_tree_cols_new_new <- filter(e6_tree_cols_new, V1 %in% askansgag_2$aaa)
e6_tree_cols_new_new$V2 <- c("IN GROUP")

ppp <- pp %<+% e6_tree_cols_new_new
pppp <- ppp + geom_text(aes(label = ppp$data$V2.y), size = 4, hjust = -0.05, color = "black", fontface = "bold") 
pppp

ggplot2::ggsave(filename = "l1_tree_hilite_diverse_TEXT.pdf" , pppp, width = 50, height = 150, units = "cm", limitsize = F)



viz <- as.data.frame(table(c(e6e7_all_all_sel_nona$id2, e6e7_all_all_sel_nona$id3)))

ggplot(viz, aes(x=Freq)) +
  geom_histogram(binwidth=25, alpha=0.5, color = "white", fill = "magenta") +
  theme_bw() +
  ggtitle("entire l1/e6 comparison")

hist(viz$Freq)

viz2 <- as.data.frame(table(c(cluster_lowe6_sep$id2, cluster_lowe6_sep$id3)))
hist(viz2$Freq)

ggplot(viz2, aes(x=Freq)) +
  geom_histogram(binwidth=25, alpha=0.5, color = "white", fill = "blue") +
  theme_bw() +
  ggtitle("cluster only")

####################################################################################################################################################################################
l1_tree_common <- read.tree("common_l1s_99_cent.aln.treefile")

#bootstrap value is saved as "label" in the tree 
l1_c <- ggtree(l1_tree_common) + geom_tiplab()


#to get the species labels for the tips, split into pave and SRA samples
l1_tree_cols_c <- as.data.frame(l1_tree_common[["tip.label"]])
colnames(l1_tree_cols_c) <- c("label")

##

ncolors <- length(unique(l1_tree_cols_c$label))
coolors_c = createPalette(ncolors,  c("#ff0000", "#00ff00", "#0000ff"))
names(coolors_c) <- NULL
l1_tree_cols_c <- as.data.frame(cbind(l1_tree_cols_c$label, coolors_c))

# Load meta
# meta=read.csv('tree_meta.csv', header=T) 

# Annotate trees
ggtree(l1_tree_common) %<+% l1_tree_cols_c + geom_tippoint(pch=16, aes(col=coolors_c))+ theme(legend.position = "none") +
  scale_color_manual(values = l1_tree_cols_c$coolors_c)

#cluster_lowe6_sep from earlier e6/e7 processing

l1_tree_cols_c$coolors_2_c <- with(l1_tree_cols_c, ifelse(V1 %in% askansgag_2$aaa, "#0000ff", l1_tree_cols_c$coolors_c)) 

l1_tree_cols_new_c <- as.data.frame(cbind(l1_tree_cols_c$V1, l1_tree_cols_c$coolors_2_c))

p <- ggtree(l1_tree_common) %<+% l1_tree_cols_new_c

p

pipi <- p + geom_tippoint(pch=16, aes(col=p$data$label)) + geom_tiplab(label = p$data$label, hjust = -1) +
  theme(legend.position = "none")

pipi

l1_tree_cols_new_new_c <- filter(l1_tree_cols_new_c, V1 %in% askansgag_2$aaa)
l1_tree_cols_new_new_c$V2 <- c("IN GROUP")

pepo <- pipi %<+% l1_tree_cols_new_new_c
pepo

popo <- pepo + geom_text(aes(label = pepo$data$V2.y), size = 4, hjust = -0.05, color = "black", fontface = "bold") 

popo

########################################
e6_tree_common <- read.tree("common_e6s_99_cent.aln.treefile")

e6_c <- ggtree(e6_tree) + geom_tiplab()


e6_tree_cols_c <- as.data.frame(e6_tree_common[["tip.label"]])
colnames(e6_tree_cols_c) <- c("label")

##

ncolors <- length(unique(e6_tree_cols_c$label))
coolors_c = createPalette(ncolors,  c("#ff0000", "#00ff00", "#0000ff"))
names(coolors_c) <- NULL
e6_tree_cols_c <- as.data.frame(cbind(e6_tree_cols_c$label, coolors_c))

# Load meta
# meta=read.csv('tree_meta.csv', header=T) 

# Annotate trees
ggtree(e6_tree_common) %<+% e6_tree_cols_c + geom_tippoint(pch=16, aes(col=coolors_c))+ theme(legend.position = "none") +
  scale_color_manual(values = e6_tree_cols_c$coolors_c)

#cluster_lowe6_sep from earlier e6/e7 processing

e6_tree_cols_c$coolors_2_c <- with(e6_tree_cols_c, ifelse(V1 %in% askansgag_2$aaa, "#0000ff", e6_tree_cols_c$coolors_c)) 

e6_tree_cols_new_c <- as.data.frame(cbind(e6_tree_cols_c$V1, e6_tree_cols_c$coolors_2_c))

p <- ggtree(e6_tree_common) %<+% e6_tree_cols_new_c

p

pipi <- p + geom_tippoint(pch=16, aes(col=p$data$label)) + geom_tiplab(label = p$data$label, hjust = -1) +
  theme(legend.position = "none")

pipi

e6_tree_cols_new_new_c <- filter(e6_tree_cols_new_c, V1 %in% askansgag_2$aaa)
e6_tree_cols_new_new_c$V2 <- c("IN GROUP")

pepo <- pipi %<+% e6_tree_cols_new_new_c
pepo

popo <- pepo + geom_text(aes(label = pepo$data$V2.y), size = 4, hjust = -0.05, color = "black", fontface = "bold") 

popo

write.table()

#investigate what the "diverged" contigs are, looking at blast results

l1_blast <- read.table("l1_blast.tsv", sep= "\t", header = F, fill = T)
#slice for the most "likely" hit using e-value?

l1_blast_eval_slice <- l1_blast %>% group_by(V1) %>% slice_min(n = 1, V11)
l1_blast_eval_slice <- l1_blast_eval_slice[!duplicated(l1_blast_eval_slice$V1),]
l1_blast_eval_slice_sel <- select(l1_blast_eval_slice, V1, V10, V11, V13)
colnames(l1_blast_eval_slice_sel) <- c("query_acc", "e_val", "l1_hit")

e6_blast <- read.table("e6_blast.tsv", sep= "\t", header = F, fill = T)
#slice for the most "likely" hit using e-value?

e6_blast_eval_slice <- e6_blast %>% group_by(V1) %>% slice_min(n = 1, V11)
e6_blast_eval_slice <- e6_blast_eval_slice[!duplicated(e6_blast_eval_slice$V1),] 
e6_blast_eval_slice_sel <- select(e6_blast_eval_slice, V1, V10, V11, V13)
colnames(e6_blast_eval_slice_sel) <- c("query_acc", "e_val", "e6_hit")

#these are the ones that are diverged: %in% askansgag_2$aaa
aaa <- c(cluster_lowe6_sep$id3, cluster_lowe6_sep$id2)
askansgag <- as.data.frame(table(as.data.frame(aaa)))
askansgag_2 <- filter(askansgag, Freq >= 600)

askansgag$aaa

l1e6_blast <- left_join(l1_blast_eval_slice_sel, e6_blast_eval_slice_sel, by = join_by(query_acc))
write.table(askansgag_2$aaa, "diverged.list", sep = "/t", col.names = F, row.names = F, quote = F)
