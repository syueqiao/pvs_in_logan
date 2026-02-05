#process unknown ORFs for mikko synthesis
lengths <- read.table("../mikko/unannotated_lengths_clean.txt", header = F, fill = T, sep = "\t")
# hist(lengths$V2)
ggplot(lengths, aes(x=V2)) +
  geom_histogram(binwidth = 10, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

lengths_partial <- filter(lengths, grepl("partial",V1))
lengths_full <- filter(lengths, !V1 %in% lengths_partial$V1)

write.table(lengths_full$V1, "../mikko/full_protein_ids.txt", sep = "\t", col.names = F, row.names = F, quote = F)

#remove FPs in bash, oops
FP_list$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", FP_list$V1)
write.table(FP_list$contig, "../mikko/FP_contigs.txt", quote = F, col.names = F, row.names = F)



nonsingleton_clusters <- read.table("../mikko/unannotated_clean_full_clusters_no_singles.uc", header = F, fill = T, sep = "\t")
nonsingleton_clusters_size <- left_join(nonsingleton_clusters, lengths_full, by = c("V9" = "V1"))
nonsingleton_clusters_size_100plus <- filter(nonsingleton_clusters_size, V2.y >= 100)
nonsingleton_clusters_size_75plus <- filter(nonsingleton_clusters_size, V2.y >= 75)
nonsingleton_clusters_size_150plus <- filter(nonsingleton_clusters_size, V2.y >= 150)


#extract cluster numbers of interest
#get IDs from full .uc file, and extract from the tab file to get all members of the cluster


ggplot(nonsingleton_clusters_size_150plus, aes(x=V3, y=V2.y)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "D") +
  xlab("cluster membership") +
  ylab("aa length")

nonsingleton_clusters_size_150plus$contig <-  str_match(nonsingleton_clusters_size_150plus$V9, "_\\s*(.*?)\\s*:")
nonsingleton_clusters_size_150plus <- filter(nonsingleton_clusters_size_150plus, contig[,2 ] %in% pv2_fasta_sto_sto$contig)




nonsingleton_clusters_98_50 <- read.table("../mikko/unannotated_clean_full_clusters_98_50_no_singles.uc", header = F, fill = T, sep = "\t")
nonsingleton_clusters_98_50_lengths <- left_join(nonsingleton_clusters_98_50, lengths_full, by = c("V9" = "V1"))
res <- as.table(str_match(nonsingleton_clusters_98_50_lengths$V9, "_\\s*(.*?)\\s*:"))
nonsingleton_clusters_98_50_lengths$contig <-  res[,2] 

nonsingleton_clusters_98_50_lengths <- filter(nonsingleton_clusters_98_50_lengths, contig %in% pv2_fasta_sto_sto$contig)

nonsingleton_clusters_98_50_lengths_100plus <- filter(nonsingleton_clusters_98_50_lengths, V2.y >= 100)

ggplot(nonsingleton_clusters_98_50_lengths_100plus, aes(x=log(V3), y=log(V2.y))) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "D") +
  xlab("cluster membership") +
  ylab("aa length")

write.table(nonsingleton_clusters_98_50_lengths_100plus$V9, "../mikko/nonsingleton_clusters_98_50_lengths_100plus.txt", quote = F, col.names = F, row.names = F)

nonsingleton_clusters_98_50_lengths_100plus_meta <- left_join(nonsingleton_clusters_98_50_lengths_100plus, concat_pv2_fil_pv_FP, by = c("contig" = "contig_clean"))
nonsingleton_clusters_98_50_lengths_100plus_meta_forward <- nonsingleton_clusters_98_50_lengths_100plus_meta

cluster_858_ids <- read.table("../mikko/cluster_858_ids.txt")
res <- as.table(str_match(cluster_858_ids$V1, "_\\s*(.*?)\\s*:"))
cluster_858_ids$contig <-  res[,2] 
concat_pv2_fil_pv_FP$contig_clean <- str_extract(concat_pv2_fil_pv_FP$contig, "[^_]*_[^_]*")
cluster_858_ids_md <- left_join(cluster_858_ids, concat_pv2_fil_pv_FP, by = c("contig" = "contig_clean"))

cluster_514_ids <- read.table("../mikko/cluster_514_ids.txt")
res <- as.table(str_match(cluster_514_ids$V1, "_\\s*(.*?)\\s*:"))
cluster_514_ids$contig <-  res[,2] 
concat_pv2_fil_pv_FP$contig_clean <- str_extract(concat_pv2_fil_pv_FP$contig, "[^_]*_[^_]*")
cluster_514_ids_md <- left_join(cluster_514_ids, concat_pv2_fil_pv_FP, by = c("contig" = "contig_clean"))

cluster_876_ids <- read.table("../mikko/cluster_876_ids.txt")
res <- as.table(str_match(cluster_876_ids$V1, "_\\s*(.*?)\\s*:"))
cluster_876_ids$contig <-  res[,2] 
concat_pv2_fil_pv_FP$contig_clean <- str_extract(concat_pv2_fil_pv_FP$contig, "[^_]*_[^_]*")
cluster_876_ids_md <- left_join(cluster_876_ids, concat_pv2_fil_pv_FP, by = c("contig" = "contig_clean"))

cluster_153_ids <- read.table("../mikko/cluster_153_ids.txt")
res <- as.table(str_match(cluster_153_ids$V1, "_\\s*(.*?)\\s*:"))
cluster_153_ids$contig <-  res[,2] 
concat_pv2_fil_pv_FP$contig_clean <- str_extract(concat_pv2_fil_pv_FP$contig, "[^_]*_[^_]*")
cluster_153_ids_md <- left_join(cluster_153_ids, concat_pv2_fil_pv_FP, by = c("contig" = "contig_clean"))


nonsingleton_clusters_98_50_lengths_100plus_meta_forward$og_id <- nonsingleton_clusters_98_50_lengths_100plus_meta_forward$V9
nonsingleton_clusters_98_50_lengths_100plus_meta_forward <- separate(nonsingleton_clusters_98_50_lengths_100plus_meta_forward, col=V9, into=c('info', 'start_orf', 'end_orf'), sep="\\:")
nonsingleton_clusters_98_50_lengths_100plus_meta_forward$end_orf <- sub(" .*", "", nonsingleton_clusters_98_50_lengths_100plus_meta_forward$end_orf)
nonsingleton_clusters_98_50_lengths_100plus_meta_forward$end_orf <- as.numeric(nonsingleton_clusters_98_50_lengths_100plus_meta_forward$end_orf)
nonsingleton_clusters_98_50_lengths_100plus_meta_forward$start_orf <- as.numeric(nonsingleton_clusters_98_50_lengths_100plus_meta_forward$start_orf)

nonsingleton_clusters_98_50_lengths_100plus_meta_forward$orf_dir <- nonsingleton_clusters_98_50_lengths_100plus_meta_forward$start_orf - nonsingleton_clusters_98_50_lengths_100plus_meta_forward$end_orf
nonsingleton_clusters_98_50_lengths_100plus_meta_forward$dir <- "fw"
nonsingleton_clusters_98_50_lengths_100plus_meta_forward$dir[nonsingleton_clusters_98_50_lengths_100plus_meta_forward$orf_dir <0] <- "rev"

nonsingleton_clusters_98_50_lengths_100plus_meta_forward$model_dir <- nonsingleton_clusters_98_50_lengths_100plus_meta_forward$start - nonsingleton_clusters_98_50_lengths_100plus_meta_forward$end
nonsingleton_clusters_98_50_lengths_100plus_meta_forward$dir_model <- "fw"
nonsingleton_clusters_98_50_lengths_100plus_meta_forward$dir_model[nonsingleton_clusters_98_50_lengths_100plus_meta_forward$model_dir <0] <- "rev"

nonsingleton_clusters_98_50_lengths_100plus_meta_forward_same_direct <- filter(nonsingleton_clusters_98_50_lengths_100plus_meta_forward, dir_model == "fw" & dir == "fw")
nonsingleton_clusters_98_50_lengths_100plus_meta_rev_same_direct <- filter(nonsingleton_clusters_98_50_lengths_100plus_meta_forward, dir_model == "rev" & dir == "rev")
nonsingleton_clusters_98_50_lengths_100plus_meta_same_direct <- rbind(nonsingleton_clusters_98_50_lengths_100plus_meta_rev_same_direct, nonsingleton_clusters_98_50_lengths_100plus_meta_forward_same_direct)

nonsingleton_clusters_98_50_lengths_100plus_meta_same_direct$og_id <- sub(" .*", "", nonsingleton_clusters_98_50_lengths_100plus_meta_same_direct$og_id)

nonsingleton_clusters_98_50_lengths_100plus_meta_same_direct_top <- nonsingleton_clusters_98_50_lengths_100plus_meta_same_direct %>% 
  arrange(desc(V3))
write.table(nonsingleton_clusters_98_50_lengths_100plus_meta_same_direct_top$og_id, "../mikko/nonsingleton_clusters_98_50_lengths_100plus_meta_same_direct_top.txt", quote = F, col.names = F, row.names = F)
