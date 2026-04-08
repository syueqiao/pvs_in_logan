#logan3
library(tidyverse)
library(ggplot2)
#function to parse the domtbl output
source("R_scripts/00_utilities/hmmsearch_utils.R")

#all clustered at 90% nuceltodie, transalted to amino acid, then searched w./ the jrHMM
#next, tabulate which are "full"
all_1100_clustered <- hmmsearch_clean("files/L1_BI_hmmer_iter2.domtbl")

all_1100_clustered$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", all_1100_clustered$query_acc)
all_1100_clustered$library <- sub("_.*", "\\1", all_1100_clustered$contig)

#record region presence in contigs
#"low fruit" pass
all_1100_clustered_mat <- all_1100_clustered %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

#filter for those that have all regions hit :)
all_1100_clustered_mat_all <- all_1100_clustered_mat %>%
  filter(if_all(ncbi_and_novel_L1_B:ncbi_and_novel_L1_I, ~ .x > 0))

all_1100_clustered_mat_all$library <- sub("_.*", "\\1", all_1100_clustered_mat_all$query_acc)

missing <- all_1100_clustered_mat[!all_1100_clustered_mat$query_acc %in% all_1100_clustered_mat_all$query_acc,]

write.table(select(all_1100_clustered_mat_all, library), "outputs/feb_7_pv_L1_BI_mat_full_library.list", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(select(all_1100_clustered_mat_all, query_acc), "outputs/feb_7_pv_L1_BI_mat_full.list", quote = FALSE, col.names = FALSE, row.names = FALSE)

#do a REALLY messy extraction of the coordinates
L1_BI_list <- all_1100_clustered_mat_all$query_acc

all_1100_clustered_full <- filter(all_1100_clustered, query_acc %in% L1_BI_list)
length(unique(all_1100_clustered_full$library))

L1_full_list <- select(all_1100_clustered_full, pfam, contig)
L1_full_list$pfam <- c("L1_full")
L1_full_list <- L1_full_list[!duplicated(L1_full_list), ]

#get the most "leftmost" coordinate
all_1100_clustered_full_envcoord_from <- all_1100_clustered_full %>% 
  group_by(query_acc) %>% 
  slice_min(n = 1, envcoord_from) %>% 
  select(query_acc, envcoord_from)

#get the most "rightmost" coordinate
all_1100_clustered_full_envcoord_to <- all_1100_clustered_full %>% 
  group_by(query_acc) %>% 
  slice_max(n = 1, envcoord_to) %>% 
  select(query_acc, envcoord_to)

all_1100_clustered_full_envcoords <- left_join(all_1100_clustered_full_envcoord_from, all_1100_clustered_full_envcoord_to, by = "query_acc")
#calculate nucleotide positions
all_1100_clustered_full_envcoords$nuc_from <- all_1100_clustered_full_envcoords$envcoord_from*2 + (all_1100_clustered_full_envcoords$envcoord_from - 3)
all_1100_clustered_full_envcoords$nuc_to <- all_1100_clustered_full_envcoords$envcoord_to*3 - 1

write.table(select(all_1100_clustered_full_envcoords, query_acc, nuc_from, nuc_to), "outputs/all_1100_clustered_full_envcoords.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
#run

#re-structured and re-clustered all the pv sequences
###for novelty search
blastn_results <- read.table("files/blastn_hits_centroids.tsv", sep = "\t", header = F)
blastn_results$library <- gsub("Score.*", "", blastn_results$V1)
blastn_results$qcov <- abs(blastn_results$V2 - blastn_results$V3)/blastn_results$V4
hist(blastn_results$qcov)

dist <- ggplot(data = blastn_results, aes(x = qcov)) + 
  geom_density(aes(y = after_stat(count)*0.01), fill = "blue", alpha = .3, col = NA)

 dist + geom_histogram(binwidth = 0.01, alpha = .5) + theme_bw()

ggsave("outputs/qcov_distribution.png", width = 10, height =10, units = "cm", limitsize = F)


#
evals fil 
blastn_results_fil <- dplyr::filter(blastn_results, V10 < 0.00001)

hist(blastn_results_fil$qcov)
blastn_results_fil <- dplyr::filter(blastn_results_fil, qcov > 0.7)

blastn_results_fil_low_conf <- dplyr::filter(blastn_results, V10 > 0.00001)
hist(blastn_results_fil_low_conf$qcov)

what_hits <- as.data.frame(table(blastn_results$V1))


#look for those with highest % iden first in confident hits
blastn_results_fil_sliced <- blastn_results_fil %>% group_by(V1) %>% slice_max(n = 1, V9)
#under 90% iden
blastn_results_fil_sliced_90 <- filter(blastn_results_fil_sliced, V9 < 90)
blastn_results_fil_sliced_90 = blastn_results_fil_sliced_90[!duplicated(blastn_results_fil_sliced_90$V1),]

blastn_results_fil_sliced_less_90_nt <- unique(blastn_results_fil_sliced_90$V1)
# all_seq_geo_qcov_less50 <- filter(all_seq_geo_qcov_less40, !V1 %in% all_seq_geo_less_90_nt)
length(unique(blastn_results_fil_sliced_less_90_nt))



hist(blastn_results_fil_sliced_90$V9)


#create list of the ones that were hit, in general
hit_list <- unique(blastn_results_fil$V1)

#look for low conf hits that were not represented in either novel already, or other filter set
blastn_results_fil_low_conf_hits <- filter(blastn_results_fil_low_conf, !V1 %in% blastn_results_fil_sliced_less_90_nt)
blastn_results_fil_low_conf_hits_all <- filter(blastn_results_fil_low_conf_hits, !V1 %in% blastn_results_fil$V1)

blastn_results_fil_low_conf_hits_all_list <- unique(blastn_results_fil_low_conf_hits_all$V1)

length(unique(blastn_results_fil_low_conf_hits_all$V1))

#group all the ones with hits that are novel
all_hit_novel <- bind_rows(blastn_results_fil_low_conf_hits_all, blastn_results_fil_sliced_90)
length(unique(all_hit_novel$V1))

write.table(unique(all_hit_novel$V1), "outputs/all_hit_novel.txt", quote = F, col.names = F, row.names = F)



#check which didn't get hit
alignments <- read.table("files/alignment_envs.table")
length(unique(alignments$V5))
length(unique(alignments$V1))

alignments$qcov <- abs(alignments$V2 - alignments$V3)/alignments$V4
hist(alignments$qcov)

dist <- ggplot(data = alignments, aes(x = qcov)) +
  geom_density(aes(y = after_stat(count)*0.01), fill = "blue", alpha = .3, col = NA)

 dist + geom_histogram(binwidth = 0.01, alpha = .5) + theme_bw()

#qcov-only fails: high identity, good evalue, but short alignment
#these are known PVs present in SRA reads but not assembled into long enough contigs
alignments_qcov_only_fail <- alignments %>%
  filter(qcov <= 0.7, V10 < 0.00001, V9 >= 90) %>%
  group_by(V1) %>%
  slice_max(n = 1, V9) %>%
  filter(!duplicated(V1))

alignments_fil <- alignments %>%
  dplyr::filter(qcov > 0.7) %>%
  dplyr::filter(V10 < 0.00001) %>%
  dplyr::filter(V9 >= 90)

length(unique(alignments_fil$V1))
length(unique(alignments_fil$V5))

#summary of 346 NCBI PVs not recovered in Logan
ncbi_missing_summary <- read.table("outputs/ncbi_pvs_not_in_logan_summary.tsv",
                                   sep = "\t", header = T, quote = "")

sero_hits <- filter(ncbi_missing_summary, miss_reason == "zero_hits")
