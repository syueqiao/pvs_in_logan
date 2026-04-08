library(tidyverse)
library(ggplot2)
library(viridisLite)
library(viridis)

hmmsearch_clean <- function(input_domtbl){
  input_domtbl_scan <- read.table(input_domtbl, sep = "",  header = FALSE, fill = TRUE) %>% na.omit() %>% .[,1:22]
  # input_domtbl_scan <- input_domtbl_scan[!(is.na(input_domtbl_scan$V21) | input_domtbl_scan$V21=="" | !(input_domtbl_scan$V2 == "-")), ]
  
  colnames(input_domtbl_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                   "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}

#script to help "choose" a set of 10 contigs for testing various graph-based contig retrieval methods
pv2_L1_sto_sto <- hmmsearch_clean("files/pv_all_L1.domtbl")
pv2_L1_sto_sto$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", pv2_L1_sto_sto$query_acc)

 
pv2_L1_sto_sto_uniq <- unique(pv2_L1_sto_sto)

pv2_L1_sto_sto_uniq$library <- sub("_.*", "", pv2_L1_sto_sto_uniq$contig)

###########

pv2_L1_sto_sto_uniq$pfam <- factor(x = pv2_L1_sto_sto_uniq$pfam,
                                   levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH"))

pv2_L1_sto_sto_uniq_tab <- pv2_L1_sto_sto_uniq %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

pv2_L1_sto_sto_uniq_tab_full <- pv2_L1_sto_sto_uniq_tab %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_GH, ~ .x > 0))

pv2_L1_sto_sto_uniq_partial <- filter(pv2_L1_sto_sto_uniq, !query_acc %in% pv2_L1_sto_sto_uniq_tab_full$query_acc)

# write.table(unique(pv2_L1_sto_sto_uniq$library), "pv2_L1_sto_sto_libs.txt", quote = FALSE, col.names = FALSE, row.names = F)

pv2_L1_sto_sto_uniq_partial_tab <- pv2_L1_sto_sto_uniq_partial %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

pv2_L1_sto_sto_uniq_partial_tab_100 <- pv2_L1_sto_sto_uniq_partial_tab[sample(nrow(pv2_L1_sto_sto_uniq_partial_tab), 1000), ]
pv2_L1_sto_sto_uniq_partial_tab_100$library <- gsub("_.*", "", pv2_L1_sto_sto_uniq_partial_tab_100$query_acc)

pv2_L1_sto_sto_uniq_partial_tab_100_sum <- select(pv2_L1_sto_sto_uniq_partial_tab_100, library, ncbi_and_pave_L1_B,
ncbi_and_pave_L1_CD, ncbi_and_pave_L1_E, ncbi_and_pave_L1_GH, ncbi_and_pave_L1_I) %>% 
    as.data.frame() %>%
    group_by(library) %>% 
    summarize_all(sum) %>%
    as.data.frame()

table(pv2_L1_sto_sto_uniq_partial_tab_100$library)

pv2_L1_sto_sto_uniq_partial_lib_tab <- pv2_L1_sto_sto_uniq_partial %>%
  count(library, pfam) %>%
  spread(pfam, n, fill = 0)

pv2_L1_sto_sto_uniq_tab$library <- sub("_.*", "", pv2_L1_sto_sto_uniq_tab$query_acc)

concat_pv2_fil_pv_FP_L1_max <- concat_pv2_fil_pv_FP_L1 %>% group_by(library) %>% slice_max(n = 1, per_iden)

##analyse output
hmmsearch_clean <- function(input_domtbl){
  input_domtbl_scan <- read.table(input_domtbl, sep = "",  header = FALSE, fill = TRUE) %>% na.omit() %>% .[,1:22]
  input_domtbl_scan <- input_domtbl_scan[!(is.na(input_domtbl_scan$V21) | input_domtbl_scan$V21=="" | !(input_domtbl_scan$V2 == "-")), ]
  
  colnames(input_domtbl_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                   "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}

mega_results <- hmmsearch_clean("files/mega_results.domtbl")
# mega_results <- separate(mega_results,query_acc,into = c("method","query_acc"),sep = "\\.",remove = FALSE,extra = "merge")
mega_results <- unique(mega_results)
mega_results_eval <- filter(mega_results, mega_results$i_eval < 0.00001)
# mega_results$library <- sub("_.*", "", mega_results$contig)

mega_results_eval$pfam <- factor(x = mega_results_eval$pfam,
                                   levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH"))


mega_results_eval_tab <- mega_results_eval %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

mega_results_eval_tab$ncbi_and_pave_L1_CD <- as.numeric(sub("2", "1", mega_results_eval_tab$ncbi_and_pave_L1_CD))
mega_results_eval_tab$ncbi_and_pave_L1_E <-  as.numeric(sub("2", "1", mega_results_eval_tab$ncbi_and_pave_L1_E))

mega_results_eval_tab <- separate(mega_results_eval_tab,query_acc,into = c("method","query_acc"),sep = "\\.",remove = FALSE,extra = "merge")
mega_results_eval_tab$sum <- rowSums(mega_results_eval_tab[, 3:8])

ggplot(mega_results_eval_tab, aes(x=sum, fill=method)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

####remove pcr amp libraries
library_metadata <- read.table("files/L1_library_metadata.txt", fill = TRUE, sep = "\t", quote = "")
#do some counting

library_amplicons <- dplyr::filter(library_metadata, V2 == "AMPLICON")
library_amplicons <- dplyr::filter(library_amplicons, V9 == "PCR")


pv2_L1_sto_sto_uniq_metadata <- left_join(pv2_L1_sto_sto_uniq, library_metadata, by = c("library" = "V1"))
# length(unique(pv2_L1_sto_sto_uniq_metadata$library))
pv2_L1_sto_sto_uniq_metadata_noamp <- filter(pv2_L1_sto_sto_uniq_metadata, !V2 == "AMPLICON")
pv2_L1_sto_sto_uniq_metadata_noamp <- filter(pv2_L1_sto_sto_uniq_metadata_noamp, !V9 == "PCR")

pv2_L1_sto_sto_uniq_metadata_noamp_RNA <-  filter(pv2_L1_sto_sto_uniq_metadata_noamp, V10 == "TRANSCRIPTOMIC")

pv2_L1_sto_sto_uniq_metadata_noamp_RNA$pfam <- factor(x = pv2_L1_sto_sto_uniq_metadata_noamp_RNA$pfam,
                                                  levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH"))

pv2_L1_sto_sto_uniq_metadata_noamp_RNA_tab <- pv2_L1_sto_sto_uniq_metadata_noamp_RNA %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)


pv2_L1_sto_sto_uniq_metadata_noamp_RNA_tab_full <- pv2_L1_sto_sto_uniq_metadata_noamp_RNA_tab %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_GH, ~ .x > 0))



length(unique(pv2_L1_sto_sto_uniq_metadata_noamp$library))

pv2_L1_sto_sto_uniq_metadata_noamp$pfam <- factor(x = pv2_L1_sto_sto_uniq_metadata_noamp$pfam,
                                   levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH"))

pv2_L1_sto_sto_uniq_metadata_noamp_tab <- pv2_L1_sto_sto_uniq_metadata_noamp %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

pv2_L1_sto_sto_uniq_metadata_noamp_tab_full <- pv2_L1_sto_sto_uniq_metadata_noamp_tab %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_GH, ~ .x > 0))

pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial <- filter(pv2_L1_sto_sto_uniq_metadata_noamp_tab, !query_acc %in% pv2_L1_sto_sto_uniq_metadata_noamp_tab_full$query_acc)
pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial$library <- sub("_.*", "", pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial$query_acc)

pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial <- subset(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial, grepl("_L", query_acc))
pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_num <- pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial
pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_num <- pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_num[,2:8]

pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num <- aggregate(. ~ library, pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_num, function(x) sum(as.numeric(x)))
pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num$sums <- rowSums(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num[,2:7])
pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_4 <- filter(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num, sums > 3)

ggplot(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_4, aes(x=sums)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0 <- filter(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num, sums < 30 )
ggplot(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0, aes(x=sums)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0 <- filter(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0, sums >= 3)

logical_df <- pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0
logical_df[-1] <- as.data.frame(+ sapply(logical_df[-1], as.logical))
logical_df$pa_sums <- rowSums(logical_df[,2:7])
logical_df_3 <- filter(logical_df, pa_sums >= 3)

pr_sample <- pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0 %>%
  group_by(sums) %>%
  sample_n(3, replace = TRUE)

logical_df_3_sample_10 <- logical_df_3 %>%
  group_by(pa_sums) %>%
  sample_n(1, replace = TRUE)

logical_df_3_pa_whole <- left_join(logical_df_3, pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0, by = "library")
logical_df_3_pa_whole <- logical_df_3_pa_whole[c(1:7, 9, 16)]
names(logical_df_3_pa_whole)[names(logical_df_3_pa_whole) == 'sums.y'] <- 'overall_hits'
logical_df_3_pa_whole_smaller <- filter(logical_df_3_pa_whole, overall_hits < 9)

logical_df_3_pa_whole %>%
  ggplot(aes(x=pa_sums, y=overall_hits)) +
  geom_point(size = 3, alpha = 0.2) +
  theme_bw()

ggplot(logical_df_3_pa_whole_smaller, aes(x=overall_hits)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

logical_df_3_pa_whole_smaller_sample10 <- logical_df_3_pa_whole_smaller %>%
  group_by(overall_hits) %>%
  sample_n(1, replace = TRUE)

write.table(pr_sample$library, "outputs/pr_sample4.txt", quote = FALSE, col.names = FALSE, row.names = F)
write.table(logical_df_3_pa_whole_smaller_sample10$library, "outputs/logical_df_3_pa_whole_smaller_sample10.txt", quote = FALSE, col.names = FALSE, row.names = F)
write.table(logical_df_3_pa_whole_smaller$library, "outputs/15k_library_for_pathracer_cutoff_9.txt", quote = FALSE, col.names = FALSE, row.names = F)

#look to see if there is pattern re:those that fail
failed <- read.table("files/acc_fail.txt", fill = TRUE, sep = "\t", quote = "")

logical_df_3_pa_whole_smaller_failed <- filter(logical_df_3_pa_whole_smaller, library %in% failed$V1)
hist(logical_df_3_pa_whole_smaller_failed$overall_hits)
threes <- filter(logical_df_3_pa_whole_smaller_failed, pa_sums =="3")
threes_in <- filter(logical_df_3_pa_whole_smaller, pa_sums =="3")
threes_in <- filter(threes_in, !library %in% threes$library)
threes

path_1k_tab_logi_no_singles_full_three <- filter(path_1k_tab_logi_no_singles_full, lib %in% threes_in$library)
       
#graph extract test
pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_ge <- pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial
lib_sums <- pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0[c(1, 8)]

pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_ge <- left_join(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_ge, lib_sums, by = "library")

write_file_func <- function(lib){
  library_filtered <- filter(pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_ge, library == lib)
  write.table(library_filtered$query_acc, paste("outputs/subgraph/", lib, "_list.txt", sep = ""), quote = FALSE, col.names = FALSE, row.names = F)
}

write_file_func("SRR25663671")

graph_ex_sample <- pr_sample$library

for (i in graph_ex_sample){
  write_file_func(i)
}

write.table(graph_ex_sample, "outputs/graph_ex_sample.txt", quote = FALSE, col.names = FALSE, row.names = F)

graph_18k <- pv2_L1_sto_sto_uniq_metadata_noamp_tab_partial_lib_num_30_0$library
for (i in graph_18k){
  write_file_func(i)
}

write.table(graph_18k, "outputs/graph_18k.txt", quote = FALSE, col.names = FALSE, row.names = F)

#empty this out 
