library(tidyverse)

hmmsearch_clean_info <- function(input_domtbl){
  input_domtbl_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:23]
  input_domtbl_scan <- input_domtbl_scan[!(is.na(input_domtbl_scan$V21) | input_domtbl_scan$V21=="" | !(input_domtbl_scan$V2 == "-")), ]
  
  colnames(input_domtbl_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                   "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc", "misc2")
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}

path_1k <- hmmsearch_clean_info("../petase_project/all_seqs_1k.domtbl")
path_1k$query_acc <- sub('_','\\|',path_1k$query_acc)
path_1k$pfam <- factor(x = path_1k$pfam,
                                   levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH"))
path_1k_tab <- path_1k %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

path_1k_tab_full <- path_1k_tab %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_GH, ~ .x > 0))

path_1k_tab_logi <- path_1k_tab
path_1k_tab_logi[-1] <- as.data.frame(+ sapply(path_1k_tab_logi[-1], as.logical))

path_1k_tab_logi$sum_p_a <- rowSums(path_1k_tab_logi[,2:7]) 

ggplot(path_1k_tab_logi, aes(x=sum_p_a)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

path_1k_tab_logi <- separate(path_1k_tab_logi, col=query_acc, into=c('lib','score', 'edges', 'pos', 'aln', 'path', 'path2'), sep="\\|")
path_1k_tab_logi_no_singles <- filter(path_1k_tab_logi, grepl("_",edges))
path_1k_tab_logi_no_singles_full <- path_1k_tab_logi_no_singles %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_GH, ~ .x > 0))

ggplot(path_1k_tab_logi_no_singles, aes(x=sum_p_a)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

ggplot(path_1k_tab_logi_no_singles_full, aes(x=sum_p_a)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

path_1k_tab_logi_no_singles_score_fil <- separate(path_1k_tab_logi_no_singles, col=score, into=c('text','score'), sep="\\=")
path_1k_tab_logi_no_singles_score_fil$score <- as.numeric(path_1k_tab_logi_no_singles_score_fil$score)
path_1k_tab_logi_no_singles_score_fil <- filter(path_1k_tab_logi_no_singles_score_fil, score > 450)

lib_table <- as.data.frame(table(path_1k_tab_logi_no_singles_full$lib))


pv2_L1_sto_sto_SRR7992617 <- pv2_L1_sto_sto[grepl("SRR7992617_", pv2_L1_sto_sto$query_acc),]
write.table(pv2_L1_sto_sto_SRR7992617$query_acc, "pv2_L1_sto_sto_SRR7992617_query_acc.txt", col.names = F, row.names = F, quote = F, sep = "\t")


concat_pv2_fil_pv_SRR7992617 <- concat_pv2_fil_pv[grepl("SRR7992617_", concat_pv2_fil_pv$V1),]
concat_pv2_fil_pv_SRR7992617_L1 <- concat_pv2_fil_pv_SRR7992617[grepl("L1.", concat_pv2_fil_pv_SRR7992617$V6),]
write.table(concat_pv2_fil_pv_SRR7992617_L1$V6, "concat_pv2_fil_pv_SRR7992617_L1_db_hits.txt", col.names = F, row.names = F, quote = F, sep = "\t")
