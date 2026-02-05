#look for novel E6 and E7 sequences?
hmmsearch_clean <- function(input_domtbl){
  input_domtbl_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:22]
  # input_domtbl_scan <- input_domtbl_scan[!(is.na(input_domtbl_scan$V21) | input_domtbl_scan$V21=="" | !(input_domtbl_scan$V2 == "-")), ]
  
  colnames(input_domtbl_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                   "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}

#from hmmer
pv2_fasta_sto_sto <- hmmsearch_clean("pv_only_sto_sto.domtbl")
pv2_fasta_sto_sto$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", pv2_fasta_sto_sto$query_acc)
pv2_fasta_sto_sto_E6 <- filter(pv2_fasta_sto_sto, pfam == "E6")
pv2_fasta_sto_sto_E7 <- filter(pv2_fasta_sto_sto, pfam == "E7")

pv2_fasta_sto_sto_E6$cov <- abs(pv2_fasta_sto_sto_E6$hmmcoord_from - pv2_fasta_sto_sto_E6$hmmcoord_to)/110
pv2_fasta_sto_sto_E6_90cov <- filter(pv2_fasta_sto_sto_E6, cov >= .9)
pv2_fasta_sto_sto_E6_90cov <- filter(pv2_fasta_sto_sto_E6_90cov, score_one >= 50)

write.table(pv2_fasta_sto_sto_E6_90cov$query_acc, "../mikko/E6_hmm_hits.txt", row.names = F, col.names = F, quote = F)

pv2_fasta_sto_sto_E7$cov <- abs(pv2_fasta_sto_sto_E7$hmmcoord_from - pv2_fasta_sto_sto_E7$hmmcoord_to)/92
pv2_fasta_sto_sto_E7_90cov <- filter(pv2_fasta_sto_sto_E7, cov >= .9)
pv2_fasta_sto_sto_E7_90cov <- filter(pv2_fasta_sto_sto_E7_90cov, score_one >= 50)

write.table(pv2_fasta_sto_sto_E7_90cov$query_acc, "../mikko/E7_hmm_hits.txt", row.names = F, col.names = F, quote = F)


#OUUGGGH
alt_sto_hmm <- hmmsearch_clean("../mikko/pv_only_min_start_alt.domtbl")
alt_sto_hmm_E6 <- filter(alt_sto_hmm, pfam == "E6")
alt_sto_hmm_E7 <- filter(alt_sto_hmm, pfam == "E7")

alt_sto_hmm_E6$cov <- abs(alt_sto_hmm_E6$hmmcoord_from - alt_sto_hmm_E6$hmmcoord_to)/110
alt_sto_hmm_E6_90 <- filter(alt_sto_hmm_E6, cov >= .9)
alt_sto_hmm_E6_90_s50 <- filter(alt_sto_hmm_E6_90, score_one >= 50)
write.table(alt_sto_hmm_E6_90_s50$query_acc, "../mikko/alt_sto_hmm_E6_90_s50.txt", row.names = F, col.names = F, quote = F)

alt_sto_hmm_E7$cov <- abs(alt_sto_hmm_E7$hmmcoord_from - alt_sto_hmm_E7$hmmcoord_to)/92
alt_sto_hmm_E7_90 <- filter(alt_sto_hmm_E7, cov >= .9)
alt_sto_hmm_E7_90_s50 <- filter(alt_sto_hmm_E7_90, score_one >= 50)
write.table(alt_sto_hmm_E7_90_s50$query_acc, "../mikko/alt_sto_hmm_E7_90_s50.txt", row.names = F, col.names = F, quote = F)

#blast results
E6_centroids_out <- read.table("../mikko/E6_centroids_out.tsv", header = F, fill = T, sep = "\t")
E6_centroids_out_max <- E6_centroids_out %>% group_by(V1) %>% slice_max(n = 1, V10)
E6_centroids_out_max_less90 <- filter(E6_centroids_out_max, V10 < 90)
E6_centroids_out_max_less90_human <- filter(E6_centroids_out_max_less90, grepl("human|Human|Gamma",V13))
write.table(E6_centroids_out_max_less90_human$V1, "../mikko/E6_centroids_out_max_less90_human.txt", row.names = F, col.names = F, quote = F)
#lcl|ORF4_SRR12903319_424646:1567:2013 was not in the output for the blast file
E6_centroids_out_max_less90_nonhuman <- filter(E6_centroids_out_max_less90, !V1 %in% E6_centroids_out_max_less90_human$V1)
write.table(E6_centroids_out_max_less90_nonhuman$V1, "../mikko/E6_centroids_out_max_less90_nonhuman.txt", row.names = F, col.names = F, quote = F)


E7_centroids_out <- read.table("../mikko/E7_centroids_out.tsv", header = F, fill = T, sep = "\t")
E7_centroids_out_max <- E7_centroids_out %>% group_by(V1) %>% slice_max(n = 1, V10)
E7_centroids_out_max_less90 <- filter(E7_centroids_out_max, V10 < 90)
write.table(E7_centroids_out_max_less90_human$V1, "../mikko/E7_centroids_out_max_less90_human.txt", row.names = F, col.names = F, quote = F)
#none not in the blast results for E7
E7_centroids_out_max_less90_nonhuman <- filter(E7_centroids_out_max_less90, !V1 %in% E7_centroids_out_max_less90_human$V1)
write.table(E7_centroids_out_max_less90_nonhuman$V1, "../mikko/E7_centroids_out_max_less90_nonhuman.txt", row.names = F, col.names = F, quote = F)
