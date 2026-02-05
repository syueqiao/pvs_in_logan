##for novelty search
feb_7_L1_nuc_hmmer_env <- read.table("feb7_all_pr_pilot_outputs_sort_centroids_huh.tsv.csv", sep = ",", header = F, quote="", comment.char = "")
# feb_7_L1_nuc_hmmer_env$V1 <- gsub("OriginScaffoldPath=.*", "", feb_7_L1_nuc_hmmer_env$V1)
length(unique(feb_7_L1_nuc_hmmer_env_less_90_nt))
#calculate qcov
feb_7_L1_nuc_hmmer_env$qcov <- abs(feb_7_L1_nuc_hmmer_env$V2 - feb_7_L1_nuc_hmmer_env$V3)/feb_7_L1_nuc_hmmer_env$V4
#evals fil 
feb_7_L1_nuc_hmmer_env_fil <- filter(feb_7_L1_nuc_hmmer_env, V10 < 0.0001)
feb_7_L1_nuc_hmmer_env_fil_low_conf <- filter(feb_7_L1_nuc_hmmer_env, V10 > 0.0001)


#look for those with highest % iden first in confident hits
feb_7_L1_nuc_hmmer_env_sliced <- feb_7_L1_nuc_hmmer_env_fil %>% group_by(V1) %>% slice_max(n = 1, V9)
#under 90% iden
feb_7_L1_nuc_hmmer_env_sliced_90 <- filter(feb_7_L1_nuc_hmmer_env_sliced, V9 < 90)
feb_7_L1_nuc_hmmer_env_sliced_90 = feb_7_L1_nuc_hmmer_env_sliced_90[!duplicated(feb_7_L1_nuc_hmmer_env_sliced_90$V1),]

feb_7_L1_nuc_hmmer_env_less_90_nt <- unique(feb_7_L1_nuc_hmmer_env_sliced_90$V1)
# feb_7_L1_nuc_hmmer_env_qcov_less50 <- filter(feb_7_L1_nuc_hmmer_env_qcov_less40, !V1 %in% feb_7_L1_nuc_hmmer_env_less_90_nt)

write.table(feb_7_L1_nuc_hmmer_env_less_90_nt, "feb_7_L1_nuc_hmmer_env_less_90_nt.txt", quote = F, col.names = F, row.names = F)

##some are not here!
#length(unique(feb_7_L1_nuc_hmmer_env_sliced$V1))
# 859 of 942, so 83 are unaccounted for

#create list of the ones that were hit, in general
hit_list <- unique(feb_7_L1_nuc_hmmer_env_sliced$V1)
write.table(hit_list, "hit_list_pr.txt", quote = F, col.names = F, row.names = F)

#look for low conf hits that were not represented in either novel already, or other filter set
feb_7_L1_nuc_hmmer_env_fil_low_conf_hits <- filter(feb_7_L1_nuc_hmmer_env_fil_low_conf, !V1 %in% feb_7_L1_nuc_hmmer_env_less_90_nt)
feb_7_L1_nuc_hmmer_env_fil_low_conf_hits <- filter(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits, !V1 %in% feb_7_L1_nuc_hmmer_env_fil$V1)

feb_7_L1_nuc_hmmer_env_fil_low_conf_hits_list <- unique(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits$V1)

write.table(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits_list, "feb_7_L1_nuc_hmmer_env_fil_low_conf_hits_list_pr.txt", quote = F, col.names = F, row.names = F)
