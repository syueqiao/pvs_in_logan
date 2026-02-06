#use output from pro_file_parsing output
july_3400_test <- read.table("july1_3400_run.pro", sep = "\t", fill = TRUE)
eval_fil_july_3400_test <- filter(july_3400_test, V11 < 0.000000001 )
jul1_test <- eval_fil_july_3400_test

#filter for pv only
##########CHANGE HERE FOR DIFFERENT INPUT#####################
pv2_test <- filter(july_3400_test, grepl("papilloma", V6))
#get clean acc
pv2_test$V14 <- sub("_.*", "", pv2_test$V1)

#get e score and choose "lowest"
pv2_test$e_score <- -log10(pv2_test$V11)

pv2_test_sliced <- pv2_test %>% group_by(V1) %>% slice_max(order_by = e_score, n = 1)
pv2_test_sliced <- pv2_test_sliced %>% group_by(V1) %>% slice_max(order_by = V10, n = 1)


pv2_test_sliced$q_cov <- abs(pv2_test_sliced$V2 - pv2_test_sliced$V3)/pv2_test_sliced$V4
pv2_test_sliced$m_cov <- abs(pv2_test_sliced$V7 - pv2_test_sliced$V8)/pv2_test_sliced$V9
pv2_test_sliced$m_cov_aa <- abs(pv2_test_sliced$V7 - pv2_test_sliced$V8)



ggplot(pv2_test_sliced, aes(x=V4)) +
  geom_histogram(, binwidth = 50, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Contig length dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(pv2_test_sliced$V4), by = 500),1))

ggplot(pv2_test_sliced, aes(x=e_score)) +
  geom_histogram(, binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("E-score dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(350), by = 20),1))

ggplot(pv2_test_sliced, aes(x=V10)) +
  geom_histogram(, binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("P-iden dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 5),1))



mean(pv2_test_sliced$V4)
  
  
ggplot(pv2_test_sliced, aes(x=m_cov, y=e_score, colour = V10)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "D")


pv2_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov, FUN = median), y=m_cov, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pv2_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov_aa, FUN = median), y=m_cov_aa, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage in aas") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pv2_test_sliced_AB331651 <- filter(pv2_test_sliced, V6 == "papilloma.PPV_E1_C.Bos_taurus_papillomavirus_10:AB331651")
write.table(pv2_test_sliced_AB331651$V1, "pv2_test_sliced_AB331651.txt", row.names = FALSE, col.names = FALSE)

pv2_test_sliced_MZ189241 <- filter(pv2_test_sliced, V6 == "papilloma.PPV_E1_DBD.Human_papillomavirus_type_18:MZ189241")
write.table(pv2_test_sliced_MZ189241$V1, "pv2_test_sliced_MZ189241", row.names = FALSE, col.names = FALSE)

write.table(paste0(pv2_test_sliced$V14, collapse = ","), "pv2_test_sliced_SRA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
names(pv2_test_sliced)[names(pv2_test_sliced) == 'V14'] <- 'Run'

sra_metadata <- read.table("3400_sra_metadata.txt", sep = ",", fill = TRUE, header = TRUE)
pv2_test_sliced_metadata <- left_join(pv2_test_sliced, sra_metadata, by = c("V14" = "Run"))
pv2_test_sliced_metadata <- pv2_test_sliced_metadata[,1:33]

pv2_test_sliced_metadata$type <- sub("[^.]+\\.([^.]+)\\..*", "\\1", pv2_test_sliced_metadata$V6)


pv2_test_sliced_metadata$type_grouped <-gsub("_N","",as.character(pv2_test_sliced_metadata$type))
pv2_test_sliced_metadata$type_grouped <-gsub("_C","",as.character(pv2_test_sliced_metadata$type_grouped))
pv2_test_sliced_metadata$type_grouped <-gsub("_DBD","",as.character(pv2_test_sliced_metadata$type_grouped))


ggplot(pv2_test_sliced_metadata, aes(fill=type, x = LibrarySource)) + 
  geom_bar() +
  scale_fill_viridis_d() +
  theme_bw()

############SEPARATE FOR ANELLO###################
july_3400_test_anello <- read.table("july1_3400_run.pro", sep = "\t", fill = TRUE)
anello1_test <- filter(july_3400_test_anello, grepl("anello", V6))
anello1_test <- filter(anello1_test, V11 < 0.0001 )

#get clean acc
anello1_test$V14 <- sub("_.*", "", anello1_test$V1)

#get e score and choose "lowest"
anello1_test$e_score <- -log10(anello1_test$V11)

anello1_test_sliced <- anello1_test %>% group_by(V1) %>% slice_max(order_by = e_score, n = 1)



anello1_test_sliced$q_cov <- abs(anello1_test_sliced$V2 - anello1_test_sliced$V3)/anello1_test_sliced$V4
anello1_test_sliced$m_cov <- abs(anello1_test_sliced$V7 - anello1_test_sliced$V8)/anello1_test_sliced$V9
anello1_test_sliced$m_cov_aa <- abs(anello1_test_sliced$V7 - anello1_test_sliced$V8)

anello1_test_sliced <- anello1_test_sliced %>% group_by(V1) %>% slice_max(order_by = m_cov, n = 1)
anello1_test_sliced <- anello1_test_sliced %>% group_by(V1) %>% slice_max(order_by = V10, n = 1)



ggplot(anello1_test_sliced, aes(x=V4)) +
  geom_histogram(, binwidth = 10, fill="#FF7354", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Contig length dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(anello1_test_sliced$V4), by = 500),1))

ggplot(anello1_test_sliced, aes(x=e_score)) +
  geom_histogram(, binwidth = 5, fill="#FF7354", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("E-score dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(350), by = 20),1))

ggplot(anello1_test_sliced, aes(x=V10)) +
  geom_histogram(, binwidth = 1, fill="#FF7354", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("P-iden dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 5),1))



mean(anello1_test_sliced$V4)


ggplot(anello1_test_sliced, aes(x=m_cov, y=e_score, colour = V10)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "C")


anello1_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov, FUN = median), y=m_cov, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c(option = "C") +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

anello1_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov_aa, FUN = median), y=m_cov_aa, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage in aas") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

anello1_test_sliced_AB331651 <- filter(anello1_test_sliced, V6 == "papilloma.PPV_E1_C.Bos_taurus_papillomavirus_10:AB331651")
# write.table(anello1_test_sliced_AB331651$V1, "anello1_test_sliced_AB331651.txt", row.names = FALSE, col.names = FALSE)

anello1_test_sliced_MZ189241 <- filter(anello1_test_sliced, V6 == "papilloma.PPV_E1_DBD.Human_papillomavirus_type_18:MZ189241")
# write.table(anello1_test_sliced_MZ189241$V1, "anello1_test_sliced_MZ189241", row.names = FALSE, col.names = FALSE)

names(anello1_test_sliced)[names(anello1_test_sliced) == 'V14'] <- 'Run'

write.table(paste0(anello1_test_sliced$Run, collapse = ","), "anello1_test_sliced_SRA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


b3400 <- read.table("3400_sra_metadata.txt", sep = ",", fill = TRUE, header = TRUE)
anello1_test_sliced_metadata <- left_join(anello1_test_sliced, sra_metadata, by = "Run")
anello1_test_sliced_metadata <- anello1_test_sliced_metadata[,1:33]

anello1_test_sliced_metadata$type <- sub("[^.]+\\.([^.]+)\\..*", "\\1", anello1_test_sliced_metadata$V6)


anello1_test_sliced_metadata$type_grouped <-gsub("_N","",as.character(anello1_test_sliced_metadata$type))
anello1_test_sliced_metadata$type_grouped <-gsub("_C","",as.character(anello1_test_sliced_metadata$type_grouped))
anello1_test_sliced_metadata$type_grouped <-gsub("_DBD","",as.character(anello1_test_sliced_metadata$type_grouped))


ggplot(anello1_test_sliced_metadata, aes(fill=type_grouped, x=LibrarySource)) + 
  geom_bar() +
  scale_fill_viridis_d() +
  theme_bw()


anello1_test_sliced_noplasmo <- filter(anello1_test_sliced, V14 != "SRR1567977")
anello1_test_sliced_noplasmo_samp <- anello1_test_sliced_noplasmo[sample(nrow(anello1_test_sliced_noplasmo), 30), ]

write.table(anello1_test_sliced_noplasmo_samp$V1, "anello1_test_sliced_noplasmo.txt", row.names = FALSE, col.names = FALSE)

##########JUL 5 run##################################################################################################################################
july5_3400_test <- read.table("july5_3400_run.pro", sep = "\t", fill = TRUE)
eval_fil_july5_3400_test <- filter(july5_3400_test, V11 < 0.000000001 )
eval_fil_july5_3400_test$V14 <- sub("_.*", "", eval_fil_july5_3400_test$V1)


############processing for "pv3"
pv3_test <- filter(eval_fil_july5_3400_test, grepl("papilloma", V6))
#get clean acc

#get e score and choose "lowest"
pv3_test$e_score <- -log10(pv3_test$V11)

pv3_test_sliced <- pv3_test %>% group_by(V1) %>% slice_max(order_by = e_score, n = 1)
pv3_test_sliced <- pv3_test_sliced %>% group_by(V1) %>% slice_max(order_by = V10, n = 1)


pv3_test_sliced$q_cov <- abs(pv3_test_sliced$V2 - pv3_test_sliced$V3)/pv3_test_sliced$V4
pv3_test_sliced$m_cov <- abs(pv3_test_sliced$V7 - pv3_test_sliced$V8)/pv3_test_sliced$V9
pv3_test_sliced$m_cov_aa <- abs(pv3_test_sliced$V7 - pv3_test_sliced$V8)



ggplot(pv3_test_sliced, aes(x=V4)) +
  geom_histogram(, binwidth = 50, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Contig length dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(pv3_test_sliced$V4), by = 500),1))

ggplot(pv3_test_sliced, aes(x=e_score)) +
  geom_histogram(, binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("E-score dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(350), by = 20),1))

ggplot(pv3_test_sliced, aes(x=V10)) +
  geom_histogram(, binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("P-iden dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 5),1))



mean(pv3_test_sliced$V4)


ggplot(pv3_test_sliced, aes(x=m_cov, y=e_score, colour = V10)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "D")


pv3_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov, FUN = median), y=m_cov, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pv3_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov_aa, FUN = median), y=m_cov_aa, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage in aas") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pv3_test_sliced_AB331651 <- filter(pv3_test_sliced, V6 == "papilloma.PPV_E1_C.Bos_taurus_papillomavirus_10:AB331651")
write.table(pv3_test_sliced_AB331651$V1, "pv3_test_sliced_AB331651.txt", row.names = FALSE, col.names = FALSE)

pv3_test_sliced_MZ189241 <- filter(pv3_test_sliced, V6 == "papilloma.PPV_E1_DBD.Human_papillomavirus_type_18:MZ189241")
write.table(pv3_test_sliced_MZ189241$V1, "pv3_test_sliced_MZ189241", row.names = FALSE, col.names = FALSE)

write.table(paste0(pv3_test_sliced$V14, collapse = ","), "pv3_test_sliced_SRA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
names(pv3_test_sliced)[names(pv3_test_sliced) == 'V14'] <- 'Run'

sra_metadata <- read.table("3400_sra_metadata.txt", sep = ",", fill = TRUE, header = TRUE)
pv3_test_sliced_metadata <- left_join(pv3_test_sliced, sra_metadata, by = "Run")
pv3_test_sliced_metadata <- pv3_test_sliced_metadata[,1:33]

pv3_test_sliced_metadata$type <- sub("[^.]+\\.([^.]+)\\..*", "\\1", pv3_test_sliced_metadata$V6)


pv3_test_sliced_metadata$type_grouped <-gsub("_N","",as.character(pv3_test_sliced_metadata$type))
pv3_test_sliced_metadata$type_grouped <-gsub("_C","",as.character(pv3_test_sliced_metadata$type_grouped))
pv3_test_sliced_metadata$type_grouped <-gsub("_DBD","",as.character(pv3_test_sliced_metadata$type_grouped))


ggplot(pv3_test_sliced_metadata, aes(fill=type_grouped, x=LibrarySource)) + 
  geom_bar() +
  scale_fill_viridis_d() +
  theme_bw()

############SEPARATE FOR ANELLO###################
july5_3400_test_anello <- read.table("july5_3400_run.pro", sep = "\t", fill = TRUE)
anello2_test <- filter(july5_3400_test_anello, grepl("anello", V6))
anello2_test <- filter(anello2_test, V11 < 0.000000001 )

#get clean acc
anello2_test$V14 <- sub("_.*", "", anello2_test$V1)

#get e score and choose "lowest"
anello2_test$e_score <- -log10(anello2_test$V11)

anello2_test_sliced <- anello2_test %>% group_by(V1) %>% slice_max(order_by = e_score, n = 1)



anello2_test_sliced$q_cov <- abs(anello2_test_sliced$V2 - anello2_test_sliced$V3)/anello2_test_sliced$V4
anello2_test_sliced$m_cov <- abs(anello2_test_sliced$V7 - anello2_test_sliced$V8)/anello2_test_sliced$V9
anello2_test_sliced$m_cov_aa <- abs(anello2_test_sliced$V7 - anello2_test_sliced$V8)

anello2_test_sliced <- anello2_test_sliced %>% group_by(V1) %>% slice_max(order_by = m_cov, n = 1)
anello2_test_sliced <- anello2_test_sliced %>% group_by(V1) %>% slice_max(order_by = V10, n = 1)



ggplot(anello2_test_sliced, aes(x=V4)) +
  geom_histogram(, binwidth = 10, fill="#FF7354", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Contig length dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(anello2_test_sliced$V4), by = 500),1))

ggplot(anello2_test_sliced, aes(x=e_score)) +
  geom_histogram(, binwidth = 5, fill="#FF7354", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("E-score dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(350), by = 20),1))

ggplot(anello2_test_sliced, aes(x=V10)) +
  geom_histogram(, binwidth = 1, fill="#FF7354", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("P-iden dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 5),1))



mean(anello2_test_sliced$V4)


ggplot(anello2_test_sliced, aes(x=m_cov, y=e_score, colour = V10)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "C")


anello2_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov, FUN = median), y=m_cov, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c(option = "C") +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

anello2_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov_aa, FUN = median), y=m_cov_aa, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage in aas") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

anello2_test_sliced_AB331651 <- filter(anello2_test_sliced, V6 == "papilloma.PPV_E1_C.Bos_taurus_papillomavirus_10:AB331651")
# write.table(anello2_test_sliced_AB331651$V1, "anello2_test_sliced_AB331651.txt", row.names = FALSE, col.names = FALSE)

anello2_test_sliced_MZ189241 <- filter(anello2_test_sliced, V6 == "papilloma.PPV_E1_DBD.Human_papillomavirus_type_18:MZ189241")
# write.table(anello2_test_sliced_MZ189241$V1, "anello2_test_sliced_MZ189241", row.names = FALSE, col.names = FALSE)

names(anello2_test_sliced)[names(anello2_test_sliced) == 'V14'] <- 'Run'

write.table(paste0(anello2_test_sliced$Run, collapse = ","), "anello2_test_sliced_SRA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


sra_metadata <- read.table("3400_sra_metadata.txt", sep = ",", fill = TRUE, header = TRUE)
anello2_test_sliced_metadata <- left_join(anello2_test_sliced, sra_metadata, by = "Run")
anello2_test_sliced_metadata <- anello2_test_sliced_metadata[,1:33]

anello2_test_sliced_metadata$type <- sub("[^.]+\\.([^.]+)\\..*", "\\1", anello2_test_sliced_metadata$V6)


anello2_test_sliced_metadata$type_grouped <-gsub("_N","",as.character(anello2_test_sliced_metadata$type))
anello2_test_sliced_metadata$type_grouped <-gsub("_C","",as.character(anello2_test_sliced_metadata$type_grouped))
anello2_test_sliced_metadata$type_grouped <-gsub("_DBD","",as.character(anello2_test_sliced_metadata$type_grouped))


ggplot(anello2_test_sliced_metadata, aes(fill=type_grouped, x=LibrarySource)) + 
  geom_bar() +
  scale_fill_viridis_d() +
  theme_bw()


anello2_test_sliced_noplasmo <- filter(anello2_test_sliced, V14 != "SRR1567977")
anello2_test_sliced_noplasmo_samp <- anello2_test_sliced_noplasmo[sample(nrow(anello2_test_sliced_noplasmo), 30), ]

write.table(anello2_test_sliced_noplasmo_samp$V1, "anello2_test_sliced_noplasmo.txt", row.names = FALSE, col.names = FALSE)

####################for logan run##############################
concat_pv2_fil_pv <- read.table("pv_only_evail_fil.pro", sep = "\t", fill = TRUE, header = FALSE)

#another FP that wasnt caught
concat_pv2_fil_pv_FP <- filter(concat_pv2_fil_pv, V6 != "papilloma.Late_protein_L2.Human_papillomavirus:LQ465342")
colnames(concat_pv2_fil_pv_FP) <- c("contig", "start", "end", "length", "strand", "model_name", "model_start", "model_end", "model_length", "per_iden", "e_val", "cigar", "aa", "nuc")
length(unique(concat_pv2_fil_pv_FP$contig))
concat_pv2_fil_pv_FP$library <- sub("_.*", "", concat_pv2_fil_pv_FP$contig)
length(unique(concat_pv2_fil_pv_FP$library))
write.table(concat_pv2_fil_pv_FP$library, "all_library_sources.txt", row.names = FALSE, col.names = FALSE, quote = TRUE)

a1 <- mean(concat_pv2_fil_pv_FP$length)
b1 <- median(concat_pv2_fil_pv_FP$length)
c1 <- sd(concat_pv2_fil_pv_FP$length)

concat_pv2_fil_pv_FP_L1 <- filter(concat_pv2_fil_pv_FP, grepl("L1",model_name))
length(unique(concat_pv2_fil_pv_FP$contig))

a <- select(concat_pv2_fil_pv_FP_L1, contig, length)
a$Set <- c("L1")

b <- select(concat_pv2_fil_pv_FP, contig, length)
b$Set <- c("all")


ab <- rbind(a, b)
colorsab <- c(all = "#83CBEB", L1 = "#104862")

p1 <- ggplot(ab, aes(x=length, fill=Set)) +
  geom_histogram(binwidth = 0.05, alpha=0.5, color="#e9ecef") +
  geom_density(aes(y=0.05*..count.., alpha=0.5)) +
  scale_x_log10() +
  annotate("text", x = 1000, y = 200000, label = paste0("Mean = ", round(a1, digits = 3), "\n", "Median = ", round(b1, digits = 3))) +
  theme_bw() +
  scale_fill_manual(values=colorsab) +
  ggtitle("Contig Length Distribution") +
  xlab("Contig Length (log scale)") +
  ylab("Count")
  

p1 + annotation_logticks(sides = "b")

length(unique(concat_pv2_fil_pv_FP$contig))

mean(concat_pv2_fil_pv_FP$V4)
median(concat_pv2_fil_pv_FP$V4)
sd(concat_pv2_fil_pv_FP$V4)


concat_pv2_fil_pv_FP$e_score <- -log(concat_pv2_fil_pv_FP$e_val)

a2 <- mean(concat_pv2_fil_pv_FP$e_score)
b2 <- median(concat_pv2_fil_pv_FP$e_score)
c2 <- sd(concat_pv2_fil_pv_FP$e_score)

ggplot(concat_pv2_fil_pv_FP, aes(x=e_score)) +
  geom_histogram(, binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=5*..count..)) +
  annotate("text", x = 100, y = 100000, label = paste0("median = ", round(b2, digits = 3))) +
  theme_bw() +
  ggtitle("E-score (-log10[eval]) dist.") +
  ylab("count")

a3 <- mean(concat_pv2_fil_pv_FP$per_iden)
b3 <- median(concat_pv2_fil_pv_FP$per_iden)
c3 <- sd(concat_pv2_fil_pv_FP$per_iden)

ggplot(concat_pv2_fil_pv_FP, aes(x=per_iden)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  annotate("text", x = 25, y = 1000000, label = paste0("mean = ", round(a3, digits = 3), "\n", "median = ", round(b3, digits = 3), "\n", "sd = ", round(c3, digits = 3))) +
  theme_bw() +
  ggtitle("P-iden dist.")

pv2_fasta <- as.data.frame(cbind(concat_pv2_fil_pv_FP$contig, concat_pv2_fil_pv_FP$nuc))
pv2_fasta <- as.data.frame(c(rbind(concat_pv2_fil_pv_FP$contig,concat_pv2_fil_pv_FP$nuc)))
pv2_fasta <- as.data.frame(matrix(rbind(concat_pv2_fil_pv_FP$contig,concat_pv2_fil_pv_FP$nuc),ncol=1))

pv2_fasta <- pv2_fasta %>%
  mutate(V1 = if_else(row_number() %% 2 == 1, paste0(">", V1), V1))
write.table(pv2_fasta, "pv2_fasta.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)


pv2_test_sliced$m_cov <- abs(pv2_test_sliced$V7 - pv2_test_sliced$V8)/pv2_test_sliced$V9

ggplot(concat_pv2_fil_pv_FP, aes(x=m_cov, y=per_iden, colour = e_score)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "D")

concat_pv2_fil_pv_FP <- sub("_.*", "", concat_pv2_fil_pv_FP$contig)


orf_stats <- read.table("pv_only_sto_sto_lengths.txt")

# 
# 
# pv2_test_sliced %>%
#   ggplot(aes(x=reorder(V6, m_cov, FUN = median), y=m_cov, color=log10(e_score))) +
#   geom_point(size = 3, alpha = 0.2) +
#   scale_color_viridis_c() +
#   # geom_jitter(color="black", size=0.4, alpha=0.9) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("model vs. coverage") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# pv2_test_sliced %>%
#   ggplot(aes(x=reorder(V6, m_cov_aa, FUN = median), y=m_cov_aa, color=log10(e_score))) +
#   geom_point(size = 3, alpha = 0.5) +
#   scale_color_viridis_c() +
#   # geom_jitter(color="black", size=0.4, alpha=0.9) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("model vs. coverage in aas") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##################################
library(tidyverse)
#visualize and look for summary stats
july_3400_test <- read.table("july1_3400_run.pro", sep = "\t", fill = TRUE)
eval_fil_july_3400_test <- filter(july_3400_test, V11 < 0.000000001 )
jul1_test <- eval_fil_july_3400_test
july_3400_test <- read.table("july1_3400_run.pro", sep = "\t", fill = TRUE)

#filter for pv only
pv2_test <- filter(jul1_test, grepl("papilloma", V6))
#get clean acc
pv2_test$V14 <- sub("_.*", "", pv2_test$V1)

#get e score and choose "lowest"
pv2_test$e_score <- -log10(pv2_test$V11)

pv2_test_sliced <- pv2_test %>% group_by(V1) %>% slice_max(order_by = e_score, n = 1)
pv2_test_sliced <- pv2_test_sliced %>% group_by(V1) %>% slice_max(order_by = V10, n = 1)


pv2_test_sliced$q_cov <- abs(pv2_test_sliced$V2 - pv2_test_sliced$V3)/pv2_test_sliced$V4
pv2_test_sliced$m_cov <- abs(pv2_test_sliced$V7 - pv2_test_sliced$V8)/pv2_test_sliced$V9
pv2_test_sliced$m_cov_aa <- abs(pv2_test_sliced$V7 - pv2_test_sliced$V8)
#visualizations follow, investigate shape of data

pv2_test_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov, FUN = median), y=m_cov, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##these look slihglty problematic
pv2_test_sliced_AB331651 <- filter(pv2_test_sliced, V6 == "papilloma.PPV_E1_C.Bos_taurus_papillomavirus_10:AB331651")
# write.table(pv2_test_sliced_AB331651$V1, "pv2_test_sliced_AB331651.txt", row.names = FALSE, col.names = FALSE)

pv2_test_sliced_MZ189241 <- filter(pv2_test_sliced, V6 == "papilloma.PPV_E1_DBD.Human_papillomavirus_type_18:MZ189241")
# write.table(pv2_test_sliced_MZ189241$V1, "pv2_test_sliced_MZ189241", row.names = FALSE, col.names = FALSE)
