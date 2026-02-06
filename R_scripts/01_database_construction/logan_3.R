#logan3
#i love tidyverse
library(tidyverse)
library(ggplot2)
#function to parse the domtbl output
hmmsearch_clean <- function(input_domtbl) {
  input_domtbl_scan <- read.table(input_domtbl, sep = "", header = FALSE, fill = TRUE) %>% na.omit() %>% .[, 1:22]  
  colnames(input_domtbl_scan) <- c(
    "query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full",
    "score_full", "bias_full", "#", "of", "c_eval", "i_eval", "score_one",
    "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to",
    "envcoord_from", "envcoord_to", "acc"
  )
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}

#import from hmmer
feb_7_pv <- hmmsearch_clean("files/feb_7_pv_fil_form.aa.domtbl")
feb_7_pv$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", feb_7_pv$query_acc)
feb_7_pv$library <- sub("_.*", "\\1", feb_7_pv$contig)
length(unique(feb_7_pv$library))

#some graphs
#for contig length
a1 <- mean(feb_7_pv$qlen)
b1 <- median(feb_7_pv$qlen)

ggplot(feb_7_pv, aes(x = qlen)) +
  geom_histogram(binwidth = 0.05, alpha = 0.5, fill = "#156082") +
  geom_density(aes(y = 0.05 * after_stat(count)), alpha = 0.5) +
  scale_x_log10() +
  annotate("text", x = 1000, y = 150000,
           label = paste0("Mean = ", round(a1, digits = 3), "\n",
                          "Median = ", round(b1, digits = 3))) +
  theme_bw() +
  ggtitle("Contig Length Distribution") +
  xlab("Contig Length (aa, log scale)") +
  ylab("Count") +
  annotation_logticks(sides = "b")

#make "superfamily" for some of the categories
feb_7_pv_oo <- unique(feb_7_pv[c("query_acc", "pfam")])
 
feb_7_pv_count <- feb_7_pv_oo %>%
  count(pfam) %>%
  group_by(pfam)

feb_7_pv_count$super <- feb_7_pv_count$pfam %>% sub(".*E1.*.", "E1", .) %>% 
  sub(".*E2.*.", "E2", .)

#look at distribution of hits
feb_7_pv_count <- filter(feb_7_pv_count, super  %in% c("E1", "E2", "Late_protein_L1", "Late_protein_L2"))

ggplot(feb_7_pv_count, aes(fill = pfam, y = n, x = super)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Set3")
  


#filter for those hit against the Pfam L1
feb_7_pv_L1 <- filter(feb_7_pv, pfam == "Late_protein_L1")

length(unique(feb_7_pv_L1$library))

#307441
length(unique(feb_7_pv_L1$query_acc)) - length(unique(feb_7_pv$query_acc))
#1058413
#generate table of these names to be used to extract sequences from the ORF file
a2 <- mean(feb_7_pv_L1$qlen)
b2 <- median(feb_7_pv_L1$qlen)
ggplot(feb_7_pv_L1, aes(x = qlen)) +
  geom_histogram(binwidth = 0.05, alpha = 0.5, fill = "#D86ECC") +
  geom_density(aes(y = 0.05 * after_stat(count)), alpha = 0.5) +
  scale_x_log10() +
  annotate("text", x = 300, y = 40000,
           label = paste0("Mean = ", round(a2, digits = 3), "\n",
                          "Median = ", round(b2, digits = 3))) +
  theme_bw() +
  ggtitle("L1 Contig Length Distribution") +
  xlab("Contig Length (aa, log scale)") +
  ylab("Count") +
  annotation_logticks(sides = "b")


write.table(feb_7_pv_L1$query_acc, "outputs/feb_7_pv_L1_orf.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


#next, tabulate which are "full"
feb_7_pv_L1_BI <- hmmsearch_clean("files/feb_7_pv_fil_form_L1_BI.domtbl")

feb_7_pv_L1_BI$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", feb_7_pv_L1_BI$query_acc)
feb_7_pv_L1_BI$library <- sub("_.*", "\\1", feb_7_pv_L1_BI$contig)

#record region presence in contigs
#"low fruit" pass
feb_7_pv_L1_BI_mat <- feb_7_pv_L1_BI %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

#filter for those that have all regions hit :)
feb_7_pv_L1_BI_mat_all <- feb_7_pv_L1_BI_mat %>%
  filter(if_all(ncbi_and_novel_L1_B:ncbi_and_novel_L1_I, ~ .x > 0))

feb_7_pv_L1_BI_mat_all$library <- sub("_.*", "\\1", feb_7_pv_L1_BI_mat_all$query_acc)
write.table(select(feb_7_pv_L1_BI_mat_all, library), "outputs/feb_7_pv_L1_BI_mat_full_library.list", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(select(feb_7_pv_L1_BI_mat_all, query_acc), "outputs/feb_7_pv_L1_BI_mat_full.list", quote = FALSE, col.names = FALSE, row.names = FALSE)

#do a REALLY messy extraction of the coordinates
L1_BI_list <- feb_7_pv_L1_BI_mat_all$query_acc

feb_7_pv_L1_full <- filter(feb_7_pv_L1_BI, query_acc %in% L1_BI_list)
length(unique(feb_7_pv_L1_full$library))

L1_full_list <- select(feb_7_pv_L1_full, pfam, contig)
L1_full_list$pfam <- c("L1_full")
L1_full_list <- L1_full_list[!duplicated(L1_full_list), ]

#get the most "leftmost" coordinate
feb_7_pv_L1_envcoord_from <- feb_7_pv_L1_full %>% 
  group_by(query_acc) %>% 
  slice_min(n = 1, envcoord_from) %>% 
  select(query_acc, envcoord_from)

#get the most "rightmost" coordinate
feb_7_pv_L1_envcoord_to <- feb_7_pv_L1_full %>% 
  group_by(query_acc) %>% 
  slice_max(n = 1, envcoord_to) %>% 
  select(query_acc, envcoord_to)

feb_7_pv_L1_envcoords <- left_join(feb_7_pv_L1_envcoord_from, feb_7_pv_L1_envcoord_to, by = "query_acc")
#calculate nucleotide positions
feb_7_pv_L1_envcoords$nuc_from <- feb_7_pv_L1_envcoords$envcoord_from*2 + (feb_7_pv_L1_envcoords$envcoord_from - 3)
feb_7_pv_L1_envcoords$nuc_to <- feb_7_pv_L1_envcoords$envcoord_to*3 - 1

write.table(select(feb_7_pv_L1_envcoords, query_acc, nuc_from, nuc_to), "outputs/feb_7_pv_L1_envcoords.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

#look for those that have all hits somewhere in libraries -1
feb_7_pv_L1_BI_lib <- filter(feb_7_pv_L1_BI, !query_acc %in% L1_BI_list)
feb_7_pv_L1_BI_lib <- subset(feb_7_pv_L1_BI_lib, grepl("_L", query_acc))

feb_7_pv_L1_BI_lib_tab <- feb_7_pv_L1_BI_lib %>%
  count(library, pfam) %>%
  spread(pfam, n, fill = 0)

#convert to yes/no matrix
feb_7_pv_L1_BI_lib_logi <- feb_7_pv_L1_BI_lib_tab
feb_7_pv_L1_BI_lib_logi[-1] <- as.data.frame(+ sapply(feb_7_pv_L1_BI_lib_logi[-1], as.logical))
feb_7_pv_L1_BI_lib_logi$pa_sums <- rowSums(feb_7_pv_L1_BI_lib_logi[,2:7])

ggplot(feb_7_pv_L1_BI_lib_logi, aes(x =pa_sums)) +
  geom_histogram(binwidth = 0.5, alpha =0.5, fill ="#D86ECC") +
  theme_bw() +
  ggtitle("Distribution of model hits") +
  xlab("number of hits") +
  ylab("freq")

feb_7_pv_L1_BI_lib_logi_4 <- filter(feb_7_pv_L1_BI_lib_logi, pa_sums > 2)

#have to extract list of expected contigs from each library
feb_7_pv_L1_BI_partials <- filter(feb_7_pv_L1_BI, !query_acc %in% feb_7_pv_L1_full)


lib_sums <- feb_7_pv_L1_BI_lib_logi_4[c(1, 8)]

feb_7_pv_L1_BI_partials_sums <- left_join(feb_7_pv_L1_BI_partials, lib_sums, by = "library")

#running this produces a lot of files!
# write_file_func_2 <- function(lib){
#   library_filtered <- filter(feb_7_pv_L1_BI_partials_sums, library == lib)
#   write.table(unique(library_filtered$query_acc), paste("outputs/2024.02.19.subgraph/", lib, "_list.txt", sep = ""), quote = FALSE, col.names = FALSE, row.names = FALSE)
# }

#test func
#write_file_func("SRR25663671")

feb_7_pv_L1_BI_lib_logi_4_list <- feb_7_pv_L1_BI_lib_logi_4$library

for (i in feb_7_pv_L1_BI_lib_logi_4_list){
  write_file_func_2(i)
}

write.table(feb_7_pv_L1_BI_lib_logi_4_list, "outputs/feb_7_pv_L1_BI_lib_logi_4_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)


###for novelty search
feb_7_L1_nuc_hmmer_env <- read.table("files/feb_7_L1_nuc_hmmer_env.tsv", sep = "\t", header = FALSE)
feb_7_L1_nuc_hmmer_env$library <- gsub("Score.*", "", feb_7_L1_nuc_hmmer_env$V1)
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
length(unique(feb_7_L1_nuc_hmmer_env_less_90_nt))


write.table(feb_7_L1_nuc_hmmer_env_less_90_nt, "outputs/feb_7_L1_nuc_hmmer_env_less_90_nt_ugh.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

hist(feb_7_L1_nuc_hmmer_env_sliced_90$V9)

##some are not here!
length(unique(feb_7_L1_nuc_hmmer_env$V1))
# 859 of 942, so 83 are unaccounted for

#create list of the ones that were hit, in general
hit_list <- unique(feb_7_L1_nuc_hmmer_env_fil$V1)
write.table(hit_list, "outputs/hit_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

#look for low conf hits that were not represented in either novel already, or other filter set
feb_7_L1_nuc_hmmer_env_fil_low_conf_hits <- filter(feb_7_L1_nuc_hmmer_env_fil_low_conf, !V1 %in% feb_7_L1_nuc_hmmer_env_less_90_nt)
feb_7_L1_nuc_hmmer_env_fil_low_conf_hits <- filter(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits, !V1 %in% feb_7_L1_nuc_hmmer_env_fil$V1)
length(unique(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits$library))

feb_7_L1_nuc_hmmer_env_fil_low_conf_hits_list <- unique(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits$V1)

length(unique(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits$V1))

write.table(feb_7_L1_nuc_hmmer_env_fil_low_conf_hits_list, "outputs/feb_7_L1_nuc_hmmer_env_fil_low_conf_hits_list_ugh.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

##
sra_metadata <- read.table("files/lib_source_big.txt", sep = "\t", fill = TRUE, header = FALSE)
#random issue manually fix
sra_metadata["144301", "V2"] <- "GENOMIC"

feb_7_pv_lib_type <- left_join(feb_7_pv, sra_metadata, by = c("library" = "V1"))

#group genes
feb_7_pv_lib_type$pfam_grouped <-gsub("_N","",as.character(feb_7_pv_lib_type$pfam))
feb_7_pv_lib_type$pfam_grouped <-gsub("_C","",as.character(feb_7_pv_lib_type$pfam_grouped))
feb_7_pv_lib_type$pfam_grouped <-gsub("_DBD","",as.character(feb_7_pv_lib_type$pfam_grouped))

#assign unknown to unk
feb_7_pv_lib_type[is.na(feb_7_pv_lib_type)] <- "UNK"


#some issues with those libraries...just get rid of for now
feb_7_pv_lib_type <- filter(feb_7_pv_lib_type, pfam != "ncbi_and_novel_L1_JR")
feb_7_pv_lib_type <- filter(feb_7_pv_lib_type, library != "SRR10752062")
feb_7_pv_lib_type <- filter(feb_7_pv_lib_type, library != "SRR10752069")

unique(feb_7_pv_lib_type$V2)

feb_7_pv_lib_type$pfam <- factor(feb_7_pv_lib_type$pfam,
levels = c("Late_protein_L1", "Late_protein_L2",
  "E7", "E6", "Papilloma_E5",
  "Pap_E4", "PPV_E1_C", "PPV_E1_DBD",
  "PPV_E1_N", "PPV_E2_C", "PPV_E2_N"))
feb_7_pv_lib_type$V2 <- factor(feb_7_pv_lib_type$V2, levels = c("GENOMIC", "GENOMIC SINGLE CELL", "METAGENOMIC", "TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "VIRAL RNA", "METATRANSCRIPTOMIC", "SYNTHETIC", "UNK", "OTHER"))

p <- feb_7_pv_lib_type %>%
  dplyr::filter(pfam %in% c("Late_protein_L1", "Late_protein_L2", "E7", "E6", "Papilloma_E5", "Pap_E4", "PPV_E2_C", "PPV_E2_N", "PPV_E1_C", "PPV_E1_DBD", "PPV_E1_N")) %>%
  ggplot(aes(fill=pfam, x =V2)) + 
  geom_bar(position="fill") +
  scale_fill_manual(values = c("#78206E",  "#D86ECC",  "#04277C", "#074C8E", "#0B76A0", "#29A2AF", "#48BDB2", "#68CBAF","#89D7B3", "#AAE4BE", "#CCEFD2"))+
  theme_bw()  + theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust =1))

p



feb_7_pv_lib_type_RNA <- filter(feb_7_pv_lib_type, V2 %in% c("TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "VIRAL RNA", "METATRANSCRIPTOMIC"))

feb_7_pv_lib_type_DNA <- filter(feb_7_pv_lib_type, V2 %in% c("GENOMIC", "GENOMIC SINGLE CELL", "METAGENOMIC"))

###do verts visualization
verts_0306 <- read.table("files/verts_info_clean_final_0306.txt", fill = TRUE)
verts_0306$index <- seq.int(nrow(verts_0306))
verts_0306[-1] <- t(apply(verts_0306[-1], 1,
                     function(x) replace(x, duplicated(x), NA)))
verts_0306$info <- verts_0306 %>% 
  rowwise() %>%
  mutate(match = str_subset(c_across(), "Saving")[1]) %>% 
  ungroup

verts_0306 <- as.data.frame(verts_0306)
verts_0306_trues <- filter(verts_0306, grepl("gfas",V15))

verts_0306_trues <- select(verts_0306_trues, V1, V2, V15)
verts_0306_trues$V15 <- gsub("gfas/", "", verts_0306_trues$V15)
verts_0306_trues$V15 <- gsub("_.*", "", verts_0306_trues$V15)
verts_0306_trues <- unique(verts_0306_trues)
verts_0306_trues$V2 <- as.numeric(verts_0306_trues$V2)
verts_0306_trues$logV2 <- log10(verts_0306_trues$V2)

mean(verts_0306_trues$logV2)
median(verts_0306_trues$logV2)
log10(50000)

ggplot(verts_0306_trues, aes(x =V2)) +
  geom_histogram(binwidth = 0.1, fill ="#bebada", color ="#e9ecef", alpha =0.9) +
  geom_density(aes(y =0.1*after_stat(count))) +
  scale_x_log10() +
  annotation_logticks(sides="b") +
  theme_bw()

##do %identity histograms
feb_7_L1_nuc_hmmer_env_sliced_90$V9

ggplot(feb_7_L1_nuc_hmmer_env_sliced, aes(x =V9)) +
  geom_histogram(binwidth = 0.5, fill ="#bebada", color ="#e9ecef", alpha =0.9) +
  geom_density(aes(y =1*after_stat(count))) +
  theme_bw()

#zoom on those less than 90

ggplot(feb_7_L1_nuc_hmmer_env_sliced_90, aes(x =V9)) +
  geom_histogram(binwidth = 0.5, fill ="#6253AD", color ="#e9ecef", alpha =0.9) +
  geom_density(aes(y =1*after_stat(count))) +
  theme_bw()

##for all prs for fullness (non bullshit ones)
all_pr_fullness <- hmmsearch_clean("files/all_pr_need_hmmer.pep.domtbl")
all_pr_fullness <- unique(all_pr_fullness)

all_pr_fullness_mat <- all_pr_fullness %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

#filter for those that have all regions hit :)
all_pr_fullness_mat_all <- all_pr_fullness_mat %>%
  filter(if_all(ncbi_and_novel_L1_B:ncbi_and_novel_L1_I, ~ .x > 0))

all_pr_full_contigs <- as.data.frame(unique(all_pr_fullness_mat_all$query_acc))
write.table(all_pr_full_contigs, "outputs/all_pr_full_contigs.list", sep = "\t", col.names = F, row.names = F, quote = F)

feb_7_pv_L1_BI_mat_all$library <- sub("_.*", "\\1", feb_7_pv_L1_BI_mat_all$query_acc)

#from the pfam models, identify accessory genes that are on each contig
all_pr_addt_orf <- hmmsearch_clean("files/all_logan.domtbl")
all_pr_addt_orf_e_5 <- filter(all_pr_addt_orf, eval_full < 0.00001)
all_pr_addt_orf_e_5$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", all_pr_addt_orf_e_5$query_acc)

all_pr_addt_orf_e_5$pfam <- all_pr_addt_orf_e_5$pfam  %>% sub(".*E1.*.", "E1", .) %>% 
  sub(".*E2.*.", "E2", .)

all_pr_addt_orf_e_5_mat <- all_pr_addt_orf_e_5 %>%
  count(contig, pfam) %>%
  spread(pfam, n, fill = 0)

all_pr_addt_orf_e_5_mat_pa <- all_pr_addt_orf_e_5_mat %>%
  mutate(across(where(is.numeric), ~ifelse(.x >0, 1, 0)))

all_pr_addt_orf_e_5_mat_pa_acc <- select(all_pr_addt_orf_e_5_mat_pa, contig, E1, E2, E6, E7, Late_protein_L2, Pap_E4, Papilloma_E5)



#for the "multi" models, if one hit is confident then it is good enough as present
hist(filter(all_pr_addt_orf_e_5, pfam == "E6")$score_full, breaks = 50)

