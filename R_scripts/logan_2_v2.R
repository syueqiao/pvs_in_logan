################################################HOLY MOLY STOP-STOP#######################################################
#analyze the HMMER output of logan2
library(tidyverse)
library(ggplot2)
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
# colnames(pv2_fasta_sto_sto) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
#                                  "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc", "misc2", "query_acc_clean", "contig")
pv2_fasta_sto_sto_dist <- pv2_fasta_sto_sto
pv2_fasta_sto_sto_dist$pfam <- sub(".*E1.*.", "E1", pv2_fasta_sto_sto$pfam)
pv2_fasta_sto_sto_dist$pfam <- sub(".*E2.*.", "E2", pv2_fasta_sto_sto_dist$pfam)
pv2_fasta_sto_sto_dist <- select(pv2_fasta_sto_sto_dist, contig, pfam)
pv2_fasta_sto_sto_dist <- unique(pv2_fasta_sto_sto_dist)

pv2_fasta_sto_sto_dist_tab <- as.data.frame(table(pv2_fasta_sto_sto_dist$pfam))
rows <- c(1, 2, 5, 6)
pv2_fasta_sto_sto_dist_tab <- pv2_fasta_sto_sto_dist_tab[rows,]

colorsbox <- c(E1 = "#E97132", E2 = "#4EA72E", Late_protein_L2 ="#A02B93", Late_protein_L1 = "#104862")
ggplot(pv2_fasta_sto_sto_dist_tab, aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  scale_fill_manual(values=colorsbox) +
  ggtitle("HMM Hit Distribution") +
  xlab("HMM") +
  ylab("Count") +
  theme(legend.position="none")

#check if there are any ORFs that are annotated as more than 1 thing
# length(pv2_fasta_sto_sto$query_acc) - length(unique(pv2_fasta_sto_sto$query_acc))
#yes, im guessing it is E1/E2 domains (probably)
#should define those that have all 3 domains of E1 and both E2 domains as having those genes correctly


#can probably rewrite the below as functions
#to do for E2, essentially filter out contigs that have E2_N, match against list that have E2_C
pv2_fasta_sto_sto_E2_N <- filter(pv2_fasta_sto_sto, pfam == "PPV_E2_N")
pv2_fasta_sto_sto_E2_C <- filter(pv2_fasta_sto_sto, pfam == "PPV_E2_C")

pv2_fasta_sto_sto_E2_full <- inner_join(pv2_fasta_sto_sto_E2_N, pv2_fasta_sto_sto_E2_C, by = c('query_acc'))
pv2_fasta_sto_sto_E2_full <- pv2_fasta_sto_sto_E2_full[,1:25]
pv2_fasta_sto_sto_E2_full$pfam.x <- c("E2_full")
names(pv2_fasta_sto_sto_E2_full) <- sub(".x", "", names(pv2_fasta_sto_sto_E2_full))
hist(pv2_fasta_sto_sto_E2_full$qlen)

# pv2_fasta_sto_sto_E2_full <- filter(pv2_fasta_sto_sto_E2_full, qlen >= 250)
#remove size for now, just use N and C?
# E2_pave <- read.table("E2_lengths.txt", sep = ",")
# min(E2_pave$V2)

ggplot(pv2_fasta_sto_sto_E2_full, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("E2_full length dist.")

E2_full_list <- select(pv2_fasta_sto_sto_E2_full, pfam, contig)


#repeat for E1
# E1_pave <- read.table("E1_lengths.txt", sep = ",")
# min(E1_pave$V2)
# ggplot(E1_pave, aes(x=V2)) +
#   geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
#   geom_density(aes(y=..count..)) +
#   theme_bw() +
#   ggtitle("E1_full length dist.")

pv2_fasta_sto_sto_E1_N <- filter(pv2_fasta_sto_sto, pfam == "PPV_E1_N")
pv2_fasta_sto_sto_E1_C <- filter(pv2_fasta_sto_sto, pfam == "PPV_E1_C")
pv2_fasta_sto_sto_E1_DBD <- filter(pv2_fasta_sto_sto, pfam == "PPV_E1_DBD")

pv2_fasta_sto_sto_E1_full <- inner_join(pv2_fasta_sto_sto_E1_N, pv2_fasta_sto_sto_E1_C, by = c('query_acc'))
pv2_fasta_sto_sto_E1_full <- inner_join(pv2_fasta_sto_sto_E1_full, pv2_fasta_sto_sto_E1_DBD, by = c('query_acc'))
pv2_fasta_sto_sto_E1_full <- pv2_fasta_sto_sto_E1_full[,1:25]
names(pv2_fasta_sto_sto_E1_full) <- sub(".x", "", names(pv2_fasta_sto_sto_E1_full))
pv2_fasta_sto_sto_E1_full$pfam <- c("E1_full")
hist(pv2_fasta_sto_sto_E1_full$qlen)
# pv2_fasta_sto_sto_E1_full <- filter(pv2_fasta_sto_sto_E1_full, qlen >= 500)


ggplot(pv2_fasta_sto_sto_E1_full, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("E1_full length dist.")

E1_full_list <- select(pv2_fasta_sto_sto_E1_full, pfam, contig)

#for L2, since there are no domains, just use P/A of this hit
pv2_fasta_sto_sto_L2 <- filter(pv2_fasta_sto_sto, pfam == "Late_protein_L2")
hist(pv2_fasta_sto_sto_L2$qlen)
ggplot(pv2_fasta_sto_sto_L2, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L2 length dist.")

##trim for L2 less than 350 aa -> PAVE = 450-550 residue long, with shortest on refseq = 358  
# L2_pave <- read.table("L2_lengths.txt", sep = ",")
# sd(L2_pave$V2)
# ggplot(L2_pave, aes(x=V2)) +
#   geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
#   geom_density(aes(y=..count..)) +
#   theme_bw() +
#   ggtitle("L2_pave length dist.")

#354 = smallest size -50 aa seems appropriate, maybe a bit lenient

#don't trim for length for now, just p/a
pv2_fasta_sto_sto_L2_full <- filter(pv2_fasta_sto_sto_L2, qlen >= 350)
pv2_fasta_sto_sto_L2_full <- pv2_fasta_sto_sto_L2



ggplot(pv2_fasta_sto_sto_L2_full, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L2_full length dist.")

L2_full_list <-select(pv2_fasta_sto_sto_L2_full, pfam, contig)


#for L1
# 
L1_pave <- read.table("L1_lengths.txt", sep = ",")
min(L1_pave$V2)
ggplot(L1_pave, aes(x=V2)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L1_pave length dist.")

#rough first view, before JR stuff
pv2_fasta_sto_sto_L1 <- filter(pv2_fasta_sto_sto, pfam == "Late_protein_L1")
write.table(pv2_fasta_sto_sto_L1$query_acc, "pv_all_fasta_sto_sto_L1_orf.txt", quote = F, col.names = F, row.names = F, sep = "\t")
pv2_fasta_sto_sto_L1_temp <- pv2_fasta_sto_sto_L1
pv2_fasta_sto_sto_L1_temp$library <- sub("_.*", "\\1", pv2_fasta_sto_sto_L1_temp$contig)
length(unique(pv2_fasta_sto_sto_L1_temp$library))
# hist(pv2_fasta_sto_sto_L1$qlen)
# 
# 
# 
# ggplot(pv2_fasta_sto_sto_L1, aes(x=qlen)) +
#   geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
#   geom_density(aes(y=..count..)) +
#   theme_bw() +
#   ggtitle("L1 length dist.")
# 
# pv2_fasta_sto_sto_L1_full <- pv2_fasta_sto_sto_L1
# ggplot(pv2_fasta_sto_sto_L1_full, aes(x=qlen)) +
#   geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
#   geom_density(aes(y=..count..)) +
#   theme_bw() +
#   ggtitle("L1_full length dist.")
# 
# 
# 
# L1_all_list <-select(pv2_fasta_sto_sto_L1_full, pfam, contig)



######################characterize L1 as full or partial
pv2_L1_sto_sto <- hmmsearch_clean("pv_all_L1.domtbl")
pv2_L1_sto_sto$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", pv2_L1_sto_sto$query_acc)
# colnames(pv2_L1_sto_sto) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
#                                  "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc", "misc2", "query_acc_clean", "contig")
# 
length(unique(pv2_L1_sto_sto$contig))
pv2_L1_sto_sto$library <- sub("_.*", "", pv2_L1_sto_sto$contig)
pv2_L1_sto_sto_uniq <- unique(pv2_L1_sto_sto)
pv2_L1_sto_sto_uniq$library  <- sub("_.*", "", pv2_L1_sto_sto_uniq$contig)
pv2_L1_sto_sto$pfam <- factor(x = pv2_L1_sto_sto$pfam,
                              levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B"))


pv2_L1_sto_sto_tab <- pv2_L1_sto_sto %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

#all L1 including partial = pv2_L1_sto_sto_tab
pv2_L1_sto_sto_tab_full <- pv2_L1_sto_sto_tab %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_B, ~ .x > 0))

L1_BI_list <- pv2_L1_sto_sto_tab_full$query_acc

pv2_L1_sto_sto_full_L1 <- filter(pv2_L1_sto_sto, query_acc_clean %in% L1_BI_list)

L1_full_list <- select(pv2_L1_sto_sto_full_L1, pfam, contig)
L1_full_list$pfam <- c("L1_full")
L1_full_list <- L1_full_list[!duplicated(L1_full_list), ]

pv2_L1_sto_sto_full_L1_envcoord_from <- pv2_L1_sto_sto_full_L1 %>% group_by(query_acc) %>% slice_min(n = 1, envcoord_from)
pv2_L1_sto_sto_full_L1_envcoord_from <- select(pv2_L1_sto_sto_full_L1_envcoord_from, query_acc, envcoord_from)

pv2_L1_sto_sto_full_L1_envcoord_to <- pv2_L1_sto_sto_full_L1 %>% group_by(query_acc) %>% slice_max(n = 1, envcoord_to)
pv2_L1_sto_sto_full_L1_envcoord_to <- select(pv2_L1_sto_sto_full_L1_envcoord_to, query_acc, envcoord_to)

pv2_L1_sto_sto_full_L1_envcoords <- left_join(pv2_L1_sto_sto_full_L1_envcoord_from, pv2_L1_sto_sto_full_L1_envcoord_to, by = "query_acc")
pv2_L1_sto_sto_full_L1_envcoords$nuc_from <- pv2_L1_sto_sto_full_L1_envcoords$envcoord_from*2 + (pv2_L1_sto_sto_full_L1_envcoords$envcoord_from - 3)
pv2_L1_sto_sto_full_L1_envcoords$nuc_to <- pv2_L1_sto_sto_full_L1_envcoords$envcoord_to*3 - 1


pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- sub("_REVERSE_SENSE_", "(-)", pv2_L1_sto_sto_full_L1_envcoords$query_acc)

pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- sub("_$","(+)",pv2_L1_sto_sto_full_L1_envcoords$query_nuc)
pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- sub("_([^_]*)$","-\\1",pv2_L1_sto_sto_full_L1_envcoords$query_nuc)
pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- sub("_([^_]*)$",":\\1",pv2_L1_sto_sto_full_L1_envcoords$query_nuc)
pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- gsub(":[^:]+:", ":", pv2_L1_sto_sto_full_L1_envcoords$query_nuc)
pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- sub("_([^_]*)$","-\\1",pv2_L1_sto_sto_full_L1_envcoords$query_nuc)


pv2_L1_sto_sto_full_L1_envcoords <- separate(pv2_L1_sto_sto_full_L1_envcoords, col=query_nuc, into=c('query_nuc', 'coords'), sep="(?<=:)")
pv2_L1_sto_sto_full_L1_envcoords <- separate(pv2_L1_sto_sto_full_L1_envcoords, col=coords, into=c('coords', 'sense'), sep="(?=\\()")
pv2_L1_sto_sto_full_L1_envcoords <- separate(pv2_L1_sto_sto_full_L1_envcoords, col=coords, into=c('coord_start', 'coord_end'), sep="-")
pv2_L1_sto_sto_full_L1_envcoords <- separate(pv2_L1_sto_sto_full_L1_envcoords, col=query_nuc, into=c('query_nuc', 'misc'), sep="-")
pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- sub("$","_:",pv2_L1_sto_sto_full_L1_envcoords$query_nuc)


pv2_L1_sto_sto_full_L1_envcoords[,8:9] <- sapply(pv2_L1_sto_sto_full_L1_envcoords[,8:9],as.numeric)


low  <- pmin(pv2_L1_sto_sto_full_L1_envcoords$coord_start, pv2_L1_sto_sto_full_L1_envcoords$coord_end)
high <- pmax(pv2_L1_sto_sto_full_L1_envcoords$coord_start, pv2_L1_sto_sto_full_L1_envcoords$coord_end)
pv2_L1_sto_sto_full_L1_envcoords$coord_start <- low
pv2_L1_sto_sto_full_L1_envcoords$coord_end  <- high
pv2_L1_sto_sto_full_L1_envcoords$coord_start <- pv2_L1_sto_sto_full_L1_envcoords$coord_start-1


pv2_L1_sto_sto_full_L1_envcoords$coord_start <- sub("$","-",pv2_L1_sto_sto_full_L1_envcoords$coord_start)
pv2_L1_sto_sto_full_L1_envcoords$final_names <- paste(pv2_L1_sto_sto_full_L1_envcoords$query_nuc, pv2_L1_sto_sto_full_L1_envcoords$coord_start, pv2_L1_sto_sto_full_L1_envcoords$coord_end, pv2_L1_sto_sto_full_L1_envcoords$sense, sep = "")
pv2_L1_sto_sto_full_L1_envcoords$query_nuc <- sub(":","_",pv2_L1_sto_sto_full_L1_envcoords$query_nuc)

pv2_L1_sto_sto_nuc_coordinates <- select(ungroup(pv2_L1_sto_sto_full_L1_envcoords), final_names, nuc_from, nuc_to)
#>DRR013318_1_ka_f_37.840_:3958-5536(-)
write.table(pv2_L1_sto_sto_nuc_coordinates, "pv_all_L1_sto_sto_nuc_coordinates.txt", quote = F, col.names = F, row.names = F, sep = "\t")

###prepare bedfile for the nt coordinate of L1s
pv2_L1_sto_sto_full_L1$query_acc_clean <-  sub("_REVERSE_SENSE_", "", pv2_L1_sto_sto_full_L1$query_acc_clean)
pv2_L1_sto_sto_full_L1$query_acc_clean <- sub("_$","",pv2_L1_sto_sto_full_L1$query_acc_clean)

pv2_L1_sto_sto_full_L1 <- pv2_L1_sto_sto_full_L1 %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_end"), "(.*)_([^_]+)$")
pv2_L1_sto_sto_full_L1 <- pv2_L1_sto_sto_full_L1 %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_start"), "(.*)_([^_]+)$")
pv2_L1_sto_sto_full_L1 <- pv2_L1_sto_sto_full_L1 %>% extract(query_acc_clean, into = c("query_acc_clean", "orf"), "(.*)_([^_]+)$")

pv2_fasta_sto_sto_L1_coord <- select(pv2_L1_sto_sto_full_L1, query_acc_clean, orf_start, orf_end)
pv2_fasta_sto_sto_L1_coord$query_acc_clean <- sub("$","_",pv2_fasta_sto_sto_L1_coord$query_acc_clean)
write.table(pv2_fasta_sto_sto_L1_coord, "pv_all_fasta_sto_sto_L1_coord.txt", quote = F, col.names = F, row.names = F, sep = "\t")

##decide here if i want to trim for full L1?
core_full_list <- rbind(L2_full_list, L1_full_list, E2_full_list, E1_full_list)
core_full_list <- add_column(core_full_list, present = 1)
core_full_list_tab <- core_full_list %>% 
  pivot_wider(names_from=pfam, values_from=present, values_fn = function(x) paste(sum(x)))
core_full_list_tab[,2:5] <- sapply(core_full_list_tab[,2:5],as.numeric)
core_full_list_tab[is.na(core_full_list_tab)] <- 0
core_full_list_tab <- arrange(core_full_list_tab, L1_full, Late_protein_L2, E2_full, E1_full)

core_full_list_tab_count <- core_full_list_tab[,2:5]
core_full_list_tab_count[core_full_list_tab_count != 0] <- 1
core_full_list_tab_count <- core_full_list_tab_count %>% group_by(across()) %>% count %>% arrange(desc(n))
#visualize this?
out <- select(core_full_list, pfam, contig) %>% with(., split(contig, pfam))
library(ggVennDiagram)
ggVennDiagram(out, label_alpha = 0, category.names = c("Full E1","Full E2","Full L1", "Full L2")) + scale_fill_gradient(low="grey90", high = "#156082") +
  guides(fill = guide_legend(title = "Hit Count"))



#things that have 1 in all categories (not dealing with repeats)
core_full_list_tab_pure <- core_full_list_tab %>%
  filter(if_all(Late_protein_L2:E1_full, ~ .x == "1"))

#things that have something in all categories (probably has repeats)
core_full_list_tab_min <- core_full_list_tab %>%
  filter(if_all(Late_protein_L2:E1_full, ~ .x > 0))
names(core_full_list_tab_pure) <- sub("contig", "contig_clean", names(core_full_list_tab_pure))


## get diamond file to put together metadata
pv_all_pro <- read.table("pv", sep = "\t", fill = T)

concat_pv2_fil_pv_FP$contig_clean <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", concat_pv2_fil_pv_FP$contig)
concat_pv2_fil_pv_FP_md <- select(concat_pv2_fil_pv_FP, contig_clean, contig, length, strand, model_name, per_iden, e_val, library_source)

char_contig <- left_join(core_full_list_tab_pure, concat_pv2_fil_pv_FP_md, by = "contig_clean")

#curiousty, what about E7E6 with HMM -> more sensitive should be pretty good



###prepare bedfile for the nt coordinate of E2s
pv2_fasta_sto_sto_E2_full$query_acc_clean <-  sub("_REVERSE_SENSE_", "", pv2_fasta_sto_sto_E2_full$query_acc_clean)
pv2_fasta_sto_sto_E2_full$query_acc_clean <- sub("_$","",pv2_fasta_sto_sto_E2_full$query_acc_clean)

pv2_fasta_sto_sto_E2_full <- pv2_fasta_sto_sto_E2_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_end"), "(.*)_([^_]+)$")
pv2_fasta_sto_sto_E2_full <- pv2_fasta_sto_sto_E2_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_start"), "(.*)_([^_]+)$")
pv2_fasta_sto_sto_E2_full <- pv2_fasta_sto_sto_E2_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf"), "(.*)_([^_]+)$")

pv2_fasta_sto_sto_E2_coord <- select(pv2_fasta_sto_sto_E2_full, query_acc_clean, orf_start, orf_end)
write.table(pv2_fasta_sto_sto_E2_coord, "pv2_fasta_sto_sto_E2_coord.bed", quote = F, col.names = F, row.names = F, sep = "\t")

#read blastnt in
L1_nt <- read.table("nt_output_centroids.tsv", sep = "\t", header = F)

L1_nt_sliced <- L1_nt %>% group_by(V1) %>% slice_max(n = 1, V9)
L1_nt_sliced_90 <- filter(L1_nt_sliced, V9 < 90)
L1_nt_sliced_90 = L1_nt_sliced_90[!duplicated(L1_nt_sliced_90$V1),]
L1_nt_sliced_90$qcov <- abs(L1_nt_sliced_90$V2 - L1_nt_sliced_90$V3)/L1_nt_sliced_90$V4


L1_nt_qcov <- L1_nt_sliced
L1_nt_qcov$qcov <- abs(L1_nt_qcov$V2 - L1_nt_qcov$V3)/L1_nt_qcov$V4
L1_nt_qcov <- L1_nt_qcov %>% group_by(V1) %>% slice_max(n = 1, qcov)
L1_nt_qcov <- L1_nt_qcov %>% group_by(V1) %>% slice_min(n = 1, V10)
L1_nt_qcov_less90 <- filter(L1_nt_qcov, V9 < 90)

L1_less_90_nt <- as.data.frame(unique(L1_nt_sliced_90$V1))
L1_nt_full_out <- as.data.frame(unique(L1_nt$V1))
colnames(L1_less_90_nt) <- c("V1")

write.table(L1_less_90_nt, "L1_less_90_nt.txt", quote = F, col.names = F, row.names = F)
write.table(L1_nt_full_out, "L1_nt_full_out", quote = F, col.names = F, row.names = F)


#compare for input that was not hit in blastn
L1_headers_in <- read.table("L1_centroids_headers.txt", sep = "\t", header = F)
L1_headers_in <- as.data.frame(gsub(">","", L1_headers_in$V1))
L1_headers_in <- as.data.frame(gsub("_$","", L1_headers_in$`gsub(">", "", L1_headers_in$V1)`))
colnames(L1_headers_in) <- c("V1")

L1_nt_in <- as.data.frame(gsub(":.*","", L1_nt$V1))
L1_nt_in <- as.data.frame(gsub("_$","", L1_nt_in$`gsub(":.*", "", L1_nt$V1)`))
colnames(L1_nt_in) <- c("V1")

not_in_nt <- setdiff(L1_headers_in, L1_less_90_nt)
not_in_nt <- unique(not_in_nt)
pv2_fasta_sto_sto_L1_coord$query_acc_clean

length(unique(L1_headers_in$V1))

novel_L1 <- rbind(not_in_nt, L1_less_90_nt)
write.table(novel_L1, "novel_L1.txt", quote = F, col.names = F, row.names = F)

########################investigate how fragmented the contigs are
pv2_L1_sto_sto_f <- hmmsearch_clean("pv_all_L1.domtbl")

pv2_L1_sto_sto_frag <- pv2_L1_sto_sto_f
pv2_L1_sto_sto_frag$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", pv2_L1_sto_sto_frag$query_acc)
pv2_L1_sto_sto_frag$library <- sub("_.*", "\\1", pv2_L1_sto_sto_frag$contig)

#record those with all regions intact
pv2_L1_sto_sto_tab_frag_mat <- pv2_L1_sto_sto_frag %>%
  count(library, pfam) %>%
  spread(pfam, n, fill = 0)

pv2_L1_sto_sto_frag_all <- pv2_L1_sto_sto_tab_frag %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_B, ~ .x > 0))

#what about just B and I?
pv2_L1_sto_sto_tab_frag_BI <- pv2_L1_sto_sto_frag
pv2_L1_sto_sto_tab_frag_BI$pfam <- factor(x = pv2_L1_sto_sto_tab_frag_BI$pfam,
                              levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B"))

pv2_L1_sto_sto_tab_frag_BI <- pv2_L1_sto_sto_tab_frag_BI %>%
  count(library, pfam) %>%
  spread(pfam, n, fill = 0)

pv2_L1_sto_sto_tab_frag_BI <- pv2_L1_sto_sto_tab_frag_BI %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_B, ~ .x > 0))

#compare against full list from before
library_with_all <- pv2_L1_sto_sto_frag_all$library

length(unique(library_with_all))

pv2_L1_sto_sto_tab_full_lib <- sub("_.*", "\\1", pv2_L1_sto_sto_tab_full$query_acc)
length(unique(pv2_L1_sto_sto_tab_full_lib))


not_full_library <- setdiff(library_with_all, pv2_L1_sto_sto_tab_full_lib)

#so there are 5547 libraries with B and I hits on same contig
#9935 libraries with all models hit somewhere

##remove the contigs that are full
pv2_L1_sto_sto_frag_remove <- filter(pv2_L1_sto_sto_frag, !contig %in% L1_full_list$contig)

pv2_L1_sto_sto_frag_remove$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", pv2_L1_sto_sto_frag_remove$query_acc)
pv2_L1_sto_sto_frag_remove$library <- sub("_.*", "\\1", pv2_L1_sto_sto_frag_remove$contig)

#record those with all regions intact
pv2_L1_sto_sto_frag_remove_mat <- pv2_L1_sto_sto_frag_remove %>%
  count(library, pfam) %>%
  spread(pfam, n, fill = 0)

pv2_L1_sto_sto_frag_remove_all <- pv2_L1_sto_sto_frag_remove_mat %>%
  filter(if_all(ncbi_and_pave_L1_I:ncbi_and_pave_L1_B, ~ .x > 0))

#9296 libraries with all fragments still

# #bug squashing
# L1_pave_coords <- read.table("L1_pave_nuc_headers.txt", sep = "\t")
# L1_pave_coords <- separate(L1_pave_coords, col=V1, into=c('V1', 'V2'), sep=" ")
# L1_pave_coords <- separate(L1_pave_coords, col=V1, into=c('V1', 'V3'), sep="/")
# L1_pave_coords <- separate(L1_pave_coords, col=V3, into=c('V3', 'V4'), sep="-")
# 
# L1_pave_coords[,2:3] <- sapply(L1_pave_coords[,2:3],as.numeric)
# L1_pave_nuc <- L1_pave_coords
# 
# L1_pave_nuc$nuc_from <- L1_pave_nuc$V3*2 + (L1_pave_nuc$V3 - 3)
# L1_pave_nuc$nuc_to <- L1_pave_nuc$V4*3 - 1
# L1_pave_nuc <- select(L1_pave_nuc, V1, nuc_from, nuc_to)
# write.table(L1_pave_nuc, "L1_pave_nuc.bed", quote = F, col.names = F, row.names = F, sep = "\t")
# 
# L1_pave_pro <- L1_pave_coords
# L1_pave_pro <- select(L1_pave_pro, V1, V3, V4)
# write.table(L1_pave_pro, "L1_pave_pro.bed", quote = F, col.names = F, row.names = F, sep = "\t")

##making figures for msa

# 
# 
# ###invesetigate missing libraries
# pv2_libraries <- as.data.frame(gsub("(.+?)(\\_.*)", "\\1", concat_pv2_fil_pv_FP$V1))
# 
# pv2_libraries <- unique(pv2_libraries)
# colnames(pv2_libraries) <- c("V14")
# 
# large_clean_fil_libraries <- as.data.frame(unique(large_clean_fil$V14))
# colnames(large_clean_fil_libraries) <- c("V14")
# 
# in_pv2_not_large <- pv2_libraries %>% 
#   filter(!pv2_libraries$V14 %in% large_cl$ean_fil_libraries$V14) #For df1 values not in df2
# 
# in_large_not_pv2 <- large_clean_fil_libraries %>% 
#   filter(!large_clean_fil_libraries$V14 %in% pv2_libraries$V14) #For df2 values not in df1
# 
# #SRR10902309 woah! white rhino papillomavirus

##for new pvdb iteration
pv2_fasta_sto_sto_ieval_fil <- filter(pv2_fasta_sto_sto, eval_full < 0.000001)

pv2_fasta_sto_sto_ieval_fil$lib <- gsub("\\_.*","",pv2_fasta_sto_sto_ieval_fil$query_acc)

pvdb2_sra <- read.table("sra_pvdb2_taxid.txt", sep = "\t")
pv2_fasta_sto_sto_ieval_fil <- left_join(pv2_fasta_sto_sto_ieval_fil, pvdb2_sra, by = c("lib" = "V1"))

accessory <- c("Papilloma_E5", "E6", "E7", "Pap_E4")

#filter on full evalue, since we are looking at only a single orf this time. Doesn't really matter
#I think I should also filter on 
pv2_fasta_sto_sto_ieval_fil <- filter(pv2_fasta_sto_sto, eval_full < 0.00001)
pv2_fasta_sto_sto_ieval_fil <- filter(pv2_fasta_sto_sto_ieval_fil, !pfam %in% accessory)

pv2_fasta_sto_sto_ieval_fil <- pv2_fasta_sto_sto_ieval_fil %>%
  mutate(V2=replace_na(V2, "unk"))

pv2_fasta_sto_sto_ieval_fil <- pv2_fasta_sto_sto_ieval_fil %>% 
  group_by(query_acc) %>% 
  slice_max(score_one, n = 1) %>%
  ungroup()
pv2_fasta_sto_sto_ieval_fil_key <- select(pv2_fasta_sto_sto_ieval_fil, query_acc, pfam, V2)
pv2_fasta_sto_sto_ieval_fil_key$info <- paste0(pv2_fasta_sto_sto_ieval_fil_key$pfam, ".",  pv2_fasta_sto_sto_ieval_fil_key$V2)
pv2_fasta_sto_sto_ieval_fil_key_final <- select(pv2_fasta_sto_sto_ieval_fil_key, query_acc, info)
write.table(pv2_fasta_sto_sto_ieval_fil_key_final, sep = " ", "hmmer_annot_pv2.txt", row.names = F, col.names = F, quote = F)
