library(tidyverse)
library(reshape2)
#analyze the HMMER output of logan2 run

hmmsearch_clean <- function(input_domtbl){
  input_domtbl_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:22]
  input_domtbl_scan <- input_domtbl_scan[!(is.na(input_domtbl_scan$V21) | input_domtbl_scan$V21=="" | !(input_domtbl_scan$V2 == "-")), ]
  
  colnames(input_domtbl_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                 "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}

pv2_fasta_star_sto <- hmmsearch_clean("pv2_fasta_star_sto.domtbl")
pv2_fasta_star_sto$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", pv2_fasta_star_sto$query_acc)



#check if there are any ORFs that are annotated as more than 1 thing
length(pv2_fasta_star_sto$query_acc) - length(unique(pv2_fasta_star_sto$query_acc))
#yes, im guessing it is E1/E2 domains (probably)
#should define those that have all 3 domains of E1 and both E2 domains as having those genes correctly


#can probably rewrite the below as functions
#to do for E2, essentially filter out contigs that have E2_N, match against list that have E2_C
pv2_fasta_star_sto_E2_N <- filter(pv2_fasta_star_sto, pfam == "PPV_E2_N")
pv2_fasta_star_sto_E2_C <- filter(pv2_fasta_star_sto, pfam == "PPV_E2_C")

pv2_fasta_star_sto_E2_full <- inner_join(pv2_fasta_star_sto_E2_N, pv2_fasta_star_sto_E2_C, by = c('query_acc'))
pv2_fasta_star_sto_E2_full <- pv2_fasta_star_sto_E2_full[,1:24]
pv2_fasta_star_sto_E2_full$pfam.x <- c("E2_full")
names(pv2_fasta_star_sto_E2_full) <- sub(".x", "", names(pv2_fasta_star_sto_E2_full))
hist(pv2_fasta_star_sto_E2_full$qlen)

pv2_fasta_star_sto_E2_full <- filter(pv2_fasta_star_sto_E2_full, qlen >= 250)
E2_pave <- read.table("E2_lengths.txt", sep = ",")
min(E2_pave$V2)

ggplot(pv2_fasta_star_sto_E2_full, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("E2_full length dist.")

E2_full_list <-select(pv2_fasta_star_sto_E2_full, pfam)

###prepare bedfile for the nt coordinate of E2s
pv2_fasta_star_sto_E2_full$query_acc_clean <-  sub("_REVERSE_SENSE_", "", pv2_fasta_star_sto_E2_full$query_acc_clean)
pv2_fasta_star_sto_E2_full$query_acc_clean <- sub("_$","",pv2_fasta_star_sto_E2_full$query_acc_clean)

pv2_fasta_star_sto_E2_full <- pv2_fasta_star_sto_E2_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_end"), "(.*)_([^_]+)$")
pv2_fasta_star_sto_E2_full <- pv2_fasta_star_sto_E2_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_start"), "(.*)_([^_]+)$")
pv2_fasta_star_sto_E2_full <- pv2_fasta_star_sto_E2_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf"), "(.*)_([^_]+)$")

pv2_fasta_star_sto_E2_coord <- select(pv2_fasta_star_sto_E2_full, query_acc_clean, orf_start, orf_end)
write.table(pv2_fasta_star_sto_E2_coord, "pv2_fasta_star_sto_E2_coord.bed", quote = F, col.names = F, row.names = F, sep = "\t")


#repeat for E1
E1_pave <- read.table("E1_lengths.txt", sep = ",")
min(E1_pave$V2)
ggplot(E1_pave, aes(x=V2)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("E1_full length dist.")

pv2_fasta_star_sto_E1_N <- filter(pv2_fasta_star_sto, pfam == "PPV_E1_N")
pv2_fasta_star_sto_E1_C <- filter(pv2_fasta_star_sto, pfam == "PPV_E1_C")
pv2_fasta_star_sto_E1_DBD <- filter(pv2_fasta_star_sto, pfam == "PPV_E1_DBD")

pv2_fasta_star_sto_E1_full <- inner_join(pv2_fasta_star_sto_E1_N, pv2_fasta_star_sto_E1_C, by = c('query_acc'))
pv2_fasta_star_sto_E1_full <- inner_join(pv2_fasta_star_sto_E1_full, pv2_fasta_star_sto_E1_DBD, by = c('query_acc'))
pv2_fasta_star_sto_E1_full <- pv2_fasta_star_sto_E1_full[,1:24]
names(pv2_fasta_star_sto_E1_full) <- sub(".x", "", names(pv2_fasta_star_sto_E1_full))
pv2_fasta_star_sto_E1_full$pfam <- c("E1_full")
hist(pv2_fasta_star_sto_E1_full$qlen)
pv2_fasta_star_sto_E1_full <- filter(pv2_fasta_star_sto_E1_full, qlen >= 500)


ggplot(pv2_fasta_star_sto_E1_full, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("E1_full length dist.")

E1_full_list <-select(pv2_fasta_star_sto_E1_full, pfam, contig)

#for L2
pv2_fasta_star_sto_L2 <- filter(pv2_fasta_star_sto, pfam == "Late_protein_L2")
hist(pv2_fasta_star_sto_L2$qlen)
ggplot(pv2_fasta_star_sto_L2, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L2 length dist.")

##trim for L2 less than 350 aa -> PAVE = 450-550 residue long, with shortest on refseq = 358  
L2_pave <- read.table("L2_lengths.txt", sep = ",")
sd(L2_pave$V2)
ggplot(L2_pave, aes(x=V2)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L2_pave length dist.")

#354 = smallest size -50 aa seems appropriate, maybe a bit lenient

pv2_fasta_star_sto_L2_full <- filter(pv2_fasta_star_sto_L2, qlen >= 350)


ggplot(pv2_fasta_star_sto_L2_full, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L2_full length dist.")

L2_full_list <-select(pv2_fasta_star_sto_L2_full, pfam, contig)


#for L1

L1_pave <- read.table("L1_lengths.txt", sep = ",")
min(L1_pave$V2)
ggplot(L1_pave, aes(x=V2)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L1_pave length dist.")

pv2_fasta_star_sto_L1 <- filter(pv2_fasta_star_sto, pfam == "Late_protein_L1")
hist(pv2_fasta_star_sto_L1$qlen)



ggplot(pv2_fasta_star_sto_L1, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L1 length dist.")

pv2_fasta_star_sto_L1_full <- filter(pv2_fasta_star_sto_L1, qlen >= 350)
ggplot(pv2_fasta_star_sto_L1_full, aes(x=qlen)) +
  geom_histogram(binwidth = 5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw() +
  ggtitle("L1_full length dist.")

L1_full_list <-select(pv2_fasta_star_sto_L1_full, pfam, contig)

###prepare bedfile for the nt coordinate of L1s
pv2_fasta_star_sto_L1_full$query_acc_clean <-  sub("_REVERSE_SENSE_", "", pv2_fasta_star_sto_L1_full$query_acc_clean)
pv2_fasta_star_sto_L1_full$query_acc_clean <- sub("_$","",pv2_fasta_star_sto_L1_full$query_acc_clean)

pv2_fasta_star_sto_L1_full <- pv2_fasta_star_sto_L1_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_end"), "(.*)_([^_]+)$")
pv2_fasta_star_sto_L1_full <- pv2_fasta_star_sto_L1_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_start"), "(.*)_([^_]+)$")
pv2_fasta_star_sto_L1_full <- pv2_fasta_star_sto_L1_full %>% extract(query_acc_clean, into = c("query_acc_clean", "orf"), "(.*)_([^_]+)$")

pv2_fasta_star_sto_L1_coord <- select(pv2_fasta_star_sto_L1_full, query_acc_clean, orf_start, orf_end)
write.table(pv2_fasta_star_sto_L1_coord, "pv2_fasta_star_sto_L1_coord.bed", quote = F, col.names = F, row.names = F, sep = "\t")
  ##decide here if i want to trim for full L1?
core_full_list <- rbind(L2_full_list, L1_full_list, E2_full_list, E1_full_list)
core_full_list <- add_column(core_full_list, present = 1)
core_full_list_tab <- core_full_list %>% 
  pivot_wider(names_from=pfam, values_from=present)

core_full_list_tab <- as.data.frame(core_full_list_tab)
core_full_list_tab[core_full_list_tab == "NULL"] <- 0
core_full_list_tab[core_full_list_tab == "c(1, 1)"] <- 2
core_full_list_tab[core_full_list_tab == "c(1, 1, 1)"] <- 3
core_full_list_tab[core_full_list_tab == "c(1, 1, 1, 1)"] <- 4
core_full_list_tab[core_full_list_tab == "c(1, 1, 1, 1, 1)"] <- 5
core_full_list_tab[core_full_list_tab == "c(1, 1, 1, 1, 1, 1)"] <- 6

core_full_list_tab[,2:5] <- sapply(core_full_list_tab[,2:5],as.numeric)

core_full_list_tab_pure <- core_full_list_tab %>%
  filter(if_all(Late_protein_L2:E1_full, ~ .x == "1"))
names(core_full_list_tab_pure) <- sub("contig", "contig_clean", names(core_full_list_tab_pure))


## get diamond file to put together metadata
concat_pv2_fil_pv_FP$contig_clean <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", concat_pv2_fil_pv_FP$contig)
concat_pv2_fil_pv_FP_md <- select(concat_pv2_fil_pv_FP, contig_clean, contig, length, strand, model_name, per_iden, e_val, library_source)

char_contig <- left_join(core_full_list_tab_pure, concat_pv2_fil_pv_FP_md, by = "contig_clean")

#curiousty, what about E7E6 with HMM -> more sensitive should be pretty good



