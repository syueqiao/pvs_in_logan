#logan parsing
library(tidyverse)

pap_pro_100 <- read.table("pap_pro_head_100k", sep = '\t', fill = T)
pap_pro_100_clean <- pap_pro_100[grepl("(\\+|-)", pap_pro_100$V5), ]
pap_pro_100_clean <- pap_pro_100_clean[grepl("*M$", pap_pro_100_clean$V12), ]

pap_pro_100_clean$V13 <- as.numeric(pap_pro_100_clean$V11)
pap_pro_100_clean$V4 <- as.numeric(pap_pro_100_clean$V4)

pap_pro_100_fil <- filter(pap_pro_100_clean, V13 < 0.0000000001)
pap_pro_100_fil$V14 <- sub("_.*", "", pap_pro_100_fil$V1)
pap_pro_100_fil <- pap_pro_100_fil[!grepl("LQ465342", pap_pro_100_fil$V6),]
pap_pro_100_fil <- pap_pro_100_fil[!grepl("OP721220", pap_pro_100_fil$V6),]
pap_pro_100_fil <- pap_pro_100_fil[!grepl("DRR027743", pap_pro_100_fil$V14),]
pap_pro_100_fil <- pap_pro_100_fil[!grepl("DRR027738", pap_pro_100_fil$V14),]
pap_pro_100_fil <- pap_pro_100_fil[!grepl("DRR027740", pap_pro_100_fil$V14),]




pap_pro_100_fil <- pap_pro_100_fil[!duplicated(pap_pro_100_fil[c(1,6)]),]
pap_pro_100_fil_top <- pap_pro_100_fil %>% group_by(V1) %>% top_n(1, V13)


pap_pro_100_fil_list <- split(pap_pro_100_fil, pap_pro_100_fil$V14)
# pap_pro_100_list_fil <- purrr::keep(pap_pro_100_fil_list, ~nrow(.) >= 5)

DRR091124
sliced <- pap_pro_100_fil %>% group_by(V1) %>% slice_max(order_by = V11, n = 1)

pap_pro_core <- read.table("../pap.core.pro", sep = '\t', fill = T)
pap_pro_core_fil <- filter(pap_pro_core, pap_pro_core$V11 < 0.0000000001)


sliced_core <- pap_pro_core_fil %>% group_by(V1) %>% slice_max(order_by = V11, n = 1)
sliced_core <- sliced_core[grepl("(\\+|-)", sliced_core$V5), ]
sliced_core <- sliced_core[grepl("*M$", sliced_core$V12), ]\
sliced_core$V4 <- as.numeric(sliced_core$V4)
sliced_core$V10 <- as.numeric(sliced_core$V10)


sliced_core_2 <- dplyr::filter(sliced_core, V4 > 400)
sliced_core_2$V11 <- as.numeric(sliced_core_2$V11)

sliced_core_3 <- dplyr::filter(sliced_core_2, V11 < 0.0000000001)

sliced_core_3 <- filter(sliced_core_3, V6 != "papilloma.Late_protein_L1.Human_papillomavirus:OP721220")
sliced_core_3 <- filter(sliced_core_3, V6 != "papilloma.Late_protein_L2.Human_papillomavirus:LQ465342")

table <- as.data.frame(table(sliced_core_3$V6))

tuni <- read.table("tuni_table.txt", sep = "\t", fill = T)
tuni_table <- as.data.frame(table(tuni$V6))

tuni_clean <- tuni[grepl("(\\+|-)", tuni$V5), ]
tuni_clean <- tuni_clean[grepl("*M$", tuni_clean$V12), ]

tuni_clean$V13 <- as.numeric(tuni_clean$V11)
tuni_clean$V4 <- as.numeric(tuni_clean$V4)

tuni_fil <- filter(tuni_clean, V13 < 0.00000001)
tuni_fil$V14 <- sub("_.*", "", tuni_fil$V1)
tuni_fil$V10 <- as.numeric(tuni_fil$V10)
tuni_fil$V7 <- as.numeric(tuni_fil$V7)
tuni_fil$V8 <- as.numeric(tuni_fil$V8)
tuni_fil$V9 <- as.numeric(tuni_fil$V9)



tuni_fil <- tuni_fil[!grepl("LQ465342", tuni_fil$V6),]
tuni_fil <- tuni_fil[!grepl("OP721220", tuni_fil$V6),]
tuni_fil <- tuni_fil[!grepl("AX020928|AX020938|OP721220|OP721198|LQ465342|JA802371|AX020936|AX020945|HE984563|MP134365|AF001600|AX020925", tuni_fil$V6),]
tuni_fil$cov <- abs(tuni_fil$V7-tuni_fil$V8)
tuni_fil$per_cov <- tuni_fil$cov/tuni_fil$V9

beet_nema <- read.table("../beet_nema_acc.core.pap.pro", sep = "\t", fill = T)

clean_inputs <- function(df){
  df_clean <- df[grepl("(\\+|-)", df$V5), ]
  df_clean <- df_clean[grepl("*M$", df_clean$V12), ]
  df_clean <- df_clean[grepl("^D|^S|^E", df_clean$V1), ]
  df_clean <- df_clean[grepl("^papilloma.", df_clean$V6), ]
  
  
  
  df_clean$V13 <- as.numeric(df_clean$V11)
  df_clean$V4 <- as.numeric(df_clean$V4)
  
  # df_fil <- filter(df_clean, V13 < 0.000001)
  df_fil <-  df_clean
  df_fil$V14 <- sub("_.*", "", df_fil$V1)
  # df_fil$V10 <- as.numeric(df_fil$V10)
  # df_fil$V7 <- as.numeric(df_fil$V7)
  # df_fil$V8 <- as.numeric(df_fil$V8)
  # df_fil$V9 <- as.numeric(df_fil$V9)
  # df_fil$V11 <- as.numeric(df_fil$V11)
  
  df_fil[,7:11] <- sapply(df_fil[,7:11],as.numeric)
  
  
  #This is remnants from pre-filtered pro file, but doesnt hurt to leave it I guess
  df_fil <- df_fil[!grepl("LQ465342", df_fil$V6),]
  df_fil <- df_fil[!grepl("OP721220", df_fil$V6),]
  df_fil <- df_fil[!grepl("AX020928|AX020938|OP721220|OP721198|LQ465342|JA802371|AX020936|AX020945|HE984563|MP134365|AF001600|AX020925|MZ189212", df_fil$V6),]
  df_fil$cov <- abs(df_fil$V7-df_fil$V8)
  df_fil$per_cov <- df_fil$cov/df_fil$V9
  # df_fil <- subset(df_fil, V4 > 200)
  df_fil <- df_fil[order(-df_fil$V10),]
  
  return(as.data.frame(df_fil))
}

beet_nema_clean <- clean_inputs(beet_nema)
beet_nema_clean_table <- as.data.frame(table(beet_nema_clean$V6))
beet_nema_clean <- beet_nema_clean[!grepl("ERR56356*", beet_nema_clean$V14),]

clean_400k <- read.table("clean.core.pap_400k.pro", sep = "\t", fill = T)
clean_400k_1 <- clean_inputs(clean_400k)


filter_inputs <- function(df, eval, contig){

df_fil <- filter(df, V13 < eval)
df_fil <- subset(df_fil, V4 > contig)
df_fil_sliced <- df_fil %>% group_by(V1) %>% slice_max(order_by = V10, n = 1)
return(df_fil_sliced)
}

clean_400k_tab <- as.data.frame(table(clean_400k_1$V6))
colnames(clean_400k_tab) <- c("V6", "freq")

sum <- clean_400k_1 %>%
  group_by(V6) %>%
  summarise_at(vars(V10), list(mean_PID = mean))

sum_2 <- clean_400k_1 %>%
  group_by(V6) %>%
  summarise_at(vars(V11), list(mean_eval = mean))

sum_3 <- left_join(sum, sum_2, by = "V6")
sum_stats_unfiltered <- left_join(sum_3, clean_400k_tab, by = "V6")


deca <- read.table("deca_clean.core.pap.pro", sep = "\t", fill = T)

clean_deca <- clean_inputs(deca)
clean_deca_fil <- filter_inputs(clean_deca, 0.00000001, 200)

###########turtles
turtle_lib <- read.table("turtle_hits.pro", sep = "\t", fill = T)
turtle_clean <- clean_inputs(turtle_lib)
turtle_clean_fil <- filter_inputs(turtle_clean, 0.00000001, 100)


####gecko
gecko_lib <- read.table("gecko_grep.pro", sep = "\t", fill = T)
gecko_clean <- clean_inputs(gecko_lib)
gecko_clean_fil <- filter_inputs(gecko_clean, 0.00000001, 100)

###try a map
library(ggplot2)
library(gggenes)

genes <- read.table("map.txt", header = T)

ggplot(genes, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, " qmm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette ="Set1") +
  theme_genes()

large <- read.table("large_contigs.pro", sep = "\t", fill = T)
large_clean <- clean_inputs(large)
large_clean_fil <- filter_inputs(large_clean, 0.00000001, 100)

salmonid <- read.table("salmonid_hits.pro", sep = "\t", fill = T)
salmonid_clean <- clean_inputs(salmonid)
salmonid_clean_fil <- filter_inputs(salmonid_clean, 0.00000001, 200)

genes_deer <- read.table("map_2.txt", header = T)

ggplot(genes_deer, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette ="Set1") +
  theme_genes()

genes_mykiss <- read.table("map_3.txt", header = T)

ggplot(genes_mykiss, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette ="Set1") +
  theme_genes()


mykiss <- read.table("seq.txt", header = F)
mykiss_urr <- substr(mykiss, 2280, 2435)


meow <- c("meow")

clean_anello <- read.table("anello2.pro", sep = "\t", fill = T)
clean_anello <- filter(clean_anello, V11 < 0.00000001 )
clean_anello <- filter(clean_anello, V10 < 90 )

##new PV
pv_2 <- read.table("updated_PV_nucs.domtbl", sep = "\t", fill = T)
salmonid_clean <- clean_inputs(salmonid)
salmonid_clean_fil <- filter_inputs(salmonid_clean, 0.00000001, 200)

##3400 accession test run
july_3400_test <- read.table("july1_3400_run.pro", sep = "\t", fill = T)
eval_fil_july_3400_test <- filter(july_3400_test, V11 < 0.000000001 )
july_3400_test$V14 <- sub("_.*", "", july_3400_test$V1)


###logan 2####
concat_pv2_fil_pv_l1 <- read.table("concat_pv2_fil_pv_l1.pro", sep = "\t", fill = T)

######################FOR ANNOTATING DIAMOND OUTPUT COLUMNS WITH NAMES######################
colnames(concat_pv2_fil_pv_l1) <- c("contig", "start", "end", "length", "strand", "model_name", "model_start", "model_end", "model_length", "per_iden", "e_val", "cigar", "aa", "nuc")
concat_pv2_fil_pv_l1$library_source <- sub("_.*", "", concat_pv2_fil_pv_l1$contig)

concat_pv2_fil_pv_l1_90 <- filter(concat_pv2_fil_pv_l1, per_iden < 90)

anello_pv2 <- read.table("anellos_fil.pro", sep = "\t", fill = T)
colnames(anello_pv2) <- c("contig", "start", "end", "length", "strand", "model_name", "model_start", "model_end", "model_length", "per_iden", "e_val", "cigar", "aa", "nuc")
anello_pv2$library_source <- sub("_.*", "", anello_pv2$contig)
anello_library <- as.data.frame(unique(anello_pv2$library_source))
write.table(paste0(anello_library$`unique(anello_pv2$library_source)`, collapse = ","), "anello_library.txt", quote = F, row.names = F, col.names = F, sep = ",")

anello_library_10k <- anello_library[1:10000,]
write.table(anello_library_10k, "anello_library_10k.txt", quote = F, row.names = F, col.names = F, sep = ",")

anello_library_20k <- anello_library[10001:20000,]
write.table(anello_library_20k, "anello_library_20k.txt", quote = F, row.names = F, col.names = F, sep = ",")

anello_library_30k <- anello_library[20001:30000,]
write.table(anello_library_30k, "anello_library_30k.txt", quote = F, row.names = F, col.names = F, sep = ",")

anello_library_40k <- anello_library[30001:40000,]
write.table(anello_library_40k, "anello_library_40k.txt", quote = F, row.names = F, col.names = F, sep = ",")

anello_library_final <- anello_library[40001:45998,]
write.table(anello_library_final, "anello_library_final", quote = F, row.names = F, col.names = F, sep = ",")

anello_pv2_ERR10301300 <- filter(anello_pv2, library_source == "ERR10301300")


# Install the SRAdb package if it is not installed.
if (!("SRAdb" %in% installed.packages())) {
  BiocManager::install("SRAdb")
}

# Declare directory to hold SRAdb file
# sra_db_dr <- <DIRECTORY NAME HERE>

# Create the directory where the SRAdb file will be downloaded to
# if (!dir.exists(sra_db_dir)) {
#   dir.create(sra_db_dir, recursive = TRUE)
# }

# Downloading this file will take some time, so we will only download it if 
# it doesn't exist.
# if (!file.exists(sqlite_file)) {
#   # Use SRAdb's function to download the file
#   SRAdb::getSRAdbFile(destdir = sra_db_dir)
# }
anello_library
anello_library_random1k <- anello_library[sample(nrow(anello_library), 1000), ]
write.table(anello_library_random1k, "anello_library_random1k.txt", quote = F, row.names = F, col.names = F, sep = ",")

concat_pv2_fil_pv <- read.table("concat_pv2_fil_pv.pro", sep = "\t", fill = T)
