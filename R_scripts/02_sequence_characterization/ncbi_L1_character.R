#characterizing the L1 sequences on NCBI
library(tidyverse)
library(gggenes)
library(ggplot2)

#maybe shouldve clustered at 90% nucleotide identity before moving foward?

#function for parsing hmmsearch outputs
source("00_utilities/hmmsearch_utils.R")

all_ncbi_pv_sto_sto <- hmmsearch_clean("files/updated_ncbi_all.domtbl")
#reparse some stuff based on new contig formatting
all_ncbi_pv_sto_sto$contig_orf <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", all_ncbi_pv_sto_sto$query_acc)
all_ncbi_pv_sto_sto$contig <- sub("_[^_]+$", "", all_ncbi_pv_sto_sto$contig_orf)
#find the start and stop positions of nucleotides
all_ncbi_pv_sto_sto$query_acc_clean <- sub("_$","",all_ncbi_pv_sto_sto$query_acc_clean)

all_ncbi_pv_sto_sto <- all_ncbi_pv_sto_sto %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_end"), "(.*)_([^_]+)$")
all_ncbi_pv_sto_sto <- all_ncbi_pv_sto_sto %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_start"), "(.*)_([^_]+)$")
all_ncbi_pv_sto_sto <- all_ncbi_pv_sto_sto %>% extract(query_acc_clean, into = c("query_acc_clean", "orf"), "(.*)_([^_]+)$")

#define direction of orf
all_ncbi_pv_sto_sto[25:26] <- lapply(all_ncbi_pv_sto_sto[,25:26], as.numeric) 
all_ncbi_pv_sto_sto$direction <- ifelse(all_ncbi_pv_sto_sto$orf_start > all_ncbi_pv_sto_sto$orf_end, "rev", "fw")

#split into different dfs for each pfam model
all_ncbi_pv_sto_sto_split <- split(all_ncbi_pv_sto_sto, all_ncbi_pv_sto_sto$pfam)

#start with those that have B and I domains
all_ncbi_pv_sto_sto_BI <- inner_join(all_ncbi_pv_sto_sto_split$L1_B_super5, all_ncbi_pv_sto_sto_split$L1_I_super5, by='query_acc_clean')
NCBI_BI_contigs <- all_ncbi_pv_sto_sto_BI$query_acc_clean
#this dataframe then has the orfs that have identifiable B and I domains
all_ncbi_pv_sto_sto_BI_orfs <- all_ncbi_pv_sto_sto_BI[all_ncbi_pv_sto_sto_BI$query_acc.x == all_ncbi_pv_sto_sto_BI$query_acc.y, ]

#look for genes that have weird patterns and visualize it

#################################first, orfs that have separated B and I (eg, on different ORFs)
all_ncbi_pv_sto_sto_BI_seporf <- all_ncbi_pv_sto_sto_BI[all_ncbi_pv_sto_sto_BI$query_acc.x != all_ncbi_pv_sto_sto_BI$query_acc.y,]

sep_orfs <- all_ncbi_pv_sto_sto_BI_seporf$query_acc_clean

all_ncbi_pv_sto_sto_B_coords <- filter(all_ncbi_pv_sto_sto_split$L1_B_super5, query_acc_clean %in% sep_orfs)
all_ncbi_pv_sto_sto_I_coords <- filter(all_ncbi_pv_sto_sto_split$L1_I_super5, query_acc_clean %in% sep_orfs)

all_ncbi_pv_sto_sto_BI_coords <- rbind(all_ncbi_pv_sto_sto_B_coords, all_ncbi_pv_sto_sto_I_coords)
# all_ncbi_pv_sto_sto_BI_coords <- select(all_ncbi_pv_sto_sto_BI_coords, query_acc_clean, orf_start, orf_end, pfam)

all_ncbi_pv_sto_sto_BI_coords$orf_start <- as.numeric(all_ncbi_pv_sto_sto_BI_coords$orf_start)
all_ncbi_pv_sto_sto_BI_coords$orf_end <- as.numeric(all_ncbi_pv_sto_sto_BI_coords$orf_end)
all_ncbi_pv_sto_sto_BI_coords$nuc_start <- all_ncbi_pv_sto_sto_BI_coords$envcoord_to*2 + (all_ncbi_pv_sto_sto_BI_coords$envcoord_from - 3) + all_ncbi_pv_sto_sto_BI_coords$orf_start
all_ncbi_pv_sto_sto_BI_coords$nuc_end <- (all_ncbi_pv_sto_sto_BI_coords$envcoord_from*3 - 1) + all_ncbi_pv_sto_sto_BI_coords$orf_start

ggplot(all_ncbi_pv_sto_sto_BI_coords, aes(xmin = orf_start, xmax = orf_end, y = query_acc_clean)) +
  geom_gene_arrow(aes(alpha = 0.5)) +
  geom_subgene_arrow(data = all_ncbi_pv_sto_sto_BI_coords,
                     aes(xmin = orf_start, xmax = orf_end, y = query_acc_clean, fill = pfam,
                         xsubmin = nuc_start, xsubmax = nuc_end), color="black", alpha=.7) +
  facet_wrap(~ query_acc_clean, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

ggsave("outputs/orfs_updated.pdf", width = 20, height = 20, limitsize = FALSE)

###############look for possible inversions
##split based on fw and rev
all_ncbi_pv_sto_sto_fr_split <- split(all_ncbi_pv_sto_sto, all_ncbi_pv_sto_sto$direction)

all_ncbi_pv_sto_sto_fr_split_contigs <- purrr::reduce(all_ncbi_pv_sto_sto_fr_split, dplyr::inner_join, by = 'contig')
#get contigs that are common between the two split dfs
fw_rev_contigs <- unique(all_ncbi_pv_sto_sto_fr_split$contig)

#start with forward contigs and start of regions
func1 <- function(x) head(x,1)   # if duplicate, use first value
all_ncbi_pv_sto_sto_order_start_fw <- aggregate(all_ncbi_pv_sto_sto_fr_split$fw[c("envcoord_from")], by=list(contig=all_ncbi_pv_sto_sto_fr_split$fw$contig,pfam=all_ncbi_pv_sto_sto_fr_split$fw$pfam), func1)

all_ncbi_pv_sto_sto_order_start_fw <- reshape(all_ncbi_pv_sto_sto_order_start_fw, idvar = "contig", timevar = "pfam", direction = "wide") %>%
  filter(., contig %in% NCBI_BI_contigs)
#warning, grab first B occurrence which is fine

#for end
func2 <- function(x) tail(x,1)   # if duplicate, use last value
all_ncbi_pv_sto_sto_order_end_fw <- aggregate(all_ncbi_pv_sto_sto_fr_split$fw[c("envcoord_to")], by=list(contig=all_ncbi_pv_sto_sto_fr_split$fw$contig,pfam=all_ncbi_pv_sto_sto_fr_split$fw$pfam), func2)

all_ncbi_pv_sto_sto_order_end_fw <- reshape(all_ncbi_pv_sto_sto_order_end_fw, idvar = "contig", timevar = "pfam", direction = "wide") %>%
  filter(., contig %in% NCBI_BI_contigs)

#look for regions that are/not in the right order
blacklist <- read.table("files/blacklist.txt", header = F)
all_ncbi_pv_sto_sto_order_start_fw_clean <- filter(all_ncbi_pv_sto_sto_order_start_fw, !contig %in% blacklist$V1)
all_ncbi_pv_sto_sto_order_CD_B <- filter(all_ncbi_pv_sto_sto_order_start_fw_clean, envcoord_from.L1_B_super5 > envcoord_from.L1_CD_super5)
all_ncbi_pv_sto_sto_order_E_CD <- filter(all_ncbi_pv_sto_sto_order_start_fw_clean, envcoord_from.L1_CD_super5 > envcoord_from.L1_E_super5)
all_ncbi_pv_sto_sto_order_F_E <- filter(all_ncbi_pv_sto_sto_order_start_fw_clean, envcoord_from.L1_E_super5 > envcoord_from.L1_F_super5)
all_ncbi_pv_sto_sto_order_GH_F <- filter(all_ncbi_pv_sto_sto_order_start_fw_clean, envcoord_from.L1_F_super5 > envcoord_from.L1_GH_super5)
all_ncbi_pv_sto_sto_order_I_GH <- filter(all_ncbi_pv_sto_sto_order_start_fw_clean, envcoord_from.L1_GH_super5 > envcoord_from.L1_I_super5)

##for reverse contigs
all_ncbi_pv_sto_sto_order_start_rev <- aggregate(all_ncbi_pv_sto_sto_fr_split$rev[c("envcoord_from")], by=list(contig=all_ncbi_pv_sto_sto_fr_split$rev$contig,pfam=all_ncbi_pv_sto_sto_fr_split$rev$pfam), func1)

all_ncbi_pv_sto_sto_order_start_rev <- reshape(all_ncbi_pv_sto_sto_order_start_rev, idvar = "contig", timevar = "pfam", direction = "wide") %>%
  filter(., contig %in% NCBI_BI_contigs)

all_ncbi_pv_sto_sto_order_end_rev <- aggregate(all_ncbi_pv_sto_sto_fr_split$rev[c("envcoord_to")], by=list(contig=all_ncbi_pv_sto_sto_fr_split$rev$contig,pfam=all_ncbi_pv_sto_sto_fr_split$rev$pfam), func2)

all_ncbi_pv_sto_sto_order_end_rev <- reshape(all_ncbi_pv_sto_sto_order_end_rev, idvar = "contig", timevar = "pfam", direction = "wide") %>%
  filter(., contig %in% NCBI_BI_contigs)

all_ncbi_pv_sto_sto_order_start_rev_clean <- filter(all_ncbi_pv_sto_sto_order_start_rev, !contig %in% blacklist$V1)
all_ncbi_pv_sto_sto_order_rev_CD_B <- filter(all_ncbi_pv_sto_sto_order_start_rev_clean, envcoord_from.L1_B_super5 > envcoord_from.L1_CD_super5)
all_ncbi_pv_sto_sto_order_rev_E_CD <- filter(all_ncbi_pv_sto_sto_order_start_rev_clean, envcoord_from.L1_CD_super5 > envcoord_from.L1_E_super5)
all_ncbi_pv_sto_sto_order_rev_F_E <- filter(all_ncbi_pv_sto_sto_order_start_rev_clean, envcoord_from.L1_E_super5 > envcoord_from.L1_F_super5)
all_ncbi_pv_sto_sto_order_rev_GH_F <- filter(all_ncbi_pv_sto_sto_order_start_rev_clean, envcoord_from.L1_F_super5 > envcoord_from.L1_GH_super5)
all_ncbi_pv_sto_sto_order_rev_I_GH <- filter(all_ncbi_pv_sto_sto_order_start_rev_clean, envcoord_from.L1_GH_super5 > envcoord_from.L1_I_super5)


##############################expand this concept to all the jellyroll regions, look for "coverage" 
all_ncbi_pv_sto_sto_bare <- select(all_ncbi_pv_sto_sto, contig, pfam)
all_ncbi_pv_sto_sto_bare <- all_ncbi_pv_sto_sto_bare %>%
  mutate(pres = 1)

func3 <- function(x) {ifelse(length(x)==1,as.character(x), length(x))}
all_ncbi_pv_sto_sto_bare_long <- aggregate(all_ncbi_pv_sto_sto_bare[c("pres")], by=list(contig=all_ncbi_pv_sto_sto_bare$contig,pfam=all_ncbi_pv_sto_sto_bare$pfam), func3)
all_ncbi_pv_sto_sto_bare_long$pres <- as.numeric(all_ncbi_pv_sto_sto_bare_long$pres)
all_ncbi_pv_sto_sto_bare_long <- filter(all_ncbi_pv_sto_sto_bare_long, !contig %in% blacklist$V1)


order_plot_hm <- function(long_df){
  long_df_peak_values <- long_df %>% 
  group_by(contig) %>% 
  slice_max(order_by = pres, n = 1, with_ties = F) %>% 
  rename(peaked_at = pfam) %>% 
  select(-pres)

number_of_peaks <- long_df_peak_values %>% 
  group_by(peaked_at) %>% 
  count() %>% 
  arrange(-n)

data_peak_values_reordered <- long_df_peak_values %>% 
  inner_join(number_of_peaks, by = "peaked_at") %>% 
  arrange(-n)  

my_data_reordered <- long_df %>% 
  inner_join(number_of_peaks, by = c("pfam" = "peaked_at")) %>% 
  mutate(pfam = reorder(pfam, -n)) %>%  # this reorders the columns
  select(-n) %>% 
  inner_join(data_peak_values_reordered, by = "contig") %>% 
  mutate(peaked_at = reorder(peaked_at, n)) %>% 
  mutate(order_rows = as.numeric(peaked_at)) %>% 
  mutate(contig = reorder(contig, order_rows)) # this reorders the rows by the "peaked_at" column. 

level_order <- c('L1_B_super5', 'L1_CD_super5', 'L1_E_super5', 'L1_F_super5', "L1_GH_super5", "L1_I_super5") 

p <- my_data_reordered %>% 
  ggplot(aes(x = factor(pfam, level=level_order), y = contig)) +
  geom_tile(aes(fill = pres, color = pfam), size = 0.5) +
  theme_bw() +
  theme(axis.ticks.y.left = element_blank())

return(p)
}

p <- order_plot_hm(all_ncbi_pv_sto_sto_bare_long)
ggsave("outputs/domains.pdf", width = 10, height = 2500, limitsize = FALSE)


#only visualize the ones with identifiable B and I domains
all_ncbi_pv_sto_sto_bare_long_BI <- filter(all_ncbi_pv_sto_sto_bare_long, contig %in% NCBI_BI_contigs)

all_ncbi_pv_sto_sto_bare_long_BI_only <- reshape(all_ncbi_pv_sto_sto_bare_long, idvar = "contig", timevar = "pfam", direction = "wide") %>%
  filter(., contig %in% NCBI_BI_contigs)
all_ncbi_pv_sto_sto_bare_long_BI_only[is.na(all_ncbi_pv_sto_sto_bare_long_BI_only)] <- 0

all_ncbi_pv_sto_sto_bare_long_BI_only_mat <- as.matrix(all_ncbi_pv_sto_sto_bare_long_BI_only)
rownames(all_ncbi_pv_sto_sto_bare_long_BI_only_mat) <- all_ncbi_pv_sto_sto_bare_long_BI_only$contig
all_ncbi_pv_sto_sto_bare_long_BI_only_mat <- all_ncbi_pv_sto_sto_bare_long_BI_only_mat[,2:7]

all_ncbi_pv_sto_sto_bare_long_BI_only_den <- as.dendrogram(hclust(d = dist(x = all_ncbi_pv_sto_sto_bare_long_BI_only_mat)))
order <- order.dendrogram(all_ncbi_pv_sto_sto_bare_long_BI_only_den)

all_ncbi_pv_sto_sto_bare_long_BI$contig <- factor(x = all_ncbi_pv_sto_sto_bare_long_BI$contig,
                               levels = all_ncbi_pv_sto_sto_bare_long_BI$contig[order], 
                               ordered = TRUE)

k <- all_ncbi_pv_sto_sto_bare_long_BI %>% 
  ggplot(aes(x = factor(pfam), y = contig)) +
  geom_tile(aes(fill = pres, color = pfam), size = 0.5) +
  theme_bw() +
  theme(axis.ticks.y.left = element_blank())

ggsave("outputs/domains_BI.pdf", width = 10, height = 1000, limitsize = FALSE)

##try this with all domains?
all_ncbi_pv_sto_sto_bare_long_all <- all_ncbi_pv_sto_sto_bare_long

all_ncbi_pv_sto_sto_bare_long_all <- reshape(all_ncbi_pv_sto_sto_bare_long_all, idvar = "contig", timevar = "pfam", direction = "wide")
all_ncbi_pv_sto_sto_bare_long_all[is.na(all_ncbi_pv_sto_sto_bare_long_all)] <- 0

all_ncbi_pv_sto_sto_bare_long_all_mat <- as.matrix(all_ncbi_pv_sto_sto_bare_long_all)
rownames(all_ncbi_pv_sto_sto_bare_long_all_mat) <- all_ncbi_pv_sto_sto_bare_long_all$contig
all_ncbi_pv_sto_sto_bare_long_all_mat <- all_ncbi_pv_sto_sto_bare_long_all_mat[,2:7]

all_ncbi_pv_sto_sto_bare_long_all_den <- as.dendrogram(hclust(d = dist(x = all_ncbi_pv_sto_sto_bare_long_all_mat)))
order_all <- order.dendrogram(all_ncbi_pv_sto_sto_bare_long_all_den)

all_ncbi_pv_sto_sto_bare_long_all$contig <- factor(x = all_ncbi_pv_sto_sto_bare_long_all$contig,
                                                  levels = all_ncbi_pv_sto_sto_bare_long_all$contig[order_all], 
                                                  ordered = TRUE)

l <- all_ncbi_pv_sto_sto_bare_long %>% 
  ggplot(aes(x = factor(pfam), y = contig)) +
  geom_tile(aes(fill = pres, color = pfam), size = 0.5) +
  theme_bw() +
  theme(axis.ticks.y.left = element_blank())

ggsave("outputs/domains_BI_all_dendro.pdf", width = 10, height = 2000, limitsize = FALSE)

#look at if B and I are missing using length based calculation
#first look at things with CD hits, and a possible B based on how much room there is in the orf/contig

all_ncbi_pv_sto_sto_30len_CD <- filter(all_ncbi_pv_sto_sto, qlen > 30 & envcoord_from > 15 & pfam == "L1_CD_super5" & !contig %in% all_ncbi_pv_sto_sto_split$L1_B_super5$contig & !contig %in% blacklist$V1)
write.table( all_ncbi_pv_sto_sto_30len_CD$contig, "outputs/ncbi_B_misses.txt", quote = F, col.names = F, row.names = F)

ggplot(all_ncbi_pv_sto_sto_30len_CD, aes(xmin = 1, xmax = qlen, y = query_acc_clean)) +
  geom_gene_arrow(aes(alpha = 0.5)) +
  geom_subgene_arrow(data = all_ncbi_pv_sto_sto_30len_CD,
                     aes(xmin = 1, xmax = qlen, y = query_acc_clean, fill = pfam,
                         xsubmin = envcoord_from, xsubmax = envcoord_to), color="black", alpha=.7) +
  facet_wrap(~ query_acc_clean, scales = "free", ncol = 1) +
  scale_x_continuous(labels = scales::comma, limits = c(1, 510)) + 
  scale_fill_brewer(palette = "Set1") +
  theme_genes()

ggsave("outputs/domains_B_fail.pdf", width = 10, height = 20, limitsize = FALSE)

#now look at GH domains for "space" after the hit
all_ncbi_pv_sto_sto_30len_GH <- all_ncbi_pv_sto_sto
all_ncbi_pv_sto_sto_30len_GH$tail <- all_ncbi_pv_sto_sto_30len_GH$qlen - all_ncbi_pv_sto_sto_30len_GH$envcoord_to
all_ncbi_pv_sto_sto_30len_GH <- filter(all_ncbi_pv_sto_sto_30len_GH, qlen > 30 & tail > 15 & pfam == "L1_GH_super5" & !contig %in% all_ncbi_pv_sto_sto_split$L1_I_super5$contig & !contig %in% blacklist$V1)

ggplot(all_ncbi_pv_sto_sto_30len_GH, aes(xmin = 1, xmax = qlen, y = query_acc_clean)) +
  geom_gene_arrow(aes(alpha = 0.5)) +
  geom_subgene_arrow(data = all_ncbi_pv_sto_sto_30len_GH,
                     aes(xmin = 1, xmax = qlen, y = query_acc_clean, fill = pfam,
                         xsubmin = envcoord_from, xsubmax = envcoord_to), color="black", alpha=.7) +
  facet_wrap(~ query_acc_clean, scales = "free", ncol = 1) +
  scale_x_continuous(labels = scales::comma, limits = c(1, 510)) + 
  scale_fill_brewer(palette = "PRGn") +
  theme_genes()

ggsave("outputs/domains_I_fail.pdf", width = 10, height = 30, limitsize = FALSE)


#finally... after trials and tribulations, we can extract the B and I domain
#tackle the easy task first of getting the ones that have B and I hits on the same orf
all_ncbi_pv_sto_sto_BI_orfs_clean <- filter(all_ncbi_pv_sto_sto_BI_orfs, !contig.x %in% blacklist$V1 & !contig.y %in% blacklist$V1)
BI_contigs <- all_ncbi_pv_sto_sto_BI_orfs_clean$query_acc.x

all_ncbi_pv_sto_sto_BI_info <- filter(all_ncbi_pv_sto_sto, query_acc %in% BI_contigs)

p2_add <- c("BK066313.1_121_4075_5487_", "BK066841.1_119_4136_5596_", "BK066884.1_125_4254_5780_", "BK066894.1_125_3999_5450_", "MZ244208.1_35_175_1623_", "PP711991.1_136_5827_7023_")
all_ncbi_pv_sto_sto_p2 <- filter(all_ncbi_pv_sto_sto, query_acc %in% p2_add)

all_ncbi_pv_sto_sto_BI_info <- rbind(all_ncbi_pv_sto_sto_BI_info, all_ncbi_pv_sto_sto_p2)

all_ncbi_pv_sto_sto_L1_envcoord_from <- all_ncbi_pv_sto_sto_BI_info %>% group_by(query_acc) %>% slice_min(n = 1, envcoord_from)
all_ncbi_pv_sto_sto_L1_envcoord_from <- select(all_ncbi_pv_sto_sto_L1_envcoord_from, query_acc, envcoord_from)

all_ncbi_pv_sto_sto_L1_envcoord_to <- all_ncbi_pv_sto_sto_BI_info %>% group_by(query_acc) %>% slice_max(n = 1, envcoord_to)
all_ncbi_pv_sto_sto_L1_envcoord_to <- select(all_ncbi_pv_sto_sto_L1_envcoord_to, query_acc, envcoord_to)

all_ncbi_pv_sto_sto_L1_envcoords <- left_join(all_ncbi_pv_sto_sto_L1_envcoord_from, all_ncbi_pv_sto_sto_L1_envcoord_to, by = "query_acc")
all_ncbi_pv_sto_sto_L1_envcoords <- unique(all_ncbi_pv_sto_sto_L1_envcoords)
all_ncbi_pv_sto_sto_L1_envcoords_envaa <- select(all_ncbi_pv_sto_sto_L1_envcoords, query_acc, envcoord_from, envcoord_to)

all_ncbi_pv_sto_sto_L1_envcoords_envaa$nuc_from <- all_ncbi_pv_sto_sto_L1_envcoords_envaa$envcoord_from*2 + (all_ncbi_pv_sto_sto_L1_envcoords_envaa$envcoord_from - 3)
all_ncbi_pv_sto_sto_L1_envcoords_envaa$nuc_to <- all_ncbi_pv_sto_sto_L1_envcoords_envaa$envcoord_to*3 - 1
write.table(select(all_ncbi_pv_sto_sto_L1_envcoords_envaa, query_acc, nuc_from, nuc_to), "outputs/all_ncbi_pv_sto_sto_L1_envcoords_p12.txt", quote = F, col.names = F, row.names = F, sep = "\t")

all_ncbi_pv_sto_sto_L1_envcoords$query_acc <- sub("_$","",all_ncbi_pv_sto_sto_L1_envcoords$query_acc)
all_ncbi_pv_sto_sto_L1_envcoords <- all_ncbi_pv_sto_sto_L1_envcoords %>% extract(query_acc, into = c("query_acc", "orf_end"), "(.*)_([^_]+)$")
all_ncbi_pv_sto_sto_L1_envcoords <- all_ncbi_pv_sto_sto_L1_envcoords %>% extract(query_acc, into = c("query_acc", "orf_start"), "(.*)_([^_]+)$")
all_ncbi_pv_sto_sto_L1_envcoords <- all_ncbi_pv_sto_sto_L1_envcoords %>% extract(query_acc, into = c("query_acc", "orf"), "(.*)_([^_]+)$")

#next to get the sequences manually for the ones that have failed B and I hits
all_ncbi_pv_sto_sto_L1_envcoords_nuc <- select(all_ncbi_pv_sto_sto_L1_envcoords, query_acc, orf_start, orf_end)
write.table(all_ncbi_pv_sto_sto_L1_envcoords_nuc, "outputs/all_ncbi_pv_sto_sto_p12_nuc.txt", quote = F, col.names = F, row.names = F, sep = "\t")


write.table(unique(all_ncbi_pv_sto_sto_p2$query_acc), "outputs/all_ncbi_pv_sto_sto_p2.txt", quote = F, col.names = F, row.names = F, sep = "\t")

#finally get the spliced sequences manually (not in R)
group3 <- hmmsearch_clean("files/group3.domtbl")

group3_envcoord_from <- group3 %>% group_by(query_acc) %>% slice_min(n = 1, envcoord_from)
group3_envcoord_from <- select(group3_envcoord_from, query_acc, envcoord_from)

group3_envcoord_to <- group3 %>% group_by(query_acc) %>% slice_max(n = 1, envcoord_to)
group3_envcoord_to <- select(group3_envcoord_to, query_acc, envcoord_to)

group3_envcoords <- left_join(group3_envcoord_from, group3_envcoord_to, by = "query_acc")
group3_envcoords <- unique(group3_envcoords)
group3_envcoords_envaa <- select(group3_envcoords, query_acc, envcoord_from, envcoord_to)

group3_envcoords_envaa$nuc_from <- group3_envcoords_envaa$envcoord_from*2 + (group3_envcoords_envaa$envcoord_from - 3)
group3_envcoords_envaa$nuc_to <- group3_envcoords_envaa$envcoord_to*3 - 1
write.table(select(group3_envcoords_envaa, query_acc, nuc_from, nuc_to), "outputs/group3_envcoords_nuc.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(select(group3_envcoords_envaa, query_acc, envcoord_from, envcoord_to), "outputs/group3_envcoords_aa.txt", quote = F, col.names = F, row.names = F, sep = "\t")


group3_envcoords$query_acc <- sub("_$","",group3_envcoords$query_acc)
group3_envcoords <- group3_envcoords %>% extract(query_acc, into = c("query_acc", "orf_end"), "(.*)_([^_]+)$")
group3_envcoords <- group3_envcoords %>% extract(query_acc, into = c("query_acc", "orf_start"), "(.*)_([^_]+)$")
group3_envcoords <- group3_envcoords %>% extract(query_acc, into = c("query_acc", "orf"), "(.*)_([^_]+)$")

#extract nucleotide and aa seq of the ncbi input

##characterize the nt diversity in L1 for the accessions to get an updated view of the known diversity of PVs
L1_cluster_nt_accs <- read.table("files/ncbi_JR_centroids_acc.nt", sep = "\t", header = F)
L1_cluster_nt_accs$V1 <- sub(">","",L1_cluster_nt_accs$V1)
L1_cluster_nt_accs$V1 <- sub("\\..*","",L1_cluster_nt_accs$V1)

ncbi_tags_manual_nt <- ncbi_tags_manual 
ncbi_tags_manual_nt$Newick_label <- sub("-.*","", ncbi_tags_manual_nt$Newick_label)

L1_cluster_nt_accs_annot <- left_join(L1_cluster_nt_accs, ncbi_tags_manual_nt, by = c("V1" = "Newick_label"))
L1_cluster_nt_accs_annot <- L1_cluster_nt_accs_annot[!duplicated(L1_cluster_nt_accs_annot), ]
L1_cluster_nt_accs_annot_NAs <- L1_cluster_nt_accs_annot[rowSums(is.na(L1_cluster_nt_accs_annot)) > 0,]
write.table(L1_cluster_nt_accs_annot_NAs$V1, "outputs/L1_cluster_nt_accs_annot_NAs.txt", quote = F, col.names = F, row.names = F, sep = "\t")

L1_cluster_nt_accs_NAs_key <- read.table("files/L1_cluster_nt_acc_annot_NAs_key.txt", sep = "\t", header = F)
L1_cluster_nt_accs_NAs_key$V1 <- sub("\\..*","",L1_cluster_nt_accs_NAs_key$V1)

L1_cluster_nt_accs_annot_all <- left_join(L1_cluster_nt_accs_annot, L1_cluster_nt_accs_NAs_key, by = "V1")
L1_cluster_nt_accs_annot_all$final <- paste(L1_cluster_nt_accs_annot_all$status,L1_cluster_nt_accs_annot_all$V2)
L1_cluster_nt_accs_annot_all$final <- sub("NA | NA","",L1_cluster_nt_accs_annot_all$final)
L1_cluster_nt_accs_annot_all_graph <- select(L1_cluster_nt_accs_annot_all, V1, final)

L1_cluster_nt_accs_annot_all_graph_table <- as.data.frame(table(L1_cluster_nt_accs_annot_all_graph$final))
sum(L1_cluster_nt_accs_annot_all_graph_table$Freq) - 485


ncbi_tax_counts <- read.table("files/host_table_ncbi.txt", sep ="\t", header = T)
ncbi_tax_counts <- ncbi_tax_counts[1:57,]
agg_tax_counts <-aggregate(ncbi_tax_counts$Recorded.PVs, by=list(Category=ncbi_tax_counts$X), FUN=sum)
sum(ncbi_tax_counts$Recorded.PVs) - 485

write.table(agg_tax_counts, "outputs/host_table_ncbi_agg.txt", quote = F, col.names = F, row.names = F, sep = "\t")

##
#essentially redo this for the new output
#im sorry its SO ugly
updated_ncbi_all <- hmmsearch_clean("files/updated_ncbi_all.domtbl")
#reparse some stuff based on new contig formatting
updated_ncbi_all$contig_orf <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", updated_ncbi_all$query_acc)
updated_ncbi_all$contig <- sub("_[^_]+$", "", updated_ncbi_all$contig_orf)
#find the start and stop positions of nucleotides
updated_ncbi_all$query_acc_clean <- sub("_$","",updated_ncbi_all$query_acc_clean)

updated_ncbi_all <- updated_ncbi_all %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_end"), "(.*)_([^_]+)$")
updated_ncbi_all <- updated_ncbi_all %>% extract(query_acc_clean, into = c("query_acc_clean", "orf_start"), "(.*)_([^_]+)$")
updated_ncbi_all <- updated_ncbi_all %>% extract(query_acc_clean, into = c("query_acc_clean", "orf"), "(.*)_([^_]+)$")

#define direction of orf
updated_ncbi_all[25:26] <- lapply(updated_ncbi_all[,25:26], as.numeric) 
updated_ncbi_all$direction <- ifelse(updated_ncbi_all$orf_start > updated_ncbi_all$orf_end, "rev", "fw")

#split into different dfs for each pfam model
updated_ncbi_all_split <- split(updated_ncbi_all, updated_ncbi_all$pfam)

#start with those that have B and I domains
updated_ncbi_all_BI <- inner_join(updated_ncbi_all_split$ncbi_and_pave_L1_B, updated_ncbi_all_split$ncbi_and_pave_L1_I, by='query_acc_clean')
NCBI_BI_contigs <- updated_ncbi_all_BI$query_acc_clean

#this dataframe then has the orfs that have identifiable B and I domains
updated_ncbi_all_BI_orfs <- updated_ncbi_all_BI[updated_ncbi_all_BI$query_acc.x == updated_ncbi_all_BI$query_acc.y, ]

#extract the domains
# getting the ones that have B and I hits on the same orf
updated_ncbi_all_BI_orfs_clean <- filter(updated_ncbi_all_BI_orfs, !contig.x %in% blacklist$V1 & !contig.y %in% blacklist$V1)
BI_contigs <- updated_ncbi_all_BI_orfs_clean$query_acc.x

updated_ncbi_all_BI_info <- filter(updated_ncbi_all, query_acc %in% BI_contigs)

p2_add <- c("BK066313.1_121_4075_5487_", "BK066841.1_119_4136_5596_", "BK066884.1_125_4254_5780_", "BK066894.1_125_3999_5450_", "MZ244208.1_35_175_1623_", "PP711991.1_136_5827_7023_")
updated_ncbi_all_p2 <- filter(updated_ncbi_all, query_acc %in% p2_add)

updated_ncbi_all_BI_info <- rbind(updated_ncbi_all_BI_info, updated_ncbi_all_p2)

updated_ncbi_all_L1_envcoord_from <- updated_ncbi_all_BI_info %>% group_by(query_acc) %>% slice_min(n = 1, envcoord_from)
updated_ncbi_all_L1_envcoord_from <- select(updated_ncbi_all_L1_envcoord_from, query_acc, envcoord_from)

updated_ncbi_all_L1_envcoord_to <- updated_ncbi_all_BI_info %>% group_by(query_acc) %>% slice_max(n = 1, envcoord_to)
updated_ncbi_all_L1_envcoord_to <- select(updated_ncbi_all_L1_envcoord_to, query_acc, envcoord_to)

updated_ncbi_all_L1_envcoords <- left_join(updated_ncbi_all_L1_envcoord_from, updated_ncbi_all_L1_envcoord_to, by = "query_acc")
updated_ncbi_all_L1_envcoords <- unique(updated_ncbi_all_L1_envcoords)
updated_ncbi_all_L1_envcoords_envaa <- select(updated_ncbi_all_L1_envcoords, query_acc, envcoord_from, envcoord_to)

updated_ncbi_all_L1_envcoords_envaa$nuc_from <- updated_ncbi_all_L1_envcoords_envaa$envcoord_from*2 + (updated_ncbi_all_L1_envcoords_envaa$envcoord_from - 3)
updated_ncbi_all_L1_envcoords_envaa$nuc_to <- updated_ncbi_all_L1_envcoords_envaa$envcoord_to*3 - 1
write.table(select(updated_ncbi_all_L1_envcoords_envaa, query_acc, nuc_from, nuc_to), "outputs/updated_ncbi_all_L1_envcoords_p12.txt", quote = F, col.names = F, row.names = F, sep = "\t")

updated_ncbi_all_L1_envcoords$query_acc <- sub("_$","",updated_ncbi_all_L1_envcoords$query_acc)
updated_ncbi_all_L1_envcoords <- updated_ncbi_all_L1_envcoords %>% extract(query_acc, into = c("query_acc", "orf_end"), "(.*)_([^_]+)$")
updated_ncbi_all_L1_envcoords <- updated_ncbi_all_L1_envcoords %>% extract(query_acc, into = c("query_acc", "orf_start"), "(.*)_([^_]+)$")
updated_ncbi_all_L1_envcoords <- updated_ncbi_all_L1_envcoords %>% extract(query_acc, into = c("query_acc", "orf"), "(.*)_([^_]+)$")

###########next to get the sequences manually for the ones that have failed B and I hits###################
updated_ncbi_all_L1_envcoords_nuc <- select(updated_ncbi_all_L1_envcoords, query_acc, orf_start, orf_end)
write.table(updated_ncbi_all_L1_envcoords_nuc, "outputs/updated_ncbi_all_p12_nuc.txt", quote = F, col.names = F, row.names = F, sep = "\t")

#finally get the spliced sequences manually (not in R)
group3 <- hmmsearch_clean("files/group3.domtbl")

group3_envcoord_from <- group3 %>% group_by(query_acc) %>% slice_min(n = 1, envcoord_from)
group3_envcoord_from <- select(group3_envcoord_from, query_acc, envcoord_from)

group3_envcoord_to <- group3 %>% group_by(query_acc) %>% slice_max(n = 1, envcoord_to)
group3_envcoord_to <- select(group3_envcoord_to, query_acc, envcoord_to)

group3_envcoords <- left_join(group3_envcoord_from, group3_envcoord_to, by = "query_acc")
group3_envcoords <- unique(group3_envcoords)
group3_envcoords_envaa <- select(group3_envcoords, query_acc, envcoord_from, envcoord_to)

group3_envcoords_envaa$nuc_from <- group3_envcoords_envaa$envcoord_from*2 + (group3_envcoords_envaa$envcoord_from - 3)
group3_envcoords_envaa$nuc_to <- group3_envcoords_envaa$envcoord_to*3 - 1
write.table(select(group3_envcoords_envaa, query_acc, nuc_from, nuc_to), "outputs/group3_envcoords_nuc.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(select(group3_envcoords_envaa, query_acc, envcoord_from, envcoord_to), "outputs/group3_envcoords_aa.txt", quote = F, col.names = F, row.names = F, sep = "\t")


group3_envcoords$query_acc <- sub("_$","",group3_envcoords$query_acc)
group3_envcoords <- group3_envcoords %>% extract(query_acc, into = c("query_acc", "orf_end"), "(.*)_([^_]+)$")
group3_envcoords <- group3_envcoords %>% extract(query_acc, into = c("query_acc", "orf_start"), "(.*)_([^_]+)$")
group3_envcoords <- group3_envcoords %>% extract(query_acc, into = c("query_acc", "orf"), "(.*)_([^_]+)$")

#extract nucleotide and aa seq of the ncbi input

##characterize the nt diversity in L1 for the accessions to get an updated view of the known diversity of PVs
L1_cluster_nt_accs <- read.table("files/ncbi_JR_centroids_acc.nt", sep = "\t", header = F)
L1_cluster_nt_accs$V1 <- sub(">","",L1_cluster_nt_accs$V1)
L1_cluster_nt_accs$V1 <- sub("\\..*","",L1_cluster_nt_accs$V1)

ncbi_tags_manual_nt <- ncbi_tags_manual 
ncbi_tags_manual_nt$Newick_label <- sub("-.*","", ncbi_tags_manual_nt$Newick_label)

L1_cluster_nt_accs_annot <- left_join(L1_cluster_nt_accs, ncbi_tags_manual_nt, by = c("V1" = "Newick_label"))
L1_cluster_nt_accs_annot <- L1_cluster_nt_accs_annot[!duplicated(L1_cluster_nt_accs_annot), ]
L1_cluster_nt_accs_annot_NAs <- L1_cluster_nt_accs_annot[rowSums(is.na(L1_cluster_nt_accs_annot)) > 0,]
write.table(L1_cluster_nt_accs_annot_NAs$V1, "outputs/L1_cluster_nt_accs_annot_NAs.txt", quote = F, col.names = F, row.names = F, sep = "\t")

L1_cluster_nt_accs_NAs_key <- read.table("files/L1_cluster_nt_acc_annot_NAs_key.txt", sep = "\t", header = F)
L1_cluster_nt_accs_NAs_key$V1 <- sub("\\..*","",L1_cluster_nt_accs_NAs_key$V1)

L1_cluster_nt_accs_annot_all <- left_join(L1_cluster_nt_accs_annot, L1_cluster_nt_accs_NAs_key, by = "V1")
L1_cluster_nt_accs_annot_all$final <- paste(L1_cluster_nt_accs_annot_all$status,L1_cluster_nt_accs_annot_all$V2)
L1_cluster_nt_accs_annot_all$final <- sub("NA | NA","",L1_cluster_nt_accs_annot_all$final)
L1_cluster_nt_accs_annot_all_graph <- select(L1_cluster_nt_accs_annot_all, V1, final)

L1_cluster_nt_accs_annot_all_graph_table <- as.data.frame(table(L1_cluster_nt_accs_annot_all_graph$final))
sum(L1_cluster_nt_accs_annot_all_graph_table$Freq) - 485


ncbi_tax_counts <- read.table("files/host_table_ncbi.txt", sep ="\t", header = T)
ncbi_tax_counts <- ncbi_tax_counts[1:57,]
agg_tax_counts <-aggregate(ncbi_tax_counts$Recorded.PVs, by=list(Category=ncbi_tax_counts$X), FUN=sum)
sum(ncbi_tax_counts$Recorded.PVs) - 485

write.table(agg_tax_counts, "outputs/host_table_ncbi_agg.txt", quote = F, col.names = F, row.names = F, sep = "\t")