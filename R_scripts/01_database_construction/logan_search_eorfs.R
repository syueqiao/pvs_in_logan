#interpret the output of logan search
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggmosaic)


map_date_df <- function(file_name, coloor){
  output_df <- read.table(file_name, sep = "\t", header = T, fill = T, quote = "\"")
  
  output_df$collection_date_sam_parsed <- parse_datetime(gsub("\\[|\\]|'", "", output_df$collection_date_sam))
  
  output_df$collection_date_sam_parsed <- date(output_df$collection_date_sam_parsed)
  world_map <- map_data("world")
  
  map <- output_df %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), data = world_map, fill = "lightgrey", color = "lightgrey") +
    geom_point(aes(x = longitude, y = latitude), color = coloor, size = 0.6, alpha = 0.4) +
    scale_color_identity() +
    coord_fixed() +
    xlab("") +
    ylab("") + theme_bw() + ggtitle(paste0(file_name))
  
  
  dates <- output_df %>%
    ggplot() +
    geom_histogram(aes(x=collection_date_sam_parsed), binwidth = 10) +
    scale_x_date(breaks="2 years") +
    ggtitle(paste0(file_name)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  png_name <-  sub('\\..*$', '', basename(file_name))


  ggplot2::ggsave(filename = paste0("../e6_e7/map_",png_name,".png"),map, width = 15, height = 10, units = "cm")
  ggplot2::ggsave(filename = paste0("../e6_e7/dates_",png_name,".png"),dates, width = 15, height = 10, units = "cm")
  
  out_list <- list(map, dates, output_df)
  return(out_list)

}

nsp7 <- map_datemap_datemap_date_df("../e6_e7/NSP7_SARSCOV2.tsv", "magenta")
nsp7_pre2020 <- filter(nsp7[[3]], collection_date_sam_parsed < '2019-01-01')
write.table(nsp7_pre2020, "../nsp7_pre2020.tsv", sep = "\t", quote = F, row.names = F)

nsp7 <- map_date_df("../e6_e7/NSP7_SARSCOV2.tsv", "magenta")
nsp7[[1]]
nsp7[[2]]
a <- nsp7[[3]]


tsvs_list <- list.files("../e6_e7", pattern = "*tsv$", full.names = T)

colors_unique <- c("blue1", "blue4", "red3", "green3", "green4", "orange4","pink4", "purple2", "aquamarine4", "darkgreen", "greenyellow", "cadetblue4", "cadetblue3", "cyan4", "royalblue")

for (j in 1:length(tsvs_list)) {
  x = paste0(tsvs_list[j])
  y = paste0(colors_unique[j])
  map_date_df(x, y)
}

datalist = list()

for (j in tsvs_list) {
  
  dat <- read.table(j, sep = "\t", header = T, fill = T, quote = "\"")
  datalist[[j]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)
big_data <- tibble::rownames_to_column(big_data, "source")
big_data$type <- with(big_data, ifelse(grepl("*C6L*|*K6L*|*VAC*", source), "Vaccinia", 
                                       ifelse(grepl("*UL145*|*US12*", source), "Herpesvirus", 
                                              ifelse(grepl("*NS*", source), "Coronavirus", "NA"))))
big_data$source <- gsub("*.?/", "", big_data$source)
big_data$source <- gsub(".e6_e", "", big_data$source)
big_data$source <- gsub("\\.tsv.*", "", big_data$source)


all <- ggplot(big_data) +
  geom_bar(aes(x = source, fill = librarysource),
           position = "stack",
           stat = "count") +
           scale_fill_brewer(palette = "Set3") +
           facet_grid(~ type, scales = "free_x") + 
           theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot2::ggsave(filename = "all_librarysource.png",all, width = 40, height = 20, units = "cm")


big_data_no_covid <- filter(big_data, !source == "NSP7_SARSCOV2")


no_covid <- ggplot(big_data_no_covid) +
  geom_bar(aes(x = source, fill = librarysource),
           position = "stack",
           stat = "count") +
  scale_fill_brewer(palette = "Set3") +
  facet_grid(~ type, scales = "free_x") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot2::ggsave(filename = "all_librarysource_nocovid.png",no_covid, width = 40, height = 20, units = "cm")

organism_table <- as.data.frame(table(big_data$organism))
tissue_table <- as.data.frame(table(big_data$tissue))
tissue_text_table <- as.data.frame(table(big_data$tissue_text))
disease_source_table <- as.data.frame(table(big_data$disease_source))
disease_text_table <- as.data.frame(table(big_data$disease_text))

#look at top for each ORF
orf_list <- unique(big_data$source)
org_list <- list()
for (i in orf_list) {
  prep <- filter(big_data, source == i)
  dat_table <- as.data.frame(table(prep$tissue))
dat_top <- top_n(dat_table, 5)
  
  org_list[[i]] <- dat_top # add it to your list
}
org_freq <- do.call(rbind, org_list)
org_freq <- tibble::rownames_to_column(org_freq, "source")
org_freq$source <- gsub("\\..*", "", org_freq$source)


P65 = createPalette(65,  c("#ff0000", "#00ff00", "#0000ff"))
names(P65) <- NULL


plot2 <- ggplot(org_freq, aes(x = source, fill = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', position = "stack") +
  scale_fill_manual(values = P65) +
  geom_text_repel(aes(label = ifelse(source == 'NSP7_SARSCOV2', as.character(Var1), "")), min.segment.length = 0, force = 15,
                  position=position_stack(vjust = 0.5), hjust = -1, vjust = -1,
                  direction="y", na.rm=TRUE, size = 5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = 'bottom') + theme(legend.key.size = unit(0.05, "cm"))
    

ggplot2::ggsave(filename = "../e6_e7/tissue_plot_covid.png",plot2, width = 80, height = 40, units = "cm")

org_freq_nocovid <- filter(org_freq, !source == 'NSP7_SARSCOV2')

plot3 <- ggplot(org_freq_nocovid, aes(x = source, fill = Var1, y = Freq)) + 
  geom_bar(stat = 'identity', position = "stack") +
  scale_fill_manual(values = P65) +
  geom_text_repel(aes(label = Var1), min.segment.length = 0, force = 1,
                  position=position_stack(vjust = 0), hjust = 0, vjust = 0,
                  direction="y", na.rm=TRUE, size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = 'bottom') + theme(legend.key.size = unit(0.05, "cm"))


ggplot2::ggsave(filename = "../e6_e7/tissue_plot_nocovid.png", plot3, width = 80, height = 40, units = "cm")


cell_lines <- big_data |> 
  filter(if_any(.cols = everything(), ~ grepl("cell line|cultured|vitro|HEK293|HeLa|Sf9|CHO|MCF-7|HL 60|Vero", .)))

write.table(big_data, "../all_logan_search_outputs.tsv", sep = "\t", quote = F, row.names = F)

####
map_date_df_deposit <- function(file_name, coloor){
  output_df <- read.table(file_name, sep = "\t", header = T, fill = T, quote = "\"")
  
  output_df$release_date_sam_parsed <- parse_datetime(gsub("\\[|\\]|'", "", output_df$releasedate))
  
  output_df$release_date_sam_parsed <- date(output_df$release_date_sam_parsed)
  world_map <- map_data("world")
  
  map <- output_df %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), data = world_map, fill = "lightgrey", color = "lightgrey") +
    geom_point(aes(x = longitude, y = latitude), color = coloor, size = 0.6, alpha = 0.4) +
    scale_color_identity() +
    coord_fixed() +
    xlab("") +
    ylab("") + theme_bw() + ggtitle(paste0(file_name))
  
  
  dates <- output_df %>%
    ggplot() +
    geom_histogram(aes(x=release_date_sam_parsed), binwidth = 10) +
    scale_x_date(breaks="2 years") +
    ggtitle(paste0(file_name)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  png_name <-  sub('\\..*$', '', basename(file_name))
  
  
  ggplot2::ggsave(filename = paste0("../e6_e7/map_",png_name,".png"),map, width = 15, height = 10, units = "cm")
  ggplot2::ggsave(filename = paste0("../e6_e7/dates_",png_name,".png"),dates, width = 15, height = 10, units = "cm")
  
  out_list <- list(map, dates, output_df)
  return(out_list)
  
}

nsp7_release <- map_date_df_deposit("../e6_e7/NSP7_SARSCOV2.tsv", "magenta")
p2 <- nsp7_release[[2]]
p2
df_release <- nsp7_release[[3]]

#import sprql query output
doid_cell_prolif <- read.table("../e6_e7/do_sparql_results.csv", sep = ",", header = T)
#parse ID to be in the same format of logan search output
doid_cell_prolif$id_parse <- gsub(":", "_", doid_cell_prolif$id)
#see which are cancer related
nsp7_cancer <- filter(df_release, do_id %in% doid_cell_prolif$id_parse)

#try for another one?
us12_hhv <- map_date_df_deposit("../e6_e7/US12_HHV1.tsv", "magenta")
us12_hhv_release <- us12_hhv[[3]]
us12_hhv_release_cancer <- filter(us12_hhv_release, do_id %in% doid_cell_prolif$id_parse)

#using the list of eORFs, randomly select those that werent included in the list of 59 ORFs that were possibly oncogenic
df_6922_0_Helicobacter_pylori <- map_date_df_deposit("../e6_e7/random_sel/6922_0_Helicobacter_pylori.tsv", "darkolivegreen")
df_6922_0_Helicobacter_pylori[[1]]
df_6922_0_Helicobacter_pylori[[2]]
df_6922 <- df_6922_0_Helicobacter_pylori[[3]]

df_6922_cancer <- filter(df_6922, do_id %in% doid_cell_prolif$id_parse)
df_6922_cancer_ratio <- (length(df_6922$ID)-length(df_6922_cancer$ID))/length(df_6922_cancer$ID)

#conduct chi squared between two input dataframes

compare_cancer_no_cancer <- function(df1, df2){
  
  oncogenic_df_cancer <- filter(df1, do_id %in% doid_cell_prolif$id_parse)
  oncogenic_df_notcancer <- filter(df1, !do_id %in% doid_cell_prolif$id_parse)
  oncogenic_df_cancer_count <- length(unique(oncogenic_df_cancer$ID))
  oncogenic_df_notcancer_count <- length(unique(oncogenic_df_notcancer$ID))
  
  nononco_df_cancer <- filter(df2, do_id %in% doid_cell_prolif$id_parse)
  nononco_df_notcancer <- filter(df2, !do_id %in% doid_cell_prolif$id_parse)
  nononco_df_cancer_count <- length(unique(nononco_df_cancer$ID))
  nononco_df_notcancer_count <- length(unique(nononco_df_notcancer$ID))
  
  

dat <- t(data.frame(
  "pred_onco" = c(oncogenic_df_notcancer_count, oncogenic_df_cancer_count),
  "notpred_onco" = c(nononco_df_notcancer_count, nononco_df_cancer_count),
  row.names = c("non_cancerous_doid", "cancerous_doid"),
  stringsAsFactors = FALSE
))


fisher_out <- fisher.test(dat)

viz <- mosaicplot(dat,
           main = paste0(deparse(substitute(df1))," and ",deparse(substitute(df2))," Mosaic plot"),
           color = 3:2)

out_list <- list(dat, fisher_out, viz)


return(out_list)

}

test_chi <- compare_cancer_no_cancer(us12_hhv_release, df_6922)
test_chi[[1]]
test_chi[[2]]$p.value

list_7074 <- map_date_df_deposit("../e6_e7/random_sel/7074_0_Streptococcus_pneumoniae.tsv", "green")
df_7074 <- list_7074[[3]]

list_4697 <- map_date_df_deposit("../e6_e7/random_sel/4697_WP_005158295.1_Yersinia_enterocolitica_subsp._palearctica_Y11.tsv", "green")
df_4697 <- list_4697[[3]]

test_chi <- compare_cancer_no_cancer(big_data_no_covid, big_non_onco)
test_chi[[1]]

full_compare_func <- function(file1, file2){
  file1_name <-  sub('\\..*$', '', basename(file1))
  file2_name <-  sub('\\..*$', '', basename(file2))
  
  
  file1_name_df <- map_date_df_deposit(file1, "blue")
  df1_file1 <- file1_name_df[[3]]
  
  file2_name_df <- map_date_df_deposit(file2, "yellow")
  df2_file2 <- file2_name_df[[3]]
  
  test_chi <- compare_cancer_no_cancer(df1_file1, df2_file2)

  ggplot_df <- test_chi[[1]]

mosaicplot(ggplot_df,
                  main = paste0(file1_name," (onco)", " and ",file2_name," Mosaic plot"),
                  color = 3:2)


return(test_chi)
}

# hhv5_ul34_vs_strep <- full_compare_func("../e6_e7/654_UL34_Human_herpesvirus_5.tsv", "../e6_e7/random_sel/7074_0_Streptococcus_pneumoniae.tsv")
# hhv5_ul34_vs_strep[[1]]
# hhv5_ul34_vs_strep[[2]]$p.value


NS3_MERSCOV_vs_pseud <- full_compare_func("../e6_e7/NS3_MERSCOV.tsv", "../e6_e7/random_sel/4815_WP_040259600.1_Pseudomonas_sp._Pseudomonas_massiliensis.tsv")
NS3_MERSCOV_vs_pseud[[1]]
NS3_MERSCOV_vs_pseud[[2]]$p.value

NS3_MERSCOV_vs_7074 <- full_compare_func("../e6_e7/random_sel/4633_WP_040232101.1_Citrobacter_pasteurii_Clermont_et_al._2015.tsv", "../e6_e7/random_sel/4697_WP_005158295.1_Yersinia_enterocolitica_subsp._palearctica_Y11.tsv")
NS3_MERSCOV_vs_7074[[1]]
NS3_MERSCOV_vs_7074[[2]]$p.value

big_non_onco <- rbind(df_7074, df_4697)

NS3_MERSCOV_vs_7074 <- full_compare_func("../e6_e7/random_sel/4633_WP_040232101.1_Citrobacter_pasteurii_Clermont_et_al._2015.tsv", "../e6_e7/random_sel/6922_0_Helicobacter_pylori.tsv")

#what if we pool them all?
tsvs_list_rand <- list.files("../e6_e7/random_sel", pattern = "*tsv$", full.names = T)


datalist_rand = list()

for (j in tsvs_list_rand) {
  
  dat <- read.table(j, sep = "\t", header = T, fill = T, quote = "\"")
  datalist_rand[[j]] <- dat # add it to your list
}

big_data_rand = do.call(rbind, datalist_rand)
big_data_rand <- tibble::rownames_to_column(big_data_rand, "source")
big_data_rand$source <- gsub("*.?/", "", big_data_rand$source)
big_data_rand$source <- gsub(".e6_e", "", big_data_rand$source)
big_data_rand$source <- gsub("random_se", "", big_data_rand$source)
big_data_rand$source <- gsub("\\.tsv.*", "", big_data_rand$source)

test_chi <- compare_cancer_no_cancer(big_data_no_covid, big_data_rand)
test_chi[[1]]
test_chi[[2]]

# chisq.test(dat)$expected

dflist <- lapply(datalist_rand, function(x) {x$do_id_type <- ifelse(x$do_id %in% doid_cell_prolif$id_parse, "cancer_doid", "noncancer_doid");return(x)})

dflist <- lapply(dflist, function(x) {x$do_id_type <- ifelse(x$do_id %in% doid_cell_prolif$id_parse, "cancer_doid", "noncancer_doid");return(x)})
dflist <- purrr::imap(dflist, ~mutate(.x, source = .y))

levels_doid <- c("cancer_doid", "noncancer_doid")
rand_counts <- lapply(dflist, function(x) {as.data.frame(table(factor(x$do_id_type, levels = levels_doid)))})
rand_counts <- lapply(rand_counts, function(x) x[1, "Freq"]/x[2, "Freq"])

###for oncogenic orfs
dflist_orfs <- lapply(datalist, function(x) {x$do_id_type <- ifelse(x$do_id %in% doid_cell_prolif$id_parse, "cancer_doid", "noncancer_doid");return(x)})

dflist_orfs <- lapply(dflist_orfs, function(x) {x$do_id_type <- ifelse(x$do_id %in% doid_cell_prolif$id_parse, "cancer_doid", "noncancer_doid");return(x)})
dflist_orfs <- purrr::imap(dflist_orfs, ~mutate(.x, source = .y))

levels_doid <- c("cancer_doid", "noncancer_doid")
rand_counts_orfs <- lapply(dflist_orfs, function(x) {as.data.frame(table(factor(x$do_id_type, levels = levels_doid)))})
rand_counts_orfs <- lapply(rand_counts_orfs, function(x) x[1, "Freq"]/x[2, "Freq"])

rand_counts_orfs_df <- as.data.frame(unlist(rand_counts_orfs))
rand_counts_orfs_df <- tibble::rownames_to_column(rand_counts_orfs_df, "source")

rand_counts_df <- as.data.frame(unlist(rand_counts))
rand_counts_df <- tibble::rownames_to_column(rand_counts_df, "source")


all_df <- as.data.frame(mapply(c, rand_counts_orfs_df,rand_counts_df))
all_df$sel <- ifelse(grepl("random_sel",all_df$source),"non_onc","onc")
all_df$`unlist(rand_counts_orfs)` <- as.numeric(all_df$`unlist(rand_counts_orfs)`)
all_df$coolor <- ifelse(grepl("non_onc", all_df$sel),"#FA8072","#00C5CD")


all_df$source <- gsub("*.?/", "", all_df$source)
all_df$source <- gsub(".e6_e", "", all_df$source)
all_df$source <- gsub("random_se", "", all_df$source)
all_df$source <- gsub("\\.tsv.*", "", all_df$source)

x <- all_df[order(all_df$`unlist(rand_counts_orfs)`, decreasing = T),]
coolors <- x$coolor
x$source <- factor(x$source, levels = x$source)


a_coolor <- ifelse(grepl("non_onc", all_df$sel),"#FA8072","#00C5CD")
# Give every color an appropriate name

p <- ggplot(x, aes(x = source, y = `unlist(rand_counts_orfs)`, fill = sel)) +
      geom_bar(stat="identity") +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45, hjust=1, colour=coolors))
p

#calculate odds ratio for each compared to full set of do_ids, and then the set of do_ids covered by the +onco_orfs set
big_data_no_covid$type <- NULL
big_data_no_covid$source <- NULL

all_all_fr <- rbind(big_data_no_covid[,1:50], big_data_rand[,2:51])

big_data_no_covid_can <- filter(big_data_no_covid, do_id %in% doid_cell_prolif$id_parse)
big_data_no_covid_noncan <- filter(big_data_no_covid, !do_id %in% doid_cell_prolif$id_parse)

big_data_no_covid_can_len <- length(unique(big_data_no_covid_can$biosample))
big_data_no_covid_noncan_len <- length(unique(big_data_no_covid_noncan$biosample))


all_all_fr_cand <- filter(all_all_fr, do_id %in% doid_cell_prolif$id_parse)
all_all_fr_noncand <- filter(all_all_fr, !do_id %in% doid_cell_prolif$id_parse)

all_cand_len <- length(unique(all_all_fr_cand$biosample))
all_noncand_len <- length(unique(all_all_fr_noncand$biosample))
baseline_len <- c(all_cand_len, all_noncand_len)

onco_orf_baseline_len <- c(big_data_no_covid_can_len, big_data_no_covid_noncan_len)
#baseline
#count cand and noncand for each tsv from oncogenic and non-oncogenic

datalist_counts= list()

for (j in tsvs_list) {
  
  dat <- map_date_df_deposit(j, "blue")
  dat <- dat[[3]]
  dat_cancer <- filter(dat, do_id %in% doid_cell_prolif$id_parse)
  dat_noncancer <- filter(dat, !do_id %in% doid_cell_prolif$id_parse)
  dat_cancer_len <- length(unique(dat_cancer$biosample))
  dat_noncancer_len <- length(unique(dat_noncancer$biosample))
  vect <- c(dat_cancer_len, dat_noncancer_len)
  datalist_counts[[j]] <- vect # add it to your list
}

datalist_counts_rand = list()
for (j in tsvs_list_rand) {
  
  dat <- map_date_df_deposit(j, "blue")
  dat <- dat[[3]]
  dat_cancer <- filter(dat, do_id %in% doid_cell_prolif$id_parse)
  dat_noncancer <- filter(dat, !do_id %in% doid_cell_prolif$id_parse)
  dat_cancer_len <- length(unique(dat_cancer$biosample))
  dat_noncancer_len <- length(unique(dat_noncancer$biosample))
  vect <- c(dat_cancer_len, dat_noncancer_len)
  datalist_counts_rand[[j]] <- vect # add it to your list
}

biglist <- c(datalist_counts_rand, datalist_counts)

dflist_counts <- lapply(biglist, function(x) {rbind(as.matrix(t(x)), as.matrix(t(onco_orf_baseline_len)))})
dflist_counts_fish <- lapply(dflist_counts, function(x) {fisher.test(x)})

chisq.test(dflist_counts$`../e6_e7/VACWR007.tsv`)$expected
dflist_counts_df <- data.frame(name = names(dflist_counts_fish), do.call(rbind, dflist_counts_fish))

dflist_counts_df <- dflist_counts_df %>% 
  unnest(conf.int) %>% 
  group_by(name) %>% 
  mutate(key = row_number()) %>% 
  spread(key, conf.int)

dflist_counts_df$name <- gsub("*.?/", "", dflist_counts_df$name)
dflist_counts_df$name <- gsub(".e6_e", "", dflist_counts_df$name)
dflist_counts_df$name <- gsub("random_se", "", dflist_counts_df$name)
dflist_counts_df$name <- gsub("\\.tsv.*", "", dflist_counts_df$name)

dflist_counts_df$type <- ifelse(grepl("^[0-9]", dflist_counts_df$name), "rand", "onco")

ggplot(dflist_counts_df, aes(y=reorder(name, as.numeric(estimate)), x=log10(as.numeric(estimate)), label=name, color = type)) +
  geom_point(size=3, shape=19) +
  geom_errorbarh(aes(xmin=log10(as.numeric(dflist_counts_df$`1`)), xmax=log10(as.numeric(dflist_counts_df$`2`))), height=.3) +
  coord_fixed(ratio=.3) +
  geom_vline(xintercept=1, linetype='longdash') + 
  scale_y_discrete(label=abbreviate) +
  geom_text_repel(nudge_y = 0.2, size = 2, aes(label = signif(as.numeric(p.value), 3))) +
  theme_bw()

dflist_counts_df$p.value <- as.numeric(dflist_counts_df$p.value)
dflist_counts_df$estimate <- as.numeric(dflist_counts_df$estimate)
dflist_counts_df$null.value <- as.numeric(dflist_counts_df$null.value)
dflist_counts_df$alternative <- as.numeric(dflist_counts_df$alternative)
dflist_counts_df$method <- as.numeric(dflist_counts_df$method)
dflist_counts_df$data.name <- as.numeric(dflist_counts_df$data.name)



write.table(dflist_counts_df, "dflist_counts_df.tsv", sep = "\t", row.names = F, col.names = T)

#all doids
bs_d <- readLines("../e6_e7/biosample_disease")
sum(str_count(bs_d, "DOID_"))

#what source?

groups_list <- as.data.frame(data.table::fread("../e6_e7/all_group_tags.list.gz", header = F))
groups_list[groups_list == "name"] <- NA
groups_list <- groups_list %>%
  mutate(V2=lag(V2)) %>%
  na.omit()
source_big_data <- left_join(all_all_fr, groups_list, by = c("ID" = "V1"))
source_big_data <- source_big_data[!(is.na(source_big_data$do_id) | source_big_data$do_id==""), ]


tabb <- as.data.frame(table(source_big_data$V2))
bar(tabb$Freq, breaks = 25)
tabb[order(tabb$Var1, tabb$Freq),]

ggplot(data=tabb, aes(x=reorder(Var1, -Freq), y=Freq)) +
  geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
