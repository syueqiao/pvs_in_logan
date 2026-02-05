#e6 e7 interrogation
library(tidyverse)

#function to parse the domtbl output
hmmsearch_clean <- function(input_domtbl){
  input_domtbl_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:22]  
  colnames(input_domtbl_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                   "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)
  return(input_domtbl_scan)
}


#next, tabulate which are "full"
feb_7_pv_L1_BI <- hmmsearch_clean("feb_7_pv_fil_form_L1_BI.domtbl")

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
# write.table(select(feb_7_pv_L1_BI_mat_all, library), "feb_7_pv_L1_BI_mat_full_library.list", quote = F, col.names = F, row.names = F)
# write.table(select(feb_7_pv_L1_BI_mat_all, query_acc), "feb_7_pv_L1_BI_mat_full.list", quote = F, col.names = F, row.names = F)
L1_full_contigs <- as.data.frame(feb_7_pv_L1_BI_mat_all$query_acc)
L1_full_contigs$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", L1_full_contigs$`feb_7_pv_L1_BI_mat_all$query_acc`)

#look at e6 e7 hits
feb_7_pv_e6_e7 <- hmmsearch_clean("table1_feb_7_pv_e6_e7.domtbl")
feb_7_pv_e6_e7$cov <- abs(feb_7_pv_e6_e7$envcoord_from - feb_7_pv_e6_e7$envcoord_to)/feb_7_pv_e6_e7$tlen
hist(feb_7_pv_e6_e7$cov)

feb_7_pv_e6_e7 <- filter(feb_7_pv_e6_e7, cov > 0.9)

feb_7_pv_e6_e7$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", feb_7_pv_e6_e7$query_acc)
feb_7_pv_e6_e7$library <- sub("_.*", "\\1", feb_7_pv_e6_e7$contig)

feb_7_pv_e6_e7 <- filter(feb_7_pv_e6_e7, score_full > 30)

feb_7_pv_e6_e7_full_L1s <- filter(feb_7_pv_e6_e7, library %in% feb_7_pv_L1_BI_mat_all$library)

feb_7_pv_e6_e7_full_L1s_slice_max_score <- feb_7_pv_e6_e7_full_L1s %>% 
  group_by(query_acc) %>% 
  slice_max(n = 1, score_one)

# write.table(select(feb_7_pv_e6_e7_full_L1s_slice_max_score, query_acc), "feb_7_pv_e6_e7_full_L1s_slice_max_score.list", quote = F, col.names = F, row.names = F)

#read in diamond pro file
l1_e6_e7_diamond <- read.table("l1_e6_e7_feb_7_pv_matches.tsv", sep = "\t", fill = T)
l1_e6_e7_diamond_eval_fil <- filter(l1_e6_e7_diamond, V11 < 0.000001 )

#get library and contig id
l1_e6_e7_diamond_eval_fil$contig <- sub("([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", l1_e6_e7_diamond_eval_fil$V1)
l1_e6_e7_diamond_eval_fil$library <- sub("_.*", "\\1", l1_e6_e7_diamond_eval_fil$contig)

l1_e6_e7_diamond_eval_fil$hit <- with(l1_e6_e7_diamond_eval_fil, ifelse(V6 == "NP_041325.1", "E6", ifelse(V6 == "NP_041326.1", "E7", "L1"))) 
l1_e6_e7_diamond_eval_fil_sub <- select(l1_e6_e7_diamond_eval_fil, V10, hit, library)

# install.packages("splitstackshape")
library(splitstackshape)

#some goddamn jank
l1_e6_e7_diamond_eval_fil_sub$accounting <- with(l1_e6_e7_diamond_eval_fil_sub, ifelse(hit == "E6", "0.1", ifelse(hit == "E7", "0.01", "0.001"))) 

l1_e6_e7_diamond_eval_fil_sub$sum <- ave(as.numeric(l1_e6_e7_diamond_eval_fil_sub$accounting), l1_e6_e7_diamond_eval_fil_sub$library, FUN=sum)

l1_e6_e7_diamond_eval_fil_sub_only1 <- filter(l1_e6_e7_diamond_eval_fil_sub, sum == 0.111)

l1_e6_e7_diamond_eval_fil_sub_only1_reshape <- reshape(l1_e6_e7_diamond_eval_fil_sub_only1, idvar="library", timevar="hit", direction="wide")

l1_e6_e7_diamond_eval_fil_sub_only1_reshape <- drop_na(l1_e6_e7_diamond_eval_fil_sub_only1_reshape)

lm_eqn <- function(df, x, y){
  m <- speedglm::speedlm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


l1_e6_e7_diamond_eval_fil_sub_only1_reshape <- filter(l1_e6_e7_diamond_eval_fil_sub_only1_reshape, library %in% feb_7_pv_e6_e7_full_L1s$library)


l1_e6_e7_diamond_eval_fil_sub_only1_reshape %>% ggplot(aes(x = V10.L1, y = V10.E6)) + 
  geom_point(size=2, alpha = 0.4)  +
  geom_smooth(method=lm , color="black", fill="#69b3a2", level=0.98) +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_bw(base_size = 20) +
  ggtitle("L1 vs E6") +
  theme(legend.position = "none") + 
  geom_text(x = 80, y = 90, label = lm_eqn(l1_e6_e7_diamond_eval_fil_sub_only1_reshape, l1_e6_e7_diamond_eval_fil_sub_only1_reshape$V10.L1, 
                                           l1_e6_e7_diamond_eval_fil_sub_only1_reshape$V10.E6), parse = TRUE)


l1_e6_e7_diamond_eval_fil_sub_only1_reshape %>% ggplot(aes(x = V10.L1, y = V10.E7)) + 
  geom_point(size=2, alpha = 0.4)  +
  geom_smooth(method=lm , color="black", fill="#69b3a2", se=TRUE) +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_bw(base_size = 20) +
  ggtitle("L1 vs E7") +
  theme(legend.position = "none") +
  geom_text(x = 80, y = 90, label = lm_eqn(l1_e6_e7_diamond_eval_fil_sub_only1_reshape, l1_e6_e7_diamond_eval_fil_sub_only1_reshape$V10.L1, 
                                           l1_e6_e7_diamond_eval_fil_sub_only1_reshape$V10.E7), parse = TRUE)


l1_e6_e7_diamond_eval_fil_sub_only1_reshape %>% ggplot(aes(x = V10.E6, y = V10.E7)) + 
  geom_point(size=2, alpha = 0.4)  +
  geom_smooth(method=lm , color="black", fill="#69b3a2", se=TRUE) +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_bw(base_size = 20) +
  ggtitle("E6 vs E7") +
  theme(legend.position = "none") +
  geom_text(x = 80, y = 90, label = lm_eqn(l1_e6_e7_diamond_eval_fil_sub_only1_reshape, l1_e6_e7_diamond_eval_fil_sub_only1_reshape$V10.E6, 
                                           l1_e6_e7_diamond_eval_fil_sub_only1_reshape$V10.E7), parse = TRUE)



l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100 <- filter(l1_e6_e7_diamond_eval_fil_sub_only1_reshape, V10.L1 < 99)
l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100 <- filter(l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100, V10.E7 < 99)
l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100 <- filter(l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100, V10.E6 < 99)

l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100 <- filter(l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100, library %in% feb_7_pv_e6_e7_full_L1s$library)



l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100 %>% ggplot(aes(x = V10.L1, y = V10.E6)) + 
  geom_point(size=2, alpha = 0.4)  +
  geom_smooth(method=lm , color="black", fill="#69b3a2", level=0.98) +
  xlim(25, 95) +
  ylim(25, 95) +
  theme_bw(base_size = 20) +
  ggtitle("L1 vs E6") +
  theme(legend.position = "none") +
  geom_text(x = 80, y = 90, label = lm_eqn(l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100, l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100$V10.L1, 
                                           l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100$V10.E6), parse = TRUE)

l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100 %>% ggplot(aes(x = V10.L1, y = V10.E7)) + 
  geom_point(size=2, alpha = 0.4)  +
  geom_smooth(method=lm , color="black", fill="#69b3a2", level=0.98) +
  xlim(25, 95) +
  ylim(25, 95) +
  theme_bw(base_size = 20) +
  ggtitle("L1 vs E7") +
  theme(legend.position = "none") +
  geom_text(x = 80, y = 90, label = lm_eqn(l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100, l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100$V10.L1, 
                                           l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100$V10.E7), parse = TRUE)


l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100 %>% ggplot(aes(x = V10.E6, y = V10.E7)) + 
  geom_point(size=2, alpha = 0.4)  +
  geom_smooth(method=lm , color="black", fill="#69b3a2", level=0.98) +
  xlim(25, 95) +
  ylim(25, 95) +
  theme_bw(base_size = 20) +
  ggtitle("E6 vs E7") +
  theme(legend.position = "none") +
  geom_text(x = 80, y = 90, label = lm_eqn(l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100, l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100$V10.E6, 
                                           l1_e6_e7_diamond_eval_fil_sub_only1_reshape_no100$V10.E7), parse = TRUE)

###redo some of this with contigs
l1_e6_e7_diamond_eval_fil_sub_contigs <- select(l1_e6_e7_diamond_eval_fil, V10, hit, contig, V1)
l1_e6_e7_diamond_eval_fil_sub_contigs_split <-split(l1_e6_e7_diamond_eval_fil_sub_contigs, l1_e6_e7_diamond_eval_fil_sub_contigs$hit)

# write.table(l1_e6_e7_diamond_eval_fil_sub_contigs_split$E6$V1, "E6s_for_evo.list", sep = "\t", quote = F, col.names = F, row.names = F)
# write.table(l1_e6_e7_diamond_eval_fil_sub_contigs_split$E7$V1, "E7s_for_evo.list", sep = "\t", quote = F, col.names = F, row.names = F)
# write.table(l1_e6_e7_diamond_eval_fil_sub_contigs_split$L1$V1, "L1s_for_evo.list", sep = "\t", quote = F, col.names = F, row.names = F)



df_list <- lapply(1:length(l1_e6_e7_diamond_eval_fil_sub_contigs_split), 
                  function(x) (pivot_wider(l1_e6_e7_diamond_eval_fil_sub_contigs_split[[x]], names_from = hit, values_from = V10 )))

df_contigs <- df_list %>% reduce(left_join, by = "contig")
df_contigs <- drop_na(df_contigs)

df_contigs[150, 2] <- c("61")
df_contigs$E6 <- as.numeric(df_contigs$E6)

evolution_scatter <- function(df, v1, v2, coloor1, coloor2){

v1name <- deparse(substitute(v1))
v2name <- deparse(substitute(v2))
dfname <- deparse(substitute(df))

v2name_text <- paste(v2name)

m <- speedglm::speedlm(df[[v2name]] ~ df[[v1name]], df, fitted = T);

pred <- fitted(m)

eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                 list(a = format(unname(coef(m)[1]), digits = 2),
                      b = format(unname(coef(m)[2]), digits = 2),
                      r2 = format(summary(m)$r.squared, digits = 3)))

p <-  df %>% ggplot(aes(x = {{ v1 }}, y = {{ v2 }})) + 
  geom_point(size=0.5, alpha = 0.2, color = coloor1) +
  geom_abline(slope = coef(m)[["df[[v1name]]"]], 
              intercept = coef(m)[["(Intercept)"]])+
  ylim(0, 110) +
  xlim(0, 110) +
  theme_bw(base_size = 20) +
  ggtitle(paste0(v1name, " vs ", v2name)) +
  theme(legend.position = "none") +
  geom_text(x = 25, y = 90, label = as.character(as.expression(eq)), parse = T)

ggplot2::ggsave(filename = paste0(v1name, "_", v2name, dfname, ".png") , p, width = 20, height = 20, units = "cm")
  
}

evolution_scatter(head(e6e7_all_all_sel_nozero_nona, 1000), pid_L1, pid_E7, "#6e3fba", "#3f8bba")

evolution_scatter(df_contigs, L1, E7, "#6e3fba", "#3f8bba")
evolution_scatter(df_contigs, L1, E6, "#ba3f4d", "#3f8bba")
evolution_scatter(df_contigs, E6, E7, "#3f4eba", "#3f8bba")



df_contigs_n100 <- filter(df_contigs, E7 < 99)
df_contigs_n100 <- filter(df_contigs_n100, E6 < 99)
df_contigs_n100 <- filter(df_contigs_n100, L1 < 99)

evolution_scatter(df_contigs_n100, L1, E7, "#6e3fba", "#3f8bba")
evolution_scatter(df_contigs_n100, L1, E6, "#ba3f4d", "#3f8bba")
evolution_scatter(df_contigs_n100, E6, E7, "#3f4eba", "#3f8bba")

#try all-all

e6_all_all <- read.table("E6_E6_all.aln", sep = "\t", fill = T)
colnames(e6_all_all) <- c("query", "target", "pid_E6", "score", "a", "b", "c", "d", "e", "f", "g", "h")
e7_all_all <- read.table("E7_E7_all.aln", sep = "\t", fill = T)
colnames(e7_all_all) <- c("query", "target", "pid_E7", "score", "a", "b", "c", "d", "e", "f", "g", "h")
l1_all_all <- read.table("L1_L1_all.aln", sep = "\t", fill = T)
colnames(l1_all_all) <- c("query", "target", "pid_L1", "score", "a", "b", "c", "d", "e", "f", "g", "h")




e6_all_all$hit <- c("E6")
e6_all_all$id <- paste0(e6_all_all$query, e6_all_all$target)
e6_all_all <- distinct(e6_all_all, id, .keep_all = T)


e7_all_all$hit <- c("E7")
e7_all_all$id <- paste0(e7_all_all$query, e7_all_all$target)

l1_all_all$hit <- c("L1")
l1_all_all$id <- paste0(l1_all_all$query, l1_all_all$target)

e6e7_all_all <- left_join(e7_all_all, e6_all_all, by = "id")
l1_all_all_sel <- select(l1_all_all, id, hit, pid_L1, target)
e6e7l1_all_all <- left_join(e6e7_all_all, l1_all_all_sel, by = "id")

e6e7_all_all_sel <- select(e6e7l1_all_all, id, pid_E6, pid_E7, pid_L1)

e6e7_all_all_sel_nona <- drop_na(e6e7_all_all_sel)

evolution_scatter(e6e7_all_all_sel_nona, pid_L1, pid_E7, "#6e3fba", "#3f8bba")

e6e7_all_all_sel_nozero <- filter(e6e7_all_all_sel, pid_L1 < 99)
e6e7_all_all_sel_nozero <- filter(e6e7_all_all_sel_nozero, pid_E7 < 99)
e6e7_all_all_sel_nozero <- filter(e6e7_all_all_sel_nozero, pid_E6 < 99)



ggplot2::ggsave(filename = "plot.png",p, width = 10, height = 10, units = "cm")

e6e7_all_all_sel_nozero_nona <- drop_na(e6e7_all_all_sel_nozero)
lm_eqn(e6e7_all_all_sel_nozero, e6e7_all_all_sel_nozero$pid_L1, e6e7_all_all_sel_nozero$pid_E7)

m <- speedglm::speedlm(e6e7_all_all_sel_nozero$pid_E7 ~ e6e7_all_all_sel_nozero$pid_L1, e6e7_all_all_sel_nozero)

summary(m)

evolution_scatter(e6e7_all_all_sel_nozero_nona, pid_L1, pid_E7, "#508004", "#3f8bba")
evolution_scatter(e6e7_all_all_sel_nozero_nona, pid_L1, pid_E6, "#807204", "#3f8bba")
evolution_scatter(e6e7_all_all_sel_nozero_nona, pid_E6, pid_E7, "#803404", "#3f8bba")

evolution_scatter(e6e7_all_all_sel_nona, pid_L1, pid_E7, "#b4cb8e", "#3f8bba")
evolution_scatter(e6e7_all_all_sel_nona, pid_L1, pid_E6, "#bcb471", "#3f8bba")
evolution_scatter(e6e7_all_all_sel_nona, pid_E6, pid_E7, "#bb9075", "#3f8bba")

cluster_lowe6 <- filter(e6e7_all_all_sel_nona, pid_L1 > 68 & pid_L1 < 75)
cluster_lowe6 <- filter(cluster_lowe6, pid_E6 > 30 & pid_E6 < 45)
cluster_lowe6_sep <- cluster_lowe6 %>%
  separate_wider_delim(id,
                       delim = stringr::regex("(?=[A-Z][R][R])"),
                       too_few = "align_start",
                       names_sep = '',
                       names_repair = ~ sub("value", "X", .x))

e6e7_all_all_sel_nona <- e6e7_all_all_sel_nona %>%
  separate_wider_delim(id,
                       delim = stringr::regex("(?=[A-Z][R][R])"),
                       too_few = "align_start",
                       names_sep = '',
                       names_repair = ~ sub("value", "X", .x))


m <- speedglm::speedlm(cluster_lowe6_sep$pid_E6 ~ cluster_lowe6_sep$pid_L1, cluster_lowe6_sep, fitted = T);

eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                 list(a = format(unname(coef(m)[1]), digits = 2),
                      b = format(unname(coef(m)[2]), digits = 2),
                      r2 = format(summary(m)$r.squared, digits = 3)))

library(Polychrome)

P800 = createPalette(800,  c("#ff0000", "#00ff00", "#0000ff"))
names(P800) <- NULL


p <-  cluster_lowe6_sep %>% ggplot(aes(x = pid_L1, y = pid_E6, color = id3)) + 
  geom_point(size=2, alpha = 0.05) +
  geom_abline(slope = coef(m)[["cluster_lowe6_sep$pid_L1"]], 
              intercept = coef(m)[["(Intercept)"]])+
  scale_color_manual(values = P800) +
  ylim(30, 50) +
  xlim(60, 80) +
  theme_bw(base_size = 20) +
  ggtitle("pid_L1 vs pid_E6 subset") +
  theme(legend.position = "none") +
  geom_text(x = 25, y = 90, label = as.character(as.expression(eq)), parse = T) 

ggplot2::ggsave(filename = "pp_id3.png" , p, width = 20, height = 20, units = "cm")


#parse and graph the output of the clusters
e6_clusters <- read.table("cluster_count_all_e6s.tsv", sep = "\t", header = T)
ggplot(e6_clusters, aes(x = iden, y = clusters)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() +
  ggtitle(label = "e6 clusters")


e7_clusters <- read.table("cluster_count_all_e7s.tsv", sep = "\t", header = T)
ggplot(e7_clusters, aes(x = iden, y = clusters)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() +
  ggtitle(label = "e7 clusters")


library(ggtree)
e6_99_tree <- read.tree("all_e6s.tree")
ggtree(e6_99_tree)
e6_99_tree_cols <- as.data.frame(e6_99_tree[["tip.label"]])
e6_99_tree_cols$tip.tabel <- e6_99_tree_cols$`e6_99_tree[["tip.label"]]`
e6_99_tree_cols <- separate(e6_99_tree_cols, col = `e6_99_tree[["tip.label"]]`, into = c("one", "two", "three", "four"), sep = "\\|")
e6_99_tree_cols <- separate(e6_99_tree_cols, col = one, into = c("five", "six", "seven", "eight"), sep = "_")
names <- select(e6_99_tree_cols, tip.tabel, five, three)
names$three <- gsub(" _.*", "", names$three)
names$five <- gsub("'", "", names$five)
write.table(names$five, "e6_sra.list", row.names = F, col.names = T, sep = "\t")

#read in
e6_sra_annot <- read.table("e6_sra_species.list")
names_e6_sra_annot <- left_join(names, e6_sra_annot, by = c("five" = "V1"))
names_e6_sra_annot[] <- t(apply(names_e6_sra_annot, 1, function(x) `length<-`(na.omit(x), length(x))))



ncolors_e6 <- length(unique(names_e6_sra_annot$three))
coolors_e6 = createPalette(ncolors_e6,  c("#ff0000", "#00ff00", "#0000ff"))
names(coolors_e6) <- NULL

key <- data.frame(unique(names_e6_sra_annot$three), coolors_e6)

names_e6_sra_annot <- left_join(names_e6_sra_annot, key, by = c("three" = "unique.names_e6_sra_annot.three."))


pwot <- ggtree(e6_99_tree) %<+% names_e6_sra_annot + geom_tippoint(aes(colour = coolors_e6.y))+ theme(legend.position = "none") + geom_tiplab(geom = "text", aes(label = three)) +
  scale_color_manual(values = names_e6_sra_annot$coolors_e6.y) + geom_tiplab(geom = "text", aes(label = five), color = "darkgrey", hjust = -3)

ggplot2::ggsave(filename = "e6_tree_99.pdf" , pwot, width = 100, height = 750, units = "cm", limitsize = F)

##############################lily dont hate me##################################
e7_99_tree <- read.tree("all_e7s.tree")
ggtree(e7_99_tree)
e7_99_tree_cols <- as.data.frame(e7_99_tree[["tip.label"]])
e7_99_tree_cols$tip.tabel <- e7_99_tree_cols$`e7_99_tree[["tip.label"]]`
e7_99_tree_cols <- separate(e7_99_tree_cols, col = `e7_99_tree[["tip.label"]]`, into = c("one", "two", "three", "four"), sep = "\\|")
e7_99_tree_cols <- separate(e7_99_tree_cols, col = one, into = c("five", "six", "seven", "eight"), sep = "_")
names <- select(e7_99_tree_cols, tip.tabel, five, three)
names$three <- gsub(" _.*", "", names$three)
names$five <- gsub("'", "", names$five)
write.table(names$five, "e7_sra.list", row.names = F, col.names = T, sep = "\t")

#read in
e7_sra_annot <- read.table("e7_sra_species.list")
names_e7_sra_annot <- left_join(names, e7_sra_annot, by = c("five" = "V1"))
names_e7_sra_annot[] <- t(apply(names_e7_sra_annot, 1, function(x) `length<-`(na.omit(x), length(x))))



ncolors_e7 <- length(unique(names_e7_sra_annot$three))
coolors_e7 = createPalette(ncolors_e7,  c("#ff0000", "#00ff00", "#0000ff"))
names(coolors_e7) <- NULL

key <- data.frame(unique(names_e7_sra_annot$three), coolors_e7)

names_e7_sra_annot <- left_join(names_e7_sra_annot, key, by = c("three" = "unique.names_e7_sra_annot.three."))


pwot <- ggtree(e7_99_tree) %<+% names_e7_sra_annot + geom_tippoint(aes(colour = coolors))+ theme(legend.position = "none") + geom_tiplab(geom = "text", aes(label = three)) +
  scale_color_manual(values = names_e7_sra_annot$coolors) + geom_tiplab(geom = "text", aes(label = five), color = "darkgrey", hjust = -3)

ggplot2::ggsave(filename = "e7_tree_99.pdf" , pwot, width = 100, height = 750, units = "cm", limitsize = F)


#additionally annotate tree with top match in blastnr
blast_dataframe <- function(df, df2){
  blast_out <- read.table(df, fill = T, sep = "\t")
  blast_out_fil <- filter(blast_out, V11 < 0.00001)
  blast_out_fil_max_e <- blast_out_fil %>% 
    group_by(V1) %>% 
    slice_min(n = 1, V11)
  blast_out_fil_max_e$tree_name <- gsub("\\|.*", "", blast_out_fil_max_e$V1)
  blast_out_fil_max_e$tree_name_2 <- gsub("_.*", "", blast_out_fil_max_e$tree_name)
  blast_out_fil_max_e$V1 <- gsub("\\/.*", "", blast_out_fil_max_e$V1)
  blast_out_fil_max_e$V1 <- gsub("-", "", blast_out_fil_max_e$V1)
  
  
  #match the names so the annotation tables can be joined
  annotation_df <- df2
  annotation_df$name_join <- sub(" .*", "", annotation_df$tip.tabel)
  annotation_df$name_join <- sub("\\'", "", annotation_df$name_join)
  annotation_df$name_join <- sub("\\/.*", "", annotation_df$name_join)
  annotation_df$name_join <- sub("\\-", "", annotation_df$name_join)
  
  
  #select relevant cols from blast output
  blast_out_sel <- select(blast_out_fil_max_e, V1, V11, V13)
  annotation_df_blast <- left_join(annotation_df, blast_out_sel, by = c("name_join" = "V1"))
  return(annotation_df_blast)
}

blast_out <- read.table("all_e6s_99_blast.tsv", fill = T, sep = "\t")

names_e6_sra_annot_blast <- blast_dataframe("all_e6s_99_blast.tsv", names_e6_sra_annot)
names_e6_sra_annot_blast <- names_e6_sra_annot_blast[!duplicated(names_e6_sra_annot_blast$tip.tabel), ]
  
pwot <- ggtree(e6_99_tree) %<+% names_e6_sra_annot_blast  + geom_tippoint(aes(colour = coolors_e6))+ theme(legend.position = "none") + geom_tiplab(geom = "text", aes(label = three)) +
  scale_color_manual(values = names_e6_sra_annot$coolors_e6) + geom_tiplab(geom = "text", aes(label = five), color = "darkgrey", hjust = -3) + geom_tiplab(geom = "text", aes(label = V13), color = "darkgrey", hjust = -4)

pwot

ggplot2::ggsave(filename = "e6_tree_99_blast.pdf" , pwot, width = 100, height = 750, units = "cm", limitsize = F)


names_e7_sra_annot_blast <- blast_dataframe("all_e7s_99_blast.tsv", names_e7_sra_annot)

names_e7_sra_annot_blast <- names_e7_sra_annot_blast[!duplicated(names_e7_sra_annot_blast$tip.tabel), ]

pwot <- ggtree(e7_99_tree) %<+% names_e7_sra_annot_blast  + geom_tippoint(aes(colour = coolors_e7))+ theme(legend.position = "none") + geom_tiplab(geom = "text", aes(label = three)) +
  scale_color_manual(values = names_e7_sra_annot$coolors_e7) + geom_tiplab(geom = "text", aes(label = five), color = "darkgrey", hjust = -3) + geom_tiplab(geom = "text", aes(label = V13), color = "darkgrey", hjust = -4)

pwot

ggplot2::ggsave(filename = "e7_tree_99_blast.pdf" , pwot, width = 100, height = 750, units = "cm", limitsize = F)

e6_tax_nam <- get_taxa_name(pwot)

write.table(e6_tax_nam, "e6_tree_tips.txt", sep = "\t", quote = F, row.names = F, col.names = F)
selected_124 <- read.table("selected.txt")
selected_124$status <- c("SYNTH")
selected_124$V1 <- paste0("'", selected_124$V1, "'")

names_e7_sra_annot_blast_synth <- left_join(names_e7_sra_annot_blast, selected_124, by = c("tip.tabel" = "V1"))


ppwot <- ggtree(e7_99_tree) %<+% names_e7_sra_annot_blast_synth  + geom_tippoint(aes(colour = coolors_e7))+ theme(legend.position = "none") + geom_tiplab(geom = "text", aes(label = label)) +
  scale_color_manual(values = names_e7_sra_annot$coolors_e7) + geom_tiplab(geom = "text", aes(label = five), color = "darkgrey", hjust = -3) + geom_tiplab(geom = "text", aes(label = V13), color = "darkgrey", hjust = -4) + geom_tiplab(geom = "text", aes(label = status), color = "red", hjust = -5, fontface = "bold", size = 5)

ggplot2::ggsave(filename = "e7_tree_99_blast_synth_test.pdf" , ppwot, width = 100, height = 750, units = "cm", limitsize = F)



##for e6
selected_e6 <- read.table("e6_sel.txt")
selected_e6$status <- c("SYNTH")
selected_e6$V1 <- paste0("'", selected_e6$V1, "'")

names_e6_sra_annot_blast_synth <- left_join(names_e6_sra_annot_blast, selected_e6, by = c("tip.tabel" = "V1"))


ppwot_e6_synth <- ggtree(e6_99_tree) %<+% names_e6_sra_annot_blast_synth  + geom_tippoint(aes(colour = coolors_e6))+ theme(legend.position = "none") + geom_tiplab(geom = "text", aes(label = label)) +
  scale_color_manual(values = names_e6_sra_annot$coolors_e6) + geom_tiplab(geom = "text", aes(label = five), color = "darkgrey", hjust = -3) + geom_tiplab(geom = "text", aes(label = V13), color = "darkgrey", hjust = -4) + geom_tiplab(geom = "text", aes(label = status), color = "red", hjust = -5, fontface = "bold", size = 5)

ggplot2::ggsave(filename = "e6_tree_99_blast_synth_test.pdf" , ppwot_e6_synth, width = 100, height = 750, units = "cm", limitsize = F)


#ughhhhhhh see if i can do a all all thing against the e6/e7 are chosen
e6_chose <- read.table("../e6_e7/e6_chosen.txt", fill = T)
e6_chose$parsed <- gsub(">", "", e6_chose$V1)
e6_chose$parsed <- gsub("__.*", "", e6_chose$parsed)

e7_chose <- read.table("../e6_e7/e7_chosen.txt", fill = T)
e7_chose$parsed <- gsub(">", "", e7_chose$V1)
e7_chose$parsed <- gsub("__.*", "", e7_chose$parsed)

e6_all_all <- read.table("E6_E6_all.aln", sep = "\t", fill = T)
colnames(e6_all_all) <- c("query", "target", "pid_E6", "score", "a", "b", "c", "d", "e", "f", "g", "h")
e7_all_all <- read.table("E7_E7_all.aln", sep = "\t", fill = T)
colnames(e7_all_all) <- c("query", "target", "pid_E7", "score", "a", "b", "c", "d", "e", "f", "g", "h")
l1_all_all <- read.table("L1_L1_all.aln", sep = "\t", fill = T)
colnames(l1_all_all) <- c("query", "target", "pid_L1", "score", "a", "b", "c", "d", "e", "f", "g", "h")

e6_all_all_chosen <- filter(e6_all_all, target %in% e6_chose$parsed)
l1_all_all_chosen_e6 <- filter(l1_all_all, query %in% e6_chose$parsed)

e7_all_all_chosen <- filter(e7_all_all, target %in% e7_chose$parsed)

l1_all_all_chosen_e7 <- filter(l1_all_all, query %in% e7_chose$parsed)

