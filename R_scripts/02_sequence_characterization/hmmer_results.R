library(tidyverse)
library(rentrez)
library(ggplot2)
library(viridis)
library(ggpubr)

#read in the table file, and then do summary statistics compared to the input accession list - for every "main genome (e.g., untranslated)" one columns plus 2nd column that describes the frame translation
comp_acc <- read.table("pv_accession.acc", sep = "",  header = F, fill = T)
herp_acc <- read.table("herp_acc.acc", sep = "",  header = F, fill = T)
poly_acc <- read.table("poly_acc.acc", sep = "", header = F, fill = T)
SRR_acc <- read.table("SRR_clean_acc", sep = "", header = F, fill = T)

##testing with tbl-out output
PVfam2_scan <- read.table("PVfam2_scan_tbl.out", sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:18]
colnames(PVfam2_scan) <- c("pfam", "pfam_acc", "query_acc", "misc", "eval_full", "score_full", "bias_full", "eval_one", "score_one", "bias_one", "exp_dom", "reg_dom", "clu_dom", "ov_dom", "env_dom", "rep_dom", "inc_dom")
PVfam2_scan$query_acc_clean <- str_split(PVfam2_scan$query_acc, "\\_", simplify=T)[,1]
in_comp_not_hit_PVfam2_scan <- data.frame(setdiff(comp_acc$V1 ,PVfam2_scan$query_acc_clean))

###Do for dom.tblout, has mor information on where the hits are etc
##Pfam full includes all of the PV-specific families on Pfam, minus less specific superfamilies
#function to make summstats table

##############################################TO CREATE WORKING DF FROM HMMER OUTPUT######################################################
hmmer_clean <- function(input_domtbl){
  PVfam_full_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:22]
  colnames(PVfam_full_scan) <- c("pfam", "pfam_acc", "tlen", "query_acc", "misc", "qlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                 "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  PVfam_full_scan$query_acc_clean <- sub("_[^_]+$", "", PVfam_full_scan$query_acc)
  
  in_comp_not_hit_PVfam_full <- data.frame(setdiff(comp_acc$V1 ,PVfam_full_scan$query_acc_clean))
  
  #about 80% hit rate for all nucleotide sequences on ncbi
  
  # write.table(in_comp_not_hit_PVfam_full, "in_comp_not_hit_PVfam_full.csv", sep=",",  col.names=FALSE, row.names = F)
  # #all ones without hit are long control regions, eg, those without protein sequences associated with them
  # 
  # #smallest qlen that gave a good eval is 10 amino acids, with 10-4 eval (decent?)
  # PVfam_full_scan_list <- split(PVfam_full_scan, PVfam_full_scan$query_acc_clean)
  # # PVfam_full_scan_ceval_min <- lapply(PVfam_full_scan_list, function(x) slice_min(x, c_eval))
  # # names(PVfam_full_scan_ceval_min) <- NULL
  # df_PVfam_full_scan_ceval_min <- do.call(rbind.data.frame, PVfam_full_scan_ceval_min)
  
  return(PVfam_full_scan)
}
# length(unique(df_PVfam_full_scan_ceval_min$query_acc_clean))
#there is one best domain for each input sequence

##############################################SUMM STATS######################################################

summ_hmm_tbl <- function(hmm_table, acc){
  in_comp_not_hit <- data.frame(setdiff(acc$V1, hmm_table$query_acc_clean))
  #about 80% hit rate for all nucleotide sequences on ncbi
  per_hits <- 100 - (length(unique(in_comp_not_hit$setdiff.acc.V1..hmm_table.query_acc_clean.))/length(acc$V1)*100)
  total_input_seq <- length(acc$V1)
  total_nohit_seq <- length(unique(in_comp_not_hit$setdiff.acc.V1..hmm_table.query_acc_clean.))
  mean_ceval <- mean(hmm_table$c_eval)
  mean_eval <- mean(hmm_table$eval_full)
  max_ceval <- max(hmm_table$c_eval)
  min_ceval <- min(hmm_table$c_eval)
  max_eval <- max(hmm_table$eval_full)
  min_eval <- min(hmm_table$eval_full)
  #10^-7 average e val
  #10^-4 tio 10^-3 max eval
  
  append_func <- function(dfNames) {
    do.call(rbind, lapply(dfNames, function(x) {
      cbind(get(x), source = x)
    }))
  }
  
  summary_table <- data.frame(append_func(c("per_hits", "mean_ceval", "mean_eval", "max_ceval", "min_ceval", "max_eval", "min_eval", 
                                            "total_input_seq", "total_nohit_seq")))
  return(summary_table)
}

####functions
df_PVfam_full_scan <- hmmer_clean_ceval_min("PVfam_full_scan_domtbl.out")
PVfam_full_scam_summ <- summ_hmm_tbl(df_PVfam_full_scan_ceval_min, comp_acc)
# df_PVfam_full_scan_list <- split(df_PVfam_full_scan, df_PVfam_full_scan$query_acc_clean)



##############################################DF TO COMPARE QLEN + HIT COUNT######################################################

###compare qlen to hit counts and plot
qlen_v_hit_fun <- function(df){
  df_Genfam_full_scan_table <-  data.frame(table(df$query_acc_clean))
  colnames(df_Genfam_full_scan_table) <- c("query_acc_clean", "hit_count")
  df_Genfam_full_scan_qlen <- select(df, query_acc_clean, qlen)
  
  ##working df
  df_Genfam_full_scan_table_qlen <- left_join(df_Genfam_full_scan_table, df_Genfam_full_scan_qlen, "query_acc_clean" )
  df_Genfam_full_scan_table_qlen$ratio <- df_Genfam_full_scan_table_qlen$qlen/df_Genfam_full_scan_table_qlen$hit_count
  
  
  # df_Genfam_full_scan_table_qlen_lm <- lm(formula = hit_count ~ qlen, data = df_Genfam_full_scan_table_qlen)
  # Gen_lm_formula <- paste0("y ~ ", round(coefficients(df_Genfam_full_scan_table_qlen_lm)[1],2), " + ", 
  #          paste(sprintf("%.2f * %s", 
  #                        coefficients(df_Genfam_full_scan_table_qlen_lm)[-1],  
  #                        names(coefficients(df_Genfam_full_scan_table_qlen_lm)[-1])), 
  #                collapse=" + ")
  return(df_Genfam_full_scan_table_qlen)
}
df_PVfam_full_scan_table_qlen_test <- qlen_v_hit_fun(df_PVfam_full_scan)

##############################################DF TO COMPARE QLEN + ##UNIQUE## HIT COUNT AND PLOT######################################################

###unique hits
qlen_v_hit_fun_uniq <- function(df){
  df <- df[!duplicated(df[c("pfam","query_acc_clean")]),]
  df_Genfam_full_scan_table <-  data.frame(table(df$query_acc_clean))
  colnames(df_Genfam_full_scan_table) <- c("query_acc_clean", "hit_count")
  df_Genfam_full_scan_qlen <- select(df, query_acc_clean, qlen)
  
  ##working df
  df_Genfam_full_scan_table_qlen <- left_join(df_Genfam_full_scan_table, df_Genfam_full_scan_qlen, "query_acc_clean" )
  df_Genfam_full_scan_table_qlen$ratio <- df_Genfam_full_scan_table_qlen$qlen/df_Genfam_full_scan_table_qlen$hit_count
  
  
  # df_Genfam_full_scan_table_qlen_lm <- lm(formula = hit_count ~ qlen, data = df_Genfam_full_scan_table_qlen)
  # Gen_lm_formula <- paste0("y ~ ", round(coefficients(df_Genfam_full_scan_table_qlen_lm)[1],2), " + ", 
  #          paste(sprintf("%.2f * %s", 
  #                        coefficients(df_Genfam_full_scan_table_qlen_lm)[-1],  
  #                        names(coefficients(df_Genfam_full_scan_table_qlen_lm)[-1])), 
  #                collapse=" + ")
  return(df_Genfam_full_scan_table_qlen)
}

df_PVfam_full_scan_table_qlen_uniq <- qlen_v_hit_fun_uniq(df_PVfam_full_scan)

plot_fun_formula <- function(df){
  sub_df_lm <- lm(formula = hit_count ~ qlen, data = df)
  PV_lm_formula <-paste0("y ~ ", round(coefficients(sub_df_lm)[1],2), " + ", 
                         paste(sprintf("%.2f * %s", 
                                       coefficients(sub_df_lm)[-1],  
                                       names(coefficients(sub_df_lm)[-1])), 
                               collapse=" + ")
  )
  p <- ggplot(df, 
              aes(x = qlen, 
                  y = hit_count)) + 
    geom_point() +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())+
    labs(color=NULL)+
    # scale_y_log10()+
    # scale_x_log10()+
    geom_smooth(method="lm", size =0.5) +
    stat_cor(method = "pearson") +
    ggtitle(deparse(substitute(df))) +
    theme_bw() +
    theme(axis.title=element_text(size=6)) +
    geom_text(data=data.frame(), aes(label = paste("formula = ", PV_lm_formula), x = 100, y = 0),
              hjust = 0, vjust = 0, color = "black")
  pp <- p + scale_colour_viridis()
  return(pp)
}

plot_fun_formula(df_PVfam_full_scan_table_qlen)  
plot_fun_formula(df_PVfam_full_scan_table_qlen_uniq)


##############################################LOOK FOR MODELS THAT HAVE FRAGMENTS (and plot)######################################################

fragment_models <- function(df){
  sub_df <- data.frame(cbind(df$query_acc_clean, df$pfam, df$query_acc))
  sub_df_counts <- data.frame(table(query_acc_clean = sub_df$X1, sub_df$X2))
  sub_df_counts <- filter(sub_df_counts, Freq > 0)
  colnames(sub_df_counts) <- c("query_acc_clean", "model", "model_hits")
  return(sub_df_counts)
}

df_PVfam_full_scan_frags <- fragment_models(df_PVfam_full_scan)

plot_frag <- function(df){
  p <- df %>%
    mutate(text = fct_reorder(model, model_hits)) %>%
    ggplot( aes(x=model_hits, color=model, fill=model)) +
    geom_histogram(alpha=0.6, binwidth = 1) +
    theme_bw() +
    scale_x_continuous(breaks = c(0:12)) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    # theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    ggtitle(deparse(substitute(df)))
  # xlab("") +
  # ylab("Assigned Probability (%)") +
  facet_wrap(~text)
  return(p)
}

plot_frag(df_PVfam_full_scan_frags)

####plot the full sequence eval versus the number of (uniq?) model hits
eval_v_hits_uniq <- function(df_og, df_frags){
  df_og_sub <- select(df_og, query_acc_clean, query_acc, eval_full, pfam)
  colnames(df_og_sub) <- c("query_acc_clean", "query_acc", "eval_full", "model")
  df_frags_eval <- left_join(df_frags, df_og_sub, by = c("query_acc_clean", "model"))
  df_frags_eval <- df_frags_eval[!duplicated(df_frags_eval[c(4,5)]),]
  df_frags_eval$eval_full_log <- -log10(df_frags_eval$eval_full)
  return(df_frags_eval)
}

df_PVfam_full_scan_eval_v_hits <- eval_v_hits_uniq(df_PVfam_full_scan, df_PVfam_full_scan_frags)

plot_eval_v_hits <- function(df){
  p <- df %>%
    # mutate(text = fct_reorder(model, model_hits)) %>%
    ggplot(aes(x=model_hits, y=eval_full_log, group = model_hits, fill = model)) +
    geom_boxplot(outlier.size = 0.5, alpha=0.6) +
    facet_wrap(~model) +
    theme_bw() +
    # scale_fill_viridis(discrete = T) +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)) +
    scale_fill_viridis(discrete = T) +
    ggtitle(deparse(substitute(df))) +
    # ) +
    #   scale_alpha_continuous(
    #     range = c(0.3, 0.9),
    #     guide = guide_legend(override.aes = list(fill = "black"))) +
    #   scale_fill_manual(values = c("goldenrod2", "forestgreen", "dodgerblue4"))
    # # # # xlab("") +
    ylab("-log10(eval)")
  
  p
  return(p)
}


plot_eval_v_hits(df_PVfam_full_scan_eval_v_hits)

######eval vs qlen per model#####
plot_eval_v_qlen <- function(df){
  df$eval_full_log <- -log10(df$eval_full)
  p <- df %>%
    ggplot(aes(x = qlen, y = eval_full_log, color = pfam))+
    geom_point(size = 1) + 
    # geom_smooth(method = "lm", se = FALSE) + 
    facet_wrap(pfam~.) +
    theme_bw() +
    scale_color_viridis(discrete = T) +
    # scale_color_viridis(discrete=TRUE) +
    theme(
      legend.position="none") +
  ggtitle(deparse(substitute(df))) +
  ylab("-log10(eval)")
return(p)
}

plot_eval_v_qlen(df_PVfam_full_scan)
##########################comparative statistics for herps SAME AS ABOVE, just compressed##################
df_herp_PVfam_full_scan <- hmmer_clean("herp_PVfam_full_scan_domtbl.out")
summ_hmm_tbl(df_herp_PVfam_full_scan, herp_acc)
df_herp_PVfam_full_scan_table_qlen_test <- qlen_v_hit_fun(df_herp_PVfam_full_scan)
df_herp_PVfam_full_scan_table_qlen_uniq <- qlen_v_hit_fun_uniq(df_herp_PVfam_full_scan)

plot_fun_formula(df_herp_PVfam_full_scan_table_qlen_test) 
plot_fun_formula(df_herp_PVfam_full_scan_table_qlen_uniq)  


df_herp_PVfam_full_scan_frags <- fragment_models(df_herp_PVfam_full_scan)
plot_frag(df_herp_PVfam_full_scan_frags)

df_herp_PVfam_full_scan_eval_v_hits <- eval_v_hits_uniq(df_herp_PVfam_full_scan, df_herp_PVfam_full_scan_frags)
plot_eval_v_hits(df_herp_PVfam_full_scan_eval_v_hits)
plot_eval_v_qlen(df_herp_PVfam_full_scan)




write.csv(df_PVfam_full_scan_ceval_min, "df_PVfam_full_scan_ceval_min.csv")
write.csv(df_herp_PVfam_full_scan_ceval_min, "df_herp_PVfam_full_scan_ceval_min.csv")


#########look for hit rate over the complete nucleotide sequences of herpes and pap#########
herp_full_acc <- read.table("complete_nucleotide_acc_herps.acc", sep = "",  header = F, fill = T)
colnames(herp_full_acc) <- c("query_acc_clean")
herp_full_acc <- filter(herp_full_acc, query_acc_clean %in% df_herp_PVfam_full_scan$query_acc_clean)
df_herp_PVfam_full_scan_comp_nuc <- filter(df_herp_PVfam_full_scan, query_acc_clean %in% herp_full_acc$query_acc_clean)
herp_full_acc_nh <- data.frame(setdiff(herp_full_acc$query_acc_clean, df_herp_PVfam_full_scan$query_acc_clean))
100 - 100*length(herp_full_acc_nh)/length(herp_full_acc$query_acc_clean)

df_PVfam_full_scan_comp_nuc <- select(df_PVfam_full_scan_comp_nuc, query_acc_clean, pfam)
comp_nuc_table <- data.frame(table(df_PVfam_full_scan_comp_nuc))
comp_nuc_spr <- spread(comp_nuc_table, key = pfam, value = Freq)

comp_nuc_spr_cov <- data.frame(colSums(comp_nuc_spr != 0))
colnames(comp_nuc_spr_cov) <- c("hits")
comp_nuc_spr_cov$per_cov <- 100*comp_nuc_spr_cov$hits/4358
PV_comp_nuc_spr_cov <- comp_nuc_spr_cov


animal_pave_genomes <- read.table("animal_pave_genomes.txt", sep = ",", header = T, fill = T)
animal_pave_genomes_tab <- data.frame(table(animal_pave_genomes$Host.Common.Name))
animal_pave_genomes_non <- filter(animal_pave_genomes_tab, !grepl('Domestic Cow|Domestic dog', Var1))
sum(animal_pave_genomes_non$Freq)

animal_pave_genomes %>%
ggplot(aes(x=Host.Common.Name, fill = Host.Common.Name)) +
  geom_histogram(stat = 'count') +
  theme_bw() +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  theme(axis.text=element_text(size=8)) +
  # theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  ggtitle(deparse(substitute(df)))

ggplot(animal_pave_genomes_tab, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="none") +
  scale_fill_viridis(discrete=TRUE)
  # 
  # geom_text(aes(y = ypos, label = Var1), color = "white", size=6) +
  # scale_fill_brewer(palette="Set1")

##########checking description of the ones with no hit - eutils is sooooooooooooooo slow and this is a list of 9k
in_comp_not_hit_PVfam_full_desc <- NULL;  
for (acc in in_comp_not_hit_PVfam_full$setdiff.comp_acc.V1..PVfam_full_scan.query_acc_clean.){
  r_search <- entrez_search(db="nucleotide", term = acc)
  entrez_db_summary("cdd")
  asumm <- entrez_summary(db = "nucleotide", id = r_search$ids)
  asumm <- asumm$title
  in_comp_not_hit_PVfam_full_desc <- rbind(in_comp_not_hit_PVfam_full_desc, asumm)
}

##################3
#testing for polyomaviruses
df_poly_PVfam_full_scan <- hmmer_clean("poly_PVfam_full_scan_domtbl.out")
summ_hmm_tbl(df_poly_PVfam_full_scan, poly_acc)
df_poly_PVfam_full_scan_table_qlen_test <- qlen_v_hit_fun(df_poly_PVfam_full_scan)
df_poly_PVfam_full_scan_table_qlen_uniq <- qlen_v_hit_fun_uniq(df_poly_PVfam_full_scan)

plot_fun_formula(df_poly_PVfam_full_scan_table_qlen_test) 
plot_fun_formula(df_poly_PVfam_full_scan_table_qlen_uniq)  


df_poly_PVfam_full_scan_frags <- fragment_models(df_poly_PVfam_full_scan)
plot_frag(df_poly_PVfam_full_scan_frags)

df_poly_PVfam_full_scan_eval_v_hits <- eval_v_hits_uniq(df_poly_PVfam_full_scan, df_poly_PVfam_full_scan_frags)
plot_eval_v_hits(df_poly_PVfam_full_scan_eval_v_hits)
plot_eval_v_qlen(df_poly_PVfam_full_scan)

sum(duplicated(df_poly_PVfam_full_scan_frags$query_acc_clean))
multi_hit_poly <- subset(df_poly_PVfam_full_scan_frags,duplicated(query_acc_clean))

##########SRR?
df_SRRtest_PVfam_full_scan <- hmmer_clean("SRR_test_PVfam_full_scan_domtbl.out")
summ_hmm_tbl(df_SRRtest_PVfam_full_scan, SRR_acc)

#plot of sequences that hit and sequences that don't hit 
#no long seuqneces that arent being hitm bascially
#asymptotic, vairable number length vs model hits
#each domain, missing model and fragmented models, size distribution of everything, what's there etc, 

#######rewrote to accomodate hmmsearch instead of hmmscan, as it seems likely that this will be what is used in later processing, and not hmmscan
hmmsearch_clean <- function(input_domtbl){
  PVfam_full_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:22]
  PVfam_full_scan <- PVfam_full_scan[!(is.na(PVfam_full_scan$V21) | PVfam_full_scan$V21=="" | !(PVfam_full_scan$V2 == "-")), ]

  colnames(PVfam_full_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                 "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  PVfam_full_scan$query_acc_clean <- sub("_[^_]+$", "", PVfam_full_scan$query_acc)
  PVfam_full_scan <- type.convert(PVfam_full_scan, as.is = TRUE)
return(PVfam_full_scan)
}

poly_search <- hmmsearch_clean("buhhh.out")
summ_hmm_tbl(poly_search, poly_acc)

poly_search_qlen <- qlen_v_hit_fun(poly_search)
poly_search_qlen_uniq <- qlen_v_hit_fun_uniq(poly_search)

plot_fun_formula(poly_search_qlen) 
plot_fun_formula(poly_search_qlen_uniq)
#####################3
pap_search <- hmmsearch_clean("bitscore_pap_test.out")
summ_hmm_tbl(pap_search, comp_acc)

pap_search_qlen <- qlen_v_hit_fun(pap_search)
pap_search_qlen_uniq <- qlen_v_hit_fun_uniq(pap_search)

plot_fun_formula(pap_search_qlen) 
plot_fun_formula(pap_search_qlen_uniq)


############select random subset of SRA runs
set.seed(69)
SRA_runs <- read.table("chordata_nomushomo.csv", sep = ",", header = T)
SRA_runs_1k <- sample_n(SRA_runs, 1000)
write.csv(SRA_runs_1k, "SRA_runs_1ksample.csv", row.names = F)

#############filter for ORFs that are hit from HMMER via bitscore for later input into DIAMOND
df_PVfam_full_scan_bitscorefil <- subset(df_PVfam_full_scan, score_one > 30)
df_PVfam_full_scan_bitscorefil <- as.data.frame(cbind(df_PVfam_full_scan_bitscorefil$query_acc, df_PVfam_full_scan_bitscorefil$envcoord_from, df_PVfam_full_scan_bitscorefil$envcoord_to))
colnames(df_PVfam_full_scan_bitscorefil) <- c("fasta", "start", "end")
write.table(df_PVfam_full_scan_bitscorefil, "df_PVfam_full_scan_bitscorefil.txt", sep = "\t", quote = F, row.names = F)

######parsing headers#####
df_PVfam_full_scan_filtered <- subset(df_PVfam_full_scan, score_one > 30)
df_PVfam_full_scan_filtered <- as.data.frame(cbind(df_PVfam_full_scan_filtered$query_acc_clean, df_PVfam_full_scan_filtered$pfam_acc,
                                                   df_PVfam_full_scan_filtered$pfam))
colnames(df_PVfam_full_scan_filtered) <- c("pfam_name", "pfam_acc", "fasta")
write.table(df_PVfam_full_scan_filtered, "df_PVfam_full_scan_filtered.txt", sep = "\t", quote = F, row.names = F)

##no more fricking with unix
all_pv_accessions <- read.table("pv_accessions_sort.txt", sep = "",  header = F, fill = T)
all_pfam_hits <- read.table("acc_pfam_pacc_sort.txt", sep = " ",  header = F, fill = T)
key <- left_join(all_pfam_hits, all_pv_accessions, by = "V1")
key$V2.y <- gsub('-', '_', key$V2.y)
key <- cbind(key$V1, key$V2.y)
write.table(key, "accession_key.txt", sep = " ", quote = F, row.names = F)
####ughhghgh

species_info <- read.csv("sequences_species.csv", sep = ",", header = T, fill = T)
species_info <- species_info
all_pfam_hits_u <- as.data.frame(unique(all_pfam_hits$V1))
colnames(all_pfam_hits_u) <- c("Accession")
all_pfam_hits_u$Accession <- gsub('>', '', all_pfam_hits_u$Accession)

key <- left_join(all_pfam_hits_u, species_info, by = "Accession")

####new results with orfs
df_ORF_scan <- hmmsearch_clean("ORFs_bugless.tbl")
df_ORF_scan <- df_ORF_scan %>% 
  group_by(query_acc) %>% 
  slice_min(i_eval, n = 1) %>%
  ungroup()

df_ORF_scan<- df_ORF_scan[!duplicated(df_ORF_scan$query_acc), ]

df_ORF_scan <- as.data.frame(cbind(df_ORF_scan$query_acc, df_ORF_scan$query_acc_clean,
                     df_ORF_scan$pfam))
colnames(df_ORF_scan) <- c("Accession_ORF", "Accession", "Family")
key$Accession <- gsub('NC_', 'NC+', key$Accession, fixed = T)
df_ORF_scan_species <- left_join(df_ORF_scan, key, by = "Accession")
df_ORF_scan_species <- df_ORF_scan_species[1:4]
df_ORF_scan_species$Species <- gsub(' ', '_', df_ORF_scan_species$Species, fixed = T)


df_ORF_scan_species<- df_ORF_scan_species[!duplicated(df_ORF_scan_species), ]

write.table(df_ORF_scan_species, "ORF_key_bugless.txt", sep = " ", quote = F, row.names = F)



df_anello_scan <- hmmsearch_clean("anello_hmmer.out")
df_anello_scan <- df_anello_scan %>% 
  group_by(query_acc) %>% 
  slice_max(score_full, n = 1) %>%
  ungroup()
df_anello_out <- as.data.frame(cbind(df_anello_scan$query_acc, df_anello_scan$pfam))
write.table(df_anello_out, "anello_map.txt", quote = F, row.names = F)

clus <- read.csv("anello_clus.csv")

ggplot(clus, aes(y=per, x=clusters)) + 
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_y_continuous(breaks = round(seq(min(clus$per), max(clus$per), by = 5),1)) +
  scale_x_continuous(breaks = round(seq(min(clus$clusters), max(clus$clusters), by = 100),1))

L1_scan <- hmmsearch_clean("../all_large_contigs_hmmer.domtblout")
L1_only <- filter(L1_scan, pfam == "Late_protein_L1")
L1 <- L1_only$query_acc
write.table(L1, "L1_contigs.txt", quote = F, row.names = F)


L1_clus <- read.csv("L1_clus.csv")

ggplot(L1_clus, aes(y=per, x=clusters)) + 
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_y_continuous(breaks = round(seq(min(clus$per), max(clus$per), by = 5),1)) +
  scale_x_continuous(breaks = round(seq(min(clus$clusters), max(clus$clusters), by = 100),1))

##new PV
pv_2 <- hmmsearch_clean("updated_PV_nucs.domtbl")
accessory <- c("Papilloma_E5", "E6", "E7", "Pap_E4")

#filter on full evalue, since we are looking at only a single orf this time. Doesn't really matter
#I think I should also filter on 
pv_2_ieval_fil <- filter(pv_2, eval_full < 0.00001)
pv_2_ieval_fil <- filter(pv_2_ieval_fil, !pfam %in% accessory)
pv_2_ieval_fil_slice <- pv_2_ieval_fil %>% 
  group_by(query_acc) %>% 
  slice_max(score_one, n = 1) %>%
  ungroup()
write.table(pv_2_ieval_fil_slice$query_acc, sep = ",", "pv_2_hmmer_hits.txt", row.names = F, col.names = F, quote = F)
pv_2_ieval_fil_key <- select(pv_2_ieval_fil_slice, query_acc, pfam)
write.table(pv_2_ieval_fil_key, sep = " ", "hmmer_key_pv2.txt", row.names = F, col.names = F, quote = F)
