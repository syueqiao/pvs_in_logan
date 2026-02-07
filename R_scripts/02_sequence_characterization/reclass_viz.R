all_pave <- read.table("files/all_pv.txt", sep = "\t", header = F)
e4 <- read.table("files/e4_pv.txt", sep = "\t", header = F)
e5 <- read.table("files/e5_pv.txt", sep = "\t", header = F)
e6 <- read.table("files/e6_pv.txt", sep = "\t", header = F)
e7 <- read.table("files/e7_pv.txt", sep = "\t", header = F)
e8 <- read.table("files/e8_pv.txt", sep = "\t", header = F)
e9 <- read.table("files/e9_pv.txt", sep = "\t", header = F)
e10 <- read.table("files/e10_pv.txt", sep = "\t", header = F)

all_pave$e4<- with(all_pave, ifelse(V1 %in% e4$V1, "1", "0"))
all_pave$e5<- with(all_pave, ifelse(V1 %in% e5$V1, "1", "0"))
all_pave$e6<- with(all_pave, ifelse(V1 %in% e6$V1, "1", "0"))
all_pave$e7<- with(all_pave, ifelse(V1 %in% e7$V1, "1", "0"))
all_pave$e8<- with(all_pave, ifelse(V1 %in% e8$V1, "1", "0"))
all_pave$e9<- with(all_pave, ifelse(V1 %in% e9$V1, "1", "0"))
all_pave$e10<- with(all_pave, ifelse(V1 %in% e10$V1, "1", "0"))

all_pave[,2:8] <- sapply(all_pave[,2:8],as.numeric)

all_pave$sum <- rowSums(all_pave[, 2:8])
all_pave <- separate(all_pave, V1, into = c("name"), sep = "\\|")

ggplot(all_pave, aes(x=sum)) +
  geom_histogram(binwidth = 0.5, alpha=0.5, fill="#D86ECC", aes(y = after_stat(count/sum(count)))) +
  geom_text(aes(label = paste0(after_stat(signif(count*100/sum(count), digits=2)), "%"), y= ..prop..), stat= "count", vjust = -.5) +
  theme_bw() +
  ggtitle("Distribution number of accessory genes in PV genomes") +
  xlab("number of accessory genes") +
  ylab("percentage of genomes") +
  scale_y_continuous(labels = scales::percent, breaks = round(seq(0, 100, by = 0.1),1))

gene_counts <- as.data.frame(colSums(all_pave[, 2:8]))
gene_counts <-tibble::rownames_to_column(gene_counts, "gene")
colnames(gene_counts) <- c("gene", "numba")
  
ggplot(gene_counts, aes(x=gene,  y=after_stat(gene_counts$numba*100/853))) +
  geom_col(fill="#8ED973", alpha=0.5, aes(y = after_stat(gene_counts$numba*100/853))) +
  geom_text(aes(label = paste0(signif(100 * after_stat(gene_counts$numba/853), 2), "%")), vjust = -0.25) +
  theme_bw() +
  ggtitle("Distribution number of accessory genes in PV genomes") +
  xlab("number of accessory genes") +
  ylab("percentage of genomes")  +
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1))



feb7_abundances <- feb_7_pv_lib_type[,c(1, 6, 4, 25, 26)]
feb7_abundances <- filter(feb7_abundances, V2 %in% c("TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "VIRAL RNA"))

feb7_abundances <- feb7_abundances %>% 
  separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
feb7_abundances$cov <- gsub("_.*","",feb7_abundances$cov)
feb7_abundances$cov <- as.numeric(feb7_abundances$cov)

feb7_abundances$cov_norm <- feb7_abundances$cov/feb7_abundances$qlen

early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
late <- c("Late_protein_L1", "Late_protein_L2")

feb7_abundances$stage <- with(feb7_abundances, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
feb7_abundances <- subset(feb7_abundances, !grepl("remove", stage))
feb7_abundances_ratio <- feb7_abundances[,c(6, 3, 9)]
feb7_abundances_ratio$cov <- as.numeric(feb7_abundances_ratio$cov)

feb7_abundances_ratio_agg <- aggregate(cov ~ library + stage, data=feb7_abundances_ratio, FUN=sum)
feb7_abundances_ratio_agg_res <- reshape(feb7_abundances_ratio_agg, idvar = "library", timevar = "stage", direction = "wide")
feb7_abundances_ratio_agg_res[is.na(feb7_abundances_ratio_agg_res)] <- 0

feb7_abundances_ratio_agg_res$le_ratio <- feb7_abundances_ratio_agg_res$cov.late/feb7_abundances_ratio_agg_res$cov.early
hist(log10(feb7_abundances_ratio_agg_res$le_ratio))

feb7_abundances_ratio_agg_res <- feb7_abundances_ratio_agg_res %>% 
  filter_at(vars(le_ratio), all_vars(!is.infinite(.)))

feb7_abundances_ratio_agg_res$le_ratio_var <- log10(feb7_abundances_ratio_agg_res$le_ratio + .Machine$double.eps) 


ggplot(feb7_abundances_ratio_agg_res, aes(x=le_ratio_var)) +
  geom_histogram( binwidth=0.5, fill="#47D45A", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Distribution of late/early gene abundance in RNA-seq libraries") +
  xlab("Log10(late/early abundance)") +
  ylab("Library count")

feb7_abundances_ratio_agg_res$le_ratio_var_nz <- log10(feb7_abundances_ratio_agg_res$le_ratio) 

ggplot(feb7_abundances_ratio_agg_res, aes(x=le_ratio_var_nz)) +
  geom_histogram( binwidth=0.1, fill="#13501B", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Distribution of late/early gene abundance in RNA-seq libraries") +
  xlab("Log10(late/early abundance)") +
  ylab("Library count")


feb7_abundances_ratio_agg_res_qq <- filter(feb7_abundances_ratio_agg_res, !le_ratio == Inf)
qqnorm(feb7_abundances_ratio_agg_res_qq$le_ratio, pch = 1, frame = FALSE)

feb7_abundances_ratio_agg_res <- subset(feb7_abundances_ratio_agg_res, !grepl("0", le_ratio))



ncbi_counts <- hmmsearch_clean("files/ncbi_input.domtbl")
ncbi_counts <- ncbi_counts[,c(1, 4)]
as.data.frame(table(ncbi_counts$pfam))


feb_7_pv_lib_type_RNA <- filter(feb_7_pv_lib_type, V2 %in% c("TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "VIRAL RNA", "METATRANSCRIPTOMIC"))
feb_7_pv_lib_type_RNA$stage <- with(feb_7_pv_lib_type_RNA, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
feb_7_pv_lib_type_RNA <- subset(feb_7_pv_lib_type_RNA, !grepl("remove", stage))

rna <- as.data.frame(table(feb_7_pv_lib_type_RNA$stage))
colnames(rna) <- c("stage", "rna")

feb_7_pv_lib_type_DNA <- filter(feb_7_pv_lib_type, V2 %in% c("GENOMIC", "GENOMIC SINGLE CELL", "METAGENOMIC"))
feb_7_pv_lib_type_DNA$stage <- with(feb_7_pv_lib_type_DNA, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
feb_7_pv_lib_type_DNA <- subset(feb_7_pv_lib_type_DNA, !grepl("remove", stage))

dna <- as.data.frame(table(feb_7_pv_lib_type_DNA$stage))
colnames(dna) <- c("stage", "dna")
chi <- cbind(rna, dna)
chi2 <- chi[,-1]
rownames(chi2) <- chi[,1]

chi2 <- chi2[,c(1,3)]

chisq.test(chi2)

length(unique(feb_7_pv_lib_type_RNA$library))

###
feb7_abundances_DNA <- filter(feb7_abundances, !V2 %in% c("TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "VIRAL RNA"))

feb7_abundances_DNA <- feb7_abundances_DNA %>% 
  separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
feb7_abundances_DNA$cov <- gsub("_.*","",feb7_abundances_DNA$cov)
feb7_abundances_DNA$cov <- as.numeric(feb7_abundances_DNA$cov)

feb7_abundances_DNA$cov_norm <- feb7_abundances_DNA$cov/feb7_abundances_DNA$tlen

early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
late <- c("Late_protein_L1", "Late_protein_L2")

feb7_abundances_DNA$stage <- with(feb7_abundances_DNA, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
feb7_abundances_DNA$gene_length <- with(feb7_abundances_DNA, ifelse(pfam == "PPV_E1_N", "650", 
                                                                    ifelse(pfam == "PPV_E1_DBD", "650",
                                                                        ifelse(pfam == "PPV_E1_C", "640", 
                                                                               ifelse(pfam == "PPV_E2_N", "365", 
                                                                                      ifelse(pfam == "PPV_E2_C", "365", 
                                                                                             ifelse(pfam == "E6", "158", 
                                                                                                    ifelse(pfam == "E7", "98", 
                                                                                                           ifelse(pfam == "Late_protein_L1", "568 ", 
                                                                                                                  ifelse(pfam == "Late_protein_L2", "473", "remove"))))))))))
feb7_abundances_DNA$cov_norm_gene <- feb7_abundances_DNA$cov/as.numeric(feb7_abundances_DNA$gene_length)

feb7_abundances_DNA <- subset(feb7_abundances_DNA, !grepl("remove", stage))
feb7_abundances_DNA_ratio <- feb7_abundances_DNA[,c(6, 3, 8, 9, 11)]

feb7_abundances_DNA_ratio_agg <- aggregate(cov ~ library + stage, data=feb7_abundances_DNA_ratio, FUN=sum)
feb7_abundances_DNA_ratio_agg_res <- reshape(feb7_abundances_DNA_ratio_agg, idvar = "library", timevar = "stage", direction = "wide")
feb7_abundances_DNA_ratio_agg_res[is.na(feb7_abundances_DNA_ratio_agg_res)] <- 0

feb7_abundances_DNA_ratio_agg_res$le_ratio <- feb7_abundances_DNA_ratio_agg_res$cov.late/feb7_abundances_DNA_ratio_agg_res$cov.early
hist(log10(feb7_abundances_DNA_ratio_agg_res$le_ratio))

feb7_abundances_DNA_ratio_agg_res <- feb7_abundances_DNA_ratio_agg_res %>% 
  filter_at(vars(le_ratio), all_vars(!is.infinite(.)))

feb7_abundances_DNA_ratio_agg_res$le_ratio_var <- log10(feb7_abundances_DNA_ratio_agg_res$le_ratio + .Machine$double.eps) 


ggplot(feb7_abundances_DNA_ratio_agg_res, aes(x=le_ratio_var)) +
  geom_histogram( binwidth=0.5, fill="#47D45A", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Distribution of late/early gene abundance in DNA-seq libraries") +
  xlab("Log10(late/early abundance)") +
  ylab("Library count")

feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz <- log10(feb7_abundances_DNA_ratio_agg_res$le_ratio) 

ggplot(feb7_abundances_DNA_ratio_agg_res, aes(x=le_ratio_var_nz)) +
  geom_histogram( binwidth=0.1, fill="#80350E", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Distribution of late/early gene abundance in DNA-seq libraries") +
  xlab("Log10(late/early abundance)") +
  ylab("Library count")

qqnorm(feb7_abundances_DNA_ratio_agg_res$le_ratio, pch = 1, frame = FALSE)


feb7_abundances_DNA_ratio_agg_res_zeros <- filter(feb7_abundances_DNA_ratio_agg_res, le_ratio == Inf)
median(feb7_abundances_DNA_ratio_agg_res_zeros$cov.late)
hist(feb7_abundances_DNA_ratio_agg_res_zeros$cov.early)

feb7_abundances_RNA_ratio_agg_res_zeros <- filter(feb7_abundances_ratio_agg_res, le_ratio == !Inf|!0)
median(feb7_abundances_ratio_agg_res$cov.early)
hist(feb7_abundances_RNA_ratio_agg_res_zeros$cov.early)

ggplot(feb7_abundances_DNA_ratio_agg_res_zeros, aes(x=le_ratio)) +
  geom_histogram( binwidth=0.005, fill="#80350E", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Distribution of late/early gene abundance in DNA-seq libraries") +
  stat_summary() +
  xlab("Log10(late/early abundance)") +
  ylab("Library count")

ggplot(feb7_abundances_DNA_ratio_agg_res, aes(x=log10(cov.early), y=log10(cov.late))) + 
  geom_point() +
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="blue")

###
library_type_ratios <- function(df, list, color1, color2){
  feb7_abundances_DNA <- filter(df, V2 %in% list)
  
  feb7_abundances_DNA <- feb7_abundances_DNA %>% 
    separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
  feb7_abundances_DNA$cov <- gsub("_.*","",feb7_abundances_DNA$cov)
  feb7_abundances_DNA$cov <- as.numeric(feb7_abundances_DNA$cov)
  
  feb7_abundances_DNA$cov_norm <- feb7_abundances_DNA$cov/feb7_abundances_DNA$tlen
  
  early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
  late <- c("Late_protein_L1", "Late_protein_L2")
  
  feb7_abundances_DNA$stage <- with(feb7_abundances_DNA, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
  feb7_abundances_DNA$gene_length <- with(feb7_abundances_DNA, ifelse(pfam == "PPV_E1_N", "650", 
                                                                      ifelse(pfam == "PPV_E1_DBD", "650",
                                                                             ifelse(pfam == "PPV_E1_C", "640", 
                                                                                    ifelse(pfam == "PPV_E2_N", "365", 
                                                                                           ifelse(pfam == "PPV_E2_C", "365", 
                                                                                                  ifelse(pfam == "E6", "158", 
                                                                                                         ifelse(pfam == "E7", "98", 
                                                                                                                ifelse(pfam == "Late_protein_L1", "568 ", 
                                                                                                                       ifelse(pfam == "Late_protein_L2", "473", "remove"))))))))))
  feb7_abundances_DNA$cov_norm_gene <- feb7_abundances_DNA$cov/as.numeric(feb7_abundances_DNA$gene_length)
  
  feb7_abundances_DNA <- subset(feb7_abundances_DNA, !grepl("remove", stage))
  feb7_abundances_DNA_ratio <- feb7_abundances_DNA[,c(6, 3, 8, 9, 11)]
  
  feb7_abundances_DNA_ratio_agg <- aggregate(cov ~ library + stage, data=feb7_abundances_DNA_ratio, FUN=sum)
  feb7_abundances_DNA_ratio_agg_res <- reshape(feb7_abundances_DNA_ratio_agg, idvar = "library", timevar = "stage", direction = "wide")
  feb7_abundances_DNA_ratio_agg_res[is.na(feb7_abundances_DNA_ratio_agg_res)] <- 0
  
  feb7_abundances_DNA_ratio_agg_res$le_ratio <- feb7_abundances_DNA_ratio_agg_res$cov.late/feb7_abundances_DNA_ratio_agg_res$cov.early
  
  feb7_abundances_DNA_ratio_agg_res <- feb7_abundances_DNA_ratio_agg_res %>% 
    filter_at(vars(le_ratio), all_vars(!is.infinite(.)))
  
  feb7_abundances_DNA_ratio_agg_res$le_ratio_var <- log10(feb7_abundances_DNA_ratio_agg_res$le_ratio + .Machine$double.eps) 
  
  
  p_full <- ggplot(feb7_abundances_DNA_ratio_agg_res, aes(x=le_ratio_var)) +
    geom_histogram(aes(y = ..density..), binwidth=0.5, fill=color1, color=color2, alpha=0.9) +
    theme_bw() +
    ggtitle(paste0("Distribution of late/early gene abundance\n", list)) +
    xlab("Log10(late/early abundance)") +
    geom_vline(aes(xintercept = median(feb7_abundances_DNA_ratio_agg_res$le_ratio_var[is.finite(feb7_abundances_DNA_ratio_agg_res$le_ratio_var)])), colour = "black", linetype = "dashed", size = 1, alpha = 0.5) +
    ylab("Library count")
  
  feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz <- log10(feb7_abundances_DNA_ratio_agg_res$le_ratio) 
  
  p_zoom <- ggplot(feb7_abundances_DNA_ratio_agg_res, aes(x=le_ratio_var_nz)) +
    geom_histogram(aes(y = ..density..), binwidth=0.05, fill=color1, color=color2, alpha=0.4) +
    geom_density(lwd = 1, colour = color1, alpha = 0.8) + 
    xlim(-5, 5) +
    theme_bw(base_size = 15) +
    ggtitle(paste0("Distribution of late/early gene abundance (inset)\n", list)) +
    xlab("Log10(late/early abundance)") +
    ylab("Library count") 
  
  p_full_annot <- annotation_custom(grid::textGrob(label = 
                                                     paste0("Median = ", 
                                                            signif(median(feb7_abundances_DNA_ratio_agg_res$le_ratio_var[is.finite(feb7_abundances_DNA_ratio_agg_res$le_ratio_var)]), 3), 
                                                            "\n", 
                                                            "Mean = ", signif(mean(feb7_abundances_DNA_ratio_agg_res$le_ratio_var[is.finite(feb7_abundances_DNA_ratio_agg_res$le_ratio_var)]), 3),
                                                            "\n", 
                                                            "SD = ", signif(sd(feb7_abundances_DNA_ratio_agg_res$le_ratio_var[is.finite(feb7_abundances_DNA_ratio_agg_res$le_ratio_var)]), 3)),
                                                  x = unit(0.75, "npc"), y = unit(0.8, "npc"),
                                                  gp = grid::gpar(cex = 1)))
  
  p_zoom_annot <- annotation_custom(grid::textGrob(label = 
                                                     paste0("Median = ", 
                                                            signif(median(feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz[is.finite(feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz)]), 3), 
                                                            "\n", 
                                                            "Mean = ", signif(mean(feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz[is.finite(feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz)]), 3), 
                                                            "\n", 
                                                            "SD = ", signif(sd(feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz[is.finite(feb7_abundances_DNA_ratio_agg_res$le_ratio_var_nz)]), 3)),
                                                   x = unit(0.75, "npc"), y = unit(0.8, "npc"),
                                                   gp = grid::gpar(cex = 1)))
  
  
  p_full_zoom_print <- p_zoom + p_zoom_annot
  
  ggplot2::ggsave(filename = paste0("outputs/plot_zoom",list,".png"),p_full_zoom_print, width = 10, height = 10, units = "cm")

}


lib_types <- c("TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "GENOMIC", "METAGENOMIC", "GENOMIC SINGLE CELL")

colors_green_brown <- c("#47D45A", "#C2F1C8", "#296030", "#C55B25", "#F9CCB5", "#80350E")

aaa <- c("GENOMIC SINGLE CELL")

library_type_ratios(feb7_abundances, aaa, "#C55B25", "#808080")

for (j in 1:length(lib_types)) {
    x = paste0(lib_types[j])
    y = paste0(colors_green_brown[j])
    library_type_ratios(feb7_abundances, x, y, "#808080")
  }

library_type_ratios(feb7_abundances, trans, "#C2F1C8", "#808080")
library_type_ratios(feb7_abundances, gen, "#80350E", "#C0C0C0")

length(feb7_abundances$library[grepl("^METATRANSCRIPTOMIC$", feb7_abundances$V2)])

risks <- read.table("files/risks.txt", sep = "\t", header = T)
risks <- risks[-21,]

risks$risk <- factor(c("high", "med", "low"))

risks$risk <- gsub("med", "low", risks$risk)

risks %>%
  mutate(risk = fct_relevel(risk, 
                            "high", "med", "low")) %>%
  ggplot( aes(x=risk, y=median, fill = risk)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set2") +
  geom_jitter(color="black", size=2, alpha=0.5, width = 0.1) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")


#test cov

trans_all <- c("METATRANSCRIPTOMIC", "TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL")

gen_all <- c("GENOMIC")

get_cov <- function(df, list, value){
feb7_abundances_trans <- filter(df, V2 %in% list)

feb7_abundances_trans <- feb7_abundances_trans %>% 
  separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
feb7_abundances_trans$cov <- gsub("_.*","",feb7_abundances_trans$cov)
feb7_abundances_trans$cov <- as.numeric(feb7_abundances_trans$cov)

early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
late <- c("Late_protein_L1", "Late_protein_L2")

feb7_abundances_trans$stage <- with(feb7_abundances_trans, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))


feb7_abundances_trans <- subset(feb7_abundances_trans, !grepl("remove", stage))
feb7_abundances_trans_ratio <- feb7_abundances_trans[,c(6, 3, 8)]

feb7_abundances_trans_ratio_agg <- aggregate(cov ~ library + stage, data=feb7_abundances_trans_ratio, FUN=sum)
feb7_abundances_trans_ratio_agg_res <- reshape(feb7_abundances_trans_ratio_agg, idvar = "library", timevar = "stage", direction = "wide")
feb7_abundances_trans_ratio_agg_res[is.na(feb7_abundances_trans_ratio_agg_res)] <- 0

feb7_abundances_trans_ratio_agg_res$le_ratio <- feb7_abundances_trans_ratio_agg_res$cov.late/feb7_abundances_trans_ratio_agg_res$cov.early

feb7_abundances_DNA_ratio_agg_res_zeros <- filter(feb7_abundances_trans_ratio_agg_res, le_ratio == value)
r_late <- median(feb7_abundances_DNA_ratio_agg_res_zeros$cov.late)
r_early <- median(feb7_abundances_DNA_ratio_agg_res_zeros$cov.early)

ee <- c(r_late, r_early)
return(ee)

}

trans <- c("TRANSCRIPTOMIC")
gen <- c("GENOMIC")
get_cov(feb7_abundances, trans, 0)

get_cov_non <- function(df, list, value){
  feb7_abundances_trans <- filter(df, V2 %in% list)
  
  feb7_abundances_trans <- feb7_abundances_trans %>% 
    separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
  feb7_abundances_trans$cov <- gsub("_.*","",feb7_abundances_trans$cov)
  feb7_abundances_trans$cov <- as.numeric(feb7_abundances_trans$cov)
  
  early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
  late <- c("Late_protein_L1", "Late_protein_L2")
  
  feb7_abundances_trans$stage <- with(feb7_abundances_trans, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
  
  
  feb7_abundances_trans <- subset(feb7_abundances_trans, !grepl("remove", stage))
  feb7_abundances_trans_ratio <- feb7_abundances_trans[,c(6, 3, 8)]
  
  feb7_abundances_trans_ratio_agg <- aggregate(cov ~ library + stage, data=feb7_abundances_trans_ratio, FUN=sum)
  feb7_abundances_trans_ratio_agg_res <- reshape(feb7_abundances_trans_ratio_agg, idvar = "library", timevar = "stage", direction = "wide")
  feb7_abundances_trans_ratio_agg_res[is.na(feb7_abundances_trans_ratio_agg_res)] <- 0

  
  feb7_abundances_trans_ratio_agg_res$le_ratio <- feb7_abundances_trans_ratio_agg_res$cov.late/feb7_abundances_trans_ratio_agg_res$cov.early
  feb7_abundances_trans_ratio_agg_res <- filter(feb7_abundances_trans_ratio_agg_res, le_ratio > 0)
  
  feb7_abundances_trans_ratio_agg_res <- feb7_abundances_trans_ratio_agg_res %>% 
    filter_at(vars(le_ratio), all_vars(!is.infinite(.)))
  r_late <- median(feb7_abundances_trans_ratio_agg_res$cov.late)
  r_early <- median(feb7_abundances_trans_ratio_agg_res$cov.early)
  
  ef <- c(r_late, r_early)
  return(ef)
}
gen <- c("TRANSCRIPTOMIC")
get_cov_non(feb7_abundances, gen)


feb7_abundances_RNA_ratio_agg_res_zeros <- filter(feb7_abundances_ratio_agg_res, le_ratio == Inf)
median(feb7_abundances_ratio_agg_res$cov.early)
hist(feb7_abundances_RNA_ratio_agg_res_zeros$cov.early)

library(ggupset)

lib_types <- c("TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL", "GENOMIC", "METAGENOMIC", "GENOMIC SINGLE CELL")

get_upset <- function(df, strategies){
  
  df_fil <- filter(df, V2 == strategies)
  df_fil <- df_fil %>% 
    separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
  df_fil$cov <- gsub("_.*","",df_fil$cov)
  df_fil$cov <- as.numeric(df_fil$cov)
  
  early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
  late <- c("Late_protein_L1", "Late_protein_L2")
  
  df_fil$stage <- with(df_fil, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
  
  
  df_fil <- subset(df_fil, !grepl("remove", stage))
  df_fil_ratio <- df_fil[,c(6, 3, 8)]
  
  df_fil_ratio_agg_mean <- aggregate(cov ~ library + stage, data=df_fil_ratio, FUN=mean)
  df_fil_ratio_agg_sum <- aggregate(cov ~ library + stage, data=df_fil_ratio, FUN=sum)
  
  df_fil_ratio_agg_mean <- reshape(df_fil_ratio_agg_mean, idvar = "library", timevar = "stage", direction = "wide")
  df_fil_ratio_agg_mean[is.na(df_fil_ratio_agg_mean)] <- 0
  
  df_fil_ratio_agg_mean$le_ratio <- df_fil_ratio_agg_mean$cov.late/df_fil_ratio_agg_mean$cov.early
  
  df_fil_ratio_agg_mean$late_status <- with(df_fil_ratio_agg_mean, ifelse(cov.late > 0, TRUE, FALSE))
  df_fil_ratio_agg_mean$early_status <- with(df_fil_ratio_agg_mean, ifelse(cov.early > 0, TRUE, FALSE))
  
  df_fil_ratio_agg_mean <- filter(df_fil_ratio_agg_mean, !library == "SRR15013819")
  
  combinations_df_fil_ratio_agg_res <- df_fil_ratio_agg_mean |>  
    mutate(
      combination = pmap(
        list(late_status, early_status),
        \(lgl1, lgl2) {
          c('Late', 'Early')[c(lgl1, lgl2)] 
        }
      )
    )
  
d <- combinations_df_fil_ratio_agg_res|> within(x <- sapply(combination, paste, collapse = "_"))
   
p <-  d %>% ggplot(aes(x = combination)) +
  geom_bar(stat='count', fill = c("#FF5050", "#8f5cb0", "#3366FF"), alpha = 0.8) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.7) +
  scale_x_upset() +
  axis_combmatrix(sep = "_", levels = c("Late", "Early")) +
  theme_bw(base_size = 15) +
  ggtitle(paste0("Upset plot of late/early gene presence\n", strategies)) +
  xlab("Combination") +
  ylab("Library count")

 ggplot2::ggsave(filename = paste0("outputs/plot_",strategies,".png"),p, width = 20, height = 20, units = "cm")
 return(p)

}

for (x in lib_types) {
  get_upset(feb7_abundances, x)
  
}

get_scatter <- function(df, aaaa){
  
  df_fil <- filter(df, V2 %in% aaaa)
  df_fil <- df_fil %>% 
    separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
  df_fil$cov <- gsub("_.*","",df_fil$cov)
  df_fil$cov <- as.numeric(df_fil$cov)
  
  early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
  late <- c("Late_protein_L1", "Late_protein_L2")
  
  df_fil$stage <- with(df_fil, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
  
  
  df_fil <- subset(df_fil, !grepl("remove", stage))
  df_fil_ratio <- df_fil[,c(6, 3, 8, 7)]
  
  df_fil_ratio_agg_mean <- aggregate(cov ~ library + stage + V2, data=df_fil_ratio, FUN=sum)
  
  df_fil_ratio_agg_mean <- reshape(df_fil_ratio_agg_mean, idvar = c("library", "V2"), timevar = "stage", direction = "wide")
  df_fil_ratio_agg_mean[is.na(df_fil_ratio_agg_mean)] <- 0
  
  df_fil_ratio_agg_mean <- filter(df_fil_ratio_agg_mean, !library == "SRR15013819")
  
  df_fil_ratio_agg_mean$late_status <- with(df_fil_ratio_agg_mean, ifelse(cov.late > 0, TRUE, FALSE))
  df_fil_ratio_agg_mean$early_status <- with(df_fil_ratio_agg_mean, ifelse(cov.early > 0, TRUE, FALSE))
  
  combinations_df_fil_ratio_agg_res <- df_fil_ratio_agg_mean |>  
    mutate(
      combination = pmap(
        list(late_status, early_status),
        \(lgl1, lgl2) {
          c('Late', 'Early')[c(lgl1, lgl2)] 
        }
      )
    )
  
  d <- combinations_df_fil_ratio_agg_res|> within(x <- sapply(combination, paste, collapse = "_"))
  
  
peep <-  d %>% ggplot(aes(x = cov.early, y = cov.late, color = as.character(x))) + 
  geom_point(size=3, alpha = 0.4)  +
  geom_smooth(data=d[d$x == "Late_Early",], method = "lm", fill = NA, linewidth = 0.5) +
  theme_bw(base_size = 20) +
  ggtitle(paste0("Scatter plot of late/early gene coverage")) +
  xlab("ka.f of early genes") +
  ylab("ka.f of late genes") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#FF5050", "#3366FF","#8f5cb0")) + expand_limits(x = 0, y = 0)

ggplot2::ggsave(filename = paste0("plot_scatter",aaaa,".png"),peep, width = 20, height = 20, units = "cm")
}
aaaaaaa <- c("GENOMIC")

df<-feb7_abundances
library <- c("GENOMIC")

for (x in lib_types) {
  get_scatter(feb7_abundances, x)
  
}

get_scatter2 <- function(df, strategies){
 
  #goddamn this step takes a long time 
  df_fil <- filter(df, V2 %in% lib_types)
  df_fil <- df_fil %>% 
    separate(query_acc,c("acc", "cov"), sep = "ka_f_", remove = FALSE)
  df_fil$cov <- gsub("_.*","",df_fil$cov)
  df_fil$cov <- as.numeric(df_fil$cov)
  
  early <- c("PPV_E1_N", "PPV_E1_DBD", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "E6", "E7")
  late <- c("Late_protein_L1", "Late_protein_L2")
  
  df_fil$stage <- with(df_fil, ifelse(pfam %in% early, "early", ifelse(pfam %in% late, "late", "remove")))
  
  
  df_fil <- subset(df_fil, !grepl("remove", stage))
  df_fil_ratio <- df_fil[,c(6, 3, 8, 7)]
  
  df_fil_ratio_agg_mean <- aggregate(cov ~ library + stage + V2, data=df_fil_ratio, FUN=sum)

  df_fil_ratio_agg_mean <- reshape(df_fil_ratio_agg_mean, idvar = c("library", "V2"), timevar = "stage", direction = "wide")
  df_fil_ratio_agg_mean[is.na(df_fil_ratio_agg_mean)] <- 0
  
  df_fil_ratio_agg_mean$le_ratio <- df_fil_ratio_agg_mean$cov.late/df_fil_ratio_agg_mean$cov.early
  df_fil_ratio_agg_mean <- filter(df_fil_ratio_agg_mean, !library == "SRR15013819")
  
  
  df_fil_ratio_agg_mean %>% ggplot(aes(x = log10(cov.early), y = log10(cov.late), color = V2)) +
    geom_point(size=0.5, alpha = 0.4)  +
    geom_smooth(method = "lm", fill = NA, linewidth = 0.5) +
    theme_bw(base_size = 15) +
    ggtitle(paste0("Scatter plot of late/early gene coverage")) +
    xlab("ka.f of early genes") +
    ylab("ka.f of late genes") +
    scale_color_brewer(palette = "Set3") + expand_limits(x = 0, y = 0)

  ggplot2::ggsave(filename = paste0("outputs/plot_",strategies,".png"),p, width = 20, height = 20, units = "cm")
  return(p)
  
}
