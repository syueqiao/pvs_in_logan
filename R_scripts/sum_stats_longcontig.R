#PV summary stats
library(tidyverse)
library(rentrez)
library(ggplot2)
library(viridis)
library(ggpubr)
library(reshape2)

#parse HMMER domtbl output
hmmsearch_clean <- function(input_domtbl){
  PVfam_full_scan <- read.table(input_domtbl, sep = "",  header = F, fill = T) %>% na.omit() %>% .[,1:22]
  PVfam_full_scan <- PVfam_full_scan[!(is.na(PVfam_full_scan$V21) | PVfam_full_scan$V21=="" | !(PVfam_full_scan$V2 == "-")), ]
  
  colnames(PVfam_full_scan) <- c("query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen", "eval_full", "score_full", "bias_full", "#", "of", 
                                 "c_eval", "i_eval", "score_one", "bias_one", "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to", "envcoord_from", "envcoord_to", "acc")
  PVfam_full_scan$query_acc_clean <- sub("_[^_]+$", "", PVfam_full_scan$query_acc)
  PVfam_full_scan <- type.convert(PVfam_full_scan, as.is = TRUE)
  return(PVfam_full_scan)
}
# 
# workflow_test <- hmmsearch_clean("large_contigs_hmmer")
# workflow_test$library_acc <- sub("_[^_]+$", "", workflow_test$query_acc_clean)
# workflow_test <- subset(workflow_test, eval_full < 0.0000000001)

#input is the 90% iden clustered ORFs, after running thru HMMER
cent_test <- hmmsearch_clean("centroid_hits.domtblout")
cent_test$library_acc <- sub("_[^_]+$", "", cent_test$query_acc_clean)
cent_test$pfam_clean <- gsub("^.*E1.*$", "E1", gsub("^.*E2.*$", "E2", cent_test$pfam))
#clean for easier parsing
cent_test <- subset(cent_test, eval_full < 0.0000000001)

#extract relevant information to see P/A of genes
cent_orfs <- as.data.frame(cbind(cent_test$query_acc_clean, cent_test$pfam))

cent_orfs_longs <- cent_orfs %>%
  mutate(yesno = 1) %>%
  distinct %>%
  spread(V2, yesno, fill = 0)

cent_orfs_longs$E1_sum <- rowSums(cent_orfs_longs[, c("PPV_E1_DBD", "PPV_E1_N", "PPV_E1_C")])
cent_orfs_longs$E2_sum <- rowSums(cent_orfs_longs[, c("PPV_E2_N", "PPV_E2_C")])
centroid_PA <- cent_orfs_longs

m <- as.matrix(centroid_PA[, -1]) # -1 to omit categories from matrix
clust <- hclust(dist(t(m)))

centroid_PA_ns <- centroid_PA[,1:12]

cent_mat <- centroid_PA %>% remove_rownames %>% column_to_rownames(var="V1")
cent_scaled <- as.matrix(scale(cent_mat))
pheatmap::pheatmap(cent_mat, cellheight=5)

cent_mat_ns <- centroid_PA_ns %>% remove_rownames %>% column_to_rownames(var="V1")
cent_scaled_ns <- as.matrix(scale(cent_mat_ns))
ns_heatmap <- pheatmap::pheatmap(cent_mat_ns, cellheight=3, cellwidth=3, fontsize = 5)

ggplot(melt(centroid_PA), aes(V1, variable, alpha = value)) + 
  geom_tile(colour = "gray50") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(limits = colnames(m)[clust$order])

##quantifying how many have all core family hits
core <- c("PPV_E1_DBD", "PPV_E1_N", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "Late_protein_L1", "Late_protein_L2")
centroid_PA_core <- select(centroid_PA, V1, all_of(core))
centroid_PA_core$sum <- rowSums(centroid_PA_core[,2:8])
#eg, 7 means we have em all etc...
ggplot(centroid_PA_core, aes(x=sum, fill = sum)) + 
  geom_histogram() +
  theme_bw() +
  scale_fill_viridis(discrete = T)

#ORF length information
cent_length <- as.data.frame(cbind(cent_test$query_acc, cent_test$qlen, cent_test$pfam_clean, cent_test$pfam))
cent_length <- cent_length[!duplicated(cent_length), ]
cent_length$V2 <- as.numeric(cent_length$V2)

ggplot(cent_length, aes(x=V3, y=V2, fill = V4)) + 
  geom_boxplot()+
  theme(legend.position="none") +
  theme_bw() +
  scale_fill_viridis(discrete = T)

ggplot(cent_length, aes(x=V4, y=V2, fill = V4)) + 
  geom_boxplot()+
  theme(legend.position="none") +
  theme_bw() +
  scale_fill_viridis(discrete = T)
#size x, bitscore on y
##genome/contig length information
#combine the diamond hits with these hits based on contig ID with join?
diamond_hits <- large_clean_fil
diamond_hits <- as.data.frame(cbind(large_clean_fil$V1, large_clean_fil$V4))
colnames(diamond_hits) <- c("query_acc_clean", "contig_length")
diamond_hits <- distinct(diamond_hits)
genome_length <- merge(x = cent_test, y = diamond_hits, by = "query_acc_clean", all.x = TRUE)
genome_length$contig_length <- as.numeric(genome_length$contig_length)

ggplot(genome_length, aes(x=contig_length)) + 
  geom_histogram(bins = 50)

#phmmer, align L1 orf to profile

#draw out for ORFs, what are the cases
#every L1, align to profile of PV L1, MSA (gappy?) phmmer
#x axis clsutering, y = # of clusters
#known reference L1 at 90% -> place "novel" on known set
#see which are novel? maybe all human PVs
#MSA of 25k L1, IQTREE, stat about clustering and level of novelty



#########################redo for the big table, hope it dont crash
curl_test <- hmmsearch_clean("curl_large_contigs_hmmer.domtblout")
curl_test$library_acc <- sub("_[^_]+$", "", curl_test$query_acc_clean)
curl_test$pfam_clean <- gsub("^.*E1.*$", "E1", gsub("^.*E2.*$", "E2", curl_test$pfam))
#clean for easier parsing
curl_test <- subset(curl_test, eval_full < 0.0000000001)

#extract relevant information to see P/A of genes
curl_orfs <- as.data.frame(cbind(curl_test$query_acc_clean, curl_test$pfam))

curl_orfs_longs <- curl_orfs %>%
  mutate(yesno = 1) %>%
  distinct %>%
  spread(V2, yesno, fill = 0)

curl_orfs_longs$E1_sum <- rowSums(curl_orfs_longs[, c("PPV_E1_DBD", "PPV_E1_N", "PPV_E1_C")])
curl_orfs_longs$E2_sum <- rowSums(curl_orfs_longs[, c("PPV_E2_N", "PPV_E2_C")])
curlroid_PA <- curl_orfs_longs

m <- as.matrix(curlroid_PA[, -1]) # -1 to omit categories from matrix
clust <- hclust(dist(t(m)))

curlroid_PA_ns <- curlroid_PA[,1:12]

curl_mat <- curlroid_PA %>% remove_rownames %>% column_to_rownames(var="V1")
curl_scaled <- as.matrix(scale(curl_mat))
pheatmap::pheatmap(curl_mat, cellheight=5)

curl_mat_ns <- curlroid_PA_ns %>% remove_rownames %>% column_to_rownames(var="V1")
curl_scaled_ns <- as.matrix(scale(curl_mat_ns))
ns_heatmap <- pheatmap::pheatmap(curl_mat_ns, cellheight=1, cellwidth=1, fontsize = 5)

ggplot(melt(curlroid_PA), aes(V1, variable, alpha = value)) + 
  geom_tile(colour = "gray50") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(limits = colnames(m)[clust$order])

##quantifying how many have all core family hits
core <- c("PPV_E1_DBD", "PPV_E1_N", "PPV_E1_C", "PPV_E2_N", "PPV_E2_C", "Late_protein_L1", "Late_protein_L2")
curlroid_PA_core <- select(curlroid_PA, V1, all_of(core))
curlroid_PA_core$sum <- rowSums(curlroid_PA_core[,2:8])
sum(curlroid_PA_core$Late_protein_L1)
#eg, 7 means we have em all etc...
ggplot(curlroid_PA_core, aes(x=sum, fill = sum)) + 
  geom_histogram() +
  theme_bw() +
  scale_fill_viridis(discrete = T)

#ORF length information
curl_length <- as.data.frame(cbind(curl_test$query_acc, curl_test$qlen, curl_test$pfam_clean, curl_test$pfam))
curl_length <- curl_length[!duplicated(curl_length), ]
curl_length$V2 <- as.numeric(curl_length$V2)

ggplot(curl_length, aes(x=V3, y=V2, fill = V4)) + 
  geom_violin()+
  theme(legend.position="none") +
  theme_bw() +
  scale_fill_viridis(discrete = T)

ggplot(curl_length, aes(x=V4, y=V2, fill = V4)) + 
  geom_boxplot()+
  theme(legend.position="none") +
  theme_bw() +
  scale_fill_viridis(discrete = T)

