#try to graph complexity with time to run the graph (this would all have a max cutoff of ~60 mins or so)
library(tidyverse)
library(ggplot2)
library(data.table)

source("00_utilities/hmmsearch_utils.R")


vert3 <- read.table("files/vert3_sum.txt", sep = "\t", header = F, fill = T, col.names = c(1:25))
vert3$index <- seq.int(nrow(vert3))
vert3[-1] <- t(apply(vert3[-1], 1,
                          function(x) replace(x, duplicated(x), NA)))
vert3$info <- vert3 %>% 
  rowwise() %>%
  mutate(match = str_subset(c_across(), "Saving")[1]) %>% 
  ungroup

vert3 <- as.data.frame(vert3)
vert3$info$match
vert3$info_clean <- vert3$info$match

vert3_nums <- vert3[c(2, 3, 28)]

vert3_nums$info_clean <- str_squish(vert3_nums$info_clean)

vert3_nums_info <- separate(vert3_nums, info_clean, into = c("time", "other"), sep = " (?=[^ ]+$)")
vert3_nums_info <- separate(vert3_nums_info, time, into = c("time", "misc"), sep = " ")
vert3_nums_info$posix <- strptime(vert3_nums_info$time, format='%H:%M:%S')
vert3_nums_info$secs <- vert3_nums_info$posix$hour * 3600 + vert3_nums_info$posix$min * 60 + vert3_nums_info$posix$sec
vert3_nums_info$X3 <- as.numeric(vert3_nums_info$X3)
vert3_nums_info$X2 <- as.numeric(vert3_nums_info$X2)

vert3_nums_info$vert_edges <- vert3_nums_info$X2/vert3_nums_info$X3

vert3_nums_info %>%
  ggplot(aes(x=log10(X3), y=secs)) +
  geom_point(size = 2, alpha = 0.2) +
  xlab("log10(vertices)") +
  theme_bw()

ggplot(vert3_nums_info, aes(x=X2)) +
  geom_histogram(binwidth = 100, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_density(aes(y=..count..)) +
  theme_bw()

vert3_nums_info_50k <- filter(vert3_nums_info, X3 < 50000)
vert3_nums_info_50k_list[,2] <- str_match(vert3_nums_info_50k$other, "/\\s*(.*?)\\s*_")
write.table(vert3_nums_info_50k_list[,2], "50k_under.txt", sep = "\t", col.names = F, row.names = F, quote = F)

pr_run1 <- hmmsearch_clean("all_pr_outputs.domtbl")
pr_run1$library <- word(pr_run1$query_acc, 2, sep="newfolder\\/")
pr_run1$library <- gsub("\\..*","",pr_run1$library)

pr_run1$pfam <- factor(x = pr_run1$pfam,
                                 levels = c("ncbi_and_pave_L1_I","ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH"))


pr_run1_tab <- pr_run1 %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

pr_run1_tab[-1] <- as.integer(pr_run1_tab[-1] != 0)

pr_run1_tab$sum <- rowSums(pr_run1_tab[, 2:7])
ggplot(pr_run1_tab, aes(x=sum)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

pr_run1_tab_reason <- filter(pr_run1_tab, sum >= 2)

pr_run1_tab_fulls <- filter(pr_run1_tab, sum == 6)
pr_run1_tab_fives <- filter(pr_run1_tab, sum == 5)

pr_run1_tab_fulls$library <- word(pr_run1_tab_fulls$query_acc, 2, sep="newfolder\\/")
pr_run1_tab_fulls$library <- gsub("\\..*","",pr_run1_tab_fulls$library)
length(unique(pr_run1_tab_fulls$library))
length(unique(pr_run1_tab_fulls$query_acc))


####run 2
pr_run2 <- hmmsearch_clean("hmmv3_test.domtbl")
pr_run2_b <- filter(pr_run2, pfam == "ncbi_and_novel_L1_B")
pr_run2$library <- word(pr_run2$query_acc, 2, sep="newfolder\\/")
pr_run2$library <- gsub("\\..*","",pr_run2$library)

pr_run2$pfam <- factor(x = pr_run2$pfam,
                       levels = c("ncbi_and_novel_L1_B", "ncbi_and_novel_L1_CD", "ncbi_and_novel_L1_E", "ncbi_and_novel_L1_F", "ncbi_and_novel_L1_GH", "ncbi_and_novel_L1_I"))


pr_run2_tab <- pr_run2 %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

pr_run2_tab[-1] <- as.integer(pr_run2_tab[-1] != 0)

pr_run2_tab$sum <- rowSums(pr_run2_tab[, 2:7])
ggplot(pr_run2_tab, aes(x=sum)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

pr_run2_tab_reason <- filter(pr_run2_tab, sum >= 2)

pr_run2_tab_fulls <- filter(pr_run2_tab, sum == 6)
pr_run2_tab_fives <- filter(pr_run2_tab, sum == 5)

pr_run2_tab_fulls$library <- word(pr_run2_tab_fulls$query_acc, 2, sep="newfolder\\/")
pr_run2_tab_fulls$library <- gsub("\\..*","",pr_run2_tab_fulls$library)
length(unique(pr_run2_tab_fulls$library))
length(unique(pr_run2_tab_fulls$query_acc))

B <- hmmsearch_clean("B_test_sense.dtbl")

samp <- sample(nrow(pr_run2_tab_fives),500)
samps <- pr_run2_tab_fives[samp,]
write.table(pr_run2_tab_fives$query_acc, "fives.txt", sep = "\t", col.names = F, row.names = F, quote = F)


##clustered
pr_run3 <- hmmsearch_clean("fold_srr11842090_14067_mstart_model_0.cif")
pr_run3_b <- filter(pr_run3, pfam == "ncbi_and_pave_L1_B")
pr_run3$library <- word(pr_run3$query_acc, 2, sep="newfolder\\/")
pr_run3$library <- gsub("\\..*","",pr_run3$library)

pr_run3$pfam <- factor(x = pr_run3$pfam,
                       levels = c("ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH", "ncbi_and_pave_L1_I"))


pr_run3_tab <- pr_run3 %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

pr_run3_tab[-1] <- as.integer(pr_run3_tab[-1] != 0)

pr_run3_tab$sum <- rowSums(pr_run3_tab[, 2:7])
ggplot(pr_run3_tab, aes(x=sum)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

pr_run3_tab_reason <- filter(pr_run3_tab, sum >= 2)

pr_run3_tab_fulls <- filter(pr_run3_tab, sum == 6)
pr_run3_tab_fives <- filter(pr_run3_tab, sum == 5)

pr_run3_tab_fulls$library <- word(pr_run3_tab_fulls$query_acc, 2, sep="newfolder\\/")
pr_run3_tab_fulls$library <- gsub("\\..*","",pr_run3_tab_fulls$library)
length(unique(pr_run3_tab_fulls$library))
length(unique(pr_run3_tab_fulls$query_acc))

write.table(pr_run3_tab_fives$query_acc, "fives_clust.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(pr_run3_tab_fives, "fives_clust_bookkeping.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(pr_run3_tab_fulls$query_acc, "full_clust_pr_pilot.txt", sep = "\t", col.names = F, row.names = F, quote = F)


#graph from B - CD
pr_run3_sliced <- pr_run3 %>% group_by(query_acc, pfam) %>% slice_max(order_by = score_one, n = 1)

pr_run3_sliced_CD <- filter(pr_run3_sliced, pfam == "ncbi_and_pave_L1_CD")
pr_run3_sliced_B <- filter(pr_run3_sliced, pfam == "ncbi_and_pave_L1_B")

hist(pr_run3_sliced_CD$score_one)
ggplot(pr_run3_sliced_CD, aes(x=score_one)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

ggplot(pr_run3_sliced_B, aes(x=score_one)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

fives_list <- pr_run3_tab_fives$query_acc
pr_run3_sliced_fives <- filter(pr_run3_sliced, query_acc %in% fives_list)
pr_run3_sliced_fives_sel <- select(pr_run3_sliced_fives, query_acc, pfam, score_one)


test <- pr_run3_sliced_fives_sel %>% .[!duplicated(.), ] %>%
  mutate(pfam = str_c(pfam)) %>%
  pivot_wider(names_from = pfam, values_from = c(score_one))

test[is.na(test)] <- 0

test %>%
  ggplot(aes(x=ncbi_and_pave_L1_GH, y=ncbi_and_pave_L1_B, color=query_acc)) +
  geom_point(size = 3, alpha = 0.5) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model1 vs. model2") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="none") + scale_x_continuous(expand = c(0.1, 0.1), limits = c(-5, NA)) + 
  scale_y_continuous(expand = c(0.1, 0.1), limits = c(-5, NA))

##for all
pr_run3_sliced_sel <- select(pr_run3_sliced, query_acc, pfam, score_one)


test2 <- pr_run3_sliced_sel %>% .[!duplicated(.), ] %>%
  mutate(pfam = str_c(pfam)) %>%
  pivot_wider(names_from = pfam, values_from = c(score_one))

test2[is.na(test2)] <- 0

test2 %>%
  ggplot(aes(x=ncbi_and_pave_L1_F, y=ncbi_and_pave_L1_B, color=query_acc)) +
  geom_point(size = 3, alpha = 0.5) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model1 vs. model2") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="none") + scale_x_continuous(expand = c(0.1, 0.1), limits = c(-5, NA)) + 
  scale_y_continuous(expand = c(0.1, 0.1), limits = c(-5, NA))


######clustered
pr_run5 <- hmmsearch_clean("mar25_all_pr_pilot_outputs_remaining_sort_centroids_nuc_pep.domtbl")
pr_run5_b <- filter(pr_run5, pfam == "ncbi_and_pave_L1_B")
pr_run5$library <- word(pr_run5$query_acc, 2, sep="newfolder\\/")
pr_run5$library <- gsub("\\..*","",pr_run5$library)

pr_run5$pfam <- factor(x = pr_run5$pfam,
                       levels = c("ncbi_and_pave_L1_B", "ncbi_and_pave_L1_CD", "ncbi_and_pave_L1_E", "ncbi_and_pave_L1_F", "ncbi_and_pave_L1_GH", "ncbi_and_pave_L1_I"))


pr_run5_tab <- pr_run5 %>%
  count(query_acc, pfam) %>%
  spread(pfam, n, fill = 0)

pr_run5_tab[-1] <- as.integer(pr_run5_tab[-1] != 0)

pr_run5_tab$sum <- rowSums(pr_run5_tab[, 2:7])
ggplot(pr_run5_tab, aes(x=sum)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

pr_run5_tab_reason <- filter(pr_run5_tab, sum >= 2)

pr_run5_tab_fulls <- filter(pr_run5_tab, sum == 6)
pr_run5_tab_fives <- filter(pr_run5_tab, sum == 5)

pr_run5_tab_fulls$library <- word(pr_run5_tab_fulls$query_acc, 2, sep="newfolder\\/")
pr_run5_tab_fulls$library <- gsub("\\..*","",pr_run5_tab_fulls$library)
length(unique(pr_run5_tab_fulls$library))
length(unique(pr_run5_tab_fulls$query_acc))

write.table(pr_run5_tab_fives$query_acc, "fives_clust.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(pr_run5_tab_fives, "fives_clust_bookkeping.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(pr_run5_tab_fulls$query_acc, "full_clust_pr_pilot.txt", sep = "\t", col.names = F, row.names = F, quote = F)


#graph from B - CD
pr_run5_sliced <- pr_run5 %>% group_by(query_acc, pfam) %>% slice_max(order_by = score_one, n = 1)

pr_run5_sliced_CD <- filter(pr_run5_sliced, pfam == "ncbi_and_pave_L1_CD")
pr_run5_sliced_B <- filter(pr_run5_sliced, pfam == "ncbi_and_pave_L1_B")

hist(pr_run5_sliced_CD$score_one)
ggplot(pr_run5_sliced_CD, aes(x=score_one)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

ggplot(pr_run5_sliced_B, aes(x=score_one)) +
  geom_histogram(binwidth = 1, alpha=0.5, position="dodge") +
  theme_bw()

fives_list <- pr_run5_tab_fives$query_acc
pr_run5_sliced_fives <- filter(pr_run5_sliced, query_acc %in% fives_list)
pr_run5_sliced_fives_sel <- select(pr_run5_sliced_fives, query_acc, pfam, score_one)


test <- pr_run5_sliced_fives_sel %>% .[!duplicated(.), ] %>%
  mutate(pfam = str_c(pfam)) %>%
  pivot_wider(names_from = pfam, values_from = c(score_one))

test[is.na(test)] <- 0

test %>%
  ggplot(aes(x=ncbi_and_pave_L1_GH, y=ncbi_and_pave_L1_B, color=query_acc)) +
  geom_point(size = 3, alpha = 0.5) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model1 vs. model2") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="none") + scale_x_continuous(expand = c(0.1, 0.1), limits = c(-5, NA)) + 
  scale_y_continuous(expand = c(0.1, 0.1), limits = c(-5, NA))

##for all
pr_run5_sliced_sel <- select(pr_run5_sliced, query_acc, pfam, score_one)


test2 <- pr_run5_sliced_sel %>% .[!duplicated(.), ] %>%
  mutate(pfam = str_c(pfam)) %>%
  pivot_wider(names_from = pfam, values_from = c(score_one))

test2[is.na(test2)] <- 0

test2 %>%
  ggplot(aes(x=ncbi_and_pave_L1_F, y=ncbi_and_pave_L1_B, color=query_acc)) +
  geom_point(size = 3, alpha = 0.5) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model1 vs. model2") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="none") + scale_x_continuous(expand = c(0.1, 0.1), limits = c(-5, NA)) + 
  scale_y_continuous(expand = c(0.1, 0.1), limits = c(-5, NA))

