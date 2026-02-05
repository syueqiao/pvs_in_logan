#mmtv weirdness
mmtv_runs <- read.table("SRR_Acc_List.pro", sep = "\t", fill = T)
mmtv_runs$V16 <- sub("_.*", "", mmtv_runs$V1)
mmtv_runs$e_score <- -log10(mmtv_runs$V11)
mmtv_runs_sliced <- mmtv_runs %>% group_by(V1) %>% slice_max(order_by = e_score, n = 1)

length(unique(mmtv_runs_sliced$V1))
length(unique(mmtv_runs_sliced$V16))

#misc coding
anello_taxid_full <- read.table("anello_taxid.txt", sep = ",", fill = T)
anello_taxid_ids <- as.data.frame(unique(anello_taxid$V4))
write.table(anello_taxid_ids, "anello_taxid_ids.txt", quote = F, row.names = F, col.names = F)

ggplot(anello_taxid, aes(x=V3)) +
  geom_histogram(stat = "count", fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Contig length dist.") +
  xlab("")

anello_taxid <- read.table("tax_id_parsed.txt", sep = ",", fill = T)

delete <- seq(0, nrow(anello_taxid), 2)
anello_taxid_cl <- anello_taxid[ delete ,]
clean <- as.data.frame(anello_taxid_cl$V1)
colnames(clean) <- c("col1")
b <- separate(clean, col1, into=c("id", "sk", "buh"), sep = "\\s+")
c <- as.data.frame(cbind(b$sk, b$buh))
d <- c[rowSums(is.na(c)) != ncol(c), ]

anello_taxid_miss <- read.table("missing_id.txt", sep = "\t", fill = T)

anello_taxid_miss <- as.data.frame(cbind(anello_taxid_miss$V1, anello_taxid_miss$V6))
e <- left_join(d, anello_taxid_miss, by = "V1")
f <- as.data.frame(cbind(e$V1, coalesce(e$V2.x, e$V2.y)))
colnames(f) <- c("V4", "sk")
f$V4 <- as.numeric(f$V4)
g <- left_join(anello_taxid_full, f, by = "V4")
g <- g %>% 
  mutate_all( ~ str_replace_all(., " ", "_"))
write.table(g, "g", quote = F, row.names = F, col.names = F)
g2 <- read.table("g", sep = "\t")
g2 <- g2 %>% 
  mutate_all( ~ str_replace_all(., "Bacteriateria", "Bacteria"))

ggplot(g2, aes(x=reorder(V7, V7, length))) +
  geom_bar(stat = "count", fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("distribution of anello hits to library superkingdom") +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5))


anello_taxid_full_plasmo <- subset(anello_taxid_full, grepl("Plasmodium ", V6))

                                   