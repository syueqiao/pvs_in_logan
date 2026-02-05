library(tidyverse)
library(gggenes)
human_sticker <- read.table("hs1.pro", sep = "\t", fill = T)
human_sticker$e_score <- -log10(human_sticker$V11)

human_sticker_sliced <- human_sticker %>% group_by(V1) %>% slice_max(order_by = e_score, n = 1)
human_sticker_sliced <- human_sticker_sliced %>% group_by(V1) %>% slice_max(order_by = V10, n = 1)


human_sticker_sliced$q_cov <- abs(human_sticker_sliced$V2 - human_sticker_sliced$V3)/human_sticker_sliced$V4
human_sticker_sliced$m_cov <- abs(human_sticker_sliced$V7 - human_sticker_sliced$V8)/human_sticker_sliced$V9
human_sticker_sliced$m_cov_aa <- abs(human_sticker_sliced$V7 - human_sticker_sliced$V8)

ggplot(human_sticker_sliced, aes(x=V4)) +
  geom_histogram(, binwidth = 50, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("Contig length dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(human_sticker_sliced$V4), by = 500),1))

ggplot(human_sticker_sliced, aes(x=e_score)) +
  geom_histogram(, binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("E-score dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(human_sticker_sliced$e_score), by = 20),1))

ggplot(human_sticker_sliced, aes(x=V10)) +
  geom_histogram(, binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  ggtitle("P-iden dist.") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 5),1))



mean(human_sticker_sliced$V4)


ggplot(human_sticker_sliced, aes(x=m_cov, y=e_score, colour = V10)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  scale_color_viridis_b(option = "D")


human_sticker_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov, FUN = median), y=m_cov, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

human_sticker_sliced %>%
  ggplot(aes(x=reorder(V6, m_cov_aa, FUN = median), y=m_cov_aa, color=log10(e_score))) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_viridis_c() +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("model vs. coverage in aas")  +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

samp <- sample(nrow(human_sticker_sliced),10)
samps <- human_sticker_sliced[samp,]


human_sticker_sliced <- human_sticker_sliced %>%
  group_by(V6) %>%
  mutate(genes_min = min(c(V7, V8))) %>%
  mutate(V7 = V7 - genes_min) %>%
  mutate(V8 = V8 - genes_min) %>%
  ungroup()

human_sticker_fil <- filter(human_sticker, V11 < 0.00001)

p <- ggplot(human_sticker_fil, aes(xmin = V7, xmax = V8, y = V6, 
                     fill = "red", label = V6)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), 
                  colour = "blue", size = 0, alpha = 0.5) +
  geom_feature( data = human_sticker_fil,
                      aes(x = V9, y = V6, forward = F)) +
  geom_feature_label( data = human_sticker_fil,
    aes(x = V9, y = V6, label = c("end"), forward = F)) +
  scale_x_continuous(labels = scales::comma, limits = c(0, 450)) + 
  theme_genes() %+replace% theme(axis.text = element_text(size = 10))


ggsave("location.pdf", width = 50, height = 25, units = "cm", limitsize = F)

### remove
df <- human_sticker_fil %>% 
  group_by(V6) %>% 
  summarise(
    minimum = min(V7),
    maximum = max(V8), 
    length = max(V9)
  )

write.table(df, "human_hits_bl.csv", sep = ",", quote = F, col.names = F, row.names = F)
