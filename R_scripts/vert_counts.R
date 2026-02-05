library(tidyverse)
gbnstall.packages("ggpmisc")
library(ggpmisc)

pv_counts <- read.table("asfasfasfasfasfaasfasfasf.csv", header = T, sep = ",")

ggplot(pv_counts, aes(x=logmb, y=logpvs)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point() +
  geom_smooth(method=lm, level=0.95, color="#3c5412", fill = "#b3de69") +
  theme_bw()


uniq_counts <- read.table("uniq.counts", sep = "\t")
vert_mbs <- read.table("vert_counts_fr.txt", sep = "\t")

#convert both to 2 names
uniq_counts <- separate(data = uniq_counts, col = V2, into = c("first", "second", "third"), sep = " ")
uniq_counts$name_2 <- paste0(uniq_counts$first, " ", uniq_counts$second)
uniq_counts <- select(uniq_counts, name_2, V1)
uniq_counts_agg <- aggregate(V1 ~ name_2, FUN = sum, data=uniq_counts)


vert_mbs <- separate(data = vert_mbs, col = V1, into = c("first", "second", "third"), sep = " ")
vert_mbs$name_2 <- paste0(vert_mbs$first, " ", vert_mbs$second)
vert_mbs <- select(vert_mbs, name_2, V2)
vert_mbs_agg <- aggregate(V2 ~ name_2, FUN = sum, data=vert_mbs)



vert_mbs_pvs <- left_join(vert_mbs_agg, uniq_counts_agg, by = "name_2")
vert_mbs_pvs[is.na(vert_mbs_pvs)] <- 0

vert_mbs_pvs$ratio <- vert_mbs_pvs$V1/(vert_mbs_pvs$V2/1000)
vert_mbs_pvs$v2_v1 <- vert_mbs_pvs$V2/vert_mbs_pvs$V1
vert_mbs_pvs$V1_v2 <- vert_mbs_pvs$V1/vert_mbs_pvs$V2


vert_mbs_pvs_fil <- filter(vert_mbs_pvs, ratio > 0)

target <- c("Homo sapiens", "Accipiter nisus")

                            
  p <- ggplot(vert_mbs_pvs_fil,aes(y='',x=v2_v1))
  
  p + geom_beeswarm(aes(y='',x=v2_v1), color = "darkgrey", alpha = 0.9, method='compactswarm', size = 3) +
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000)) +
    geom_point(data=vert_mbs_pvs_fil %>%filter( name_2 %in% target), color="#3c5412", size = 3) +
    theme_bw()

  ggplot(vert_mbs_pvs_fil, aes(x =log10(v2_v1) , y="")) +
    geom_boxplot(width=0.1, outlier.colour=rgb(.5,.5,.5, 0)) +
    geom_violin(trim=FALSE, fill="#3c5412", alpha = 0.3, color = rgb(.5,.5,.5, 0.5)) + 
    geom_jitter(alpha = 0.7, position = position_jitter(seed = 1, width = 1, height = 0.05), aes(color = ifelse(v2_v1 >= 13537.00|v2_v1 <= 24787625.00, "black", NA))) +
    geom_jitter(alpha = 0.7, position = position_jitter(seed = 1, width = 1, height = 0.05), aes(color = ifelse(v2_v1 < 4.034271e-08|v2_v1 > 24787625.00, "#699320", NA))) +
    geom_text_repel(aes(label = ifelse(v2_v1 < 13537.00|v2_v1 > 24787625.00, name_2, "")), position = position_jitter(seed = 1, width = 1, height = 0.05)) +
    scale_color_identity() +
    theme_bw(base_size = 15)
  
  
  ##############
  p <- ggplot(vert_mbs_pvs_fil,aes(y='',x=V1_v2))
  
  p + geom_beeswarm(aes(y='',x=V1_v2), color = "darkgrey", alpha = 0.9, method='compactswarm', size = 3) +
    scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000)) +
    geom_point(data=vert_mbs_pvs_fil %>%filter( name_2 %in% target), color="#3c5412", size = 3) +
    theme_bw()
  
  ggplot(vert_mbs_pvs_fil, aes(x =log10(V1_v2) , y="")) +
    geom_boxplot(width=0.1, outlier.colour=rgb(.5,.5,.5, 0)) +
    geom_violin(trim=FALSE, fill="#3c5412", alpha = 0.3, color = rgb(.5,.5,.5, 0.5)) + 
    geom_jitter(alpha = 0.7, position = position_jitter(seed = 2, width = 1, height = 0.05), aes(color = ifelse(V1_v2 >= 4.034271e-08|V1_v2 <= 7.387161e-05, "black", NA))) +
    geom_jitter(alpha = 0.7, position = position_jitter(seed = 2, width = 1, height = 0.05), aes(color = ifelse(V1_v2 < 4.034271e-08|V1_v2 > 7.387161e-05, "#699320", NA))) +
    geom_text_repel(aes(label = ifelse(V1_v2 < 4.034271e-08|V1_v2 > 7.387161e-05, name_2, "")), position = position_jitter(seed = 2, width = 1, height = 0.05), max.overlaps = 10) +
    scale_color_identity() +
    theme_bw(base_size = 15)
  
  

ggplot(vert_mbs_pvs_fil, aes(x=v2_v1))+
  geom_histogram(binwidth = 0.05, fill="#bebada", color="#e9ecef", alpha=0.9) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  # geom_density(aes(y=0.1*..count..)) +
  # annotation_logticks(sides="b") +
  theme_bw()                                       



ggplot(vert_mbs_pvs_agg_no_zero, aes(x=logmb, y=logpv)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point() +
  geom_smooth(method=lm, level=0.95, color="#3c5412", fill = "#b3de69", formula=y~x+0) +
  theme_bw()

vert_mbs_pvs_nz <- filter(vert_mbs_pvs, !V1 == 0)

ggplot(vert_mbs_pvs, aes(x=V2, y=V1)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  stat_poly_line(color = "darkgreen") +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point() +
  theme_bw()

vert_mbs_pvs_genus <- separate(data = vert_mbs_pvs, col = name_2, into = c("genus", "species"), sep = " ")
vert_mbs_pvs_genus <- select(vert_mbs_pvs_genus, genus, V2, V1)
vert_mbs_pvs_genus$genus <- gsub("\\(", "", vert_mbs_pvs_genus$genus)
vert_mbs_pvs_genus_agg <- aggregate(.~genus, vert_mbs_pvs_genus, sum)
 
vert_mbs_pvs_genus_agg_nz <- filter(vert_mbs_pvs_genus_agg, !V1 == 0)
vert_mbs_pvs_genus_agg_nz2 <- filter(vert_mbs_pvs_genus_agg, !V2 == 0)n_occur <- data.frame(table(vocabulary$id))




ggplot(vert_mbs_pvs_nz, aes(x=V2/1000, y=V1)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  stat_poly_line(color = "darkgreen") +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point(color = "darkgreen") +
  geom_smooth(method=lm, level=0.95, color="#3c5412", fill = "#b3de69", formula=y~x+0) +
  theme_bw()

vert_mbs_pvs_nz$ratio <- vert_mbs_pvs_nz$V1/vert_mbs_pvs_nz$V2
