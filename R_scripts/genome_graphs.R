#genome vis
library(tidyverse)
 library(ggplot2)
library(gggenes)
library(viridis)
 
 rhino <- read.table("rhino.txt", sep ="\t", header = T)
 
 ggplot(rhino, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
   geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
   facet_wrap(~ molecule, scales = "free", ncol = 1) +
   scale_fill_brewer(palette = "Set3") +
   ylab("Contig") +
   guides(fill=guide_legend(title="Gene")) +
   theme_genes()
 
 ggsave("rhino.png", width = 20, height = 10, units = "cm", limitsize = F)
 

manidae <- read.table("manidae_2.txt", sep ="\t", header = T)

ggplot(manidae, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("manidae_2.png", width = 20, height = 10, units = "cm", limitsize = F)


ggsave("manidae_2.png")


hawk <- read.table("hawk.txt", sep ="\t", header = T)

ggplot(hawk, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_viridis(discrete = T)+
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("hawk.png")

fox <- read.table("vulpes.txt", sep ="\t", header = T)

ggplot(fox, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_viridis(discrete = T)+
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("fox.png")

rat <- read.table("rat.txt", sep ="\t", header = T)

ggplot(manidae, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_viridis(discrete = T)+
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("rat.svg")

rat <- read.table("rat.txt", sep ="\t", header = T)

ggplot(rat, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_viridis(discrete = T)+
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("rat.svg")

tick <- read.table("tick.txt", sep ="\t", header = T)

tick <- ggplot(tick, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_viridis(discrete = T)+
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("tick.png")

###########weird mmtv

SRR830327_9 <- read.table("SRR830327_9.txt", sep ="\t", header = T)

ggplot(SRR830327_9, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

SRR830321_224 

SRR830321_224 <- read.table("SRR830321_224.txt", sep ="\t", header = T)

ggplot(SRR830321_224, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()


NODE_1_length_8802_cov_9.743642 <- read.table("NODE_1_length_8802_cov_9.743642.txt", sep ="\t", header = T)

ggplot(NODE_1_length_8802_cov_9.743642, aes(xmin = start, xmax = end, y = molecule, fill = origin, label = gene, alpha = 0.1)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  ylab("Contig") +
  guides(label=guide_legend(title="Gene")) +
  scale_fill_viridis(discrete = T)+
  theme_genes()


NODE_1_length_8802_cov_9.743642_SRR830310 <- read.table("NODE_1_length_8802_cov_9.743642_SRR830310.txt", sep ="\t", header = T)

ggplot(NODE_1_length_8802_cov_9.743642_SRR830310, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  ylab("Contig") +
  guides(label=guide_legend(title="Gene")) +
  scale_fill_brewer(palette = "Set3")+
  theme_genes()

NODE_1_length_8802_cov_9.743642_SRR830298 <- read.table("NODE_1_length_8802_cov_9.743642_SRR830298.txt", sep ="\t", header = T)

ggplot(manidae, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  ylab("Contig") +
  guides(label=guide_legend(title="Gene")) +
  scale_fill_brewer(palette = "Set3")+
  theme_genes()


fish <- read.table("fish.txt", sep ="\t", header = T)

ggplot(fish, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("fish_genes.png", width = 20, height = 10, units = "cm", limitsize = F)

bear <- read.table("bear.txt", sep ="\t", header = T)

ggplot(bear, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set2") +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

ggsave("bear_genes.svg", width = 20, height = 10, units = "cm", limitsize = F)


plot_from_table <- function(gene_table){
gene_table <- read.table(gene_table, sep ="\t", header = T)

p <- ggplot(gene_table, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, alpha=0.5)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  ylab("Contig") +
  guides(fill=guide_legend(title="Gene")) +
  theme_genes()

return(p)}

plot_from_table("SRR22213301_62214.txt")
