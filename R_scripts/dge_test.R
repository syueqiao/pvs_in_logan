BiocManager::install("DESeq2")
BiocManager::install("apeglm")


library("DESeq2")
library("apeglm")


pasCts <- read.table(gzfile("GSE191090_raw_counts_GRCh38.p13_NCBI.tsv.gz"), sep = "\t", header = T, row.names = "GeneID")


cts <- as.matrix(pasCts)

coldata <- read.table("conditions2.txt", sep = "\t", header = F, row.names = 1, fill = T)
coldata$gene <- str_sub(coldata$V3, -2, -1)
coldata$pv <- str_sub(coldata$V3, -4, -3)

coldata <- filter(coldata, gene == "E7")
  
colnames(coldata) <- c("cells", "condition", "rep", "gene", "pv")
coldata$condition <- gsub("-", "_", coldata$condition)

coldata$condition <- factor(coldata$condition)
coldata$rep <- factor(coldata$rep)

cts <- cts[, colnames(cts) %in%  as.character(rownames(coldata))]


all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ rep + pv)
dds

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
# res_all <- results(dds, contrast=c("gene", "E5", "E6"))
res_all <- results(dds, contrast=c("pv", "84", "16"))

res_all
sum(res_all$padj < 0.05, na.rm=TRUE)
summary(res_all)



sum(res_all$padj < 0.0001, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="pv_pu_vs_16", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))



resOrdered <- res_all[order(res_all$pvalue),]
resSig <- subset(resOrdered, padj < 0.1)
write.csv(resSig, "resSig_E7.csv")

plotCounts(dds, gene=which(rownames(res_all) == "8638") , intgroup="pv")

install.packages("pheatmap")
library("pheatmap")

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]

df <- as.data.frame(colData(dds)[,c("pv", "gene")])
rownames(df) <- colnames(cts)


ntd <- normTransform(dds)

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("rep"))

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$pv, vsd$rep, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

plotPCA(vsd, intgroup=c("rep", ""))
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$rep)
plotPCA(vsd, "rep")
plotPCA(vsd, "pv")

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$pv, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


dds2 <- DESeq(dds, test="LRT", reduced=~rep)
res2 <- results(dds2)
summary(res2)
res2

plotMA(res2, ylim=c(-2,2))


alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res_all$log2FoldChange, -log10(res_all$pvalue))
plot(res_all$log2FoldChange, -log10(res_all$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res_all$log2FoldChange) > 2.5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.4)


par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- as.data.frame(res_all)
http://127.0.0.1:32037/graphics/plot_zoom_png?width=654&height=625
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<0.05 & log2FoldChange>=1), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))
with(subset(topT, padj<0.05 & log2FoldChange< -1), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

sign.genes=which(topT$padj<0.05)

text(x=topT$log2FoldChange[sign.genes] , y=-log10(topT$pvalue[sign.genes]), label=row.names(topT)[sign.genes], cex=0.5)

#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1, col="black", lty=4, lwd=2.0)
abline(v=1, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

library(EnhancedVolcano)

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL

EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5e-2,
                title = '84 vs 16',
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 2.0, drawConnectors = TRUE, max.overlaps = 10)

write.csv(resSig, "resSig_E7.csv")

resOrdered <- res_all[order(res_all$pvalue),]
resSig <- subset(resOrdered, padj < 0.1)
plotMA(res2, ylim=c(-2,2))
