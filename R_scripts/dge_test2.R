library("DESeq2")
library("apeglm")


pasCts <- read.table(gzfile("GSE191090_raw_counts_GRCh38.p13_NCBI.tsv.gz"), sep = "\t", header = T, row.names = "GeneID")


cts <- as.matrix(pasCts)

coldata <- read.table("conditions3.csv", sep = ",", header = F, row.names = 1, fill = T)
coldata <- coldata %>% 
  extract(V2, into = c("condition", "rep"), "(.*)_([^_]+)$")

# colnames(coldata) <- c("cells", "condition", "rep", "gene", "pv")
# coldata$condition <- gsub("-", "_", coldata$condition)



coldata$condition <- factor(coldata$condition)
coldata$rep <- factor(coldata$rep)
coldata$risk <- with(coldata, ifelse(condition %in% c("FBE6E7_18", "FBE6E7_16"), "high", ifelse(condition == "FBpLVX", "control", "low")))

cts <- cts[, colnames(cts) %in%  as.character(rownames(coldata))]


all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ rep + condition)
dds

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
# res_all <- results(dds, contrast=c("gene", "E5", "E6"))
res_all <- results(dds, contrast=c("condition", "FBE6E7_16", "FBE6E7_122"))

res_all
sum(res_all$padj < 0.05, na.rm=TRUE)
sum(abs(res_all$log2FoldChange)> 1 & res_all$padj < 0.05, na.rm=TRUE)

summary(res_all)

dds$condition <- relevel(dds$condition, ref = "FBE6E7_16")

dds$condition
dds <- nbinomWaldTest(dds)


resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_FBE6E7_107_vs_FBE6E7_16", type="apeglm")
plotMA(dds, ylim=c(-2,2))



resOrdered <- res_all[order(res_all$pvalue),]
resSig <- subset(resOrdered, padj < 0.1)
write.csv(resSig, "resSig_E7_3.csv")

vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("condition"))
pcaData <- plotPCA(vsd, intgroup=c("condition", "risk", "rep"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


mat <- assay(vsd)
mm <- model.matrix(~condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$rep, design=mm)
assay(vsd) <- mat

plotPCA(vsd, intgroup=c("condition"))
pcaData <- plotPCA(vsd, intgroup=c("condition", "risk", "rep"), returnData=TRUE)

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=risk)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_brewer(palette = "Set2") +
  theme_bw()

plotCounts(dds, gene=which(rownames(res_all) == "1026") , intgroup="condition")

EnhancedVolcano(res_all,
                lab = rownames(res_all),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5e-2,
                title = 'FBE6E7_16 vs FBE6E7_122',
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 2.0, drawConnectors = TRUE, max.overlaps = 10)


