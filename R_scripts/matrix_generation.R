library(tidyverse)
library(ggplot2)
#import dataset with fasta headers
genomes_tsv_df <- read.table("C:/Users/Jessica Shen/Documents/rnalab_2024/20_test_set_run1_folder/ncbi_dataset/cat_20_test_set_run1_headers.tsv", sep = "\t" , header = F, quote=NULL, comment='')
#clean because ncbi..why
matrix_maker <- function(genomes_tsv_df){
genomes_extracted <- genomes_tsv_df
genomes_extracted$protein <- sub(".*(protein=.*)", "\\1", genomes_extracted$V1)
genomes_extracted$gene <- sub(".*(gene=.*)", "\\1", genomes_extracted$V1)

# str_extract(genomes_extracted, "(?<=22k_).*(?=_res)")
# genomes_sep <- separate(genomes_tsv_df, 1, into=c("ids", "gene", "locus_tag", "db"), sep= " ")
genomes_extracted$V1 <- gsub( " .*$", "", genomes_extracted$V1 )
genomes_sep <- genomes_extracted
colnames(genomes_sep) <- c("ids", "protein", "gene")
genomes_sep$ids <- gsub(".*\\|","",genomes_sep$ids)%>%str_remove(., "_cds")
genomes_sep_g_p <- separate(genomes_sep, ids, into = c("genome_id", "protein_id"), sep = "(?<=[0-9])_")
genomes_sep_g_p$protein <- sub("].*", "", genomes_sep_g_p$protein)
genomes_sep_g_p$gene[] <- lapply(genomes_sep_g_p$gene, function(x) replace(x, grep("[|]", x), NA))
genomes_sep_g_p$gene <- sub("].*", "", genomes_sep_g_p$gene)
# genomes_sep_g_p[] <- lapply(genomes_sep_g_p, gsub, pattern='\\[', replacement='')
# genomes_sep_g_p[] <- lapply(genomes_sep_g_p, gsub, pattern='\\]', replacement='')
# genomes_list <- group_split(genomes_sep_g_p, genome_id)

#dealing with edge cases...such as weird annotating from ncbi

pv_genes_joins <- c("E1\\^E4", "E2\\^E8")
indx1 <- paste(pv_genes_joins, collapse="|")
genomes_sep_joins <- genomes_sep_g_p[grep(indx1, genomes_sep_g_p$protein),]
genomes_sep_joins$protein_extract <- str_extract(genomes_sep_joins$protein, indx1)
genomes_sep_joins$gene_extract <- str_extract(genomes_sep_joins$gene, indx1)

pv_genes <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "L1", "L2", "L3")
indx2 <- paste(pv_genes, collapse="|")
genomes_sep_g_p_nj <- genomes_sep_g_p[genomes_sep_g_p$protein_id != "YP_009268708.1", ]
genomes_sep_g_p_nj <- genomes_sep_g_p_nj[grep(indx2, genomes_sep_g_p_nj$protein),]
genomes_sep_g_p_nj$protein_extract <- str_extract(genomes_sep_g_p_nj$protein, indx2)
genomes_sep_g_p_nj <- genomes_sep_g_p_nj[grep(indx2, genomes_sep_g_p_nj$gene),]
genomes_sep_g_p_nj$gene_extract <- str_extract(genomes_sep_g_p_nj$gene, indx2)


# genomes_sep_g_p$test <- regmatches(genomes_sep_g_p$gene, 
#                                    gregexpr(pv_genes, genomes_sep_g_p$gene))
#characterize which genes are present in each genome
genomes_sep_final <- rbind(genomes_sep_g_p_nj, genomes_sep_joins)
genomes_sep_clean <- genomes_sep_final[!duplicated(genomes_sep_final$protein_id),]
genomes_sep_clean
# matrix <- table(select(genomes_sep_clean, genome_id, gene))
# matrix
}
test <- matrix_maker(genomes_tsv_df)
LARGE_genomes_tsv_df <- read.table("C:/Users/Jessica Shen/Documents/rnalab_2024/bweee_folder/ncbi_dataset/cat_bweee_headers.tsv", sep = "\t" , header = F, quote=NULL, comment='')
large_table_lol <- matrix_maker(LARGE_genomes_tsv_df)
write.csv(large_table_lol, "large_table.csv")
write.csv(test$genome, "rs_accession_js.csv")
