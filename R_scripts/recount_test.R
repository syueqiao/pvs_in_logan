#inferferon
# BiocManager::install("recount")
# BiocManager::install("edgeR")
# install.packages("remotes")
# remotes::install_github("ss-lab-cancerunit/SpiderSeqR")
# BiocManager::install("GEOquery")
# install.packages("XML")

library(XML)
library("SpiderSeqR")
library(recount)
library(DESeq2)
library(limma)
library("edgeR")
library(GEOquery)
library(tidyverse)

#go from sra run acession to the bioproject
link <- find_geo("SRR14508844")
#extract information from the link and associated gse
link_gse <- geo_info(link)
link_gse_in <- as.character(link_gse$series)

#get list of all processed rna seq counts 
geo_lines <- as.data.frame(readLines("gds_result (1).txt"))
sum(str_count(geo_lines$`readLines("gds_result (1).txt")`, "Accession"))
geo_fil <- geo_lines %>% filter(rowSums(across(everything(), ~ str_detect(., "Accession"))) > 0)
colnames(geo_fil) <- c("V1")
geo_fil$V2 <- base::gsub("Series\t\tAccession: ", "", geo_fil$V1)
geo_fil$V2<-gsub("\tID: .*", "", geo_fil$V2)

library_infg_list <- list()

for (j in head(geo_fil$V2, 5)) {
  #extract the links
  link_gse_url <- getGSEDownloadPageURLs(j)
  read_in_matrix <- getRNAQuantRawCountsURL(link_gse_url)
  annot_link <- getRNAQuantAnnotationURL(link_gse_url)
  
  #read in from the url, without reading to disk
  df_in <- as.data.frame(data.table::fread(read_in_matrix))
  df_in_annot <- as.data.frame(data.table::fread(annot_link))
  
  #convert to correct matrix format for deseq2
  mat_data <- data.matrix(df_in[,2:ncol(df_in)])
  rnames <- df_in[, 1]
  rownames(mat_data) <- rnames 
  mat_data <- as.matrix(mat_data)
  
  #get annotation data
  sum_exp <- getRNASeqData(j)
  #this could potentially work, but, given that i have to "add" my own metadata field
  #wont be as nice
  # ddsSE <- DESeqDataSet(sum_exp, design = ~ cell + dex)
  
  #get dataframe
  mat <- as.data.frame(sum_exp@assays@data)
  mat <- subset(mat, select=-c(group, group_name))
  
  col <- as.data.frame(sum_exp@colData@listData)
  
  #convert to appropriate matrix format
  col_data <- subset(col, select=-c(geo_accession))

  col_data_rnames <- as_vector(select(col, geo_accession))
  rownames(col_data) <- col_data_rnames
  
  #combine count data with metadata
  # ~1 is for no metadata, since I just want to look at counts right now
  dds <- DESeqDataSetFromMatrix(countData = mat_data,
                                colData = col_data,
                                design = ~ 1)
  
  dds <- DESeq(dds)
  
  deseq_norm_counts <- counts(dds, normalized=T)
  #gene id for infG
  # infg_row <- deseq_norm_counts[c("3458"),]
  infg_row <- deseq_norm_counts
  library_infg_list[[j]] <- infg_row # add it to your list
  
}
  
#deseq process
  res <- results(dds)
  return(res)
  
  #visualize if it worked
  resLFC <- lfcShrink(dds, coef="source_name_ch1_Skin_vs_Blood", type="apeglm")
  vsd <- vst(dds, blind=FALSE)
  plotPCA(vsd, intgroup=c("source_name_ch1"))
  #plotMA overlaps with another package, needs to be specified
  DESeq2::plotMA(resLFC, ylim=c(-2,2))
  
  plotCounts(dds, gene=which.min(res$padj), intgroup="title_rep")
  resOrdered <- res[order(res$pvalue),]
  
  resSig <- subset(resOrdered, padj < 0.1)
  change_rank <- resSig$log2FoldChange
  

#nice

#so, now to annotate PV+ libraries


######################################RUN FOR FUNCTIONS. THANK YOU SEAN.######################
#' get all download links from a GEO accession
#'
#' This function gets all download links from a GEO accession number.
#'
#' @param gse GEO accession number
#'
#' @return A character vector with all download links
#'
#' @keywords internal
getGSEDownloadPageURLs <- function(gse) {
  url <- "https://ncbi.nlm.nih.gov/geo/download/"
  links <- httr2::request(url) |>
    httr2::req_timeout(15) |>
    httr2::req_retry(3) |>
    httr2::req_url_query(acc = gse) |>
    httr2::req_perform() |>
    httr2::resp_body_string() |>
    rvest::read_html() |>
    rvest::html_nodes("a") |>
    rvest::html_attr("href") |>
    stringr::str_replace("^/geo/", "https://www.ncbi.nlm.nih.gov/geo/") |>
    stringr::str_replace("^ftp://", "https://")
  class(links) <- c("geoDownloadLinks", class(links))
  return(links)
}


#' Get the link to the raw counts file from GEO
#'
#' This function extracts the link to the raw counts file from a
#' geoDownloadLinks object.
#'
#' @param links A geoDownloadLinks object
#'
#' @return A character vector with the link to the raw counts file
#'
#' @keywords internal
getRNAQuantRawCountsURL <- function(links) {
  if (!inherits(links, "geoDownloadLinks")) {
    stop("Input must be a geoDownloadLinks object")
  }
  link <- stringr::str_subset(links, "_raw_counts_GRCh38")
  return(link)
}

#' Get the RNA-seq quantification annotation link
#'
#' This function extracts the link to the RNA-seq quantification annotation
#' file from a geoDownloadLinks object.
#'
#' @param links A geoDownloadLinks object
#'
#' @return A character vector with the link to the annotation file
#'
#' @keywords internal
getRNAQuantAnnotationURL <- function(links) {
  if (!inherits(links, "geoDownloadLinks")) {
    stop("Input must be a geoDownloadLinks object")
  }
  link <- stringr::str_subset(links, "annot.tsv.gz")
  return(link)
}

#' Extract filename from a GEO download URL
#'
#' This function extracts the filename from a GEO download URL.
#' The filename is expected to be a query parameter called "file".
#' If the query parameter is not found, the function returns NULL.
#'
#' The idea is to use this function to extract filenames that
#' contain important metadata from the GEO RNA-seq quantification.
#'
#' In particular, the filename is expected to contain the genome build
#' and species information that we can attach to the SummarizedExperiment.
#'
#' @param url A GEO download URL
#'
#' @return A character vector with the filename
#'
#' @keywords internal
extractFilenameFromDownloadURL <- function(url) {
  # get a query parameter called "file" and its value from
  # the URL
  
  # example URL: https://www.ncbi.nlm.nih.gov/geo/download/\
  # ?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz
  parsed <- httr2::url_parse(url)
  fname <- NULL
  if ("query" %in% names(parsed)) {
    if ("file" %in% names(parsed$query)) {
      fname <- parsed$query$file
    }
  }
  fname
}

#' Extract genome build and species from a GEO download URL
#'
#' This function extracts the genome build and species information
#' from a GEO download URL. The genome build and species information
#' is expected to be in the filename of the download URL.
#'
#' @param url A GEO annotation file download URL
#'
#' @return A character vector with the genome build and species information
#'
#' @keywords internal
urlExtractRNASeqQuantGenomeInfo <- function(url) {
  fname <- extractFilenameFromDownloadURL(url)
  if (is.null(fname)) {
    return(NULL)
  }
  splits <- stringr::str_split(fname, "\\.")[[1]]
  genome_build <- paste0(splits[2], ".", splits[3])
  species <- splits[1]
  return(c(genome_build = genome_build, species = species, fname = fname))
}

#' Extract genome build and species for GEO RNA-seq quantification
#'
#' This function extracts the genome build and species information
#' for a GEO RNA-seq quantification.
#'
#' @param gse GEO accession number
#'
#' @return A character vector with the genome build and species information
#'
#' @examples
#' getRNASeqQuantGenomeInfo("GSE164073")
#'
#' @export
getRNASeqQuantGenomeInfo <- function(gse) {
  links <- getGSEDownloadPageURLs(gse)
  annotation_link <- getRNAQuantAnnotationURL(links)
  metadata <- urlExtractRNASeqQuantGenomeInfo(annotation_link)
  return(metadata)
}


#' Read RNA-seq quantification annotation from GEO
#'
#' This function reads the annotation file from a GEO link. The annotation
#' file is expected to be a tab-separated file with the first column
#' containing the gene IDs and the remaining columns containing the
#' annotation information.
#'
#' @param link A link to the annotation file
#'
#' @return A data frame of annotation information with gene IDs as row names
#'
#' @keywords internal
readRNAQuantAnnotation <- function(link) {
  annotation <- as.data.frame(readr::read_tsv(link, show_col_types = FALSE))
  rownames(annotation) <- as.character(annotation$GeneID)
  return(annotation)
}


#' Read raw counts from GEO
#'
#' This function reads the raw counts from a GEO link. The raw counts
#' are expected to be in a tab-separated file with the first column
#' containing the gene IDs and the remaining columns containing the
#' raw counts.
#'
#' This function reads the raw counts and returns a matrix with the
#' gene IDs as the row names, ready for use in creating a
#' SummarizedExperiment.
#'
#' @param link A link to the raw counts file
#'
#' @return A matrix of raw counts with gene IDs as row names
#'
#' @keywords internal
readRNAQuantRawCounts <- function(link) {
  quants <- readr::read_tsv(link, show_col_types = FALSE)
  gene_ids <- as.character(quants$GeneID)
  quants <- as.matrix(quants[, -1])
  rownames(quants) <- gene_ids
  return(quants)
}


#' Get RNA-seq quantification and annotation from GEO
#'
#' This function downloads the raw counts and annotation files from GEO
#' for a given GEO accession number.
#'
#' @param gse GEO accession number
#'
#' @return A list with two elements: quants (a matrix of raw counts) and
#' annotation (a data frame of annotation information).
#'
#' @keywords internal
getRNASeqQuantResults <- function(gse) {
  links <- getGSEDownloadPageURLs(gse)
  raw_counts_link <- getRNAQuantRawCountsURL(links)
  if (length(raw_counts_link) == 0) {
    stop(
      "No raw counts file found.\n",
      "Navigate to: \n  https://ncbi.nlm.nih.gov/geo/download/?acc=",
      gse,
      "\nand check if the 'RNA-Seq raw counts' link is available."
    )
  }
  annotation_link <- getRNAQuantAnnotationURL(links)
  quants <- readRNAQuantRawCounts(raw_counts_link)
  annotation <- readRNAQuantAnnotation(annotation_link)
  return(list(quants = quants, annotation = annotation))
}

#' Browse GEO search website for RNA-seq datasets
#'
#' This function opens a browser window to the NCBI GEO website
#' with a search for RNA-seq datasets. It is included as a convenience
#' function to remind users of how to search for RNA-seq datasets using
#' the NCBI GEO website and an "rnaseq counts" filter.
#'
#' @examples
#' \dontrun{
#' browseWebsiteRNASeqSearch()
#' }
#' @export
browseWebsiteRNASeqSearch <- function() {
  browseURL(
    "https://ncbi.nlm.nih.gov/gds?term=%22rnaseq%20counts%22%5BFilter%5D"
  )
}


#' Does a GEO accession have RNA-seq quantifications?
#'
#' This function checks if a GEO accession number has RNA-seq quantifications
#' available. It does this by checking if the GEO accession number has a
#' "RNA-Seq raw counts" link available on the GEO download page.
#'
#' @param accession GEO accession number
#'
#' @return TRUE if the GEO accession number has RNA-seq quantifications
#' available, FALSE otherwise.
#'
#' @examples
#' hasRNASeqQuantifications("GSE164073")
#'
#' @export
hasRNASeqQuantifications <- function(accession) {
  links <- getGSEDownloadPageURLs(accession)
  raw_counts_link <- getRNAQuantRawCountsURL(links)
  if (length(raw_counts_link) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Get GEO RNA-seq quantifications as a SummarizedExperiment object
#'
#' For human and mouse GEO datasets, NCBI GEO attempts to process
#' the raw data and provide quantifications in the form of raw counts
#' and an annotation file. This function downloads the raw counts and
#' annotation files from GEO and merges that with the metadata from the GEO
#' object to create a SummarizedExperiment.
#'
#' @details
#' A major barrier to fully exploiting and reanalyzing the massive volumes
#' of public RNA-seq data archived by SRA is the cost and effort required to
#' consistently process raw RNA-seq reads into concise formats that summarize
#' the expression results. To help address this need, the NCBI SRA and GEO
#' teams have built a pipeline that precomputes RNA-seq gene expression counts
#' and delivers them as count matrices that may be incorporated into commonly
#' used differential expression analysis and visualization software.
#'
#' The pipeline processes RNA-seq data from SRA using the HISAT2 aligner and
#' and then generates gene expression counts using the featureCounts program.
#'
#' See the
#' [GEO documentation](https://ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html)
#' for more details.
#'
#' @import SummarizedExperiment
#'
#'
#' @param accession GEO accession number
#'
#' @return A SummarizedExperiment object with the raw counts as the counts
#' assay, the annotation as the rowData, and the metadata from GEO as
#' the colData.
#'
#'
#' @examples
#' se <- getRNASeqData("GSE164073")
#' se
#'
#' @export
getRNASeqData <- function(accession) {
  quantifications <- getRNASeqQuantResults(accession)
  se <- as(
    getGEO(accession)[[1]],
    "SummarizedExperiment"
  )
  se <- se[,se@colData$geo_accession %in% colnames(quantifications$quants)]
  old_metadata <- S4Vectors::metadata(se)
  old_metadata$genomeInfo <- getRNASeqQuantGenomeInfo(accession)
  old_metadata$created_at <- Sys.time()
  new_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = quantifications$quants),
    rowData = quantifications$annotation,
    colData = SummarizedExperiment::colData(se),
    metadata = old_metadata
  )
  new_se
}

