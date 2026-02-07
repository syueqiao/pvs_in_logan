# hmmsearch_utils.R
# Shared utility functions for parsing HMMER output files
#
# Usage: source("00_utilities/hmmsearch_utils.R")

library(tidyverse)

#' Parse HMMER domtblout file
#'
#' Reads and cleans a HMMER domain table output file (.domtbl),
#' extracting relevant columns and converting to proper types.
#'
#' @param input_domtbl Path to the HMMER domtblout file
#' @return A data frame with cleaned HMMER results
#'
#' @examples
#' results <- hmmsearch_clean("output.domtbl")
#'
hmmsearch_clean <- function(input_domtbl) {
  input_domtbl_scan <- read.table(input_domtbl, sep = "", header = FALSE, fill = TRUE) %>%
    na.omit() %>%
    .[, 1:22]

  # Filter out header lines and incomplete rows
  input_domtbl_scan <- input_domtbl_scan[
    !(is.na(input_domtbl_scan$V21) |
      input_domtbl_scan$V21 == "" |
      !(input_domtbl_scan$V2 == "-")),
  ]

  colnames(input_domtbl_scan) <- c(
    "query_acc",      # Query sequence accession
    "misc",           # Miscellaneous (usually "-")
    "qlen",           # Query sequence length
    "pfam",           # Target HMM name
    "pfam_acc",       # Target HMM accession
    "tlen",           # Target HMM length
    "eval_full",      # Full sequence E-value
    "score_full",     # Full sequence score
    "bias_full",      # Full sequence bias
    "#",              # Domain number (this domain)
    "of",             # Total domains (of this many)
    "c_eval",         # Conditional E-value
    "i_eval",         # Independent E-value
    "score_one",      # Domain score
    "bias_one",       # Domain bias
    "hmmcoord_from",  # HMM coordinate start
    "hmmcoord_to",    # HMM coordinate end
    "alicoord_from",  # Alignment coordinate start
    "alicoord_to",    # Alignment coordinate end
    "envcoord_from",  # Envelope coordinate start
    "envcoord_to",    # Envelope coordinate end
    "acc"             # Accuracy
  )

  # Clean query accession (remove trailing ORF identifier)
  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)

  # Convert columns to appropriate types
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)

  return(input_domtbl_scan)
}


#' Parse HMMER domtblout file with additional info column
#'
#' Similar to hmmsearch_clean() but reads 23 columns instead of 22,
#' useful when the domtbl contains an additional description column.
#'
#' @param input_domtbl Path to the HMMER domtblout file
#' @return A data frame with cleaned HMMER results including misc2 column
#'
hmmsearch_clean_info <- function(input_domtbl) {
  input_domtbl_scan <- read.table(input_domtbl, sep = "", header = FALSE, fill = TRUE) %>%
    na.omit() %>%
    .[, 1:23]

  input_domtbl_scan <- input_domtbl_scan[
    !(is.na(input_domtbl_scan$V21) |
      input_domtbl_scan$V21 == "" |
      !(input_domtbl_scan$V2 == "-")),
  ]

  colnames(input_domtbl_scan) <- c(
    "query_acc", "misc", "qlen", "pfam", "pfam_acc", "tlen",
    "eval_full", "score_full", "bias_full", "#", "of",
    "c_eval", "i_eval", "score_one", "bias_one",
    "hmmcoord_from", "hmmcoord_to", "alicoord_from", "alicoord_to",
    "envcoord_from", "envcoord_to", "acc", "misc2"
  )

  input_domtbl_scan$query_acc_clean <- sub("_[^_]+$", "", input_domtbl_scan$query_acc)
  input_domtbl_scan <- type.convert(input_domtbl_scan, as.is = TRUE)

  return(input_domtbl_scan)
}


#' Parse Diamond .pro output file
#'
#' Reads Diamond protein alignment output in tabular format.
#'
#' @param input_file Path to the Diamond .pro file
#' @return A data frame with Diamond alignment results
#'
clean_inputs <- function(input_file) {
  diamond_output <- read.table(input_file, sep = "\t", header = FALSE, fill = TRUE)

  colnames(diamond_output) <- c(
    "qseqid",    # Query sequence ID
    "qstart",    # Query start
    "qend",      # Query end
    "qlen",      # Query length
    "sseqid",    # Subject sequence ID
    "sstart",    # Subject start
    "send",      # Subject end
    "slen",      # Subject length
    "pident",    # Percent identity
    "evalue"     # E-value
  )

  return(diamond_output)
}
