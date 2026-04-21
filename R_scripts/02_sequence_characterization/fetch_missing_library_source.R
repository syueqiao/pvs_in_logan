library(tidyverse)

# =====================================================================
# Fetch library_source for accessions missing from lib_source_big.txt
# Uses the ENA Portal API (mirrors all SRA data, no API key needed)
# =====================================================================

# --- 1. Identify missing accessions ---
domtbl_raw <- read.table("files/feb_7_pv_fil_form.aa.domtbl",
                         comment.char = "#", fill = TRUE,
                         stringsAsFactors = FALSE)

domtbl_srrs <- unique(sub("^([SED]RR[0-9]+)_.*", "\\1", domtbl_raw$V1))
cat("Unique accessions in domtbl:", length(domtbl_srrs), "\n")

# Load all existing library source files
lib1 <- read.table("files/lib_source_big.txt", sep = "\t",
                    header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("srr", "library_source"),
                    fill = TRUE, quote = "")

lib2 <- read.table("files/L1_library_metadata.txt", sep = "\t",
                    header = FALSE, stringsAsFactors = FALSE,
                    fill = TRUE, quote = "") %>%
  select(srr = V1, library_source = V10) %>%
  filter(!is.na(library_source), library_source != "")

lib3 <- read.table("files/hits_library_biosample.list", sep = "\t",
                    header = TRUE, stringsAsFactors = FALSE) %>%
  select(srr = Run, library_source = LibrarySource) %>%
  filter(!is.na(library_source), library_source != "")

# Combine existing sources (take first non-NA match per accession)
existing <- bind_rows(lib1, lib2, lib3) %>%
  filter(!is.na(library_source), library_source != "") %>%
  distinct(srr, .keep_all = TRUE)

cat("Existing library source entries:", nrow(existing), "\n")

missing_srrs <- setdiff(domtbl_srrs, existing$srr)
cat("Missing accessions to fetch:", length(missing_srrs), "\n")

# --- 2. Batch query ENA Portal API ---
# ENA mirrors all SRA/DRA/ERA data and supports bulk queries
# Endpoint: https://www.ebi.ac.uk/ena/portal/api/search
# We query in batches of 200 accessions

fetch_ena_batch <- function(accessions) {
  # ENA query: accession="X" OR accession="Y" ...
  query_str <- paste0('run_accession="',
                       paste(accessions, collapse = '" OR run_accession="'),
                       '"')
  url <- paste0("https://www.ebi.ac.uk/ena/portal/api/search?",
                "result=read_run&",
                "fields=run_accession,library_source&",
                "format=tsv&",
                "limit=0")

  resp <- tryCatch({
    httr::POST(url,
               body = list(query = query_str),
               encode = "form",
               httr::timeout(120))
  }, error = function(e) {
    message("  Request failed: ", e$message)
    return(NULL)
  })

  if (is.null(resp) || httr::status_code(resp) != 200) {
    status <- ifelse(is.null(resp), "error", httr::status_code(resp))
    message("  HTTP ", status)
    return(tibble(run_accession = character(), library_source = character()))
  }

  content <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (nchar(trimws(content)) < 10) {
    return(tibble(run_accession = character(), library_source = character()))
  }

  tryCatch(
    read_tsv(content, show_col_types = FALSE),
    error = function(e) {
      message("  Parse error: ", e$message)
      tibble(run_accession = character(), library_source = character())
    }
  )
}

# Check httr is available
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")

batch_size <- 100
batches <- split(missing_srrs, ceiling(seq_along(missing_srrs) / batch_size))
cat("Total batches:", length(batches), "\n")

results <- list()
fetched_so_far <- 0
for (i in seq_along(batches)) {
  if (i %% 50 == 0 || i == 1) {
    cat(sprintf("  Batch %d / %d (%d accessions fetched so far)\n",
                i, length(batches), fetched_so_far))
  }
  res <- fetch_ena_batch(batches[[i]])
  results[[i]] <- res
  fetched_so_far <- fetched_so_far + nrow(res)

  # Rate limit: be polite to ENA
  Sys.sleep(0.5)
}

ena_results <- bind_rows(results)
cat("\nENA returned results for", nrow(ena_results), "accessions\n")

# --- 3. Check what's still missing after ENA ---
still_missing <- setdiff(missing_srrs, ena_results$run_accession)
cat("Still missing after ENA query:", length(still_missing), "\n")

# --- 4. Combine everything and save ---
ena_clean <- ena_results %>%
  select(srr = run_accession, library_source) %>%
  filter(!is.na(library_source), library_source != "")

all_lib_source <- bind_rows(existing, ena_clean) %>%
  distinct(srr, .keep_all = TRUE)

cat("\nFinal combined library source entries:", nrow(all_lib_source), "\n")
cat("Coverage of domtbl accessions:",
    sum(domtbl_srrs %in% all_lib_source$srr), "/",
    length(domtbl_srrs), "\n")

# Save the combined file
write.table(all_lib_source, "files/lib_source_combined.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("Saved: files/lib_source_combined.txt\n")

# Also save the still-missing list for reference
if (length(still_missing) > 0) {
  writeLines(still_missing, "outputs/lib_source_still_missing.list")
  cat("Saved missing list:", length(still_missing),
      "accessions -> outputs/lib_source_still_missing.list\n")
}
