library(tidyverse)
library(Polychrome)
library(ggtree)
library(ggplot2)
library(magrittr)
library(ggnewscale)
library(treeio)
library(ape)

#need to set this manually, some update to ggtree and ggplot that made this not play nice anymore??
is.waive <- function(x) {
  return(inherits(x, "waiver"))
}



#' Encapsulate path
#'
#' This function is used to encapsulate paths to allow for behaviour switches.
#' One use is for example when plotting. The plot_all method will encapsulate a
#' path so that plots may be saved to a directory structure. Other plot methods,
#' e.g. plot_model_performance do not encapsulate a path, and if the user calls
#' these functions directly, the plot may be written to the provided path
#' instead of a directory structure.
#'
#' @return encapsulated_path object
#' @md
#' @keywords internal
encapsulate_path <- function(path) {
  structure(path, class = "encapsulated_path")
}

#read in final tree file with correct counts...

final_tree <- read.tree("files/2026.04.01.tree_seqs_1375.treefile")

# Remove SRR6976994 tips (simple repeats)
drop_tips <- final_tree$tip.label[grepl("^SRR6976994", final_tree$tip.label)]
final_tree <- drop.tip(final_tree, drop_tips)

tree <- ggtree(final_tree)
# view(tree$data)
tree_tips <- filter(tree$data, isTip == T)
tree$data$label <- gsub("Edges.*", "", tree$data$label)
tree$data$label <- gsub("ka_f.*", "", tree$data$label)

# tree$data$x10 <- iconv(tree$data$x10, from = "UTF-8", to = "ASCII", sub = "")
# tree$data[1092, 10] <- "Zhangixalus_dugritei_associated"

tree$data$status <- ifelse(grepl("[E|S|D][R][R]", tree$data$label), "novel", "ncbi")

novel_tips <- filter(tree$data, status == "novel" & isTip == T)
check_tree <- tree$data
tree_plot <- tree + geom_aline(aes(color = status), linetype = "solid", linewidth = 0.5, position = position_nudge(x = -0.035)) + 
  scale_color_manual(values = c("grey","#4646d4")) +
  geom_tiplab(color = "black",size = 2, hjust = 0, align= T, linetype = "blank") +
  theme(legend.position= "none")

#add library source data
sra_metadata <- read.table("files/all_hits_info.list", sep = ",", fill = T, header = F, quote = "", comment.char = "#")
#fix single broken cell
# sra_metadata["144301", "V2"] <- "GENOMIC"
sra_metadata <- unique(sra_metadata)

#read in the first version of the annotations
addtl_ano <- read.table("files/2026.02.12all_genus_addtl_anno_fixed.csv", sep = ",", header = T)
addtl_ano$lib <- gsub("_.*", "", addtl_ano$label)
addtl_ano$lib <- gsub("Score.*", "", addtl_ano$lib)

addtl_ano <- left_join(addtl_ano, sra_metadata, by = c("lib" = "V1"))
addtl_tip <- filter(addtl_ano, isTip == "TRUE")
addtl_tip$broader_gen <- gsub("Cattle", "Bovine", addtl_tip$broader_gen)

addtl_tip[is.na(addtl_tip$broader_gen)] <- "Other"

addtl_tip$broader_gen <- factor(addtl_tip$broader_gen, levels = c("Human", "Bat", "Bovine",
"Canine", "Ray-finned Fish",
"Rodent", "Avian", 
"Cetacean", "Amphibian", 
"Cervine", "Feline", "Pangolin", 
"Non-human Primate", 
"Reptile", "Equine", 
"Porcine", "Other"))

addtl_tip[is.na(addtl_tip)] <- "Other"


summary_table <- addtl_tip %>%
  group_by(generalization, status) %>%
  summarize(
    Count = n() # n() counts the number of rows in each group
  )

df_wide <- pivot_wider(
  summary_table,
  names_from = status,
  values_from = Count
)

df_wide[is.na(df_wide)] <- 0

#add back the additional information from final accounting pass
addtl_tip_plus_new_manual <- read.table("outputs/2026.02.13.addtl_tip_plus_new.tsv.txt", sep = '\t', header = TRUE)
addtl_tip_plus_new_manual$library <- gsub("_.*", "", addtl_tip_plus_new_manual$V1)
addtl_tip_plus_new_manual$library <- gsub("Score.*", "", addtl_tip_plus_new_manual$library)

addtl_tip_plus_new_manual <- left_join(addtl_tip_plus_new_manual, sra_metadata, by = c("library" = "V1"))

addtl_tip_plus_new_manual[is.na(addtl_tip_plus_new_manual$broader_gen)] <- "Other"
check_table <- as.data.frame(table(addtl_tip_plus_new_manual$broader_gen, addtl_tip_plus_new_manual$status))


summary_table_updated <- addtl_tip_plus_new_manual %>%
  group_by(generalization, status) %>%
  summarize(
    Count = n() # n() counts the number of rows in each group
  )

df_wide_updated <- pivot_wider(
  summary_table_updated,
  names_from = status,
  values_from = Count
)

df_wide_updated[is.na(df_wide_updated)] <- 0

write.table(df_wide_updated, "outputs/2025.02.20.all_organism_group_assignments_updated.tsv", sep = '\t', quote = F, row.names = F, col.names = F)

check_table <- as.data.frame(table(addtl_tip_plus_new_manual$broader_gen))

addtl_tip_plus_new_manual$broader_gen <- factor(addtl_tip_plus_new_manual$broader_gen, levels = c("Human", "Bat", "Bovine",
"Canine", "Rodent","Ray-finned Fish",
 "Avian", 
"Amphibian","Feline","Cetacean",  
"Non-human Primate","Cervine", "Reptile",  "Equine", 
"Pangolin", 
"Pinnipeds", "Other"))

addtl_tip_plus_new_manual[is.na(addtl_tip_plus_new_manual)] <- "Other"

summary_table_updated <- addtl_tip_plus_new_manual %>%
  group_by(broader_gen, status) %>%
  summarize(
    Count = n() # n() counts the number of rows in each group
  )

df_wide_updated <- pivot_wider(
  summary_table_updated,
  names_from = status,
  values_from = Count
)

                         
write.table(df_wide, "outputs/all_organism_group_assignments_broad.tsv", sep = '\t', quote = F, row.names = F, col.names = F)


addtl_tip_novel <- subset(addtl_tip_plus_new_manual, grepl("[SED][R]{2}", addtl_tip_plus_new_manual$V1))

addtl_tip_novel_tb <- as.data.frame(table(addtl_tip_novel$broader_gen))

addtl_tip_novel_tb_ord <- addtl_tip_novel_tb %>%
  mutate(Var1 = forcats::fct_reorder(factor(Var1), Freq))

levels(factor(addtl_tip_novel_tb_ord$Var1)) 


#see if lollipop chart looks nice?
ggplot(addtl_tip_novel_tb_ord, aes(x=Var1, y=Freq)) +
  geom_segment( aes(x=Var1, xend=Var1, y=0, yend=Freq), color="grey20") +
  geom_point( color="#4646d4", size=2, alpha=0.9) +
  theme_bw() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+
  labs(title = "Novel Papillomaviruses by Associated Species",
       x = "Associated Species",
       y = "Count")

#it doesnt.

ggplot(addtl_tip_plus_new_manual, aes(x=specific_association, fill = status)) + 
  geom_bar(stat="count") + scale_fill_manual(values = c("grey", "#4646d4")) + theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  labs(title = "Known and Novel PVs\nAssociated with Specific Species",
       x = "Category",
       y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(text = element_text(family="Noto Sans", size = 20))

#thats not nice either

#count how many novel types in novel associated hosts
novel_in_novel_hosts <- addtl_tip_plus_new_manual %>% 
  group_by(specific_species) %>% 
  mutate(status_count = !("ncbi" %in% status)) %>% filter(status_count == TRUE) %>% filter(grepl("[SED][R]{2}",V1))

novel_in_novel_hosts$library <- gsub("Score.*", "", novel_in_novel_hosts$V1)
novel_in_novel_hosts$library <- gsub("_.*", "", novel_in_novel_hosts$library)

write.table(novel_in_novel_hosts, "outputs/novel_in_novel_hosts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

#get library for all the novel ones
novel_types <- addtl_tip_plus_new_manual %>% 
filter(grepl("[SED][R]{2}",V1))

novel_types$library <- gsub("Score.*", "", novel_types$V1)
novel_types$library <- gsub("_.*", "", novel_types$library)

novel_headers_type <- left_join(novel_types, sra_metadata, by = c("library" = "V1"))

write.table(novel_headers_type$V1, "outputs/all_novel_libraries.tsv", sep = "\t", col.names = F, row.names = F, quote = F)


type_table <- as.data.frame(table(novel_headers_type$V2))
type_table = type_table[-1,]

#use all_pvs mapping from map_gen_for_pub.R
all_pvs_mapping$library <- gsub("Score.*", "", all_pvs_mapping$V1.x)
all_pvs_mapping$library <- gsub("_.*", "", all_pvs_mapping$library)


novel_headers_type <- left_join(all_pvs_mapping, sra_metadata, by = c("library" = "V1"))
novel_headers_type$collapsed_type <-coalesce(novel_headers_type$LibrarySource, novel_headers_type$V15)

novel_headers_type_missing_type <- 
  novel_headers_type %>% 
  filter(is.na(collapsed_type))

# Get library source for SRA accessions missing type info using rentrez/efetch
library(rentrez)

missing_libraries <- as.data.frame(unique(novel_headers_type_missing_type$library))
colnames(missing_libraries) <- c("V1")
missing_libraries <- missing_libraries %>% filter(!grepl('*SRR2707*', V1))
#manuallyfix parsing errors
missing_libraries$V1 <- gsub("SRR23971524SRR23971524", "SRR23971524", missing_libraries$V1)
missing_libraries$V1 <- gsub("SRR24184635SRR24184635", "SRR24184635", missing_libraries$V1)


results <- data.frame(library = character(), librarysource = character(), stringsAsFactors = FALSE)

for (acc in missing_libraries$V1) {
  search_results <- entrez_search(db = "sra", term = paste0(acc, "[ACCN]"))
  
  if (search_results$count > 0) {
    fetch_results <- entrez_fetch(db = "sra", id = search_results$ids, rettype = "xml")
    xml_doc <- xml2::read_xml(fetch_results)
    lib_source <- xml2::xml_text(xml2::xml_find_first(xml_doc, ".//LIBRARY_SOURCE"))
    results <- rbind(results, data.frame(library = acc, librarysource = lib_source, stringsAsFactors = FALSE))
  }
  
  Sys.sleep(0.4)
}
 
results


novel_headers_type <- left_join(novel_headers_type, results, by = c("library" = "library"))

novel_headers_type$collapsed_type_2 <- coalesce(novel_headers_type$collapsed_type, novel_headers_type$librarysource)

# write.table(novel_headers_type, "type_missing_manual.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
# novel_headers_type_MANUAL <- read.table("files/2025.11.28.type_missing_manual.tsv.txt", sep = "\t", header = T)

#okay deal with missing ENA stuff
# --- ENA library source lookup for ERR/ERA accessions missing from NCBI ---
library(httr)
library(xml2)

# Find accessions still missing collapsed_type_2 (NA) or that are ERR/ERA prefixed
ena_missing <- novel_headers_type %>%
  filter(is.na(collapsed_type_2) | collapsed_type_2 == "") %>%  pull(library) %>%
  unique()

# Also catch any ERR accessions that the rentrez loop silently skipped
ena_accs <- missing_libraries$V1[grepl("^(ERR|ERA|DRR)", missing_libraries$V1)]
ena_missing <- unique(c(ena_missing, ena_accs))
ena_missing <- ena_missing[!is.na(ena_missing) & nchar(ena_missing) > 0]

cat("Querying ENA for", length(ena_missing), "accessions\n")

ena_results <- data.frame(library = character(), librarysource = character(), stringsAsFactors = FALSE)

for (acc in ena_missing) {
  tryCatch({
    url <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", acc)
    resp <- httr::GET(url)
    
    if (httr::status_code(resp) == 200) {
      xml_doc <- xml2::read_xml(httr::content(resp, as = "text", encoding = "UTF-8"))
      
      # Try to find LIBRARY_SOURCE in the run XML
      lib_source <- xml2::xml_text(xml2::xml_find_first(xml_doc, ".//LIBRARY_SOURCE"))
      
      # If not in run XML, follow the experiment reference
      if (is.na(lib_source)) {
        exp_acc <- xml2::xml_attr(xml2::xml_find_first(xml_doc, ".//EXPERIMENT_REF"), "accession")
        if (!is.na(exp_acc)) {
          exp_url <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", exp_acc)
          exp_resp <- httr::GET(exp_url)
          if (httr::status_code(exp_resp) == 200) {
            exp_xml <- xml2::read_xml(httr::content(exp_resp, as = "text", encoding = "UTF-8"))
            lib_source <- xml2::xml_text(xml2::xml_find_first(exp_xml, ".//LIBRARY_SOURCE"))
          }
        }
      }
      
      if (!is.na(lib_source)) {
        ena_results <- rbind(ena_results, data.frame(library = acc, librarysource = lib_source, stringsAsFactors = FALSE))
        cat("  ", acc, "->", lib_source, "\n")
      } else {
        cat("  ", acc, "-> LIBRARY_SOURCE not found in XML\n")
      }
    } else {
      cat("  ", acc, "-> HTTP", httr::status_code(resp), "\n")
    }
    
    Sys.sleep(0.3)  # rate limit
  }, error = function(e) {
    cat("  ", acc, "-> Error:", e$message, "\n")
  })
}

# Merge ENA results into the existing results and re-coalesce
results <- rbind(results, ena_results)
results <- results[!duplicated(results$library), ]

# Re-join and fill collapsed_type_2
novel_headers_type <- novel_headers_type %>%
  left_join(results, by = c("library" = "library"))

novel_headers_type <- novel_headers_type %>%
  mutate(collapsed_type_2 = replace_na(collapsed_type_2, "UNKNOWN")) %>%
  group_by(library) %>%
  fill(collapsed_type_2, .direction = "downup") %>%
  ungroup()


novel_headers_type_MANUAL_table <- as.data.frame(table(novel_headers_type$collapsed_type_2, novel_headers_type$status))

library(extrafont)

#clean up metadata
`%notin%` <- Negate(`%in%`)

status <- c("METAGENOMIC", "GENOMIC", "TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "GENOMIC SINGLE CELL", "TRANSCRIPTOMIC SINGLE CELL")

i <- novel_headers_type_MANUAL_table$Var1 %notin% status
novel_headers_type_MANUAL_table$Var1[i] <- "OTHER"

novel_headers_type_MANUAL_table <- unique(novel_headers_type_MANUAL_table %>% group_by(Var1, Var2) %>% mutate(Freq = sum(Freq)))


ggplot(novel_headers_type_MANUAL_table, aes(x=reorder(Var1, Freq), y=Freq, fill = fct_rev(Var2))) + 
  geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = c("#4646d4", "grey")) +
  coord_flip()+ 
  theme(text = element_text(family = "Noto Sans", size = 25))+ 
  theme(axis.text.y = element_text(hjust=1))

ggsave("outputs/2026.02.17.breakdown_new_small.jpg", width = 30, height = 15, units = "cm", limitsize = F)


biosamp_data <- read.table("files/all_runs.csv", sep = ",", header = T, fill = T)

biosamp_novel_headers <- left_join(novel_headers_type, biosamp_data, by = c("library" = "Run"))

write.table(biosamp_novel_headers$BioSample, "outputs/all_novel_biosamples.tsv", sep = "\t", col.names = F, row.names = F, quote = F)


################################END LIB TYPE############################

addtl_ano <- left_join(addtl_ano, sra_metadata, by = c("lib" = "V1"))
addtl_ano$generalization <- iconv(addtl_ano$generalization, from = "UTF-8", to = "ASCII", sub = "")
addtl_ano <- filter(addtl_ano, isTip == "TRUE")


################################ TREE ANNOTATION & METADATA ################################

#manually add version of gheatmap to resolve ggplot conflicts 
gheatmap <- function(p, data, offset=0, width=1, low="green", high="red", color="white",
                     colnames=TRUE, colnames_position="bottom", colnames_angle=0, colnames_level=NULL,
                     colnames_offset_x = 0, colnames_offset_y = 0, font.size=4, family="",
                     hjust=0.5, legend_title = "value", custom_column_labels = NULL) {
  
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  
  ## if (is.null(width)) {
  ##     width <- (p$data$x %>% range %>% diff)/30
  ## }
  
  ## convert width to width of each cell
  width <- width * (p$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)
  
  isTip <- x <- from <- to <- custom_labels <- NULL
  
  ## handle the display of heatmap on collapsed nodes
  ## https://github.com/GuangchuangYu/ggtree/issues/242
  ## extract data on leaves (& on collapsed internal nodes) 
  ## (the latter is extracted only when the input data has data on collapsed
  ## internal nodes)
  df <- p$data
  nodeCo <- intersect(df %>% dplyr::filter(is.na(x)) %>% 
                        select(.data$parent, .data$node) %>% unlist(), 
                      df %>% dplyr::filter(!is.na(x)) %>% 
                        select(.data$parent, .data$node) %>% unlist())
  labCo <- df %>% dplyr::filter(.data$node %in% nodeCo) %>% 
    select(.data$label) %>% unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- df$label %in% selCo
  
  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm=TRUE) + offset
  
  dd <- as.data.frame(data)
  ## dd$lab <- rownames(dd)
  i <- order(df$y)
  
  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  ## dd <- dd[lab, , drop=FALSE]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  
  
  dd$y <- sort(df$y)
  dd$lab <- lab
  ## dd <- melt(dd, id=c("lab", "y"))
  dd <- gather(dd, variable, value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels=colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels=colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  data_axis <- data.frame(from=dd$variable, to=V2)
  data_axis <- unique(data_axis)
  
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  if (is.null(color)) {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value), width=width, inherit.aes=FALSE)
  } else {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value), width=width, color=color, inherit.aes=FALSE)
  }
  if (is(dd$value,"numeric")) {
    p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
  } else {
    p2 <- p2 + scale_fill_discrete(na.value=NA, name = legend_title) #"white")
  }
  
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    data_axis$y <- y
    data_axis[[".panel"]] <- factor("Tree")
    # if custom column annotations are provided
    if (!is.null(custom_column_labels)) {
      # assess the type of input for the custom column annotation
      # either a vector or a named vector with positions for specific names
      if (is.null(names(custom_column_labels))) {
        if (length(custom_column_labels) > nrow(data_axis)) {
          warning(paste("Input column label vector has more elements than there are columns.",
                        "\n", "Using the first ", nrow(data_axis)," elements as labels", sep=""))
          data_axis[["custom_labels"]] <- as.character(custom_column_labels[1:nrow(data_axis)])
        } else if (length(custom_column_labels) < nrow(data_axis)) {
          warning(paste("Input column label vector has fewer elements than there are columns.",
                        "\n", "Using all available labels, n = ",
                        length(custom_column_labels), sep=""))
          data_axis[["custom_labels"]] <- as.character(c(custom_column_labels,
                                                         rep("", nrow(data_axis) - length(custom_column_labels))))
        } else {
          data_axis[["custom_labels"]] <- custom_column_labels
        }
      } else {
        if (!is.null(colnames_level)) {
          # use the colnames levels if available
          # otherwise use the default order provided by the data frame
          vector_order <- colnames_level
          
        } else {
          vector_order <- as.character(data_axis$from)
        }
        for (elem in custom_column_labels) {
          vector_order[which(vector_order == elem)] = names(which(custom_column_labels == elem))
        }
        data_axis[["custom_labels"]] <- vector_order
      }
      p2 <- p2 + geom_text(data=data_axis, aes(x=to, y = y, label=custom_labels),
                           size=font.size, family=family, inherit.aes = FALSE, angle=colnames_angle,
                           nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
    } else {
      p2 <- p2 + geom_text(data=data_axis, aes(x=to, y = y, label=from), size=font.size, family=family,
                           inherit.aes = FALSE, angle=colnames_angle,
                           nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
    }
  }
  p2 <- p2 + theme(legend.position="right")
  ## p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))
  
  if (!colnames) {
    ## https://github.com/GuangchuangYu/ggtree/issues/204
    p2 <- p2 + scale_y_continuous(expand = c(0,0))
  }
  
  attr(p2, "data_axis") <- data_axis
  return(p2)
}

# row.names(addtl_ano_tree_broad) <- gsub("\\.", "-", row.names(addtl_ano_tree_broad))

tree_plot$data <- tree_plot$data %>%
  mutate(across(where(is.character), str_trim))
# addtl_ano_tree_broad$label <- gsub("__1", "_", addtl_ano_tree_broad$label)

tree_plot$data$label <- gsub("__1", "_", tree_plot$data$label)

select <- dplyr::select

# Read missing_manual once (used for both heatmap and metadata)
missing_manual_full <- read.table("outputs/2026.03.24.missing_manual.tsv.txt", sep = "\t", header = T)

# --- Build all_tree_metadata: one row per tree tip with full annotation ---
combined_metadata <- bind_rows(
  addtl_tip_plus_new_manual %>%
    dplyr::select(label = V1, status, generalization, broader_gen, specific_association),
  missing_manual_full %>%
    mutate(specific_association = coalesce(specific_association, specific_species)) %>%
    dplyr::select(label, status, generalization, broader_gen, specific_association)
) %>%
  mutate(label = gsub("\\.", "_", label)) %>%
  distinct(label, .keep_all = TRUE)

tree_tips <- tree_plot$data %>% filter(isTip == TRUE) %>%
  mutate(label = gsub("\\.", "_", label)) %>%
  dplyr::select(label)
all_tree_metadata <- left_join(tree_tips, combined_metadata, by = "label")
all_tree_metadata$status <- ifelse(grepl("[E|S|D][R][R]", all_tree_metadata$label), "novel", "ncbi")

# --- Build addtl_ano_tree_broad: heatmap annotation (rownames = tip labels) ---
addtl_ano_tree_broad <- dplyr::select(addtl_tip_plus_new_manual, V1, broader_gen)
rownames(addtl_ano_tree_broad) <- gsub("\\.", "_", addtl_ano_tree_broad$V1)
addtl_ano_tree_broad <- dplyr::select(addtl_ano_tree_broad, broader_gen)

# Add missing sequences
missing_heatmap <- select(missing_manual_full, label, broader_gen) %>%
  remove_rownames() %>% column_to_rownames(var = "label")
addtl_ano_tree_broad <- rbind(addtl_ano_tree_broad, missing_heatmap)

# Collapse generalization groups with < 10 representatives into "Other",
# and promote large groups currently labelled "Other" to their generalization name
gen_counts <- all_tree_metadata %>% count(generalization)
small_groups <- gen_counts %>% filter(n < 10) %>% pull(generalization)
large_groups <- gen_counts %>% filter(n >= 10) %>% pull(generalization)
all_tree_metadata$broader_gen <- case_when(
  all_tree_metadata$generalization %in% small_groups ~ "Other",
  all_tree_metadata$broader_gen == "Other" & all_tree_metadata$generalization %in% large_groups ~ all_tree_metadata$generalization,
  TRUE ~ all_tree_metadata$broader_gen
)

# --- Color key ---
ordered_gen_colors <- c(
  "Human" = "#D44746", "Non-human Primate" = "#F0B2B2",
  "Rodent" = "#D4B046", "Cetacean" = "#175F67",
  "Bovine" = "#8C46D4", "Pinnipeds" = "#E146D4",
  "Cervine" = "#882020", "Canine" = "#44C8D6",
  "Feline" = "#A3E4EB", "Equine" = "#941888",
  "Pangolin" = "#F0E0B2", "Bat" = "#C099E7",
  "Avian" = "#3092D4", "Reptile" = "#C9E349",
  "Amphibian" = "#53741D", "Ray-finned Fish" = "#9BBED1",
  "Other" = "#CCCCCC"
)
active_colors <- ordered_gen_colors[intersect(names(ordered_gen_colors), unique(addtl_ano_tree_broad$broader_gen))]

# --- Bar plots ---
ggplot(all_tree_metadata, aes(x=fct_relevel(fct_infreq(broader_gen), "Other", after = Inf), fill = fct_rev(status))) +
  geom_bar(stat="count") + scale_fill_manual(values = c("#4646d4", "grey")) + theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  labs(title = "Known and Novel PVs\nAssociated with Generalized Host Categories",
       x = "Category",
       y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(text = element_text(family="Noto Sans", size = 30))

ggsave("outputs/2026.03.30.host_group_chart.png", width = 30, height = 20, units = "cm", limitsize = F)

ggplot(filter(all_tree_metadata, broader_gen == "Other"), aes(x=fct_infreq(generalization) , fill =  fct_rev(status))) +
  geom_bar(stat="count") + scale_fill_manual(values = c("#4646d4", "grey")) + theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  labs(title = "Known and Novel PVs Associated with Generalized Host Categories",
       x = "Category",
       y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(text = element_text(family="Noto Sans", size = 30))

ggsave("outputs/2026.03.30.host_group_chart_ZOOM.png", width = 35, height = 20, units = "cm", limitsize = F)

# --- Tree heatmap (using all_tree_metadata) ---
heatmap_df <- all_tree_metadata %>%
  dplyr::select(label, broader_gen) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = "label")

# Rebuild active_colors from all_tree_metadata groups (includes promoted groups)
active_colors <- ordered_gen_colors[intersect(names(ordered_gen_colors), unique(heatmap_df$broader_gen))]
# Add colors for any new groups not in ordered_gen_colors (promoted from generalization)
new_groups <- setdiff(unique(heatmap_df$broader_gen), names(ordered_gen_colors))
if (length(new_groups) > 0) {
  new_colors <- setNames(rep("#ffb44a", length(new_groups)), new_groups)
  active_colors <- c(active_colors, new_colors)
}


gheatmap(tree_plot, heatmap_df, width=0.05, font.size=0,  offset = 0.4,  color = NA) +
  scale_fill_manual(breaks=names(active_colors),
                    values=active_colors, name="host_group") +
                    theme(legend.position = 'right')


ggsave("outputs/2026.02.23.fixed_new_broad_colors_medium_FINAL.png", width = 50, height =75, units = "cm", limitsize = F)
ggsave("outputs/2026.03.24.fixed_new_broad_colors_LARGE_FINAL.pdf", width = 50, height =200, units = "cm", limitsize = F)
ggsave("outputs/2026.02.23.fixed_new_broad_colors_LARGE_final.png", width = 50, height =200, units = "cm", limitsize = F)

tree_plot_joined <- left_join(tree_plot$data, rownames_to_column(addtl_ano_tree_broad), by = c("label" = "rowname"))
tree_plot_joined <- filter(tree_plot_joined, isTip == TRUE)
per_host_group <- as.data.frame(table(tree_plot_joined$broader_gen, tree_plot_joined$status))

per_host_group_wide <- pivot_wider(
  per_host_group,
  names_from = Var2,
  values_from = Freq
)
write.table(per_host_group_wide, "outputs/per_host_group_wide.tsv", sep = "\t", quote = F, row.names = F)
#generate circular tree for figure

tree_unroot <- ggtree(final_tree, layout = "circular")
  
tree_unroot$data$label <- gsub("Edges.*", "", tree_unroot$data$label)
tree_unroot$data$label <- gsub("ka_f.*", "", tree_unroot$data$label)
  

# tree_unroot$data$x10 <- iconv(tree_unroot$data$x10, from = "UTF-8", to = "ASCII", sub = "")
# tree_unroot$data[1092, 10] <- "Zhangixalus_dugritei_associated"

tree_unroot$data$status <- ifelse(grepl("[E|S|D][R][R]", tree_unroot$data$label), "novel", "ncbi")

  
tree_unroot$data$color_values <- ifelse(grepl("ncbi", tree_unroot$data$status), "grey", "#3333E7")
  
circular_tree_plot <- tree_unroot %<+% select(tree_unroot$data, label, status) + aes(colour=I(tree_unroot$data$color_values))+
    # geom_tippoint() +  # moved to p3 as layered novel-on-top

    # geom_tiplab(color = "black",size = 2, hjust = 0, linetype = "blank") +
    # geom_text(aes(label=label), hjust=-0.5, size = 2, color = "black") +
    # geom_nodelab(label = ncbi_and_novel_tree$node.label, geom = 'text', size = 1.5) +
    theme(legend.position= "none")

    circular_tree_plot


circular_tree_plot$data <- circular_tree_plot$data %>%
  mutate(across(where(is.character), str_trim))
# addtl_ano_tree_broad$label <- gsub("__1", "_", addtl_ano_tree_broad$label)

circular_tree_plot$data$label <- gsub("__1", "_", circular_tree_plot$data$label)

  p1 <- gheatmap(circular_tree_plot, addtl_ano_tree_broad, offset=-0.5, width=0.2,
           colnames_angle=-45, color = NA, font.size=0) + 
    scale_fill_manual(breaks=addtl_ano_tree_broad$broader_gen, 
                      values=addtl_ano_tree_broad_vals$color_values, 
                      name="status") + 
    theme(legend.position = 'none' )

  p1
  
  
  p2 <- p1 + new_scale_fill()

addtl_ano_tree_broad_stats <- addtl_ano_tree_broad
addtl_ano_tree_broad_stats$status <- ifelse(grepl("[ESD]RR", rownames(addtl_ano_tree_broad_stats)), "novel", "ncbi")
addtl_ano_tree_broad_stats <- select(addtl_ano_tree_broad_stats, status)

  # Tip points: grey (ncbi) drawn first, blue (novel) on top
  p3 <- p2 +
    geom_tippoint(data = . %>% filter(isTip, !grepl("[ESD]RR", label)),
                  color = "grey", size = 0.6) +
    geom_tippoint(data = . %>% filter(isTip, grepl("[ESD]RR", label)),
                  color = "#3333E7", size = 0.6) +
    theme(legend.position = 'none')

  p3

ggsave("outputs/2026.03.24.circular_tree.png", width = 15, height =15, units = "cm", limitsize = F)





check_tree <- tree_plot$data
check_tree <- filter(check_tree, isTip == T)

ggplot(all_tree_metadata, aes(x=fct_rev(fct_relevel(fct_infreq(broader_gen), "Other", after = Inf)), fill = fct_rev(status))) +
  geom_bar(stat="count") + scale_fill_manual(values = c("#4646d4", "grey")) + theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), panel.grid.minor.y = element_line(colour = "grey95")) +
  labs(title = "Known and Novel PVs\nAssociated with Generalized Host Categories",
       x = "Category",
       y = "Count") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(text = element_text(family="Noto Sans", size = 30)) +
      coord_flip()

ggsave("outputs/2026.03.24.poster_host.png", width = 30, height =20, units = "cm", limitsize = F)


##test 'ratios'
# ---- Fisher's exact test: Human proportion in Logan vs NCBI ----

# Set 1: Logan clusters (~1100)
logan_set <- all_clusters_summary %>%
  filter(centroid_library != "SRR6976994") %>%
  transmute(source = "Logan", host_group = as.character(host_group))

add <- c("Logan", "Human")
logan_set <- rbind(logan_set, add)
# Set 2: NCBI reference PVs from tree (~992)
ncbi_set <- all_tree_metadata %>%
  filter(status == "ncbi") %>%
  transmute(source = "NCBI", host_group = broader_gen)

cat("Logan clusters:", nrow(logan_set), "\n")
cat("NCBI sequences:", nrow(ncbi_set), "\n")

# Combine and build 2x2 table
combined <- bind_rows(logan_set, ncbi_set) %>%
  mutate(is_human = if_else(host_group == "Human", "Human", "Non-Human"))

contingency <- table(combined$source, combined$is_human)
print(contingency)
fisher.test(contingency)

combined_clean <- combined %>% filter(!is.na(host_group))

fisher_results <- map_dfr(hosts_in_both, function(hg) {
  combined_hg <- combined_clean %>%
    mutate(is_target = if_else(host_group == hg, "Yes", "No"))
  ct <- table(combined_hg$source, combined_hg$is_target)
  if (nrow(ct) < 2 || ncol(ct) < 2) return(NULL)
  ft <- fisher.test(ct)
  tibble(
    host_group = hg,
    logan_yes  = ct["Logan", "Yes"],
    logan_no   = ct["Logan", "No"],
    ncbi_yes   = ct["NCBI", "Yes"],
    ncbi_no    = ct["NCBI", "No"],
    p_value    = ft$p.value,
    odds_ratio = ft$estimate
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

print(fisher_results, n = nrow(fisher_results))

write.table(fisher_results, "outputs/2026.04.02.fisher_results.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

