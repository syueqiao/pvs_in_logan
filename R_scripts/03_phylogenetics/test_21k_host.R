library(tidyverse)
###look at all 21k
all_21323_libraries_metadata <- read.table("files/all_21323_libraries_metadata.list", fill = TRUE, quote = "", sep = ",")

# Classify V3 (organism/host) into generalized host groups
all_21323_libraries_metadata <- all_21323_libraries_metadata %>%
  mutate(V3_clean = str_trim(gsub('"', '', V3))) %>%
  mutate(host_group = case_when(
    # Human
    str_detect(V3_clean, regex("^Homo sapiens|^human |^Human p|human feces|human gut|human skin|human vaginal|human oral|human nasopharyngeal|human reproductive|human metagenome|human viral|human lung|human urinary|human saliva|human blood|human milk|uncultured human|^vaginal metagenome", ignore_case = TRUE)) ~ "Human",
    str_detect(V3_clean, regex("^Human rhinovirus|^Human respiratory|^Human alpha|^Human gamma|human papillomavirus", ignore_case = TRUE)) ~ "Human",

    # Bat
    str_detect(V3_clean, regex("^bat metagenome|^Chiroptera|^Desmodus rotundus|^Eidolon helvum|^Monophyllus redmani", ignore_case = TRUE)) ~ "Bat",

    # Bovine
    str_detect(V3_clean, regex("^Bos taurus|^Bos indicus|^Bos grunniens|^bovine|^Bubalus bubalis|^Syncerus caffer|^Budorcas taxicolor|^milk metagenome", ignore_case = TRUE)) ~ "Bovine",

    # Non-human Primate
    str_detect(V3_clean, regex("^Macaca|^Gorilla|^Saimiri|^Cebus|^Prolemur|^Perodicticus|^Tarsius|^primate metagenome", ignore_case = TRUE)) ~ "Non-human Primate",

    # Avian
    str_detect(V3_clean, regex("^Gallus|^Nipponia|^Thamnophilus|^Hieraaetus|^Geranoaetus|^Falco|^Aquila|^Accipiter|^Circaetus|^Clanga|^Anser|^Alipiopsitta|^Dendrortyx|^Conuropsis|^bird metagenome", ignore_case = TRUE)) ~ "Avian",

    # Cetacean
    str_detect(V3_clean, regex("^Sousa chinensis", ignore_case = TRUE)) ~ "Cetacean",

    # Ray-finned Fish
    str_detect(V3_clean, regex("^Gadus morhua|^Cynoglossus|^Coilia nasus|^Paralichthys|^Oncorhynchus|^Cyclopterus|^Danio rerio|^Dicentrarchus|^Hypophthalmichthys|^Carassius|^Coryphaenoides|^Misgurnus|^Malapterurus|^Astyanax|^Triplophysa|^Tachysurus|^Chionodraco|^Parachondrostoma|Ray-finned", ignore_case = TRUE)) ~ "Ray-finned Fish",

    # Cervine (deer)
    str_detect(V3_clean, regex("^Rangifer", ignore_case = TRUE)) ~ "Cervine",

    # Equine
    str_detect(V3_clean, regex("^Equus caballus|^horse metagenome", ignore_case = TRUE)) ~ "Equine",

    # Rodent
    str_detect(V3_clean, regex("^Mus musculus|^Peromyscus|^Rattus|^Mesocricetus|^Myodes|^Arvicanthis|^Cavia porcellus|^Marmota|^rodent metagenome|^mouse gut metagenome|^rat gut metagenome", ignore_case = TRUE)) ~ "Rodent",

    # Feline
    str_detect(V3_clean, regex("^Felis catus|^Felis silvestris|^Panthera tigris", ignore_case = TRUE)) ~ "Feline",

    # Insectivore (Eulipotyphla: shrews, moles, solenodons, hedgehogs-adjacent)
    str_detect(V3_clean, regex("^Sorex|^Crocidura|^Suncus|^Talpa|^Condylura|^Scalopus|^Solenodon|^Tenrec|^Hemiechinus|^Atelerix|eulipotyphla", ignore_case = TRUE)) ~ "Insectivore",

    # Canine
    str_detect(V3_clean, regex("^Canis lupus|^Canis latrans|^Canis familiaris|^Vulpes|^Urocyon|^canine metagenome", ignore_case = TRUE)) ~ "Canine",

    # Pinniped
    str_detect(V3_clean, regex("^Arctocephalus|Phocidae|phocicerebrale|phocidae", ignore_case = TRUE)) ~ "Pinniped",

    # Amphibian
    str_detect(V3_clean, regex("^Ambystoma|^Amolops|^Ptychadena|^Leptobrachium|^Lithobates|^Bufo|^Rana |^Zhangixalus|^Phyllomedusa|^Oreolalax|^Oophaga|^Limnonectes|^Desmognathus|^Triturus", ignore_case = TRUE)) ~ "Amphibian",

    # Hedgehog (Erinaceidae) - none apparent in data, but keep category
    str_detect(V3_clean, regex("^Erinaceus|hedgehog", ignore_case = TRUE)) ~ "Hedgehog",

    # Mustelid
    str_detect(V3_clean, regex("^Mustela|^Gulo gulo", ignore_case = TRUE)) ~ "Mustelid",

    # Reptile
    str_detect(V3_clean, regex("^Phrynocephalus|^Takydromus|^Podarcis|^Dibamus", ignore_case = TRUE)) ~ "Reptile",

    # Camelid
    str_detect(V3_clean, regex("^Camelus|^Lama|^Vicugna|camelid", ignore_case = TRUE)) ~ "Camelid",

    # Caprine
    str_detect(V3_clean, regex("^Capra|^caprine", ignore_case = TRUE)) ~ "Caprine",

    # Manatee (Sirenia)
    str_detect(V3_clean, regex("^Trichechus|manatee", ignore_case = TRUE)) ~ "Manatee",

    # Marsupial
    str_detect(V3_clean, regex("^Antechinus|^Sarcophilus", ignore_case = TRUE)) ~ "Marsupial",

    # Ovine
    str_detect(V3_clean, regex("^Ovis aries", ignore_case = TRUE)) ~ "Ovine",

    # Giraffid
    str_detect(V3_clean, regex("^Giraffa", ignore_case = TRUE)) ~ "Giraffid",

    # Ursine
    str_detect(V3_clean, regex("^Ursus|^Ailuropoda", ignore_case = TRUE)) ~ "Ursine",

    # Hyena
    str_detect(V3_clean, regex("^Crocuta|hyena|Hyaenidae", ignore_case = TRUE)) ~ "Hyena",

    # Lapine
    str_detect(V3_clean, regex("^Oryctolagus", ignore_case = TRUE)) ~ "Lapine",

    # Pangolin
    str_detect(V3_clean, regex("^Manis|^Manidae", ignore_case = TRUE)) ~ "Pangolin",

    # Porcine
    str_detect(V3_clean, regex("^Sus scrofa|^Phacochoerus|^Potamochoerus", ignore_case = TRUE)) ~ "Porcine",

    # Raccoon (Procyonidae)
    str_detect(V3_clean, regex("^Procyon|raccoon", ignore_case = TRUE)) ~ "Raccoon",

    # PV species named by host association
    str_detect(V3_clean, regex("^Deltapapillomavirus", ignore_case = TRUE)) ~ "Bovine",
    str_detect(V3_clean, regex("^Macaca.*papillomavirus", ignore_case = TRUE)) ~ "Non-human Primate",
    str_detect(V3_clean, regex("^Canis.*papillomavirus", ignore_case = TRUE)) ~ "Canine",
    str_detect(V3_clean, regex("^Ceratotherium", ignore_case = TRUE)) ~ "Other",

    # Everything else -> Other
    TRUE ~ "Other"
  )) %>%
  mutate(host_group = factor(host_group, levels = c(
    "Human", "Bat", "Bovine", "Non-human Primate", "Avian", "Cetacean",
    "Ray-finned Fish", "Cervine", "Equine", "Rodent", "Feline", "Insectivore",
    "Canine", "Pinniped", "Amphibian", "Hedgehog", "Mustelid", "Reptile",
    "Camelid", "Caprine", "Manatee", "Marsupial", "Ovine", "Giraffid",
    "Other", "Ursine", "Hyena", "Lapine", "Pangolin", "Porcine", "Raccoon"
  ))) %>%
  # Confidence ranking: high = exact species match, medium = host-associated metagenome/PV name,
  # low = ambiguous metagenome or organism with unclear host relevance, flag = needs manual inspection
  mutate(host_confidence = case_when(
    # HIGH: direct species binomial -> unambiguous host
    str_detect(V3_clean, regex("^Homo sapiens$|^Bos taurus$|^Bos indicus|^Bos grunniens|^Bubalus bubalis|^Syncerus caffer|^Budorcas taxicolor", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Mus musculus$|^Felis catus$|^Felis silvestris|^Canis lupus|^Canis latrans|^Equus caballus", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Ovis aries|^Capra |^Sus scrofa|^Oryctolagus|^Giraffa|^Ursus |^Ailuropoda", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Rattus |^Rattus$|^Peromyscus|^Mesocricetus|^Cavia porcellus|^Marmota|^Myodes|^Arvicanthis", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Macaca [a-z]|^Gorilla|^Saimiri|^Cebus|^Prolemur|^Perodicticus|^Tarsius", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Desmodus rotundus|^Eidolon helvum|^Monophyllus", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Gadus morhua|^Danio rerio|^Oncorhynchus|^Paralichthys|^Dicentrarchus|^Cyclopterus", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Cynoglossus|^Coilia|^Carassius|^Hypophthalmichthys|^Misgurnus|^Malapterurus|^Astyanax", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Triplophysa|^Tachysurus|^Chionodraco|^Parachondrostoma|^Coryphaenoides", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Rangifer|^Sousa chinensis|^Arctocephalus|^Panthera tigris|^Vulpes", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Mustela|^Gulo gulo|^Phacochoerus|^Potamochoerus|^Manis |^Manidae", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Gallus|^Falco |^Anser |^Accipiter|^Aquila|^Nipponia", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Ambystoma|^Amolops|^Ptychadena|^Lithobates|^Bufo |^Zhangixalus|^Rana ", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Leptobrachium|^Phyllomedusa|^Oreolalax|^Oophaga|^Limnonectes|^Desmognathus|^Triturus", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Phrynocephalus|^Takydromus|^Podarcis|^Dibamus", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Antechinus|^Sarcophilus|^Urocyon|^Canis familiaris$", ignore_case = TRUE)) ~ "high",
    str_detect(V3_clean, regex("^Spilogale|^Ceratotherium", ignore_case = TRUE)) ~ "high",

    # MEDIUM: host-specific metagenomes or PV species named by host
    str_detect(V3_clean, regex("^human |human feces|human gut|human skin|human vaginal|human oral|human nasopha|human reproductive|human viral|human lung|human urinary|human saliva|human blood|human milk|uncultured human|vaginal metagenome", ignore_case = TRUE)) ~ "medium",
    str_detect(V3_clean, regex("^bat metagenome|^bovine |^canine metagenome|^horse metagenome|^rodent metagenome|^bird metagenome|^primate metagenome", ignore_case = TRUE)) ~ "medium",
    str_detect(V3_clean, regex("^mouse gut metagenome|^rat gut metagenome|^bovine gut metagenome|^milk metagenome", ignore_case = TRUE)) ~ "medium",
    str_detect(V3_clean, regex("papillomavirus", ignore_case = TRUE)) ~ "medium",
    str_detect(V3_clean, regex("^Human rhinovirus|^Human respiratory|^Human alpha|^Human gamma", ignore_case = TRUE)) ~ "medium",
    str_detect(V3_clean, regex("^Chiroptera$|^Deltapapillomavirus|Phocidae|phocicerebrale", ignore_case = TRUE)) ~ "medium",

    # LOW: ambiguous metagenomes where host is unclear - flag for manual inspection
    str_detect(V3_clean, regex("^metagenome$|^viral metagenome|^skin metagenome$|^vaginal metagenome$|^gut metagenome$|^eye metagenome", ignore_case = TRUE)) ~ "low",
    str_detect(V3_clean, regex("^soil metagenome|^air metagenome|^dust metagenome|^wastewater|^marine |^cold seep|^indoor|^whole organism", ignore_case = TRUE)) ~ "low",
    str_detect(V3_clean, regex("^sediment|^aquatic|^food |^hospital|^urban|^brine|^riverine|^groundwater|^oil production|^respiratory tract metagenome", ignore_case = TRUE)) ~ "low",
    str_detect(V3_clean, regex("^synthetic|^medical device|^museum specimen|^fossil|^organismal|^tick metagenome|^insect gut|^ciliate|^algae|^urine metagenome", ignore_case = TRUE)) ~ "low",
    str_detect(V3_clean, regex("^mixed sample|^unidentified|^Papillomavirus sp|^Papillomaviridae", ignore_case = TRUE)) ~ "low",

    # LOW: non-vertebrate organisms (plants, fungi, bacteria, invertebrates, protists) - unlikely true PV hosts
    TRUE ~ "low"
  )) %>%
  mutate(needs_review = host_confidence == "low")

# --- BioSample lookup for ambiguous metagenome entries ---
# Load the larger SRA metadata that has BioSample accessions
sra_full <- read.table("files/all_hits_info.list", sep = ",", fill = TRUE,
                        header = TRUE, quote = "", comment.char = "#")

# Identify rows needing BioSample resolution
ambiguous_rows <- all_21323_libraries_metadata %>%
  filter(needs_review == TRUE,
         str_detect(V3_clean, regex("metagenome|^viral metagenome|^skin metagenome$|^vaginal metagenome$|^gut metagenome$|^eye metagenome$|^mixed sample|^unidentified", ignore_case = TRUE)))

# Get BioSample accessions by joining with the full SRA metadata
# V1 has quotes from read.table; strip them to match Run in sra_full
ambiguous_with_biosample <- ambiguous_rows %>%
  mutate(V1_clean = gsub('"', '', V1)) %>%
  left_join(sra_full %>% dplyr::select(Run, BioSample) %>% distinct(),
            by = c("V1_clean" = "Run")) %>%
  filter(!is.na(BioSample), BioSample != "")

biosample_ids <- unique(ambiguous_with_biosample$BioSample)
cat("Querying NCBI BioSample for", length(biosample_ids), "unique BioSamples...\n")

# Batch query NCBI BioSample for host and isolation_source attributes
# rentrez efetch in batches of 200
library(rentrez)
library(xml2)

batch_size <- 200
biosample_results <- data.frame(BioSample = character(), host = character(),
                                 isolation_source = character(), organism = character(),
                                 title = character(), stringsAsFactors = FALSE)

if (length(biosample_ids) == 0) {
  cat("No BioSample IDs found - check that V1 values match Run column in all_hits_info.list\n")
}

# Step 1: Convert BioSample accessions to NCBI numeric UIDs via batched OR queries
# One API call per batch instead of one per accession - avoids hanging on individual timeouts
# Accession->BioSample mapping is recovered from XML in Step 2 (acc_returned from xml_attr)
cat("Converting BioSample accessions to NCBI UIDs...\n")
acc_to_uid <- data.frame(BioSample = character(), uid = character(), stringsAsFactors = FALSE)

search_batch_size <- 50  # smaller batches keep OR queries short and reliable
for (i in seq_len(ceiling(length(biosample_ids) / search_batch_size))) {
  start_idx <- (i - 1) * search_batch_size + 1
  batch <- biosample_ids[start_idx:min(start_idx + search_batch_size - 1, length(biosample_ids))]

  term <- paste0('"', batch, '"[Accession]', collapse = " OR ")

  tryCatch({
    search_res <- entrez_search(db = "biosample", term = term, retmax = length(batch))
    if (length(search_res$ids) > 0) {
      acc_to_uid <- rbind(acc_to_uid,
                          data.frame(BioSample = NA_character_, uid = search_res$ids,
                                     stringsAsFactors = FALSE))
    }
  }, error = function(e) {
    cat("  Warning: search batch", i, "failed -", conditionMessage(e), "\n")
  })
  Sys.sleep(0.34)  # respect ~3 req/sec without API key
  if (i %% 10 == 0) cat("  Searched", min(start_idx + search_batch_size - 1, length(biosample_ids)),
                          "of", length(biosample_ids), "accessions\n")
}

cat("Retrieved", nrow(acc_to_uid), "UIDs for", length(biosample_ids), "accessions\n")

# Step 2: Fetch BioSample XML in batches using numeric UIDs
uid_list <- unique(acc_to_uid$uid)
cat("Fetching BioSample details for", length(uid_list), "UIDs...\n")

for (i in seq_len(ceiling(length(uid_list) / batch_size))) {
  start_idx <- (i - 1) * batch_size + 1
  batch_uids <- uid_list[start_idx:min(start_idx + batch_size - 1, length(uid_list))]
  cat("  Fetching batch", i, "of", ceiling(length(uid_list) / batch_size), "...\n")

  tryCatch({
    fetch_result <- entrez_fetch(db = "biosample", id = batch_uids, rettype = "xml")
    xml_doc <- read_xml(fetch_result)
    samples <- xml_find_all(xml_doc, "//BioSample")

    for (s in samples) {
      uid_returned <- xml_attr(s, "id")  # numeric UID
      acc_returned <- xml_attr(s, "accession")  # SAMN accession

      # Map back to original queried accession(s) via UID
      # Note: BioSample column may be NA when using batched OR queries in Step 1
      original_accs <- acc_to_uid$BioSample[acc_to_uid$uid == uid_returned]
      if (length(original_accs) == 0 || all(is.na(original_accs))) original_accs <- acc_returned

      # Extract host attribute
      host_node <- xml_find_first(s, ".//Attribute[@attribute_name='host']")
      host_val <- if (!is.na(host_node)) xml_text(host_node) else NA_character_
      # Extract isolation_source attribute
      iso_node <- xml_find_first(s, ".//Attribute[@attribute_name='isolation_source']")
      iso_val <- if (!is.na(iso_node)) xml_text(iso_node) else NA_character_
      # Extract Organism from Description (e.g. <Organism taxonomy_name="Homo sapiens"/>)
      org_node <- xml_find_first(s, ".//Description/Organism")
      org_val <- if (!is.na(org_node)) xml_attr(org_node, "taxonomy_name") else NA_character_
      # Extract Title from Description (sometimes contains host info)
      title_node <- xml_find_first(s, ".//Description/Title")
      title_val <- if (!is.na(title_node)) xml_text(title_node) else NA_character_

      # Store one row per original accession (handles SAMD/SAME both mapping to same UID)
      for (orig_acc in original_accs) {
        biosample_results <- rbind(biosample_results,
                                    data.frame(BioSample = orig_acc, host = host_val,
                                               isolation_source = iso_val,
                                               organism = org_val, title = title_val,
                                               stringsAsFactors = FALSE))
      }
    }
    Sys.sleep(0.4)  # rate limit
  }, error = function(e) {
    cat("  Warning: batch failed -", conditionMessage(e), "\n")
  })
}

cat("Retrieved host info for", nrow(biosample_results), "BioSamples\n")
cat("Non-NA hosts:", sum(!is.na(biosample_results$host)), "\n")
cat("Non-NA isolation_source:", sum(!is.na(biosample_results$isolation_source)), "\n")
cat("Non-NA organism:", sum(!is.na(biosample_results$organism)), "\n")
saveRDS(biosample_results, "outputs/biosample_results_cache.rds")
# To reload without re-running the API: biosample_results <- readRDS("outputs/biosample_results_cache.rds")

biosample_results <- readRDS("outputs/biosample_results_cache.rds")

# Merge biosample host info back to the ambiguous rows
# Priority: host attribute > organism field > isolation_source > title
ambiguous_resolved <- ambiguous_with_biosample %>%
  left_join(biosample_results, by = "BioSample") %>%
  mutate(biosample_organism = coalesce(host, organism, isolation_source, title)) %>%
  filter(!is.na(biosample_organism), biosample_organism != "", biosample_organism != "not collected",
         biosample_organism != "not applicable", biosample_organism != "missing",
         biosample_organism != "metagenome", biosample_organism != "viral metagenome") %>%
  # Reclassify based on biosample host info using the same logic
  mutate(host_group_resolved = case_when(
    str_detect(biosample_organism, regex("homo sapiens|human|vaginal metagenome", ignore_case = TRUE)) ~ "Human",
    str_detect(biosample_organism, regex("chiroptera|bat|Desmodus|Eidolon|Monophyllus|Rhinolophus|Myotis|Pteropus|Rousettus|Miniopterus|Pipistrellus|Carollia|Artibeus|Molossus|Tadarida|Hipposideros|Noctilio|Phyllostomus|Glossophaga", ignore_case = TRUE)) ~ "Bat",
    str_detect(biosample_organism, regex("^Bos |bovine|cow|cattle|buffalo|Bubalus|Syncerus|yak", ignore_case = TRUE)) ~ "Bovine",
    str_detect(biosample_organism, regex("macaca|gorilla|saimiri|cebus|primate|monkey|chimp|Pan troglo|baboon|Papio|lemur|marmoset|Callithrix", ignore_case = TRUE)) ~ "Non-human Primate",
    str_detect(biosample_organism, regex("chicken|gallus|avian|bird|duck|goose|eagle|falcon|hawk|parrot|penguin|pigeon|turkey|Meleagris|Anas |Anser", ignore_case = TRUE)) ~ "Avian",
    str_detect(biosample_organism, regex("cetacea|whale|dolphin|porpoise|Tursiops|Sousa|Balaen", ignore_case = TRUE)) ~ "Cetacean",
    str_detect(biosample_organism, regex("fish|salmon|trout|cod|Gadus|Danio|zebrafish|carp|bass|tilapia|Oncorhynchus|Cyprinus|Salmo ", ignore_case = TRUE)) ~ "Ray-finned Fish",
    str_detect(biosample_organism, regex("deer|cervid|Cervus|Rangifer|Odocoileus|elk|moose|Alces", ignore_case = TRUE)) ~ "Cervine",
    str_detect(biosample_organism, regex("horse|equine|Equus|donkey|mule", ignore_case = TRUE)) ~ "Equine",
    str_detect(biosample_organism, regex("mouse|^Mus |rat|Rattus|hamster|Mesocricetus|rodent|Peromyscus|vole|Myodes|gerbil|guinea pig|Cavia", ignore_case = TRUE)) ~ "Rodent",
    str_detect(biosample_organism, regex("^cat$|feline|Felis|Panthera|lion|tiger|leopard|cheetah", ignore_case = TRUE)) ~ "Feline",
    str_detect(biosample_organism, regex("shrew|mole|Sorex|Crocidura|Suncus|Talpa|Scalopus|Condylura|Solenodon|eulipotyphla", ignore_case = TRUE)) ~ "Insectivore",
    str_detect(biosample_organism, regex("dog|canine|Canis|wolf|fox|Vulpes|coyote|jackal", ignore_case = TRUE)) ~ "Canine",
    str_detect(biosample_organism, regex("seal|sea lion|pinniped|Phoca|Arctocephalus|Zalophus|walrus|Mirounga|Halichoerus", ignore_case = TRUE)) ~ "Pinniped",
    str_detect(biosample_organism, regex("frog|toad|salamander|newt|amphibian|Ambystoma|Rana |Bufo |Xenopus|axolotl|caecilian", ignore_case = TRUE)) ~ "Amphibian",
    str_detect(biosample_organism, regex("hedgehog|Erinaceus|Atelerix", ignore_case = TRUE)) ~ "Hedgehog",
    str_detect(biosample_organism, regex("weasel|ferret|otter|badger|marten|wolverine|Mustela|Gulo|Lutra|Meles", ignore_case = TRUE)) ~ "Mustelid",
    str_detect(biosample_organism, regex("lizard|snake|turtle|tortoise|crocodil|alligator|gecko|iguana|reptil|Podarcis|chameleon", ignore_case = TRUE)) ~ "Reptile",
    str_detect(biosample_organism, regex("camel|llama|alpaca|Camelus|Lama |Vicugna|dromedary", ignore_case = TRUE)) ~ "Camelid",
    str_detect(biosample_organism, regex("goat|Capra ", ignore_case = TRUE)) ~ "Caprine",
    str_detect(biosample_organism, regex("manatee|Trichechus|dugong|sirenia", ignore_case = TRUE)) ~ "Manatee",
    str_detect(biosample_organism, regex("marsupial|kangaroo|koala|opossum|wombat|Antechinus|Sarcophilus|wallaby|Didelphis|Macropus", ignore_case = TRUE)) ~ "Marsupial",
    str_detect(biosample_organism, regex("sheep|ovine|Ovis|lamb", ignore_case = TRUE)) ~ "Ovine",
    str_detect(biosample_organism, regex("giraffe|Giraffa|okapi", ignore_case = TRUE)) ~ "Giraffid",
    str_detect(biosample_organism, regex("bear|Ursus|panda|Ailuropoda", ignore_case = TRUE)) ~ "Ursine",
    str_detect(biosample_organism, regex("hyena|Crocuta|Hyaena", ignore_case = TRUE)) ~ "Hyena",
    str_detect(biosample_organism, regex("rabbit|hare|Oryctolagus|Lepus|lagomorph", ignore_case = TRUE)) ~ "Lapine",
    str_detect(biosample_organism, regex("pangolin|Manis |Manidae", ignore_case = TRUE)) ~ "Pangolin",
    str_detect(biosample_organism, regex("pig|swine|porcine|Sus scrofa|boar|warthog|Phacochoerus|Potamochoerus", ignore_case = TRUE)) ~ "Porcine",
    str_detect(biosample_organism, regex("raccoon|Procyon", ignore_case = TRUE)) ~ "Raccoon",
    TRUE ~ NA_character_
  ))

# Summary of BioSample resolution
cat("\n--- BioSample resolution summary ---\n")
cat("Ambiguous entries queried:", nrow(ambiguous_rows), "\n")
cat("Resolved to a host group:", sum(!is.na(ambiguous_resolved$host_group_resolved)), "\n")
cat("Still unresolved:", sum(is.na(ambiguous_resolved$host_group_resolved)), "\n")
print(table(ambiguous_resolved$host_group_resolved, useNA = "ifany"))

ambiguous_resolved <- as.data.frame(ambiguous_resolved)
test <- c("test")

# Apply resolved host groups back to the main table
resolved_lookup <- ambiguous_resolved %>%
  filter(!is.na(host_group_resolved)) %>%
  dplyr::select(V1_clean, host_group_resolved, biosample_organism) %>%
  distinct()

cat("Resolved lookup has", nrow(resolved_lookup), "entries to merge back\n")

# Drop columns from any previous run to avoid conflicts on re-run
all_21323_libraries_metadata <- all_21323_libraries_metadata %>%
  dplyr::select(-any_of(c("V1_clean", "host_group_resolved", "biosample_organism")))

all_21323_libraries_metadata <- all_21323_libraries_metadata %>%
  mutate(V1_clean = gsub('"', '', V1)) %>%
  left_join(resolved_lookup, by = "V1_clean") %>%
  mutate(
    # Update host_group where BioSample resolved it
    host_group = if_else(!is.na(host_group_resolved),
                         factor(host_group_resolved, levels = levels(host_group)),
                         host_group),
    # Update confidence: biosample-resolved entries become "medium"
    host_confidence = if_else(!is.na(host_group_resolved), "medium", host_confidence),
    # Update needs_review flag
    needs_review = host_confidence == "low"
  )

cat("\n--- Final host_group distribution ---\n")
print(table(all_21323_libraries_metadata$host_group, useNA = "ifany"))
cat("\nRemaining entries needing review:", sum(all_21323_libraries_metadata$needs_review), "\n")

# --- Diagnostic: inspect unresolved BioSample values to guide case_when expansion ---
# These are entries where BioSample was fetched but didn't match any host regex

unresolved_biosample <- ambiguous_resolved %>%
  filter(is.na(host_group_resolved)) %>%
  dplyr::select(V1_clean, BioSample, biosample_organism, host, organism, isolation_source, title) %>%
  mutate(across(everything(), as.character)) %>%  # flatten any list columns
  distinct() %>%
  as.data.frame()

cat("\n--- Top unresolved biosample_organism values (candidates to add to case_when) ---\n")
print(
  sort(table(unresolved_biosample$biosample_organism), decreasing = TRUE) |> head(40)
)

# Entries where BioSample returned nothing useful at all (not in biosample_results)
no_biosample_info <- ambiguous_with_biosample %>%
  anti_join(biosample_results, by = "BioSample") %>%
  dplyr::select(V1_clean, BioSample, V3_clean) %>%
  distinct() %>%
  as.data.frame()

cat("\n--- Entries with no usable BioSample info returned (", nrow(no_biosample_info), "runs) ---\n")
if (nrow(no_biosample_info) > 0) print(head(no_biosample_info, 20))

# --- Combined resolved + still-ambiguous dataframe for inspection ---
ambiguous_final <- all_21323_libraries_metadata %>%
  mutate(V1_clean = gsub('"', '', V1)) %>%
  filter(V1_clean %in% ambiguous_with_biosample$V1_clean) %>%
  dplyr::select(V1_clean, V3_clean, host_group, host_confidence, biosample_organism, needs_review) %>%
  mutate(
    biosample_organism = as.character(biosample_organism),
    resolution_status  = if_else(host_confidence == "medium", "resolved", "unresolved")
  ) %>%
  as.data.frame()

cat("\n--- ambiguous_final summary ---\n")
cat("Resolved via BioSample:", sum(ambiguous_final$resolution_status == "resolved"), "\n")
cat("Still unresolved:      ", sum(ambiguous_final$resolution_status == "unresolved"), "\n")
print(table(ambiguous_final$resolution_status, ambiguous_final$host_group, useNA = "ifany"))

ambiguous_final_low <- filter(ambiguous_final, host_confidence == 'low')
write.table(ambiguous_final_low, "outputs/ambiguous_final_low.tsv", sep = "\t", row.names = F)

#manual annotation phew

ambiguous_final_low_manual <- read.table("outputs/2026.03.06.ambiguous_final_low_manual.tsv.txt",
                                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Build lookup: V1_clean -> manually assigned host group
manual_lookup <- ambiguous_final_low_manual %>%
  filter(!is.na(host_group), host_group != "") %>%
  dplyr::select(V1_clean, host_group) %>%
  rename(manual_add = host_group) %>%
  distinct()

# Apply manual annotations to all_21323_libraries_metadata
all_21323_libraries_metadata <- all_21323_libraries_metadata %>%
  select(-any_of("manual_host")) %>%
  mutate(V1_clean = gsub('"', '', V1)) %>%
  left_join(manual_lookup, by = "V1_clean") %>%
  mutate(
    host_group      = as.character(host_group),
    host_group      = if_else(!is.na(manual_add), manual_add, host_group),
    host_confidence = if_else(!is.na(manual_add), "manual", host_confidence),
    needs_review = host_confidence == "low"
  ) %>%
  dplyr::rename(manual_host = manual_add)  # keep for granular reference in libraries_final

# --- Final comprehensive output dataframe ---
libraries_final <- all_21323_libraries_metadata %>%
  dplyr::select(V1_clean, V3_clean, host_group, host_confidence,
                biosample_organism, manual_host, needs_review) %>%
  mutate(biosample_organism = as.character(biosample_organism)) %>%
  # Join individual BioSample fields for full provenance
  left_join(
    biosample_results %>%
      mutate(across(everything(), as.character)) %>%
      distinct(BioSample, .keep_all = TRUE),
    by = c("V1_clean" = "BioSample")
  ) %>%
  rename(
    original_metadata     = V3_clean,
    biosample_host        = host,
    biosample_organism_field = organism,
    biosample_isolation_source = isolation_source,
    biosample_title       = title
  ) %>%
  as.data.frame()

cat("\n--- Final library annotation summary ---\n")
cat("Total libraries:        ", nrow(libraries_final), "\n")
cat("High confidence:        ", sum(libraries_final$host_confidence == "high",   na.rm = TRUE), "\n")
cat("Medium (BioSample):     ", sum(libraries_final$host_confidence == "medium", na.rm = TRUE), "\n")
cat("Manual annotation:      ", sum(libraries_final$host_confidence == "manual", na.rm = TRUE), "\n")
cat("Low (unresolved):       ", sum(libraries_final$host_confidence == "low",    na.rm = TRUE), "\n")
cat("NA host_group remaining:", sum(is.na(libraries_final$host_group)), "\n")
print(table(libraries_final$host_group, useNA = "ifany"))

libraries_final$host_group[grepl("vaginal metagenome", libraries_final$V3_clean, ignore.case = TRUE)] <- "Human"
libraries_final$host_group[libraries_final$V1_clean %in% c("SRR27078610", "SRR27078611", "SRR27078612", "SRR27078618")] <- "Human"

# --- Host metadata resolution summary ---
# Categorize by SOURCE of host assignment, not confidence level:
#   1) Input file: host_group assigned directly from V3 organism/host metadata
#      (includes both "high" confidence exact species AND original "medium" host-associated metagenomes)
#   2) BioSample API: originally ambiguous (low), resolved via NCBI BioSample lookup
#      (these were re-labeled "medium" during resolution, but originated as "low")
#   3) Manual: remaining ambiguous entries curated by hand after all API searches
# Distinguish input-file "medium" from BioSample-resolved "medium" using biosample_organism:
#   - BioSample-resolved entries have a non-NA biosample_organism value
#   - Input-file "medium" entries have NA biosample_organism
n_from_input <- sum(
  libraries_final$host_confidence == "high" |
  (libraries_final$host_confidence == "medium" & is.na(libraries_final$biosample_organism)) |
  (libraries_final$host_confidence == "low"    & is.na(libraries_final$biosample_organism)),
  na.rm = TRUE
)
n_from_biosample <- sum(
  (libraries_final$host_confidence == "medium" & !is.na(libraries_final$biosample_organism)) |
  (libraries_final$host_confidence == "low"    & !is.na(libraries_final$biosample_organism)),
  na.rm = TRUE
)

n_manual <- sum(libraries_final$host_confidence == "manual", na.rm = TRUE)


# view(libraries_final)
libraries_final$host_group_broad <- libraries_final$host_group

table(libraries_final$host_group)

write.table(libraries_final, "outputs/libraries_final.tsv", sep = '\t', quote = F, row.names = F, col.names = T)
libraries_final <- read.table("outputs/2026.03.17.libraries_final.tsv.txt", sep = '\t', header = T)
###read in clustering information

all_full_l1s_logan_and_pr_90_uc <- read.table("files/all_full_l1s_logan_and_pr_90.uc", sep = "\t")



cat("\n========== Host metadata resolution summary ==========\n")
cat("Total libraries:                              ", nrow(libraries_final), "\n")
cat("1) Assigned from SRA library metadata:        ", n_from_input, sprintf("(%.1f%%)\n", 100 * n_from_input / nrow(libraries_final)))
cat("2) Resolved via BioSample API:                ", n_from_biosample, sprintf("(%.1f%%)\n", 100 * n_from_biosample / nrow(libraries_final)))
cat("3) Manually curated post-API:                 ", n_manual, sprintf("(%.1f%%)\n", 100 * n_manual / nrow(libraries_final)))
cat("=====================================================\n")


all_full_l1s_logan_and_pr_90_uc_key <- dplyr::select(all_full_l1s_logan_and_pr_90_uc, 
                                                    V1, V2, V9, V10) %>% filter(V1 == "H")
colnames(all_full_l1s_logan_and_pr_90_uc_key) <- c("record_type", "cluster_number", "sequence", "centroid")
all_full_l1s_logan_and_pr_90_uc_key$centroid_library <- gsub("Score.*", "", all_full_l1s_logan_and_pr_90_uc_key$centroid)
all_full_l1s_logan_and_pr_90_uc_key$centroid_library <- gsub("_.*", "", all_full_l1s_logan_and_pr_90_uc_key$centroid_library)

all_full_l1s_logan_and_pr_90_uc_key$sequence_library <- gsub("Score.*", "", all_full_l1s_logan_and_pr_90_uc_key$sequence)
all_full_l1s_logan_and_pr_90_uc_key$sequence_library <- gsub("_.*", "", all_full_l1s_logan_and_pr_90_uc_key$sequence_library)

all_full_l1s_logan_and_pr_90_uc_key_centroid_host <- left_join(
  all_full_l1s_logan_and_pr_90_uc_key,
  libraries_final %>% rename_with(~ paste0("centroid_", .), -V1_clean),
  by = c("centroid_library" = "V1_clean")
)

all_full_l1s_logan_and_pr_90_uc_key_centroid_host_sequence_host <- left_join(
  all_full_l1s_logan_and_pr_90_uc_key_centroid_host,
  libraries_final %>% rename_with(~ paste0("sequence_", .), -V1_clean),
  by = c("sequence_library" = "V1_clean")
)

# ---- Expand new_table to all cluster members via library-level UC join ----
# new_table (alignment_ncbi_centroids_and_pr.table) was searched against the same
# sequences as all_full_l1s_logan_and_pr_90.uc -> 100% library-level overlap confirmed
# Library extractor: handles both Score= (Logan) and _NNN_ka_f_ (pathracer) formats
get_lib <- function(x) gsub("Score.*|_[0-9].*", "", x)
# Normalize V5 to UC format: underscores->dots in floats, strip trailing :start-end
# Handles both Score=609_111 -> Score=609.111 and ka_f_7_631 -> ka_f_7.631
normalize_v5 <- function(x) {
  x <- gsub("(Score=[0-9]+)_([0-9]+)", "\\1.\\2", x)       # Score float: 609_111 -> 609.111
  x <- gsub("(ka_f_[0-9]+)_([0-9]+)", "\\1.\\2", x)        # pathracer float: ka_f_7_631 -> ka_f_7.631
  x <- gsub("_([0-9]+)_([0-9]+):.*$", "_:\\1-\\2", x)       # pathracer ORF coords: _57_1235:0-1181 -> _:57-1235
  x <- gsub("(ScaffoldPath=):([0-9]+-[0-9]+)$", "\\1", x)   # Score= entries: strip :33-1166 suffix
  x
}
new_table <- read.table("files/alignment_envs.table")
# Filter new_table with same thresholds as alignments_fil
new_table$qcov <- abs(new_table$V2 - new_table$V3) / new_table$V4
new_table_fil <- new_table %>%
  filter(qcov > 0.7, V10 < 1e-5, V9 >= 90) %>%
  distinct(V1, V5, .keep_all = TRUE) %>%   # keep all (NCBI ref, centroid) pairs
  mutate(centroid_library = get_lib(V5),
         V5_norm = normalize_v5(V5))

# Step 1: Logan sequence (V5) -> cluster number
sequence_cluster_map <- all_full_l1s_logan_and_pr_90_uc %>%
  filter(V1 %in% c("S", "H")) %>%
  transmute(uc_sequence = V9, cluster_number = V2) %>%
  distinct(uc_sequence, .keep_all = TRUE)

# Step 2: cluster number -> all member sequences
cluster_members <- all_full_l1s_logan_and_pr_90_uc %>%
  filter(V1 %in% c("H", "S")) %>%
  transmute(cluster_number = V2, member_sequence = V9,
            member_library = get_lib(V9))

cat("new_table_fil V5_norm matching UC:",
    sum(unique(new_table_fil$V5_norm) %in% sequence_cluster_map$uc_sequence),
    "of", length(unique(new_table_fil$V5_norm)), "\n")

# Step 3: new_table_fil -> cluster number (via V5_norm) -> all members -> member host info
new_table_expanded <- new_table_fil %>%
  left_join(sequence_cluster_map, by = c("V5_norm" = "uc_sequence")) %>%
  left_join(cluster_members, by = "cluster_number") %>%
  mutate(
    member_sequence = if_else(is.na(member_sequence), V5, member_sequence),
    member_library  = if_else(is.na(member_library), centroid_library, member_library)
  ) %>%
  left_join(
    libraries_final %>% rename_with(~ paste0("member_", .), -V1_clean),
    by = c("member_library" = "V1_clean")
  )

cat("Rows in new_table_expanded:", nrow(new_table_expanded), "\n")
cat("Unique NCBI references (V1):", length(unique(new_table_expanded$V1)), "\n")

# Join NCBI ground-truth host from addtl_tip_plus_new_manual (status == "ncbi")
# V1 in addtl_tip_plus_new_manual has a trailing _1 suffix -> strip before joining
ncbi_host_truth <- addtl_tip_plus_new_manual %>%
  filter(status == "ncbi") %>%
  mutate(V1 = sub("_1$", "", V1)) %>%
  dplyr::select(V1, ncbi_host_broad = broader_gen, ncbi_host_specific = generalization) %>%
  distinct(V1, .keep_all = TRUE)

new_table_expanded <- new_table_expanded %>% dplyr::select(-any_of(c("ncbi_host_broad", "ncbi_host_specific")))
new_table_expanded <- left_join(new_table_expanded, ncbi_host_truth, by = "V1")

cat("Rows with ncbi_host_broad assigned:", sum(!is.na(new_table_expanded$ncbi_host_broad)),
    "of", nrow(new_table_expanded), "\n")
cat("ncbi_host_truth rows:", nrow(ncbi_host_truth), "\n")
cat("V1 overlap:", sum(unique(new_table_expanded$V1) %in% ncbi_host_truth$V1),
    "of", length(unique(new_table_expanded$V1)), "unique V1s\n")
# Show a few unmatched V1s for debugging
unmatched_v1 <- setdiff(unique(new_table_expanded$V1), ncbi_host_truth$V1)
if (length(unmatched_v1) > 0) cat("First 5 unmatched V1s:", head(unmatched_v1, 5), "\n")

# view(head(new_table_expanded, 100))

new_table_expanded$member_host_group[grepl("^SRR270786", new_table_expanded$member_library)] <- "Human"

# Add BioProject from L1_library_metadata.txt (column 17)
lib_bioproject <- read.table("files/L1_library_metadata.txt", sep = "\t", quote = "", fill = TRUE) %>%
  transmute(member_library = V1, member_bioproject = V17) %>%
  distinct(member_library, .keep_all = TRUE)
# Drop any pre-existing member_bioproject column (from re-runs) to avoid .x/.y suffixes
new_table_expanded <- new_table_expanded %>% dplyr::select(-any_of("member_bioproject"))
new_table_expanded <- left_join(new_table_expanded, lib_bioproject, by = "member_library")
new_table_expanded$member_bioproject[new_table_expanded$member_bioproject == ""] <- NA

# Fetch BioProject for libraries missing from metadata file via NCBI API
missing_libs <- new_table_expanded %>%
  filter(is.na(member_bioproject)) %>%
  pull(member_library) %>%
  unique()
cat("Fetching BioProject for", length(missing_libs), "missing libraries\n")

fetch_bioproject <- function(srr_ids) {
  results <- data.frame(member_library = character(), member_bioproject = character(), stringsAsFactors = FALSE)
  for (i in seq_along(srr_ids)) {
    srr <- srr_ids[i]
    if (i %% 50 == 0) cat(i, "of", length(srr_ids), "\n")
    tryCatch({
      sr <- rentrez::entrez_search(db = "sra", term = paste0(srr, "[ACCN]"))
      if (sr$count > 0) {
        lnk <- rentrez::entrez_link(dbfrom = "sra", db = "bioproject", id = sr$ids[1])
        bp_ids <- lnk$links$sra_bioproject
        if (!is.null(bp_ids) && length(bp_ids) > 0) {
          bp_summary <- rentrez::entrez_summary(db = "bioproject", id = bp_ids[1])
          results <- rbind(results, data.frame(
            member_library = srr,
            member_bioproject = bp_summary$project_acc,
            stringsAsFactors = FALSE
          ))
        }
      }
      Sys.sleep(0.35)
    }, error = function(e) cat("Error for", srr, ":", e$message, "\n"))
  }
  results
}

if (length(missing_libs) > 0) {
  bioproject_api <- fetch_bioproject(missing_libs)
  cat("Retrieved BioProject for", nrow(bioproject_api), "of", length(missing_libs), "libraries\n")
  # Fill in the missing values
  bp_lookup <- setNames(bioproject_api$member_bioproject, bioproject_api$member_library)
  idx <- is.na(new_table_expanded$member_bioproject) &
         new_table_expanded$member_library %in% names(bp_lookup)
  new_table_expanded$member_bioproject[idx] <- bp_lookup[new_table_expanded$member_library[idx]]
  # Cache results
  saveRDS(bioproject_api, "outputs/bioproject_api_cache.rds")
}

# Re-apply cached API results on re-runs
bioproject_api <- readRDS("outputs/bioproject_api_cache.rds")
bp_lookup <- setNames(bioproject_api$member_bioproject, bioproject_api$member_library)
idx <- is.na(new_table_expanded$member_bioproject) &
       new_table_expanded$member_library %in% names(bp_lookup)
new_table_expanded$member_bioproject[idx] <- bp_lookup[new_table_expanded$member_library[idx]]

# Fetch BioProject for any STILL missing libraries (e.g., not in first API run)
still_missing <- new_table_expanded %>%
  filter(is.na(member_bioproject)) %>%
  pull(member_library) %>%
  unique()
cat("Libraries still missing BioProject after cache:", length(still_missing), "\n")

if (length(still_missing) > 0) {
  bioproject_api_2 <- fetch_bioproject(still_missing)
  cat("Retrieved BioProject for", nrow(bioproject_api_2), "of", length(still_missing), "libraries\n")
  if (nrow(bioproject_api_2) > 0) {
    bp_lookup_2 <- setNames(bioproject_api_2$member_bioproject, bioproject_api_2$member_library)
    idx2 <- is.na(new_table_expanded$member_bioproject) &
            new_table_expanded$member_library %in% names(bp_lookup_2)
    new_table_expanded$member_bioproject[idx2] <- bp_lookup_2[new_table_expanded$member_library[idx2]]
    # Update cache with new results
    bioproject_api <- rbind(bioproject_api, bioproject_api_2) %>% distinct()
    saveRDS(bioproject_api, "outputs/bioproject_api_cache.rds")
  }
}

bioproject_api_2 <- readRDS("outputs/bioproject_api_cache.rds")

check_table <- filter(new_table_expanded, V1 == "X70829.1_102_1310")
check_table <- headnew_table_expanded, V1 == "X70829.1_102_1310")


accuracy_per_ncbi <- new_table_expanded %>%
  filter(!is.na(ncbi_host_broad), !is.na(member_host_group)) %>%
  dplyr::select(V1, member_library, ncbi_host_broad, ncbi_host_specific, member_host_group) %>%
  distinct(V1, member_library, .keep_all = TRUE) %>%   # dedupe members across centroids
  group_by(V1, ncbi_host_broad, ncbi_host_specific) %>%
  summarize(
    n_members      = n(),
    n_concordant   = sum(member_host_group == ncbi_host_specific, na.rm = TRUE),
    pct_concordant = n_concordant / n_members,
    .groups = "drop"
  ) %>%
  arrange(ncbi_host_broad, desc(pct_concordant))

hist(accuracy_per_ncbi$pct_concordant)

write.table(accuracy_per_ncbi, "outputs/host_concordance.tsv", sep = "\t", row.names = F, col.names = T)
write.table(new_table_expanded, "outputs/host_concordance_full_table.tsv", sep = "\t", row.names = F, col.names = T)


accuracy_per_ncbi$flag <- ifelse(accuracy_per_ncbi$n_members > 10 & accuracy_per_ncbi$pct_concordant < 0.5, "low_accuracy", "ok")

accuracy_per_ncbi_plot <- ggplot(data = accuracy_per_ncbi, aes(x = n_members, y = pct_concordant)) +
  geom_count(aes(color = after_stat(n)), alpha = 0.8) +
  scale_color_viridis_b(option = 'plasma') +
  geom_count(data = filter(accuracy_per_ncbi, flag == "low_accuracy"), , shape = 5,
             color = "red", alpha = 0.9, show.legend = FALSE, stroke = 0.5, size = 1.5) +
  scale_x_log10() + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'grey40') +
  geom_vline(xintercept = 10, linetype = 'dotted', color = 'grey40')


ggsave("outputs/2026.03.13.accuracy_per_ncbi_plot.png", accuracy_per_ncbi_plot, width = 10, height = 10, units = "cm", limitsize = F,bg='transparent')

# Quadrant breakdown: NCBI PV types by n_members (x=10) and pct_concordant (y=0.5)
accuracy_per_ncbi <- accuracy_per_ncbi %>%
  mutate(quadrant = case_when(
    n_members <= 10 & pct_concordant >= 0.5 ~ "UL (few, high concordance)",
    n_members >  10 & pct_concordant >= 0.5 ~ "UR (many, high concordance)",
    n_members <= 10 & pct_concordant <  0.5 ~ "LL (few, low concordance)",
    n_members >  10 & pct_concordant <  0.5 ~ "LR (many, low concordance)"
  ))

quadrant_summary <- accuracy_per_ncbi %>%
  group_by(quadrant) %>%
  summarize(
    n_ncbi_types = n(),
    n_libraries  = sum(n_members),
    .groups = "drop"
  ) %>%
  mutate(pct_types     = round(100 * n_ncbi_types / sum(n_ncbi_types), 1),
         pct_libraries = round(100 * n_libraries / sum(n_libraries), 1))

cat("\n========== Accuracy plot quadrant summary ==========\n")
print(quadrant_summary, n = 4)
cat("====================================================\n")

accuracy_per_ncbi_lower_right <- filter(accuracy_per_ncbi, n_members > 10 & pct_concordant < 0.5)
print(accuracy_per_ncbi_lower_right, max = NULL)
write.table(accuracy_per_ncbi_lower_right, "outputs/accuracy_per_ncbi_lower_right.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# Stacked bar: what Logan metadata says for lower-right quadrant NCBI PV types
lr_host_data <- new_table_expanded %>%
  filter(V1 %in% accuracy_per_ncbi_lower_right$V1,
         !is.na(member_host_group)) %>%
  distinct(V1, member_library, .keep_all = TRUE) %>%
  mutate(broad_host = if_else(member_host_group %in% names(tree_host_colors),
                              member_host_group, "Other")) %>%
  count(V1, ncbi_host_specific, broad_host)

# Use same color palette as phylogenetic tree
tree_host_colors <- c(
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

# Label: "NCBI_host (accession)" ordered by total observations
lr_labels <- lr_host_data %>%
  group_by(V1, ncbi_host_specific) %>%
  summarize(total = sum(n), .groups = "drop") %>%
  left_join(accuracy_per_ncbi_lower_right %>%
              dplyr::select(V1, pct_concordant), by = "V1") %>%
  arrange(desc(total)) %>%
  mutate(label = paste0(ncbi_host_specific, " (", V1, ")"),
         annotation = paste0("n=", total, ", ", round(pct_concordant * 100, 1), "% conc."))

lr_host_data <- lr_host_data %>%
  left_join(lr_labels %>% dplyr::select(V1, label, annotation, total), by = "V1") %>%
  mutate(label = factor(label, levels = rev(lr_labels$label)))

lr_bar_plot <- ggplot(lr_host_data, aes(x = label, y = n, fill = broad_host)) +
  geom_col(width = 0.7) +
  geom_text(data = lr_labels %>% mutate(label = factor(label, levels = rev(label))),
            aes(x = label, y = total, label = annotation, fill = NULL),
            hjust = -0.05, size = 2.8, family = "Noto Sans", color = "grey30") +
  scale_fill_manual(values = tree_host_colors, name = "Logan host metadata") +
  scale_x_discrete() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_flip() +
  labs(x = NULL, y = "Number of SRA libraries",
       title = "Lower-right quadrant: Logan host calls\nfor discordant NCBI PVs (per-library)") +
  theme_bw() +
  theme(text = element_text(family = "Noto Sans", size = 10),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("outputs/accuracy_lower_right_host_bar_library.png", lr_bar_plot,
       width = 22, height = 14, units = "cm", dpi = 300, bg = "white")

###different vis test####
accuracy_per_ncbi <- accuracy_per_ncbi %>%
  mutate(obs_bin = cut(n_members, 
                       breaks = c(0, 1, 5, 10, 50, 100, Inf),
                       labels = c("1", "2–5", "6–10", "11–50", "51–100", "100+"),
                       right = TRUE))

accuracy_per_ncbi_violin <- ggplot(data = accuracy_per_ncbi, 
                                    aes(x = obs_bin, y = pct_concordant)) +
  geom_violin(fill = "grey85", color = "grey50", scale = "width", trim = FALSE) +
  geom_jitter(data = filter(accuracy_per_ncbi, flag == "low_accuracy"), aes(color = flag), 
              width = 0.15, height = 0, alpha = 0.4, size = 1.2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey30"), guide = "none") +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'grey40') +
  theme_bw() +
  theme(text = element_text(family = "Noto Sans"),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Number of SRA libraries", y = "Concordance with NCBI host")



# ---- Same NCBI concordance analysis using BioProject as observation unit ----
# Show members with no BioProject assigned
no_bioproject_members <- new_table_expanded %>%
  filter(!is.na(ncbi_host_broad), !is.na(member_host_group), is.na(member_bioproject)) %>%
  dplyr::select(V1, member_library, ncbi_host_broad, ncbi_host_specific, member_host_group) %>%
  distinct(V1, member_library, .keep_all = TRUE)
cat("Members with no BioProject:", nrow(no_bioproject_members), "across",
    length(unique(no_bioproject_members$V1)), "NCBI refs\n")
cat("Unique libraries missing BioProject:", length(unique(no_bioproject_members$member_library)), "\n")

accuracy_per_ncbi_bp <- new_table_expanded %>%
  filter(!is.na(ncbi_host_broad), !is.na(member_host_group), !is.na(member_bioproject)) %>%
  dplyr::select(V1, member_bioproject, ncbi_host_broad, ncbi_host_specific, member_host_group) %>%
  distinct(V1, member_bioproject, member_host_group, .keep_all = TRUE) %>%   # keep unique (V1, bioproject, host_group) combos
  group_by(V1, ncbi_host_broad, ncbi_host_specific) %>%
  summarize(
    n_bioprojects    = n(),
    n_concordant     = sum(member_host_group == ncbi_host_specific, na.rm = TRUE),
    pct_concordant   = n_concordant / n_bioprojects,
    .groups = "drop"
  ) %>%
  arrange(ncbi_host_broad, desc(pct_concordant))

accuracy_per_ncbi_bp$flag <- ifelse(accuracy_per_ncbi_bp$n_bioprojects > 10 & accuracy_per_ncbi_bp$pct_concordant < 0.5, "low_accuracy", "ok")

accuracy_per_ncbi_bp_plot <- ggplot(data = accuracy_per_ncbi_bp, aes(x = n_bioprojects, y = pct_concordant)) +
  geom_count(aes(color = after_stat(n)), alpha = 0.8) +
  scale_color_viridis_b(option = 'plasma') +
  geom_count(data = filter(accuracy_per_ncbi_bp, flag == "low_accuracy"), shape = 5,
             color = "red", alpha = 0.9, show.legend = FALSE, stroke = 0.5, size = 1.5) +
  scale_x_log10() + theme_bw() + theme(text = element_text(family = "Noto Sans")) +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'grey40') +
  geom_vline(xintercept = 10, linetype = 'dotted', color = 'grey40')

accuracy_per_ncbi_bp_plot
accuracy_per_ncbi_bp_lower_right <- filter(accuracy_per_ncbi_bp, n_bioprojects > 10 & pct_concordant < 0.5)
write.table(accuracy_per_ncbi_bp_lower_right, "outputs/accuracy_per_ncbi_bp_lower_right.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


ggsave("outputs/2026.03.16.accuracy_per_ncbi_bioproject_plot.png", accuracy_per_ncbi_bp_plot, width = 10, height = 10, units = "cm", limitsize = F, bg = 'transparent')

accuracy_per_ncbi_bp <- accuracy_per_ncbi_bp %>%
  mutate(quadrant = case_when(
    n_bioprojects <= 10 & pct_concordant >= 0.5 ~ "UL (few, high concordance)",
    n_bioprojects >  10 & pct_concordant >= 0.5 ~ "UR (many, high concordance)",
    n_bioprojects <= 10 & pct_concordant <  0.5 ~ "LL (few, low concordance)",
    n_bioprojects >  10 & pct_concordant <  0.5 ~ "LR (many, low concordance)"
  ))

quadrant_summary <- accuracy_per_ncbi_bp %>%
  group_by(quadrant) %>%
  summarize(
    n_ncbi_types = n(),
    n_bioprojects  = sum(n_bioprojects),
    .groups = "drop"
  ) %>%
  mutate(pct_types     = round(100 * n_ncbi_types / sum(n_ncbi_types), 1),
         pct_libraries = round(100 * n_bioprojects / sum(n_bioprojects), 1))

cat("\n========== Accuracy plot quadrant summary ==========\n")
print(quadrant_summary, n = 4)
cat("====================================================\n")

write.table(accuracy_per_ncbi_bp, "outputs/host_concordance_bioproject.tsv", sep = "\t", row.names = F, col.names = T)

# Stacked bar: what Logan metadata says for lower-right quadrant (per-BioProject)
lr_bp_host_data <- new_table_expanded %>%
  filter(V1 %in% accuracy_per_ncbi_bp_lower_right$V1,
         !is.na(member_host_group), !is.na(member_bioproject)) %>%
  distinct(V1, member_bioproject, member_host_group, .keep_all = TRUE) %>%
  mutate(broad_host = if_else(member_host_group %in% names(tree_host_colors),
                              member_host_group, "Other")) %>%
  count(V1, ncbi_host_specific, broad_host)

lr_bp_labels <- lr_bp_host_data %>%
  group_by(V1, ncbi_host_specific) %>%
  summarize(total = sum(n), .groups = "drop") %>%
  left_join(accuracy_per_ncbi_bp_lower_right %>%
              dplyr::select(V1, pct_concordant), by = "V1") %>%
  arrange(desc(total)) %>%
  mutate(label = paste0(ncbi_host_specific, " (", V1, ")"),
         annotation = paste0("n=", total, ", ", round(pct_concordant * 100, 1), "% conc."))

lr_bp_host_data <- lr_bp_host_data %>%
  left_join(lr_bp_labels %>% dplyr::select(V1, label, annotation, total), by = "V1") %>%
  mutate(label = factor(label, levels = rev(lr_bp_labels$label)))

lr_bp_bar_plot <- ggplot(lr_bp_host_data, aes(x = label, y = n, fill = broad_host)) +
  geom_col(width = 0.7) +
  geom_text(data = lr_bp_labels %>% mutate(label = factor(label, levels = rev(label))),
            aes(x = label, y = total, label = annotation, fill = NULL),
            hjust = -0.05, size = 2.8, family = "Noto Sans", color = "grey30") +
  scale_fill_manual(values = tree_host_colors, name = "Logan host metadata") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_flip() +
  labs(x = NULL, y = "Number of BioProjects",
       title = "Lower-right quadrant: Logan host calls\nfor discordant NCBI PVs (per-BioProject)") +
  theme_bw() +
  theme(text = element_text(family = "Noto Sans", size = 10),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("outputs/accuracy_lower_right_host_bar_bioproject.png", lr_bp_bar_plot,
       width = 22, height = 14, units = "cm", dpi = 300, bg = "white")

# Combined two-panel figure: library + bioproject lower-right barcharts
lr_combined <- ggarrange(
  lr_bar_plot + theme(legend.position = "none"),
  lr_bp_bar_plot + theme(legend.position = "none"),
  ncol = 2, nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 14, family = "Noto Sans")
)

ggsave("outputs/accuracy_lower_right_combined.png", lr_combined,
       width = 36, height = 16, units = "cm", dpi = 300, bg = "white")

# ---- Novel sequences: UC members NOT in high-confidence NCBI clusters ----
# High-confidence clusters: only those where the V5 hit IS the cluster centroid (S record)
uc_centroids <- all_full_l1s_logan_and_pr_90_uc %>%
  filter(V1 == "S") %>%
  transmute(uc_sequence = V9, cluster_number = V2)

high_conf_clusters <- new_table_fil %>%
  distinct(V5_norm) %>%
  inner_join(uc_centroids, by = c("V5_norm" = "uc_sequence")) %>%
  pull(cluster_number) %>%
  unique()
cat("High-confidence NCBI clusters (centroid is the hit):", length(high_conf_clusters), "\n")

# All sequences in the UC file
all_uc_sequences <- all_full_l1s_logan_and_pr_90_uc %>%
  filter(V1 %in% c("S", "H")) %>%
  transmute(sequence = V9, cluster_number = V2,
            library = get_lib(V9))

# Split into 3 categories:
# 1) High-confidence NCBI: members of the 646 centroid clusters
# 2) Novel: members of clusters with NO NCBI match at all
# 3) Low-confidence NCBI: members of the remaining 42 clusters (matched via other V1+V5 pairs)
all_cluster_numbers <- unique(all_uc_sequences$cluster_number)
ncbi_all_clusters <- unique(new_table_expanded$cluster_number)
low_conf_clusters <- setdiff(ncbi_all_clusters, high_conf_clusters)
novel_cluster_numbers <- setdiff(all_cluster_numbers, ncbi_all_clusters)

cat("High-confidence NCBI clusters:", length(high_conf_clusters), "\n")
cat("Low-confidence NCBI clusters:", length(low_conf_clusters), "\n")
cat("Novel clusters (no NCBI match):", length(novel_cluster_numbers), "\n")
cat("Total:", length(high_conf_clusters) + length(low_conf_clusters) + length(novel_cluster_numbers), "of", length(all_cluster_numbers), "\n")

novel_members <- all_uc_sequences %>%
  filter(cluster_number %in% novel_cluster_numbers) %>%
  distinct(sequence, .keep_all = TRUE)

low_conf_members <- all_uc_sequences %>%
  filter(cluster_number %in% low_conf_clusters) %>%
  distinct(sequence, .keep_all = TRUE)

cat("Novel sequences:", nrow(novel_members), "\n")
cat("Low-confidence NCBI sequences:", nrow(low_conf_members), "\n")
cat("High-confidence NCBI sequences:", sum(all_uc_sequences$cluster_number %in% high_conf_clusters), "\n")

# Add centroid info from UC (S records = centroids, H records point to centroid in V10)
uc_centroid_lookup <- all_full_l1s_logan_and_pr_90_uc %>%
  filter(V1 == "H") %>%
  transmute(sequence = V9, centroid = V10) %>%
  bind_rows(
    # S records are their own centroid
    all_full_l1s_logan_and_pr_90_uc %>%
      filter(V1 == "S") %>%
      transmute(sequence = V9, centroid = V9)
  ) %>%
  distinct(sequence, .keep_all = TRUE) %>%
  mutate(centroid_library = get_lib(centroid))

# Helper: add host metadata + centroid info to a members table
add_member_metadata <- function(members_df) {
  members_df %>%
    left_join(uc_centroid_lookup, by = "sequence") %>%
    left_join(libraries_final, by = c("library" = "V1_clean")) %>%
    left_join(lib_bioproject, by = c("library" = "member_library")) %>%
    left_join(
      libraries_final %>% dplyr::select(V1_clean, host_group) %>% rename(centroid_host_group = host_group),
      by = c("centroid_library" = "V1_clean")
    ) %>%
    mutate(host_group = if_else(grepl("^SRR270786", library), "Human", host_group))
}

novel_members <- add_member_metadata(novel_members)
low_conf_members <- add_member_metadata(low_conf_members)

# Fix duplicated library name bug (e.g., SRR24184635SRR24184635 -> SRR24184635)
novel_members$library <- gsub("^(SRR24184635)\\1$", "\\1", novel_members$library)
novel_members$host_group[novel_members$library == "SRR24184635"] <- "Amphibian"
novel_members$original_metadata[novel_members$library == "SRR24184635"] <- "Zhangixalus dugritei"

cat("\n--- Novel sequence host group counts ---\n")
print(table(novel_members$host_group, useNA = "ifany"))
cat("\n--- Low-confidence NCBI sequence host group counts ---\n")
print(table(low_conf_members$host_group, useNA = "ifany"))

cat("\n--- Novel sequence BioProject counts (top 20) ---\n")
print(head(sort(table(novel_members$member_bioproject), decreasing = TRUE), 20))


test_table <- as.data.frame(table(novel_members$cluster_number, novel_members$host_group))

summary_observations <- test_table %>% 
                        pivot_wider(names_from = Var2, values_from = Freq, names_prefix = 'inst_')

summary_observations$observation_n <- rowSums(summary_observations[-1], na.rm = TRUE)
hist(summary_observations$observation_n)

# ---- Same analysis but using BioProject as the observation unit ----
# Deduplicate: one row per (cluster, bioproject, host_group) instead of per library
novel_by_bioproject <- novel_members %>%
  filter(!is.na(member_bioproject)) %>%
  distinct(cluster_number, member_bioproject, host_group)

test_table_bp <- as.data.frame(table(novel_by_bioproject$cluster_number, novel_by_bioproject$host_group))

summary_observations_bp <- test_table_bp %>%
  pivot_wider(names_from = Var2, values_from = Freq, names_prefix = 'inst_')

summary_observations_bp$observation_n <- rowSums(summary_observations_bp[-1], na.rm = TRUE)
hist(summary_observations_bp$observation_n)

# ---- Summary: all clusters by novelty status and host group ----
# Three-way classification:
#   1) ncbi_virus (688): matched via alignment_envs.table (NCBI virus refs aligned to Logan centroids)
#   2) blastn 26): matched via blastn_hits_centroids.tsv (Logan centroids vs NCBI nt, >=90% identity)
#      but NOT already captured by alignment_envs.table
#   3) novel (386): no match in either search

# Read blastn results and identify clusters with >=90% identity hits to NCBI nt
blastn_results <- read.table("files/blastn_hits_centroids.tsv", sep = "\t", header = FALSE)
blastn_results$qcov <- abs(blastn_results$V2 - blastn_results$V3) / blastn_results$V4

# Filter: same thresholds as logan_3_novelty.R (e-value < 1e-5, qcov >= 0.7, identity >= 90%)
blastn_fil <- blastn_results %>%
  filter(V10 < 1e-5, qcov >= 0.7) %>%
  group_by(V1) %>% slice_max(n = 1, V9) %>% ungroup() %>%
  distinct(V1, .keep_all = TRUE) %>%
  filter(V9 >= 90)

# Map blastn query sequences to UC cluster numbers
# blastn V1 names match UC V9 names exactly — no normalization needed
blastn_clusters <- blastn_fil %>%
  left_join(sequence_cluster_map, by = c("V1" = "uc_sequence")) %>%
  filter(!is.na(cluster_number)) %>%
  pull(cluster_number) %>%
  unique()

# blastn-only clusters: in blastn but NOT in alignment_envs.table (ncbi_all_clusters)
blastn_only_clusters <- setdiff(blastn_clusters, ncbi_all_clusters)
cat("  Blastn clusters (total):   ", length(blastn_clusters), "\n")
cat("  Overlap with ncbi_virus:   ", length(intersect(blastn_clusters, ncbi_all_clusters)), "\n")
cat("  Blastn-only clusters:      ", length(blastn_only_clusters), "\n")

# Build one row per cluster using the centroid's host group
all_clusters_summary <- all_full_l1s_logan_and_pr_90_uc %>%
  filter(V1 == "S") %>%
  transmute(cluster_number = V2, centroid = V9,
            centroid_library = get_lib(V9)) %>%
  left_join(libraries_final %>% dplyr::select(V1_clean, host_group),
            by = c("centroid_library" = "V1_clean")) %>%
  mutate(
    # Apply the same manual fixes used elsewhere
    host_group = if_else(grepl("^SRR270786", centroid_library), "Human", host_group),
    host_group = if_else(grepl("^SRR24184635", centroid_library), "Amphibian", host_group),
    # Three-way novelty classification
    status = case_when(
      cluster_number %in% ncbi_all_clusters    ~ "ncbi_virus",
      cluster_number %in% blastn_only_clusters ~ "blastn",
      TRUE ~ "novel"
    )
  )

# For known clusters (ncbi_virus), also grab NCBI host from the concordance data
known_ncbi_host <- new_table_expanded %>%
  distinct(cluster_number, ncbi_host_specific) %>%
  filter(!is.na(ncbi_host_specific)) %>%
  group_by(cluster_number) %>%
  count(ncbi_host_specific, sort = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(cluster_number, ncbi_host = ncbi_host_specific)

all_clusters_summary <- left_join(all_clusters_summary, known_ncbi_host, by = "cluster_number")

# Count members per cluster
members_per_cluster <- all_uc_sequences %>%
  count(cluster_number, name = "n_members")
all_clusters_summary <- left_join(all_clusters_summary, members_per_cluster, by = "cluster_number")

cat("\n========== All clusters: novelty status by host group ==========\n")
cat("Total clusters:", nrow(all_clusters_summary), "\n")
cat("NCBI virus:", sum(all_clusters_summary$status == "ncbi_virus"), "\n")
cat("Blastn:    ", sum(all_clusters_summary$status == "blastn"), "\n")
cat("Novel:     ", sum(all_clusters_summary$status == "novel"), "\n\n")

cluster_host_table <- table(all_clusters_summary$status, all_clusters_summary$host_group, useNA = "ifany")
print(addmargins(cluster_host_table))

cat("\n--- As a dataframe (cluster counts) ---\n")
cluster_host_df <- all_clusters_summary %>%
  count(host_group, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  mutate(total = ncbi_virus + blastn + novel,
         host_group = if_else(total < 9, "Other", as.character(host_group))) %>%
  group_by(host_group) %>%
  summarize(across(c(blastn, ncbi_virus, novel), sum), .groups = "drop") %>%
  mutate(total_known = blastn + ncbi_virus,
         total = ncbi_virus + blastn + novel) %>%
  select(host_group, blastn, ncbi_virus, total_known, novel, total) %>%
  arrange(desc(total))
print(cluster_host_df, n = nrow(cluster_host_df))
cat("=================================================================\n")
all_clusters_summary %>% filter(is.na(host_group)) %>% select(cluster_number, centroid, centroid_library)
write.table(select(cluster_host_df, host_group, total_known, novel), "outputs/table_3.tsv", sep = "\t", quote = F, row.names = F)

# SRR6976994: mostly simple repeat artifacts. Only keep members of clusters
# where the centroid is NOT from SRR6976994 (e.g., cluster 282 Zea mays PV).
srr6976994_centroid_clusters <- all_clusters_summary %>%
  filter(centroid_library == "SRR6976994") %>%
  pull(cluster_number)

# ========== Sequences by library type and novelty status ==========
# Read SRA metadata for library source (METAGENOMIC, GENOMIC, etc.)
sra_metadata <- read.table("files/all_hits_info.list", sep = ",", fill = TRUE,
                           header = TRUE, quote = "", comment.char = "#")
sra_metadata <- unique(sra_metadata)

# Map cluster -> status from all_clusters_summary
cluster_status <- all_clusters_summary %>% dplyr::select(cluster_number, status)

# Join each sequence to its cluster status and library source
seq_by_type <- all_uc_sequences %>%
  filter(!(library == "SRR6976994" & cluster_number %in% srr6976994_centroid_clusters)) %>%
  left_join(cluster_status, by = "cluster_number") %>%
  left_join(sra_metadata %>% dplyr::select(Run, LibrarySource),
            by = c("library" = "Run")) %>%
  mutate(
    # Collapse rare categories into OTHER
    lib_type = case_when(
      LibrarySource %in% c("METAGENOMIC", "GENOMIC", "TRANSCRIPTOMIC",
                           "METATRANSCRIPTOMIC", "GENOMIC SINGLE CELL",
                           "TRANSCRIPTOMIC SINGLE CELL") ~ LibrarySource,
      TRUE ~ "OTHER"
    )
  )

cat("\n========== Full L1 sequences by library type and novelty status ==========\n")
cat("Total sequences:", nrow(seq_by_type), "\n\n")

all_uc_sequences

seq_type_table <- seq_by_type %>%
  count(lib_type, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  mutate(total_known = blastn + ncbi_virus,
         total = rowSums(across(c(blastn, ncbi_virus, novel)))) %>%
  arrange(desc(total))
print(seq_type_table, n = nrow(seq_type_table))
cat("=========================================================================\n")

seq_type_table <- seq_type_table %>% 
    rownames_to_column("id") %>%
    gather(status, values, -lib_type)


seq_type_table <- filter(seq_type_table, grepl('total_known|novel', status))
seq_type_table$values <- as.numeric(seq_type_table$values)

ggplot(seq_type_table, aes(reorder(lib_type, values), y=values, fill = status)) + 
  geom_bar(stat = "identity") + theme_classic() + scale_fill_manual(values = c("#4646d4", "grey")) +
  coord_flip()+ 
  theme(text = element_text(family = "Noto Sans", size = 25))+ 
  theme(axis.text.y = element_text(hjust=1))

ggsave("outputs/2026.03.30.breakdown_new_small.jpg", width = 25, height = 15, units = "cm", limitsize = F)

#remove clusters where SRR6976994 is the centroid (simple repeat artifacts)
all_clusters_novel <- filter(all_clusters_summary, !cluster_number %in% srr6976994_centroid_clusters)
#get list of novel ones to put together a new tree

all_clusters_novel <- filter(all_clusters_novel, status == "novel")$centroid
write.table(all_clusters_novel, "outputs/all_novel_centroids.tsv", sep = "\t", quote = F, row.names = F, col.names = F)

# =====================================================================
# =================== SUMMARY FIGURES =================================
# =====================================================================

library(ggpubr)

# Consistent theme for all plots
theme_summary <- theme_bw() +
  theme(text = element_text(family = "Noto Sans", size = 11),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank())

# Shared filters: remove SRR6976994 and collapse small host groups into "Other"
small_host_groups <- all_clusters_summary %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters, !is.na(host_group)) %>%
  count(host_group) %>%
  filter(n < 10) %>%
  pull(host_group)

all_clusters_summary_fig <- all_clusters_summary %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters) %>%
  mutate(host_group = as.character(host_group),
         host_group = if_else(host_group %in% small_host_groups, "Other", host_group))

seq_by_type_fig <- seq_by_type %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters)

novel_members_fig <- novel_members %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters)

new_table_expanded_fig <- new_table_expanded %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters)

libraries_final_fig <- libraries_final %>%
  filter(V1_clean != "SRR6976994")

# ---- Summary: individual full-length PVs per host group ----
pv_per_host <- all_uc_sequences %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters) %>%
  left_join(all_clusters_summary %>% dplyr::select(cluster_number, status, host_group),
            by = "cluster_number") %>%
  left_join(lib_bioproject, by = c("library" = "member_library")) %>%
  mutate(host_group = as.character(host_group),
         host_group = if_else(host_group %in% small_host_groups, "Other", host_group)) %>%
  filter(!is.na(host_group)) %>%
  group_by(host_group) %>%
  summarize(
    n_sequences  = n(),
    n_novel      = sum(status == "novel", na.rm = TRUE),
    n_known      = sum(status != "novel", na.rm = TRUE),
    n_clusters   = n_distinct(cluster_number),
    n_libraries  = n_distinct(library),
    n_bioprojects = n_distinct(member_bioproject[!is.na(member_bioproject)]),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sequences))

cat("\n========== Full-length PVs per host group ==========\n")
print(pv_per_host, n = nrow(pv_per_host))
cat("Total sequences:", sum(pv_per_host$n_sequences), "\n")
cat("====================================================\n")

# ---- Fig 0: PV per host summary (horizontal stacked bar + annotations) ----
fig0_data <- pv_per_host %>%
  pivot_longer(cols = c(n_novel, n_known), names_to = "type", values_to = "count") %>%
  mutate(type = factor(type, levels = c("n_known", "n_novel"),
                       labels = c("Known", "Novel")))

fig0_data$host_group <- factor(fig0_data$host_group,
                                levels = rev(pv_per_host$host_group))  # bottom-to-top by n_sequences

fig0_annotations <- pv_per_host %>%
  mutate(host_group = factor(host_group, levels = rev(host_group)),
         label = paste0(n_clusters, "c / ", n_libraries, "L / ", n_bioprojects, "BP"))

fig0 <- ggplot(fig0_data, aes(x = count, y = host_group, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(data = fig0_annotations,
            aes(x = n_sequences, y = host_group, label = label, fill = NULL),
            hjust = -0.05, size = 2.8, family = "Noto Sans", color = "grey30") +
  scale_fill_manual(values = c("Known" = "grey", "Novel" = "#4646d4"), name = NULL) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Number of full-length L1 sequences", y = NULL,
       title = "Full-length PV sequences per host group",
       subtitle = "Annotations: clusters / libraries / BioProjects") +
  theme_summary +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"))

ggsave("outputs/summary_fig0_pv_per_host.png", fig0,
       width = 22, height = 16, units = "cm", dpi = 300, bg = "white")

# ---- Fig 1: Library type bar chart (sequences by library source + novelty) ----
fig1_data <- seq_by_type_fig %>%
  filter(!is.na(lib_type)) %>%
  count(lib_type, status) %>%
  mutate(status = factor(status, levels = c("novel", "blastn", "ncbi_virus"),
                         labels = c("Novel", "Blastn only", "NCBI virus")))

fig1 <- ggplot(fig1_data, aes(x = reorder(lib_type, -n, sum), y = n, fill = status)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Novel" = "#4646d4", "Blastn only" = "grey60", "NCBI virus" = "grey80"),
                    name = "Status") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Number of full-length L1 sequences",
       title = "Full-length L1 sequences by library source") +
  theme_summary +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("outputs/summary_fig1_library_type_bar.png", fig1,
       width = 18, height = 12, units = "cm", dpi = 300, bg = "white")

# ---- Fig 2: Host group bar chart — clusters by novelty status ----
fig2_data <- all_clusters_summary_fig %>%
  filter(!is.na(host_group)) %>%
  count(host_group, status) %>%
  mutate(status = factor(status, levels = c("novel", "blastn", "ncbi_virus"),
                         labels = c("Novel", "Blastn only", "NCBI virus")))

# Order host groups by total count
host_order <- fig2_data %>% group_by(host_group) %>% summarize(tot = sum(n)) %>% arrange(desc(tot)) %>% pull(host_group)
fig2_data$host_group <- factor(fig2_data$host_group, levels = host_order)

fig2 <- ggplot(fig2_data, aes(x = host_group, y = n, fill = status)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = c("Novel" = "#4646d4", "Blastn only" = "grey60", "NCBI virus" = "grey80"),
                    name = "Status") +
  labs(x = NULL, y = "Number of 90% clusters",
       title = "PV clusters by host group and novelty status") +
  theme_summary +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outputs/summary_fig2_host_novelty_bar.png", fig2,
       width = 22, height = 14, units = "cm", dpi = 300, bg = "white")

# ---- Fig 3: Proportion stacked bar — novel vs known per host group ----
fig3_data <- all_clusters_summary_fig %>%
  filter(!is.na(host_group)) %>%
  mutate(known_novel = if_else(status == "novel", "Novel", "Known")) %>%
  count(host_group, known_novel) %>%
  group_by(host_group) %>%
  mutate(prop = n / sum(n), total = sum(n)) %>%
  ungroup()

fig3_data$host_group <- factor(fig3_data$host_group, levels = host_order)

fig3 <- ggplot(fig3_data, aes(x = host_group, y = prop, fill = known_novel)) +
  geom_col(width = 0.75) +
  geom_text(data = fig3_data %>% filter(known_novel == "Novel"),
            aes(label = paste0("n=", total)), y = 1.03, size = 2.8, family = "Noto Sans") +
  scale_fill_manual(values = c("Known" = "grey", "Novel" = "#4646d4"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Proportion of clusters",
       title = "Proportion of novel vs known PV clusters per host") +
  theme_summary +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outputs/summary_fig3_host_proportion.png", fig3,
       width = 22, height = 14, units = "cm", dpi = 300, bg = "white")

# ---- Fig 4: Host confidence resolution pie/donut ----
fig4_data <- libraries_final_fig %>%
  mutate(source = case_when(
    host_confidence == "high" ~ "SRA metadata (high)",
    host_confidence == "medium" & is.na(biosample_organism) ~ "SRA metadata (medium)",
    host_confidence == "medium" & !is.na(biosample_organism) ~ "BioSample API",
    host_confidence == "manual" ~ "Manual curation",
    host_confidence == "low" ~ "Unresolved (low)",
    TRUE ~ "Unknown"
  )) %>%
  count(source) %>%
  mutate(pct = n / sum(n),
         label = paste0(source, "\n", scales::comma(n), " (", scales::percent(pct, 0.1), ")"))

fig4 <- ggplot(fig4_data, aes(x = 2, y = n, fill = source)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = c("SRA metadata (high)" = "#2A9D8F",
                                "SRA metadata (medium)" = "#8AB17D",
                                "BioSample API" = "#457B9D",
                                "Manual curation" = "#F4A261",
                                "Unresolved (low)" = "#E76F51",
                                "Unknown" = "grey70"),
                    name = "Resolution source") +
  labs(title = "Host metadata resolution sources") +
  theme_void() +
  theme(text = element_text(family = "Noto Sans", size = 11))

ggsave("outputs/summary_fig4_host_resolution_donut.png", fig4,
       width = 16, height = 14, units = "cm", dpi = 300, bg = "white")

# ---- Fig 5: Cluster size distribution — novel vs known (violin + boxplot) ----
fig5_data <- all_clusters_summary_fig %>%
  mutate(known_novel = if_else(status == "novel", "Novel", "Known"))

fig5 <- ggplot(fig5_data, aes(x = known_novel, y = n_members, fill = known_novel)) +
  geom_violin(alpha = 0.4, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_y_log10() +
  scale_fill_manual(values = c("Known" = "grey", "Novel" = "#4646d4"), guide = "none") +
  labs(x = NULL, y = "Cluster members (log scale)",
       title = "Cluster size: novel vs known PV clusters") +
  theme_summary

ggsave("outputs/summary_fig5_cluster_size_violin.png", fig5,
       width = 12, height = 12, units = "cm", dpi = 300, bg = "white")

# ---- Fig 6: Top host groups — grouped bar (known vs novel counts) ----
cluster_host_df_fig <- all_clusters_summary_fig %>%
  filter(!is.na(host_group)) %>%
  count(host_group, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  mutate(total_known = blastn + ncbi_virus,
         total = ncbi_virus + blastn + novel) %>%
  arrange(desc(total))

fig6_data <- cluster_host_df_fig %>%
  filter(host_group != "Other") %>%
  dplyr::select(host_group, total_known, novel) %>%
  pivot_longer(cols = c(total_known, novel), names_to = "type", values_to = "count") %>%
  mutate(type = factor(type, levels = c("total_known", "novel"),
                       labels = c("Known", "Novel")))

fig6_data$host_group <- factor(fig6_data$host_group,
                                levels = cluster_host_df_fig$host_group[cluster_host_df_fig$host_group != "Other"])

fig6 <- ggplot(fig6_data, aes(x = host_group, y = count, fill = type)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Known" = "grey", "Novel" = "#4646d4"), name = NULL) +
  labs(x = NULL, y = "Number of 90% clusters",
       title = "Known vs novel PV clusters per host group") +
  theme_summary +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outputs/summary_fig6_known_vs_novel_dodge.png", fig6,
       width = 22, height = 14, units = "cm", dpi = 300, bg = "white")

# ---- Fig 7: Concordance heatmap — NCBI host vs SRA host ----
fig7_data <- new_table_expanded_fig %>%
  filter(!is.na(ncbi_host_broad), !is.na(member_host_group)) %>%
  distinct(V1, member_library, .keep_all = TRUE) %>%
  count(ncbi_host_specific, member_host_group) %>%
  group_by(ncbi_host_specific) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Keep only host groups with enough observations
big_hosts <- fig7_data %>%
  group_by(ncbi_host_specific) %>%
  summarize(total = sum(n)) %>%
  filter(total >= 20) %>%
  pull(ncbi_host_specific)

fig7_data_filt <- fig7_data %>% filter(ncbi_host_specific %in% big_hosts)

fig7 <- ggplot(fig7_data_filt, aes(x = member_host_group, y = ncbi_host_specific, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n), size = 2.5, family = "Noto Sans") +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Proportion",
                       labels = scales::percent) +
  labs(x = "SRA library host group", y = "NCBI reference host",
       title = "Host concordance: NCBI reference vs SRA metadata") +
  theme_minimal() +
  theme(text = element_text(family = "Noto Sans", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("outputs/summary_fig7_concordance_heatmap.png", fig7,
       width = 22, height = 16, units = "cm", dpi = 300, bg = "white")

# ---- Fig 8: Novel cluster observations — histogram of how many libraries per novel cluster ----
fig8_data <- novel_members_fig %>%
  distinct(cluster_number, library) %>%
  count(cluster_number, name = "n_libraries")

fig8 <- ggplot(fig8_data, aes(x = n_libraries)) +
  geom_histogram(binwidth = 1, fill = "#4646d4", color = "white", alpha = 0.85) +
  scale_x_continuous(breaks = seq(0, max(fig8_data$n_libraries, na.rm = TRUE), by = 5)) +
  labs(x = "Number of SRA libraries per novel cluster",
       y = "Number of clusters",
       title = "Distribution of library support for novel PV clusters") +
  theme_summary

ggsave("outputs/summary_fig8_novel_library_support.png", fig8,
       width = 16, height = 12, units = "cm", dpi = 300, bg = "white")

# ---- Fig 9: Library type proportions for novel vs known (percent stacked bar) ----
fig9_data <- seq_by_type_fig %>%
  filter(!is.na(lib_type)) %>%
  mutate(known_novel = if_else(status == "novel", "Novel", "Known")) %>%
  count(known_novel, lib_type) %>%
  group_by(known_novel) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

fig9 <- ggplot(fig9_data, aes(x = known_novel, y = prop, fill = lib_type)) +
  geom_col(width = 0.6) +
  scale_fill_brewer(palette = "Set2", name = "Library source") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of sequences",
       title = "Library source composition: novel vs known sequences") +
  theme_summary

ggsave("outputs/summary_fig9_libtype_proportion.png", fig9,
       width = 14, height = 12, units = "cm", dpi = 300, bg = "white")

# ---- Fig 10: Multi-panel summary (ggarrange) ----
fig_combined <- ggarrange(fig1, fig2, fig5, fig3,
                          ncol = 2, nrow = 2,
                          labels = c("A", "B", "C", "D"),
                          font.label = list(size = 14, family = "Noto Sans"))

ggsave("outputs/summary_fig10_combined_panel.png", fig_combined,
       width = 36, height = 28, units = "cm", dpi = 300, bg = "white")

cat("\n========== Summary figures saved to outputs/ ==========\n")
cat("  Fig 1: summary_fig1_library_type_bar.png\n")
cat("  Fig 2: summary_fig2_host_novelty_bar.png\n")
cat("  Fig 3: summary_fig3_host_proportion.png\n")
cat("  Fig 4: summary_fig4_host_resolution_donut.png\n")
cat("  Fig 5: summary_fig5_cluster_size_violin.png\n")
cat("  Fig 6: summary_fig6_known_vs_novel_dodge.png\n")
cat("  Fig 7: summary_fig7_concordance_heatmap.png\n")
cat("  Fig 8: summary_fig8_novel_library_support.png\n")
cat("  Fig 9: summary_fig9_libtype_proportion.png\n")
cat("  Fig 10: summary_fig10_combined_panel.png\n")
cat("=======================================================\n")

# How many Ray-finned Fish clusters exist before the summary?
all_clusters_summary %>% filter(host_group == "Ray-finned Fish") %>% nrow()

# How many sequences lack a libraries_final match?
all_uc_sequences %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters) %>%
  left_join(libraries_final %>% dplyr::select(V1_clean, host_group),
            by = c("library" = "V1_clean")) %>%
  left_join(all_clusters_summary %>% dplyr::select(cluster_number, host_group) %>%
              rename(centroid_host = host_group), by = "cluster_number") %>%
  filter(centroid_host == "Ray-finned Fish") %>%
  count(is.na(host_group))

# =====================================================================
# ========= Fisher's exact test: Logan vs NCBI host proportions =======
# =====================================================================
# Compare host group proportions between Logan clusters (~1100) and
# NCBI reference PVs (~992) from the phylogenetic tree.

# Set 1: Logan clusters
logan_set <- all_clusters_summary %>%
  filter(!cluster_number %in% srr6976994_centroid_clusters) %>%
  transmute(source = "Logan", host_group = as.character(host_group))

# Set 2: NCBI reference PVs from tree
ncbi_set <- all_tree_metadata %>%
  filter(status == "ncbi") %>%
  transmute(source = "NCBI", host_group = broader_gen)

combined_fisher <- bind_rows(logan_set, ncbi_set) %>%
  filter(!is.na(host_group)) %>%
  mutate(source = factor(source, levels = c("NCBI", "Logan")))

# Only test host groups present in both sets
hosts_in_both <- combined_fisher %>%
  distinct(source, host_group) %>%
  count(host_group) %>%
  filter(n == 2) %>%
  pull(host_group)

fisher_results <- map_dfr(hosts_in_both, function(hg) {
  combined_hg <- combined_fisher %>%
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
    odds_ratio = ft$estimate,
    ci_low     = ft$conf.int[1],
    ci_high    = ft$conf.int[2]
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("\n========== Fisher's exact test: Logan vs NCBI host proportions ==========\n")
print(fisher_results, n = nrow(fisher_results))
cat("OR > 1 = Logan-enriched, OR < 1 = NCBI-enriched\n")
cat("=========================================================================\n")

# Forest plot of odds ratios
fisher_plot_data <- fisher_results %>%
  mutate(
    host_group = factor(host_group, levels = rev(host_group)),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    enrichment = if_else(odds_ratio > 1, "Logan-enriched", "NCBI-enriched")
  )

forest_plot <- ggplot(fisher_plot_data, aes(x = odds_ratio, y = host_group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, color = "grey40") +
  geom_point(aes(color = enrichment, size = logan_yes + ncbi_yes)) +
  geom_text(aes(x = ci_high, label = sig), hjust = -0.3, size = 3.5, family = "Noto Sans") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_manual(values = c("Logan-enriched" = "#4646d4", "NCBI-enriched" = "grey40"),
                     name = "Enrichment") +
  scale_size_continuous(name = "Total count", range = c(2, 6)) +
  labs(x = "Odds ratio (log scale)", y = NULL,
       title = "Host group enrichment: Logan vs NCBI clusters",
       subtitle = "Fisher's exact test, BH-adjusted (* p<0.05, ** p<0.01, *** p<0.001)") +
  theme_bw() +
  theme(text = element_text(family = "Noto Sans", size = 11),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("outputs/fisher_host_forest_plot.png", forest_plot,
       width = 22, height = 14, units = "cm", dpi = 300, bg = "white")

write.table(fisher_results, "outputs/fisher_host_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# =====================================================================
# =================== MASTER TABLE (Figure 1 data) ====================
# =====================================================================
# Combines: novelty status, host group, NCBI/Logan/nt classification,
#           library source, and geographic information — one row per sequence.

# --- Load geo data (same source files as map_gen_for_pub.R) ---
geo_biosamps <- read.csv("files/geo_data_annotation_for_all_biosamps.txt", header = FALSE)
geo_lib_biosamp <- read.table("files/hits_library_biosample.list", sep = "\t", header = TRUE, fill = TRUE)

geo_joined <- left_join(geo_biosamps, geo_lib_biosamp, by = c("V1" = "BioSample"))
geo_joined <- dplyr::select(geo_joined, V1, V2, V4, Run, LibrarySource)

# Prioritize geo source: lat_lon > geo_loc_name > country_calc > region > center
geo_joined$V2 <- factor(geo_joined$V2, levels = c(
  "lat_lon",
  "geographic location (latitude),geographic location (longitude)",
  "geo_loc_name_sam", "geo_loc_name",
  "geo_loc_name_country_calc",
  "geographic location (region and locality)",
  "region", "birth_location", "INSDC center name"
))
geo_joined <- geo_joined[order(geo_joined$V2), ]
geo_per_library <- geo_joined[!duplicated(geo_joined[, "Run"]), ]

# Parse WKB geometry to lat/lon coordinates
library(sf)
wkb_obj <- structure(as.list(geo_per_library$V4), class = "WKB")
geo_per_library$geometry <- st_as_sfc(wkb_obj, EWKB = TRUE)
geo_sf <- st_as_sf(geo_per_library, geometry = geo_per_library$geometry, crs = 4326)
geo_coords <- st_coordinates(geo_sf)
geo_per_library$longitude <- geo_coords[, 1]
geo_per_library$latitude  <- geo_coords[, 2]

geo_for_join <- geo_per_library %>%
  dplyr::select(Run, geo_source = V2, longitude, latitude) %>%
  mutate(geo_source = as.character(geo_source))

# --- Build master table ---
master_table <- all_uc_sequences %>%
  filter(!(library == "SRR6976994" & cluster_number %in% srr6976994_centroid_clusters)) %>%
  # Cluster-level: novelty status, centroid host group, cluster size
  left_join(
    all_clusters_summary %>%
      dplyr::select(cluster_number, centroid, centroid_library, host_group, status,
                    any_of("ncbi_host"), n_members),
    by = "cluster_number"
  ) %>%
  # Library-level: SRA metadata (library source)
  left_join(
    sra_metadata %>% dplyr::select(Run, LibrarySource) %>% distinct(),
    by = c("library" = "Run")
  ) %>%
  # Library-level: geo data (lat/lon)
  left_join(geo_for_join, by = c("library" = "Run")) %>%
  # Clean up library source labels
  mutate(
    lib_source = case_when(
      LibrarySource %in% c("METAGENOMIC", "GENOMIC", "TRANSCRIPTOMIC",
                           "METATRANSCRIPTOMIC", "GENOMIC SINGLE CELL",
                           "TRANSCRIPTOMIC SINGLE CELL") ~ LibrarySource,
      TRUE ~ "OTHER"
    )
  ) %>%
  dplyr::select(
    sequence, cluster_number, centroid, centroid_library, library,
    status,          # ncbi_virus / blastn / novel
    host_group,      # generalized host group (from centroid library metadata)
    any_of("ncbi_host"),  # NCBI reference host (for ncbi_virus clusters, if available)
    n_members,       # cluster size
    lib_source,      # library source category
    LibrarySource,   # raw library source
    geo_source,      # type of geo annotation used
    longitude, latitude
  )

cat("\n========== Master table summary ==========\n")
cat("Rows (sequences):", nrow(master_table), "\n")
cat("Unique clusters: ", n_distinct(master_table$cluster_number), "\n")
cat("Unique libraries:", n_distinct(master_table$library), "\n")
cat("With geo coords: ", sum(!is.na(master_table$latitude)), "\n")
cat("Status breakdown (sequences):\n")
print(table(master_table$status, useNA = "ifany"))
cat("\nStatus breakdown (clusters):\n")
print(table(master_table %>% distinct(cluster_number, .keep_all = TRUE) %>% pull(status), useNA = "ifany"))
cat("\nHost group breakdown:\n")
print(table(master_table$host_group, useNA = "ifany"))
cat("==========================================\n")

write.table(master_table, "outputs/2026.04.06.master_table_all_pvs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: outputs/master_table_all_pvs.tsv\n")


check <- new_table_expanded %>%
  filter(V1 == "M32305.1_60_1277") %>%
  distinct(member_library, member_host_group, .keep_all = TRUE) %>%
  left_join(libraries_final %>% dplyr::select(V1_clean, original_metadata),
            by = c("member_library" = "V1_clean")) %>%
  select(member_library, member_host_group, original_metadata)
table(check$original_metadata)
