###look at all 21k
all_21323_libraries_metadata <- read.table("files/all_21323_libraries_metadata.list", fill = TRUE, quote = "", sep = ",")

# Classify V3 (organism/host) into generalized host groups
all_21323_libraries_metadata <- all_21323_libraries_metadata %>%
  mutate(V3_clean = str_trim(gsub('"', '', V3))) %>%
  mutate(host_group = case_when(
    # Human
    str_detect(V3_clean, regex("^Homo sapiens|^human |^Human p|human feces|human gut|human skin|human vaginal|human oral|human nasopharyngeal|human reproductive|human metagenome|human viral|human lung|human urinary|human saliva|human blood|human milk|uncultured human", ignore_case = TRUE)) ~ "Human",
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

# Step 1: Convert all BioSample accessions to NCBI numeric UIDs via esearch
# This handles SAMD/SAME/SAMN cross-references reliably
cat("Converting BioSample accessions to NCBI UIDs...\n")
acc_to_uid <- data.frame(BioSample = character(), uid = character(), stringsAsFactors = FALSE)

for (i in seq_len(ceiling(length(biosample_ids) / batch_size))) {
  start_idx <- (i - 1) * batch_size + 1
  batch <- biosample_ids[start_idx:min(start_idx + batch_size - 1, length(biosample_ids))]

  tryCatch({
    # Search each accession individually to maintain mapping
    for (acc in batch) {
      search_res <- entrez_search(db = "biosample", term = acc)
      if (length(search_res$ids) > 0) {
        acc_to_uid <- rbind(acc_to_uid,
                            data.frame(BioSample = acc, uid = search_res$ids[1],
                                       stringsAsFactors = FALSE))
      }
    }
    Sys.sleep(0.1)
  }, error = function(e) {
    cat("  Warning: search batch failed -", conditionMessage(e), "\n")
  })
  if (i %% 5 == 0) cat("  Searched", min(start_idx + batch_size - 1, length(biosample_ids)),
                         "of", length(biosample_ids), "accessions\n")
}

cat("Mapped", nrow(acc_to_uid), "of", length(biosample_ids), "accessions to UIDs\n")

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
      original_accs <- acc_to_uid$BioSample[acc_to_uid$uid == uid_returned]
      if (length(original_accs) == 0) original_accs <- acc_returned

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
    str_detect(biosample_organism, regex("homo sapiens|human", ignore_case = TRUE)) ~ "Human",
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
