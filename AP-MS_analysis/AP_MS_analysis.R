
# Load functions and packages from Z:/Viviane/Manuscripts/supplement/github/AP_MS_functions.R

rm(list = ls())

source("AP_MS_functions.R")

###################################################################################
######################### Info to change for new analysis #########################
###################################################################################

# Load data

DIA_raw_EGFP <- read_protti("20241113_145624_VR_Ex75_EGFP_only_CAMKK2_wt_fasta_Report.tsv") %>% 
  mutate(bait_label = "EGFP")

DIA_raw_CAMKK2_WT <- DIA_raw_EGFP %>% 
  rbind(read_protti("20241113_142710_VR_Ex75_CAMKK2_wt_Report.tsv") %>% 
          mutate(bait_label = "CAMKK2_WT"))

DIA_raw_CAMKK2_R311C <- DIA_raw_EGFP %>% 
  rbind(read_protti("20241108_072530_VR_Ex75_CAMKK2_R311C_Report.tsv") %>% 
          mutate(bait_label = "CAMKK2_R311C"))

DIA_raw_CAMKK2_D312A <- DIA_raw_EGFP %>%
  rbind(read_protti("20241114_152135_VR_Ex75_CAMKK2_D312A_rem_Report.tsv") %>%
          mutate(bait_label = "CAMKK2_D312A"))

DIA_raw_CAMKK2_T483D <- DIA_raw_EGFP %>%
  rbind(read_protti("20241115_083241_VR_Ex75_CAMKK2_T483D_Report.tsv") %>%
          mutate(bait_label = "CAMKK2_T483D"))

DIA_raw_CHEK1_WT <- read_protti("VR_Ex78_EGFP_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "EGFP") %>%
  rbind(read_protti("VR_Ex78_CHEK1_WT_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "CHEK1_WT"))

DIA_raw_CHEK1_D130A <- read_protti("VR_Ex78_EGFP_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "EGFP") %>%
  rbind(read_protti("VR_Ex78_CHEK1_D130A_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "CHEK1_D130A"))

DIA_raw_CHEK1_L449R <- read_protti("VR_Ex78_EGFP_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "EGFP") %>%
  rbind(read_protti("VR_Ex78_CHEK1_L449R_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "CHEK1_L449R"))

DIA_raw_PRKCA_WT <- read_protti("VR_Ex81_EGFP_Report_VR_QC (Normal).tsv") %>% 
  mutate(bait_label = "EGFP") %>% 
  rbind(read_protti("VR_Ex81_PRKCA_WT_Report_VR_QC (Normal).tsv") %>% 
          mutate(bait_label = "PRKCA_WT"))

DIA_raw_Extra <- read_protti("VR_Ex70_EGFP_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "EGFP") %>% 
  rbind(read_protti("VR_Ex70_bait_extra_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "Extra"))


# Initialize list to iterate over

bait_gene_list <- c("CAMKK2_WT",
                    "CHEK1_WT",
                    "PRKCA_WT", 
                    "Extra")

bait_accession_list <- c("Q96RR4",
                         "O14757", 
                         "P17252",
                         "O14936")
names(bait_accession_list) <- c(bait_gene_list)


biogrid_interactors_CAMKK2 <-
  read_delim(
    "230926_CAMKK2_BioGrid_data_tab3.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  janitor::clean_names()

biogrid_interactors_CHEK1 <-
  read_delim(
    "241205_CHEK1_BioGrid_data_tab3.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  janitor::clean_names()

biogrid_interactors_PRKCA <-
  read_delim(
    "240516_PRKCA_BioGrid_data_tab3.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  janitor::clean_names()


biogrid_methods <- c("Affinity Capture-Western", "Affinity Capture-MS", "Two-hybrid", "Reconstituted Complex", 
                     "Affinity Capture-Luminescence", "Co-purification", "PCA")


###################################################################################
###################################################################################
###################################################################################



# Determine interactors using adapted WD score


prot_intensity_data_all_baits <- data.frame()

for(bait_gene in c(bait_gene_list)) {
  
  # Read in DIA_raw file
  
  DIA_raw <- sprintf("DIA_raw_%s", bait_gene) %>% 
    get() %>% 
    filter(r_condition %in% c("probe", "DMSO")) %>% 
    mutate(pg_protein_accessions = ifelse(pg_protein_accessions == "P0DP23;P0DP24;P0DP25", "P0DP23", pg_protein_accessions))
  
  bait_accession <- bait_accession_list[bait_gene]
  
  if(bait_gene != "Extra"){
  biogrid_interactors <- sprintf("biogrid_interactors_%s", str_replace(bait_gene, "\\_.*", "")) %>% 
    get()
  
  biogrid_interactors_vector <- biogrid_interactors %>% 
    filter(experimental_system %in% biogrid_methods) %>% 
    dplyr::select(c(swiss_prot_accessions_interactor_a, official_symbol_interactor_a)) %>% 
    dplyr::rename(protein_id = swiss_prot_accessions_interactor_a,
                  gene = official_symbol_interactor_a) %>% 
    rbind(biogrid_interactors %>% 
            filter(experimental_system %in% biogrid_methods) %>% 
            dplyr::select(c(swiss_prot_accessions_interactor_b, official_symbol_interactor_b)) %>% 
            dplyr::rename(protein_id = swiss_prot_accessions_interactor_b,
                          gene = official_symbol_interactor_b)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull(protein_id)
  } else{biogrid_interactors_vector = c() }

  
  # Filter, normalise and impute data
  
  unis <- DIA_raw %>%
    dplyr::pull(pg_protein_accessions) %>%
    strsplit(";") %>%
    unlist() %>%
    unique()
  
  
  uniprot <-
    read.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = "")) %>%
    dplyr::select(-c(X))
    
  DIA_raw_uniprot <- DIA_raw %>% 
    left_join(uniprot, by = c("pg_protein_accessions")) %>% 
    find_peptide(protein_sequence, pep_stripped_sequence) %>%
    assign_peptide_type(aa_before, last_aa)
  
  DIA_clean <- filter_data_and_log2_transform(data = DIA_raw_uniprot,
                                              proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"))
  
  DIA_clean_normalised <- normalise_data(data = DIA_clean)
  
  
  protein_abundance <- calculate_protein_abundance(
    data = DIA_clean_normalised,
    sample = new_sample_id,
    protein_id = pg_protein_accessions,
    precursor = eg_precursor_id,
    peptide = pep_stripped_sequence,
    intensity_log2 = normalised_intensity_log2,
    method = "iq",
    for_plot = FALSE,
    retain_columns = c(new_condition_id)
  )
  
  
  prot_intensity_data <- protein_abundance %>% 
    dplyr::rename(protein_intensity_log2 = normalised_intensity_log2) %>% 
    left_join(DIA_clean_normalised, by = c("new_condition_id", "new_sample_id", "pg_protein_accessions")) %>% 
    mutate(is_interactor = ifelse(pg_protein_accessions %in% biogrid_interactors_vector, TRUE, FALSE)) %>% 
    mutate(kinase_label = bait_gene)
  
  
  prot_intensity_data_all_baits <- prot_intensity_data_all_baits %>% 
    rbind(prot_intensity_data)
  
}

prot_intensity_data_all_baits <- prot_intensity_data_all_baits %>% 
  mutate(bait_label = ifelse(bait_label == "EGFP", paste(kinase_label, bait_label, sep = "_"), bait_label)) %>%
  mutate(new_condition_id = ifelse(grepl("EGFP", new_condition_id), paste(kinase_label, new_condition_id, sep = "_"), new_condition_id))%>%
  mutate(new_sample_id = ifelse(grepl("EGFP", new_sample_id), paste(kinase_label, new_sample_id, sep = "_"), new_sample_id)) %>% 
  mutate(kinase_label = ifelse(grepl("EGFP", new_sample_id), "EGFP", kinase_label))

proteins_in_condition_count <- prot_intensity_data_all_baits %>% 
  filter(kinase_label %in% c("CAMKK2_WT", "CHEK1_WT", "PRKCA_WT", "Extra", "EGFP")) %>% 
  dplyr::distinct(kinase_label, bait_label, new_condition_id, new_sample_id, pg_protein_accessions, protein_intensity_log2, is_interactor) %>%
  group_by(new_condition_id, pg_protein_accessions) %>% 
  mutate(condition_sum_intensity_for_protein = sum(protein_intensity_log2)) %>% 
  ungroup() %>% 
  distinct(kinase_label, new_condition_id, pg_protein_accessions, condition_sum_intensity_for_protein) %>% 
  group_by(pg_protein_accessions) %>% 
  mutate(condition_sum_intensity_for_protein_EGFP = ifelse(!is.null(condition_sum_intensity_for_protein[which(kinase_label == "EGFP")]), condition_sum_intensity_for_protein[which(kinase_label == "EGFP")], 0)) %>%
  mutate(condition_sum_intensity_for_protein_EGFP = ifelse(is.na(condition_sum_intensity_for_protein_EGFP), 0, condition_sum_intensity_for_protein_EGFP)) %>% 
  ungroup() %>% 
  mutate(condition_sum_intensity_for_protein_final = condition_sum_intensity_for_protein + condition_sum_intensity_for_protein_EGFP) %>% 
  filter(kinase_label != "EGFP")

wd_score_per_bait_and_treatment <- prot_intensity_data_all_baits %>% 
  filter(kinase_label %in% c("CAMKK2_WT", "CHEK1_WT", "PRKCA_WT")) %>% 
  dplyr::distinct(kinase_label, bait_label, new_condition_id, new_sample_id, pg_protein_accessions, protein_intensity_log2, is_interactor) %>%
  left_join(proteins_in_condition_count, by = c("pg_protein_accessions", "kinase_label", "new_condition_id")) %>% 
  mutate(D_score_R = sqrt((708/condition_sum_intensity_for_protein_final) ^ (condition_sum_intensity_for_protein / 170))) %>% #708 is the max total intensity, 170 is the max in a single condition
  mutate(D_score_scaled = ifelse(bait_label != "EGFP" & is_interactor == TRUE, D_score_R * 1.1, D_score_R)) %>% 
  group_by(kinase_label) %>% 
  mutate(is_bait = ifelse(bait_accession_list[[kinase_label[1]]] == pg_protein_accessions, TRUE, FALSE)) %>% 
  ungroup()


roc_1 <- wd_score_per_bait_and_treatment %>% 
  filter(!grepl("EGFP", bait_label)) %>% 
  group_by(kinase_label, pg_protein_accessions) %>% 
  mutate(D_score_R = max(D_score_R)) %>% 
  ungroup() %>% 
  roc(is_interactor, D_score_scaled)

print(plot.roc(roc_1, print.auc = TRUE))

# Make ROC plot

# pdf(paste("roc_curve_mod_wd_score.pdf", sep = ""), 
#     width = 6,
#     height = 6,
#     pointsize = 20)
# print(plot.roc(roc_1, print.auc = TRUE))
# dev.off()

ranked_protein_kinases <- wd_score_per_bait_and_treatment %>%
  mutate(D_score_R = D_score_scaled) %>% 
  filter(!grepl("EGFP", bait_label)) %>% 
  #filter(bait_label == "CHEK1_WT") %>% 
  distinct(bait_label, pg_protein_accessions, D_score_R, is_interactor) %>% 
  arrange(D_score_R, desc = FALSE) %>% 
  mutate(rank = row_number()) %>% 
  ggplot(aes(x = rank, y = D_score_R, color = is_interactor)) +
  geom_point() +
  theme_bw() +
  labs(color = "Known interactor", y = "Modified D-score", x = "Rank") +
  geom_hline(yintercept = 1.425, linetype = "dashed") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

ranked_protein_kinases

# Make ranked kinase plot

# pdf(paste("mod_d_score_ranked_proteins_kinases.pdf", sep = ""), 
#     width = 4,
#     height = 3,
#     pointsize = 20)
# ranked_protein_kinases
# dev.off()


unis <- wd_score_per_bait_and_treatment %>%
  dplyr::pull(pg_protein_accessions) %>%
  strsplit(";") %>%
  unlist() %>%
  unique()

uniprot <-
  fetch_uniprot(
    unis,
    columns = c(
      "protein_name",
      "gene_primary"
    ),
    batchsize = 100
  ) %>%
  dplyr::rename(
    pg_protein_accessions = accession
  )

interactors <- wd_score_per_bait_and_treatment %>% 
  filter(!grepl("EGFP", bait_label)) %>% 
  filter(D_score_scaled > 1.425) %>% 
  distinct(bait_label, pg_protein_accessions, is_interactor) %>% 
  group_by(bait_label) %>% 
  mutate(no_interactors = n()) %>% 
  ungroup() %>% 
  left_join(uniprot, by = "pg_protein_accessions")

interactors %>% write.csv(paste(getwd(), "/interactors_all_baits_per_treatment_condition.csv", sep = ""))

# Use determined interactors to find probe-induced changes

AP_eval_list <- c("CAMKK2_WT", "CHEK1_WT", "PRKCA_WT")


for(bait_gene in c(AP_eval_list)) {
  
  # Read in DIA_raw file
  
  DIA_raw <- sprintf("DIA_raw_%s", bait_gene) %>% 
    get() %>% 
    mutate(pg_protein_accessions = ifelse(pg_protein_accessions == "P0DP23;P0DP24;P0DP25", "P0DP23", pg_protein_accessions))
  
  bait_accession <- bait_accession_list[bait_gene]
  
  biogrid_interactors <- sprintf("biogrid_interactors_%s", str_replace(bait_gene, "\\_.*", "")) %>% 
    get()
  
  biogrid_interactors_vector <- biogrid_interactors %>% 
    filter(experimental_system %in% biogrid_methods) %>% 
    dplyr::select(c(swiss_prot_accessions_interactor_a, official_symbol_interactor_a)) %>% 
    dplyr::rename(protein_id = swiss_prot_accessions_interactor_a,
                  gene = official_symbol_interactor_a) %>% 
    rbind(biogrid_interactors %>% 
            dplyr::select(c(swiss_prot_accessions_interactor_b, official_symbol_interactor_b)) %>% 
            dplyr::rename(protein_id = swiss_prot_accessions_interactor_b,
                          gene = official_symbol_interactor_b)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull(protein_id)
  
  # Filter, normalise and impute data
  
  unis <- DIA_raw %>%
    dplyr::pull(pg_protein_accessions) %>%
    strsplit(";") %>%
    unlist() %>%
    unique()
  
  uniprot <- read.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = "")) %>% 
      dplyr::select(-c(X))
  
  DIA_raw_uniprot <- DIA_raw %>% 
    left_join(uniprot, by = c("pg_protein_accessions")) %>% 
    find_peptide(protein_sequence, pep_stripped_sequence) %>%
    assign_peptide_type(aa_before, last_aa)
  
  DIA_clean <- filter_data_and_log2_transform(data = DIA_raw_uniprot,
                                              proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"))
  
  DIA_clean_normalised <- normalise_data(data = DIA_clean,
                                         imputation_method = "msimpute")
  
  DIA_clean_normalised <- bait_normalise(DIA_clean_normalised, bait_accession = bait_accession)
  
  
  protein_abundance <- calculate_protein_abundance(
    data = DIA_clean_normalised,
    sample = new_sample_id,
    protein_id = pg_protein_accessions,
    precursor = eg_precursor_id,
    peptide = pep_stripped_sequence,
    intensity_log2 = normalised_intensity_log2,
    method = "iq",
    for_plot = FALSE,
    retain_columns = c(new_condition_id)
  )
  
  diff_abundance_tmp <- protein_abundance %>%
    filter(grepl(bait_gene, new_condition_id)) %>% 
    filter(pg_protein_accessions %in% (interactors %>% filter(bait_label == bait_gene) %>% pull(pg_protein_accessions))) %>% 
    assign_missingness(
      new_sample_id,
      new_condition_id,
      pg_protein_accessions,
      normalised_intensity_log2,
      ref_condition = paste(bait_gene, "DMSO", sep = "_")
    ) %>%
    calculate_diff_abundance(
      new_sample_id,
      new_condition_id,
      pg_protein_accessions,
      normalised_intensity_log2,
      missingness,
      comparison,
      ref_condition = paste(bait_gene, "DMSO", sep = "_"),
      method = "moderated_t-test"
    ) %>%
    left_join(uniprot, by = "pg_protein_accessions") %>% 
    dplyr::mutate(significant = abs(diff) > 1 & adj_pval < 0.05)
  
  assign(sprintf("diff_abundance_%s", bait_gene), diff_abundance_tmp)
  
  
}

for(bait_gene in c(AP_eval_list)) {
  
  diff_abundance_tmp <- sprintf("diff_abundance_%s", bait_gene) %>% 
    get() %>%
    filter(grepl("probe", comparison)) %>% 
    dplyr::mutate(significant = abs(diff) > 1 & adj_pval < 0.05)
  
  
  # diff_abundance_tmp %>% write.csv(paste(getwd(), "/diff_abundance_probe_vs_DMSO_", bait_gene, ".csv", sep = ""))
  
  diff_abundance_tmp <- diff_abundance_tmp %>% 
    left_join((interactors %>% 
                 filter(bait_label == bait_gene) %>% 
                 distinct(pg_protein_accessions, is_interactor)), 
              by = "pg_protein_accessions")
  
  uniprot <-
    fetch_uniprot(
      uniprot_ids = (diff_abundance_tmp %>% pull(pg_protein_accessions)),
      columns = c(
        "protein_name",
        "length",
        "sequence",
        "protein_name",
        "gene_primary"),
      batchsize = 100
    ) %>%
    dplyr::rename(
      protein_sequence = sequence,
      length_protein = length,
      pg_protein_accessions = accession
    )
  
  volcano_plot <- make_volcano_plot(diff_abundance_tmp)
  
  
  # Make volcano plot
  
  # pdf(paste(getwd(), "/FINAL_volcano_plot_probe_vs_DMSO_", bait_gene, ".pdf", sep = ""),
  #     width = 4,
  #     height = 4,
  #     pointsize = 20)
  # print(volcano_plot)
  # dev.off()

  
}


AP_eval_list_CAMKK2 <- c("CAMKK2_WT", "CAMKK2_R311C", "CAMKK2_D312A", "CAMKK2_T483D")


for(bait_gene in c(AP_eval_list_CAMKK2)) {
  
  # Read in DIA_raw file
  
  DIA_raw <- sprintf("DIA_raw_%s", bait_gene) %>% 
    get() %>% 
    mutate(pg_protein_accessions = ifelse(pg_protein_accessions == "P0DP23;P0DP24;P0DP25", "P0DP23", pg_protein_accessions))
  
  bait_accession <- bait_accession_list[bait_gene]
  
  biogrid_interactors <- sprintf("biogrid_interactors_%s", str_replace(bait_gene, "\\_.*", "")) %>% 
    get()
  
  biogrid_interactors_vector <- biogrid_interactors %>% 
    filter(experimental_system %in% biogrid_methods) %>% 
    dplyr::select(c(swiss_prot_accessions_interactor_a, official_symbol_interactor_a)) %>% 
    dplyr::rename(protein_id = swiss_prot_accessions_interactor_a,
                  gene = official_symbol_interactor_a) %>% 
    rbind(biogrid_interactors %>% 
            dplyr::select(c(swiss_prot_accessions_interactor_b, official_symbol_interactor_b)) %>% 
            dplyr::rename(protein_id = swiss_prot_accessions_interactor_b,
                          gene = official_symbol_interactor_b)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull(protein_id)
  
  # Filter, normalise and impute data
  
  unis <- DIA_raw %>%
    dplyr::pull(pg_protein_accessions) %>%
    strsplit(";") %>%
    unlist() %>%
    unique()
  
  uniprot <- read.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = "")) %>% 
      dplyr::select(-c(X))

  
  DIA_raw_uniprot <- DIA_raw %>% 
    left_join(uniprot, by = c("pg_protein_accessions")) %>% 
    find_peptide(protein_sequence, pep_stripped_sequence) %>%
    assign_peptide_type(aa_before, last_aa)
  
  DIA_clean <- filter_data_and_log2_transform(data = DIA_raw_uniprot,
                                              proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"))
  
  DIA_clean_normalised <- normalise_data(data = DIA_clean,
                                         imputation_method = "msimpute")
  
  DIA_clean_normalised <- bait_normalise(DIA_clean_normalised, bait_accession = bait_accession)
  
  
  protein_abundance <- calculate_protein_abundance(
    data = DIA_clean_normalised,
    sample = new_sample_id,
    protein_id = pg_protein_accessions,
    precursor = eg_precursor_id,
    peptide = pep_stripped_sequence,
    intensity_log2 = normalised_intensity_log2,
    method = "iq",
    for_plot = FALSE,
    retain_columns = c(new_condition_id)
  )
  
  assign(sprintf("AP_protein_abundance_%s", bait_gene), protein_abundance)
  
  
  diff_abundance_tmp <- protein_abundance %>%
    filter(grepl(bait_gene, new_condition_id)) %>% 
    filter(pg_protein_accessions %in% (interactors %>% filter(bait_label == "CAMKK2_WT") %>% pull(pg_protein_accessions))) %>% 
    assign_missingness(
      new_sample_id,
      new_condition_id,
      pg_protein_accessions,
      normalised_intensity_log2,
      ref_condition = paste(bait_gene, "DMSO", sep = "_")
    ) %>%
    calculate_diff_abundance(
      new_sample_id,
      new_condition_id,
      pg_protein_accessions,
      normalised_intensity_log2,
      missingness,
      comparison,
      ref_condition = paste(bait_gene, "DMSO", sep = "_"),
      method = "moderated_t-test"
    ) %>%
    left_join(uniprot, by = "pg_protein_accessions") %>% 
    dplyr::mutate(significant = abs(diff) > 1 & adj_pval < 0.05)
  
  assign(sprintf("diff_abundance_%s", bait_gene), diff_abundance_tmp)
  
  #diff_abundance_tmp %>% write.csv(paste(getwd(), "/diff_abundance_probe_vs_DMSO_", bait_gene, "_mut.csv", sep = ""))
  
}

all_CAMKK2_changes <- diff_abundance_CAMKK2_WT %>%
  rbind(diff_abundance_CAMKK2_D312A) %>%
  rbind(diff_abundance_CAMKK2_R311C) %>%
  rbind(diff_abundance_CAMKK2_T483D) %>% 
  mutate(current_bait_id = str_extract(comparison, "CAMKK2_[A-Z0-9]*"))


levels_vector = c("WT", "D312A", "T483D", "R311C")
significance_vector_PRKAA1 <- c(paste0("p == 1.1", "%*%", "10^{-9}"), paste0("p == 2.9", "%*%", "10^{-4}"), paste0("p == 5.5", "%*%", "10^{-6}"), paste0("ND"))
color = c("#7E69AF", "#7E69AF", "#7E69AF", "#7E69AF")
all_CAMKK2_data_sem <- read.csv("all_CAMKK2_data_sem.csv")


CAMKK2_PRKAA1_across_mutants_plot <- plot_individual_protein_change_CAMKK2(all_CAMKK2_changes, precursor = "Q13131", grouping = pg_protein_accessions, 
                                       gene_id = "PRKAA1", levels_vector = levels_vector, label_vector = significance_vector_PRKAA1,
                                       color = color, sem_data = all_CAMKK2_data_sem)

# pdf(paste(getwd(), "/mutant_interaction_plot_CAMKK2_PRKAA1_errorbars.pdf", sep = ""), 
#     width = 5,
#     height = 3.5,
#     pointsize = 20)
# print(CAMKK2_PRKAA1_across_mutants_plot)
# dev.off()





AP_eval_list_CHEK1 <- c("CHEK1_WT", "CHEK1_D130A", "CHEK1_L449R")

for(bait_gene in c(AP_eval_list_CHEK1)) {
  
  # Read in DIA_raw file
  
  DIA_raw <- sprintf("DIA_raw_%s", bait_gene) %>% 
    get() %>% 
    mutate(pg_protein_accessions = ifelse(pg_protein_accessions == "P0DP23;P0DP24;P0DP25", "P0DP23", pg_protein_accessions))
  
  bait_accession <- bait_accession_list[bait_gene]
  
  biogrid_interactors <- sprintf("biogrid_interactors_%s", str_replace(bait_gene, "\\_.*", "")) %>% 
    get()
  
  biogrid_interactors_vector <- biogrid_interactors %>% 
    filter(experimental_system %in% biogrid_methods) %>% 
    dplyr::select(c(swiss_prot_accessions_interactor_a, official_symbol_interactor_a)) %>% 
    dplyr::rename(protein_id = swiss_prot_accessions_interactor_a,
                  gene = official_symbol_interactor_a) %>% 
    rbind(biogrid_interactors %>% 
            dplyr::select(c(swiss_prot_accessions_interactor_b, official_symbol_interactor_b)) %>% 
            dplyr::rename(protein_id = swiss_prot_accessions_interactor_b,
                          gene = official_symbol_interactor_b)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull(protein_id)
  
  # Filter, normalise and impute data
  
  unis <- DIA_raw %>%
    dplyr::pull(pg_protein_accessions) %>%
    strsplit(";") %>%
    unlist() %>%
    unique()
  
  uniprot <- read.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = "")) %>% 
      dplyr::select(-c(X))

  
  DIA_raw_uniprot <- DIA_raw %>% 
    left_join(uniprot, by = c("pg_protein_accessions")) %>% 
    find_peptide(protein_sequence, pep_stripped_sequence) %>%
    assign_peptide_type(aa_before, last_aa)
  
  DIA_clean <- filter_data_and_log2_transform(data = DIA_raw_uniprot,
                                              proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"))
  
  DIA_clean_normalised <- normalise_data(data = DIA_clean,
                                         imputation_method = "msimpute")
  
  DIA_clean_normalised <- bait_normalise(DIA_clean_normalised, bait_accession = bait_accession)
  
  
  protein_abundance <- calculate_protein_abundance(
    data = DIA_clean_normalised,
    sample = new_sample_id,
    protein_id = pg_protein_accessions,
    precursor = eg_precursor_id,
    peptide = pep_stripped_sequence,
    intensity_log2 = normalised_intensity_log2,
    method = "iq",
    for_plot = FALSE,
    retain_columns = c(new_condition_id)
  )
  
  assign(sprintf("AP_protein_abundance_%s", bait_gene), protein_abundance)
  
  diff_abundance_tmp <- protein_abundance %>%
    filter(grepl(bait_gene, new_condition_id)) %>% 
    filter(pg_protein_accessions %in% (interactors %>% filter(bait_label == "CHEK1_WT") %>% pull(pg_protein_accessions))) %>% 
    assign_missingness(
      new_sample_id,
      new_condition_id,
      pg_protein_accessions,
      normalised_intensity_log2,
      ref_condition = paste(bait_gene, "DMSO", sep = "_")
    ) %>%
    calculate_diff_abundance(
      new_sample_id,
      new_condition_id,
      pg_protein_accessions,
      normalised_intensity_log2,
      missingness,
      comparison,
      ref_condition = paste(bait_gene, "DMSO", sep = "_"),
      method = "moderated_t-test"
    ) %>%
    left_join(uniprot, by = "pg_protein_accessions") %>% 
    dplyr::mutate(significant = abs(diff) > 1 & adj_pval < 0.05)
  
  diff_abundance_tmp %>% write.csv(paste(getwd(), "/diff_abundance_probe_vs_DMSO_", bait_gene, "_mut.csv", sep = ""))
  
  assign(sprintf("diff_abundance_%s", bait_gene), diff_abundance_tmp)
  
  
}

all_CHEK1_changes <- diff_abundance_CHEK1_WT %>%
  left_join((uniprot %>% dplyr::select(-c(input_id, length_protein, protein_sequence, protein_name))), by = "pg_protein_accessions") %>% 
  rbind(diff_abundance_CHEK1_D130A) %>%
  rbind(diff_abundance_CHEK1_L449R) %>% 
  mutate(current_bait_id = str_extract(comparison, "CHEK1_[A-Z0-9]*"))



levels_vector = c("WT", "D130A", "L449R")
significance_vector_CLPB <- c(paste0("p == 6.8", "%*%", "10^{-4}"), paste0("p == 5.8", "%*%", "10^{-3}"), paste0("NS"))
color = c("#8bbd26", "#8bbd26", "#c5aef2")
all_CHEK1_data_sem <- read.csv("all_CHEK1_data_sem.csv")

CHEK1_CLPB_across_mutants_plot <- plot_individual_protein_change_CHEK1(all_CHEK1_changes, precursor = "Q9H078", grouping = pg_protein_accessions, 
                                       gene_id = "CLPB", levels_vector = levels_vector, label_vector = significance_vector_CLPB,
                                       color = color, sem_data = all_CHEK1_data_sem)


# pdf(paste(getwd(), "/mutant_interaction_plot_CHEK1_CLPB.pdf", sep = ""), 
#     width = 4,
#     height = 3.5,
#     pointsize = 20)
# print(CHEK1_CLPB_across_mutants_plot)
# dev.off()
