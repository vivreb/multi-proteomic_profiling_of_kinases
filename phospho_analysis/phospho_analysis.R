
# Load functions and packages 

rm(list = ls())

source("phospho_functions.R")

###################################################################################
######################### Info to change for new analysis #########################
###################################################################################

# Load data

DIA_raw_CHEK1_WT <- read_protti("VR_Ex57_EGFP_Report_VR_QC (Normal).tsv") %>% 
  mutate(bait_label = "EGFP") %>% 
  rbind(read_protti("VR_Ex57_CHEK1_Report_VR_QC (Normal).tsv") %>% 
          mutate(bait_label = "CHEK1_WT"))


###################################################################################
###################################################################################
###################################################################################

interactors <- read.csv("AP-MS_interactors.csv") %>% 
  dplyr::select(-c(X))

crapome_vector_final <- read.csv("final_crapome_proteins.csv")


for(bait_gene in c("CHEK1_WT")) {
  
  # Read in DIA_raw file
  
  DIA_raw <- sprintf("DIA_raw_%s", bait_gene) %>% 
    get() %>% 
    mutate(pg_protein_accessions = ifelse(pg_protein_accessions == "P0DP23;P0DP24;P0DP25", "P0DP23", pg_protein_accessions))
  
  bait_accession <- bait_accession_list[bait_gene]
  
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
                                              proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"),
                                              cutoff_abs_diff_rt = 0.1,
                                              cutoff_fg_quantity_raw = 1024,
                                              cutoff_min_number_of_obs_per_peptide = 3)
  
  DIA_clean_normalised <- normalise_data(data = DIA_clean,
                                         imputation_method = "msimpute")
  
  DIA_clean_normalised <- bait_normalise(DIA_clean_normalised, bait_accession = bait_accession)

  
  diff_abundance_tmp <- DIA_clean_normalised %>%
    filter(!pg_protein_accessions %in% (crapome_vector_final %>% pull(entry))) %>% 
    filter(grepl("Phospho", eg_precursor_id)) %>% 
    filter(grepl(bait_gene, new_condition_id)) %>% 
    assign_missingness(
      new_sample_id,
      new_condition_id,
      eg_precursor_id,
      normalised_intensity_log2,
      ref_condition = paste(bait_gene, "DMSO", sep = "_")
    ) %>%
    calculate_diff_abundance(
      new_sample_id,
      new_condition_id,
      eg_precursor_id,
      normalised_intensity_log2,
      missingness,
      comparison,
      ref_condition = paste(bait_gene, "DMSO", sep = "_"),
      method = "moderated_t-test"
    ) %>%
    left_join((DIA_clean_normalised %>% distinct(eg_precursor_id, pg_protein_accessions)), by = "eg_precursor_id") %>% 
    left_join(uniprot, by = "pg_protein_accessions") %>% 
    dplyr::mutate(significant = abs(diff) > 0.585 & adj_pval < 0.05)
  
  assign(sprintf("diff_abundance_phospho_%s", bait_gene), diff_abundance_tmp)
  
  
  
}

peptide_label_CHEK1 <- c("CHEK1 pSer296", "CHEK1 pSer317")
names(peptide_label_CHEK1) <- c("_HIQS[Phospho (STY)]NLDFSPVNSASSEENVK_.3", "_YSSS[Phospho (STY)]QPEPR_.2")

for(bait_gene in c("CHEK1_WT")) { #bait_gene_list
  
  diff_abundance_tmp <- sprintf("diff_abundance_phospho_%s", bait_gene) %>% 
    get() %>%
    filter(grepl("probe", comparison))
  
  diff_abundance_tmp %>% write.csv(paste(getwd(), "/diff_abundance_phospho_probe_vs_DMSO_", bait_gene, ".csv", sep = ""))
  

  uniprot <- read.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = "")) %>% 
    dplyr::select(-c(X))
  
  volcano_plot <- make_volcano_plot(diff_abundance_tmp, peptide_labels = peptide_label_CHEK1)

  pdf(paste(getwd(), "/volcano_plot_phospho_probe_vs_DMSO_", bait_gene, ".pdf", sep = ""),
      width = 4,
      height = 4,
      pointsize = 20)
  print(volcano_plot)
  dev.off()
  
  
}
