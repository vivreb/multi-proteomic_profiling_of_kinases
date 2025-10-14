rm(list = ls())

source("BioID_functions.R")

###################################################################################
######################### Info to change for new analysis #########################
###################################################################################

# Load data

DIA_raw_EGFP <- read_protti("VR_Ex82_EGFP_for_PRKCA_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "EGFP") %>% 
  dplyr::select(-c(pep_digest_type_trypsin_p))

DIA_raw_PRKCA_WT <- read_protti("VR_Ex82_PRKCA_WT_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_WT") %>%
  rbind(DIA_raw_EGFP)

DIA_raw_PRKCA_R22A_A25E <- read_protti("VR_Ex82_PRKCA_R22A_A25E_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_R22A_A25E") %>%
  rbind(DIA_raw_EGFP)

DIA_raw_PRKCA_D463A <- read_protti("VR_Ex82_PRKCA_D463A_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_D463A") %>%
  rbind(DIA_raw_EGFP)

DIA_raw_CAMKK2_Nt <- read_protti("SL_Ex1_Nt_CAMKK2_Report_VR_QC (Normal).tsv") %>% 
  mutate(bait_label = "CAMKK2_Nt") %>% 
  rbind(read_protti("SL_Ex1_EGFP_Report_VR_QC (Normal).tsv") %>% 
          mutate(bait_label = "EGFP"))

DIA_raw_PRKAA1 <- read_protti("SL_Ex1_PRKAA1_Report_VR_QC (Normal).tsv") %>% 
  mutate(bait_label = "PRKAA1") %>% 
  rbind(read_protti("SL_Ex1_EGFP_Report_VR_QC (Normal).tsv") %>% 
          mutate(bait_label = "EGFP"))

DIA_raw_CHEK1 <- read_protti("VR_Ex90_CHEK1_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "CHEK1") %>%
  rbind(read_protti("VR_Ex90_EGFP_for_CHEK1_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "EGFP"))

label_list = list()
label_list[["CAMKK2_Nt"]] <- c("Q13131", "P54646", "Q9Y478", "O43741", "P54619")
label_list[["PRKAA1"]] <- c("Q96RR4") 
label_list[["CHEK1"]] = c("Q13535", "Q92547", "Q7L590")
label_list[["PRKCA_WT"]] = c("P05556", "P23229", "P98172", "Q14126", "P43121", "P17301", "P78310", "Q16625", "Q8NFZ8")
label_list[["PRKCA_R22A_A25E"]] = c("P05556", "P23229")
label_list[["PRKCA_D463A"]] = c("P05556", "P23229")


# Initialize list to iterate over

bait_gene_list <- c("CAMKK2_Nt",
                    "PRKAA1",
                    "CHEK1",
                    "PRKCA_WT",
                    "PRKCA_R22A_A25E",
                    "PRKCA_D463A")

bait_accession_list <- c("Q96RR4",
                         "Q13131",
                         "O14757",
                         "P17252", "P17252", "P17252")
names(bait_accession_list) <- c(bait_gene_list)


bait_gene_list <- c("PRKAA1")

bait_accession_list <- c("Q13131")

names(bait_accession_list) <- c(bait_gene_list)


biogrid_interactors_CAMKK2 <-
  read_delim(
    "Z:/Viviane/Experiments/VR_Ex42/BIOGRID-GENE-115889-4.4.225.tab3.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  janitor::clean_names()

biogrid_interactors_CHEK1 <-
  read_delim(
    "Z:/Viviane/Experiments/VR_Ex57/BIOGRID-GENE-107536-4.4.233.tab3.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  janitor::clean_names()

biogrid_interactors_PRKCA <-
  read_delim(
    "Z:/Viviane/Experiments/VR_Ex59/BIOGRID-GENE-111564-4.4.233.tab3.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  janitor::clean_names()

biogrid_interactors_PRKAA1 <-
  read_delim(
    "Z:/Viviane/Experiments/SL_Ex1/BIOGRID-GENE-111549-4.4.241.tab3.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  janitor::clean_names()


biogrid_methods <- c("Affinity Capture-Western", "Affinity Capture-MS", "Two-hybrid", "Reconstituted Complex", 
                     "Affinity Capture-Luminescence", "Co-purification", "PCA")

biogrid_interactors_vector <- biogrid_interactors_PRKCA %>% 
  filter(experimental_system %in% biogrid_methods) %>% 
  dplyr::select(c(swiss_prot_accessions_interactor_a, official_symbol_interactor_a)) %>% 
  dplyr::rename(protein_id = swiss_prot_accessions_interactor_a,
                gene = official_symbol_interactor_a) %>% 
  rbind(biogrid_interactors_PRKCA %>% 
          filter(experimental_system %in% biogrid_methods) %>% 
          dplyr::select(c(swiss_prot_accessions_interactor_b, official_symbol_interactor_b)) %>% 
          dplyr::rename(protein_id = swiss_prot_accessions_interactor_b,
                        gene = official_symbol_interactor_b)) %>% 
  dplyr::distinct() %>% 
  dplyr::pull(protein_id)

crapome_vector_final <- read.csv("final_crapome_proteins.csv")

###################################################################################
###################################################################################
###################################################################################


for(bait_gene in c(bait_gene_list)) {
  
    # Read in DIA_raw file
    
    DIA_raw <- sprintf("DIA_raw_%s", bait_gene) %>% 
      get() %>% 
      mutate(pg_protein_accessions = ifelse(pg_protein_accessions == "P0DP23;P0DP24;P0DP25", "P0DP23", pg_protein_accessions))
    
    bait_accession <- bait_accession_list[bait_gene]
    
    # Filter, normalise and impute data

    if(file.exists(paste(getwd(), "/protein_abundance_imputed_", bait_gene, ".csv", sep = ""))){
      protein_abundance_imputed <- read.csv(paste(getwd(), "/protein_abundance_imputed_", bait_gene, ".csv", sep = "")) %>%
        dplyr::select(-c(X))

      print("Imputed data loaded from saved csv")

    } else {


      
      DIA_clean <- filter_data_and_log2_transform(data = DIA_raw,
                                                  proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"))
      
      DIA_clean_normalised <- normalise_data(data = DIA_clean, imputation_method = "msimpute")
      
      DIA_clean_normalised %>% write.csv(paste(getwd(), "/DIA_clean_imputed_res_only_", bait_gene, ".csv", sep = ""))
      
      
      protein_abundance_imputed <- calculate_protein_abundance(
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
      
      protein_abundance_imputed %>% write.csv(paste(getwd(), "/protein_abundance_imputed_", bait_gene, ".csv", sep = ""))
      
    }
    
    
    unis <- protein_abundance_imputed %>%
      dplyr::pull(pg_protein_accessions) %>%
      strsplit(";") %>%
      unlist() %>%
      unique()
    
    if(file.exists(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = ""))){
      uniprot <- read.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = "")) %>% 
        dplyr::select(-c(X))
      
      print("Uniprot data loaded from saved csv")
      
      
    } else{
      uniprot <-
        fetch_uniprot(
          unis,
          columns = c(
            "protein_name",
            "length",
            "sequence",
            "gene_primary",
            "gene_names",
            "go_f",
            "go_p",
            "go_c"
          ),
          batchsize = 100
        ) %>%
        dplyr::rename(
          protein_sequence = sequence,
          length_protein = length,
          pg_protein_accessions = accession
        )
      
      uniprot %>% write.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = ""))
      
    }
    
    protein_abundance_uniprot <- protein_abundance_imputed %>%
      left_join(uniprot, by = "pg_protein_accessions") %>% 
      filter(!pg_protein_accessions %in% (crapome_vector_final %>% pull(entry)))
    
    
    diff_abundance_treatment_vs_DMSO <- protein_abundance_uniprot %>%
      filter(grepl(bait_gene, new_condition_id) ) %>% 
      filter(grepl("probe", new_condition_id) | grepl("DMSO", new_condition_id)) %>% 
      assign_missingness(new_sample_id,
                         new_condition_id,
                         pg_protein_accessions,
                         normalised_intensity_log2,
                         ref_condition = paste(bait_gene, "DMSO", sep = "_")) %>%
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
      left_join(protein_abundance_uniprot, by = "pg_protein_accessions") %>%
      dplyr::mutate(significant = abs(diff) > 1 & adj_pval < 0.05) 
    
    
    # diff_abundance_treatment_vs_DMSO %>% write.csv(paste(getwd(), "/diff_abundance_probe_vs_DMSO_", bait_gene, ".csv", sep = ""))
    
    
    vol_plot <- make_volcano_plot(diff_abundance_treatment_vs_DMSO, label_vector = label_list[[bait_gene]])
    
    # pdf(paste("volcano_plot_probe_vs_DMSO_", bait_gene, ".pdf", sep = ""), 
    #     width = 5,
    #     height = 5,
    #     pointsize = 20)
    # print(vol_plot)
    # dev.off()
    
    diff_abundance_treatment_vs_DMSO_EGFP <- protein_abundance_uniprot %>%
      filter(grepl("EGFP", new_condition_id)) %>% 
      filter(grepl("probe", new_condition_id) | grepl("DMSO", new_condition_id)) %>% 
      assign_missingness(new_sample_id,
                         new_condition_id,
                         pg_protein_accessions,
                         normalised_intensity_log2,
                         ref_condition = paste("EGFP", "DMSO", sep = "_")) %>%
      calculate_diff_abundance(
        new_sample_id,
        new_condition_id,
        pg_protein_accessions,
        normalised_intensity_log2,
        missingness,
        comparison,
        ref_condition = paste("EGFP", "DMSO", sep = "_"),
        method = "moderated_t-test"
      ) %>%
      left_join(protein_abundance_uniprot, by = "pg_protein_accessions") %>%
      dplyr::mutate(significant = abs(diff) > 1 & adj_pval < 0.05) 
    
    # diff_abundance_treatment_vs_DMSO_EGFP %>% write.csv(paste(getwd(), "/diff_abundance_probe_vs_DMSO_", bait_gene, "_EGFP.csv", sep = ""))
    
    
    vol_plot <- make_volcano_plot(diff_abundance_treatment_vs_DMSO_EGFP, label_vector = label_list[[bait_gene]])
    
    # pdf(paste("volcano_plot_probe_vs_DMSO_", bait_gene, "_EGFP.pdf", sep = ""), 
    #     width = 5,
    #     height = 5,
    #     pointsize = 20)
    # print(vol_plot)
    # dev.off()
    

    kinase_upregulated_proteins <- diff_abundance_treatment_vs_DMSO %>% 
      filter(significant == TRUE) %>% 
      filter(diff > 1)
        
    EGFP_upregulated_proteins <- diff_abundance_treatment_vs_DMSO_EGFP %>% 
      filter(significant == TRUE) %>% 
      filter(diff > 1)
    
    
    kinase_downregulated_proteins <- diff_abundance_treatment_vs_DMSO %>% 
      filter(significant == TRUE) %>% 
      filter(diff < -1)
    
    EGFP_downregulated_proteins <- diff_abundance_treatment_vs_DMSO_EGFP %>% 
      filter(significant == TRUE) %>% 
      filter(diff < -1)
    
    proteins_to_exclude <- kinase_upregulated_proteins %>% 
      filter(pg_protein_accessions %in% (EGFP_upregulated_proteins %>% pull(pg_protein_accessions))) %>% 
      rbind(kinase_downregulated_proteins %>% 
              filter(pg_protein_accessions %in% (EGFP_downregulated_proteins %>% pull(pg_protein_accessions)))) %>% 
      pull(pg_protein_accessions) %>% 
      unique()
    
    diff_abundance_treatment_vs_DMSO_EGFP_filtered <- protein_abundance_uniprot %>%      
      filter(!pg_protein_accessions %in% c(proteins_to_exclude)) %>% 
      filter(grepl(bait_gene, new_condition_id) ) %>% 
      filter(grepl("probe", new_condition_id) | grepl("DMSO", new_condition_id)) %>% 
      assign_missingness(new_sample_id,
                         new_condition_id,
                         pg_protein_accessions,
                         normalised_intensity_log2,
                         ref_condition = paste(bait_gene, "DMSO", sep = "_")) %>%
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
      left_join(protein_abundance_uniprot, by = "pg_protein_accessions") %>%
      dplyr::mutate(significant = abs(diff) > 1 & adj_pval < 0.05)
      
    
    # diff_abundance_treatment_vs_DMSO_EGFP_filtered %>% write.csv(paste(getwd(), "/diff_abundance_probe_vs_DMSO_", bait_gene, "_EGFP_corr.csv", sep = ""))

    
    #volplot
    
    vol_plot <- make_volcano_plot(diff_abundance_treatment_vs_DMSO_EGFP_filtered, label_vector = label_list[[bait_gene]])
    
    pdf(paste("volcano_plot_probe_vs_DMSO_", bait_gene, "_EGFP_corr.pdf", sep = ""),
        width = 4,
        height = 4,
        pointsize = 20)
    print(vol_plot)
    dev.off()
    
    
    assign(sprintf("diff_abundance_%s", bait_gene), diff_abundance_treatment_vs_DMSO_EGFP_filtered)
  
}

all_PRKCA_changes <- diff_abundance_PRKCA_WT %>%
  rbind(diff_abundance_PRKCA_D463A) %>%
  rbind(diff_abundance_PRKCA_R22A_A25E) %>% 
  mutate(current_bait_id = str_extract(comparison, "PRKCA_[A-Z0-9]*")) %>% 
  mutate(current_bait_id = ifelse(current_bait_id == "PRKCA_R22A", "PRKCA_R22A_A25E", current_bait_id))

all_changes <- plot_changes_in_heatmap(data = all_PRKCA_changes, 
                                conditions = c("probe"),
                                significance_column = adj_pval,
                                levels_vector = c("PRKCA WT", "PRKCA D463A", "PRKCA R22A A25E"),
                                wt_bait_id = "PRKCA_WT",
                                fold_change_cutoff = 1,
                                include_annotations = TRUE)

# pdf(paste(getwd(), "/interaction_plot_PRKCA_mutants_only_1_mut.pdf", sep = ""), 
#     width = 3.5,
#     height = 5,
#     pointsize = 20)
# print(all_changes)
# dev.off()

