rm(list = ls())

source("BioID_functions.R")


###################################################################################
######################### Info to change for new analysis #########################
###################################################################################


DIA_raw_PRKCA_WT <- read_protti("Z:/Viviane/Experiments/VR_Ex114/20260502_215103_VR_Ex114_PRKCA_dose_response_BioID/VR_Ex114_PRKCA_dose_response_BioID_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_WT") %>% 
  mutate(r_condition = str_replace_all(r_condition, " ", "_"))


# Define treatments

condition_vector <- c("PRKCA_WT_DMSO", "PRKCA_WT_100_pM", "PRKCA_WT_1_nM", "PRKCA_WT_10_nM", "PRKCA_WT_100_nM", "PRKCA_WT_1_uM", "PRKCA_WT_10_uM")
numeric_condition_vector <- c(0, 100, 1000, 10000, 100000, 1000000, 10000000)
names(numeric_condition_vector) <- condition_vector

treatment_list = list()
treatment_list[["PRKCA_WT"]] = condition_vector[-1]


# Initialize list to iterate over

bait_gene_list <- c("PRKCA_WT")

bait_accession_list <- c("P17252")
names(bait_accession_list) <- c(bait_gene_list)


biogrid_interactors_PRKCA <-
  read_delim(
    "Z:/Viviane/Experiments/VR_Ex59/BIOGRID-GENE-111564-4.4.233.tab3.txt",
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

bait_gene = "PRKCA_WT"

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
  
  
  fit <- protein_abundance_uniprot %>% 
    mutate(doses = numeric_condition_vector[new_condition_id]) %>% 
    fit_drc_4p(sample = new_sample_id, 
               grouping = pg_protein_accessions,
               response = normalised_intensity_log2,
               dose = doses,
               filter = "post", 
               retain_columns = c(protein_name, gene_primary)) 
  
  fit %>% 
    dplyr::select(-c(plot_curve, plot_points)) %>% 
    write.csv("fit_df_raw.csv")
  
}
  
  hits <- fit %>% 
    filter(passed_filter == TRUE) %>% 
    filter(abs(max_model - min_model) > 0.585) %>% 
    filter(abs(hill_coefficient) > 0.35) %>% 
    filter(correlation > 0.85) %>% 
    arrange((abs(ec_50))) %>% 
    pull(gene_primary)
  
  sign_drc_plot <-   fit %>% 
    drc_4p_plot(
      grouping = gene_primary,
      dose = doses,
      response = normalised_intensity_log2,
      targets = cell_junction_plot_gene_vector,
      unit = "pM",
      export = TRUE
    )
  
  cell_junction_plot_protein_vector <- cell_junction_plot_proteins %>% pull(pg_protein_accessions)
  cell_junction_plot_gene_vector <- cell_junction_plot_proteins %>% pull(gene_primary)
  
  significant_proteins_single_dose <- diff_abundance_probe_vs_DMSO_PRKCA_WT_EGFP_corr %>% filter(significant == TRUE) %>% pull(pg_protein_accessions) %>% unique()
  detected_proteins_single_dose <- diff_abundance_probe_vs_DMSO_PRKCA_WT_EGFP_corr %>% pull(pg_protein_accessions) %>% unique()
  
  # Turn curves into data frame for plot
  
  curve_df <- data.frame(pg_protein_accessions = character(0), dose = numeric(0), prediction = numeric(0))
  
  for(protein in fit$pg_protein_accessions){
    
    tmp_curve <- fit$plot_curve[which(fit$pg_protein_accessions == protein)]
    
    tmp_df <- data.frame(pg_protein_accessions = rep(protein, 100), dose = tmp_curve[[1]]$dose, prediction = tmp_curve[[1]]$Prediction)
    
    curve_df <- curve_df %>% rbind(tmp_df)
    
  }
  
  curve_df_annot <- curve_df %>% 
    filter(!is.na(dose)) %>% 
    mutate(is_target = ifelse(pg_protein_accessions %in% significant_proteins_single_dose, "Yes", "No")) %>% 
    group_by(pg_protein_accessions) %>% 
    mutate(min_dose = min(dose)) %>% 
    mutate(
      dose_100 = prediction[abs(dose - 100) < 1e-5][1], 
      normalised_prediction = abs(prediction - dose_100),
      scaled_prediction = (normalised_prediction - min(normalised_prediction)) / 
        (max(normalised_prediction) - min(normalised_prediction))) %>% 
  ungroup()
  
  filter_vector <- fit %>% 
    filter(pg_protein_accessions %in% detected_proteins_single_dose) %>% 
    filter(passed_filter == TRUE) %>% 
    filter(abs(max_model - min_model) > 0.585) %>% 
    filter(abs(hill_coefficient) > 0.35) %>% 
    filter((correlation > 0.8)) %>% 
    pull(pg_protein_accessions)
  
  for_plot <- curve_df_annot %>% 
    filter(pg_protein_accessions %in% filter_vector)

 # Test if groups are different
  
  wil_test_1 <- for_plot %>% 
    filter(is_target == "Yes") %>% 
    pull(ec_50) %>% 
    unique()
  
  median(wil_test_1)
  
  wil_test_2 <- for_plot %>% 
    filter(is_target == "No") %>% 
    pull(ec_50) %>% 
    unique()
  
  median(wil_test_2)
  
  wilcox.test(wil_test_1, wil_test_2, alternative = "two.sided")
  
  
  bars.df <- data.frame(x = c(4.5), y = c(0.9), xend = c(4.5), yend = c(2.1))
  
  label.df <- data.frame(ec_50 = c(9000), 
                         y = c(1.5),
                         label = paste0("p == 1.8", "%*%", "10^{-11}"))
  
  pdf(paste("box_plot_filtered_PRKCA.pdf", sep = ""),
      width = 5,
      height = 1.5,
      pointsize = 20)
  for_plot %>% 
    distinct(pg_protein_accessions, ec_50, is_target) %>% 
    ggplot(aes(y = is_target, x = log10(ec_50), fill = is_target)) +
    geom_boxplot(outliers = FALSE) +
    geom_text(data = label.df, aes(y = y, x = log10(ec_50), label = label), inherit.aes = FALSE, size = 3, parse = TRUE) +
    geom_segment(data = bars.df, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, size = 0.5) +
    scale_x_continuous(limits = c(2, 7.5), expand = c(0.01, 0.01)) +
    scale_fill_manual(values = c("Yes" = "mediumpurple", "No" = "#BBBBBB")) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  dev.off()

  pdf(paste("curve_plot_filtered_PRKCA.pdf", sep = ""),
      width = 5,
      height = 3,
      pointsize = 20)
  for_plot %>% 
    ggplot(aes(x = log10(dose), y = abs(scaled_prediction), group = pg_protein_accessions, color = is_target, size = is_target, alpha = is_target)) +
    geom_line() +
    scale_size_manual(values = c("Yes" = 1.5, "No" = 1.5, "lower_ec50" = 1)) +
    scale_alpha_manual(values = c("Yes" = 1, "No" = 0.2, "lower_ec50" = 0.5)) +
    scale_color_manual(values = c("Yes" = "mediumpurple", "No" = "#BBBBBB")) +
    #geom_vline(xintercept = 6, linetype = "dashed")+
    geom_hline(yintercept = 0.5, linetype = "dashed")+
    theme_bw() +
    geom_line(data = filter(for_plot, is_target == "Yes"), 
              color = "mediumpurple", size = 1.5) +
    labs(x = "log10(concentration (pM))", y = "Absolute relative abundance change") +
    scale_x_continuous(limits = c(2, 7.5), expand = c(0.01, 0.01)) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  dev.off()
