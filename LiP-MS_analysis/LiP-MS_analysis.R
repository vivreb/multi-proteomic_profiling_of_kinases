
# Load functions and packages 

rm(list = ls())

source("LiP-MS_functions.R")

###################################################################################
######################### Info to change for new analysis #########################
###################################################################################

# Load data

DIA_raw_CAMKK2_WT <- read_protti("VR_Ex76_CAMKK2_wt_Report_VR_QC (Normal).tsv") %>% 
                          mutate(bait_label = "CAMKK2_WT")

DIA_raw_CAMKK2_R311C <- read_protti("VR_Ex76_CAMKK2_R311C_Report_VR_QC (Normal).tsv") %>% 
                                mutate(bait_label = "CAMKK2_R311C")

DIA_raw_CAMKK2_D312A <- read_protti("VR_Ex76_CAMKK2_D312A_Report_VR_QC (Normal).tsv") %>%
                                mutate(bait_label = "CAMKK2_D312A")

DIA_raw_CAMKK2_T483D <- read_protti("VR_Ex76_CAMKK2_T483D_Report_VR_QC (Normal).tsv") %>%
                                mutate(bait_label = "CAMKK2_T483D")

DIA_raw_CHEK1_WT <- read_protti("VR_Ex56_CHEK1_Report_VR_QC (Normal).tsv") %>%
                          mutate(bait_label = "CHEK1_WT")

DIA_raw_PRKCA_WT <- read_protti("VR_Ex63_PRKCA_Report_VR_QC (Normal).tsv") %>% 
                            mutate(bait_label = "PRKCA_WT")

DIA_raw_DCLK1 <- read_protti("VR_Ex46_DCLK1_SN19_Report_VR_QC (Normal).tsv") %>%
                                mutate(bait_label = "DCLK1")


# Initialize list to iterate over

bait_gene_list <- c("CAMKK2_WT", "CAMKK2_R311C", "CAMKK2_D312A", "CAMKK2_T483D", 
                    "CHEK1_WT", 
                    "PRKCA_WT",
                    "DCLK1")

bait_accession_list <- c("Q96RR4", "Q96RR4", "Q96RR4", "Q96RR4",
                         "O14757",
                         "P17252", 
                         "O15075")
names(bait_accession_list) <- c(bait_gene_list)

crapome_vector_final <- read.csv("final_crapome_proteins.csv")

###################################################################################
###################################################################################
###################################################################################



for(bait_gene in bait_gene_list) {

    set.seed(123)
    
    DIA_raw <- sprintf("DIA_raw_%s", bait_gene) %>% 
      get() 
    
    bait_accession <- bait_accession_list[bait_gene]
    
    DIA_raw <- DIA_raw %>%
      filter(!pg_protein_accessions %in% crapome_vector_final$entry)
    
    
    DIA_clean <- DIA_raw  %>%       
      filter_data_and_log2_transform(
        cutoff_eg_qvalue = 0.001,
        cutoff_abs_diff_rt = 0.1,
        cutoff_fg_quantity_raw = 1024)
    
    DIA_clean_normalised <- DIA_clean %>% 
      normalise_data(imputation_method = "none") %>% 
      dplyr::select(new_sample_id, eg_precursor_id, normalised_intensity_log2) %>% 
      left_join((DIA_clean %>%
                   mutate(new_sample_id = paste(bait_label, r_condition, r_replicate, sep = "_"))),
                by = c("new_sample_id", "eg_precursor_id"))
    
    
    unis <- DIA_clean %>%
      pull(pg_protein_accessions) %>%
      strsplit(";") %>%
      unlist() %>%
      unique()
    
    
    uniprot <- read.csv(paste(getwd(), "/uniprot_", bait_gene, ".csv", sep = "")) %>% 
        dplyr::select(-c(X))
    
    low_cv_peptides <- DIA_clean_normalised %>%
      filter(r_condition %in% c("DMSO", "probe")) %>%
      filter(!is.na(normalised_intensity_log2)) %>%
      mutate(normalised_intensity = 2 ^ normalised_intensity_log2) %>%
      mutate(peptide = paste(eg_precursor_id, r_condition, sep = "+")) %>%
      qc_cvs(grouping = eg_precursor_id,
             condition = peptide,
             intensity = normalised_intensity,
             plot = FALSE) %>%
      mutate(peptide = str_replace(peptide, "[\\+].+", "")) %>%
      group_by(peptide) %>%
      mutate(max_cv = max(median_cv)) %>%
      ungroup() %>%
      filter(max_cv < 1000) %>%
      pull(peptide)
    
    DIA_clean_uniprot <- DIA_clean_normalised %>%
      filter(r_condition %in% c("DMSO", "probe")) %>% 
      filter(eg_precursor_id %in% low_cv_peptides)%>%
      left_join(uniprot, by = c("pg_protein_accessions" = "uniprot_id")) %>%
      find_peptide(protein_sequence, pep_stripped_sequence) %>%
      assign_peptide_type(aa_before, last_aa)  %>% 
      filter(n() >= 3) %>%
      ungroup() %>%
      group_by(eg_precursor_id) %>%
      filter(n() >= 6) %>%
      ungroup()
    
    
    diff_abundance <- DIA_clean_uniprot %>%
      mutate(r_condition = str_replace(r_condition, "/", "_")) %>% 
      assign_missingness(r_file_name,
                         r_condition,
                         eg_precursor_id,
                         normalised_intensity_log2,
                         ref_condition = "DMSO") %>%
      calculate_diff_abundance(
        r_file_name,
        r_condition,
        eg_precursor_id,
        normalised_intensity_log2,
        missingness,
        comparison,
        ref_condition = "DMSO",
        method = "moderated_t-test"
      ) %>%
      left_join(DIA_clean_uniprot, by = c("eg_precursor_id" = "eg_precursor_id")) %>%
      mutate(significant = abs(diff) > 0.585 & adj_pval < 0.05) 
    
    diff_abundance %>% write.csv(paste(getwd(), "/diff_abundance_all_", bait_gene, "_probe.csv", sep = ""))
    
    diff_abundance %>%
      filter(pg_protein_accessions == bait_accession) %>%
      filter(significant == TRUE) %>%
      write.csv(paste(getwd(), "/diff_abundance_bait_significant_", bait_gene, "_probe.csv", sep = ""))
    
    
    volcano_plot <- diff_abundance %>%
      mutate(hit_type = "no") %>%
      mutate(hit_type = ifelse((pg_protein_accessions == bait_accession) & (significant == TRUE), "on_target_significant", hit_type)) %>%
      mutate(hit_type = ifelse((pg_protein_accessions == bait_accession) & (significant == FALSE), "on_target_insignificant", hit_type)) %>%
      volcano_plot(
        eg_precursor_id,
        diff,
        pval,
        colour = c("grey60", "#5680C1", "#d9363c", "#B96DAD", "#64CACA", "#81ABE9", "#F6B8D1", "#99F1E4", "#9AD1FF", "#548BDF", "#A55098", "#3EB6B6",
                   "#87AEE8", "#CA91C1", "#A4E0E0", "#1D4F9A", "#D7ACD2", "#49C1C1"),
        method = "target",
        target = c("on_target_significant", "on_target_insignificant"),
        target_column = hit_type,
        x_axis_label = "log2(fold change)",
        title = "",
        log2FC_cutoff = 0.585,
        significance_cutoff = c(0.05, "adj_pval"),
        interactive = FALSE
      ) + theme(legend.position = "none",
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                panel.border = element_blank(),
                axis.line = element_line(colour = "black"))
    
    
    barcode_plot <- diff_abundance %>% 
      filter(pg_protein_accessions == bait_accession) %>% 
      calculate_sequence_coverage(protein_sequence, pep_stripped_sequence) %>% 
      barcode_plot_new(
        start_position = start,
        end_position = end,
        protein_length = length_protein,
        protein_id = gene_primary,
        coverage = coverage,
        cutoffs = c(diff = 0.585, adj_pval = 0.05)
      ) +
      labs(fill = "", title = "", x = "") +
      theme(legend.position = "none")
    
    pdf(paste("barcode_plot_protti_", bait_gene, "_probe.pdf", sep = ""), 
        width = 10,
        height = 2,
        pointsize = 30)
    print(barcode_plot)
    dev.off()

}


