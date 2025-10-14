library(tidyverse)
library(protti)
library(readxl)
library(MBQN)
library(ggrepel)
library(doParallel)
library(msImpute)


####### Functions ####### 



filter_data_and_log2_transform <- function(data, 
                                           filter_for_proteotypicity = TRUE,
                                           proteins_exempt_from_proteotypicity_filter = bait_accession,
                                           protein_group_column = pg_protein_accessions,
                                           cutoff_eg_qvalue = 0.0001,
                                           qvalue_column = eg_qvalue,
                                           cutoff_fg_quantity_raw = 512,
                                           intensity_column = fg_quantity,
                                           cutoff_abs_diff_rt = 0.5,
                                           cutoff_min_number_of_obs_per_peptide = 3,
                                           apex_rt_column = eg_apex_rt,
                                           predicted_rt_column = eg_rt_predicted,
                                           digest_type_column = pep_type) {
  
  
  if(filter_for_proteotypicity == TRUE){
    data <- data %>% 
      filter(pep_is_proteotypic | ({{ protein_group_column }} %in% c(proteins_exempt_from_proteotypicity_filter)))
  }
  
  data <- data %>% 
    dplyr::mutate(new_condition_id = paste(bait_label, r_condition, sep = "_")) %>% 
    dplyr::mutate(new_sample_id = paste(new_condition_id, r_replicate, sep = "_")) %>% 
    filter({{ qvalue_column }} < 1e-3) %>%
    filter({{ intensity_column }} > 512) %>% 
    dplyr::mutate(eg_diff_rt = {{ apex_rt_column }} - {{ predicted_rt_column }}) %>% 
    filter(abs(eg_diff_rt) < 1) %>% 
    filter(!is.na(pg_organism_id)) %>% 
    filter({{ digest_type_column }} == "fully-tryptic") %>% 
    group_by(eg_precursor_id) %>% 
    mutate(n_obs = n()) %>% 
    ungroup() %>% 
    filter(n_obs > cutoff_min_number_of_obs_per_peptide) %>% 
    dplyr::mutate(intensity_log2 = log2(fg_quantity))
  
  return(data)
  
}


normalise_data <- function(data,
                           bait_accession = bait_accession,
                           precursor_column = eg_precursor_id,
                           measured_intensity = intensity_log2,
                           imputation_method = "none"
){
  
  
  DIA_matrix <- data  %>% 
    distinct({{ precursor_column }}, new_sample_id, {{ measured_intensity }}) %>% 
    pivot_wider(names_from = "new_sample_id", values_from = {{ measured_intensity }}) %>% 
    as.matrix(labels = TRUE)
  
  
  rownames(DIA_matrix) <- DIA_matrix[,rlang::as_name(rlang::enquo(precursor_column))]
  DIA_matrix <- DIA_matrix[, colnames(DIA_matrix) != rlang::as_name(rlang::enquo(precursor_column))]
  storage.mode(DIA_matrix) <- "numeric"
  
  normalised_matrix <- mbqnNRI(DIA_matrix, median, na.rm = TRUE)
  
  groups <- data %>% distinct(new_condition_id, new_sample_id) %>% pull(new_condition_id)
  
  if(imputation_method == "msimpute"){
    imputed_matrix <- msImpute(normalised_matrix, method = "v2-mnar", group = groups)
  } else{
    imputed_matrix <- normalised_matrix
  }
  
  
  DIA_matrix_normalised <- data.frame(imputed_matrix) %>% 
    rownames_to_column(var = rlang::as_name(rlang::enquo(precursor_column))) %>% 
    pivot_longer(cols = c(2:(normalised_matrix %>% colnames() %>% length() + 1)), names_to = "new_sample_id", values_to = "normalised_intensity_log2") %>% 
    filter(!is.na(normalised_intensity_log2)) %>% 
    left_join(data %>% distinct(eg_precursor_id, pg_protein_accessions, pep_stripped_sequence), 
              by = c(rlang::as_name(rlang::enquo(precursor_column)))) %>% 
    left_join(data %>% distinct(new_sample_id, eg_precursor_id, intensity_log2), 
              by = c(rlang::as_name(rlang::enquo(precursor_column)), "new_sample_id")) %>% 
    left_join(data %>% distinct(new_sample_id, new_condition_id, bait_label), 
              by = c(rlang::as_name("new_sample_id")))
  
  
  return(DIA_matrix_normalised)
  
}

bait_normalise <- function(data, bait_accession){
  
  DIA_matrix_normalised <- data %>%
    dplyr::filter(pg_protein_accessions %in% c(paste(bait_accession), "P42212")) %>% 
    group_by(eg_precursor_id) %>% 
    filter(n() >= (data$new_sample_id %>% unique() %>% length() / 2)) %>% 
    dplyr::mutate(med_bait_intensity = median(normalised_intensity_log2)) %>% 
    ungroup() %>%     
    dplyr::mutate(diff_bait_intensity = normalised_intensity_log2 - med_bait_intensity) %>% 
    group_by(new_sample_id) %>% 
    dplyr::mutate(med_diff_bait_intensity = mean(diff_bait_intensity)) %>% 
    ungroup() %>% 
    dplyr::distinct(new_sample_id, med_diff_bait_intensity) %>% 
    left_join(data, by = "new_sample_id", relationship = "many-to-many") %>% 
    dplyr::mutate(normalised_intensity_log2 = normalised_intensity_log2 - med_diff_bait_intensity) %>% 
    dplyr::mutate(uniprot_id = str_extract(pg_protein_accessions, "\\w+\\-?\\w?"))
  return(DIA_matrix_normalised)
  
}

find_cutoff_pval <- function(diff_abundance, cutoff = 0.05){
  diff_abundance <- diff_abundance %>% 
    filter(!is.na(adj_pval))
  
  if((pull(diff_abundance, adj_pval) %>% min()) > cutoff){ return(0) }
  
  largest_significant_pval <- diff_abundance %>% 
    filter(adj_pval < cutoff) %>% 
    pull(pval) %>% 
    max()
  
  smallest_insignificant_pval <- diff_abundance %>% 
    filter(adj_pval > cutoff) %>% 
    pull(pval) %>% 
    min()
  
  pval_matrix <- diff_abundance %>% 
    filter(pval %in% c(largest_significant_pval, smallest_insignificant_pval)) %>% 
    distinct(pval, adj_pval)
  
  regression_of_pvals <- lm(adj_pval ~ pval, pval_matrix)
  
  cutoff_pval <- (cutoff - regression_of_pvals$coefficients["(Intercept)"])/regression_of_pvals$coefficients["pval"]
  
  return(cutoff_pval)
}

make_volcano_plot <- function(data, log2fc_cutoff = 0.585, pval_cutoff = 0.05, peptide_labels){
  
  target_up <- data %>% 
    dplyr::filter(diff > 0) %>% 
    dplyr::filter(significant == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions, eg_precursor_id) %>% 
    left_join(uniprot %>% dplyr::select(c(pg_protein_accessions, gene_primary)), by = c("pg_protein_accessions" = "pg_protein_accessions"))
  
  target_down <- data %>% 
    dplyr::filter(diff < 0) %>% 
    dplyr::filter(significant == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions, eg_precursor_id) %>% 
    left_join(uniprot %>% dplyr::select(c(pg_protein_accessions, gene_primary)), by = c("pg_protein_accessions" = "pg_protein_accessions"))
  
  all_other_proteins <- data %>% 
    dplyr::filter(!pg_protein_accessions %in% (target_up %>% dplyr::pull(pg_protein_accessions))) %>% 
    dplyr::filter(!pg_protein_accessions %in% (target_down %>% dplyr::pull(pg_protein_accessions))) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions)
  
  if(!is.null(peptide_labels)){
      target_label <- data %>% 
    filter(significant == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions, eg_precursor_id) %>% 
    group_by(eg_precursor_id) %>% 
    mutate(label = peptide_labels[which(names(peptide_labels) == eg_precursor_id)[1]]) %>% 
    ungroup()
  } else{
    
    target_label <- data.frame(diff = numeric(0), pval = numeric(0), eg_precursor_id = character(0), label = character(0))
    
  }



  cutoff_p_value <- find_cutoff_pval(data, cutoff = pval_cutoff)
  
  
  volcano_plot_treatment_vs_DMSO <- data %>% 
    distinct(diff, pval, pg_protein_accessions, eg_precursor_id) %>% 
    ggplot(aes(x = diff, y = -log10(pval))) +
    geom_point(color = "white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = -log10(cutoff_p_value), linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "grey80") +
    # geom_label_repel(data = target_label,
    #                  aes(label = gene_primary, alpha = 1),
    #                  force = 2,
    #                  min.segment.length = unit(0, 'lines'),
    #                  #nudge_y = -0.5,
    #                  size = 2.5) +
    geom_point(data = all_other_proteins, 
               shape = 21,
               size = 2, 
               fill = "grey70") +
    geom_point(data = target_up,
               shape = 21,
               size = 2, 
               fill = "goldenrod", 
               colour = "black") + 
    geom_point(data = target_down,
               shape = 21,
               size = 2, 
               fill = "lightblue", 
               colour = "black") +
    geom_label_repel(data = target_label,
                     aes(label = label, alpha = 1),
                     max.overlaps = 4,
                     force = 2,
                     force_pull = 0.8,
                     direction = "both", 
                     nudge_x = -0.1,
                     nudge_y = -0.1,
                     min.segment.length = unit(0, 'lines')) +
    labs(x = expression(paste("log"[2]*"(fold change)")), y = expression(paste("-log"[10]*"(p-value)"))) +
    annotate("text", x = -max(abs(data %>% filter(!is.na(diff)) %>% pull(diff))) * 0.75, 
             y = (-log10(cutoff_p_value) - 0.15), label = paste("adjusted p = ", pval_cutoff, sep = ""), size = 3, alpha = 0.5) +
    xlim(c(-max(abs(data %>% filter(!is.na(diff)) %>% pull(diff))) - 0.2,
           max(abs(data %>% filter(!is.na(diff)) %>% pull(diff)))) + 0.2) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  
  
  return(volcano_plot_treatment_vs_DMSO)
  
}


