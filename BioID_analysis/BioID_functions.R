library(tidyverse)
library(protti)
library(readxl)
library(MBQN)
library(ggrepel)
library(doParallel)
library(msImpute)
library(pheatmap)





filter_data_and_log2_transform <- function(data, 
                                           filter_for_proteotypicity = TRUE,
                                           proteins_exempt_from_proteotypicity_filter = bait_accession,
                                           protein_group_column = pg_protein_accessions,
                                           cutoff_eg_qvalue = 0.0001,
                                           qvalue_column = eg_qvalue,
                                           cutoff_fg_quantity_raw = 512,
                                           intensity_column = fg_quantity,
                                           cutoff_abs_diff_rt = 0.5,
                                           apex_rt_column = eg_apex_rt,
                                           predicted_rt_column = eg_rt_predicted) {
  
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
    group_by(eg_precursor_id) %>% 
    mutate(n_obs = n()) %>% 
    ungroup() %>% 
    filter(n_obs > 3) %>% 
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

make_volcano_plot <- function(data, label_vector, log2fc_cutoff = 1, pval_cutoff = 0.05){
  
  target_up <- data %>% 
    dplyr::filter(diff > 0) %>% 
    dplyr::filter(significant == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions) %>% 
    left_join(uniprot %>% dplyr::select(c(pg_protein_accessions, gene_primary)), by = c("pg_protein_accessions" = "pg_protein_accessions"))
  
  target_down <- data %>% 
    dplyr::filter(diff < 0) %>% 
    dplyr::filter(significant == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions) %>% 
    left_join(uniprot %>% dplyr::select(c(pg_protein_accessions, gene_primary)), by = c("pg_protein_accessions" = "pg_protein_accessions"))
  
  all_other_proteins <- data %>% 
    dplyr::filter(!pg_protein_accessions %in% (target_up %>% dplyr::pull(pg_protein_accessions))) %>% 
    dplyr::filter(!pg_protein_accessions %in% (target_down %>% dplyr::pull(pg_protein_accessions))) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions)
  
  target_label <- data %>% 
    #filter(significant == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions) %>% 
    left_join(uniprot %>% dplyr::select(c(pg_protein_accessions, gene_primary)), by = c("pg_protein_accessions" = "pg_protein_accessions")) %>% 
    filter(pg_protein_accessions %in% label_vector)
  
  cutoff_p_value <- find_cutoff_pval(data, cutoff = pval_cutoff)
  
  
  volcano_plot_treatment_vs_DMSO <- data %>% 
    distinct(diff, pval, pg_protein_accessions) %>% 
    ggplot(aes(x = diff, y = -log10(pval))) +
    geom_point(color = "white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = -log10(cutoff_p_value), linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "grey80") +
    geom_point(data = all_other_proteins, 
               shape = 21,
               size = 2, 
               fill = "grey70") +
    geom_point(data = target_up,
               shape = 21,
               size = 2, 
               fill = "mediumpurple", 
               colour = "black") + 
    geom_point(data = target_down,
               shape = 21,
               size = 2, 
               fill = "yellowgreen", 
               colour = "black") +
    geom_label_repel(data = target_label,
                     aes(label = gene_primary, alpha = 1),
                     max.overlaps = 20,
                     force = 1.5,
                     force_pull = 2,
                     direction = "y", 
                     #nudge_x = -0.1,
                     #nudge_y = 0.1, #remove for AMPK
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

plot_changes_in_heatmap <- function(data, 
                                    conditions,
                                    bait_accession = bait_accession,
                                    protein_column = gene_primary,
                                    diff_column = diff,
                                    significance_column = adj_pval,
                                    levels_vector = NULL,
                                    wt_bait_id = "CAMKK2_WT",
                                    fold_change_cutoff = 1,
                                    include_annotations = TRUE) {
  
  plot_data <- data.frame()
  
  for(condition_index in conditions){
    
    plot_data <- plot_data %>% 
      rbind(data %>% 
              filter(!grepl("EGFP", comparison)) %>% 
              filter(grepl(condition_index, comparison)) %>%
              mutate(treatment_condition = condition_index) %>% 
              distinct({{ protein_column }}, pg_protein_accessions, {{ diff_column }}, {{ significance_column }}, current_bait_id, treatment_condition) )
  }
  
  max_abs_value <- plot_data %>% 
    pull({{diff_column}}) %>% 
    abs() %>% 
    max()
  
  significant_proteins_WT <- plot_data %>% 
    mutate(is_significant = ifelse((( abs({{ diff_column }}) > fold_change_cutoff) & ( {{ significance_column }} < 0.05)), TRUE, FALSE)) %>% 
    filter(current_bait_id == wt_bait_id) %>% 
    filter(is_significant == TRUE) %>% 
    pull({{ protein_column }})
  
  annotation_data <- plot_data %>% 
    mutate(is_significant = ifelse((( abs({{ diff_column }}) > fold_change_cutoff) & ( {{ significance_column }} < 0.05)), TRUE, FALSE)) %>% 
    group_by({{ protein_column }}) %>%
    mutate(protein_significance = ifelse({{ protein_column }} %in% significant_proteins_WT, 1, 0)) %>% 
    #mutate(protein_significance = ifelse(current_bait_id == wt_bait_id, is_significant[which(current_bait_id == wt_bait_id)], FALSE)) %>%                    #max(is_significant)) %>%
    ungroup() %>%
    #mutate(protein_significance = ifelse( ({{ protein_column }} == "Q96RR4"), 1, protein_significance)) %>%
    filter(protein_significance == 1) %>%
    mutate(current_bait_id = str_replace_all(current_bait_id, "_", " ")) %>%
    mutate(label = paste(current_bait_id)) %>% #, treatment_condition
    mutate(sign_col = ifelse(({{ significance_column }} < 0.01) & include_annotations == TRUE, "**", "*")) %>%
    mutate(sign_col = ifelse(({{ significance_column }} > 0.05) | include_annotations == FALSE, " ", sign_col)) %>%
    #mutate(sign_col = ifelse(include_annotations == FALSE, " ", sign_col)) %>%
    dplyr::select(-c({{ diff_column }}, is_significant, protein_significance, {{ significance_column }}, pg_protein_accessions)) %>% #, current_bait_id, treatment_condition)) %>%
    pivot_wider(names_from = {{ protein_column }}, values_from = sign_col, values_fill = "ND")
    
  
  if(!is.null(levels_vector)){
    annotation_data <- annotation_data   %>%
      mutate(label = factor(label, levels = levels_vector)) %>%
      arrange(label)
  }
  
  annotation_data <- as.matrix(annotation_data, labels = TRUE)
  rownames(annotation_data) <- annotation_data[,"label"]
  annotation_data <- annotation_data[, -c(1:3)]
  #storage.mode(annotation_data) <- "numeric"
  

  annotation_row_mat <- plot_data %>% 
    mutate(Interactor = ifelse( pg_protein_accessions %in% biogrid_interactors_vector, "Known", "Novel")) %>% 
    mutate(label = {{ protein_column }}) %>% 
    distinct(label, Interactor) %>% 
    filter(!is.na(label)) %>% 
    column_to_rownames(var="label")
    
  #return(annotation_row_mat)

  
  plot_data <- plot_data %>% 
    mutate(is_significant = ifelse((( abs({{ diff_column }}) > 1) & ( {{ significance_column }} < 0.05)), TRUE, FALSE)) %>% 
    group_by({{ protein_column }}) %>% 
    mutate(protein_significance = ifelse({{ protein_column }} %in% significant_proteins_WT, 1, 0)) %>% 
    ungroup() %>% 
    #mutate(protein_significance = ifelse( ({{ protein_column }} == "Q96RR4"), 1, protein_significance)) %>% 
    filter(protein_significance == 1)  %>% 
    mutate(current_bait_id = str_replace_all(current_bait_id, "_", " ")) %>% 
    mutate(label = paste(current_bait_id)) %>% #, treatment_condition
    dplyr::select(-c({{ significance_column }}, is_significant, protein_significance, pg_protein_accessions)) %>% #, current_bait_id, treatment_condition)) %>% 
    pivot_wider(names_from = {{ protein_column }}, values_from = {{ diff_column }}, values_fill = 0) 
  
  
  if(!is.null(levels_vector)){
    plot_data <- plot_data   %>%
      mutate(label = factor(label, levels = levels_vector)) %>%
      arrange(label)
  }
  

  plot_data <- as.matrix(plot_data, labels = TRUE)
  rownames(plot_data) <- plot_data[,"label"]
  plot_data <- plot_data[, -c(1:3)]
  storage.mode(plot_data) <- "numeric"
  plot_data <- t(plot_data)
  
  ann_colors = list(
    Interactor = c(Known="#3DC1FF", Novel="#ebeb26"))
  
  my_breaks <- seq(from = (-max_abs_value - 0.25), to = (max_abs_value + 0.25), length.out = 100)
  
  plot <-
    pheatmap(
      plot_data,
      color = colorRampPalette(
        c(
          "#8bbd26",
          "#9cc532",
          "#a6d448",
          "#b7e063",
          "#c6e880",
          "#d7f2a0",
          "#f7f7f7",
          "#e0d2fc",
          "#c5aef2",
          "#b69aed",
          "#a07fe3",
          "#8f75CB",
          "#7E69AF"
        )
      )(100), 
      breaks = my_breaks,
      #colorspace::diverge_hsv(n = 11),
      #annotation_row = annotation_col_mat,
      annotation_row = annotation_row_mat,
      annotation_colors = ann_colors,
      cluster_cols = FALSE,
      #cluster_rows = FALSE,
      #clustering_distance_rows = "manhattan",
      clustering_distance_rows = "euclidean",
      display_numbers = t(annotation_data),
      fontsize_number = 5,
      fontsize = 5,
      cutree_rows = 3,
      number_color = "#3b3100",
      border_color = "white",
      angle_col = 45
    ) 
  
  
  
  # clustered_data <- plot_data %>% 
  #   do(clusters = hclust(dist(.[colnames(plot_data %>% dplyr::select( {{ diff_column }}))])))
  
  
  return(plot)
  
}
