library(tidyverse)
library(protti)
library(readxl)
library(MBQN)
library(ggrepel)
library(doParallel)
library(msImpute)
library(pheatmap)
library(dichromat)
library(pROC)


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

make_volcano_plot <- function(data, log2fc_cutoff = 1, pval_cutoff = 0.05){
  
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
  
  known_interactors <- data %>% 
    filter(is_interactor == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions)
  
  all_other_proteins <- data %>% 
    dplyr::filter(!pg_protein_accessions %in% (target_up %>% dplyr::pull(pg_protein_accessions))) %>% 
    dplyr::filter(!pg_protein_accessions %in% (target_down %>% dplyr::pull(pg_protein_accessions))) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions)
  
  target_label <- data %>% 
    filter(significant == TRUE) %>% 
    dplyr::distinct(diff, pval, pg_protein_accessions) %>% 
    left_join(uniprot %>% dplyr::select(c(pg_protein_accessions, gene_primary)), by = c("pg_protein_accessions" = "pg_protein_accessions")) %>% 
    filter(pg_protein_accessions %in% c("Q13131",
                                        "P08237",
                                        "A0FGR8",
                                        "O75147",
                                        "P55786",
                                        "P27348",
                                        "P25205",
                                        "P33992",
                                        "Q9H078",
                                        "P12004",
                                        "Q93009",
                                        "P11802"
                                        
                                        ))
  
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
    geom_point(data = known_interactors,
               shape = 21,
               size = 2, 
               stroke = 1,
               fill = "grey70") +
    geom_point(data = target_up,
               shape = 21,
               size = 2, 
               fill = "mediumpurple2", 
               colour = "black") + 
    geom_point(data = target_down,
               shape = 21,
               size = 2, 
               fill = "yellowgreen", 
               colour = "black") +
    geom_point(data = known_interactors,
               shape = 21,
               size = 2, 
               stroke = 1,
               colour = "blue4") +
    geom_label_repel(data = target_label %>% filter(diff < 0),
                     aes(label = gene_primary, alpha = 1),
                     max.overlaps = 15,
                     force = 1,
                     force_pull = 2,
                     direction = "both", 
                     nudge_x = -0.1,
                     #nudge_y = -0.1,
                     min.segment.length = unit(0, 'lines'),
                     box.padding = 0.15) +
    geom_label_repel(data = target_label %>% filter(diff > 0),
                     aes(label = gene_primary, alpha = 1),
                     max.overlaps = 15,
                     force = 2,
                     force_pull = 2,
                     direction = "both", 
                     nudge_x = 0.2,
                     #nudge_y = -0.1,
                     min.segment.length = unit(0, 'lines'),
                     box.padding = 0.15) +
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




plot_individual_protein_change_CAMKK2 <- function(diff_abundance, precursor, colours, current_bait_id = current_bait_id,
                                                  grouping, diff = diff, pval = adj_pval, gene_id, levels_vector, label_vector, color, sem_data){
  
  data <- diff_abundance %>%
    filter(grepl("probe", comparison)) %>% 
    dplyr::distinct({{ grouping }}, {{ diff }}, {{ pval }}, gene_primary, current_bait_id) %>% 
    complete(current_bait_id, {{ grouping }}) %>% 
    filter({{ grouping }} == precursor) %>%
    dplyr::mutate(gene_primary = gene_id) %>% 
    dplyr::mutate(significant = ifelse(adj_pval < 0.05, TRUE, FALSE)) %>% 
    dplyr::mutate(current_bait_id = str_remove(current_bait_id, "CAMKK2_")) %>% 
    left_join((sem_data %>% 
                 dplyr::mutate(current_bait_id = str_remove(current_bait_id, "CAMKK2_")) %>% 
                 rbind(data.frame(current_bait_id = "R311C", total_sem = 0))), by = "current_bait_id") %>%  
    dplyr::mutate(current_bait_id = factor(current_bait_id, levels = levels_vector)) %>% 
    arrange(current_bait_id)
  
  data$label <- label_vector
  
  plot <- data %>%
    dplyr::mutate(diff = ifelse( is.na({{ diff }}), 0, {{ diff }})) %>% 
    mutate(label_pos = {{ diff }} + sign({{ diff }}) * (total_sem)) %>% 
    ggplot(aes(x = current_bait_id, y = {{diff}}, fill = current_bait_id)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_text(aes(y = label_pos, label = label), vjust = -0.5, size = 4, parse=TRUE) +
    geom_errorbar(aes(ymin={{diff}} - total_sem, ymax= {{diff}} + total_sem), width=.2,
                  position=position_dodge(.9)) +
    #geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "solid") +
    scale_y_continuous(limits = c(-1, 4), expand = c(0,0.1)) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = "black") +
    labs(title = paste(unique(dplyr::pull(data, gene_primary))),
         y = expression(log[2] (fold~change))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 12))
  
  return(plot)
}

plot_individual_protein_change_CHEK1 <- function(diff_abundance, precursor, colours, current_bait_id = current_bait_id,
                                                 grouping, diff = diff, pval = adj_pval, gene_id, levels_vector, label_vector, color, sem_data){
  
  data <- diff_abundance %>%
    filter(grepl("probe", comparison)) %>% 
    dplyr::distinct({{ grouping }}, {{ diff }}, {{ pval }}, gene_primary, current_bait_id) %>% 
    complete(current_bait_id, {{ grouping }}) %>% 
    filter({{ grouping }} == precursor) %>%
    dplyr::mutate(gene_primary = gene_id) %>% 
    dplyr::mutate(significant = ifelse(adj_pval < 0.05, TRUE, FALSE)) %>% 
    dplyr::mutate(current_bait_id = str_remove(current_bait_id, "CHEK1_")) %>% 
    left_join((sem_data %>% 
                 dplyr::mutate(current_bait_id = str_remove(current_bait_id, "CHEK1_"))), by = "current_bait_id") %>%  
    dplyr::mutate(current_bait_id = factor(current_bait_id, levels = levels_vector)) %>% 
    arrange(current_bait_id)
  
  data$label <- label_vector
  
  plot <- data %>%
    dplyr::mutate(diff = ifelse( is.na({{ diff }}), 0, {{ diff }})) %>% 
    mutate(label_pos = {{ diff }} + sign({{ diff }}) * (total_sem + 0.1)) %>% 
    ggplot(aes(x = current_bait_id, y = {{diff}}, fill = current_bait_id)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_text(aes(y = label_pos, label = label), size = 4, parse=TRUE) +
    geom_errorbar(aes(ymin={{diff}} - total_sem, ymax= {{diff}} + total_sem), width=.2,
                  position=position_dodge(.9)) +
    #geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "solid") +
    scale_y_continuous(limits = c(-2, 1), expand = c(0,0.1)) +
    scale_fill_manual(values = color) +
    scale_color_manual(values = "black") +
    labs(title = paste(unique(dplyr::pull(data, gene_primary))),
         y = expression(log[2] (fold~change))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 12))
  
  return(plot)
}


