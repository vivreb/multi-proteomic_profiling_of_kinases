library(tidyverse)
library(protti)
library(readxl)
library(MBQN)
library(ggrepel)

barcode_plot_new <- function(data,
                             start_position,
                             end_position,
                             protein_length,
                             coverage = NULL,
                             colouring = NULL,
                             protein_id = NULL,
                             facet = NULL,
                             cutoffs = NULL) {
  # Check if there is more than one protein even though protein_id was specified.
  if (!missing(protein_id)) {
    if (length(unique(dplyr::pull(data, {{ protein_id }}))) > 1) {
      stop("If data contains information of multiple proteins use the facet argument, not the protein_id argument")
    }
  }
  # Check if there are more than 20 proteins for faceting.
  if (!missing(facet)) {
    if (length(unique(dplyr::pull(data, {{ facet }}))) > 20) {
      n_proteins <- length(unique(dplyr::pull(data, {{ facet }})))
      twenty_proteins <- unique(dplyr::pull(data, {{ facet }}))[1:20]
      data <- data %>%
        dplyr::filter({{ facet }} %in% twenty_proteins)
      warning(paste(
        "Only the first 20 proteins from", rlang::as_name(enquo(facet)),
        "have been used for plotting since there are", n_proteins,
        "proteins. Consider mapping over subsetted datasets."
      ))
    }
  }
  # Apply fold change  and significance cutoff if fold change is provided
  if (!missing(cutoffs)) {
    fc_name <- names(cutoffs)[1]
    sig_name <- names(cutoffs)[2]
    fc <- cutoffs[1]
    sig <- cutoffs[2]
    
    colouring <- sym("change")
    
    data <- data %>%
      dplyr::mutate({{ colouring }} := ifelse(((!!ensym(fc_name) >= fc | !!ensym(fc_name) <= -fc) & !!ensym(sig_name) <= sig), "Changed", "Unchanged")) %>%
      dplyr::mutate({{ colouring }} := forcats::fct_rev(ifelse(is.na({{ colouring }}), "Unchanged", {{ colouring }}))) %>%
      dplyr::arrange({{ colouring }})
  }
  # Add coverage to protein ID name if present.
  if (!missing(coverage) & !missing(facet)) {
    data <- data %>%
      mutate({{ facet }} := paste0({{ facet }}, " (", round({{ coverage }}, digits = 1), "%)"))
  }
  if (!missing(coverage) & !missing(protein_id)) {
    data <- data %>%
      mutate({{ protein_id }} := paste0({{ protein_id }}, " (", round({{ coverage }}, digits = 1), "%)"))
  }
  # Create plot
  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(
      ymin = -2.5,
      ymax = 2.5,
      xmax = {{ end_position }} / {{ protein_length }} * 100,
      xmin = ({{ start_position }} - 1) / {{ protein_length }} * 100,
      fill = {{ colouring }}
    ),
    size = 0.7
    ) +
    ggplot2::scale_fill_manual(values = c(
      "#5680C1", "#d9363c", "#B96DAD", "#64CACA", "#81ABE9", "#F6B8D1", "#99F1E4", "#9AD1FF", "#548BDF", "#A55098", "#3EB6B6",
      "#87AEE8", "#CA91C1", "#A4E0E0", "#1D4F9A", "#D7ACD2", "#49C1C1"
    )) +
    ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = "Protein Sequence", title = {
      if (!missing(protein_id)) unique(dplyr::pull(data, {{ protein_id }}))
    }) +
    {
      if (!missing(facet)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet)))
    } +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.ticks.x = element_blank(),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 15),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA)
    )
}



filter_data_and_log2_transform <- function(data, 
                                           filter_for_proteotypicity = TRUE,
                                           proteins_exempt_from_proteotypicity_filter = bait_accession,
                                           protein_group_column = pg_protein_accessions,
                                           cutoff_eg_qvalue = 0.001,
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
    filter({{ qvalue_column }} < cutoff_eg_qvalue) %>%
    filter({{ intensity_column }} > cutoff_fg_quantity_raw) %>% 
    dplyr::mutate(eg_diff_rt = {{ apex_rt_column }} - {{ predicted_rt_column }}) %>% 
    filter(abs(eg_diff_rt) < cutoff_abs_diff_rt) %>% 
    filter(!is.na(pg_organism_id)) %>% 
    dplyr::mutate(intensity_log2 = log2(fg_quantity))
  
  return(data)
  
}


normalise_data <- function(data,
                           precursor_column = eg_precursor_id,
                           measured_intensity = intensity_log2,
                           imputation_method
){
  
  if(imputation_method == "msimpute"){
    data <- data %>% 
      group_by({{precursor_column}}) %>% 
      mutate(n_obs = n()) %>% 
      ungroup() %>% 
      filter(n_obs > 3) 
  }
  
  DIA_matrix <- data  %>% 
    distinct({{ precursor_column }}, new_sample_id, {{ measured_intensity }}) %>% 
    pivot_wider(names_from = "new_sample_id", values_from = {{ measured_intensity }}) %>% 
    as.matrix(labels = TRUE)
  
  
  rownames(DIA_matrix) <- DIA_matrix[,rlang::as_name(rlang::enquo(precursor_column))]
  DIA_matrix <- DIA_matrix[, colnames(DIA_matrix) != rlang::as_name(rlang::enquo(precursor_column))]
  storage.mode(DIA_matrix) <- "numeric"
  
  normalised_matrix <- mbqnNRI(DIA_matrix, median, na.rm = TRUE)
  
  if(imputation_method == "msimpute"){
    groups <- data %>% distinct(new_condition_id, new_sample_id) %>% pull(new_condition_id)
  
    imputed_matrix <- msImpute(normalised_matrix, method = "v2-mnar", group = groups)
    
  } else {
    imputed_matrix = normalised_matrix
  }
  
  
  
  DIA_matrix_normalised <- data.frame(imputed_matrix) %>% 
    rownames_to_column(var = rlang::as_name(rlang::enquo(precursor_column))) %>% 
    pivot_longer(cols = c(2:(normalised_matrix %>% colnames() %>% length() + 1)), names_to = "new_sample_id", values_to = "normalised_intensity_log2") %>% 
    filter(!is.na(normalised_intensity_log2)) %>% 
    left_join(data %>% distinct(eg_precursor_id, pg_protein_accessions, pep_stripped_sequence), 
              by = c(rlang::as_name(rlang::enquo(precursor_column)))) %>% 
    left_join(data %>% distinct(new_sample_id, eg_precursor_id, intensity_log2), 
              by = c(rlang::as_name(rlang::enquo(precursor_column)), "new_sample_id")) %>% 
    left_join(data %>% distinct(new_sample_id, new_condition_id), 
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


