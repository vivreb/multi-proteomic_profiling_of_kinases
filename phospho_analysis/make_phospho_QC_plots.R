library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(protti)


make_phospho_peptide_protein_count_plot <- function(diff_abundance_data, raw_data, bait = bait_accession){

  number_of_proteins_detected <- raw_data %>% pull(pg_protein_accessions) %>% unique() %>% length()
  
  number_of_proteins_with_phosphosite <- diff_abundance_data %>% pull(pg_protein_accessions) %>% unique() %>% length()
  
  number_of_proteins_significant <- diff_abundance_data %>% filter(significant == TRUE) %>% 
    pull(pg_protein_accessions) %>% unique() %>% length()
  
  
  number_of_peptides_detected <- raw_data %>% pull(eg_precursor_id) %>% unique() %>% length()
  
  number_of_peptides_with_phosphosite <- diff_abundance_data %>% pull(eg_precursor_id) %>% unique() %>% length()
  
  number_of_peptides_bait <- raw_data %>% filter(pg_protein_accessions == bait) %>% 
    pull(eg_precursor_id) %>% unique() %>% length()
  
  number_of_peptides_bait_with_phosphosite <- diff_abundance_data %>% filter(pg_protein_accessions == bait) %>% 
    pull(eg_precursor_id) %>% unique() %>% length()
  
  number_of_peptides_significant <- diff_abundance_data %>% filter(significant == TRUE) %>% 
    pull(eg_precursor_id) %>% unique() %>% length()
  
  number_of_peptides_bait_significant <- diff_abundance_data %>% filter(significant == TRUE) %>% 
    filter(pg_protein_accessions == bait) %>% 
    pull(eg_precursor_id) %>% unique() %>% length()
  
  names_vector <- c("Proteins detected", "Phosphoproteins detected", "Significant proteins",
                    "Peptides detected", "Phosphopeptides detected", "Target peptides detected", "Target phosphopep. detected", 
                    "Significant peptides", "Significant target peptides")
  
  
  
  plot_data <- data.frame("statistic" = names_vector,
                          "value" = c(number_of_proteins_detected, number_of_proteins_with_phosphosite, number_of_proteins_significant,
                                      number_of_peptides_detected, number_of_peptides_with_phosphosite, number_of_peptides_bait, 
                                      number_of_peptides_bait_with_phosphosite, number_of_peptides_significant, number_of_peptides_bait_significant)) %>% 
    dplyr::mutate(statistic = factor(statistic, levels = names_vector))
  
  plot <- plot_data %>% 
    ggplot(aes(x = statistic, y = value, fill = statistic)) +
    geom_bar(stat="identity") +
    geom_text(aes(label = value, vjust = -0.15, size = 2), size = 2) +
    scale_fill_manual(values = c("grey90", "grey90", "grey90", "grey90", "grey90", "#5680C1", "lightblue", "grey90", "goldenrod")) +
    labs(y = "Count") +
    scale_y_continuous(limits = c(0, 25500)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=50, hjust = 1))
  
  return(plot)
}

make_phospho_intensity_plot <- function(raw_data, inhibitor){
  
  if(inhibitor == "Staurosporine"){filter_vector = c("DMSO", "STS")} 
  else filter_vector = c("DMSO", "probe")
  
  data <- raw_data %>% 
    filter(r_condition %in% filter_vector) %>% 
    dplyr::mutate(intensity_log2 = log2(fg_quantity)) %>% 
    dplyr::distinct(r_condition, r_replicate, intensity_log2) %>% 
    dplyr::mutate(r_condition = ifelse(r_condition != "DMSO", inhibitor, r_condition)) %>% 
    dplyr::mutate(x_label = paste(r_condition, r_replicate)) %>% 
    dplyr::mutate(x_label = factor(x_label, levels = c("DMSO 1", "DMSO 2", "DMSO 3",
                                                       paste(inhibitor, "1"), paste(inhibitor, "2"), paste(inhibitor, "3"))))
  
  plot <- data %>% 
    ggplot(aes(x = x_label, y = intensity_log2, fill = r_condition)) +
    geom_violin() +
    scale_fill_manual(values = c("#c4d3e9", "#fee55d")) +
    labs(y = expression(paste("Normalized log "[2]*"intensity"))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=50, hjust = 1))
  
  return(plot)
  
}


diff_abundance_all_CHEK1_WT <- read_csv("diff_abundance_phospho_probe_vs_DMSO_CHEK1_WT.csv")

DIA_raw_CHEK1_WT <- read_protti("VR_Ex57_CHEK1_Report_VR_QC (Normal).tsv") %>% 
  mutate(bait_label = "CHEK1_WT")

uniprot_data_CHEK1_WT <- read_csv("uniprot_CHEK1_WT.csv")



bait_gene_list <- c("CHEK1_WT")

bait_accession_list <- c("O14757")

inhibitor_list <- c("Rabusertib")

names(bait_accession_list) <- c(bait_gene_list)
names(inhibitor_list) <- c(bait_gene_list)


for(bait_gene in bait_gene_list) {
  
  set.seed(123)
  
  raw_data <- sprintf("DIA_raw_%s", bait_gene) %>%
    get()
  
  diff_abundance_data <- sprintf("diff_abundance_all_%s", bait_gene) %>%
    get()
  
  uniprot_data <- sprintf("uniprot_data_%s", bait_gene) %>%
    get()
  
  bait_accession <- bait_accession_list[bait_gene]
  inhibitor <- inhibitor_list[bait_gene]
  
  protein_count_plot <- make_phospho_peptide_protein_count_plot(diff_abundance_data = diff_abundance_data, raw_data = raw_data) #2 by 3
  
  phospho_intensity_plot <- make_phospho_intensity_plot(raw_data = raw_data, inhibitor = inhibitor) # 2 by 3
  

  pdf(paste("phospho_protein_count_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(protein_count_plot)
  dev.off()
  
  pdf(paste("phospho_intensity_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(phospho_intensity_plot)
  dev.off()
  
}

