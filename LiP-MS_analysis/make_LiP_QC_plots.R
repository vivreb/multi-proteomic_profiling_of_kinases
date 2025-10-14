library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(protti)


make_lip_peptide_protein_count_plot <- function(diff_abundance_data, bait = bait_accession){
  
  data <- diff_abundance_data %>% 
    dplyr::distinct(eg_precursor_id, pg_protein_accessions, significant)
  
  number_of_proteins <- data %>% pull(pg_protein_accessions) %>% unique() %>% length()
  
  number_of_proteins_significant <- data %>% filter(significant == TRUE) %>% 
    pull(pg_protein_accessions) %>% unique() %>% length()
  
  
  number_of_peptides <- data %>% pull(eg_precursor_id) %>% unique() %>% length()

  number_of_peptides_bait <- data %>% filter(pg_protein_accessions == bait) %>% 
    pull(eg_precursor_id) %>% unique() %>% length()
  
  number_of_peptides_significant <- data %>% filter(significant == TRUE) %>% 
    pull(eg_precursor_id) %>% unique() %>% length()

  number_of_peptides_bait_significant <- data %>% filter(significant == TRUE) %>% 
    filter(pg_protein_accessions == bait) %>% 
    pull(eg_precursor_id) %>% unique() %>% length()
  
  names_vector <- c("Proteins detected", "Significant proteins", 
                    "Peptides detected", "Target peptides detected", "Significant peptides", "Significant target peptides")
  
  
  
  plot_data <- data.frame("statistic" = names_vector,
                          "value" = c(number_of_proteins, number_of_proteins_significant,
                                      number_of_peptides, number_of_peptides_bait, number_of_peptides_significant, number_of_peptides_bait_significant)) %>% 
    dplyr::mutate(statistic = factor(statistic, levels = names_vector))
  
  plot <- plot_data %>% 
    ggplot(aes(x = statistic, y = value, fill = statistic)) +
    geom_bar(stat="identity") +
    geom_text(aes(label = value, vjust = -0.15, size = 2), size = 2) +
    scale_fill_manual(values = c("grey90", "grey90", "grey90", "#5680C1", "grey90", "#D9353C")) +
    labs(y = "Count") +
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

make_lip_intensity_plot <- function(raw_data, inhibitor){
  
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
  
make_lip_peptype_plot <- function(raw_data, uniprot_data, inhibitor){
  
  if(inhibitor == "Staurosporine"){filter_vector = c("DMSO", "STS")} 
  else filter_vector = c("DMSO", "probe")
  
  data <- raw_data %>% 
    filter(r_condition %in% filter_vector) %>% 
    left_join(uniprot_data, by = c("pg_protein_accessions" = "uniprot_id")) %>% 
    find_peptide(protein_sequence = protein_sequence, peptide_sequence = pep_stripped_sequence) %>% 
    assign_peptide_type() %>% 
    dplyr::distinct(r_condition, r_replicate, pep_type, eg_precursor_id) %>% 
    group_by(r_condition, r_replicate) %>%
    count(pep_type) %>%
    ungroup() %>%
    dplyr::mutate(r_condition = ifelse(r_condition != "DMSO", inhibitor, r_condition)) %>% 
    dplyr::mutate(x_label = paste(r_condition, r_replicate)) %>% 
    dplyr::mutate(x_label = factor(x_label, levels = c("DMSO 1", "DMSO 2", "DMSO 3",
                                                       paste(inhibitor, "1"), paste(inhibitor, "2"), paste(inhibitor, "3")))) %>% 
    filter(!is.na(pep_type))
  
  plot <- data %>% 
    ggplot(aes(x = x_label, y = n, fill = pep_type)) +
    geom_bar(position = "fill", stat="identity") +
    scale_fill_manual(values = c("#1a5276", "grey90", "#d4ac0d")) +
    labs(y = expression(paste("Proportion of peptide type")), fill = "Peptide type") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=50, hjust = 1))
  
  return(plot)
  
}

diff_abundance_all_CAMKK2_WT <- read_csv("diff_abundance_all_CAMKK2_WT_probe.csv")
diff_abundance_all_CAMKK2_R311C <- read_csv("diff_abundance_all_CAMKK2_R311C_probe.csv")
diff_abundance_all_CAMKK2_D312A <- read_csv("diff_abundance_all_CAMKK2_D312A_probe.csv")
diff_abundance_all_CAMKK2_T483D <- read_csv("diff_abundance_all_CAMKK2_T483D_probe.csv")
diff_abundance_all_DCLK1 <- read_csv("diff_abundance_all_DCLK1_probe.csv")
diff_abundance_all_CHEK1_WT <- read_csv("diff_abundance_all_CHEK1_WT_probe.csv")
diff_abundance_all_PRKCA_WT <- read_csv("diff_abundance_all_PRKCA_WT_probe.csv")


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

uniprot_data_CAMKK2_WT <- read_csv("uniprot_CAMKK2_WT.csv")
uniprot_data_CAMKK2_R311C <- read_csv("uniprot_CAMKK2_R311C.csv")
uniprot_data_CAMKK2_D312A <- read_csv("uniprot_CAMKK2_D312A.csv")
uniprot_data_CAMKK2_T483D <- read_csv("uniprot_CAMKK2_T483D.csv")

uniprot_data_CHEK1_WT <- read_csv("uniprot_CHEK1_WT.csv")
uniprot_data_PRKCA_WT <- read_csv("uniprot_PRKCA_WT.csv")
uniprot_data_DCLK1 <- read_csv("uniprot_DCLK1.csv")



bait_gene_list <- c("CAMKK2_WT", "CAMKK2_R311C", "CAMKK2_D312A", "CAMKK2_T483D",
                    "CHEK1_WT", 
                    "PRKCA_WT",
                    "DCLK1")

bait_accession_list <- c("Q96RR4", "Q96RR4", "Q96RR4", "Q96RR4",
                         "O14757",
                         "P17252", 
                         "O15075")

inhibitor_list <- c("SGC-CAMKK2-1", "SGC-CAMKK2-1", "SGC-CAMKK2-1", "SGC-CAMKK2-1", 
               "Rabusertib",
               "Gö 6983",
               "DCLK1-IN-1")

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
  
  protein_count_plot <- make_lip_peptide_protein_count_plot(diff_abundance_data = diff_abundance_data) #2 by 3
  
  lip_intensity_plot <- make_lip_intensity_plot(raw_data = raw_data, inhibitor = inhibitor) # 2 by 3
  
  lip_peptype_plot <- make_lip_peptype_plot(raw_data = raw_data, uniprot_data = uniprot_data, inhibitor = inhibitor) # 3 by 2
  
  pdf(paste("lip_protein_count_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(protein_count_plot)
  dev.off()
  
  pdf(paste("lip_intensity_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(lip_intensity_plot)
  dev.off()
  
  pdf(paste("lip_peptype_plot_", bait_gene, ".pdf", sep = ""), 
      width = 4,
      height = 3,
      pointsize = 20)
  print(lip_peptype_plot)
  dev.off()
  
}

