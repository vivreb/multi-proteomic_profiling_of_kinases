library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(protti)

make_bioid_protein_count_plot <- function(diff_abundance_data, raw_data, bait = bait_accession){
  
  number_of_proteins_detected <- raw_data %>% pull(pg_protein_accessions) %>% unique() %>% length()
  
  number_of_proteins_quantified <- diff_abundance_data %>% pull(pg_protein_accessions) %>% unique() %>% length()
  
  number_of_proteins_significant <- diff_abundance_data %>% filter(significant == TRUE) %>% 
    pull(pg_protein_accessions) %>% unique() %>% length()
  
  names_vector <- c("Proteins detected", "Proteins quantified", "Significant proteins")
  
  
  
  plot_data <- data.frame("statistic" = names_vector,
                          "value" = c(number_of_proteins_detected, number_of_proteins_quantified, number_of_proteins_significant)) %>% 
    dplyr::mutate(statistic = factor(statistic, levels = names_vector))
  
  plot <- plot_data %>% 
    ggplot(aes(x = statistic, y = value, fill = statistic)) +
    geom_bar(stat="identity") +
    geom_text(aes(label = value, vjust = -0.15, size = 2), size = 2) +
    scale_fill_manual(values = c("grey90", "grey70", "#8bbd26")) +
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

make_bioid_intensity_plot <- function(raw_data, inhibitor){
  
  if(inhibitor == "Staurosporine"){filter_vector = c("DMSO", "STS")} 
  else filter_vector = c("DMSO", "probe")
  
  data <- raw_data %>% 
    filter(r_condition %in% filter_vector) %>% 
    dplyr::mutate(intensity_log2 = log2(fg_quantity)) %>% 
    dplyr::distinct(r_condition, r_replicate, intensity_log2) %>% 
    dplyr::mutate(r_condition = ifelse(r_condition != "DMSO", inhibitor, r_condition)) %>% 
    dplyr::mutate(x_label = paste(r_condition, r_replicate)) %>% 
    dplyr::mutate(x_label = factor(x_label, levels = c("DMSO 1", "DMSO 2", "DMSO 3", "DMSO 4",
                                                       paste(inhibitor, "1"), paste(inhibitor, "2"), paste(inhibitor, "3"), paste(inhibitor, "4"))))
  
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


diff_abundance_all_CAMKK2_WT <- read_csv("diff_abundance_probe_vs_DMSO_CAMKK2_Nt_EGFP_corr.csv")
diff_abundance_all_PRKAA1 <- read_csv("diff_abundance_probe_vs_DMSO_PRKAA1_EGFP_corr.csv")
diff_abundance_all_CHEK1_WT <- read_csv("diff_abundance_probe_vs_DMSO_CHEK1_EGFP_corr.csv")
diff_abundance_all_PRKCA_WT <- read_csv("diff_abundance_probe_vs_DMSO_PRKCA_WT_EGFP_corr.csv")
diff_abundance_all_PRKCA_D463A <- read_csv("diff_abundance_probe_vs_DMSO_PRKCA_D463A_EGFP_corr.csv")
diff_abundance_all_PRKCA_R22A_A25E <- read_csv("diff_abundance_probe_vs_DMSO_PRKCA_R22A_A25E_EGFP_corr.csv")



DIA_raw_CAMKK2_WT <- read_protti("SL_Ex1_Nt_CAMKK2_Report_VR_QC (Normal).tsv") %>% 
  mutate(bait_label = "CAMKK2_Nt")

DIA_raw_PRKAA1 <- read_protti("SL_Ex1_PRKAA1_Report_VR_QC (Normal).tsv") %>% 
  mutate(bait_label = "PRKAA1")

DIA_raw_CHEK1_WT <- read_protti("VR_Ex90_CHEK1_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "CHEK1")

DIA_raw_PRKCA_WT <- read_protti("VR_Ex82_PRKCA_WT_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_WT")

DIA_raw_PRKCA_R22A_A25E <- read_protti("VR_Ex82_PRKCA_R22A_A25E_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_R22A_A25E") 

DIA_raw_PRKCA_D463A <- read_protti("VR_Ex82_PRKCA_D463A_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_D463A")


uniprot_data_CAMKK2_WT <- read_csv("uniprot_CAMKK2_Nt.csv")
uniprot_data_PRKAA1 <- read_csv("uniprot_PRKAA1.csv")
uniprot_data_CHEK1_WT <- read_csv("uniprot_CHEK1.csv")
uniprot_data_PRKCA_WT <- read_csv("uniprot_PRKCA_WT.csv")
uniprot_data_PRKCA_R22A_A25E <- read_csv("uniprot_PRKCA_R22A_A25E.csv")
uniprot_data_PRKCA_D463A <- read_csv("uniprot_PRKCA_D463A.csv")



bait_gene_list <- c("CAMKK2_WT",
                    "PRKAA1",
                    "CHEK1_WT", 
                    "PRKCA_WT", "PRKCA_D463A", "PRKCA_R22A_A25E")

bait_accession_list <- c("Q96RR4",
                         "Q13131",
                         "O14757",
                         "P17252", "P17252", "P17252")

inhibitor_list <- c("SGC-CAMKK2-1",
                    "SGC-CAMKK2-1",
                    "Rabusertib",
                    "Gö 6983", "Gö 6983", "Gö 6983")

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
  
  protein_count_plot <- make_bioid_protein_count_plot(diff_abundance_data = diff_abundance_data, raw_data = raw_data) #2 by 3
  
  bioid_intensity_plot <- make_bioid_intensity_plot(raw_data = raw_data, inhibitor = inhibitor) # 2 by 3
  
  pdf(paste("bioid_protein_count_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(protein_count_plot)
  dev.off()
  
  pdf(paste("bioid_intensity_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(bioid_intensity_plot)
  dev.off()
  
  
}

