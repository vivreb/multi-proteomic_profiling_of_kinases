library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(protti)


make_ap_protein_count_plot <- function(diff_abundance_data, raw_data, bait = bait_accession){

  number_of_proteins_detected <- raw_data %>% pull(pg_protein_accessions) %>% unique() %>% length()
  
  number_of_interactors <- diff_abundance_data %>% pull(pg_protein_accessions) %>% unique() %>% length()
  
  number_of_proteins_significant <- diff_abundance_data %>% 
    filter(significant == TRUE) %>% 
    filter(abs(diff) > 1) %>% 
    filter(grepl("probe", comparison)) %>% 
    pull(pg_protein_accessions) %>% unique() %>% length()
  
  names_vector <- c("Proteins detected", "Number of interactors", "Significant proteins")
  
  
  
  plot_data <- data.frame("statistic" = names_vector,
                          "value" = c(number_of_proteins_detected, number_of_interactors, number_of_proteins_significant)) %>% 
    dplyr::mutate(statistic = factor(statistic, levels = names_vector))
  
  plot <- plot_data %>% 
    ggplot(aes(x = statistic, y = value, fill = statistic)) +
    geom_bar(stat="identity") +
    geom_text(aes(label = value, vjust = -0.15, size = 2), size = 2) +
    scale_fill_manual(values = c("grey90", "blue4", "#8bbd26")) +
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

make_ap_intensity_plot <- function(raw_data, inhibitor){
  
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


diff_abundance_all_CAMKK2_WT <- read_csv("diff_abundance_probe_vs_DMSO_CAMKK2_WT.csv")
diff_abundance_all_CAMKK2_R311C <- read_csv("diff_abundance_probe_vs_DMSO_CAMKK2_R311C_mut.csv")
diff_abundance_all_CAMKK2_D312A <- read_csv("diff_abundance_probe_vs_DMSO_CAMKK2_D312A_mut.csv")
diff_abundance_all_CAMKK2_T483D <- read_csv("diff_abundance_probe_vs_DMSO_CAMKK2_T483D_mut.csv")

diff_abundance_all_CHEK1_WT <- read_csv("diff_abundance_probe_vs_DMSO_CHEK1_WT.csv")
diff_abundance_all_CHEK1_D130A <- read_csv("diff_abundance_probe_vs_DMSO_CHEK1_D130A_mut.csv")
diff_abundance_all_CHEK1_L449R <- read_csv("diff_abundance_probe_vs_DMSO_CHEK1_L449R_mut.csv")

diff_abundance_all_PRKCA_WT <- read_csv("diff_abundance_probe_vs_DMSO_PRKCA_WT.csv")



DIA_raw_CAMKK2_WT <- read_protti("20241113_142710_VR_Ex75_CAMKK2_wt_Report.tsv") %>% 
          mutate(bait_label = "CAMKK2_WT")
DIA_raw_CAMKK2_R311C <- read_protti("20241108_072530_VR_Ex75_CAMKK2_R311C_Report.tsv") %>% 
          mutate(bait_label = "CAMKK2_R311C")
DIA_raw_CAMKK2_D312A <- read_protti("20241114_152135_VR_Ex75_CAMKK2_D312A_rem_Report.tsv") %>%
          mutate(bait_label = "CAMKK2_D312A")
DIA_raw_CAMKK2_T483D <- read_protti("20241115_083241_VR_Ex75_CAMKK2_T483D_Report.tsv") %>%
          mutate(bait_label = "CAMKK2_T483D")

DIA_raw_CHEK1_WT <- read_protti("VR_Ex78_CHEK1_WT_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "CHEK1_WT")
DIA_raw_CHEK1_D130A <- read_protti("VR_Ex78_CHEK1_D130A_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "CHEK1_D130A")
DIA_raw_CHEK1_L449R <- read_protti("VR_Ex78_CHEK1_L449R_Report_VR_QC (Normal).tsv") %>%
          mutate(bait_label = "CHEK1_L449R")

DIA_raw_PRKCA_WT <- read_protti("VR_Ex81_PRKCA_WT_Report_VR_QC (Normal).tsv") %>% 
          mutate(bait_label = "PRKCA_WT")


uniprot_data_CAMKK2_WT <- read_csv("uniprot_CAMKK2_WT.csv")
uniprot_data_CAMKK2_R311C <- read_csv("uniprot_CAMKK2_R311C.csv")
uniprot_data_CAMKK2_D312A <- read_csv("uniprot_CAMKK2_D312A.csv")
uniprot_data_CAMKK2_T483D <- read_csv("uniprot_CAMKK2_T483D.csv")

uniprot_data_CHEK1_WT <- read_csv("uniprot_CHEK1_WT.csv")
uniprot_data_CHEK1_D130A <- read_csv("uniprot_CHEK1_D130A.csv")
uniprot_data_CHEK1_L449R <- read_csv("uniprot_CHEK1_L449R.csv")

uniprot_data_PRKCA_WT <- read_csv("uniprot_PRKCA_WT.csv")


bait_gene_list <- c("CAMKK2_WT", "CAMKK2_R311C", "CAMKK2_D312A", "CAMKK2_T483D",
                    "CHEK1_WT", "CHEK1_D130A", "CHEK1_L449R",
                    "PRKCA_WT")

bait_accession_list <- c("Q96RR4", "Q96RR4", "Q96RR4", "Q96RR4",
                         "O14757", "O14757", "O14757",
                         "P17252")

inhibitor_list <- c("SGC-CAMKK2-1", "SGC-CAMKK2-1", "SGC-CAMKK2-1", "SGC-CAMKK2-1", 
                    "Rabusertib", "Rabusertib", "Rabusertib", 
                    "Gö 6983")

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
  
  protein_count_plot <- make_ap_protein_count_plot(diff_abundance_data = diff_abundance_data, raw_data = raw_data) #2 by 3
  
  ap_intensity_plot <- make_ap_intensity_plot(raw_data = raw_data, inhibitor = inhibitor) # 2 by 3
  
  pdf(paste("ap_protein_count_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(protein_count_plot)
  dev.off()
  
  pdf(paste("ap_intensity_plot_", bait_gene, ".pdf", sep = ""), 
      width = 2,
      height = 3,
      pointsize = 20)
  print(ap_intensity_plot)
  dev.off()
  
  
}

