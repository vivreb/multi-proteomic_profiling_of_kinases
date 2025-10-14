library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(protti)

source("BioID_functions.R")

DIA_raw_PRKCA_AP <- read_protti("VR_Ex81_PRKCA_WT_Report_VR_QC (Normal).tsv") %>% 
          mutate(bait_label = "PRKCA_WT")


unis <- DIA_raw_PRKCA_AP %>%
  dplyr::pull(pg_protein_accessions) %>%
  strsplit(";") %>%
  unlist() %>%
  unique()

  uniprot <-
    fetch_uniprot(
      unis,
      columns = c(
        "protein_name",
        "length",
        "sequence",
        "protein_name",
        "gene_primary")
    ) %>%
    dplyr::rename(
      protein_sequence = sequence,
      length_protein = length,
      pg_protein_accessions = accession
    )

DIA_raw_uniprot <- DIA_raw_PRKCA_AP %>% 
  left_join(uniprot, by = c("pg_protein_accessions")) %>% 
  find_peptide(protein_sequence, pep_stripped_sequence) %>%
  assign_peptide_type(aa_before, last_aa)

DIA_clean <- filter_data_and_log2_transform(data = DIA_raw_uniprot,
                                            proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"))

DIA_clean_normalised <- normalise_data(data = DIA_clean,
                                       imputation_method = "none")

DIA_final_AP <- DIA_clean_normalised %>% 
  filter(pg_protein_accessions == "P17252")

protein_abundance_AP <- calculate_protein_abundance(
  data = DIA_final_AP,
  sample = new_sample_id,
  protein_id = pg_protein_accessions,
  precursor = eg_precursor_id,
  peptide = pep_stripped_sequence,
  intensity_log2 = intensity_log2,
  method = "iq",
  for_plot = FALSE,
  retain_columns = c(new_condition_id)
)


DIA_raw_PRKCA_BioID <- read_protti("VR_Ex82_PRKCA_WT_Report_VR_QC (Normal).tsv") %>%
  mutate(bait_label = "PRKCA_WT") %>% 
  filter(pg_protein_accessions == "P17252") %>% 
  mutate(intensity_log2 = log2(fg_quantity))

unis <- DIA_raw_PRKCA_BioID %>%
  dplyr::pull(pg_protein_accessions) %>%
  strsplit(";") %>%
  unlist() %>%
  unique()

uniprot <-
  fetch_uniprot(
    unis,
    columns = c(
      "protein_name",
      "length",
      "sequence",
      "protein_name",
      "gene_primary")
  ) %>%
  dplyr::rename(
    protein_sequence = sequence,
    length_protein = length,
    pg_protein_accessions = accession
  )

DIA_raw_uniprot <- DIA_raw_PRKCA_BioID %>% 
  left_join(uniprot, by = c("pg_protein_accessions")) %>% 
  find_peptide(protein_sequence, pep_stripped_sequence) %>%
  assign_peptide_type(aa_before, last_aa)

DIA_clean <- filter_data_and_log2_transform(data = DIA_raw_uniprot,
                                            proteins_exempt_from_proteotypicity_filter = c(unname(bait_accession), "P0DP23"))

DIA_clean_normalised <- normalise_data(data = DIA_clean,
                                       imputation_method = "none")

DIA_final_BioID <- DIA_clean_normalised %>% 
  filter(pg_protein_accessions == "P17252")


protein_abundance_BioID <- calculate_protein_abundance(
  data = DIA_final_BioID,
  sample = new_sample_id,
  protein_id = pg_protein_accessions,
  precursor = eg_precursor_id,
  peptide = pep_stripped_sequence,
  intensity_log2 = intensity_log2,
  method = "iq",
  for_plot = FALSE,
  retain_columns = c(new_condition_id)
)

protein_abundance_all <- protein_abundance_AP %>% 
  mutate(experiment = "Affinity purification") %>% 
  rbind(protein_abundance_BioID %>% 
          mutate(experiment = "Proximity labeling"))

pval_AP <- t.test(protein_abundance_AP$intensity_log2[which(protein_abundance_AP$new_condition_id == "PRKCA_WT_DMSO")],
                  protein_abundance_AP$intensity_log2[which(protein_abundance_AP$new_condition_id == "PRKCA_WT_probe")],
                  alternative = "greater",
                  paired = FALSE)

pval_BioID <- t.test(protein_abundance_BioID$intensity_log2[which(protein_abundance_BioID$new_condition_id == "PRKCA_WT_DMSO")],
                     protein_abundance_BioID$intensity_log2[which(protein_abundance_BioID$new_condition_id == "PRKCA_WT_probe")],
                     alternative = "greater",
                     paired = FALSE)


bars.df <- data.frame(x = c(1), y = c(22), xend = c(2), yend = c(22))

label.df <- data.frame(experiment = c("Affinity purification", "Proximity labeling"),
                       Intensity = c(22.2, 22.2), 
                       x = c(1.5, 1.5),
                       label = c("p = 0.03875", "p = 0.9939"))

protein_abundance_all %>% 
  mutate(new_condition_id = ifelse(new_condition_id == "PRKCA_WT_DMSO", "DMSO", "Gö 6983")) %>% 
  group_by(pg_protein_accessions, new_condition_id) %>% 
  mutate(Intensity = mean(intensity_log2)) %>% 
  mutate(SEM = sd(intensity_log2) / sqrt(n())) %>% 
  ungroup() %>% 
  distinct(pg_protein_accessions, new_condition_id, experiment, Intensity, SEM, intensity_log2) %>% 
  ggplot(aes(x = new_condition_id, y = intensity_log2, fill = new_condition_id)) +
  geom_boxplot(position = position_dodge()) +
  #geom_errorbar(stat = "identity", aes(ymin=Intensity-SEM, ymax=Intensity+SEM), width=.2, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("indianred", "cadetblue")) +
  theme_bw() +
  #geom_hline(yintercept = 0) +
  labs(y = bquote("log"[2]~"(Intensity)")) + #title = paste("Changes in PKC\u03b1 Intensity in Soluble Fraction upon PKC\u03b1 Inhibition")
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8)) +
  scale_y_continuous(limits = c(17, 23)) +
  facet_grid(~experiment)

pdf("prkca_change_in_MS_intensity_AP_vs_BioID.pdf", 
    width = 4,
    height = 3,
    pointsize = 20)

dev.off()
