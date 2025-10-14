library(protti)
library(dplyr)
library(tidyr)


# Load data

whole_proteome_for_GO_background <- read.csv("whole_proteome_uniprot_for_GO_background.csv")

diff_abundance_CAMKK2_AP <- read.csv("diff_abundance_probe_vs_DMSO_CAMKK2_WT.csv") %>% 
  left_join((whole_proteome_for_GO_background %>% dplyr::select(accession, go_f, go_p, go_c)), by = c("pg_protein_accessions" = "accession"))

diff_abundance_CHEK1_AP <- read.csv("diff_abundance_probe_vs_DMSO_CHEK1_WT.csv") %>% 
  left_join((whole_proteome_for_GO_background %>% dplyr::select(accession, go_f, go_p, go_c)), by = c("pg_protein_accessions" = "accession"))

diff_abundance_CAMKK2_BioID <- read.csv("diff_abundance_probe_vs_DMSO_CAMKK2_Nt_EGFP_corr.csv")

diff_abundance_CHEK1_BioID <- read.csv("diff_abundance_probe_vs_DMSO_CHEK1_EGFP_corr.csv")

diff_abundance_PRKAA1_BioID <- read.csv("diff_abundance_probe_vs_DMSO_PRKAA1_EGFP_corr.csv")

diff_abundance_PRKCA_BioID <- read.csv("diff_abundance_probe_vs_DMSO_PRKCA_WT_EGFP_corr.csv")



# Gene names

gene_names <- c("CAMKK2", "CHEK1", "PRKAA1", "PRKCA")

# Experiment types

experiment_types <- list()
experiment_types[["CAMKK2"]] <- c("AP", "BioID")
experiment_types[["CHEK1"]] <- c("AP", "BioID")
experiment_types[["PRKAA1"]] <- c("BioID")
experiment_types[["PRKCA"]] <- c("BioID")
# go_terms

go_types <- c("go_f", "go_c", "go_p")

go_titles <- c("MF", "CC", "BP")
names(go_titles) <- go_types


for(gene in gene_names){
  for(experiment in experiment_types[[gene]]){
    for(go_type in go_types){
      
      df_ending <- paste(gene, experiment, sep = "_")
      
      diff_abundance <- sprintf("diff_abundance_%s", df_ending) %>% 
        get()
      
      # go enrichment compared to all proteins
      
      go_analysis_interactome <- whole_proteome_for_GO_background %>% 
        mutate(go_significance = (accession %in% (diff_abundance %>% pull(pg_protein_accessions)))) %>% 
        calculate_go_enrichment(protein_id = accession,
                                is_significant = go_significance,
                                go_annotations_uniprot = !!sym(go_type),
                                plot = FALSE)
      
      go_analysis_interactome %>% write.csv(paste(getwd(), "/GO_enrichment_against_all_proteins_", gene, "_", experiment, "_", go_type, ".csv", sep = ""))
      
      if(length(go_analysis_interactome %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 3){height = 2}
      else if(length(go_analysis_interactome %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 5){height = 3}
      else if(length(go_analysis_interactome %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 8){height = 4}
      else if(length(go_analysis_interactome %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 15){height = 5}
      else(height = 6)
      
      if(length(go_analysis_interactome %>% filter(adj_pval < 0.05) %>% pull(go_id)) > 0){
        if(length(go_analysis_interactome %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 15){
          
        go_analysis_interactome <-
          whole_proteome_for_GO_background %>%
          mutate(go_significance = (accession %in% (
            diff_abundance %>% pull(pg_protein_accessions)
          ))) %>%
          calculate_go_enrichment(
            protein_id = accession,
            is_significant = go_significance,
            go_annotations_uniprot = !!sym(go_type),
            plot_cutoff = "adj_pval 0.05",
            barplot_fill_colour = c("#c4d3e9", "#c4d3e9"),
            plot_title = paste(gene, " GO enrichment (", go_titles[[go_type]], ")", sep = "")
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none",
            title = element_blank()
          )
        
        } else {
          go_analysis_interactome <-
            whole_proteome_for_GO_background %>%
            mutate(go_significance = (accession %in% (
              diff_abundance %>% pull(pg_protein_accessions)
            ))) %>%
            calculate_go_enrichment(
              protein_id = accession,
              is_significant = go_significance,
              go_annotations_uniprot = !!sym(go_type),
              plot_cutoff = "adj_pval top10",
              barplot_fill_colour = c("#c4d3e9", "#c4d3e9"),
              plot_title = paste(gene, " GO enrichment (", go_titles[[go_type]], ")", sep = ""),
              label_move_frac = 0.25,
              enrichment_type = "enriched"
            ) +
            theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position = "none",
              title = element_blank()
            )
          
        }
        
        pdf(
          paste(
            getwd(),
            "/GO_enrichment_against_all_proteins_",
            gene,
            "_",
            experiment,
            "_",
            go_type,
            ".pdf",
            sep = ""
          ),
          width = 6,
          height = 4,
          pointsize = 20
        )
        print(go_analysis_interactome)
        dev.off()
      }
      
   
      
      # go enrichment of significant proteins
      go_cellcomp_plot <- diff_abundance %>% 
        mutate(go_significance = ifelse(diff > 1 & adj_pval < 0.05, TRUE, FALSE)) %>% 
        calculate_go_enrichment(protein_id = pg_protein_accessions,
                                is_significant = go_significance,
                                go_annotations_uniprot = !!sym(go_type),
                                plot = FALSE)
      
      go_cellcomp_plot %>% write.csv(paste(getwd(), "/GO_enrichment_signficant_proteins_", gene, "_", experiment, "_", go_type, ".csv", sep = ""))
      
      if(length(go_cellcomp_plot %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 3){height = 2}
      else if(length(go_cellcomp_plot %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 5){height = 3}
      else if(length(go_cellcomp_plot %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 8){height = 4}
      else if(length(go_cellcomp_plot %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 11){height = 5}
      else(height = 6)
      
      
      
      if(length(go_cellcomp_plot %>% filter(adj_pval < 0.05) %>% pull(go_id)) > 0){
        if(length(go_cellcomp_plot %>% filter(adj_pval < 0.05) %>% pull(go_id)) < 15){
          
        go_cellcomp_plot <- diff_abundance %>%
          mutate(go_significance = ifelse(diff > 1 &
                                            adj_pval < 0.05, TRUE, FALSE)) %>%
          calculate_go_enrichment(
            protein_id = pg_protein_accessions,
            is_significant = go_significance,
            go_annotations_uniprot = !!sym(go_type),
            plot_cutoff = "adj_pval 0.05",
            barplot_fill_colour = c("#f6e7a9", "#f6e7a9"),
            plot_title = paste(gene, " GO enrichment (", go_titles[[go_type]], ")", sep = "")
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none"
          )
        } else{
          
          go_cellcomp_plot <- diff_abundance %>%
            mutate(go_significance = ifelse(diff > 1 &
                                              adj_pval < 0.05, TRUE, FALSE)) %>%
            calculate_go_enrichment(
              protein_id = pg_protein_accessions,
              is_significant = go_significance,
              go_annotations_uniprot = !!sym(go_type),
              plot_cutoff = "adj_pval top15",
              barplot_fill_colour = c("#f6e7a9", "#f6e7a9"),
              plot_title = paste(gene, " GO enrichment (", go_titles[[go_type]], ")", sep = "")
            ) +
            theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position = "none"
            )
          
        }
        
        pdf(
          paste(
            getwd(),
            "/GO_enrichment_signficant_proteins_",
            gene,
            "_",
            experiment,
            "_",
            go_type,
            ".pdf",
            sep = ""
          ),
          width = 6,
          height = height,
          pointsize = 20
        )
        print(go_cellcomp_plot)
        dev.off()
      }
 
    }
  }
}


# Venn diagrams

gene_names <- c("CAMKK2", "CHEK1")

for(gene in gene_names){
      
  df_ending <- paste(gene, "AP", sep = "_")
  
  diff_abundance_AP <- sprintf("diff_abundance_%s", df_ending) %>%
    get()
  
  df_ending <- paste(gene, "BioID", sep = "_")
  
  diff_abundance_BioID <- sprintf("diff_abundance_%s", df_ending) %>%
    get()
  
  # All proteins overlap
  
  venn_list <- list()
  venn_list[["Affinity purification"]] <- diff_abundance_AP %>% pull(pg_protein_accessions)
  venn_list[["Proximity labeling"]] <- diff_abundance_BioID %>% pull(pg_protein_accessions)
  
  venn_diagram <- ggvenn(venn_list, show_percentage = FALSE, fill_color = c("#c4d3e9", "#7788aa"), set_name_size = 5)
  
  pdf(
    paste(getwd(),  "/venn_diagram_quantified_proteins_", gene, ".pdf", sep = ""), 
    width = 4,
    height = 3,
    pointsize = 20
  )
  print(venn_diagram)
  dev.off()
  
  protein_list <- diff_abundance_AP %>% 
    dplyr::select(pg_protein_accessions, gene_primary) %>% 
    mutate(detected_in_BioID = ifelse(pg_protein_accessions %in% (diff_abundance_BioID %>% pull(pg_protein_accessions)), TRUE, FALSE))
  
  protein_list %>%  write_csv(paste(getwd(),  "/venn_diagram_quantified_proteins_", gene, ".csv", sep = ""))
  
  # Significant proteins overlap
  
  venn_list <- list()
  venn_list[["Affinity purification"]] <- diff_abundance_AP %>% filter(significant == TRUE) %>% pull(pg_protein_accessions)
  venn_list[["Proximity labeling"]] <- diff_abundance_BioID %>% filter(significant == TRUE) %>% pull(pg_protein_accessions)
  
  venn_diagram <- ggvenn(venn_list, show_percentage = FALSE, fill_color = c("#f6e7a9", "#ddbb77"), set_name_size = 5)
  
  pdf(
    paste(getwd(),  "/venn_diagram_significant_proteins_", gene, ".pdf", sep = ""), 
    width = 4,
    height = 3,
    pointsize = 20
  )
  print(venn_diagram)
  dev.off()
  
  protein_list <- diff_abundance_AP %>% 
    filter(significant == TRUE) %>% 
    dplyr::distinct(pg_protein_accessions, gene_primary) %>% 
    mutate(significant_in = "AP") %>% 
    rbind(diff_abundance_BioID %>% 
            filter(significant == TRUE) %>% 
            dplyr::distinct(pg_protein_accessions, gene_primary) %>% 
            mutate(significant_in = "BioID")) %>% 
    group_by(pg_protein_accessions) %>% 
    mutate(significant_in = ifelse(n() == 2, "AP, BioID", significant_in)) %>% 
    distinct()
  
  protein_list %>%  write_csv(paste(getwd(),  "/venn_diagram_significant_proteins_", gene, ".csv", sep = ""))
  
  # BioID significant proteins in AP interactors overlap
  
  venn_list <- list()
  venn_list[["Affinity purification"]] <- diff_abundance_AP %>% pull(pg_protein_accessions)
  venn_list[["Proximity labeling"]] <- diff_abundance_BioID %>% filter(significant == TRUE) %>% pull(pg_protein_accessions)
  
  venn_diagram <- ggvenn(venn_list, show_percentage = FALSE, fill_color = c("#c4d3e9", "#f6e7a9"), set_name_size = 5)
  
  pdf(
    paste(getwd(),  "/venn_diagram_significant_BioID_proteins_quantified_AP_interactors_", gene, ".pdf", sep = ""), 
    width = 4,
    height = 3,
    pointsize = 20
  )
  print(venn_diagram)
  dev.off()
  
  protein_list <- diff_abundance_BioID %>%
    filter(significant == TRUE) %>% 
    dplyr::distinct(pg_protein_accessions, gene_primary) %>% 
    mutate(detected_in_AP = ifelse(pg_protein_accessions %in% (diff_abundance_AP %>% pull(pg_protein_accessions)), TRUE, FALSE))
  
  protein_list %>%  write_csv(paste(getwd(),  "/venn_diagram_significant_BioID_proteins_quantified_AP_interactors_", gene, ".csv", sep = ""))
  
}

## Compare PRKAA1 and CAMKK2



# All proteins overlap

venn_list <- list()
venn_list[["CAMKK2"]] <- diff_abundance_CAMKK2_BioID %>% pull(pg_protein_accessions)
venn_list[["PRKAA1"]] <- diff_abundance_PRKAA1_BioID %>% pull(pg_protein_accessions)

venn_diagram <- ggvenn(venn_list, show_percentage = FALSE, fill_color = c("#c4d3e9", "#7788aa"), set_name_size = 5)

pdf(
  paste(getwd(),  "/venn_diagram_quantified_proteins_CAMKK2_PRKAA1.pdf", sep = ""), 
  width = 4,
  height = 3,
  pointsize = 20
)
print(venn_diagram)
dev.off()


# Significant proteins overlap

venn_list <- list()
venn_list[["CAMKK2"]] <- diff_abundance_CAMKK2_BioID %>% filter(significant == TRUE) %>% pull(pg_protein_accessions)
venn_list[["PRKAA1"]] <- diff_abundance_PRKAA1_BioID %>% filter(significant == TRUE) %>% pull(pg_protein_accessions)

venn_diagram <- ggvenn(venn_list, show_percentage = FALSE, fill_color = c("#f6e7a9", "#ddbb77"), set_name_size = 5)

pdf(
  paste(getwd(),  "/venn_diagram_significant_proteins_CAMKK2_PRKAA1.pdf", sep = ""), 
  width = 4,
  height = 3,
  pointsize = 20
)
print(venn_diagram)
dev.off()

