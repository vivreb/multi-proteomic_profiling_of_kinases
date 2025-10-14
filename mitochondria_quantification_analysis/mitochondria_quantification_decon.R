rm(list = ls())

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)


all_samples_raw_data_cmle_table <- read.csv("raw_mitochondrial_quant_data_Rab_RO.csv")

# QCs

all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class) %>% 
  ggplot(aes(x = condition, group = predicted_class, fill = predicted_class)) +
  geom_bar(position = "stack")

test_total_size_in_pixels <- all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class, size_in_pixels) %>% 
  group_by(condition, predicted_class) %>% 
  mutate(total_size_in_pixels = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  distinct(condition, predicted_class, total_size_in_pixels)

test_total_size_in_pixels_all_morphologies <- all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class, size_in_pixels) %>% 
  group_by(condition) %>% 
  mutate(total_size_in_pixels = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  distinct(condition, total_size_in_pixels)

test_mean_size_in_pixels <- all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class, size_in_pixels) %>% 
  group_by(condition) %>% 
  mutate(mean_size_in_pixels = mean(size_in_pixels)) %>% 
  ungroup() %>% 
  distinct(condition, mean_size_in_pixels)

# ANOVA

test <- all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class, size_in_pixels) %>% 
  group_by(condition, replicate) %>% 
  mutate(total_size_in_pixels_condition = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  group_by(condition, predicted_class, replicate) %>% 
  mutate(total_size_in_pixels_morphology = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  mutate(fraction_of_pixels_per_morphology = total_size_in_pixels_morphology/total_size_in_pixels_condition) %>% 
  distinct(condition, replicate, predicted_class, fraction_of_pixels_per_morphology)


anova_test <- all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class, size_in_pixels) %>% 
  group_by(condition, replicate) %>% 
  mutate(total_size_in_pixels_condition = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  group_by(condition, predicted_class, replicate) %>% 
  mutate(total_size_in_pixels_morphology = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  mutate(fraction_of_pixels_per_morphology = total_size_in_pixels_morphology/total_size_in_pixels_condition) %>% 
  distinct(condition, replicate, predicted_class, fraction_of_pixels_per_morphology) %>% 
  mutate(condition_1 = str_replace(condition, "_[A-Za-z]+", "")) %>% 
  mutate(condition_2 = str_replace(condition, "[A-Za-z]+_", "")) %>% 
  dplyr::select(-c(condition)) %>% 
  filter(predicted_class == "fragmented") %>% 
  group_by(predicted_class) %>% 
  do({
    model_results <-
      aov(fraction_of_pixels_per_morphology ~ condition_1 * condition_2,
          data = .)
    tukey_result <- TukeyHSD(model_results)
    
    tukey_Rab <- as.data.frame(tukey_result$condition_1)
    tukey_RO <- as.data.frame(tukey_result$condition_2)
    tukey_interaction <-
      as.data.frame(tukey_result$`condition_1:condition_2`)
    
    tukey_Rab$comparison <- rownames(tukey_Rab)
    tukey_RO$comparison <- rownames(tukey_RO)
    tukey_interaction$comparison <- rownames(tukey_interaction)
    
    tibble(
      Rab_comparisons = list(tukey_Rab),
      RO_comparisons = list(tukey_RO),
      interaction_comparisons = list(tukey_interaction)
    )
  }) %>% 
  pivot_longer(cols = c(Rab_comparisons:interaction_comparisons), names_to = "comparison_type", values_to = "tukey_result") %>% 
  unnest(tukey_result) %>% 
  filter(!is.na(diff))

bars.df <- data.frame(x = c(1, 2, 0.8, 2.8, 1.5, 3.5, 1.5), y = c(0.75, 0.85, 0.95, 0.95, 0.95, 0.95, 1), 
                      xend = c(3, 4, 2.2, 4.2, 1.5, 3.5, 3.5), yend = c(0.75, 0.85, 0.95, 0.95, 1, 1, 1))

label.df <- data.frame(y = c(0.8, 0.9, 1.05), 
                       x = c(2, 3, 2.5),
                       label = c("p = 4.3e-3", "p = 1.4e-6", "p = 4.8e-9"))

label.df <- data.frame(y = c(0.8, 0.9, 1.05), 
                       x = c(2, 3, 2.5),
                       label = c(paste0("p == 4.3", "%*%", "10^{-3}"), paste0("p == 1.4", "%*%", "10^{-6}"), paste0("p == 4.8", "%*%", "10^{-9}")))


quanti_image <- all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class, size_in_pixels) %>% 
  group_by(condition, replicate) %>% 
  mutate(total_size_in_pixels_condition = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  group_by(condition, predicted_class, replicate) %>% 
  mutate(total_size_in_pixels_morphology = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  mutate(fraction_of_pixels_per_morphology = total_size_in_pixels_morphology/total_size_in_pixels_condition) %>% 
  distinct(condition, replicate, predicted_class, fraction_of_pixels_per_morphology) %>% 
  filter(predicted_class == "fragmented") %>% 
  group_by(condition) %>%
  mutate(
    mean = mean(fraction_of_pixels_per_morphology),
    sd = sd(fraction_of_pixels_per_morphology),
    n = n(),
    se = sd / sqrt(n)
  ) %>% 
  mutate(labels = ifelse(condition == "DMSO_DMSO", "Control", "RO-3306")) %>% 
  mutate(labels = ifelse(condition == "Rab_DMSO", "Rabusertib", labels)) %>% 
  mutate(labels = ifelse(condition == "Rab_RO", "Rabusertib + RO-3306", labels)) %>% 
  mutate(labels = factor(labels, levels = c("Control", "RO-3306", "Rabusertib", "Rabusertib + RO-3306"))) %>% 
  ggplot(aes(x = labels, y = fraction_of_pixels_per_morphology, group = condition)) +
  geom_crossbar(aes(ymin = mean - se, ymax = mean + se, y = mean),
                width = 0.5) +
  geom_point() +
  labs(x = "", y = "Fraction of fragmented mitochondria") +
  scale_y_continuous(limits = c(0, 1.1)) +
  geom_segment(data = bars.df, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, size = 0.5) +
  geom_text(data = label.df, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.5, parse=TRUE)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=50, hjust = 1),
        axis.title.y = element_text(size = 10))

quanti_image

# pdf("Z:/Viviane/Experiments/data_for_manuscript/mitochondria_quantification/CHEK1_inhibition_microscopy_mitochondria_quantification_deconvoluted.pdf",
#     height = 4,
#     width = 2)
# quanti_image
# dev.off()






