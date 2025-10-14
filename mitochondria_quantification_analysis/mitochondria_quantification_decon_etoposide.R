rm(list = ls())

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)


all_samples_raw_data_cmle_table <- read.csv("raw_mitochondrial_quant_data_etoposide.csv")

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

# T-test (equivalent to ANOVA on single sample)

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

test_input_df <- all_samples_raw_data_cmle_table %>% 
  distinct(condition, replicate, object_id, predicted_class, size_in_pixels) %>% 
  group_by(condition, replicate) %>% 
  mutate(total_size_in_pixels_condition = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  group_by(condition, predicted_class, replicate) %>% 
  mutate(total_size_in_pixels_morphology = sum(size_in_pixels)) %>% 
  ungroup() %>% 
  mutate(fraction_of_pixels_per_morphology = total_size_in_pixels_morphology/total_size_in_pixels_condition) %>% 
  distinct(condition, replicate, predicted_class, fraction_of_pixels_per_morphology) %>% 
  filter(predicted_class == "fragmented")

t_test_1 <- test_input_df %>% 
  filter(condition == "DMSO") %>% 
  pull(fraction_of_pixels_per_morphology)

t_test_2 <- test_input_df %>% 
  filter(condition == "Etoposide") %>% 
  pull(fraction_of_pixels_per_morphology)

t.test(t_test_1, t_test_2, paired = F)

anova_tes_decon <- test_input_df %>% do({
  model_results <-
    aov(fraction_of_pixels_per_morphology ~ condition,
        data = .)
  tukey_result <- TukeyHSD(model_results)
  
  tukey_Eto <- as.data.frame(tukey_result$condition)
  
  tukey_Eto$comparison <- rownames(tukey_Eto)
  
  tibble(
    Eto_comparisons = list(tukey_Eto)
  )
}) %>% 
  pivot_longer(cols = c(Eto_comparisons), names_to = "comparison_type", values_to = "tukey_result") %>% 
  unnest(tukey_result) %>% 
  filter(!is.na(diff))


bars.df <- data.frame(x = c(1), y = c(0.85), 
                      xend = c(2), yend = c(0.85))

label.df <- data.frame(y = c(0.9), 
                       x = c(1.5),
                       label = c(paste0("p == 0.151")))


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
  mutate(labels = factor(condition, levels = c("DMSO", "Etoposide"))) %>% 
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

# pdf("Etoposide_microscopy_mitochondria_quantification_deconvoluted_t-test.pdf",
#     height = 4,
#     width = 2)
# quanti_image
# dev.off()




