### Peatland Phenolics Project - Diversity 

## library: 
suppressPackageStartupMessages({  
  library(tidyverse)
  library(vegan)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
}) 

# =========================================================
# 0) FONT SIZE STANDARDS (Springer/Biogeochemistry: 8-12 pt final size)
# =========================================================

BASE_FONT_PT  <- 9    # axis text, axis titles, legend text
TITLE_FONT_PT <- 10   # plot titles
TAG_FONT_PT   <- 12   # bold panel tags (a, b)

## ---- Consistent Type colors ----
type_colors <- c(
  "Bog"           = "#8B4513",
  "Shrub Wetland" = "#2F8B8B",
  "Fen"           = "#F5DEB3"
)

## ---- Common site order + extended names for axis display ----
site_levels <- c("S1", "BLF", "S3", "US2", "LC")

site_labels <- c(
  "S1"  = "S1-bog",
  "BLF" = "BLF-poor fen",
  "S3"  = "S3-fen",
  "US2" = "US2-fen",
  "LC"  = "LC-shrub wetland"
)

## ---- Shared theme (transparent, standardized fonts) 
my_theme_div <- theme_bw(base_size = BASE_FONT_PT) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = TITLE_FONT_PT),
    axis.text.x   = element_text(face = "bold", size = BASE_FONT_PT, color = "black",
                                 angle = 45, hjust = 1),
    axis.text.y   = element_text(size = BASE_FONT_PT, color = "black"),
    axis.title.y  = element_text(size = BASE_FONT_PT, face = "bold"),
    axis.title.x  = element_text(size = BASE_FONT_PT, face = "bold"),
    strip.text    = element_text(face = "bold", size = BASE_FONT_PT),
    legend.title  = element_text(face = "bold", size = BASE_FONT_PT),
    legend.text   = element_text(size = BASE_FONT_PT),
    panel.grid    = element_blank(),
    
    # transparent backgrounds
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    
    legend.key = element_rect(fill = "transparent", color = NA)
  )

# RP:  --------------------------------------------------------------------
compounds_table_RP <- read.csv("Tables/compounds_table_non_real_gap_filled_RP.csv")
compounds_table_RP <- compounds_table_RP %>% 
  select(- gap_status.x, -gap_status.y, -Type, - Site, - Date, - Depth, - Subsite, - rep, - pH, -WT, -WT_Level) %>%
  pivot_wider(names_from = 'SampleID', values_from = 'AUC')

# Import metadata and fix names
metadata <- read.csv('Tables/fixed_metadata_rp.csv') %>%
  mutate(Site = recode(Site, "T1" = "S1"))     # T1 -> S1, consistent with other figures

## ---- Site -> Type 
site_type <- metadata %>%            
  select(Site, Type) %>%
  distinct()

## getting intensity matrix
intensity_matrix <- compounds_table_RP %>%
  select(FeatureID, all_of(metadata$SampleID)) %>%
  pivot_longer(!FeatureID, names_to = 'SampleID', values_to = 'intensity') %>%
  arrange(desc(intensity)) %>%
  distinct(FeatureID, SampleID, .keep_all = TRUE)%>%
  pivot_wider(names_from = 'FeatureID', values_from = 'intensity')%>%
  column_to_rownames(var = "SampleID")

intensity_matrix[is.na(intensity_matrix)] <- 0

# Sum normalize intensities
sample_sum = rowSums(intensity_matrix)
norm_intensity_matrix <- intensity_matrix / sample_sum

# Diversity Index

Shannon_table_RP <- tibble(SampleID = rownames(norm_intensity_matrix),
                           Shannon = diversity(norm_intensity_matrix, index = 'shannon'))

Shannon_table_RP <- tibble(
  SampleID = rownames(norm_intensity_matrix),
  Shannon = diversity(norm_intensity_matrix, index = "shannon")
) %>%
  pivot_longer(!SampleID, names_to = "index", values_to = "values") %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

## get significance
stat_table_RP <- Shannon_table_RP %>%
  filter(!is.na(Site)) %>%
  select(values, Site) %>%
  wilcox_test(values ~ Site) %>%
  add_xy_position(step.increase = 0.1) %>%
  add_significance()

Shannon_plot_RP <- ggplot(Shannon_table_RP, aes(x = Site, y = values, fill = Type)) +
  geom_boxplot(outlier.size = 1, width = 0.7) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  scale_x_discrete(labels = site_labels) +
  stat_pvalue_manual(stat_table_RP, inherit.aes = FALSE, hide.ns = TRUE,
                     size = (BASE_FONT_PT - 1) / .pt) +
  labs(title = "Metabolite Diversity - RP", y = "Shannon diversity", x = NULL) +
  my_theme_div +
  guides(fill = "none")

Shannon_plot_RP


# HILIC -------------------------------------------------------------------
compounds_table_HILIC <- read.csv("Tables/compounds_table_non_real_gap_filled_HILIC.csv")
compounds_table_HILIC <- compounds_table_HILIC %>% 
  select(- gap_status.x, -gap_status.y, -Type, - Site, - Date, - Depth, - Subsite, - rep, - pH, -WT, -WT_Level) %>%
  pivot_wider(names_from = 'SampleID', values_from = 'AUC')

# Import metadata and fix names
metadata <- read.csv('Tables/fixed_metadata_hilic.csv') %>%
  mutate(Site = recode(Site, "T1" = "S1"))     # T1 -> S1

## ---- Site -> Type 
site_type <- metadata %>%            
  select(Site, Type) %>%
  distinct()

## getting intensity matrix
intensity_matrix <- compounds_table_HILIC %>%
  select(FeatureID, all_of(metadata$SampleID)) %>%
  pivot_longer(!FeatureID, names_to = 'SampleID', values_to = 'intensity') %>%
  arrange(desc(intensity)) %>%
  distinct(FeatureID, SampleID, .keep_all = TRUE)%>%
  pivot_wider(names_from = 'FeatureID', values_from = 'intensity')%>%
  column_to_rownames(var = "SampleID")

intensity_matrix[is.na(intensity_matrix)] <- 0

# Sum normalize intensities
sample_sum = rowSums(intensity_matrix)
norm_intensity_matrix <- intensity_matrix / sample_sum

# Diversity Index

Shannon_table_HILIC <- tibble(SampleID = rownames(norm_intensity_matrix),
                              Shannon = diversity(norm_intensity_matrix, index = 'shannon'))

Shannon_table_HILIC <- tibble(
  SampleID = rownames(norm_intensity_matrix),
  Shannon = diversity(norm_intensity_matrix, index = "shannon")
) %>%
  pivot_longer(!SampleID, names_to = "index", values_to = "values") %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

## get significance
stat_table_HILIC <- Shannon_table_HILIC %>%
  filter(!is.na(Site)) %>%
  select(values, Site) %>%
  wilcox_test(values ~ Site) %>%
  add_xy_position(step.increase = 0.1) %>%
  add_significance()

Shannon_plot_HILIC <- ggplot(Shannon_table_HILIC, aes(x = Site, y = values, fill = Type)) +
  geom_boxplot(outlier.size = 1, width = 0.7) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  scale_x_discrete(labels = site_labels) +
  stat_pvalue_manual(stat_table_HILIC, inherit.aes = FALSE, hide.ns = TRUE,
                     size = (BASE_FONT_PT - 1) / .pt) +
  labs(title = "Metabolite Diversity - HILIC", y = "Shannon diversity", x = NULL) +
  my_theme_div +
  guides(fill = "none")

Shannon_plot_HILIC


# Merged Figure -----------------------------------------------------------
Combined_diversity <- (Shannon_plot_HILIC | Shannon_plot_RP) +
  plot_annotation(tag_levels = "a") &          # lowercase panel tags
  theme(plot.tag = element_text(face = "bold", size = TAG_FONT_PT))

Combined_diversity

ggsave("Plots/New/SupFigure2_Diversity_HILIC_RP.png", Combined_diversity,
       dpi = 300, width = 174, height = 100, units = "mm", bg = "transparent")


Combined_diversity <- (Shannon_plot_HILIC | Shannon_plot_RP) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = TAG_FONT_PT)
  )

ggsave("Plots/New/Diversity_HILIC_RP.png", Combined_diversity,
       dpi = 300, width = 174, height = 100, units = "mm", bg = "transparent")


# 1. bind both tables and label method
Shannon_table_all <- bind_rows(
  Shannon_table_HILIC %>% mutate(Method = "HILIC"),
  Shannon_table_RP    %>% mutate(Method = "RP")
) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Method = factor(Method, levels = c("HILIC", "RP")),
    Type = factor(Type, levels = names(type_colors))
  )

y_range <- range(Shannon_table_all$values, na.rm = TRUE)
pad <- diff(y_range) * 0.25
y_limits <- c(y_range[1], y_range[2] + pad)

stat_table_all <- Shannon_table_all %>%
  filter(!is.na(Site)) %>%
  group_by(Method) %>%
  wilcox_test(values ~ Site, p.adjust.method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "Site", step.increase = 0.06)

combined_plot <- ggplot(Shannon_table_all, aes(x = Site, y = values, fill = Type)) +
  geom_boxplot(outlier.size = 1, width = 0.7) +
  facet_wrap(~ Method, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  scale_x_discrete(labels = site_labels) +
  stat_pvalue_manual(
    stat_table_all,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE,
    inherit.aes = FALSE,
    size = (BASE_FONT_PT - 1) / .pt
  ) +
  coord_cartesian(ylim = y_limits) +
  labs(title = "Metabolite Shannon diversity (HILIC vs RP)",
       y = "Shannon diversity",
       x = NULL) +
  my_theme_div +
  theme(legend.position = "bottom") 

combined_plot

ggsave("Plots/Diversity_Site_Combined.png", combined_plot,
       dpi = 300, width = 174, height = 110, units = "mm", bg = "transparent")