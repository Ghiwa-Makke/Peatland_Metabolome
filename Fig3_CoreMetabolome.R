#### Peatlands Project - 2023 Surface Samples 
### Figure 3. Core metabolome overview
## Core metabolome - VK - Mass Distribution
# Notes: "Core metabolome" = Universal (100% prevalence) + Near-universal (>=95% prevalence)

### Setting Libraries - bakcgrounds 
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
  library(MASS)
  library(scales)
})

source("functions_cdis_exploration_1.R")
separate_formula2 <- function(df, formula_col = "Formula", keep_formula = TRUE) {
  # Validate
  if (!formula_col %in% colnames(df)) stop("Formula column not found in df")
  df <- df %>% mutate(!!formula_col := as.character(.data[[formula_col]]))
  
  # remove whitespace inside formulas (just in case)
  df[[formula_col]] <- str_remove_all(df[[formula_col]], "\\s+")
  
  # 1) extract all element symbols across all formulas: uppercase letter + optional lowercase
  all_elems <- df %>%
    pull(!!sym(formula_col)) %>%
    # extract all element tokens per formula and flatten
    map(~ str_extract_all(.x, "[A-Z][a-z]?")[[1]]) %>%
    unlist() %>%
    unique() %>%
    sort()
  
  # if no elements found, return original
  if (length(all_elems) == 0) return(df)
  
  # 2) for each element, compute counts (use map_dfc to bind columns)
  counts_df <- map_dfc(all_elems, function(el) {
    # regex to capture digits immediately after the element symbol
    # (?<=El)\d+  matches digits following element
    pattern_digits <- paste0("(?<=", el, ")\\d+")
    # detect presence of element (either followed by digits or not)
    # use lookarounds to ensure we match the element symbol as an element token
    present_pattern <- paste0("(?<![A-Za-z])", el, "(?![a-z])|(?<![A-Za-z])", el, "(?=\\d)")
    
    formulas <- df[[formula_col]]
    digits <- str_extract(formulas, pattern_digits)          # NA if no digits
    present <- str_detect(formulas, present_pattern)         # TRUE if element symbol appears
    
    counts <- as.integer(digits)
    counts[is.na(counts) & present] <- 1L
    counts[is.na(counts) & !present] <- 0L
    tibble(!!el := counts)
  })
  
  # 3) combine
  if (keep_formula) {
    out <- bind_cols(df, counts_df)
  } else {
    out <- bind_cols(df %>% select(-all_of(formula_col)), counts_df)
  }
  
  out <- distinct(out)
  return(out)
}

prev_cut <- 0.95
cv_cut   <- 10
presence_thresh <- 0

# Colors
pie_colors <- c(
  "Remaining features" = "#E0E0E0",
  "95% prevalence"     = "#2F8B8B",
  "100% prevalence"    = "#8B4513"
)

mass_colors <- c(
  "Core metabolome"  = "#2F8B8B",
  "Distinctive (RF)" = "#7A7A7A"
)

# Shared figure theme
figure_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.margin = margin(5, 5, 5, 5)
  )


# 1) CORE CATEGORY CALCULATIONS -----
##  HILIC ------
HILIC_raw <- read.csv(
  "Tables/compounds_table_non_real_gap_filled_HILIC.csv",
  check.names = FALSE
)

HILIC_matrix <- read.csv(
  "Tables/HILIC_norm_mean.csv",
  row.names = 1,
  check.names = FALSE
)

HILIC_matrix_long <- HILIC_matrix %>%
  rownames_to_column("FeatureID") %>%
  pivot_longer(
    cols = -FeatureID,
    names_to = "SampleID",
    values_to = "Intensity_norm"
  )

HILIC_with_norm <- HILIC_raw %>%
  left_join(HILIC_matrix_long, by = c("FeatureID", "SampleID")) %>%
  filter(!is.na(Date))

n_samples_hilic <- HILIC_with_norm %>%
  distinct(SampleID) %>%
  nrow()

feature_prevalence_hilic <- HILIC_with_norm %>%
  group_by(FeatureID) %>%
  summarise(
    n_present = n_distinct(SampleID[Intensity_norm > presence_thresh]),
    prop_present = n_present / n_samples_hilic,
    .groups = "drop"
  )

feature_cv_hilic <- HILIC_with_norm %>%
  group_by(FeatureID) %>%
  summarise(
    mean_intensity = mean(Intensity_norm, na.rm = TRUE),
    sd_intensity   = sd(Intensity_norm, na.rm = TRUE),
    cv_percent     = (sd_intensity / mean_intensity) * 100,
    .groups = "drop"
  ) %>%
  mutate(cv_percent = ifelse(is.finite(cv_percent), cv_percent, NA_real_))

feature_stats_hilic <- feature_prevalence_hilic %>%
  left_join(feature_cv_hilic, by = "FeatureID") %>%
  mutate(
    DonutCategory = case_when(
      prop_present == 1 & cv_percent <= cv_cut ~ "Universal",
      prop_present >= prev_cut & prop_present < 1 & cv_percent <= cv_cut ~ "Near-universal",
      TRUE ~ "Remaining"
    ),
    CoreGroup = if_else(
      DonutCategory %in% c("Universal", "Near-universal"),
      "Core metabolome",
      "Remaining"
    )
  )

## RP -----------------------------
RP_raw <- read.csv(
  "Tables/compounds_table_non_real_gap_filled_RP.csv",
  check.names = FALSE
)

RP_matrix <- read.csv(
  "Tables/RP_norm_mean.csv",
  row.names = 1,
  check.names = FALSE
)

# Restore sample ID style if needed
colnames(RP_matrix) <- str_replace_all(colnames(RP_matrix), "\\.", "-")

RP_matrix_long <- RP_matrix %>%
  rownames_to_column("FeatureID") %>%
  pivot_longer(
    cols = -FeatureID,
    names_to = "SampleID",
    values_to = "Intensity_norm"
  )

RP_with_norm <- RP_raw %>%
  left_join(RP_matrix_long, by = c("FeatureID", "SampleID")) %>%
  filter(!is.na(Date))

n_samples_rp <- RP_with_norm %>%
  distinct(SampleID) %>%
  nrow()

feature_prevalence_rp <- RP_with_norm %>%
  group_by(FeatureID) %>%
  summarise(
    n_present = n_distinct(SampleID[Intensity_norm > presence_thresh]),
    prop_present = n_present / n_samples_rp,
    .groups = "drop"
  )

feature_cv_rp <- RP_with_norm %>%
  group_by(FeatureID) %>%
  summarise(
    mean_intensity = mean(Intensity_norm, na.rm = TRUE),
    sd_intensity   = sd(Intensity_norm, na.rm = TRUE),
    cv_percent     = (sd_intensity / mean_intensity) * 100,
    .groups = "drop"
  ) %>%
  mutate(cv_percent = ifelse(is.finite(cv_percent), cv_percent, NA_real_))

feature_stats_rp <- feature_prevalence_rp %>%
  left_join(feature_cv_rp, by = "FeatureID") %>%
  mutate(
    DonutCategory = case_when(
      prop_present == 1 & cv_percent <= cv_cut ~ "Universal",
      prop_present >= prev_cut & prop_present < 1 & cv_percent <= cv_cut ~ "Near-universal",
      TRUE ~ "Remaining"
    ),
    CoreGroup = if_else(
      DonutCategory %in% c("Universal", "Near-universal"),
      "Core metabolome",
      "Remaining"
    )
  )

# 2) Pie Data -----
# HILIC pie data
library(ggrepel)

hilic_pie_df <- feature_stats_hilic %>%
  count(DonutCategory, name = "Count") %>%
  mutate(
    DonutCategory = factor(
      DonutCategory,
      levels = c("Remaining", "Near-universal", "Universal"),
      labels = c("Remaining features", "95% prevalence", "100% prevalence")
    )
  ) %>%
  arrange(DonutCategory) %>%
  mutate(
    Fraction = Count / sum(Count),
    ymax = cumsum(Fraction),
    ymin = lag(ymax, default = 0),
    mid = (ymin + ymax) / 2,
    Label = paste0(round(Fraction * 100, 1), "%"),
    
    # anchor point on slice
    x_anchor = c(2.00, 2.12, 2.24),
    y_anchor = mid,
    
    # label positions:
    # keep remaining where it is, move the two small slices to the top
    x_label = c(1.75, 2.55, 2.55),
    y_label = c(mid[1],0.89, 1),
    # line end BEFORE label (buffer space)
    x_line_end = x_label - 0.15
  )

# Inset exploded pie: HILIC
pie_hilic <- ggplot(hilic_pie_df) +
  geom_rect(
    aes(
      xmin = c(1.00, 1.12, 1.24),
      xmax = c(2.00, 2.12, 2.24),
      ymin = ymin,
      ymax = ymax,
      fill = DonutCategory
    ),
    color = "white",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      x = x_anchor,
      xend = x_line_end,
      y = y_anchor,
      yend = y_label
    ),
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      x = x_label,
      y = y_label,
      label = Label
    ),
    size = 2.5,
    fontface = "bold",
    hjust = c(1, 0.5, 0.5)
  ) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.7) +
  scale_fill_manual(values = pie_colors, drop = FALSE) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "none",
    plot.background  = element_rect(fill = alpha("white", 0.1), color = NA),
    panel.background = element_rect(fill = alpha("white", 0.1), color = NA)
  )
pie_hilic
# Inset exploded pie: RP
# RP pie data
rp_pie_df <- feature_stats_rp %>%
  count(DonutCategory, name = "Count") %>%
  mutate(
    DonutCategory = factor(
      DonutCategory,
      levels = c("Remaining", "Near-universal", "Universal"),
      labels = c("Remaining features", "95% prevalence", "100% prevalence")
    )
  ) %>%
  arrange(DonutCategory) %>%
  mutate(
    Fraction = Count / sum(Count),
    ymax = cumsum(Fraction),
    ymin = lag(ymax, default = 0),
    mid = (ymin + ymax) / 2,
    Label = paste0(round(Fraction * 100, 1), "%"),
    
    # anchor point on slice
    x_anchor = c(2.00, 2.12, 2.24),
    y_anchor = mid,
    
    # label positions:
    # keep remaining where it is, move the two small slices to the top
    x_label = c(1.75, 2.55, 2.55),
    y_label = c(mid[1],0.87, 0.99),
    # line end BEFORE label (buffer space)
    x_line_end = x_label - 0.15
  )
pie_rp <- ggplot(rp_pie_df) +
  geom_rect(
    aes(
      xmin = c(1.00, 1.12, 1.24),
      xmax = c(2.00, 2.12, 2.24),
      ymin = ymin,
      ymax = ymax,
      fill = DonutCategory
    ),
    color = "white",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      x = x_anchor,
      xend = x_line_end,
      y = y_anchor,
      yend = y_label
    ),
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      x = x_label,
      y = y_label,
      label = Label
    ),
    size = 2.5,
    fontface = "bold",
    hjust = c(1, 0.5, 0.5)
  ) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.7) +
  scale_fill_manual(values = pie_colors, drop = FALSE) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "none",
    plot.background  = element_rect(fill = alpha("white", 0.1), color = NA),
    panel.background = element_rect(fill = alpha("white", 0.1), color = NA)
  )
pie_rp
# Legend-only plot for pie categories
pie_legend_dummy <- ggplot(
  tibble(
    x = c(1, 2, 3),
    y = c(1, 1, 1),
    DonutCategory = factor(
      c("Remaining features", "95% prevalence", "100% prevalence"),
      levels = c("Remaining features", "95% prevalence", "100% prevalence")
    )
  ),
  aes(x = x, y = y, fill = DonutCategory)
) +
  geom_col(width = 0.8, alpha = 0) +   # invisible plot, legend still collected
  scale_fill_manual(
    values = pie_colors,
    name = NULL
  ) +
  theme_void() +
  theme(legend.position = "right")


# 3) VAN KREVELEN DATA FOR CORE METABOLOME ============================================

# VK classification boxes
classification <- tribble(
  ~Class, ~OC_low, ~OC_high, ~HC_low, ~HC_high,
  "Lipid",        0.000, 0.300, 1.5, 2.5,
  "Unsat. HC",    0.000, 0.125, 0.8, 1.5,
  "Cond. HC",     0.000, 0.950, 0.2, 0.8,
  "Protein",      0.300, 0.550, 1.5, 2.3,
  "Amino sugar",  0.550, 0.700, 1.5, 2.2,
  "Carbohydrate", 0.700, 1.500, 1.5, 2.5,
  "Lignin",       0.125, 0.650, 0.8, 1.5,
  "Tannin",       0.650, 1.100, 0.8, 1.5
) %>%
  mutate(
    label_x = (OC_low + OC_high) / 2,
    label_y = HC_high - 0.10
  )


## HILIC core-metabolome VK --------------------------
HILIC_core_ids <- read.csv("Tables/Near_Univeral_Core_HILIC.csv") %>%
  mutate(CoreGroup = "Core metabolome") %>%
  dplyr::select(FeatureID, CoreGroup)

HILIC_vk <- read_xlsx("Tables/HILIC_annot_features_final.xlsx") %>%
  dplyr::select(FeatureID, Annotation_Confidence, Final_formula) %>%
  left_join(HILIC_core_ids, by = "FeatureID") %>%
  filter(!is.na(CoreGroup))

HILIC_vk <- separate_formula2(HILIC_vk, formula_col = "Final_formula", keep_formula = TRUE)
HILIC_vk <- calc_ratios_n_idxs(HILIC_vk)
HILIC_vk <- calc_classes(HILIC_vk)

VK_df_HILIC <- HILIC_vk %>%
  transmute(
    FeatureID,
    O_to_C = as.numeric(O_to_C),
    H_to_C = as.numeric(H_to_C),
    CoreGroup
  ) %>%
  filter(is.finite(O_to_C), is.finite(H_to_C))

## RP core-metabolome VK ----------------------------
RP_core_ids <- read.csv("Tables/Near_Univeral_Core_RP.csv") %>%
  mutate(CoreGroup = "Core metabolome") %>%
  dplyr::select(FeatureID, CoreGroup)

RP_vk <- read_xlsx("Tables/RP_annot_features_final.xlsx") %>%
  dplyr::select(FeatureID, Annotation_Confidence, Final_formula) %>%
  left_join(RP_core_ids, by = "FeatureID") %>%
  filter(!is.na(CoreGroup))

RP_vk <- separate_formula2(RP_vk, formula_col = "Final_formula", keep_formula = TRUE)

# Safeguard if these are absent in RP formulas
if (!"P" %in% names(RP_vk)) RP_vk$P <- 0
if (!"S" %in% names(RP_vk)) RP_vk$S <- 0

RP_vk <- calc_ratios_n_idxs(RP_vk)
RP_vk <- calc_classes(RP_vk)

VK_df_RP <- RP_vk %>%
  transmute(
    FeatureID,
    O_to_C = as.numeric(O_to_C),
    H_to_C = as.numeric(H_to_C),
    CoreGroup
  ) %>%
  filter(is.finite(O_to_C), is.finite(H_to_C))


# Shared continuous density scale across HILIC and RP

hilic_kde <- kde2d(
  x = VK_df_HILIC$O_to_C,
  y = VK_df_HILIC$H_to_C,
  n = 200,
  lims = c(0, 1.5, 0, 2.5)
)

rp_kde <- kde2d(
  x = VK_df_RP$O_to_C,
  y = VK_df_RP$H_to_C,
  n = 200,
  lims = c(0, 1.5, 0, 2.5)
)

shared_density_max <- max(c(hilic_kde$z, rp_kde$z), na.rm = TRUE)

# Base HILIC VK plot
vk_plot_hilic_base <- ggplot(VK_df_HILIC, aes(x = O_to_C, y = H_to_C)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE,
    n = 200
  ) +
  geom_density_2d(
    color = "black",
    linewidth = 0.35,
    bins = 6
  ) +
  geom_rect(
    data = classification,
    aes(xmin = OC_low, xmax = OC_high, ymin = HC_low, ymax = HC_high),
    inherit.aes = FALSE,
    fill = NA,
    color = "grey25",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_text(
    data = classification,
    aes(x = label_x, y = label_y, label = Class),
    inherit.aes = FALSE,
    size = 2.8,
    color = "grey20"
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2171B5",
    limits = c(0, shared_density_max),
    name = "Density"
  ) +
  labs(
    title = "HILIC",
    x = "O:C",
    y = "H:C"
  ) +
  figure_theme +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 2.5))

# Base RP VK plot
vk_plot_rp_base <- ggplot(VK_df_RP, aes(x = O_to_C, y = H_to_C)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE,
    n = 200
  ) +
  geom_density_2d(
    color = "black",
    linewidth = 0.35,
    bins = 6
  ) +
  geom_rect(
    data = classification,
    aes(xmin = OC_low, xmax = OC_high, ymin = HC_low, ymax = HC_high),
    inherit.aes = FALSE,
    fill = NA,
    color = "grey25",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_text(
    data = classification,
    aes(x = label_x, y = label_y, label = Class),
    inherit.aes = FALSE,
    size = 2.8,
    color = "grey20"
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2171B5",
    limits = c(0, shared_density_max),
    name = "Density"
  ) +
  labs(
    title = "RP",
    x = "O:C",
    y = "H:C"
  ) +
  figure_theme +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 2.5))

### Add inset exploded pies into bottom-right corner -----
vk_plot_hilic <- vk_plot_hilic_base +
  inset_element(pie_hilic, left = 0.62, bottom = 0.01, right = 0.99, top = 0.54)

vk_plot_rp <- vk_plot_rp_base +
  inset_element(pie_rp, left = 0.62, bottom = 0.01, right = 0.99, top = 0.54)


# 4) MASS DISTRIBUTION DATA ================================================
## HILIC mass: Core metabolome vs Distinctive (RF plateau set) -----

HILIC_Core_mass <- read.csv("Tables/Near_Univeral_Core_HILIC.csv") %>%
  mutate(
    Category = "Core metabolome",
    Mass = as.numeric(str_extract(FeatureID, "(?<!\\d)(?:[5-9]\\d|\\d{3,4})\\.\\d+"))
  ) %>%
  dplyr::select(FeatureID, Mass, Category)

HILIC_RF_mass <- read.csv("Tables/RF_selected_features_plateau_HILIC.csv") %>%
  rename(FeatureID = feature) %>%
  mutate(
    Mass = as.numeric(str_extract(FeatureID, "(?<!\\d)(?:[5-9]\\d|\\d{3,4})\\.\\d+")),
    Category = "Distinctive (RF)"
  ) %>%
  dplyr::select(FeatureID, Mass, Category)

mass_hilic_df <- bind_rows(HILIC_Core_mass, HILIC_RF_mass) %>%
  filter(!is.na(Mass)) %>%
  mutate(
    Category = factor(Category, levels = c("Core metabolome", "Distinctive (RF)"))
  )

## RP mass: Core metabolome vs top 30 RF features -----
RP_Core_mass <- read.csv("Tables/Near_Univeral_Core_RP.csv") %>%
  mutate(
    Category = "Core metabolome",
    Mass = as.numeric(str_extract(FeatureID, "(?<!\\d)(?:[5-9]\\d|\\d{3,4})\\.\\d+"))
  ) %>%
  dplyr::select(FeatureID, Mass, Category)

RP_RF_mass <- read.csv("Tables/Site_Features_RF_importance_fullmodel_RP.csv") %>%
  mutate(feature = gsub("`", "", feature)) %>%
  rename(FeatureID = feature) %>%
  arrange(desc(importance)) %>%
  slice(1:30) %>%
  mutate(
    Mass = as.numeric(str_extract(FeatureID, "(?<!\\d)(?:[5-9]\\d|\\d{3,4})\\.\\d+")),
    Category = "Distinctive (RF)"
  ) %>%
  dplyr::select(FeatureID, Mass, Category)

mass_rp_df <- bind_rows(RP_Core_mass, RP_RF_mass) %>%
  filter(!is.na(Mass)) %>%
  mutate(
    Category = factor(Category, levels = c("Core metabolome", "Distinctive (RF)"))
  )

# Shared x-axis range for comparison
shared_mass_limits <- range(c(mass_hilic_df$Mass, mass_rp_df$Mass), na.rm = TRUE)

# HILIC mass plot
mass_plot_hilic <- ggplot(mass_hilic_df, aes(x = Mass, fill = Category, color = Category)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(
    data = mass_hilic_df %>%
      group_by(Category) %>%
      summarise(med = median(Mass), .groups = "drop"),
    aes(xintercept = med, color = Category),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = mass_colors, drop = FALSE, name = NULL) +
  scale_color_manual(values = mass_colors, drop = FALSE, name = NULL) +
  labs(
    title = "HILIC",
    x = "m/z",
    y = "Density"
  ) +
  figure_theme +
  coord_cartesian(xlim = shared_mass_limits)

# RP mass plot
mass_plot_rp <- ggplot(mass_rp_df, aes(x = Mass, fill = Category, color = Category)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(
    data = mass_rp_df %>%
      group_by(Category) %>%
      summarise(med = median(Mass), .groups = "drop"),
    aes(xintercept = med, color = Category),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = mass_colors, drop = FALSE, name = NULL) +
  scale_color_manual(values = mass_colors, drop = FALSE, name = NULL) +
  labs(
    title = "RP",
    x = "m/z",
    y = "Density"
  ) +
  figure_theme +
  coord_cartesian(xlim = shared_mass_limits)


# 5) FINAL PUBLICATION-READY MULTIPANEL FIGURE ==============================================

# Main 2x2 panel figure
figure3_final <- (
  (vk_plot_hilic   | vk_plot_rp) /
    (mass_plot_hilic | mass_plot_rp)
) +
  plot_layout(
    heights = c(1.15, 0.85, 0.001),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.tag = element_text(face = "bold", size = 13)
  )

# Display
figure3_final

ggsave(
  filename = "Plots/Figure3_core_metabolome_final.png",
  plot = figure3_final,
  width = 9.2,
  height = 6,
  units = "in",
  bg = "white"
)



