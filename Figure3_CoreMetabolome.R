#### Peatlands Project - 2023 Surface Samples 
### Figure 3. Core metabolome overview
## Core metabolome - VK - Mass Distribution
# Notes: "Core metabolome" = Ubiquitous (100% prevalence) + Near-universal (>=95% prevalence)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
  library(MASS)
  library(scales)
  library(ggrepel)
  library(ggforce)   # required for geom_arc_bar() in make_pie()
  library(magick)
})

source("functions_cdis_exploration_1.R")
separate_formula2 <- function(df, formula_col = "Formula", keep_formula = TRUE) {
  if (!formula_col %in% colnames(df)) stop("Formula column not found in df")
  df <- df %>% mutate(!!formula_col := as.character(.data[[formula_col]]))
  df[[formula_col]] <- str_remove_all(df[[formula_col]], "\\s+")
  
  all_elems <- df %>%
    pull(!!sym(formula_col)) %>%
    map(~ str_extract_all(.x, "[A-Z][a-z]?")[[1]]) %>%
    unlist() %>%
    unique() %>%
    sort()
  
  if (length(all_elems) == 0) return(df)
  
  counts_df <- map_dfc(all_elems, function(el) {
    pattern_digits <- paste0("(?<=", el, ")\\d+")
    present_pattern <- paste0("(?<![A-Za-z])", el, "(?![a-z])|(?<![A-Za-z])", el, "(?=\\d)")
    
    formulas <- df[[formula_col]]
    digits <- str_extract(formulas, pattern_digits)
    present <- str_detect(formulas, present_pattern)
    
    counts <- as.integer(digits)
    counts[is.na(counts) & present] <- 1L
    counts[is.na(counts) & !present] <- 0L
    tibble(!!el := counts)
  })
  
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

## size 

BASE_FONT_PT  <- 9    # axis text, axis titles, legend text, pie/mass labels
TITLE_FONT_PT <- 10   # plot titles
TAG_FONT_PT   <- 12   # bold panel tags (a, b, c...)

pt_to_mm <- function(pt) pt / .pt   # ggplot2 geom_text/annotate 'size' is in mm

# Colors
pie_colors <- c(
  "Remaining features" = "#7A7A7A",
  "95% prevalence"     = "#2F8B8B",
  "100% prevalence"    = "#8B4513"
)

mass_colors <- c(
  "Core metabolome"  = "#2F8B8B",
  "Distinctive (RF)" = "#7A7A7A"
)

figure_theme <- theme_bw(base_size = BASE_FONT_PT) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = TITLE_FONT_PT),
    axis.title = element_text(face = "bold", size = BASE_FONT_PT),
    axis.text = element_text(color = "black", size = BASE_FONT_PT),
    legend.title = element_text(face = "bold", size = BASE_FONT_PT),
    legend.text = element_text(size = BASE_FONT_PT),
    legend.position = "right",
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.margin = margin(4, 4, 4, 4)
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

# 2) PIE CHART DATA + PLOTS ==================================================
make_pie_plot <- function(feature_stats_df, method_label) {
  pie_df <- feature_stats_df %>%
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
      ymax     = cumsum(Fraction),
      ymin     = lag(ymax, default = 0),
      mid      = (ymin + ymax) / 2,
      Label    = paste0(round(Fraction * 100, 1), "%"),
      x_anchor    = c(2.00, 2.12, 2.24),
      y_anchor    = mid,
      x_label     = c(1.75, 2.55, 2.55),
      y_label     = c(mid[1], 0.89, 1.00),
      x_line_end  = x_label - 0.15
    )
  
  ggplot(pie_df) +
    geom_rect(
      aes(
        xmin = c(1.00, 1.12, 1.24),
        xmax = c(2.00, 2.12, 2.24),
        ymin = ymin,
        ymax = ymax,
        fill = DonutCategory
      ),
      color = "white", linewidth = 0.5
    ) +
    geom_segment(
      aes(x = x_anchor, xend = x_line_end,
          y = y_anchor, yend = y_label),
      color = "black", linewidth = 0.35
    ) +
    geom_text(
      aes(x = x_label, y = y_label, label = Label),
      size = 3, fontface = "bold",
      hjust = c(1, 0.5, 0.5)
    ) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.9) +
    scale_fill_manual(
      values = pie_colors,
      name   = NULL,
      drop   = FALSE
    ) +
    labs(title = method_label) +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "right",
      plot.margin     = margin(1, 1, 1, 1),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

pie_hilic <- make_pie_plot(feature_stats_hilic, "HILIC")
pie_hilic
pie_rp    <- make_pie_plot(feature_stats_rp,    "RP")
pie_rp

pie_hilic | pie_rp


# 3) VAN KREVELEN PLOTS (no insets) ===========================================

classification <- tribble(
  ~Class,          ~OC_low, ~OC_high, ~HC_low, ~HC_high,
  "Lipid",         0.000,   0.300,    1.5,     2.5,
  "Unsat. HC",     0.000,   0.125,    0.8,     1.5,
  "Cond. HC",      0.000,   0.950,    0.2,     0.8,
  "Protein",       0.300,   0.550,    1.5,     2.3,
  "Amino sugar",   0.550,   0.700,    1.5,     2.2,
  "Carbohydrate",  0.700,   1.500,    1.5,     2.5,
  "Lignin",        0.125,   0.650,    0.8,     1.5,
  "Tannin",        0.650,   1.100,    0.8,     1.5
) %>%
  mutate(label_x = (OC_low + OC_high) / 2,
         label_y = HC_high - 0.10)

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

hilic_kde <- kde2d(VK_df_HILIC$O_to_C, VK_df_HILIC$H_to_C,
                   n = 200, lims = c(0, 1.5, 0, 2.5))
rp_kde    <- kde2d(VK_df_RP$O_to_C,    VK_df_RP$H_to_C,
                   n = 200, lims = c(0, 1.5, 0, 2.5))
shared_density_max <- max(c(hilic_kde$z, rp_kde$z), na.rm = TRUE)

make_vk_plot <- function(vk_df, method_label) {
  ggplot(vk_df, aes(x = O_to_C, y = H_to_C)) +
    stat_density_2d(
      aes(fill = after_stat(density)),
      geom = "raster", contour = FALSE, n = 200
    ) +
    geom_density_2d(color = "black", linewidth = 0.35, bins = 6) +
    geom_rect(
      data = classification,
      aes(xmin = OC_low, xmax = OC_high, ymin = HC_low, ymax = HC_high),
      inherit.aes = FALSE,
      fill = NA, color = "grey25", linetype = "dashed", linewidth = 0.4
    ) +
    geom_text(
      data = classification,
      aes(x = label_x, y = label_y, label = Class),
      inherit.aes = FALSE, size = pt_to_mm(BASE_FONT_PT - 1), color = "grey20"
    ) +
    scale_fill_gradient(
      low = "#F7FBFF", high = "#2171B5",
      limits = c(0, shared_density_max),
      name = "Density"
    ) +
    labs(title = method_label, x = "O:C", y = "H:C") +
    figure_theme +
    coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 2.5))
}

vk_plot_hilic <- make_vk_plot(VK_df_HILIC, "HILIC")
vk_plot_rp    <- make_vk_plot(VK_df_RP,    "RP")

vk_plot_hilic | vk_plot_rp

# 4) MASS DISTRIBUTION PLOTS ==================================================

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

shared_mass_limits <- range(c(mass_hilic_df$Mass, mass_rp_df$Mass), na.rm = TRUE)

make_mass_plot <- function(mass_df, method_label) {
  ggplot(mass_df, aes(x = Mass, fill = Category, color = Category)) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    geom_vline(
      data = mass_df %>%
        group_by(Category) %>%
        summarise(med = median(Mass), .groups = "drop"),
      aes(xintercept = med, color = Category),
      linetype = "dashed", linewidth = 0.7, show.legend = FALSE
    ) +
    scale_fill_manual(values  = mass_colors, drop = FALSE, name = NULL) +
    scale_color_manual(values = mass_colors, drop = FALSE, name = NULL) +
    labs(title = method_label, x = "m/z", y = "Density") +
    figure_theme +
    coord_cartesian(xlim = shared_mass_limits)
}

mass_plot_hilic <- make_mass_plot(mass_hilic_df, "HILIC")
mass_plot_rp    <- make_mass_plot(mass_rp_df,    "RP")

mass_plot_hilic | mass_plot_rp

# 5) ASSEMBLE: pie → VK → mass ================================================

figure3_final <- (
  (pie_hilic | pie_rp) /
    (vk_plot_hilic | vk_plot_rp) /
    (mass_plot_hilic | mass_plot_rp)
) +
  plot_layout(
    heights = c(2, 1.1, 0.85),
    guides  = "collect"
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position  = "right",
    legend.box       = "vertical",
    plot.tag         = element_text(face = "bold", size = TAG_FONT_PT)
  )

figure3_final

ggsave(
  "Plots/Figure3_core_metabolome_final_upated.svg",
  figure3_final,
  width  = 174,
  height = 200,
  units  = "mm",
  bg     = "white"
)


row3 <- (mass_plot_hilic | mass_plot_rp) +
  plot_layout(guides  = "collect") &
  theme(legend.position  = "bottom")

ggsave(
  "Plots/New/Figure3_row3.svg",
  row3,
  width  = 174,
  height = 90,
  units  = "mm",
  bg     = "white"
)


row2 <- (vk_plot_hilic | vk_plot_rp) +
  plot_layout(guides  = "collect") &
  theme(legend.position  = "bottom")

ggsave(
  "Plots/New/Figure3_row2.svg",
  row2,
  width  = 174,
  height = 110,
  units  = "mm",
  bg     = "white"
)


row1 <- (pie_hilic | pie_rp) +
  plot_layout(guides  = "collect")

ggsave(
  "Plots/New/Figure3_row1.png",
  row1,
  width  = 174,
  height = 90,
  units  = "mm",
  bg     = "white"
)

