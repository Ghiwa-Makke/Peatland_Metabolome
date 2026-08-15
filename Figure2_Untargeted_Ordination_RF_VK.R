#### Peatlands Project - 2023 Surface Samples 
### Figure 2 - Untargeted Metabolomics - overview
## Ordination - RF - VK of RF features 

suppressPackageStartupMessages({
  library(tidyverse)
  library(tibble)
  library(readr)
  library(readxl)
  library(factoextra)
  library(caret)
  library(ranger)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Setting theme: ------
set.seed(1)

# font size 
BASE_FONT_PT  <- 9    # axis text, axis titles, legend text
TITLE_FONT_PT <- 10   # plot titles
TAG_FONT_PT   <- 12   # bold panel tags (a, b, c...)

pt_to_mm <- function(pt) pt / .pt   # ggplot2 geom_text/annotate 'size' is in mm

type_colors <- c(
  "Bog" = "#8B4513",
  "Shrub Wetland" = "#2F8B8B",
  "Fen" = "#F5DEB3"
)

# shapes ----
site_shapes <- c(
  "S1"  = 21,   
  "BLF" = 24,
  "S3"  = 22,
  "US2" = 23,
  "LC"  = 25
)

# ---- legend display labels (full descriptive site names) ----
site_labels <- c(
  "S1"  = "S1-bog",
  "BLF" = "BLF-poor fen",
  "S3"  = "S3-fen",
  "US2" = "US2-fen",
  "LC"  = "LC-shrub wetland"
)

figure_theme <- theme_bw(base_size = BASE_FONT_PT) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = TITLE_FONT_PT),
    axis.title = element_text(face = "bold", size = BASE_FONT_PT),
    axis.text = element_text(color = "black", size = BASE_FONT_PT),
    legend.title = element_text(face = "bold", size = BASE_FONT_PT),
    legend.text = element_text(size = BASE_FONT_PT),
    legend.position = "bottom",
    legend.box = "vertical"
  )


# HILIC PCA ------
hilic_feat <- read_csv("Tables/HILIC_norm_mean.csv") %>%
  column_to_rownames(var = "...1")
hilic_meta <- read_csv("Tables/fixed_metadata_hilic.csv")
hilic_meta$SampleID <- gsub("-", "_", hilic_meta$SampleID)

hilic_pca <- prcomp(t(hilic_feat), scale. = FALSE)
hilic_eig <- get_eigenvalue(hilic_pca)

hilic_pca_df <- as_tibble(hilic_pca$x, rownames = "SampleID") %>%
  left_join(hilic_meta, by = "SampleID") %>%
  mutate(
    Type = factor(Type, levels = c("Bog", "Fen", "Shrub Wetland")),
    Site = factor(Site, levels = c("S1", "BLF", "S3", "US2", "LC"))
  )
hilic_pc1 <- paste0("PC1 (", round(hilic_eig$variance.percent[1], 1), "%)")
hilic_pc2 <- paste0("PC2 (", round(hilic_eig$variance.percent[2], 1), "%)")

pca_plot_hilic <- ggplot(hilic_pca_df, aes(PC1, PC2)) +
  stat_ellipse(aes(color = Type), linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  geom_point(
    aes(fill = Type, shape = Site),
    color = "black",       # fixed black outline
    size = 2.8, stroke = 0.4
  ) +
  scale_color_manual(values = type_colors, drop = FALSE) +
  scale_fill_manual(values = type_colors, drop = FALSE) +
  scale_shape_manual(values = site_shapes, labels = site_labels, drop = FALSE) +
  labs(
    title = "HILIC mode \u2013 Polar compounds",
    x = hilic_pc1,
    y = hilic_pc2,
    fill = "Peatland type",
    shape = "Site"
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  figure_theme 

pca_plot_hilic


# RP PCA ---- 
rp_feat <- read_csv("Tables/RP_norm_mean.csv") %>%
  column_to_rownames(var = "...1")

rp_meta <- read_csv("Tables/fixed_metadata_rp.csv")
rp_meta$SampleID <- gsub("-", "_", rp_meta$SampleID)

rp_pca <- prcomp(t(rp_feat), scale. = FALSE)
rp_eig <- get_eigenvalue(rp_pca)

rp_pca_df <- as_tibble(rp_pca$x, rownames = "SampleID") %>%
  left_join(rp_meta, by = "SampleID") %>%
  mutate(
    Type = factor(Type, levels = c("Bog", "Fen", "Shrub Wetland")),
    Site = factor(Site, levels = c("S1", "BLF", "S3", "US2", "LC"))
  )

rp_pc1 <- paste0("PC1 (", round(rp_eig$variance.percent[1], 1), "%)")
rp_pc2 <- paste0("PC2 (", round(rp_eig$variance.percent[2], 1), "%)")

pca_plot_rp <- ggplot(rp_pca_df, aes(PC1, PC2)) +
  stat_ellipse(aes(color = Type), linewidth = 0.6, alpha = 0.9, show.legend = FALSE) +
  geom_point(
    aes(fill = Type, shape = Site),
    color = "black",
    size = 2.8, stroke = 0.4
  ) +
  scale_color_manual(values = type_colors, drop = FALSE) +
  scale_fill_manual(values = type_colors, drop = FALSE) +
  scale_shape_manual(values = site_shapes, labels = site_labels, drop = FALSE) +
  labs(
    title = "RP mode \u2013 Semi-polar compounds",
    x = rp_pc1,
    y = rp_pc2,
    fill = "Peatland type",
    shape = "Site"
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  figure_theme 

pca_plot_rp

# HILIC RANDOM FOREST + CONFUSION MATRIX -----
hilic_X <- t(hilic_feat) %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(hilic_meta %>% select(SampleID, Site), by = "SampleID") %>%
  filter(!is.na(Site))

hilic_X$Site <- factor(hilic_X$Site, levels = c("S1", "BLF", "S3", "US2", "LC"))
hilic_X[is.na(hilic_X)] <- 0
hilic_predictors <- hilic_X %>% select(-SampleID, -Site)
hilic_nzv <- nearZeroVar(hilic_predictors)
if (length(hilic_nzv) > 0) {
  hilic_predictors <- hilic_predictors[, -hilic_nzv, drop = FALSE]
}
hilic_rf_df <- bind_cols(Site = hilic_X$Site, as.data.frame(hilic_predictors))

ctrl_hilic <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = multiClassSummary
)

rf_hilic <- train(
  Site ~ .,
  data = hilic_rf_df,
  method = "ranger",
  trControl = ctrl_hilic,
  metric = "Accuracy",
  importance = "permutation",
  num.trees = 2000
)

hilic_imp <- varImp(rf_hilic, scale = FALSE)$importance %>%
  rownames_to_column("feature") %>%
  rename(importance = Overall) %>%
  mutate(feature = gsub("`", "", feature)) %>%
  arrange(desc(importance))

write_csv(hilic_imp, "Tables/Site_Features_RF_importance_fullmodel_HILIC.csv")

hilic_pred <- rf_hilic$pred
if (!is.null(rf_hilic$bestTune) && ncol(rf_hilic$bestTune) > 0) {
  for (nm in names(rf_hilic$bestTune)) {
    hilic_pred <- hilic_pred[hilic_pred[[nm]] == rf_hilic$bestTune[[nm]], , drop = FALSE]
  }
}
cm_hilic <- confusionMatrix(data = hilic_pred$pred, reference = hilic_pred$obs)
cm_hilic_df <- as.data.frame(cm_hilic$table) %>%
  mutate(Prediction = factor(Prediction, levels = rev(levels(Prediction))))

cm_hilic_df <- cm_hilic_df %>%
  group_by(Reference) %>%
  mutate(
    Percent = 100 * Freq / sum(Freq),
    Label = sprintf("%.1f%%", Percent)
  ) %>%
  ungroup()

hilic_acc <- round(cm_hilic$overall["Accuracy"] * 100, 1)

cm_plot_hilic <- ggplot(cm_hilic_df, aes(x = Reference, y = Prediction, fill = Percent)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = pt_to_mm(BASE_FONT_PT), color = "black", fontface = "bold") +
  scale_fill_gradient(
    low = "#F5F5F5",
    high = "#8EC5E8",
    limits = c(0, 100),
    name = "Percent"
  ) +
  #scale_x_discrete(labels = site_labels) +
  #scale_y_discrete(labels = site_labels) +
  labs(
    title = paste0("HILIC Random Forest\nOverall CV accuracy = ", hilic_acc, "%"),
    x = "Observed site",
    y = "Predicted site"
  ) +
  figure_theme +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

cm_plot_hilic


# RP RANDOM FOREST + CONFUSION MATRIX ----
rp_X <- t(rp_feat) %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(rp_meta %>% select(SampleID, Site), by = "SampleID") %>%
  filter(!is.na(Site))

rp_X$Site <- factor(rp_X$Site, levels = c("S1", "BLF", "S3", "US2", "LC"))
rp_X[is.na(rp_X)] <- 0

rp_predictors <- rp_X %>% select(-SampleID, -Site)
rp_nzv <- nearZeroVar(rp_predictors)
if (length(rp_nzv) > 0) {
  rp_predictors <- rp_predictors[, -rp_nzv, drop = FALSE]
}
rp_rf_df <- bind_cols(Site = rp_X$Site, as.data.frame(rp_predictors))

ctrl_rp <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = multiClassSummary
)

rf_rp <- train(
  Site ~ .,
  data = rp_rf_df,
  method = "ranger",
  trControl = ctrl_rp,
  metric = "Accuracy",
  importance = "permutation",
  num.trees = 2000
)

rp_imp <- varImp(rf_rp, scale = FALSE)$importance %>%
  rownames_to_column("feature") %>%
  rename(importance = Overall) %>%
  mutate(feature = gsub("`", "", feature)) %>%
  arrange(desc(importance))

write_csv(rp_imp, "Tables/Site_Features_RF_importance_fullmodel_RP.csv")

rp_pred <- rf_rp$pred
if (!is.null(rf_rp$bestTune) && ncol(rf_rp$bestTune) > 0) {
  for (nm in names(rf_rp$bestTune)) {
    rp_pred <- rp_pred[rp_pred[[nm]] == rf_rp$bestTune[[nm]], , drop = FALSE]
  }
}
cm_rp <- confusionMatrix(data = rp_pred$pred, reference = rp_pred$obs)
cm_rp_df <- as.data.frame(cm_rp$table) %>%
  mutate(Prediction = factor(Prediction, levels = rev(levels(Prediction))))

cm_rp_df <- cm_rp_df %>%
  group_by(Reference) %>%
  mutate(
    Percent = 100 * Freq / sum(Freq),
    Label = sprintf("%.1f%%", Percent)
  ) %>%
  ungroup()

rp_acc <- round(cm_rp$overall["Accuracy"] * 100, 1)

cm_plot_rp <- ggplot(cm_rp_df, aes(x = Reference, y = Prediction, fill = Percent)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = pt_to_mm(BASE_FONT_PT), color = "black", fontface = "bold") +
  scale_fill_gradient(
    low = "#F5F5F5",
    high = "#8EC5E8",
    limits = c(0, 100),
    name = "Percent"
  ) +
  #scale_x_discrete(labels = site_labels) +
  #scale_y_discrete(labels = site_labels) +
  labs(
    title = paste0("RP Random Forest\nOverall CV accuracy = ", rp_acc, "%"),
    x = "Observed site",
    y = "Predicted site"
  ) +
  figure_theme +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
cm_plot_rp


# SELECT PLATEAU FEATURES - HILIC ------
panel_sizes_hilic <- c(5, 10, 20, 30, 50, 75, 100, 150, 200, 300)
panel_sizes_hilic <- panel_sizes_hilic[panel_sizes_hilic <= nrow(hilic_imp)]

hilic_perf <- tibble(N = panel_sizes_hilic, Accuracy = NA_real_)

for (i in seq_along(panel_sizes_hilic)) {
  N <- panel_sizes_hilic[i]
  top_feats <- hilic_imp$feature[1:N]
  dfN <- hilic_rf_df[, c("Site", top_feats), drop = FALSE]
  set.seed(1)
  fitN <- train(
    Site ~ .,
    data = dfN,
    method = "ranger",
    trControl = ctrl_hilic,
    metric = "Accuracy",
    importance = "permutation",
    num.trees = 2000
  )
  hilic_perf$Accuracy[i] <- max(fitN$results$Accuracy)
}

hilic_max_acc <- max(hilic_perf$Accuracy, na.rm = TRUE)
hilic_target_acc <- hilic_max_acc - 0.01
hilic_N_plateau <- hilic_perf$N[min(which(hilic_perf$Accuracy >= hilic_target_acc))]
hilic_top_feats_plateau <- hilic_imp$feature[1:hilic_N_plateau]

write_csv(
  tibble(feature = hilic_top_feats_plateau),
  "Tables/RF_selected_features_plateau_HILIC.csv"
)

plateau_plot_hilic <- ggplot(hilic_perf, aes(x = N, y = Accuracy)) +
  geom_line(linewidth = 0.7, color = "#333333") +
  geom_point(size = 2.5, color = "#333333") +
  geom_hline(yintercept = hilic_max_acc, linetype = "dashed", color = "#888888") +
  geom_hline(yintercept = hilic_target_acc, linetype = "dotted", color = "#888888") +
  geom_vline(xintercept = hilic_N_plateau, linetype = "dashed", color = "#2F8B8B") +
  annotate(
    "label",
    x = hilic_N_plateau,
    y = min(hilic_perf$Accuracy) + 0.02,
    label = paste0("Chosen N = ", hilic_N_plateau),
    size = pt_to_mm(BASE_FONT_PT)
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "HILIC feature-panel plateau",
    x = "Number of top-ranked features",
    y = "Repeated CV accuracy"
  ) +
  figure_theme

plateau_plot_hilic

# SELECT PLATEAU FEATURES - RP -----

panel_sizes_rp <- c(5, 10, 20, 30, 50, 75, 100, 150, 200, 300)
panel_sizes_rp <- panel_sizes_rp[panel_sizes_rp <= nrow(rp_imp)]

rp_perf <- tibble(N = panel_sizes_rp, Accuracy = NA_real_)

for (i in seq_along(panel_sizes_rp)) {
  N <- panel_sizes_rp[i]
  top_feats <- rp_imp$feature[1:N]
  dfN <- rp_rf_df[, c("Site", top_feats), drop = FALSE]
  set.seed(1)
  fitN <- train(
    Site ~ .,
    data = dfN,
    method = "ranger",
    trControl = ctrl_rp,
    metric = "Accuracy",
    importance = "permutation",
    num.trees = 2000
  )
  rp_perf$Accuracy[i] <- max(fitN$results$Accuracy)
}

rp_max_acc <- max(rp_perf$Accuracy, na.rm = TRUE)
rp_target_acc <- rp_max_acc - 0.01
rp_N_plateau <- rp_perf$N[min(which(rp_perf$Accuracy >= rp_target_acc))]
rp_top_feats_plateau <- rp_imp$feature[1:rp_N_plateau]

write_csv(
  tibble(feature = rp_top_feats_plateau),
  "Tables/RF_selected_features_plateau_RP.csv"
)

plateau_plot_rp <- ggplot(rp_perf, aes(x = N, y = Accuracy)) +
  geom_line(linewidth = 0.7, color = "#333333") +
  geom_point(size = 2.5, color = "#333333") +
  geom_hline(yintercept = rp_max_acc, linetype = "dashed", color = "#888888") +
  geom_hline(yintercept = rp_target_acc, linetype = "dotted", color = "#888888") +
  geom_vline(xintercept = rp_N_plateau, linetype = "dashed", color = "#2F8B8B") +
  annotate(
    "label",
    x = rp_N_plateau,
    y = min(rp_perf$Accuracy) + 0.02,
    label = paste0("Chosen N = ", rp_N_plateau),
    size = pt_to_mm(BASE_FONT_PT)
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "RP feature-panel plateau",
    x = "Number of top-ranked features",
    y = "Repeated CV accuracy"
  ) +
  figure_theme
plateau_plot_rp

supp_plateau <- plateau_plot_hilic | plateau_plot_rp
supp_plateau 


# E-F VAN KREVELEN SETUP -----
classification_boxes <- tribble(
  ~Class,          ~OC_low, ~OC_high, ~HC_low, ~HC_high,
  "Lipid",          0.00,    0.30,     1.5,     2.5,
  "Unsat. HC",      0.00,    0.125,    0.8,     1.5,
  "Cond. HC",       0.00,    0.95,     0.2,     0.8,
  "Protein",        0.30,    0.55,     1.5,     2.3,
  "Amino sugar",    0.55,    0.70,     1.5,     2.2,
  "Carbohydrate",   0.70,    1.50,     1.5,     2.5,
  "Lignin",         0.125,   0.65,     0.8,     1.5,
  "Tannin",         0.65,    1.10,     0.8,     1.5
) %>%
  mutate(
    label_x = (OC_low + OC_high) / 2,
    label_y = HC_high - 0.08
  )

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

library(MASS)

# E) HILIC VAN KREVELEN -------
hilic_vk_raw <- read_xlsx("Tables/HILIC_annot_features_final.xlsx") %>%
  filter(FeatureID %in% hilic_top_feats_plateau) %>%
  dplyr::select(FeatureID, Final_formula)

hilic_vk_raw <- separate_formula2(hilic_vk_raw, formula_col = "Final_formula", keep_formula = TRUE)
hilic_vk_raw <- calc_ratios_n_idxs(hilic_vk_raw)
hilic_vk_raw <- calc_classes(hilic_vk_raw)

VK_df_HILIC <- hilic_vk_raw %>%
  transmute(
    FeatureID,
    O_to_C = as.numeric(O_to_C),
    H_to_C = as.numeric(H_to_C),
    Class
  ) %>%
  filter(!is.na(O_to_C), !is.na(H_to_C))


# F) RP VAN KREVELEN ------
rp_vk_raw <- read_xlsx("Tables/RP_annot_features_final.xlsx") %>%
  filter(FeatureID %in% rp_top_feats_plateau) %>%
  dplyr::select(FeatureID, Final_formula)

rp_vk_raw <- separate_formula2(rp_vk_raw, formula_col = "Final_formula", keep_formula = TRUE) %>%
  mutate(C13 = 0, S = 0, P = 0, Na = 0)
rp_vk_raw <- calc_ratios_n_idxs(rp_vk_raw)
rp_vk_raw <- calc_classes(rp_vk_raw)

VK_df_RP <- rp_vk_raw %>%
  transmute(
    FeatureID,
    O_to_C = as.numeric(O_to_C),
    H_to_C = as.numeric(H_to_C),
    Class
  ) %>%
  filter(!is.na(O_to_C), !is.na(H_to_C))

hilic_kde <- kde2d(
  x = VK_df_HILIC$O_to_C, y = VK_df_HILIC$H_to_C,
  n = 200, lims = c(0, 1.5, 0, 2.5)
)
rp_kde <- kde2d(
  x = VK_df_RP$O_to_C, y = VK_df_RP$H_to_C,
  n = 200, lims = c(0, 1.5, 0, 2.5)
)
shared_density_max <- max(c(hilic_kde$z, rp_kde$z), na.rm = TRUE)

vk_plot_hilic <- ggplot(VK_df_HILIC, aes(x = O_to_C, y = H_to_C)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE, n = 200) +
  geom_density_2d(color = "black", linewidth = 0.35, bins = 6) +
  geom_rect(
    data = classification_boxes,
    aes(xmin = OC_low, xmax = OC_high, ymin = HC_low, ymax = HC_high),
    inherit.aes = FALSE, fill = NA, color = "grey20", linetype = "dashed", linewidth = 0.4
  ) +
  geom_text(
    data = classification_boxes,
    aes(x = label_x, y = label_y, label = Class),
    inherit.aes = FALSE, size = pt_to_mm(BASE_FONT_PT - 1), color = "grey15"
  ) +
  scale_fill_gradient(
    low = "#F7FBFF", high = "#2171B5",
    limits = c(0, shared_density_max), name = "Density"
  ) +
  labs(title = "HILIC RF-selected features", x = "O:C", y = "H:C") +
  figure_theme +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 2.5)) +
  theme(legend.position = "right")
vk_plot_hilic


vk_plot_rp <- ggplot(VK_df_RP, aes(x = O_to_C, y = H_to_C)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE, n = 200) +
  geom_density_2d(color = "black", linewidth = 0.35, bins = 6) +
  geom_rect(
    data = classification_boxes,
    aes(xmin = OC_low, xmax = OC_high, ymin = HC_low, ymax = HC_high),
    inherit.aes = FALSE, fill = NA, color = "grey20", linetype = "dashed", linewidth = 0.4
  ) +
  geom_text(
    data = classification_boxes,
    aes(x = label_x, y = label_y, label = Class),
    inherit.aes = FALSE, size = pt_to_mm(BASE_FONT_PT - 1), color = "grey15"
  ) +
  scale_fill_gradient(
    low = "#F7FBFF", high = "#2171B5",
    limits = c(0, shared_density_max), name = "Density"
  ) +
  labs(title = "RP RF-selected features", x = "O:C", y = "H:C") +
  figure_theme +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 2.5)) +
  theme(legend.position = "right")
vk_plot_rp


# # FINAL MULTIPANEL FIGURE -----------------------------------------------

figure2 <- (
  (pca_plot_hilic | pca_plot_rp) /
    (cm_plot_hilic  | cm_plot_rp) /
    (vk_plot_hilic  | vk_plot_rp)
) +
  plot_layout(guides = "collect", heights = c(1, 1, 1.15)) +
  plot_annotation(tag_levels = "a") &     # lowercase panel tags
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = TAG_FONT_PT)
  )

figure2

ggsave(
  filename = "Plots/New/Figure2_multimodal_overview.png",
  plot = figure2,
  width = 200,
  height = 180,
  units = "mm",
  bg = "white"
)
