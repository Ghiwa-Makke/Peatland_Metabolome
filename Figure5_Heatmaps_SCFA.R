suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
  library(scales)
  library(ggh4x)
  library(rstatix)
  library(ggpubr)
})

## font size 

BASE_FONT_PT  <- 9    # axis text, axis titles, legend text, site labels
TITLE_FONT_PT <- 10   # plot titles
TAG_FONT_PT   <- 12   # bold panel tags (a, b, c)

pt_to_mm <- function(pt) pt / .pt   # ggplot2 geom_text/annotate 'size' is in mm

#  COMMON SETTINGS

site_levels <- c("S1", "BLF", "S3", "US2", "LC")   

site_labels <- c(
  "S1"  = "S1-bog",
  "BLF" = "BLF-poor fen",
  "S3"  = "S3-fen",
  "US2" = "US2-fen",
  "LC"  = "LC-shrub wetland"
)

type_colors <- c(
  "Bog"          = "#8B4513",
  "Shrub Wetland"= "#2F8B8B",
  "Fen"          = "#F5DEB3"
)

site_colors <- c(
  "S1"  = "#8B4513",  
  "BLF" = "#F5DEB3",   
  "S3"  = "#F5DEB3",   
  "US2" = "#F5DEB3",   
  "LC"  = "#2F8B8B"  
)

scfa_cols <- c(
  "Acetate"    = "#7570b3",
  "Propionate" = "#d95f02",
  "Butyrate"   = "#1b9e77"
)

my_theme <- theme_classic(base_size = BASE_FONT_PT) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = TITLE_FONT_PT),
    strip.background  = element_rect(fill = "white", colour = NA),
    strip.text        = element_text(face = "bold", size = BASE_FONT_PT),
    axis.title        = element_text(face = "bold", size = BASE_FONT_PT),
    axis.text         = element_text(size = BASE_FONT_PT, color = "black"),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key        = element_rect(fill = "transparent", color = NA),
    legend.title      = element_text(face = "bold", size = BASE_FONT_PT),
    legend.text       = element_text(size = BASE_FONT_PT),
    panel.grid        = element_blank()
  )


# HELPERS
make_sample_order <- function(metadata) {
  metadata %>%
    mutate(Site = factor(Site, levels = site_levels)) %>%
    arrange(Site, SampleID) %>%
    pull(SampleID)
}

make_site_bar_df <- function(metadata) {
  sample_order <- make_sample_order(metadata)
  
  meta_ord <- metadata %>%
    mutate(
      Site     = factor(Site, levels = site_levels),
      SampleID = factor(SampleID, levels = sample_order)
    ) %>%
    arrange(SampleID) %>%
    mutate(x = seq_len(n()))
  
  site_centers <- meta_ord %>%
    group_by(Site) %>%
    summarise(
      xmin = min(x) - 0.5,
      xmax = max(x) + 0.5,
      xmid = mean(x),
      .groups = "drop"
    )
  
  site_breaks <- site_centers$xmax[-nrow(site_centers)]
  
  list(
    meta_ord     = meta_ord,
    site_centers = site_centers,
    site_breaks  = site_breaks
  )
}

prep_one_category <- function(annot, intensity_mat, metadata,
                              category = c("Polyphenols", "Sugars"),
                              mode     = c("HILIC", "RP")) {
  category <- match.arg(category)
  mode     <- match.arg(mode)
  
  if (category == "Polyphenols") {
    guide <- annot %>%
      select(FeatureID, Final_name,
             ClassyFire_superclass, ClassyFire_class,
             ClassyFire_subclass, `ClassyFire_level 5`) %>%
      filter(ClassyFire_superclass %in% c(
        "Benzenoids",
        "Phenylpropanoids and polyketides",
        "Lignans, neolignans and related compounds"
      )) %>%
      distinct(FeatureID, .keep_all = TRUE) %>%
      filter(!is.na(ClassyFire_class))
    
    if (mode == "RP")    guide <- guide %>% filter(ClassyFire_class != "Tetralins")
    if (mode == "HILIC") guide <- guide %>% filter(!ClassyFire_class %in% c(
      "Naphthalenes", "Fluorenes", "Macrolides and analogues", "NA"
    ))
    
    group_var <- "ClassyFire_class"
  }
  
  if (category == "Sugars") {
    guide <- annot %>%
      select(FeatureID, Final_name,
             ClassyFire_superclass, ClassyFire_class,
             ClassyFire_subclass, `ClassyFire_level 5`) %>%
      filter(`ClassyFire_level 5` %in% c(
        "Monosaccharides", "Disaccharides",
        "Sugar acids and derivatives", "Aminosaccharides"
      )) %>%
      distinct(FeatureID, .keep_all = TRUE) %>%
      filter(!is.na(`ClassyFire_level 5`))
    
    group_var <- "ClassyFire_level 5"
  }
  
  sample_order <- make_sample_order(metadata)
  
  hm_df <- intensity_mat %>%
    semi_join(guide, by = "FeatureID") %>%
    pivot_longer(cols = -FeatureID, names_to = "SampleID", values_to = "Intensity") %>%
    left_join(guide %>% select(FeatureID, all_of(group_var)), by = "FeatureID") %>%
    left_join(metadata %>% select(SampleID, Site, Type), by = "SampleID") %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(.data[[group_var]], SampleID, Site, Type) %>%
    summarise(Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Group          = .data[[group_var]],
      Category       = category,
      Mode           = mode,
      SampleID       = factor(SampleID, levels = sample_order),
      ## ── log10 transform ──────────────────────────────────
      Intensity_plot = log10(Intensity + 1)
    )
  
  row_order <- hm_df %>%
    group_by(Group) %>%
    summarise(avg = mean(Intensity_plot, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg)) %>%
    pull(Group)
  
  hm_df %>% mutate(Group = factor(Group, levels = rev(row_order)))
}

prep_mode_heatmap <- function(annotation_file, intensity_file,
                              metadata_file, mode = c("HILIC", "RP")) {
  mode <- match.arg(mode)
  
  annot <- read_xlsx(annotation_file) %>%
    select(FeatureID, Final_name, starts_with("ClassyFire_"))
  
  intensity_mat <- read.csv(intensity_file, row.names = 1) %>%
    rownames_to_column(var = "FeatureID")
  
  metadata <- read.csv(metadata_file) %>%
    mutate(
      Site = recode(Site, "T1" = "S1"),           # T1 -> S1
      Site = factor(Site, levels = site_levels),
      Type = factor(Type, levels = names(type_colors))
    )
  
  polyphenols_df <- prep_one_category(annot, intensity_mat, metadata, "Polyphenols", mode)
  sugars_df      <- prep_one_category(annot, intensity_mat, metadata, "Sugars",      mode)
  
  hm_df <- bind_rows(polyphenols_df, sugars_df) %>%
    mutate(Category = factor(Category, levels = c("Polyphenols", "Sugars")))
  
  bar_info <- make_site_bar_df(metadata)
  
  list(
    hm_df        = hm_df,
    metadata     = metadata,
    meta_ord     = bar_info$meta_ord,
    site_centers = bar_info$site_centers,
    site_breaks  = bar_info$site_breaks
  )
}

make_mode_plot <- function(prepped, mode_title, fill_limits, show_legend = FALSE) {
  
  # extended site names for in-plot bar labels
  site_centers_labeled <- prepped$site_centers %>%
    mutate(SiteLabel = site_labels[as.character(Site)])
  
  p_bar <- ggplot(prepped$meta_ord, aes(x = SampleID, y = 1, fill = Site)) +
    geom_tile(height = 1) +
    geom_text(
      data = site_centers_labeled,
      aes(x = xmid, y = 1, label = SiteLabel),      # extended names in-plot
      inherit.aes = FALSE,
      size = pt_to_mm(BASE_FONT_PT - 2),             # slightly smaller: labels are longer now
      fontface = "bold",
      angle = 0
    ) +
    scale_fill_manual(values = site_colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    theme_void() +
    theme(legend.position = "none")
  
  p_hm <- ggplot(prepped$hm_df, aes(x = SampleID, y = Group, fill = Intensity_plot)) +
    geom_tile() +
    facet_grid(rows = vars(Category), scales = "free_y",
               space = "free_y", switch = "y") +
    geom_vline(xintercept = prepped$site_breaks,
               color = "white", linewidth = 0.7) +
    scale_fill_viridis_c(
      option    = "magma",
      direction = -1,
      limits    = fill_limits,
      oob       = scales::squish,
      name      = expression(log[10] ~ "(mean intensity)")   
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(title = mode_title, x = NULL, y = NULL) +
    my_theme +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = BASE_FONT_PT - 1),
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 90, face = "bold", size = BASE_FONT_PT - 1),
      legend.position    = if (show_legend) "bottom" else "none",
      legend.text        = element_text(size = BASE_FONT_PT - 1),
      plot.margin        = margin(0, 0, 0, 0)
    )
  
  p_bar / p_hm + plot_layout(heights = c(0.12, 1))
}


# LOAD & PLOT HEATMAPS


hilic_dat <- prep_mode_heatmap(
  annotation_file = "Input_Data/LC/annot/Annotation_HILIC.xlsx",
  intensity_file  = "Tables/HILIC_norm_mean.csv",
  metadata_file   = "Tables/fixed_metadata_hilic.csv",
  mode            = "HILIC"
)

rp_dat <- prep_mode_heatmap(
  annotation_file = "Input_Data/LC/annot/Annotation_RP.xlsx",
  intensity_file  = "Tables/RP_norm_mean.csv",
  metadata_file   = "Tables/fixed_metadata_rp.csv",
  mode            = "RP"
)

## Stats - supplementary tables ---------

sig_class_stats <- function(prepped, mode_label) {
  
  pw <- prepped$hm_df %>%
    group_by(Category, Group) %>%
    wilcox_test(Intensity ~ Site, p.adjust.method = "BH") %>%
    add_significance("p.adj") %>%
    ungroup()
  
  medians <- prepped$hm_df %>%
    group_by(Category, Group, Site) %>%
    summarise(median_val = median(Intensity, na.rm = TRUE), .groups = "drop")
  
  pw %>%
    left_join(
      medians %>% rename(group1 = Site, median1 = median_val),
      by = c("Category", "Group", "group1")
    ) %>%
    left_join(
      medians %>% rename(group2 = Site, median2 = median_val),
      by = c("Category", "Group", "group2")
    ) %>%
    mutate(
      Mode = mode_label,
      Direction = case_when(
        p.adj >= 0.05 ~ "ns",
        median1 > median2 ~ paste0(group1, " > ", group2),
        median1 < median2 ~ paste0(group2, " > ", group1),
        TRUE ~ "ns"
      )
    ) %>%
    select(Mode, Category, Variable = Group, Site1 = group1, Site2 = group2,
           n1, n2, median1, median2, statistic, p, p.adj, p.adj.signif, Direction)
}

kw_class_stats <- function(prepped, mode_label) {
  prepped$hm_df %>%
    group_by(Category, Group) %>%
    kruskal_test(Intensity ~ Site) %>%
    add_significance("p") %>%
    ungroup() %>%
    mutate(Mode = mode_label) %>%
    select(Mode, Category, Variable = Group, n, statistic, df, p, p.signif)
}

kw_hilic_class <- kw_class_stats(hilic_dat, "HILIC")
kw_rp_class    <- kw_class_stats(rp_dat,    "RP")

sig_hilic_class <- sig_class_stats(hilic_dat, "HILIC")
sig_rp_class    <- sig_class_stats(rp_dat,    "RP")

# Table A: ALL pairwise comparisons ----------
sup_table3_all <- bind_rows(sig_hilic_class, sig_rp_class) %>%
  select(Mode, Category, Variable, Site1, Site2, median1, median2, p, p.adj, p.adj.signif, Direction) %>%
  arrange(Category, Variable, Site1, Site2)

# Table B: SIGNIFICANT pairwise comparisons only ----
sup_table3_sig <- sup_table3_all %>%
  filter(p.adj < 0.05, p.adj.signif != "ns")

# Omnibus Kruskal-Wallis summary, one row per class
sup_table3_omnibus <- bind_rows(kw_hilic_class, kw_rp_class) %>%
  arrange(Category, p)

write_csv(sup_table3_all,     "Tables/SupTable3_class_pairwise_ALL.csv")
write_csv(sup_table3_sig,     "Tables/SupTable3_class_pairwise_SIGNIFICANT.csv")
write_csv(sup_table3_omnibus, "Tables/SupTable3_class_kruskalwallis_omnibus.csv")

# SCFA WITH SIGNIFICANCE ---------------

fill_limits <- range(
  c(hilic_dat$hm_df$Intensity_plot, rp_dat$hm_df$Intensity_plot),
  na.rm = TRUE
)

p_hilic <- make_mode_plot(hilic_dat, "HILIC", fill_limits, show_legend = TRUE)
p_rp    <- make_mode_plot(rp_dat,    "RP",    fill_limits, show_legend = FALSE)

scfa <- read_csv("Input_Data/SCFA/SCFA_Results.csv", show_col_types = FALSE)

scfa_long <- scfa %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors)),
    across(c(Acetate, Propionate, Butyrate),
           ~ na_if(., "<LOQ") %>% as.numeric())
  ) %>%
  pivot_longer(
    cols      = c(Acetate, Propionate, Butyrate),
    names_to  = "SCFA",
    values_to = "Conc_mM"
  ) %>%
  mutate(SCFA = factor(SCFA, levels = c("Acetate", "Propionate", "Butyrate"))) %>%
  filter(!is.na(Conc_mM), Conc_mM > 0)

## pairwise significance per SCFA (across Sites)---------
sig_scfa <- scfa_long %>%
  group_by(SCFA) %>%
  wilcox_test(Conc_mM ~ Site, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05, p.adj.signif != "ns") %>%
  add_xy_position(x = "Site") %>%
  ungroup()

## coloring the strips
site_strip_df <- tibble(Site = factor(site_levels, levels = site_levels)) %>%
  mutate(SiteLabel = site_labels[as.character(Site)])
strip_fills   <- unname(site_colors[as.character(site_strip_df$Site)])

scfa_strip <- strip_themed(
  background_x = elem_list_rect(fill = strip_fills, colour = NA),
  text_x       = elem_list_text(face = "bold", colour = "black", size = BASE_FONT_PT - 1)
)

p_scfa <- ggplot(
  scfa_long %>% mutate(SiteLabel = factor(site_labels[as.character(Site)],
                                          levels = site_labels[site_levels])),
  aes(x = SCFA, y = Conc_mM)
) +
  geom_boxplot(
    aes(fill = SCFA),
    outlier.shape = NA, alpha = 0.75, width = 0.65
  ) +
  geom_jitter(
    aes(color = Type),
    width = 0.12, size = 1.8, alpha = 0.9
  ) +
  { if (nrow(sig_scfa) > 0)
    stat_pvalue_manual(
      sig_scfa %>% mutate(SiteLabel = factor(site_labels[as.character(Site)],
                                             levels = site_labels[site_levels])),
      label         = "p.adj.signif",
      tip.length    = 0.01,
      size          = pt_to_mm(BASE_FONT_PT - 1),
      step.increase = 0.08
    )
  } +
  facet_wrap2(~ SiteLabel, nrow = 1) +               # extended names on facet strip
  scale_fill_manual(values  = scfa_cols,   drop = FALSE, name = "SCFA") +
  scale_color_manual(values = type_colors, drop = FALSE, name = "Type") +
  scale_y_log10(labels = label_number(), expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "Concentration (mM, log10 scale)") +
  theme_bw(base_size = BASE_FONT_PT) +
  theme(
    legend.position   = "left",
    legend.text       = element_text(size = BASE_FONT_PT),
    legend.title      = element_text(face = "bold", size = BASE_FONT_PT),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = BASE_FONT_PT),
    axis.text.y       = element_text(size = BASE_FONT_PT),
    axis.title        = element_text(face = "bold", size = BASE_FONT_PT),
    strip.text        = element_text(face = "bold", size = BASE_FONT_PT - 1),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key        = element_rect(fill = "transparent", color = NA),
    panel.grid        = element_blank()
  )


# final 

p_hilic <- p_hilic + theme(plot.margin = margin(0, 0, 0, 0))
p_rp    <- p_rp    + theme(plot.margin = margin(0, 0, 0, 0))
p_scfa  <- p_scfa  + theme(plot.margin = margin(0, 0, 0, 0))

final_figure5 <- wrap_plots(
  wrap_elements(p_hilic),
  wrap_elements(p_rp),
  wrap_elements(p_scfa),
  ncol    = 1,
  heights = c(1, 1, 1)
) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = TAG_FONT_PT))

final_figure5

ggsave("Plots/Figure5_final_with_SCFA.png", final_figure5,
       width = 174, height = 220, units = "mm", dpi = 300, bg = "white")
ggsave("Plots/Figure5_final_with_SCFA.svg", final_figure5,
       width = 174, height = 220, units = "mm", bg = "white")
