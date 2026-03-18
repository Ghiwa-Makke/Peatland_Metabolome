suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
  library(scales)
  library(ggh4x)
})

# =========================================================
# 1) COMMON SETTINGS
# =========================================================

site_levels <- c("T1", "BLF", "S3", "US2", "LC")

type_colors <- c(
  "Bog" = "#8B4513",
  "Shrub Wetland" = "#2F8B8B",
  "Fen" = "#F5DEB3"
)

site_colors <- c(
  "T1"  = "#8B4513",
  "BLF" = "#D8C3A5",
  "S3"  = "#E8D8B8",
  "US2" = "#CDBA96",
  "LC"  = "#2F8B8B"
)
scfa_cols <- c(
  "Acetate"    = "#7570b3",
  "Propionate" = "#d95f02",
  "Butyrate"   = "#1b9e77"
)


my_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold", size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

# =========================================================
# 2) HELPERS
# =========================================================

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
      Site = factor(Site, levels = site_levels),
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
    meta_ord = meta_ord,
    site_centers = site_centers,
    site_breaks = site_breaks
  )
}

prep_one_category <- function(annot,
                              intensity_mat,
                              metadata,
                              category = c("Polyphenols", "Sugars"),
                              mode = c("HILIC", "RP")) {
  
  category <- match.arg(category)
  mode <- match.arg(mode)
  
  if (category == "Polyphenols") {
    guide <- annot %>%
      select(
        FeatureID, Final_name,
        ClassyFire_superclass, ClassyFire_class,
        ClassyFire_subclass, `ClassyFire_level 5`
      ) %>%
      filter(ClassyFire_superclass %in% c(
        "Benzenoids",
        "Phenylpropanoids and polyketides",
        "Lignans, neolignans and related compounds"
      )) %>%
      distinct(FeatureID, .keep_all = TRUE) %>%
      filter(!is.na(ClassyFire_class))
    
    if (mode == "RP") {
      guide <- guide %>%
        filter(ClassyFire_class != "Tetralins")
    }
    
    if (mode == "HILIC") {
      guide <- guide %>%
        filter(!ClassyFire_class %in% c(
          "Naphthalenes",
          "Fluorenes",
          "Macrolides and analogues",
          "NA"
        ))
    }
    
    group_var <- "ClassyFire_class"
  }
  
  if (category == "Sugars") {
    guide <- annot %>%
      select(
        FeatureID, Final_name,
        ClassyFire_superclass, ClassyFire_class,
        ClassyFire_subclass, `ClassyFire_level 5`
      ) %>%
      filter(`ClassyFire_level 5` %in% c(
        "Monosaccharides",
        "Disaccharides",
        "Sugar acids and derivatives",
        "Aminosaccharides"
      )) %>%
      distinct(FeatureID, .keep_all = TRUE) %>%
      filter(!is.na(`ClassyFire_level 5`))
    
    group_var <- "ClassyFire_level 5"
  }
  
  sample_order <- make_sample_order(metadata)
  
  hm_df <- intensity_mat %>%
    semi_join(guide, by = "FeatureID") %>%
    pivot_longer(
      cols = -FeatureID,
      names_to = "SampleID",
      values_to = "Intensity"
    ) %>%
    left_join(
      guide %>% select(FeatureID, all_of(group_var)),
      by = "FeatureID"
    ) %>%
    left_join(
      metadata %>% select(SampleID, Site, Type),
      by = "SampleID"
    ) %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(.data[[group_var]], SampleID, Site, Type) %>%
    summarise(Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Group = .data[[group_var]],
      Category = category,
      Mode = mode,
      SampleID = factor(SampleID, levels = sample_order),
      Intensity_plot = log1p(Intensity)
    )
  
  row_order <- hm_df %>%
    group_by(Group) %>%
    summarise(avg = mean(Intensity_plot, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg)) %>%
    pull(Group)
  
  hm_df %>%
    mutate(Group = factor(Group, levels = rev(row_order)))
}

prep_mode_heatmap <- function(annotation_file,
                              intensity_file,
                              metadata_file,
                              mode = c("HILIC", "RP")) {
  
  mode <- match.arg(mode)
  
  annot <- read_xlsx(annotation_file) %>%
    select(FeatureID, Final_name, starts_with("ClassyFire_"))
  
  intensity_mat <- read.csv(intensity_file, row.names = 1) %>%
    rownames_to_column(var = "FeatureID")
  
  metadata <- read.csv(metadata_file) %>%
    mutate(
      Site = factor(Site, levels = site_levels),
      Type = factor(Type, levels = names(type_colors))
    )
  
  polyphenols_df <- prep_one_category(
    annot = annot,
    intensity_mat = intensity_mat,
    metadata = metadata,
    category = "Polyphenols",
    mode = mode
  )
  
  sugars_df <- prep_one_category(
    annot = annot,
    intensity_mat = intensity_mat,
    metadata = metadata,
    category = "Sugars",
    mode = mode
  )
  
  hm_df <- bind_rows(polyphenols_df, sugars_df) %>%
    mutate(Category = factor(Category, levels = c("Polyphenols", "Sugars")))
  
  bar_info <- make_site_bar_df(metadata)
  
  list(
    hm_df = hm_df,
    metadata = metadata,
    meta_ord = bar_info$meta_ord,
    site_centers = bar_info$site_centers,
    site_breaks = bar_info$site_breaks
  )
}

make_mode_plot <- function(prepped,
                           mode_title,
                           fill_limits,
                           show_legend = FALSE) {
  
  p_bar <- ggplot(prepped$meta_ord, aes(x = SampleID, y = 1, fill = Site)) +
    geom_tile(height = 1) +
    geom_text(
      data = prepped$site_centers,
      aes(x = xmid, y = 1, label = Site),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(values = site_colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "none"#,
      #plot.margin = margin(2, 5.5, 0, 5.5)
    )
  
  p_hm <- ggplot(prepped$hm_df, aes(x = SampleID, y = Group, fill = Intensity_plot)) +
    geom_tile() +
    facet_grid(
      rows = vars(Category),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    geom_vline(
      xintercept = prepped$site_breaks,
      color = "white",
      linewidth = 0.7
    ) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,   # reversed
      limits = fill_limits,
      oob = scales::squish,
      name = "log1p(mean intensity)"
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = mode_title,
      x = NULL,
      y = NULL
    ) +
    my_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        angle = 90,
        face = "bold",
        size = 8
      ),
      legend.position = if (show_legend) "bottom" else "none",
      plot.margin = margin(0, 0, 0, 0)
    )
  
  p_bar / p_hm + plot_layout(heights = c(0.12, 1))
}

hilic_dat <- prep_mode_heatmap(
  annotation_file = "Data/LC/annot/Annotation_Filt_HILIC.xlsx",
  intensity_file  = "Tables/HILIC_norm_mean.csv",
  metadata_file   = "Tables/fixed_metadata_hilic.csv",
  mode            = "HILIC"
)

rp_dat <- prep_mode_heatmap(
  annotation_file = "Data/LC/annot/Annotation_Filt_updated_oct2025_RP.xlsx",
  intensity_file  = "Tables/RP_norm_mean.csv",
  metadata_file   = "Tables/fixed_metadata_rp.csv",
  mode            = "RP"
)

fill_limits <- range(
  c(hilic_dat$hm_df$Intensity_plot, rp_dat$hm_df$Intensity_plot),
  na.rm = TRUE
)

p_hilic <- make_mode_plot(
  prepped = hilic_dat,
  mode_title = "HILIC",
  fill_limits = fill_limits,
  show_legend = TRUE
)

p_rp <- make_mode_plot(
  prepped = rp_dat,
  mode_title = "RP",
  fill_limits = fill_limits,
  show_legend = FALSE
)




# =========================================================
# 5) SCFA PANEL
# =========================================================

scfa <- read_csv("Data/SCFA/Results_Edited_Ghiwa.csv", show_col_types = FALSE)

scfa_long <- scfa %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors)),
    across(c(Acetate, Propionate, Butyrate),
           ~ na_if(., "<LOQ") %>% as.numeric())
  ) %>%
  pivot_longer(
    cols = c(Acetate, Propionate, Butyrate),
    names_to = "SCFA",
    values_to = "Conc_mM"
  ) %>%
  mutate(
    SCFA = factor(SCFA, levels = c("Acetate", "Propionate", "Butyrate"))
  ) %>%
  filter(!is.na(Conc_mM), Conc_mM > 0)

site_strip_df <- tibble(
  Site = factor(site_levels, levels = site_levels)
)

strip_fills <- unname(site_colors[as.character(site_strip_df$Site)])

scfa_strip <- strip_themed(
  background_x = elem_list_rect(
    fill = strip_fills,
    colour = NA
  ),
  text_x = elem_list_text(
    face = "bold",
    colour = "black"
  )
)

p_scfa <- ggplot(scfa_long, aes(x = SCFA, y = Conc_mM)) +
  geom_boxplot(
    aes(fill = SCFA),
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.65
  ) +
  geom_jitter(
    aes(color = Type),
    width = 0.12,
    size = 1.8,
    alpha = 0.9
  ) +
  facet_wrap2(
    ~ Site,
    nrow = 1,
    strip = scfa_strip
  ) +
  scale_fill_manual(values = scfa_cols, drop = FALSE, name = "SCFA") +
  scale_color_manual(values = type_colors, drop = FALSE, name = "Type") +
  scale_y_log10(labels = label_number()) +
  labs(
    x = NULL,
    y = "Concentration (mM, log10 scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "left",
    axis.text.x = element_text(angle = 45, hjust = 1),
    #plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
    #plot.title = element_text(hjust = 0.5, face = "bold"),
    #strip.background = element_rect(fill = "white", colour = NA),
    #strip.text = element_text(face = "bold", size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(face = "bold"),
    #legend.direction = "vertical",
    panel.grid = element_blank()
  )#+
  #guides(fill = guide_legend(nrow = 2))
p_scfa

library(patchwork)
p_hilic <- p_hilic +
  theme(plot.margin = margin(0, 0, 0, 0), 
        legend.margin = margin(0, 0, 0, 0))

p_rp <- p_rp +
  theme(plot.margin = margin(0, 0, 0, 0))

p_scfa <- p_scfa +
  theme(plot.margin = margin(0, 0, 0, 0))

final_figure5 <- wrap_plots(
  wrap_elements(p_hilic),
  wrap_elements(p_rp),
  wrap_elements(p_scfa),
  ncol = 1, 
  heights = c(1, 1, 1)
) +
  plot_annotation(
    tag_levels = "a"
  )

final_figure5

ggsave(
  "Plots/Figure5_final_with_SCFA.png",
  final_figure5,
  width = 7,
  height = 9,
  dpi = 300,
  bg = "transparent"
)

ggsave(
  "Plots/Figure5_final_with_SCFA.svg",
  final_figure5,
  width = 7,
  height = 9,
  dpi = 300,
  bg = "transparent"
)

