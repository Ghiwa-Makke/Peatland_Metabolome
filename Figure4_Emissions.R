## Peatland - Metabolome Project - 2023 Samples
## Emissions + Whiticar multi-panel with Panel D

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggforce)
  library(rstatix)
  library(ggpubr)
  library(ggsignif)
  library(scales)
})

## font size

BASE_FONT_PT  <- 9    # axis text, axis titles, legend text
TITLE_FONT_PT <- 10   # plot titles
TAG_FONT_PT   <- 12   # bold panel tags (a-f)

pt_to_mm <- function(pt) pt / .pt   # ggplot2 geom_text/annotate 'size' is in mm


## Data -------

emission <- read.csv("Input_Data/Emissions/Emission_Data.csv")

site_levels <- c("S1", "BLF", "S3", "US2", "LC")

# short codes for axes; full names for legends only
site_labels <- c(
  "S1"  = "S1-bog",
  "BLF" = "BLF-poor fen",
  "S3"  = "S3-fen",
  "US2" = "US2-fen",
  "LC"  = "LC-shrub wetland"
)

type_colors <- c(
  "Bog" = "#8B4513",
  "Shrub Wetland" = "#2F8B8B",
  "Fen" = "#CC9B7AFF"
)

site_shapes <- c(
  "S1"  = 21,   
  "BLF" = 24,
  "S3"  = 22,
  "US2" = 23,
  "LC"  = 25
)


site_type <- emission %>%
  select(Site, Type) %>%
  distinct() %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )


## Common theme ------

my_theme <- theme_classic(base_size = BASE_FONT_PT) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = TITLE_FONT_PT),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold", size = BASE_FONT_PT),
    axis.title = element_text(size = BASE_FONT_PT, face = "bold"),
    axis.text = element_text(size = BASE_FONT_PT, color = "black"),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(face = "bold", size = BASE_FONT_PT),
    legend.text = element_text(size = BASE_FONT_PT),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = margin(1, 1, 1, 1)
  )



## Panel a: DIC concentration --------------

dic_data <- emission %>%
  rename(DIC = CO2mM) %>%
  select(SampleID, Site, Type, DIC, Subsite, rep) %>%
  filter(!is.na(DIC)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

sig_dic <- dic_data %>%
  wilcox_test(DIC ~ Site, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05, p.adj.signif != "ns") %>%
  add_xy_position(x = "Site", step.increase = 0.1)

DIC_plot <- ggplot(dic_data, aes(x = Site, y = DIC, fill = Type)) +
  geom_boxplot(linewidth = 0.5) +
  stat_pvalue_manual(sig_dic, label = "p.adj.signif", tip.length = 0.01,
                     size = pt_to_mm(BASE_FONT_PT - 1)) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  labs(x = "Site", y = "DIC (mM)", title = "DIC concentration") +
  my_theme

DIC_plot


## Panel b: CH4 concentration ------------------

ch4_data <- emission %>%
  rename(CH4 = CH4mM) %>%
  select(SampleID, Site, Type, CH4, Subsite, rep) %>%
  filter(!is.na(CH4)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

sig_ch4 <- ch4_data %>%
  wilcox_test(CH4 ~ Site, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05, p.adj.signif != "ns") %>%
  add_xy_position(x = "Site", step.increase = 0.1)

CH4_plot <- ggplot(ch4_data, aes(x = Site, y = CH4, fill = Type)) +
  geom_boxplot(linewidth = 0.5) +
  stat_pvalue_manual(sig_ch4, label = "p.adj.signif", tip.length = 0.01,
                     size = pt_to_mm(BASE_FONT_PT - 1)) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  labs(x = "Site", y = expression(bold(CH[4] ~ "(mM)")), title = expression(bold(CH[4] ~ "concentration"))) +
  my_theme

CH4_plot


## Panel c: d13C-CO2 -------------------

co2_iso_data <- emission %>%
  select(SampleID, Site, Type, d13CO2, Subsite, rep) %>%
  filter(!is.na(d13CO2)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

sig_co2_iso <- co2_iso_data %>%
  wilcox_test(d13CO2 ~ Site, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05, p.adj.signif != "ns") %>%
  add_xy_position(x = "Site", step.increase = 0.1)

co2_isotope_plot <- ggplot(co2_iso_data, aes(x = Site, y = d13CO2, fill = Type)) +
  geom_boxplot(linewidth = 0.5) +
  stat_pvalue_manual(sig_co2_iso, label = "p.adj.signif", tip.length = 0.01,
                     size = pt_to_mm(BASE_FONT_PT - 1)) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  labs(
    x = "Site",
    y = expression(bold(delta^13 * C ~ "(\u2030)")),
    title = "\u03B4\u00B9\u00B3C-CO\u2082"     # δ¹³C-CO₂
  ) +
  my_theme

co2_isotope_plot


## 6) Panel d: d13C-CH4 ---------------

ch4_iso_data <- emission %>%
  select(SampleID, Site, Type, d13CH4, Subsite, rep) %>%
  filter(!is.na(d13CH4)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

sig_ch4_iso <- ch4_iso_data %>%
  wilcox_test(d13CH4 ~ Site, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05, p.adj.signif != "ns") %>%
  add_xy_position(x = "Site", step.increase = 0.1)

ch4_isotope_plot <- ggplot(ch4_iso_data, aes(x = Site, y = d13CH4, fill = Type)) +
  geom_boxplot(linewidth = 0.5) +
  stat_pvalue_manual(sig_ch4_iso, label = "p.adj.signif", tip.length = 0.01,
                     size = pt_to_mm(BASE_FONT_PT - 1)) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  labs(
    x = "Site",
    y = expression(bold(delta^13 * C ~ "(\u2030)")),
    title = "\u03B4\u00B9\u00B3C-CH\u2084"     # δ¹³C-CH₄
  ) +
  my_theme

ch4_isotope_plot

## Panel e: alpha ---------------
alpha_data <- emission %>%
  select(SampleID, Site, Type, alpha, Subsite, rep) %>%
  filter(!is.na(alpha)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

sig_alpha <- alpha_data %>%
  wilcox_test(alpha ~ Site, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  filter(p.adj < 0.05, p.adj.signif != "ns") %>%
  add_xy_position(x = "Site", step.increase = 0.05)

alpha_plot <- ggplot(alpha_data, aes(x = Site, y = alpha, fill = Type)) +
  geom_boxplot(linewidth = 0.5) +
  stat_pvalue_manual(sig_alpha, label = "p.adj.signif", tip.length = 0.01,
                     size = pt_to_mm(BASE_FONT_PT - 1)) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  scale_y_continuous(
    n.breaks = 10,                                  # finer gridlines than default
    labels = scales::number_format(accuracy = 0.001)  # more decimal precision
  ) +
  labs(x = "Site", y = expression(bold(alpha))) +
  my_theme +
  theme(legend.position = "right")

alpha_plot


## Panel f: Whiticar plot --------------

df_whiticar <- emission %>%
  select(SampleID, Site, Type, d13CO2, d13CH4, alpha) %>%
  filter(!is.na(d13CO2), !is.na(d13CH4)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

## Original Whiticar fields
carbonate_reduction <- tibble(
  field = "Carbonate\nReduction",
  x = c(-101, -95, -82, -70, -60, -55, -50, -55, -70, -101),
  y = c(-28,  -29, -30, -20, -10,   -3,  10,  15,  10,  -10)
)

methyl_fermentation <- tibble(
  field = "Methyl\nFermentation",
  x = c(-67, -66, -57, -49, -48, -51, -54, -56, -61, -65, -67, -67),
  y = c(-28, -28, -25, -19, -11, -5, -2.5, -2.5, -8, -16.5, -25, -28)
)

methane_oxidation_field <- tibble(
  field = "Methane\nOxidation",
  x = c(-65, -60, -54, -50, -22, -32, -65),
  y = c(-37, -43, -45, -45, -17,  -3, -37)
)

field_labels <- tribble(
  ~field,                    ~x,   ~y,
  "Carbonate\nReduction",   -80,  -15,
  "Methyl\nFermentation",   -45, -2,
  "Methane\nOxidation",     -32, -15
)

alphas <- c(1.005, 1.03, 1.055)

xlim_use <- c(-100, -20)
ylim_use <- c(-50, 20)

alpha_df <- tidyr::crossing(
  alpha = alphas,
  x = seq(xlim_use[1], xlim_use[2], by = 0.5)
) %>%
  mutate(
    y = alpha * (x + 1000) - 1000,
    alpha_lab = paste0("\u03B1 = ", format(alpha, trim = TRUE))
  ) %>%
  filter(y >= ylim_use[1], y <= ylim_use[2])

alpha_labels <- alpha_df %>%
  group_by(alpha, alpha_lab) %>%
  slice_min(order_by = x, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    x = x + 1.2,
    y = y + 1.2
  )

## Helper to smooth polygon edges
smooth_polygon <- function(df, n = 300) {
  df_closed <- df
  
  if (df_closed$x[1] != df_closed$x[nrow(df_closed)] ||
      df_closed$y[1] != df_closed$y[nrow(df_closed)]) {
    df_closed <- bind_rows(df_closed, df_closed[1, ])
  }
  
  tibble(
    x = spline(seq_along(df_closed$x), df_closed$x, n = n, method = "natural")$y,
    y = spline(seq_along(df_closed$y), df_closed$y, n = n, method = "natural")$y
  )
}

carbonate_s <- carbonate_reduction %>%
  select(x, y) %>%
  smooth_polygon()

methyl_s <- methyl_fermentation %>%
  select(x, y) %>%
  smooth_polygon()

methox_s <- methane_oxidation_field %>%
  select(x, y)

## Production arrow
production_arrow <- tibble(
  x = c(-100, -90, -80, -70),
  y = c(-18,  -11,  -2,   8)
)

## Oxidation arrow
oxidation_arrow <- tibble(
  x = c(-100, -82, -55, -30),
  y = c(-24,  -24, -24, -22)
)

process_labels <- tibble(
  label = c("Production", "Oxidation"),
  x = c(-88, -68),
  y = c(-8, -28),
  angle = c(25, 0)
)

## Plot theme
whiticar_theme <- theme_classic(base_size = BASE_FONT_PT) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = BASE_FONT_PT),
    axis.text = element_text(color = "black", size = BASE_FONT_PT),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = BASE_FONT_PT),
    legend.text = element_text(size = BASE_FONT_PT),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )

## panel f - whiticar_panel
whiticar_panel <- ggplot() +
  ## Alpha lines
  geom_line(
    data = alpha_df,
    aes(x = x, y = y, group = alpha),
    linewidth = 0.8,
    color = "black",
    alpha = 0.3
  ) +
  geom_text(
    data = alpha_labels,
    aes(x = x, y = y, label = alpha_lab),
    size = pt_to_mm(BASE_FONT_PT - 1),
    hjust = 0.2,
    vjust = 0
  ) +
  ## Background mechanistic fields
  geom_polygon(
    data = carbonate_s,
    aes(x = x, y = y),
    fill = "grey85",
    color = "black",
    linewidth = 0.45,
    alpha = 1
  ) +
  geom_polygon(
    data = methyl_s,
    aes(x = x, y = y),
    fill = "grey85",
    color = "black",
    linewidth = 0.45,
    alpha = 1
  ) +
  geom_polygon(
    data = methox_s,
    aes(x = x, y = y),
    fill = "grey85",
    color = "black",
    linewidth = 0.45,
    alpha = 1
  ) +
  
  ## Production arrow
  geom_bezier(
    data = production_arrow,
    aes(x = x, y = y),
    linewidth = 6,
    color = "#9fd3f2",
    alpha = 0.9,
    lineend = "round"
  ) +
  geom_segment(
    aes(x = -75, y = 5, xend = -65, yend = 10),
    arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
    linewidth = 1.5,
    color = "#9fd3f2"
  ) +
  
  ## Oxidation arrow
  geom_curve(
    aes(x = -100, y = -20, xend = -40, yend = -20),
    curvature = 0.25,
    linewidth = 6,
    color = "#f4b6c2",
    lineend = "round"
  ) +
  geom_segment(
    aes(x = -40, y = -20, xend = -35, yend = -18),
    arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
    linewidth = 1.5,
    color = "#f4b6c2"
  ) +
  
  ## Region labels
  geom_text(
    data = field_labels,
    aes(x = x, y = y, label = field),
    size = pt_to_mm(BASE_FONT_PT),
    lineheight = 0.95
  ) +
  
  ## Process labels
  geom_text(
    data = process_labels,
    aes(x = x, y = y, label = label, angle = angle),
    fontface = "bold",
    size = pt_to_mm(TITLE_FONT_PT)
  ) +
  
  ## Sample points
  geom_point(
    data = df_whiticar,
    aes(x = d13CH4, y = d13CO2, fill = Type, shape = Site),
    size = 1.5,
    alpha = 1
  ) +
  
  ## Axes and scales
  scale_fill_manual(
    values = type_colors,
    name = "Classification",
    drop = FALSE
  ) +
  scale_shape_manual(
    values = site_shapes,
    name = "Site",
    labels = site_labels,     
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(-101, -20),
    breaks = seq(-100, -20, 10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-50, 20),
    breaks = seq(-50, 20, 10),
    expand = c(0, 0)
  ) +
  labs(
    x = "\u03B4\u00B9\u00B3C-CH\u2084 (\u2030)",
    y = "\u03B4\u00B9\u00B3C-CO\u2082 (\u2030)"
  ) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  whiticar_theme

## Show plot
whiticar_panel


## Final multi-panel figure a-f --------------------
##    (a) DIC (b) CH4 (c) d13C-CO2 (d) d13C-CH4 (e) alpha (f) whiticar


overall <- (
  (DIC_plot | CH4_plot) /
    (co2_isotope_plot | ch4_isotope_plot) /
    (alpha_plot | whiticar_panel)
) +
  plot_layout(heights = c(1, 1, 1.2), guides = "collect") +
  plot_annotation(tag_levels = "a") &   
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = TAG_FONT_PT)
  )

overall



top_row <- (DIC_plot | CH4_plot)
mid_row <- (co2_isotope_plot | ch4_isotope_plot) 

bottom_row <- alpha_plot + whiticar_panel +
  plot_layout(widths = c(1.5, 2.5))



overall_2 <- (top_row / mid_row /bottom_row)+
  plot_layout(heights = c(1, 1, 1.2),guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right")


ggsave(
  "Plots/New/Figure4_Emissions_with_Whiticar_stat_vertical.png",
  overall_2,
  width  = 174,
  height = 200,
  units  = "mm",
  dpi    = 300,
  bg     = "transparent"
)



## Export supplementary statistics tables (all pairwise comparisons)

run_site_pairwise_stats <- function(data, value_col, variable_name) {
  
  pw <- data %>%
    wilcox_test(as.formula(paste(value_col, "~ Site")), p.adjust.method = "BH") %>%
    add_significance("p.adj")
  
  medians <- data %>%
    group_by(Site) %>%
    summarise(median_val = median(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  
  pw %>%
    left_join(medians %>% rename(group1 = Site, median1 = median_val), by = "group1") %>%
    left_join(medians %>% rename(group2 = Site, median2 = median_val), by = "group2") %>%
    mutate(
      Variable = variable_name,
      Direction = case_when(
        p.adj >= 0.05 ~ "ns",
        median1 > median2 ~ paste0(group1, " > ", group2),
        median1 < median2 ~ paste0(group2, " > ", group1),
        TRUE ~ "ns"
      )
    ) %>%
    select(Variable, group1, group2, n1, n2, median1, median2,
           statistic, p, p.adj, p.adj.signif, Direction)
}

# Run for all 
stats_dic     <- run_site_pairwise_stats(dic_data,     "DIC",    "DIC (mM)")
stats_ch4     <- run_site_pairwise_stats(ch4_data,     "CH4",    "CH4 (mM)")
stats_co2_iso <- run_site_pairwise_stats(co2_iso_data, "d13CO2", "delta13C-CO2 (per mil)")
stats_ch4_iso <- run_site_pairwise_stats(ch4_iso_data, "d13CH4", "delta13C-CH4 (per mil)")
stats_alpha   <- run_site_pairwise_stats(alpha_data,   "alpha",  "alpha")

# Combine into one supplementary table
figure4_supp_stats <- bind_rows(
  stats_dic, stats_ch4, stats_co2_iso, stats_ch4_iso, stats_alpha
)

write_csv(
  figure4_supp_stats,
  "Tables/SupplementaryTable_Figure4_EmissionStats.csv"
)
