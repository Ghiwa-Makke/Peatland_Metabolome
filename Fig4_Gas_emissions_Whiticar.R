## Peatland - Metabolome Project - 2023 Samples
## Emissions + Whiticar multi-panel with Panel D

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

## ------------------------------------------------------------
## 1) Data
## ------------------------------------------------------------
emission <- read.csv("Data/Emissions/Emission_Data.csv")

site_levels <- c("T1", "BLF", "S3", "US2", "LC")

type_colors <- c(
  "Bog" = "#8B4513",
  "Shrub Wetland" = "#2F8B8B",
  "Fen" = "#CC9B7AFF"
)

site_shapes <- c(
  "T1"  = 16,
  "BLF" = 17,
  "S3"  = 15,
  "US2" = 9,
  "LC"  = 8
)

site_type <- emission %>%
  select(Site, Type) %>%
  distinct() %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

## ------------------------------------------------------------
## 2) Common theme
## ------------------------------------------------------------
my_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold", size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5)
  )

## ------------------------------------------------------------
## 3) Panel A: CO2 & CH4 concentrations
## ------------------------------------------------------------
gas_data <- emission %>%
  rename(CH4 = CH4mM, CO2 = CO2mM) %>%
  select(SampleID, Site, Type, CO2, CH4, Subsite, rep) %>%
  pivot_longer(cols = c(CO2, CH4), names_to = "gas", values_to = "mM") %>%
  filter(!is.na(mM)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors)),
    gas = factor(
      gas,
      levels = c("CO2", "CH4"),
      labels = c(paste0("CO", "\u2082"), paste0("CH", "\u2084"))
    )
  )

CO2_CH4 <- ggplot(gas_data, aes(x = Site, y = mM, fill = Type)) +
  geom_boxplot(size = 0.5) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  facet_wrap(~ gas, scales = "free_y") +
  labs(
    x = "Site",
    y = "Concentration (mM)",
    title = expression(CO[2] ~ "and" ~ CH[4] ~ "concentrations")
  ) +
  my_theme

## ------------------------------------------------------------
## 4) Panel B: δ13C
## ------------------------------------------------------------
isotope_data <- emission %>%
  select(SampleID, Site, Type, d13CH4, d13CO2, Subsite, rep) %>%
  pivot_longer(
    cols = c(d13CH4, d13CO2),
    names_to = "gas",
    values_to = "d13C"
  ) %>%
  filter(!is.na(d13C)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

isotope_plot <- ggplot(isotope_data, aes(x = Site, y = d13C, fill = Type)) +
  geom_boxplot(size = 0.5) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  facet_wrap(~ gas, scales = "free_y") +
  labs(
    x = "Site",
    y = expression(delta^13 * C),
    title = expression(delta^13 * C ~ "of CO"[2] * " and CH"[4])
  ) +
  my_theme

## ------------------------------------------------------------
## 5) Panel C: alpha
## ------------------------------------------------------------
alpha_data <- emission %>%
  select(SampleID, Site, Type, alpha, Subsite, rep) %>%
  filter(!is.na(alpha)) %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    Type = factor(Type, levels = names(type_colors))
  )

alpha_plot <- ggplot(alpha_data, aes(x = Site, y = alpha, fill = Type)) +
  geom_boxplot(size = 0.5) +
  scale_fill_manual(values = type_colors, name = "Classification", drop = FALSE) +
  labs(
    x = "Site",
    y = expression(alpha),
    #title = expression(alpha ~ "(" * alpha * ")")
  ) +
  my_theme +
  theme(legend.position = "right")

## ------------------------------------------------------------
## 6) Panel D: Whiticar plot
## ------------------------------------------------------------
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
  "Methyl\nFermentation",   -45, -3,
  "Methane\nOxidation",     -32, -15
)

alphas <- c(1.005, 1.03, 1.04, 1.055)

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


## Production arrow: oblique smooth arrow
production_arrow <- tibble(
  x = c(-100, -90, -80, -70),
  y = c(-18,  -11,  -2,   8)
)

## Oxidation arrow: mostly horizontal with slight curve
oxidation_arrow <- tibble(
  x = c(-100, -82, -55, -30),
  y = c(-24,  -24, -24, -22)
)

process_labels <- tibble(
  label = c("Production", "Oxidation"),
  x = c(-90, -68),
  y = c(-8, -32),
  angle = c(25, 0)
)

## Plot theme
whiticar_theme <- theme_classic(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    #plot.margin = margin(8, 8, 8, 8),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )
## panel D - whiticar_panel 
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
    size = 2,
    hjust = 0,
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
    curvature = 0.25,  # subtle bend near the end
    linewidth = 6,
    color = "#f4b6c2",
    lineend = "round"
  )+
  geom_segment(
    aes(x = -40, y =-20, xend = -35, yend = -18),
    arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
    linewidth = 1.5,
    color = "#f4b6c2"
  ) +

  ## Region labels
  geom_text(
    data = field_labels,
    aes(x = x, y = y, label = field),
    size = 2.5,
    lineheight = 0.95
  ) +
  
  ## Process labels
  geom_text(
    data = process_labels,
    aes(x = x, y = y, label = label, angle = angle),
    fontface = "bold",
    size = 3
  ) +
  
  ## Sample points
  geom_point(
    data = df,
    aes(x = d13CH4, y = d13CO2, color = Type, shape = Site),
    size = 1.5,
    alpha = 1
  ) +
  
  ## Axes and scales
  scale_color_manual(
    values = type_colors,
    name = "Classification",
    drop = FALSE
  ) +
  scale_shape_manual(
    values = site_shapes,
    name = "Site",
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
    x = expression(delta^{13} * C[CH[4]] ~ "(\u2030)"),
    y = expression(delta^{13} * C[CO[2]] ~ "(\u2030)")
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  whiticar_theme

## Show plot
whiticar_panel

## ------------------------------------------------------------
## 7) Save individual panels
## ------------------------------------------------------------
# dir.create("Plots", showWarnings = FALSE)
# 
# ggsave("Plots/CO2_CH4_Emissions.png", CO2_CH4, width = 4.5, height = 3.5, dpi = 300, bg = "transparent")
# ggsave("Plots/isotope_plot_Emissions.png", isotope_plot, width = 4.5, height = 3.5, dpi = 300, bg = "transparent")
# ggsave("Plots/alpha_plot_Emissions.png", alpha_plot, width = 4.5, height = 3.0, dpi = 300, bg = "transparent")
# ggsave("Plots/Figure4D_Whiticar_panel.png", whiticar_panel, width = 7.4, height = 4.6, units = "in", bg = "white")

## ------------------------------------------------------------
## 8) Final multi-panel figure with A–D
## ------------------------------------------------------------
top_row <- CO2_CH4 | isotope_plot

bottom_row <- alpha_plot + whiticar_panel +
  plot_layout(widths = c(1, 2.5))

overall <- (top_row / bottom_row) +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right")

overall

ggsave(
  "Plots/Fig4_Emissions_with_Whiticar_2.png",
  overall,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "transparent"
)
