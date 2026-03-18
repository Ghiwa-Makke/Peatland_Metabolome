suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(cowplot)
  library(patchwork)
})

# =========================================================
# 1) READ METADATA
# =========================================================

data <- read.csv("Data/Sample_Metadata.csv")

# one row per site
sites_all <- data %>%
  select(Site, Type, latitude, longitude) %>%
  distinct() %>%
  mutate(
    Site = factor(Site, levels = c("T1", "BLF", "S3", "US2", "LC"))
  ) %>%
  arrange(Site)

# points to show on main map only
sites_map <- sites_all %>%
  filter(Site %in% c("T1", "US2", "LC")) %>%
  mutate(
    Label = c("MEF", "US2", "LC"),
    PlotType = c("MEF", as.character(Type[Site == "US2"]), as.character(Type[Site == "LC"]))
  )

# safer explicit version for PlotType
sites_map <- sites_map %>%
  mutate(
    PlotType = case_when(
      Label == "MEF" ~ "MEF",
      Label == "US2" ~ "Fen",
      Label == "LC"  ~ "Shrub Wetland"
    )
  )

# =========================================================
# 2) COLORS
# =========================================================

point_colors <- c(
  "MEF" = "#8B4513",
  "Fen" = "#F5DEB3",
  "Shrub Wetland" = "#2F8B8B"
)

# =========================================================
# 3) IMAGE PATHS
# =========================================================

# panel B: uploaded MEF figure
mef_figure_path <- "Data/Pictures/MEF_marked.png"

# panel C: local photos for all 5 sites
site_photos <- tibble(
  Site = c("T1", "BLF", "S3", "US2", "LC"),
  Title = c("T1", "BLF", "S3", "US2", "LC"),
  image = c(
    "Data/Pictures/T1_cropped.jpg",
    "Data/Pictures/BLF_Cropped.png",
    "Data/Pictures/S3_Cropped.jpg",
    "Data/Pictures/US2_Cropped.jpg",
    "Data/Pictures/LC_cropped.jpg"
  )
)

# optional check
print(file.exists(site_photos$image))
print(file.exists(mef_figure_path))

# =========================================================
# 4) MAIN MAP (PANEL A)
# =========================================================

states <- ne_states(
  country = "United States of America",
  returnclass = "sf"
) %>%
  filter(name %in% c("Minnesota", "Wisconsin"))

sites_sf <- st_as_sf(
  sites_map,
  coords = c("longitude", "latitude"),
  crs = 4326
)

coords_df <- st_coordinates(sites_sf) %>%
  as_tibble() %>%
  bind_cols(sites_map)

pA_map <- ggplot() +
  geom_sf(
    data = states,
    fill = "grey96",
    color = "grey35",
    linewidth = 0.4
  ) +
  geom_sf(
    data = sites_sf,
    aes(fill = PlotType),
    shape = 21,
    color = "black",
    size = 4.5,
    stroke = 0.8
  ) +
  geom_text(
    data = coords_df,
    aes(x = X, y = Y, label = Label),
    nudge_y = 0.22,
    fontface = "bold",
    size = 4.2
  ) +
  scale_fill_manual(values = point_colors) +
  coord_sf(
    xlim = c(-97.8, -86.8),
    ylim = c(43.5, 49.8),
    expand = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    #plot.margin = margin(0, 8, 8, 8)
  )

# =========================================================
# 5) PANEL B: MEF FIGURE
# =========================================================

pB_mef <- ggdraw() +
  draw_image(mef_figure_path, x = 0, y = 0, width = 1, height = 1)

# =========================================================
# 6) PANEL C: FIVE SITE PHOTOS IN ONE ROW
# =========================================================

make_photo_panel <- function(img_path, title_text) {
  ggdraw() +
    draw_image(
      img_path,
      x = 0.05, y = 0.18,
      width = 0.90, height = 0.72,
      hjust = 0, vjust = 0
    ) +
    draw_label(
      title_text,
      x = 0.5, y = 0.06,
      fontface = "bold",
      size = 11
    )
}

photo_panels <- lapply(seq_len(nrow(site_photos)), function(i) {
  make_photo_panel(site_photos$image[i], site_photos$Title[i])
})

pC_row <- wrap_plots(photo_panels, nrow = 1)
pC_row_labeled <- ggdraw(pC_row) +
  draw_label(
    "C",
    x = 0.01, y = 0.98,
    hjust = 0, vjust = 1,
    fontface = "bold",
    size = 16
  )
# =========================================================
# 7) COMBINE PANELS
# =========================================================

pA_map_labeled <- pA_map +
  annotate("text", x = -97.2, y = 49.5, label = "a", fontface = "bold", size = 6)

pB_mef_labeled <- ggdraw(pB_mef) +
  draw_label(
    "b",
    x = 0.03, y = 0.97,
    hjust = 0, vjust = 1,
    fontface = "bold",
    size = 16
  )

pC_row_labeled <- ggdraw(pC_row) +
  draw_label(
    "c",
    x = 0.01, y = 0.98,
    hjust = 0, vjust = 1,
    fontface = "bold",
    size = 16
  )

top_row <- pA_map_labeled + pB_mef_labeled +
  plot_layout(ncol = 2, widths = c(1.15, 1))

final_plot <- top_row / pC_row_labeled +
  plot_layout(heights = c(1, 1))

final_plot

# =========================================================
# 8) SAVE
# =========================================================

ggsave(
  "Plots/Figure1_site_map_MEF_and_photos.png",
  final_plot,
  width = 6,
  height = 4,
  dpi = 600,
  bg = "white"
)
