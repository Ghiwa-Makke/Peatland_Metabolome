#### Northern Peatland - Metabolomics Paper ------
### Author: Ghiwa Makke
### Figure 1 - Sites 

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(cowplot)
  library(patchwork)
  library(grid)
  library(magick)   # for MEF image aspect ratio
})

# Setting font sizes -------

BASE_FONT_PT  <- 9    # site labels, axis text, scale bar text
TAG_FONT_PT   <- 12   # bold panel tags (a, b, c)
TITLE_FONT_PT <- 10   # site photo titles

pt_to_mm <- function(pt) pt / .pt   # ggplot2 geom_text/annotate 'size' is in mm

# 1) Read data --------
data <- read.csv("Data/Sample_Metadata.csv")

sites_all <- data %>%
  select(Site, Type, latitude, longitude) %>%
  distinct() %>%
  mutate(Site = factor(Site, levels = c("T1", "BLF", "S3", "US2", "LC"))) %>%
  arrange(Site)

sites_map <- sites_all %>%
  filter(Site %in% c("T1", "US2", "LC")) %>%
  mutate(Label = c("MEF", "US2", "LC")) %>%
  mutate(
    PlotType = case_when(
      Label == "MEF" ~ "MEF",
      Label == "US2" ~ "Fen",
      Label == "LC"  ~ "Shrub Wetland"
    )
  )

# 2) set colors --------
point_colors <- c(
  "MEF" = "#8B4513",
  "Fen" = "#F5DEB3",
  "Shrub Wetland" = "#2F8B8B"
)

# 3) picture paths -----------
mef_figure_path <- "Data/Pictures/MEF_marked.png"

site_photos <- tibble(
  Site  = c("T1", "BLF", "US2", "S3", "LC"),
  Title = c("S1-bog", "BLF-poor fen", "US2-fen", "S3-fen", "LC-shrub wetland"),
  image = c(
    "Data/Pictures/T1_cropped.jpg",
    "Data/Pictures/BLF_Cropped.png",
    "Data/Pictures/US2_Cropped.jpg",
    "Data/Pictures/S3_Cropped.jpg",
    "Data/Pictures/LC_cropped.jpg"
  )
)

print(file.exists(site_photos$image))
print(file.exists(mef_figure_path))

# 4) Main map - panel A  -------
 states <- ne_states(
  country = "United States of America",
  returnclass = "sf"
) %>%
  filter(name %in% c("Minnesota", "Wisconsin"))

sites_sf <- st_as_sf(sites_map, coords = c("longitude", "latitude"), crs = 4326)

coords_df <- st_coordinates(sites_sf) %>%
  as_tibble() %>%
  bind_cols(sites_map)

# ---- scale bar  ----
scale_lat        <- 43.9
scale_lon_start  <- -97.2
scale_length_km  <- 200

lon_distance_km <- function(delta_lon, lat, lon0 = 0) {
  p1 <- st_sfc(st_point(c(lon0, lat)), crs = 4326)
  p2 <- st_sfc(st_point(c(lon0 + delta_lon, lat)), crs = 4326)
  as.numeric(st_distance(p1, p2)) / 1000
}

delta_lon <- uniroot(
  f = function(d) lon_distance_km(d, lat = scale_lat) - scale_length_km,
  interval = c(0, 10)
)$root

scale_lon_end  <- scale_lon_start + delta_lon
scale_lon_mid  <- scale_lon_start + delta_lon / 2
scale_bar_y    <- scale_lat
scale_tick_y0  <- scale_lat - 0.05
scale_tick_y1  <- scale_lat + 0.05
scale_label_y  <- scale_lat + 0.12

pA_map <- ggplot() +
  geom_sf(data = states, fill = "grey96", color = "grey35", linewidth = 0.4) +
  geom_sf(
    data = sites_sf, aes(fill = PlotType),
    shape = 21, color = "black", size = 3.5, stroke = 0.8
  ) +
  geom_text(
    data = coords_df,
    aes(x = X, y = Y, label = Label),
    nudge_y = 0.4, nudge_x = 0.6,
    fontface = "bold",
    size = pt_to_mm(BASE_FONT_PT)
  ) +
  annotate("segment", x = scale_lon_start, xend = scale_lon_end,
           y = scale_bar_y, yend = scale_bar_y, linewidth = 0.5, color = "black") +
  annotate("segment", x = scale_lon_start, xend = scale_lon_start,
           y = scale_tick_y0, yend = scale_tick_y1, linewidth = 1.0, color = "black") +
  annotate("segment", x = scale_lon_mid, xend = scale_lon_mid,
           y = scale_tick_y0, yend = scale_tick_y1, linewidth = 1.0, color = "black") +
  annotate("segment", x = scale_lon_end, xend = scale_lon_end,
           y = scale_tick_y0, yend = scale_tick_y1, linewidth = 1.0, color = "black") +
  annotate("text", x = scale_lon_start, y = scale_label_y, label = "0",
           size = pt_to_mm(BASE_FONT_PT), vjust = -0.3) +
  annotate("text", x = scale_lon_mid, y = scale_label_y,
           label = paste0(scale_length_km / 2),
           size = pt_to_mm(BASE_FONT_PT), vjust = -0.3, hjust = 0.5) +
  annotate("text", x = scale_lon_end, y = scale_label_y,
           label = paste0(scale_length_km, " km"),
           size = pt_to_mm(BASE_FONT_PT), hjust = 0.1, vjust = -0.3) +
  scale_fill_manual(values = point_colors) +
  coord_sf(
    xlim = c(-97.8, -86.8),
    ylim = c(43.5, 49.8),
    expand = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = BASE_FONT_PT) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text  = element_text(color = "black", size = BASE_FONT_PT),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1)     # trimmed margin, reduces whitespace
  )

mean_lat_rad <- mean(c(43.5, 49.8)) * pi / 180
lon_range    <- 86.8 - 97.8  # negative sign doesn't matter, using abs below
lat_range    <- 49.8 - 43.5
map_aspect   <- (abs(lon_range) * cos(mean_lat_rad)) / lat_range


# 5) Panel B - MEF map -------------
mef_info   <- image_info(image_read(mef_figure_path))
mef_aspect <- mef_info$width / mef_info$height

pB_mef <- ggdraw() +
  draw_image(mef_figure_path, x = 0, y = 0, width = 1, height = 1)

# 6) Panel C: site pictures -------
make_photo_panel <- function(img_path, title_text) {
  ggdraw() +
    draw_image(
      img_path,
      x = 0.02, y = 0.14,      # tighter crop margins = less whitespace
      width = 0.96, height = 0.80,
      hjust = 0, vjust = 0
    ) +
    draw_label(
      title_text,
      x = 0.5, y = 0.05,
      fontface = "bold",
      size = TITLE_FONT_PT
    )
}

photo_panels <- lapply(seq_len(nrow(site_photos)), function(i) {
  make_photo_panel(site_photos$image[i], site_photos$Title[i])
})

pC_row <- wrap_plots(photo_panels, nrow = 1) &
  theme(plot.margin = margin(1, 1, 1, 1))

# 7) Combining panels --------------
pA_map_labeled <- ggdraw(pA_map) +
  draw_label("a", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = TAG_FONT_PT)

pB_mef_labeled <- ggdraw(pB_mef) +
  draw_label("b", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = TAG_FONT_PT)

pC_row_labeled <- ggdraw(pC_row) +
  draw_label("c", x = 0.005, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = TAG_FONT_PT)



top_row <- pA_map_labeled + pB_mef_labeled +
  plot_layout(ncol = 2, widths = c(map_aspect, mef_aspect))

final_plot <- (top_row / pC_row_labeled) +
  plot_layout(heights = c(1, 0.75)) & 
  theme(plot.margin = margin(1, 1, 1, 1))

final_plot


# save -----

ggsave(
  "Plots/New/Figure1_site_map_MEF_and_photos_.png",
  final_plot,
  width  = 174 / 25.4,   # 174 mm Springer single-column width
  height = 120 / 25.4,   # reduced from 4.5in — less dead vertical space
  dpi    = 600,
  bg     = "white"
)
