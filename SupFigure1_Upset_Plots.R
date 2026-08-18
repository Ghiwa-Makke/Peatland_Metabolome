suppressPackageStartupMessages({
  library(tidyverse)
  library(UpSetR)
  library(VennDiagram)
  library(RColorBrewer)
  library(grid)
})


site_levels <- c("S1", "BLF", "S3", "US2", "LC")

site_labels_venn <- c(
  "S1"  = "S1-bog",
  "BLF" = "BLF-poor\nfen",
  "S3"  = "S3-fen",
  "US2" = "US2-fen",
  "LC"  = "LC-shrub\nwetland"
)

site_colors <- c(
  "S1"  = "#8B4513",
  "BLF" = "#E6A15C",
  "S3"  = "#4C9A6B",
  "US2" = "#8B6FB3",
  "LC"  = "#2F8B8B"
)

upset_set_names <- c(
  "S1_bog"      = "S1",
  "BLF_poorfen" = "BLF",
  "S3_fen"      = "S3",
  "US2_fen"     = "US2",
  "LC_shrub"    = "LC"
)

add_alpha <- function(col, alpha = 0.6) grDevices::adjustcolor(col, alpha.f = alpha)


save_upset_svg <- function(compounds_df, mode_label, outfile,
                           width = 7, height = 4.5) {
  
  site_list <- lapply(names(upset_set_names), function(nm) {
    unique(compounds_df$FeatureID[compounds_df$Site == upset_set_names[[nm]]])
  })
  names(site_list) <- names(upset_set_names)
  
  site_queries <- lapply(names(upset_set_names), function(nm) {
    list(
      query  = intersects,
      params = list(nm),
      color  = unname(site_colors[upset_set_names[[nm]]]),
      active = TRUE
    )
  })
  
  full_intersection_query <- list(
    query  = intersects,
    params = as.list(names(upset_set_names)),
    color  = "#903495",
    active = TRUE
  )
  
  svg(outfile, width = width, height = height)
  print(
    upset(
      fromList(site_list),
      order.by = "freq",
      mainbar.y.label = paste0("Intersection Size \u2013 ", mode_label),
      sets.x.label = "Set Size",
      text.scale = c(1.3, 1.1, 1.3, 1.1, 1.2, 1.0),
      point.size = 2.6,
      line.size  = 0.9,
      queries = c(site_queries, list(full_intersection_query))
    )
  )
  dev.off()
  
  message("Saved: ", outfile)
}

save_venn_svg <- function(compounds_df, mode_label, outfile,
                          width = 5, height = 5) {
  
  sets <- lapply(site_levels, function(s) {
    unique(compounds_df$FeatureID[compounds_df$Site == s])
  })
  names(sets) <- site_levels
  
  comb_size <- function(idx) length(Reduce(intersect, sets[idx]))
  
  area1 <- comb_size(1); area2 <- comb_size(2); area3 <- comb_size(3)
  area4 <- comb_size(4); area5 <- comb_size(5)
  
  n12 <- comb_size(c(1,2)); n13 <- comb_size(c(1,3)); n14 <- comb_size(c(1,4)); n15 <- comb_size(c(1,5))
  n23 <- comb_size(c(2,3)); n24 <- comb_size(c(2,4)); n25 <- comb_size(c(2,5))
  n34 <- comb_size(c(3,4)); n35 <- comb_size(c(3,5)); n45 <- comb_size(c(4,5))
  
  n123 <- comb_size(c(1,2,3)); n124 <- comb_size(c(1,2,4)); n125 <- comb_size(c(1,2,5))
  n134 <- comb_size(c(1,3,4)); n135 <- comb_size(c(1,3,5)); n145 <- comb_size(c(1,4,5))
  n234 <- comb_size(c(2,3,4)); n235 <- comb_size(c(2,3,5)); n245 <- comb_size(c(2,4,5))
  n345 <- comb_size(c(3,4,5))
  
  n1234 <- comb_size(c(1,2,3,4)); n1235 <- comb_size(c(1,2,3,5))
  n1245 <- comb_size(c(1,2,4,5)); n1345 <- comb_size(c(1,3,4,5)); n2345 <- comb_size(c(2,3,4,5))
  
  n12345 <- comb_size(1:5)
  
  fills    <- add_alpha(unname(site_colors[site_levels]), alpha = 0.6)
  cat_cols <- unname(site_colors[site_levels])
  
  svg(outfile, width = width, height = height)
  grid.newpage()
  draw.quintuple.venn(
    area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5,
    n12 = n12, n13 = n13, n14 = n14, n15 = n15,
    n23 = n23, n24 = n24, n25 = n25,
    n34 = n34, n35 = n35, n45 = n45,
    n123 = n123, n124 = n124, n125 = n125,
    n134 = n134, n135 = n135, n145 = n145,
    n234 = n234, n235 = n235, n245 = n245,
    n345 = n345,
    n1234 = n1234, n1235 = n1235, n1245 = n1245, n1345 = n1345, n2345 = n2345,
    n12345 = n12345,
    category   = site_labels_venn[site_levels],
    fill       = fills,
    lty        = "solid", lwd = 1,
    col        = "black",
    cat.col    = cat_cols,
    cat.cex    = 0.9,
    cex        = 0.8,
    cat.pos    = c(-20, 20, 120, -20, 20),
    cat.dist   = c(0.075, 0.075, 0.075, 0.09, 0.11),
    margin     = 0.1,
    scaled     = FALSE
  )
  dev.off()
  
  message("Saved: ", outfile)
}

# DATA -------------

HILIC_compounds <- read.csv("Tables/compounds_table_non_real_gap_filled_HILIC.csv") %>%
  filter(!is.na(Date))

RP_compounds <- read.csv("Tables/compounds_table_non_real_gap_filled_RP.csv") %>%
  filter(!is.na(Date))


# save panels ---------


dir.create("Plots/SVG", showWarnings = FALSE, recursive = TRUE)

save_upset_svg(HILIC_compounds, "HILIC", "Plots/SVG/SupFig1a_upset_HILIC.svg")
save_venn_svg(HILIC_compounds,  "HILIC", "Plots/SVG/SupFig1a_venn_HILIC.svg")

save_upset_svg(RP_compounds, "RP", "Plots/SVG/SupFig1b_upset_RP.svg")
save_venn_svg(RP_compounds,  "RP", "Plots/SVG/SupFig1b_venn_RP.svg")
