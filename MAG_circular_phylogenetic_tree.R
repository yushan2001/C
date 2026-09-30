# MAG circular phylogenetic tree with genome-quality annotations
#
# Annotation tracks:
#   Completeness, contamination, GC content, species-level assignment,
#   genome size, abundance, and phylum.
#
# Original paths, plotting parameters, color scales, dimensions,
# and PNG/PDF output settings are retained.
# ======================================================================

# Install and load necessary packages
if (!requireNamespace("ggtree", quietly = TRUE)) {
  install.packages("ggtree")
}
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}

library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(readxl)
library(viridis)
library(cowplot)
library(ggnewscale)

# ===================== Plot and value-range parameters =====================
# Set numeric limits manually if needed; NA uses the observed Excel range.
completeness_min <- NA  # 完整性最小值（Excel minimum）
completeness_max <- NA  # 完整性最大值（Excel maximum）
contamination_min <- NA # 污染度最小值（Excel minimum）
contamination_max <- NA # 污染度最大值（Excel maximum）
gc_min <- NA            # GC含量最小值（Excel minimum）
gc_max <- NA            # GC含量最大值（Excel maximum）
size_min <- NA          # 大小最小值（Excel minimum）
size_max <- NA          # 大小最大值（Excel maximum）
abundance_min <- NA     # 丰度最小值（Excel minimum）
abundance_max <- NA     # 丰度最大值（Excel maximum）

# Plot layout parameters
strip_offset <- 0.1
strip_width <- 0.3
tree_line_size <- 0.8
tip_point_size <- 2.5
hilight_alpha <- 0.3
hilight_extend <- 0.05

# ================================================================================

# Read Excel data and extract target bins
checkm_path <- "C:/Users/余山小可爱/Desktop/环形树数据.xlsx"
checkm_df <- read_excel(checkm_path) %>%
  rename(Bin = `Bin Id`) %>%
  mutate(Bin = trimws(Bin))
excel_bins <- checkm_df$Bin
cat("Excel中提取到", length(excel_bins), "个目标bin\n")

# Read and filter the phylogenetic tree
tree_path <- "C:/Users/余山小可爱/Desktop/gtdb_results.backbone.bac120.classify.tree"
tree <- read.tree(tree_path)
keep_index <- which(tree$tip.label %in% excel_bins)
filtered_tree <- drop.tip(tree, tree$tip.label[-keep_index])
cat("进化树过滤后保留", length(filtered_tree$tip.label), "个目标bin\n")

# Retain bins present in the filtered tree
checkm_df <- checkm_df %>%
  filter(Bin %in% filtered_tree$tip.label)

# 3. Determine numeric ranges
# 完整性范围
if (is.na(completeness_min)) completeness_min <- min(checkm_df$`Completeness（%）`, na.rm = TRUE)
if (is.na(completeness_max)) completeness_max <- max(checkm_df$`Completeness（%）`, na.rm = TRUE)
# 污染度范围
if (is.na(contamination_min)) contamination_min <- min(checkm_df$`Contamination（%）`, na.rm = TRUE)
if (is.na(contamination_max)) contamination_max <- max(checkm_df$`Contamination（%）`, na.rm = TRUE)
# GC含量范围
if (is.na(gc_min)) gc_min <- min(checkm_df$`GC含量（%）`, na.rm = TRUE)
if (is.na(gc_max)) gc_max <- max(checkm_df$`GC含量（%）`, na.rm = TRUE)
# 大小范围
if (is.na(size_min)) size_min <- min(checkm_df$`大小（Mbp）`, na.rm = TRUE)
if (is.na(size_max)) size_max <- max(checkm_df$`大小（Mbp）`, na.rm = TRUE)
# 丰度范围
if (is.na(abundance_min)) abundance_min <- min(checkm_df$丰度, na.rm = TRUE)
if (is.na(abundance_max)) abundance_max <- max(checkm_df$丰度, na.rm = TRUE)

# Print numeric ranges
cat("\nNumeric ranges：\n")
cat("完整性：", completeness_min, "-", completeness_max, "%\n")
cat("污染度：", contamination_min, "-", contamination_max, "%\n")
cat("GC含量：", gc_min, "-", gc_max, "%\n")
cat("大小：", size_min, "-", size_max, "Mbp\n")
cat("丰度：", abundance_min, "-", abundance_max, "\n")

# 4. Data preparation
checkm_df <- checkm_df %>%
  mutate(
    # 完整性：填充缺失值为最小值，限制在适配范围
    `Completeness（%）` = ifelse(is.na(`Completeness（%）`), completeness_min, `Completeness（%）`),
    `Completeness（%）` = pmin(pmax(`Completeness（%）`, completeness_min), completeness_max),
    
    # 污染度：填充缺失值为最小值，限制在适配范围
    `Contamination（%）` = ifelse(is.na(`Contamination（%）`), contamination_min, `Contamination（%）`),
    `Contamination（%）` = pmin(pmax(`Contamination（%）`, contamination_min), contamination_max),
    
    # GC含量：填充缺失值为最小值，限制在适配范围
    `GC含量（%）` = ifelse(is.na(`GC含量（%）`), gc_min, `GC含量（%）`),
    `GC含量（%）` = pmin(pmax(`GC含量（%）`, gc_min), gc_max),
    
    # 种水平：填充缺失值为NO
    `种水平` = ifelse(is.na(`种水平`), "NO", `种水平`),
    
    # 大小：填充缺失值为最小值，限制在适配范围
    `大小（Mbp）` = ifelse(is.na(`大小（Mbp）`), size_min, `大小（Mbp）`),
    `大小（Mbp）` = pmin(pmax(`大小（Mbp）`, size_min), size_max),
    
    # 丰度：填充缺失值为最小值，限制在适配范围
    丰度 = ifelse(is.na(丰度), abundance_min, 丰度),
    丰度 = pmin(pmax(丰度, abundance_min), abundance_max)
  )

# 5. Taxonomic information
taxonomy_for_tree <- checkm_df %>%
  mutate(label = Bin) %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum))

# 6. Phylum color mapping
unique_phyla <- unique(taxonomy_for_tree$Phylum)
phylum_colors <- scales::hue_pal()(length(unique_phyla))
phylum_color_map <- setNames(phylum_colors, unique_phyla)

# 7. Circular phylogenetic tree
circular_plot_base <- ggtree(filtered_tree, layout = "circular", color = "black", size = tree_line_size)

# Merge taxonomy data
tree_data <- circular_plot_base$data %>%
  left_join(taxonomy_for_tree %>% select(label, Phylum), by = "label") %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum))
circular_plot_base$data <- tree_data

# Add phylum blocks
for (phylum in unique_phyla) {
  tips_in_phylum <- taxonomy_for_tree %>%
    filter(Phylum == phylum) %>%
    pull(label)
  
  valid_tips <- tips_in_phylum[tips_in_phylum %in% filtered_tree$tip.label]
  
  if (length(valid_tips) >= 1) {
    if (length(valid_tips) == 1) {
      mrca_node <- which(filtered_tree$tip.label == valid_tips)
    } else {
      mrca_node <- getMRCA(filtered_tree, tip = valid_tips)
    }
    
    if (!is.null(mrca_node) && mrca_node > 0) {
      circular_plot_base <- circular_plot_base +
        geom_hilight(
          node = mrca_node,
          fill = phylum_color_map[phylum],
          alpha = hilight_alpha,
          extend = hilight_extend
        )
    }
  }
}

# 8. Calculate maximum radius
max_radius <- max(circular_plot_base$data$x)

# 9. Prepare annotation-ring data
strip_data <- circular_plot_base$data %>%
  filter(isTip) %>%
  left_join(checkm_df, by = c("label" = "Bin")) %>%
  mutate(
    comp_xmin = max_radius + strip_offset,
    comp_xmax = max_radius + strip_offset + strip_width,
    cont_xmin = max_radius + strip_offset*2 + strip_width,
    cont_xmax = max_radius + strip_offset*2 + strip_width*2,
    gc_xmin = max_radius + strip_offset*3 + strip_width*2,
    gc_xmax = max_radius + strip_offset*3 + strip_width*3,
    species_xmin = max_radius + strip_offset*4 + strip_width*3,
    species_xmax = max_radius + strip_offset*4 + strip_width*4,
    size_xmin = max_radius + strip_offset*5 + strip_width*4,
    size_xmax = max_radius + strip_offset*5 + strip_width*5,
    abundance_xmin = max_radius + strip_offset*6 + strip_width*5,
    abundance_xmax = max_radius + strip_offset*6 + strip_width*6
  )

# 10. Plot with legends
plot_with_legends <- circular_plot_base +
  # 10.1 Completeness
  geom_rect(data = strip_data, aes(xmin = comp_xmin + 0.3, xmax = comp_xmax - 0.3, 
                                  ymin = y - 0.5, ymax = y + 0.5, fill = `Completeness（%）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("Completeness (%) [", round(completeness_min,1), "-", round(completeness_max,1), "]"), 
                       colours = c("#E6E6FA", "#C8A2C8", "#9966CC", "#8A2BE2", "#4B0082"), 
                       limits = c(completeness_min, completeness_max),
                       guide = guide_colorbar(order = 1, title.position = "top", title.hjust = 0.5)) +
  new_scale_fill() +
  # 10.2 Contamination
  geom_rect(data = strip_data, aes(xmin = cont_xmin + 0.3, xmax = cont_xmax - 0.3, 
                                  ymin = y - 0.5, ymax = y + 0.5, fill = `Contamination（%）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("Contamination (%) [", round(contamination_min,1), "-", round(contamination_max,1), "]"), 
                       colours = c("#FFE4E1", "#FFB6C1", "#FF69B4", "#DC143C", "#8B0000"), 
                       limits = c(contamination_min, contamination_max),
                       guide = guide_colorbar(order = 2, title.position = "top", title.hjust = 0.5)) +
  new_scale_fill() +
  # 10.3 GC Content
  geom_rect(data = strip_data, aes(xmin = gc_xmin + 0.3, xmax = gc_xmax - 0.3, 
                                  ymin = y - 0.5, ymax = y + 0.5, fill = `GC含量（%）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("GC Content (%) [", round(gc_min,1), "-", round(gc_max,1), "]"), 
                       colours = c("#E6F7FF", "#BAE7FF", "#91D5FF", "#1890FF", "#0050B3"), 
                       limits = c(gc_min, gc_max),
                       guide = guide_colorbar(order = 3, title.position = "top", title.hjust = 0.5)) +
  new_scale_fill() +
  # 10.4 Species Level
  geom_point(data = strip_data %>% filter(!is.na(`种水平`)), 
             aes(x = species_xmin + (species_xmax - species_xmin)/2, y = y, fill = `种水平`), 
             shape = 21, size = 3, color = "white", stroke = 0.2, inherit.aes = FALSE) +
  scale_fill_manual(name = "Species Level", 
                    values = c("NO" = "gray70", "YES" = "red2"), 
                    guide = guide_legend(order = 4, title.position = "top", title.hjust = 0.5)) +
  new_scale_fill() +
  # 10.5 Size
  geom_rect(data = strip_data, aes(xmin = size_xmin, 
                                  xmax = size_xmin + ((`大小（Mbp）` - size_min)/(size_max - size_min)) * (size_xmax - size_xmin), 
                                  ymin = y - 0.4, ymax = y + 0.4, fill = `大小（Mbp）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("Size (Mbp) [", round(size_min,1), "-", round(size_max,1), "]"), 
                       colours = c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"), 
                       limits = c(size_min, size_max),
                       guide = guide_colorbar(order = 5, title.position = "top", title.hjust = 0.5)) +
  new_scale_fill() +
  # 10.6 Abundance
  geom_rect(data = strip_data, aes(xmin = abundance_xmin, 
                                  xmax = abundance_xmin + ((丰度 - abundance_min)/(abundance_max - abundance_min)) * (abundance_xmax - abundance_xmin), 
                                  ymin = y - 0.4, ymax = y + 0.4, fill = 丰度), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("Abundance [", round(abundance_min,1), "-", round(abundance_max,1), "]"), 
                       colours = c("#FFFDE7", "#FFF9C4", "#FFF59D", "#FFEE58", "#FFD600", "#FFC107", "#FFA000", "#FF8F00", "#FF6F00"), 
                       limits = c(abundance_min, abundance_max),
                       guide = guide_colorbar(order = 6, title.position = "top", title.hjust = 0.5)) +
  # 10.7 Phylum
  geom_tippoint(aes(color = Phylum), size = tip_point_size, alpha = 0.9) +
  scale_color_manual(values = phylum_color_map, 
                     name = "Phylum", 
                     guide = guide_legend(order = 7, title.position = "top", title.hjust = 0.5)) +
  # Theme
  theme_tree2() +
  ggtitle("Circular Phylogenetic Tree with Six Annotation Rings") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  )

# 11. Extract and save individual legends
legend <- get_legend(plot_with_legends)
legend_components <- list(
  Completeness = legend$grobs[[1]],
  Contamination = legend$grobs[[2]],
  GC_Content = legend$grobs[[3]],
  Species_Level = legend$grobs[[4]],
  Size = legend$grobs[[5]],
  Abundance = legend$grobs[[6]],
  Phylum = legend$grobs[[7]]
)

# 保存图例
output_dir <- "C:/Users/余山小可爱/Desktop/legends/"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (name in names(legend_components)) {
  single_legend <- ggplot() +
    annotation_custom(grob = legend_components[[name]]) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  # Save PNG
  filename_png <- paste0(output_dir, name, "_legend.png")
  ggsave(
    filename = filename_png,
    plot = single_legend,
    width = 5,
    height = 3,
    dpi = 600,
    bg = "transparent",
    device = png(type = "cairo")
  )
  
  # Save PDF
  filename_pdf <- paste0(output_dir, name, "_legend.pdf")
  ggsave(
    filename = filename_pdf,
    plot = single_legend,
    width = 5,
    height = 3,
    bg = "transparent",
    device = "pdf"
  )
  
  cat("已保存", name, "图例: ", filename_png, " + ", filename_pdf, "\n")
}

# 12. Main plot without legends
main_plot <- circular_plot_base +
  # 12.1 Completeness
  geom_rect(data = strip_data, aes(xmin = comp_xmin + 0.3, xmax = comp_xmax - 0.3, 
                                  ymin = y - 0.5, ymax = y + 0.5, fill = `Completeness（%）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = c("#E6E6FA", "#C8A2C8", "#9966CC", "#8A2BE2", "#4B0082"), 
                       limits = c(completeness_min, completeness_max), guide = "none") +
  new_scale_fill() +
  # 12.2 Contamination
  geom_rect(data = strip_data, aes(xmin = cont_xmin + 0.3, xmax = cont_xmax - 0.3, 
                                  ymin = y - 0.5, ymax = y + 0.5, fill = `Contamination（%）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = c("#FFE4E1", "#FFB6C1", "#FF69B4", "#DC143C", "#8B0000"), 
                       limits = c(contamination_min, contamination_max), guide = "none") +
  new_scale_fill() +
  # 12.3 GC Content
  geom_rect(data = strip_data, aes(xmin = gc_xmin + 0.3, xmax = gc_xmax - 0.3, 
                                  ymin = y - 0.5, ymax = y + 0.5, fill = `GC含量（%）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = c("#E6F7FF", "#BAE7FF", "#91D5FF", "#1890FF", "#0050B3"), 
                       limits = c(gc_min, gc_max), guide = "none") +
  new_scale_fill() +
  # 12.4 Species Level
  geom_point(data = strip_data %>% filter(!is.na(`种水平`)), 
             aes(x = species_xmin + (species_xmax - species_xmin)/2, y = y, fill = `种水平`), 
             shape = 21, size = 3, color = "white", stroke = 0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = c("NO" = "gray70", "YES" = "red2"), guide = "none") +
  new_scale_fill() +
  # 12.5 Size
  geom_rect(data = strip_data, aes(xmin = size_xmin, 
                                  xmax = size_xmin + ((`大小（Mbp）` - size_min)/(size_max - size_min)) * (size_xmax - size_xmin), 
                                  ymin = y - 0.4, ymax = y + 0.4, fill = `大小（Mbp）`), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"), 
                       limits = c(size_min, size_max), guide = "none") +
  new_scale_fill() +
  # 12.6 Abundance
  geom_rect(data = strip_data, aes(xmin = abundance_xmin, 
                                  xmax = abundance_xmin + ((丰度 - abundance_min)/(abundance_max - abundance_min)) * (abundance_xmax - abundance_xmin), 
                                  ymin = y - 0.4, ymax = y + 0.4, fill = 丰度), 
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = c("#FFFDE7", "#FFF9C4", "#FFF59D", "#FFEE58", "#FFD600", "#FFC107", "#FFA000", "#FF8F00", "#FF6F00"), 
                       limits = c(abundance_min, abundance_max), guide = "none") +
  # 12.7 Phylum
  geom_tippoint(aes(color = Phylum), size = tip_point_size, alpha = 0.9) +
  scale_color_manual(values = phylum_color_map, guide = "none") +
  # Theme
  theme_tree2() +
  ggtitle("Circular Phylogenetic Tree with Six Annotation Rings") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# 13. Save main plot
# Save PNG
ggsave(
  filename = "C:/Users/余山小可爱/Desktop/main_plot.png",
  plot = main_plot,
  width = 22,
  height = 22,
  dpi = 600,
  bg = "transparent",
  device = png(type = "cairo")
)

# Save PDF
ggsave(
  filename = "C:/Users/余山小可爱/Desktop/main_plot.pdf",
  plot = main_plot,
  width = 22,
  height = 22,
  bg = "transparent",
  device = "pdf"
)

cat("主图已保存至: C:/Users/余山小可爱/Desktop/main_plot.png + main_plot.pdf\n")
cat("\n All PNG and PDF files have been saved.\n")
