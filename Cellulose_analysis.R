# Cellulose degradation analysis
# NMDS/PERMANOVA, standardized heatmap, MAG phylogeny,
# and seasonal order-level cellulase composition
# ======================================================================

# ======================================================================
# Part I: NMDS and PERMANOVA
# ======================================================================
if (!require(readxl)) install.packages("readxl")
if (!require(vegan)) install.packages("vegan")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")

library(readxl)    
library(vegan)     
library(ggplot2)   
library(dplyr)     
library(tidyr)     

# ===================== 核心：基于Sheet5纤维素酶数据做 NMDS =====================
# Input data
file_path <- "C:/Users/余山小可爱/Desktop/纤维素原始代码数据.xlsx"
df <- read_excel(file_path, sheet = "Sheet1")

df_long <- df %>%
  pivot_longer(cols = -cellulase_type, names_to = "Sample", values_to = "Abundance") %>%
  pivot_wider(names_from = cellulase_type, values_from = Abundance)

df_selected <- df_long %>% select(-Sample)
rownames(df_selected) <- df_long$Sample
df_selected[is.na(df_selected)] <- 0
df_transposed <- as.matrix(df_selected)

set.seed(42)
nmds_result <- metaMDS(df_transposed, distance = "bray", k = 2, trymax = 100)

nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_coords$Time <- ifelse(grepl("^04", rownames(nmds_coords)), "April 2025", "October 2025")
group_factor <- factor(nmds_coords$Time)

adonis_result <- adonis2(df_transposed ~ group_factor, 
                         data = nmds_coords, 
                         method = "bray", 
                         permutations = 999)
adonis_r2 <- round(adonis_result$R2[1], 4)
adonis_p <- round(adonis_result$`Pr(>F)`[1], 4)

calculate_hull <- function(data) {
  data[chull(data$NMDS1, data$NMDS2), ]
}
hull_data <- nmds_coords %>%
  group_by(Time) %>%
  do(calculate_hull(.))

nmds_plot <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Time, shape = Time)) +
  geom_polygon(data = hull_data, aes(fill = Time), alpha = 0.2, show.legend = FALSE) +
  geom_point(size = 5, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = c("April 2025" = "#FFB3DD", "October 2025" = "#ACFFFF")) +
  scale_fill_manual(values = c("April 2025" = "#FFB3DD", "October 2025" = "#ACFFFF")) +
  scale_shape_manual(values = c("April 2025" = 16, "October 2025" = 17)) +
  labs(title = "NMDS Analysis of Cellulase TPM (April 2025 vs October 2025)",
       x = "NMDS1", y = "NMDS2", 
       color = "Time Point",
       shape = "Time Point") +
  annotate("text", x = min(nmds_coords$NMDS1)*0.9, y = min(nmds_coords$NMDS2)*0.9,
           label = paste0("Stress = ", round(nmds_result$stress, 4), "\n",
                          "Adonis: R² = ", adonis_r2, ", p = ", adonis_p),
           size = 4, fontface = "bold", hjust = 0,
           bbox = list(boxstyle = "round,pad=0.5", fill = "white", alpha = 0.9)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

# Display plot
print(nmds_plot)

ggsave("C:/Users/余山小可爱/Desktop/纤维素酶_NMDS图.png", 
       plot = nmds_plot, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("C:/Users/余山小可爱/Desktop/纤维素酶_NMDS图.pdf", 
       plot = nmds_plot, width = 10, height = 8, bg = "white")

cat("=== Cellulase TPM NMDS Analysis Results ===\n")
cat("分析数据：Sheet5 所有纤维素酶TPM值\n")
cat("NMDS Stress：", round(nmds_result$stress, 4), "\n")
cat("Adonis R²：", adonis_r2, "\n")
cat("Adonis P-value：", adonis_p, "\n")


# ======================================================================
# Part II: Standardized cellulase heatmap
# ======================================================================
if (!require(readxl)) install.packages("readxl")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(scales)) install.packages("scales")

library(readxl)    
library(ggplot2)   
library(dplyr)     
library(tidyr)     
library(scales)

# Input data
file_path <- "C:/Users/余山小可爱/Desktop/纤维素原始代码数据.xlsx"
df_raw <- read_excel(file_path, sheet = "Sheet1")

# Check data structure
cat("=== 🔍 原始数据核对 ===\n")
cat("数据维度：", nrow(df_raw), " 种纤维素酶 × ", ncol(df_raw)-1, " 个样本\n")
cat("样本列表：", paste(colnames(df_raw)[-1], collapse = ", "), "\n")
cat("纤维素酶原始顺序：", paste(df_raw$cellulase_type, collapse = ", "), "\n")
cat("\n")

# Z-score normalization
# Convert to numeric matrix
df_matrix <- df_raw %>%
  column_to_rownames("cellulase_type") %>%
  as.matrix()

# Replace missing values with zero
df_matrix[is.na(df_matrix)] <- 0

# Row-wise Z-score normalization
# Z = (x - mean) / SD
df_zscore <- t(apply(df_matrix, 1, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

# Convert to long format
df_zscore_long <- as.data.frame(df_zscore) %>%
  rownames_to_column("cellulase_type") %>%
  pivot_longer(
    cols = -cellulase_type,
    names_to = "Sample",
    values_to = "Z_score"
  )

# Preserve sample order
sample_order <- colnames(df_raw)[-1]
df_zscore_long$Sample <- factor(df_zscore_long$Sample, levels = sample_order)

enzyme_original_order <- df_raw$cellulase_type
df_zscore_long$cellulase_type <- factor(df_zscore_long$cellulase_type, levels = enzyme_original_order)

heatmap_plot <- ggplot(
  df_zscore_long,
  aes(x = Sample, y = cellulase_type, fill = Z_score)
) +
  # Heatmap tiles
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-3, 3),  # 限制Z值范围，避免极端值影响配色
    oob = squish,
    name = "Z-score"
  ) +
  # Labels
  labs(
    title = "Standardized Heatmap of Cellulase TPM",
    x = "Sample",
    y = "Cellulase Type"
  ) +
  theme_bw() +
  theme(
    # Title
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    # Axes
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    # Legend
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    # Remove panel grid
    panel.grid = element_blank()
  )

# Display heatmap
print(heatmap_plot)

ggsave(
  "C:/Users/余山小可爱/Desktop/纤维素酶标准化热图_原始顺序.png",
  plot = heatmap_plot,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)
ggsave(
  "C:/Users/余山小可爱/Desktop/纤维素酶标准化热图_原始顺序.pdf",
  plot = heatmap_plot,
  width = 14,
  height = 10,
  bg = "white"
)

# Heatmap summary


# ======================================================================
# Part III: Cellulose-related genes across MAGs
# ======================================================================
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

# Packages
if (!requireNamespace("ggtree", quietly = TRUE)) { install.packages("BiocManager"); BiocManager::install("ggtree") }
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(readxl)
library(viridis)
library(cowplot)
library(ggnewscale)

# Input and output paths
file_path  <- "C:/Users/余山小可爱/Desktop/纤维素原始代码数据.xlsx"
tree_path  <- "C:/Users/余山小可爱/Desktop/gtdb_results.backbone.bac120.classify.tree"
output_dir <- "C:/Users/余山小可爱/Desktop/Output_Phylogeny_Pro/"
if (!dir.exists(output_dir)) dir.create(output_dir)

tree_line_size <- 0.8     # 树枝线条粗细
tip_point_size <- 2.5     # 门分类 Tip 点大小
hilight_alpha  <- 0.3     # 背景高亮透明度

# Ring layout
strip_offset   <- 0.08    # 第一圈热图距离树顶端的起始距离
ring_width     <- 0.045   # 每圈热图的物理宽度
ring_gap       <- 0.006   # 圈与圈之间的空隙
border_size    <- 0.05    # 热图小方块的黑色边框粗细

gene_cols <- c(
  "1,4-beta-D-glucan glucohydrolase",
  "Beta-N-acetylglucosaminidase/beta-glucosidase",
  "Beta-glucosidase",
  "Beta-mannanase/endoglucanase",
  "Bifunctional beta-D-glucosidase/beta-D-fucosidase",
  "Cellodextrinase",
  "Cellulase",
  "Cellulase CelDZ1",
  "Cellulase/esterase CelE",
  "Cellulose 1,4-beta-cellobiosidase",
  "Endoglucanase",
  "Endoglucanase/exoglucanase",
  "Exoglucanase",
  "Exoglucanase/xylanase",
  "Major extracellular endoglucanase",
  "Minor endoglucanase",
  "Periplasmic beta-glucosidase",
  "Periplasmic beta-glucosidase/beta-xylosidase",
  "Probable endoglucanase",
  "Thermostable beta-glucosidase",
  "Thermostable celloxylanase"
)

# Prepare MAG data
# Read data and trim column/Bin names
full_df <- read_excel(file_path, sheet = "Sheet2")
colnames(full_df) <- trimws(colnames(full_df)) 

full_df <- full_df %>%
  rename(Bin = `bin_name`) %>%  
  mutate(Bin = trimws(Bin))

# Match tree tips to bins in the data
tree <- read.tree(tree_path)
keep_bins <- intersect(tree$tip.label, full_df$Bin)
filtered_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep_bins])

checkm_df <- full_df %>% 
  filter(Bin %in% filtered_tree$tip.label) %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum))

# Convert gene counts to numeric values
for(col in gene_cols) {
  if(col %in% colnames(checkm_df)) {
    checkm_df[[col]] <- as.numeric(as.character(checkm_df[[col]]))
    checkm_df[[col]][is.na(checkm_df[[col]])] <- 0
  } else {
    warning(paste0("注意：在 Excel 中未找到列名：'", col, "'，请检查是否字打错了！"))
  }
}

# Maximum gene copy number
max_actual_val <- checkm_df %>% select(any_of(gene_cols)) %>% max(na.rm = TRUE)
if(is.na(max_actual_val) || max_actual_val == 0) max_actual_val <- 10

# Circular phylogeny and phylum sectors
# Phylum color mapping
unique_phyla <- unique(checkm_df$Phylum)
phylum_color_map <- setNames(scales::hue_pal()(length(unique_phyla)), unique_phyla)

p_base <- ggtree(filtered_tree, layout = "circular", color = "black", size = tree_line_size)

# Tree radius
max_radius <- max(p_base$data$x)
tip_coords <- p_base$data %>% filter(isTip) %>% select(label, y)

current_plot <- p_base

# 3.2 动态计算高亮背景延伸长度（基础半径 + 起始偏移 + 21圈的总物理宽度）
heatmap_total_width <- length(gene_cols) * (ring_width + ring_gap)
extend_length <- strip_offset + heatmap_total_width - ring_gap

# 沿用参考脚本的单/多节点门混合兼容高亮逻辑
for (phylum_name in unique_phyla) {
  tips_in_phylum <- checkm_df %>% filter(Phylum == phylum_name) %>% pull(Bin)
  valid_tips <- tips_in_phylum[tips_in_phylum %in% filtered_tree$tip.label]
  
  if (length(valid_tips) >= 1) {
    if (length(valid_tips) == 1) {
      mrca_node <- which(filtered_tree$tip.label == valid_tips)
    } else {
      mrca_node <- getMRCA(filtered_tree, tip = valid_tips)
    }
    
    if (!is.null(mrca_node) && mrca_node > 0) {
      current_plot <- current_plot + 
        geom_hilight(node = mrca_node, fill = phylum_color_map[phylum_name], 
                     alpha = hilight_alpha, extend = extend_length)
    }
  }
}

# Gene copy-number rings
for (i in seq_along(gene_cols)) {
  gene_name <- gene_cols[i]
  if(!gene_name %in% colnames(checkm_df)) next
  
  ring_data <- checkm_df %>% 
    select(Bin, val = !!sym(gene_name)) %>%
    inner_join(tip_coords, by = c("Bin" = "label")) %>%
    # Ring position
    mutate(x = max_radius + strip_offset + (i-1) * (ring_width + ring_gap) + (ring_width/2))
  
  current_plot <- current_plot +
    new_scale_fill() +
    geom_tile(data = ring_data, aes(x = x, y = y, fill = val), 
              width = ring_width, height = 1, color = "black", 
              size = border_size, inherit.aes = FALSE) +
    scale_fill_gradientn(colors = c("white", "#deebf7", "#3182ce", "#084594"),
                         values = scales::rescale(c(0, 0.1, max_actual_val * 0.5, max_actual_val)), 
                         limits = c(0, max_actual_val), oob = scales::squish, guide = "none")
}

# Phylum tip points
final_plot <- current_plot + 
  new_scale_color() + 
  geom_tippoint(data = p_base$data %>% filter(isTip) %>% left_join(checkm_df, by = c("label" = "Bin")), 
                aes(x = x, y = y, color = Phylum), size = tip_point_size, alpha = 0.9, inherit.aes = FALSE) +
  scale_color_manual(values = phylum_color_map, name = "Phylum")

# Export phylogeny and legend
first_valid_gene <- gene_cols[1]

legend_plot <- ggplot(checkm_df) +
  geom_point(aes(x = 1, y = 1, color = Phylum)) +
  scale_color_manual(values = phylum_color_map) +
  new_scale_fill() +
  geom_tile(aes(x = 1, y = 1, fill = !!sym(first_valid_gene))) + 
  scale_fill_gradientn(colors = c("white", "#deebf7", "#3182ce", "#084594"),
                       values = scales::rescale(c(0, 0.1, max_actual_val * 0.5, max_actual_val)), 
                       limits = c(0, max_actual_val), 
                       breaks = seq(0, max_actual_val, by = max(1, round(max_actual_val / 5))),
                       name = "Gene Copy Number\n(基因拷贝绝对数量)") +
  theme_bw()

all_legend <- get_legend(legend_plot)

# Remove main legend and extend x-axis for outer rings
final_main <- final_plot + 
  theme(legend.position = "none") + 
  xlim(0, max_radius + strip_offset + heatmap_total_width + 0.4)

# Save PNG and PDF
ggsave(paste0(output_dir, "Main_Tree_Cellulose_v3.png"), final_main, width = 18, height = 18, dpi = 600, device = png(type = "cairo"))
ggsave(paste0(output_dir, "Main_Tree_Cellulose_v3.pdf"), final_main, width = 18, height = 18, device = "pdf")

# Save legend
ggsave(paste0(output_dir, "Legends_Cellulose_v3.pdf"), ggdraw(all_legend), width = 5, height = 10, device = "pdf")


# ======================================================================
# Part IV: Seasonal order-level cellulase composition
# ======================================================================
# Packages
if (!require(readxl)) install.packages("readxl")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggpubr)) install.packages("ggpubr")
library(readxl)
library(tidyverse)
library(ggpubr)

# Output directory
output_dir <- "C:/Users/余山小可爱/Desktop/纤维素酶家族_季节对比饼图"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cellulase data by bin
df_enz <- read_excel("C:/Users/余山小可爱/Desktop/纤维素原始代码数据.xlsx", sheet = "Sheet3")
# Genome taxonomy and sample TPM
df_tpm <- read_excel("C:/Users/余山小可爱/Desktop/纤维素原始代码数据.xlsx", sheet = "Sheet3")

# Harmonize genome identifiers
df_enz_renamed <- df_enz %>% rename(Genome = bin_name)

df_merge <- df_tpm %>%
  inner_join(df_enz_renamed, by = "Genome") %>%
  filter(!is.na(Order)) # 过滤掉无分类学信息的行

# Cellulase columns
all_enz_cols <- setdiff(colnames(df_enz_renamed), "Genome")

df_merge <- df_merge %>%
  mutate(
    Month04 = rowMeans(select(., starts_with("04")), na.rm = TRUE),
    Month10 = rowMeans(select(., starts_with("10")), na.rm = TRUE)
  )

# Order-level pie-chart function
plot_cellulase_pie <- function(enz_name = "Total_Cellulase", enz_cols = all_enz_cols, topN = 10) {
  # Total target-enzyme count per genome
  df_tmp <- df_merge %>%
    mutate(Total_enz = rowSums(select(., all_of(enz_cols)), na.rm = TRUE))
  
  # Aggregate weighted TPM by order
  df_sum <- df_tmp %>%
    group_by(Order) %>%
    summarise(
      `04月` = sum(Total_enz * Month04, na.rm = TRUE),
      `10月` = sum(Total_enz * Month10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(`04月` > 0 | `10月` > 0) # 过滤掉无丰度的目
  
  # Retain top 10 orders; combine the remainder as Others
  top_orders <- df_sum %>%
    mutate(Total = `04月` + `10月`) %>%
    arrange(desc(Total)) %>%
    slice_head(n = topN) %>%
    pull(Order)
  
  df_plot <- df_sum %>%
    mutate(Order = ifelse(Order %in% top_orders, Order, "Others")) %>%
    group_by(Order) %>%
    summarise(
      `04月` = sum(`04月`),
      `10月` = sum(`10月`),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = -Order, names_to = "Month", values_to = "TPM") %>%
    group_by(Month) %>%
    mutate(Rel_Abund = TPM / sum(TPM) * 100) %>%
    ungroup() %>%
    mutate(Month = factor(Month, levels = c("04月", "10月")))
  
  macaron_colors <- c(
    "#FFB3BA","#FFDFBA","#FFFFBA","#BAFFC9","#BAE1FF",
    "#D6BAFF","#FFBAF0","#BAB0FF","#FFC8BA","#BAFFFD"
  )
  color_map <- setNames(macaron_colors, top_orders)
  color_map["Others"] <- "#888888"
  
  # Plot each season
  for (mon in levels(df_plot$Month)) {
    plot_data <- df_plot %>% filter(Month == mon)
    
    p <- ggplot(plot_data, aes(x = "", y = Rel_Abund, fill = Order)) +
      geom_col(color = "white", size = 0.2) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = color_map) +
      labs(
        title = paste0(enz_name, " 目水平相对丰度（", mon, "）"),
        fill = "Order (目)"
      ) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold")
      )
    
    # Save plots
    ggsave(
      paste0(output_dir, "/", enz_name, "_", mon, "_目水平饼图.png"),
      plot = p, width = 8, height = 6, dpi = 300, bg = "white"
    )
    ggsave(
      paste0(output_dir, "/", enz_name, "_", mon, "_目水平饼图.pdf"),
      plot = p, width = 8, height = 6, bg = "white"
    )
  }
}

# Generate cellulase plots
cat("\n========== 绘制 所有纤维素酶总类 饼图 ==========\n")
plot_cellulase_pie(enz_name = "Total_Cellulase", enz_cols = all_enz_cols, topN = 10)

# Optional individual enzyme plots
# cat("\n========== 绘制 Cellulase（纤维素酶） 饼图 ==========\n")
# plot_cellulase_pie(enz_name = "Cellulase", enz_cols = "Cellulase", topN = 10)
# 
# cat("\n========== 绘制 Endoglucanase（内切葡聚糖酶） 饼图 ==========\n")
# plot_cellulase_pie(enz_name = "Endoglucanase", enz_cols = "Endoglucanase", topN = 10)
# 
# cat("\n========== 绘制 Exoglucanase（外切葡聚糖酶） 饼图 ==========\n")
# plot_cellulase_pie(enz_name = "Exoglucanase", enz_cols = "Exoglucanase", topN = 10)

# Shared legend
# Top 10 orders plus Others
top_orders_all <- df_merge %>%
  mutate(Total_enz = rowSums(select(., all_of(all_enz_cols)), na.rm = TRUE)) %>%
  group_by(Order) %>%
  summarise(Total = sum(Total_enz * (Month04 + Month10), na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10) %>%
  pull(Order)

macaron_colors <- c(
  "#FFB3BA","#FFDFBA","#FFFFBA","#BAFFC9","#BAE1FF",
  "#D6BAFF","#FFBAF0","#BAB0FF","#FFC8BA","#BAFFFD"
)
color_map <- setNames(macaron_colors, top_orders_all)
color_map["Others"] <- "#888888"

# Export legend
leg_p <- ggplot(data.frame(Order = names(color_map)), aes(x = "", y = 1, fill = Order)) +
  geom_col() +
  scale_fill_manual(values = color_map) +
  labs(fill = "Order (目)") +
  theme_void() +
  theme(
    legend.position = "center",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold")
  )

only_leg <- get_legend(leg_p)
ggsave(paste0(output_dir, "/统一图例.png"), only_leg, width = 6, height = 8, dpi = 300, bg = "white")
ggsave(paste0(output_dir, "/统一图例.pdf"), only_leg, width = 6, height = 8, bg = "white")

cat("1. 自动过滤了无酶数据的基因组，只保留bin_name和Genome匹配的有效数据\n")
cat("2. 04月/10月TPM为对应月份所有样本的均值，和原代码逻辑完全一致\n")
cat("3. 饼图为目水平的相对丰度，TOP10目单独展示，其余合并为Others\n")
