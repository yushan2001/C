# Lignin degradation analysis
# Order-level abundance, MAG phylogeny, and seasonal enzyme composition
# ======================================================================

# ======================================================================
# Part I: Order-level abundance
# ======================================================================
# ----------------------------------------------------------------------1. 安装并加载必要包 =====================
if (!require("readxl")) install.packages("readxl")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("scales")) install.packages("scales")

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd("C:/Users/余山小可爱/Desktop")
df_raw <- read_excel("木质素原始代码数据.xlsx", sheet = 1)  

# 数据清洗
df <- df_raw %>%
  select(`包含的Prc`, Order, `04总丰度`, `10总丰度`) %>%  
  filter(
    !is.na(`包含的Prc`), `包含的Prc` != "",
    !is.na(Order), Order != "",
    !is.na(`04总丰度`),
    !is.na(`10总丰度`)
  ) %>%
  mutate(
    `04总丰度` = as.numeric(`04总丰度`) %>% replace_na(0),
    `10总丰度` = as.numeric(`10总丰度`) %>% replace_na(0)
  ) %>%
  filter(`04总丰度` > 0 | `10总丰度` > 0)

# ----------------------------------------------------------------------3. 数据预处理 =====================
plot_data <- df %>%
  pivot_longer(
    cols = c(`04总丰度`, `10总丰度`),
    names_to = "Time",
    values_to = "Abundance"
  ) %>%
  mutate(
    Time = factor(Time,
                  levels = c("04总丰度", "10总丰度"),
                  labels = c("202504", "202510")),
    Order = factor(Order)
  )

# ----------------------------------------------------------------------4. 精选浅色系 + 大跨度打乱色板 =====================
# Order levels
unique_orders <- sort(unique(plot_data$Order))
n_orders <- length(unique_orders)

shuffled_light_palette <- c(
  "#FFB6C1", "#87CEFA", "#98FB98", "#DDA0DD", "#FFE4B5", "#7FFFD4", # 1-6:   浅粉、天蓝、薄荷绿、淡紫、奶油桔、碧绿
  "#F0E68C", "#F08080", "#B0E0E6", "#90EE90", "#EE82EE", "#FFD700", # 7-12:  淡金黄、浅珊瑚红、粉蓝、浅绿、兰花紫、明黄
  "#FFA07A", "#AFEEEE", "#E6E6FA", "#F5DEB3", "#9370DB", "#20B2AA", # 13-18: 浅鲑红、淡青绿、薰衣草紫、麦香色、中紫、浅海蓝
  "#FFC0CB", "#B0C4DE", "#8FBC8F", "#D8BFD8", "#FFF8DC", "#FA8072", # 19-24: 纯粉红、冰山蓝、暗海绿、蓟紫、米绢色、鲜鲑红
  "#E0FFFF", "#F4A460", "#9ACD32", "#F0F8FF", "#DEB887", "#E6A8D7", # 25-30: 浅青色、沙褐色、黄绿、爱丽丝蓝、硬木色、浅芭蕾粉
  "#7B68EE", "#FFDEAD", "#ADFF2F", "#BA55D3", "#FFE4E1", "#48D1CC", # 31-36: 暗灰蓝、白杏色、绿黄、中兰花紫、玫瑰白、中青绿
  "#F5F5DC", "#BC8F8F", "#ADD8E6", "#E0FFFE"                        # 37-40: 米色、褐玫瑰红、浅蓝、极淡青
)

# Select colors for observed orders
final_colors <- shuffled_light_palette[1:n_orders]

# Map colors to orders
order_color_map <- setNames(final_colors, unique_orders)

# ----------------------------------------------------------------------5. 创建保存文件夹 =====================
if (!dir.exists("Protein_Order_高辨识浅纯色图")) {
  dir.create("Protein_Order_高辨识浅纯色图")
}

# ----------------------------------------------------------------------6. 批量绘图函数 =====================
draw_plot <- function(prot_name) {
  sub_data <- plot_data %>% filter(`包含的Prc` == prot_name)
  
  if (nrow(sub_data) == 0 || all(sub_data$Abundance == 0)) {
    return(NULL)
  }
  
  p <- ggplot(sub_data, aes(x = Time, y = Abundance, fill = Order)) +
    geom_col(position = "stack", width = 0.8, color = "black", linewidth = 0.4) +
    scale_fill_manual(values = order_color_map, drop = FALSE) +
    labs(title = prot_name, x = "", y = "") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 18, color = "black", face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "none"
    )
  
  safe_name <- gsub("[^0-9a-zA-Z_]", "_", prot_name)
  save_path = "Protein_Order_高辨识浅纯色图/"
  
  ggsave(filename = paste0(save_path, safe_name, "_原始丰度图.png"),
         plot = p, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(filename = paste0(save_path, safe_name, "_原始丰度图.pdf"),
         plot = p, width = 8, height = 6, device = "pdf", bg = "white")
  
  return(p)
}

# ----------------------------------------------------------------------7. 批量生成所有图 =====================
all_proteins <- unique(plot_data$`包含的Prc`)
for (prot in all_proteins) {
  draw_plot(prot)
}

# ----------------------------------------------------------------------8. 统一图例 =====================
legend_df <- data.frame(
  Order = factor(unique_orders, levels = unique_orders),
  y = seq_along(unique_orders)
)

legend_plot <- ggplot(legend_df, aes(x = 1, y = y, fill = Order)) +
  geom_tile(color = "black", linewidth = 0.4) + # Black border for legend keys
  scale_fill_manual(values = order_color_map) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.key.size = unit(1.4, "cm"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(title = "Order（目水平）", ncol = 1, byrow = TRUE))

# Legend height
legend_height <- max(6, length(unique_orders) * 0.5 + 2)

ggsave("Protein_Order_高辨识浅纯色图/Order_统一图例.png",
       legend_plot, width = 6, height = legend_height, dpi = 300, bg = "white")
ggsave("Protein_Order_高辨识浅纯色图/Order_统一图例.pdf",
       legend_plot, width = 6, height = legend_height, device = "pdf", bg = "white")


# ======================================================================
# Part II: Lignin-related genes across MAGs
# ======================================================================
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

# Packages
if (!require("ggtree")) { install.packages("BiocManager"); BiocManager::install("ggtree") }
if (!require("ape")) install.packages("ape")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("viridis")) install.packages("viridis")
if (!require("cowplot")) install.packages("cowplot")
if (!require("ggnewscale")) install.packages("ggnewscale")
if (!require("tidytree")) install.packages("tidytree")
if (!require("readxl")) install.packages("readxl")

library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(readxl)
library(viridis)
library(cowplot)
library(ggnewscale)
library(tidytree)

# ----------------------------------------------------------------------1. 参数与路径配置 =====================
file_path  <- "C:/Users/余山小可爱/Desktop/木质素原始代码数据.xlsx"
tree_path  <- "C:/Users/余山小可爱/Desktop/gtdb_results.backbone.bac120.classify.tree"
output_dir <- "C:/Users/余山小可爱/Desktop/Output_Phylogeny_Pro/"
if (!dir.exists(output_dir)) dir.create(output_dir)

tree_line_size <- 0.6       # 树分枝线条粗细
tip_point_size <- 1.5       # 树末端 Phylum 节点大小
strip_offset   <- 0.05      # 【关键】树最外层尖端到第一圈热图的空白间距 (等同于最内环色带配置)
ring_width     <- 0.03      # 每圈热图的物理宽度
ring_gap       <- 0.005     # 热图圈与圈之间的微小空隙
border_size    <- 0.05      # 热图单元格的黑色边框粗细
hilight_alpha  <- 0.15      # 背景高亮扇区的透明度

gene_cols <- c(
  "DesA", "DesB", "DesZ", "DypB", "DyPs", "FerA", "FerB/FerB2", "LACs", "LDS", 
  "LigA", "LigB", "LigC", "LigD", "LigEFG", "LigH", "LigI", "LigJ", "LigK", 
  "LigM", "LigV", "LigW", "LigX", "LigY", "LigZ", "PcaB", "PcaC", "pcaD", 
  "PcaG", "pcaH", "PcaI", "PcaJ", "PraA", "PraB", "PraC", "PraD", "PraE", 
  "PraF", "PraG", "PraH", "VanA", "VanB"
)

# ----------------------------------------------------------------------2. 数据处理 =====================
full_df <- read_excel(file_path, sheet = "Sheet2") %>%
  rename(Bin = `Bin ID`) %>%
  mutate(Bin = trimws(Bin))

tree <- read.tree(tree_path)
keep_bins <- intersect(tree$tip.label, full_df$Bin)
filtered_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep_bins])
checkm_df <- full_df %>% filter(Bin %in% filtered_tree$tip.label)

for(col in gene_cols) {
  if(col %in% colnames(checkm_df)) {
    checkm_df[[col]] <- as.numeric(as.character(checkm_df[[col]]))
    checkm_df[[col]][is.na(checkm_df[[col]])] <- 0
  }
}

# Maximum gene copy number
max_actual_val <- checkm_df %>% 
  select(any_of(gene_cols)) %>% 
  max(na.rm = TRUE)

if(is.na(max_actual_val) || max_actual_val == 0) max_actual_val <- 10

# ----------------------------------------------------------------------3. 基础树绘制与动态半径获取 =====================
unique_phyla <- unique(checkm_df$Phylum)
phylum_color_map <- setNames(scales::hue_pal()(length(unique_phyla)), unique_phyla)

p_base <- ggtree(filtered_tree, layout = "circular", color = "black", size = tree_line_size)

max_radius <- max(p_base$data$x)
tip_coords <- p_base$data %>% filter(isTip) %>% select(label, y)

current_plot <- p_base

heatmap_total_width <- length(gene_cols) * (ring_width + ring_gap)
extend_length <- (max_radius + strip_offset + heatmap_total_width) - max_radius

for (phylum_name in unique_phyla) {
  nodes <- checkm_df %>% filter(Phylum == phylum_name) %>% pull(Bin)
  if (length(nodes) > 1) {
    mrca_node <- getMRCA(filtered_tree, nodes)
    if (!is.null(mrca_node)) {
      current_plot <- current_plot + 
        geom_hilight(node = mrca_node, fill = phylum_color_map[phylum_name], 
                     alpha = hilight_alpha, extend = extend_length)
    }
  }
}

# ----------------------------------------------------------------------5. 绘制绝对数量热图环 (基于自适应坐标体系) =====================
for (i in seq_along(gene_cols)) {
  gene_name <- gene_cols[i]
  if(!gene_name %in% colnames(checkm_df)) next
  
  ring_data <- checkm_df %>% 
    select(Bin, val = !!sym(gene_name)) %>%
    inner_join(tip_coords, by = c("Bin" = "label")) %>%
    mutate(x = max_radius + strip_offset + (i-1) * (ring_width + ring_gap))
  
  current_plot <- current_plot +
    new_scale_fill() +
    geom_tile(data = ring_data, aes(x = x, y = y, fill = val), 
              width = ring_width, height = 1, color = "black", 
              size = border_size, inherit.aes = FALSE) +
    scale_fill_gradientn(colors = c("white", "#deebf7", "#3182ce", "#084594"),
                         values = scales::rescale(c(0, 0.1, max_actual_val * 0.5, max_actual_val)), 
                         limits = c(0, max_actual_val), oob = scales::squish, guide = "none")
}

# ----------------------------------------------------------------------6. 添加  树尖端末端点 (Phylum) =====================
final_plot <- current_plot + 
  new_scale_color() + 
  geom_tippoint(data = p_base$data %>% filter(isTip) %>% left_join(checkm_df, by = c("label" = "Bin")), 
                aes(x = x, y = y, color = Phylum), size = tip_point_size, alpha = 0.9, inherit.aes = FALSE) +
  scale_color_manual(values = phylum_color_map, name = "Phylum")

# ----------------------------------------------------------------------7. 导出最终结果与真实数量图例 =====================
legend_plot <- ggplot(checkm_df) +
  geom_point(aes(x = 1, y = 1, color = Phylum)) +
  scale_color_manual(values = phylum_color_map) +
  new_scale_fill() +
  geom_tile(aes(x = 1, y = 1, fill = DesA)) + 
  scale_fill_gradientn(colors = c("white", "#deebf7", "#3182ce", "#084594"),
                       values = scales::rescale(c(0, 0.1, max_actual_val * 0.5, max_actual_val)), 
                       limits = c(0, max_actual_val), 
                       breaks = seq(0, max_actual_val, by = max(1, round(max_actual_val / 5))),
                       name = "Gene Copy Number\n(基因拷贝绝对数量)") +
  theme_bw()

all_legend <- get_legend(legend_plot)

# Extend plotting range for outer rings
final_main <- final_plot + theme(legend.position = "none") + 
  xlim(0, max_radius + strip_offset + heatmap_total_width + 0.1)

ggsave(paste0(output_dir, "Tree_Absolute_Counts_Adaptive.png"), final_main, 
       width = 18, height = 18, dpi = 600, device = png(type = "cairo"))
ggsave(paste0(output_dir, "Tree_Absolute_Counts_Adaptive.pdf"), final_main, 
       width = 18, height = 18, device = "pdf")
ggsave(paste0(output_dir, "Legends_Final.pdf"), ggdraw(all_legend), width = 5, height = 10)

# ----------------------------------------------------------------------8. 运行反馈 =====================


# ======================================================================
# Part III: Seasonal order-level lignin enzyme composition
# ======================================================================
if (!require(readxl)) install.packages("readxl")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggpubr)) install.packages("ggpubr")
library(readxl)
library(tidyverse)
library(ggpubr)

# ----------------------------------------------------------------------1. 输出文件夹（适配木质素酶主题） =====================
output_dir <- "C:/Users/余山小可爱/Desktop/木质素酶家族_季节对比饼图"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Sheet3: lignin enzyme counts by bin
df_enz <- read_excel("C:/Users/余山小可爱/Desktop/木质素原始代码数据.xlsx", sheet = "Sheet3")
# Sheet4: taxonomy and sample TPM by genome
df_tpm <- read_excel("C:/Users/余山小可爱/Desktop/木质素原始代码数据.xlsx", sheet = "Sheet4")

df_enz_renamed <- df_enz %>% rename(Genome = bin_name)

df_merge <- df_tpm %>%
  inner_join(df_enz_renamed, by = "Genome") %>%
  filter(!is.na(Order)) # Remove rows without order-level taxonomy

# Lignin enzyme columns
all_enz_cols <- setdiff(colnames(df_enz_renamed), "Genome")

df_merge <- df_merge %>%
  mutate(
    Month04 = rowMeans(select(., starts_with("04")), na.rm = TRUE),
    Month10 = rowMeans(select(., starts_with("10")), na.rm = TRUE)
  )

# ----------------------------------------------------------------------5. 通用绘图函数（适配木质素酶，可做总类/单个酶） =====================
plot_lignin_pie <- function(enz_name = "Total_Ligninase", enz_cols = all_enz_cols, topN = 10) {
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
    filter(`04月` > 0 | `10月` > 0) # Remove orders with zero abundance
  
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

# ----------------------------------------------------------------------6. 批量绘图（总类+关键单酶，可按需调整） =====================
cat("\n========== 绘制 所有木质素酶总类 饼图 ==========\n")
plot_lignin_pie(enz_name = "Total_Ligninase", enz_cols = all_enz_cols, topN = 10)

# cat("\n========== 绘制 LigB（木质素β-醚酶） 饼图 ==========\n")
# plot_lignin_pie(enz_name = "LigB", enz_cols = "LigB", topN = 10)
# 
# cat("\n========== 绘制 pcaD（原儿茶酸脱羧酶） 饼图 ==========\n")
# plot_lignin_pie(enz_name = "pcaD", enz_cols = "pcaD", topN = 10)
# 
# cat("\n========== 绘制 DyPs（过氧化物酶） 饼图 ==========\n")
# plot_lignin_pie(enz_name = "DyPs", enz_cols = "DyPs", topN = 10)
# 
# cat("\n========== 绘制 DesA（芳香还原酶） 饼图 ==========\n")
# plot_lignin_pie(enz_name = "DesA", enz_cols = "DesA", topN = 10)

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

cat("1. 自动过滤了无酶数据的基因组，只保留bin_name和Genome完全匹配的有效数据\n")
cat("2. 04月/10月TPM为对应月份所有样本的均值，和原代码逻辑100%一致\n")
cat("3. 饼图为目水平的相对丰度，TOP10目单独展示，其余合并为Others\n")
cat("4. 已预留你LEfSe分析核心差异基因的单酶绘图入口，按需取消注释即可运行\n")
