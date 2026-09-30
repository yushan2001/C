# CAZyme analysis workflow
#
# This script combines the phylogenetic annotation, seasonal comparison,
# heatmap, NMDS/PERMANOVA, and taxonomic contribution analyses used in the study.
# Input/output paths and analysis parameters are retained from the working scripts.

酶家族环形树
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

# ===================== 核心参数：可手动微调数值范围 =====================
# 酶家族数值范围（NA=自动适配Excel最大值/最小值）
aas_min <- NA; aas_max <- NA
cbms_min <- NA; cbms_max <- NA
ces_min <- NA; ces_max <- NA
ghs_min <- NA; ghs_max <- NA
gts_min <- NA; gts_max <- NA
pls_min <- NA; pls_max <- NA
cazy_num_min <- NA; cazy_num_max <- NA

# 绘图布局参数
strip_offset <- 0.1    # 条带间距
strip_width <- 0.8    # 条带宽度
bar_width <- 0.8      # 柱状图宽度
tree_line_size <- 0.8
tip_point_size <- 2.5
hilight_alpha <- 0.3
hilight_extend <- 0.05

# ================================================================================

# 读取 酶家族.xlsx → Sheet1
checkm_path <- "C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx"
checkm_df <- read_excel(checkm_path, sheet = "Sheet1") %>%  # 读取Sheet1
  mutate(
    Genome = trimws(Genome),          
    Phylum = trimws(Phylum)     
  )

excel_bins <- checkm_df$Genome  # 这里改成 Genome
cat("Excel中提取到", length(excel_bins), "个目标bin\n")

# 读取进化树并过滤
tree_path <- "C:/Users/余山小可爱/Desktop/gtdb_results.backbone.bac120.classify.tree"
tree <- read.tree(tree_path)
keep_index <- which(tree$tip.label %in% excel_bins)
filtered_tree <- drop.tip(tree, tree$tip.label[-keep_index])
cat("进化树过滤后保留", length(filtered_tree$tip.label), "个目标bin\n")

# 清洗数据+确保酶家族列是数值型
checkm_df <- checkm_df %>%
  filter(Genome %in% filtered_tree$tip.label) %>%  # 这里改成 Genome
  mutate(
    AAs = as.numeric(AAs), AAs = ifelse(is.na(AAs), 0, AAs),
    CBMs = as.numeric(CBMs), CBMs = ifelse(is.na(CBMs), 0, CBMs),
    CEs = as.numeric(CEs), CEs = ifelse(is.na(CEs), 0, CEs),
    GHs = as.numeric(GHs), GHs = ifelse(is.na(GHs), 0, GHs),
    GTs = as.numeric(GTs), GTs = ifelse(is.na(GTs), 0, GTs),
    PLs = as.numeric(PLs), PLs = ifelse(is.na(PLs), 0, PLs),
    `number of CAZyme families` = as.numeric(`number of CAZyme families`),
    `number of CAZyme families` = ifelse(is.na(`number of CAZyme families`), 0, `number of CAZyme families`),
    Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum)
  )

# 自动计算数值范围
if (is.na(aas_min)) aas_min <- min(checkm_df$AAs, na.rm = TRUE); if (is.na(aas_max)) aas_max <- max(checkm_df$AAs, na.rm = TRUE)
if (is.na(cbms_min)) cbms_min <- min(checkm_df$CBMs, na.rm = TRUE); if (is.na(cbms_max)) cbms_max <- max(checkm_df$CBMs, na.rm = TRUE)
if (is.na(ces_min)) ces_min <- min(checkm_df$CEs, na.rm = TRUE); if (is.na(ces_max)) ces_max <- max(checkm_df$CEs, na.rm = TRUE)
if (is.na(ghs_min)) ghs_min <- min(checkm_df$GHs, na.rm = TRUE); if (is.na(ghs_max)) ghs_max <- max(checkm_df$GHs, na.rm = TRUE)
if (is.na(gts_min)) gts_min <- min(checkm_df$GTs, na.rm = TRUE); if (is.na(gts_max)) gts_max <- max(checkm_df$GTs, na.rm = TRUE)
if (is.na(pls_min)) pls_min <- min(checkm_df$PLs, na.rm = TRUE); if (is.na(pls_max)) pls_max <- max(checkm_df$PLs, na.rm = TRUE)
if (is.na(cazy_num_min)) cazy_num_min <- min(checkm_df$`number of CAZyme families`, na.rm = TRUE)
if (is.na(cazy_num_max)) cazy_num_max <- max(checkm_df$`number of CAZyme families`, na.rm = TRUE)

cat("\n自动适配的数值范围：\n")
cat("AAs：", aas_min, "-", aas_max, "\n")
cat("CBMs：", cbms_min, "-", cbms_max, "\n")
cat("CEs：", ces_min, "-", ces_max, "\n")
cat("GHs：", ghs_min, "-", ghs_max, "\n")
cat("GTs：", gts_min, "-", gts_max, "\n")
cat("PLs：", pls_min, "-", pls_max, "\n")
cat("CAZyme家族数量：", cazy_num_min, "-", cazy_num_max, "\n")

taxonomy_for_tree <- checkm_df %>%
  mutate(label = Genome)  # 这里改成 Genome

unique_phyla <- unique(taxonomy_for_tree$Phylum)
if (length(unique_phyla) == 0) unique_phyla <- "Unclassified"
phylum_colors <- scales::hue_pal()(length(unique_phyla))
phylum_color_map <- setNames(phylum_colors, unique_phyla)

circular_plot_base <- ggtree(filtered_tree, layout = "circular", color = "black", size = tree_line_size)

tree_data <- circular_plot_base$data %>%
  left_join(taxonomy_for_tree %>% select(label, Phylum), by = "label") %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum))
circular_plot_base$data <- tree_data

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

max_radius <- max(circular_plot_base$data$x)

strip_data <- circular_plot_base$data %>%
  filter(isTip) %>%
  left_join(checkm_df, by = c("label" = "Genome")) %>%  # 这里改成 Genome
  mutate(
    aas_xmin = max_radius + strip_offset,
    aas_xmax = max_radius + strip_offset + strip_width,
    cbms_xmin = max_radius + strip_offset*2 + strip_width,
    cbms_xmax = max_radius + strip_offset*2 + strip_width*2,
    ces_xmin = max_radius + strip_offset*3 + strip_width*2,
    ces_xmax = max_radius + strip_offset*3 + strip_width*3,
    ghs_xmin = max_radius + strip_offset*4 + strip_width*3,
    ghs_xmax = max_radius + strip_offset*4 + strip_width*4,
    gts_xmin = max_radius + strip_offset*5 + strip_width*4,
    gts_xmax = max_radius + strip_offset*5 + strip_width*5,
    pls_xmin = max_radius + strip_offset*6 + strip_width*5,
    pls_xmax = max_radius + strip_offset*6 + strip_width*6,
    cazy_xmin = max_radius + strip_offset*7 + strip_width*6,
    cazy_xmax = max_radius + strip_offset*7 + strip_width*6 + 
      ((`number of CAZyme families` - cazy_num_min)/(cazy_num_max - cazy_num_min)) * strip_width*2
  )

enzyme_palettes <- list(
  AAs = c("#FFE4E1", "#FFB6C1", "#FF69B4", "#DC143C"),
  CBMs = c("#BBD9F0", "#7FB8E0", "#3D88C0", "#005A99"),
  CEs = c("#B8E6C8", "#7FCC99", "#40B366", "#00802B"),
  GHs = c("#FFF380", "#FFEA40", "#FFE000", "#CCB300"),
  GTs = c("#F0E6FA", "#DDA0DD", "#C9A0DC", "#B37FCC"),
  PLs = c("#FFB3B3", "#FF7070", "#FF2828", "#CC0000"),
  CAZy = c("#A8E6C8", "#66CC99", "#2DB366", "#008039")
)

plot_with_legends <- circular_plot_base +
  geom_rect(data = strip_data, aes(xmin = aas_xmin + 0.1, xmax = aas_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = AAs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("AAs [", round(aas_min,1), "-", round(aas_max,1), "]"), colours = enzyme_palettes$AAs, limits = c(aas_min, aas_max), guide = guide_colorbar(order = 1, title.position = "top")) + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = cbms_xmin + 0.1, xmax = cbms_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = CBMs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("CBMs [", round(cbms_min,1), "-", round(cbms_max,1), "]"), colours = enzyme_palettes$CBMs, limits = c(cbms_min, cbms_max), guide = guide_colorbar(order = 2, title.position = "top")) + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = ces_xmin + 0.1, xmax = ces_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = CEs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("CEs [", round(ces_min,1), "-", round(ces_max,1), "]"), colours = enzyme_palettes$CEs, limits = c(ces_min, ces_max), guide = guide_colorbar(order = 3, title.position = "top")) + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = ghs_xmin + 0.1, xmax = ghs_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = GHs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("GHs [", round(ghs_min,1), "-", round(ghs_max,1), "]"), colours = enzyme_palettes$GHs, limits = c(ghs_min, ghs_max), guide = guide_colorbar(order = 4, title.position = "top")) + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = gts_xmin + 0.1, xmax = gts_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = GTs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("GTs [", round(gts_min,1), "-", round(gts_max,1), "]"), colours = enzyme_palettes$GTs, limits = c(gts_min, gts_max), guide = guide_colorbar(order = 5, title.position = "top")) + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = pls_xmin + 0.1, xmax = pls_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = PLs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("PLs [", round(pls_min,1), "-", round(pls_max,1), "]"), colours = enzyme_palettes$PLs, limits = c(pls_min, pls_max), guide = guide_colorbar(order = 6, title.position = "top")) + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = cazy_xmin, xmax = cazy_xmax, ymin = y - bar_width/2, ymax = y + bar_width/2, fill = `number of CAZyme families`), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(name = paste0("Number of CAZyme Families [", round(cazy_num_min,1), "-", round(cazy_num_max,1), "]"), colours = enzyme_palettes$CAZy, limits = c(cazy_num_min, cazy_num_max), guide = guide_colorbar(order = 7, title.position = "top")) + new_scale_fill() +
  geom_tippoint(aes(color = Phylum), size = tip_point_size, alpha = 0.9) +
  scale_color_manual(values = phylum_color_map, name = "Phylum", guide = guide_legend(order = 8, title.position = "top")) +
  theme_tree2() + ggtitle("Circular Phylogenetic Tree with CAZyme Family Annotations") +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 12, face = "bold"))

legend <- get_legend(plot_with_legends)
legend_components <- list(AAs = legend$grobs[[1]], CBMs = legend$grobs[[2]], CEs = legend$grobs[[3]], GHs = legend$grobs[[4]], GTs = legend$grobs[[5]], PLs = legend$grobs[[6]], CAZyme_Families = legend$grobs[[7]], Phylum = legend$grobs[[8]])

output_dir <- "C:/Users/余山小可爱/Desktop/legends-cazys/"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (name in names(legend_components)) {
  single_legend <- ggplot() + annotation_custom(grob = legend_components[[name]]) + theme_void() +
    theme(plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA))
  
  filename_png <- paste0(output_dir, name, "_legend-cazys.png")
  ggsave(filename_png, plot = single_legend, width = 5, height = 3, dpi = 600, bg = "transparent", device = png(type = "cairo"))
  
  filename_pdf <- paste0(output_dir, name, "_legend-cazys.pdf")
  ggsave(filename = filename_pdf, plot = single_legend, width = 5, height = 3, bg = "transparent", device = "pdf", useDingbats = FALSE)
  
  cat("已保存", name, "图例至:", filename_png, " + ", filename_pdf, "\n")
}

main_plot <- circular_plot_base +
  geom_rect(data = strip_data, aes(xmin = aas_xmin + 0.1, xmax = aas_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = AAs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = enzyme_palettes$AAs, limits = c(aas_min, aas_max), guide = "none") + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = cbms_xmin + 0.1, xmax = cbms_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = CBMs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = enzyme_palettes$CBMs, limits = c(cbms_min, cbms_max), guide = "none") + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = ces_xmin + 0.1, xmax = ces_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = CEs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = enzyme_palettes$CEs, limits = c(ces_min, ces_max), guide = "none") + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = ghs_xmin + 0.1, xmax = ghs_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = GHs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = enzyme_palettes$GHs, limits = c(ghs_min, ghs_max), guide = "none") + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = gts_xmin + 0.1, xmax = gts_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = GTs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = enzyme_palettes$GTs, limits = c(gts_min, gts_max), guide = "none") + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = pls_xmin + 0.1, xmax = pls_xmax - 0.1, ymin = y - 0.5, ymax = y + 0.5, fill = PLs), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = enzyme_palettes$PLs, limits = c(pls_min, pls_max), guide = "none") + new_scale_fill() +
  geom_rect(data = strip_data, aes(xmin = cazy_xmin, xmax = cazy_xmax, ymin = y - bar_width/2, ymax = y + bar_width/2, fill = `number of CAZyme families`), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_gradientn(colours = enzyme_palettes$CAZy, limits = c(cazy_num_min, cazy_num_max), guide = "none") + new_scale_fill() +
  geom_tippoint(aes(color = Phylum), size = tip_point_size, alpha = 0.9) +
  scale_color_manual(values = phylum_color_map, guide = "none") +
  theme_tree2() + ggtitle("Circular Phylogenetic Tree with CAZyme Family Annotations") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA))

ggsave(filename = "C:/Users/余山小可爱/Desktop/main_plot_cazy.png", plot = main_plot, width = 22, height = 22, dpi = 600, bg = "transparent", device = png(type = "cairo"))
ggsave(filename = "C:/Users/余山小可爱/Desktop/main_plot_cazy.pdf", plot = main_plot, width = 22, height = 22, bg = "transparent", device = "pdf", useDingbats = FALSE)

cat("主图已保存至: C:/Users/余山小可爱/Desktop/main_plot_cazy.png + main_plot_cazy.pdf\n")
cat("\n所有图例(PNG+PDF)和主图(PNG+PDF)均保存完成。\n")


# 酶家族显著性结果.xlsx
# ===================== 1. 加载包 =====================
if (!require("readxl")) install.packages("readxl")
if (!require("writexl")) install.packages("writexl")
if (!require("dplyr")) install.packages("dplyr")

library(readxl)
library(writexl)
library(dplyr)

# ===================== 2. 读取数据 =====================
file_path <- "C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx"
df <- read_excel(file_path, sheet = "Sheet2")

# ===================== 3. 查看列名 =====================
cat("原始数据列名：\n")
print(colnames(df))

# ===================== 4. 定义分组 =====================
cols_04 <- c("04DC","04JH","04MH","04RC","04SX","04TC","04TJ")
cols_10 <- c("10DC","10JH","10MH","10RC","10SX","10TC","10TJ")

# ===================== 5. Wilcoxon检验 =====================
result <- df %>%
  rowwise() %>%
  mutate(
    # 提取两组数据
    group04 = list(c_across(all_of(cols_04))),
    group10 = list(c_across(all_of(cols_10))),
    
    # 计算均值
    mean_04 = mean(unlist(group04), na.rm = TRUE),
    mean_10 = mean(unlist(group10), na.rm = TRUE),
    
    # fold change（防止除0）
    fold_change = ifelse(mean_04 == 0, NA, mean_10 / mean_04),
    
    # Wilcoxon检验
    p_value = tryCatch(
      wilcox.test(unlist(group04), unlist(group10))$p.value,
      error = function(e) NA
    ),
    
    # 显著性（原始p值）
    差异酶家族 = ifelse(!is.na(p_value) & p_value < 0.05, "是", "否")
  ) %>%
  ungroup()

# ===================== 6. 多重检验校正 =====================
result <- result %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    差异酶家族_FDR = ifelse(!is.na(p_adj) & p_adj < 0.05, "是", "否")
  )

# ===================== 7. 整理输出列 =====================
result <- result %>%
  select(CAZy_Family,
         mean_04, mean_10, fold_change,
         p_value, p_adj,
         差异酶家族, 差异酶家族_FDR,
         everything())

# ===================== 8. 输出结果 =====================
output_path <- "C:/Users/余山小可爱/Desktop/酶家族显著性结果_Wilcoxon.xlsx"
write_xlsx(result, output_path)

# ===================== 9. 输出统计 =====================
cat("Wilcoxon分析完成！\n")
cat("结果已保存至：", output_path, "\n\n")

cat("原始p值显著统计：\n")
print(table(result$差异酶家族))

cat("\nFDR校正后显著统计（推荐看这个）：\n")
print(table(result$差异酶家族_FDR))

# CAZy_heatmap_Zscore.png
# 安装加载包
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("ggh4x")) install.packages("ggh4x")

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(ggh4x) # 用于分栏独立着色

# ===================== 1. 读取数据 =====================
file_path <- "C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx"
raw_data <- read_excel(file_path, sheet = "Sheet3")

# 定义样本列
cols_04 <- c("04DC","04JH","04MH","04RC","04SX","04TC","04TJ")
cols_10 <- c("10DC","10JH","10MH","10RC","10SX","10TC","10TJ")

# ===================== 2. 数据整理 + Z-score 标准化 =====================
plot_data <- raw_data %>%
  select(CAZy_Family, all_of(cols_04), all_of(cols_10)) %>%
  pivot_longer(-CAZy_Family, names_to="Sample", values_to="Value") %>%
  mutate(
    Time = ifelse(grepl("^04", Sample), "04", "10"),
    Value = as.numeric(Value)
  ) %>%
  group_by(CAZy_Family) %>%
  mutate(Relative = scale(Value)) %>%
  ungroup() %>%
  mutate(CAZy_Family = factor(CAZy_Family, levels = unique(raw_data$CAZy_Family)))

# ===================== 3. 绘制热图（分栏颜色独立） =====================
p <- ggplot(plot_data, aes(x=Sample, y=CAZy_Family, fill=Relative)) +
  geom_tile(color="white", size=0.15) +
  facet_wrap2(~Time, scales="free_x", nrow=1,
              strip = strip_themed(
                background_x = list(
                  element_rect(fill = "#FFB3DD"),  # 04 粉色
                  element_rect(fill = "#ACFFFF")   # 10 青色
                ),
                text_x = element_text(color = "white", size = 12, face = "bold")
              )) +
  scale_fill_gradient2(
    low = "#2E86AB",
    mid = "white",
    high = "#E63946",
    midpoint = 0,
    name = "Z-score\n(Normalized Expression)"
  ) +
  labs(title="CAZy family abundance heatmap (Z-score normalized)",
       x="Sample", y="CAZy Family") +
  theme_bw() +
  theme(
    plot.title=element_text(hjust=0.5,size=15,face="bold"),
    axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=9),
    axis.text.y=element_text(size=8),
    panel.spacing = unit(0.2, "lines")
  ) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0))

# ===================== 4. 保存图片 =====================
ggsave("C:/Users/余山小可爱/Desktop/CAZy_heatmap_Zscore.png", p, width=10, height=16, dpi=300)
ggsave("C:/Users/余山小可爱/Desktop/CAZy_heatmap_Zscore.pdf", p, width=10, height=16, device="pdf", useDingbats=FALSE)

write.xlsx(plot_data, "C:/Users/余山小可爱/Desktop/CAZy_heatmap_data_Zscore.xlsx", overwrite=T)

cat("完成！04=#FFB3DD | 10=#ACFFFF \n")
print(p)


# 酶家族NMDS图_英文标题图例.png
# ===================== 仅加载基础包 =====================
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

# ===================== 核心分析代码（已替换为你的新数据） =====================
# 1. 读取数据：酶家族.xlsx → Sheet7
file_path <- "C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx.xlsx"
df <- read_excel(file_path, sheet = "Sheet4")

# 2. 数据预处理（适配你的样本格式：04DC / 10DC）
df_long <- df %>%
  pivot_longer(cols = -EnzymeFamily, names_to = "Sample", values_to = "Abundance") %>%
  pivot_wider(names_from = EnzymeFamily, values_from = Abundance)

# 构建丰度矩阵
df_selected <- df_long %>% select(-Sample)
rownames(df_selected) <- df_long$Sample
df_selected[is.na(df_selected)] <- 0
df_transposed <- as.matrix(df_selected)

# 3. NMDS分析
set.seed(42)
nmds_result <- metaMDS(df_transposed, distance = "bray", k = 2, trymax = 100)

# 4. 提取坐标 + 自动分组（04=April 2025，10=October 2025）
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_coords$Time <- ifelse(grepl("^04", rownames(nmds_coords)), "April 2025", "October 2025")
group_factor <- factor(nmds_coords$Time)

# 5. Adonis 统计检验
adonis_result <- adonis2(df_transposed ~ group_factor, 
                         data = nmds_coords, 
                         method = "bray", 
                         permutations = 999)
adonis_r2 <- round(adonis_result$R2[1], 4)
adonis_p <- round(adonis_result$`Pr(>F)`[1], 4)

# 6. 分组凸包（聚类圈）
calculate_hull <- function(data) {
  data[chull(data$NMDS1, data$NMDS2), ]
}
hull_data <- nmds_coords %>%
  group_by(Time) %>%
  do(calculate_hull(.))

# 7. 绘制NMDS图
nmds_plot <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Time, shape = Time)) +
  geom_polygon(data = hull_data, aes(fill = Time), alpha = 0.2, show.legend = FALSE) +
  geom_point(size = 5, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = c("April 2025" = "#FFB3DD", "October 2025" = "#ACFFFF")) +
  scale_fill_manual(values = c("April 2025" = "#FFB3DD", "October 2025" = "#ACFFFF")) +
  scale_shape_manual(values = c("April 2025" = 16, "October 2025" = 17)) +
  labs(title = "NMDS Analysis of Enzyme Family (April 2025 vs October 2025)",
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

# 显示 + 保存图片到桌面
print(nmds_plot)
ggsave("C:/Users/余山小可爱/Desktop/酶家族NMDS图_英文标题图例.png", 
       plot = nmds_plot, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("C:/Users/余山小可爱/Desktop/酶家族NMDS图_英文标题图例.pdf", 
       plot = nmds_plot, width = 10, height = 8, bg = "white")

# 输出结果
cat("=== Enzyme Family NMDS Analysis Results ===\n")
cat("NMDS Stress value：", round(nmds_result$stress, 4), "\n")
cat("Adonis analysis R²：", adonis_r2, "\n")
cat("Adonis analysis p-value：", adonis_p, "\n")
cat("酶家族NMDS图 PNG+PDF 已保存至桌面。\n")








# 目水平饼图

# ===================== 加载包 =====================
if (!require(readxl)) install.packages("readxl")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggpubr)) install.packages("ggpubr")
library(readxl)
library(tidyverse)
library(ggpubr)

# ===================== 1. 输出文件夹 =====================
output_dir <- "C:/Users/余山小可爱/Desktop/三大酶家族_季节对比饼图"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("输出文件夹已创建：", output_dir, "\n")

# ===================== 2. 读取数据 =====================
df_tpm <- read_excel("C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx", sheet = "Sheet6")
df_enz <- read_excel("C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx", sheet = "Sheet4")

# ===================== 3. 去重分类列，按Genome合并 =====================
df_tpm_clean <- df_tpm %>% select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species))
df_merge <- df_tpm_clean %>%
  inner_join(df_enz, by = "Genome") %>%
  filter(!is.na(Order))

# ===================== 4. 计算04月/10月平均TPM =====================
df_merge <- df_merge %>%
  mutate(
    Month04 = rowMeans(select(., starts_with("04")), na.rm = TRUE),
    Month10 = rowMeans(select(., starts_with("10")), na.rm = TRUE)
  )

# ===================== 5. 定义绘图函数（通用AA/GH/CE） =====================
plot_season_pie <- function(prefix, topN = 10) {
  # 提取对应前缀的所有酶家族，汇总为该大类的总酶数
  all_cols <- colnames(df_merge)
  fam_cols <- all_cols[str_starts(all_cols, prefix)]
  cat(prefix, "家族共", length(fam_cols), "个\n")
  
  # 计算每个bin的该大类总酶数
  df_tmp <- df_merge %>%
    mutate(Total_enz = rowSums(select(., all_of(fam_cols)), na.rm = TRUE))
  
  # 按目汇总两个月份的总TPM
  df_sum <- df_tmp %>%
    group_by(Order) %>%
    summarise(
      `04月` = sum(Total_enz * Month04, na.rm = TRUE),
      `10月` = sum(Total_enz * Month10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(`04月` > 0 | `10月` > 0)
  
  # 统计所有目总丰度，取TOP10
  top_orders <- df_sum %>%
    mutate(Total = `04月` + `10月`) %>%
    arrange(desc(Total)) %>%
    slice_head(n = topN) %>%
    pull(Order)
  
  # 非TOP10目合并为Others
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
  
  # 统一马卡龙配色
  macaron_colors <- c(
    "#FFB3BA","#FFDFBA","#FFFFBA","#BAFFC9","#BAE1FF",
    "#D6BAFF","#FFBAF0","#BAB0FF","#FFC8BA","#BAFFFD"
  )
  color_map <- setNames(macaron_colors, top_orders)
  color_map["Others"] <- "#888888"
  
  # 分季节绘制饼图
  for (mon in levels(df_plot$Month)) {
    plot_data <- df_plot %>% filter(Month == mon)
    
    p <- ggplot(plot_data, aes(x = "", y = Rel_Abund, fill = Order)) +
      geom_col(color = "white", size = 0.2) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = color_map) +
      labs(
        title = paste0(prefix, "家族 目水平相对丰度（", mon, "）"),
        fill = "Order (目)"
      ) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold")
      )
    
    # 保存文件
    ggsave(
      paste0(output_dir, "/", prefix, "_", mon, "_目水平饼图.png"),
      plot = p, width = 8, height = 6, dpi = 300, bg = "white"
    )
    ggsave(
      paste0(output_dir, "/", prefix, "_", mon, "_目水平饼图.pdf"),
      plot = p, width = 8, height = 6, bg = "white"
    )
  }
  cat("已完成", prefix, "的两个季节饼图\n")
}

# ===================== 6. 绘制三大类饼图 =====================
cat("\n========== 绘制 AA 家族 ==========\n")
plot_season_pie("AA", 10)

cat("\n========== 绘制 GH 家族 ==========\n")
plot_season_pie("GH", 10)

cat("\n========== 绘制 CE 家族 ==========\n")
plot_season_pie("CE", 10)

# ===================== 7. 单独输出统一图例 =====================
# 提取TOP10目+Others的统一图例
all_orders <- unique(c(
  df_merge %>% filter(str_starts(colnames(df_merge), "AA") %>% any()) %>% pull(Order),
  df_merge %>% filter(str_starts(colnames(df_merge), "GH") %>% any()) %>% pull(Order),
  df_merge %>% filter(str_starts(colnames(df_merge), "CE") %>% any()) %>% pull(Order)
)) %>% unique() %>% .[1:10]

macaron_colors <- c(
  "#FFB3BA","#FFDFBA","#FFFFBA","#BAFFC9","#BAE1FF",
  "#D6BAFF","#FFBAF0","#BAB0FF","#FFC8BA","#BAFFFD"
)
color_map <- setNames(macaron_colors, all_orders)
color_map["Others"] <- "#888888"

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

cat("\n全部完成。6张饼图已保存到：\n", output_dir)

# 指定核心酶家族_统一图例平面饼图
# ==============================================================================
# 指定核心酶家族“完全独立出图” —— 3D立体 + 侧边栏全局统一图例 ()
# ==============================================================================

# ===================== 加载与验证包 =====================
if (!require(readxl)) install.packages("readxl")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(plotrix)) install.packages("plotrix") 
library(readxl)
library(tidyverse)
library(plotrix)

# ===================== 1. 输出文件夹 =====================
output_dir <- "C:/Users/余山小可爱/Desktop/指定核心酶家族_全局统一3D饼图"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("输出文件夹已创建：", output_dir, "\n")

# ===================== 2. 读取数据 =====================
df_tpm <- read_excel("C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx", sheet = "Sheet6")
df_enz <- read_excel("C:/Users/余山小可爱/Desktop/酶家族代码原始数据.xlsx", sheet = "Sheet4")

# ===================== 3. 去重分类列，按Genome合并 =====================
df_tpm_clean <- df_tpm %>% select(-c(Kingdom, Phylum, Class, Order, Family, Genus, Species))
df_merge <- df_tpm_clean %>%
  inner_join(df_enz, by = "Genome") %>%
  filter(!is.na(Order))

# ===================== 4. 计算04月/10月平均TPM =====================
df_merge <- df_merge %>%
  mutate(
    Month04 = rowMeans(select(., starts_with("04")), na.rm = TRUE),
    Month10 = rowMeans(select(., starts_with("10")), na.rm = TRUE)
  )

# ===================== 5. 确定并验证目标酶家族 =====================
target_enzymes <- c("GH5", "GH1", "CE1", "GH3", "GH43", "GH51", "GH20", "GH53", "GH2", "GH4", "AA1", "AA3") %>% unique()
all_cols <- colnames(df_merge)
valid_enzymes <- target_enzymes[target_enzymes %in% all_cols]

# ===================== 6. 计算指定基因群的全局 Top 10 个目 =====================
cat("\n计算指定酶家族的全局 Top 10 目...\n")

df_global_ranking <- df_merge %>%
  select(Order, Month04, Month10, all_of(valid_enzymes)) %>%
  pivot_longer(cols = all_of(valid_enzymes), names_to = "Enzyme", values_to = "Count") %>%
  filter(Count > 0) %>%
  group_by(Order) %>%
  summarise(
    Total_Global_TPM = sum(Count * Month04, na.rm = TRUE) + sum(Count * Month10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Global_TPM))

# 锁定前 10 个核心目名字并确定全局因子的层级
global_top_10_orders <- df_global_ranking %>% slice_head(n = 10) %>% pull(Order)
final_order_levels <- c(global_top_10_orders, "Others")

# 锁定全局统一调色盘（学术马卡龙 10 色 + 1个灰色给 Others）
macaron_base <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF",
                  "#D6BAFF", "#FFBAF0", "#BAB0FF", "#FFC8BA", "#BAFFFD")
global_color_map <- setNames(macaron_base[1:length(global_top_10_orders)], global_top_10_orders)
global_color_map["Others"] <- "#888888"

# ===================== 7. 循环独立清洗并绘制每个酶家族的 3D 立体饼图 =====================
for (enz in valid_enzymes) {
  cat("处理酶家族：", enz, "\n")
  
  # 提取当前单个酶家族的数据
  df_single <- df_merge %>%
    select(Order, Month04, Month10, current_enz = all_of(enz)) %>%
    filter(current_enz > 0)
  
  if (nrow(df_single) == 0) next
  
  # 按目汇总两个月份下的表达量
  df_sum <- df_single %>%
    group_by(Order) %>%
    summarise(
      `04月` = sum(current_enz * Month04, na.rm = TRUE),
      `10月` = sum(current_enz * Month10, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 全员严格套用全局 Top 10 规则
  df_plot <- df_sum %>%
    mutate(Order = ifelse(Order %in% global_top_10_orders, Order, "Others")) %>%
    group_by(Order) %>%
    summarise(
      `04月` = sum(`04月`),
      `10月` = sum(`10月`),
      .groups = "drop"
    )
  
  seasons <- c("04月", "10月")
  for (mon in seasons) {
    # 提取单季节有效数据
    season_data <- df_plot %>% select(Order, TPM = all_of(mon)) %>% filter(TPM > 0)
    if (nrow(season_data) == 0) next
    
    # 严格按照全局因子的级别排序
    season_data <- season_data %>%
      mutate(Order = factor(Order, levels = final_order_levels)) %>%
      arrange(Order)
    
    # 计算相对百分比
    slices <- season_data$TPM
    pct <- round(slices / sum(slices) * 100, 1)
    
    # 饼图外侧只显示纯百分比符号，低于 1.5% 的隐藏
    labels_only_pct <- ifelse(pct >= 1.5, paste0(pct, "%"), "")
    plot_colors <- global_color_map[as.character(season_data$Order)]
    
    # 定义保存路径
    file_base <- paste0(output_dir, "/", enz, "_", mon, "_独立3D内嵌图例饼图")
    
    # 使用两阶段设备写入（PNG & PDF均应用 layout 分栏布局）
    for (f_type in c("png", "pdf")) {
      if (f_type == "png") {
        png(paste0(file_base, ".png"), width = 3200, height = 1800, res = 300, bg = "white")
      } else {
        pdf(paste0(file_base, ".pdf"), width = 11, height = 6.2, bg = "white")
      }
      
      # 将画布强制切分为左右两块（1号画饼图，2号画图例），权重比例为 7.5 : 2.5
      layout(matrix(c(1, 2), nrow = 1), widths = c(7.5, 2.5))
      
      # ---- 【第一步：在左侧区域1画 3D 饼图】 ----
      par(mar = c(2, 4, 5, 2)) # 收紧边距
      pie3D(slices, 
            labels = labels_only_pct, 
            edges = 350, 
            radius = 0.9, 
            height = 0.11, 
            theta = pi/4, 
            explode = 0.04, 
            col = plot_colors, 
            main = paste0(enz, " 核心物种演替 (", mon, ")"), 
            cex.main = 1.3, font.main = 2,
            labelcex = 0.9)
      
      # ---- 【第二步：在右侧区域2专门画图例】 ----
      # 清空右侧单独区域的坐标轴系，创建一个干净的透明背景层
      par(mar = c(0, 0, 0, 0)) 
      plot.new() 
      
      # 在这个完全干净不受干扰的独立空间里绘制大图例，用于绘制独立图例
      legend(x = "center",                  # 图例置于右侧区域中央
             legend = final_order_levels,   # 全局统一 10 个目 + Others
             fill = global_color_map,       # 全局统一颜色
             bty = "n",                     # 不显示边框
             title = "全局统一核心目\n(Order Level)", 
             title.font = 2,                # 标题加粗
             cex = 0.9,                     # 文字大小调节
             y.intersp = 1.4)               # 行间距加宽，调整行间距
      
      dev.off() # 关闭绘图设备
    }
  }
}

cat("\n========================================================")
cat("\n3D饼图及图例输出完成。")
cat("\n每张3D饼图右侧均已稳固、清晰地内嵌了【全局统一图例】。")
cat("\n结果请见：\n ", output_dir)
cat("\n========================================================\n")








