# Lignin analysis
#
# Modules:
#   1. Order-level lignin-related protein abundance plots
#   2. Circular phylogenetic tree with 41 lignin-related gene rings
#   3. Seasonal Order-gene chord diagrams
#
# Original paths, sheet names, gene order, analysis settings, colors,
# plotting parameters, and output filenames are retained.
# ==============================================================================

# ===================== 1. 安装并加载必要包 =====================
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

# ===================== 2. 读取数据 =====================
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

# ===================== 3. 数据预处理 =====================
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

# ===================== 4. Order 配色 =====================
# 获取并排序 Order
unique_orders <- sort(unique(plot_data$Order))
n_orders <- length(unique_orders)

# 40 种浅色系
# 颜色顺序用于提高相邻类别辨识度
shuffled_light_palette <- c(
  "#FFB6C1", "#87CEFA", "#98FB98", "#DDA0DD", "#FFE4B5", "#7FFFD4", # 1-6:   浅粉、天蓝、薄荷绿、淡紫、奶油桔、碧绿
  "#F0E68C", "#F08080", "#B0E0E6", "#90EE90", "#EE82EE", "#FFD700", # 7-12:  淡金黄、浅珊瑚红、粉蓝、浅绿、兰花紫、明黄
  "#FFA07A", "#AFEEEE", "#E6E6FA", "#F5DEB3", "#9370DB", "#20B2AA", # 13-18: 浅鲑红、淡青绿、薰衣草紫、麦香色、中紫、浅海蓝
  "#FFC0CB", "#B0C4DE", "#8FBC8F", "#D8BFD8", "#FFF8DC", "#FA8072", # 19-24: 纯粉红、冰山蓝、暗海绿、蓟紫、米绢色、鲜鲑红
  "#E0FFFF", "#F4A460", "#9ACD32", "#F0F8FF", "#DEB887", "#E6A8D7", # 25-30: 浅青色、沙褐色、黄绿、爱丽丝蓝、硬木色、浅芭蕾粉
  "#7B68EE", "#FFDEAD", "#ADFF2F", "#BA55D3", "#FFE4E1", "#48D1CC", # 31-36: 暗灰蓝、白杏色、绿黄、中兰花紫、玫瑰白、中青绿
  "#F5F5DC", "#BC8F8F", "#ADD8E6", "#E0FFFE"                        # 37-40: 米色、褐玫瑰红、浅蓝、极淡青
)

# 根据 Order 数量截取颜色
final_colors <- shuffled_light_palette[1:n_orders]

# Order 颜色映射
order_color_map <- setNames(final_colors, unique_orders)

# ===================== 5. 创建保存文件夹 =====================
if (!dir.exists("Protein_Order_高辨识浅纯色图")) {
  dir.create("Protein_Order_高辨识浅纯色图")
}

# ===================== 6. 批量绘图函数 =====================
draw_plot <- function(prot_name) {
  sub_data <- plot_data %>% filter(`包含的Prc` == prot_name)
  
  if (nrow(sub_data) == 0 || all(sub_data$Abundance == 0)) {
    cat("跳过：Protein", prot_name, "无有效数据\n")
    return(NULL)
  }
  
  p <- ggplot(sub_data, aes(x = Time, y = Abundance, fill = Order)) +
    # 添加黑色边框以区分堆叠类别
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
  
  cat("已生成：", prot_name, "\n")
  return(p)
}

# ===================== 7. 批量生成所有图 =====================
all_proteins <- unique(plot_data$`包含的Prc`)
for (prot in all_proteins) {
  draw_plot(prot)
}

# ===================== 8. 统一图例 =====================
legend_df <- data.frame(
  Order = factor(unique_orders, levels = unique_orders),
  y = seq_along(unique_orders)
)

legend_plot <- ggplot(legend_df, aes(x = 1, y = y, fill = Order)) +
  geom_tile(color = "black", linewidth = 0.4) + # 图例添加黑色边框
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

# 根据 Order 数量调整图例高度
legend_height <- max(6, length(unique_orders) * 0.5 + 2)

ggsave("Protein_Order_高辨识浅纯色图/Order_统一图例.png",
       legend_plot, width = 6, height = legend_height, dpi = 300, bg = "white")
ggsave("Protein_Order_高辨识浅纯色图/Order_统一图例.pdf",
       legend_plot, width = 6, height = legend_height, device = "pdf", bg = "white")


# ==============================================================================
# 木质素相关基因环形进化树
# ==============================================================================

# 1. 加载依赖包
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

# ===================== 1. 参数与路径配置 =====================
file_path  <- "C:/Users/余山小可爱/Desktop/木质素原始代码数据.xlsx"
tree_path  <- "C:/Users/余山小可爱/Desktop/gtdb_results.backbone.bac120.classify.tree"
output_dir <- "C:/Users/余山小可爱/Desktop/Output_Phylogeny_Pro/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ----------------- 环形树布局参数 -----------------
tree_line_size <- 0.6       # 树分枝线条粗细
tip_point_size <- 1.5       # 树末端 Phylum 节点大小
strip_offset   <- 0.05      # 树最外层尖端到第一圈热图的空白间距
ring_width     <- 0.03      # 每圈热图的物理宽度
ring_gap       <- 0.005     # 热图圈与圈之间的微小空隙
border_size    <- 0.05      # 热图单元格的黑色边框粗细
hilight_alpha  <- 0.15      # 背景高亮扇区的透明度

# 41 个蛋白/基因列，按当前顺序排列
gene_cols <- c(
  "DesA", "DesB", "DesZ", "DypB", "DyPs", "FerA", "FerB/FerB2", "LACs", "LDS", 
  "LigA", "LigB", "LigC", "LigD", "LigEFG", "LigH", "LigI", "LigJ", "LigK", 
  "LigM", "LigV", "LigW", "LigX", "LigY", "LigZ", "PcaB", "PcaC", "pcaD", 
  "PcaG", "pcaH", "PcaI", "PcaJ", "PraA", "PraB", "PraC", "PraD", "PraE", 
  "PraF", "PraG", "PraH", "VanA", "VanB"
)

# ===================== 2. 数据处理 =====================
full_df <- read_excel(file_path, sheet = "Sheet2") %>%
  rename(Bin = `Bin ID`) %>%
  mutate(Bin = trimws(Bin))

tree <- read.tree(tree_path)
keep_bins <- intersect(tree$tip.label, full_df$Bin)
filtered_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep_bins])
checkm_df <- full_df %>% filter(Bin %in% filtered_tree$tip.label)

# 基因拷贝数转换为数值
for(col in gene_cols) {
  if(col %in% colnames(checkm_df)) {
    checkm_df[[col]] <- as.numeric(as.character(checkm_df[[col]]))
    checkm_df[[col]][is.na(checkm_df[[col]])] <- 0
  }
}

# 获取基因拷贝数最大值
max_actual_val <- checkm_df %>% 
  select(any_of(gene_cols)) %>% 
  max(na.rm = TRUE)

if(is.na(max_actual_val) || max_actual_val == 0) max_actual_val <- 10

# ===================== 3. 基础树绘制与动态半径获取 =====================
unique_phyla <- unique(checkm_df$Phylum)
phylum_color_map <- setNames(scales::hue_pal()(length(unique_phyla)), unique_phyla)

# 保留进化树原始比例
p_base <- ggtree(filtered_tree, layout = "circular", color = "black", size = tree_line_size)

# 获取进化树最大半径
max_radius <- max(p_base$data$x)
tip_coords <- p_base$data %>% filter(isTip) %>% select(label, y)

current_plot <- p_base

# ===================== 4. 背景高亮区域 =====================
# 计算 41 层热图总延伸长度
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

# ===================== 5. 绘制绝对数量热图环 =====================
for (i in seq_along(gene_cols)) {
  gene_name <- gene_cols[i]
  if(!gene_name %in% colnames(checkm_df)) next
  
  # 根据 max_radius 和 strip_offset 定位各热图环
  ring_data <- checkm_df %>% 
    select(Bin, val = !!sym(gene_name)) %>%
    inner_join(tip_coords, by = c("Bin" = "label")) %>%
    mutate(x = max_radius + strip_offset + (i-1) * (ring_width + ring_gap))
  
  current_plot <- current_plot +
    new_scale_fill() +
    geom_tile(data = ring_data, aes(x = x, y = y, fill = val), 
              width = ring_width, height = 1, color = "black", 
              size = border_size, inherit.aes = FALSE) +
    # 基因拷贝数渐变映射
    scale_fill_gradientn(colors = c("white", "#deebf7", "#3182ce", "#084594"),
                         values = scales::rescale(c(0, 0.1, max_actual_val * 0.5, max_actual_val)), 
                         limits = c(0, max_actual_val), oob = scales::squish, guide = "none")
}

# ===================== 6. 添加树尖端 Phylum 点 =====================
# 在树尖端添加 Phylum 分类颜色点
final_plot <- current_plot + 
  new_scale_color() + 
  geom_tippoint(data = p_base$data %>% filter(isTip) %>% left_join(checkm_df, by = c("label" = "Bin")), 
                aes(x = x, y = y, color = Phylum), size = tip_point_size, alpha = 0.9, inherit.aes = FALSE) +
  scale_color_manual(values = phylum_color_map, name = "Phylum")

# ===================== 7. 导出结果与图例 =====================
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

# 调整画布范围以避免最外圈截断
final_main <- final_plot + theme(legend.position = "none") + 
  xlim(0, max_radius + strip_offset + heatmap_total_width + 0.1)

# 保存 PNG 和 PDF 主图
ggsave(paste0(output_dir, "Tree_Absolute_Counts_Adaptive.png"), final_main, 
       width = 18, height = 18, dpi = 600, device = png(type = "cairo"))
ggsave(paste0(output_dir, "Tree_Absolute_Counts_Adaptive.pdf"), final_main, 
       width = 18, height = 18, device = "pdf")
ggsave(paste0(output_dir, "Legends_Final.pdf"), ggdraw(all_legend), width = 5, height = 10)


# ============================================================
# 木质素降解潜在贡献 Chord Diagram
# TOP20 lignin-related genes + Other genes
# ×
# TOP10 Orders + Others
#
# 数据：
# Sheet3 = bin_name + lignin-related gene copy number
# Sheet4 = Genome + taxonomy + April/October TPM
#
# Potential contribution =
# MAG mean TPM × gene copy number
#
# 最终显示：
#
# Order：
# TOP10 Orders + Others
#
# Gene：
# TOP20 Genes + Other genes
#
# 排名依据：
# April + October 总 potential contribution
#
# 注意：
# 弦表示 potential contribution relationship，
# 不代表真实代谢流、物质流或生态相互作用。
# ============================================================


# ============================================================
# 0. 加载包
# ============================================================

packages <- c(
  "readxl",
  "tidyverse",
  "circlize"
)

for (pkg in packages) {

  if (!require(pkg, character.only = TRUE)) {

    install.packages(
      pkg,
      dependencies = TRUE
    )

    library(
      pkg,
      character.only = TRUE
    )
  }
}

library(readxl)
library(tidyverse)
library(circlize)


# ============================================================
# 1. 输入 / 输出路径
# ============================================================

input_file <-
  "C:/Users/余山小可爱/Desktop/MAGs-lingin.xlsx"


output_dir <-
  "C:/Users/余山小可爱/Desktop/木质素_Chord_TOP20基因_TOP10目"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat(
  " 输出目录：",
  output_dir,
  "\n"
)


# ============================================================
# 2. 参数
# ============================================================

# ------------------------------------------------------------
# Order：
# TOP10 + Others
# ------------------------------------------------------------

TOP_ORDER_N <- 10


# ------------------------------------------------------------
# Gene：
# TOP20 + Other genes
# ------------------------------------------------------------

TOP_GENE_N <- 20


# ------------------------------------------------------------
# TRUE：
# 对最终绘图数据在每个季节重新归一化为 100%
# ------------------------------------------------------------

RENORMALIZE_SELECTED <- TRUE


# ------------------------------------------------------------
# 连线透明度
# ------------------------------------------------------------

LINK_TRANSPARENCY <- 0.38


# ============================================================
# 3. 配色
# ============================================================

# ------------------------------------------------------------
# TOP10 Order 配色
# ------------------------------------------------------------

order_palette <- c(

  "#E99AAA",
  "#F2B880",
  "#E6CF72",
  "#8CC7A1",
  "#77B8D7",
  "#A99AD9",
  "#D69CCB",
  "#8FA7D6",
  "#E59B83",
  "#76C9C3"

)


# ------------------------------------------------------------
# Others
# ------------------------------------------------------------

OTHERS_COLOR <- "#9E9E9E"


# ------------------------------------------------------------
# TOP20 Gene sector
# 统一浅灰蓝
# ------------------------------------------------------------

GENE_COLOR <- "#DDE6EB"


# ------------------------------------------------------------
# Other genes
# 稍深灰色，便于和TOP20区分
# ------------------------------------------------------------

OTHER_GENE_COLOR <- "#B8C2C8"


# ============================================================
# 4. 读取数据
# ============================================================

# ------------------------------------------------------------
# Sheet3：
# bin_name + 木质素相关基因
# ------------------------------------------------------------

df_gene <- read_excel(
  input_file,
  sheet = "Sheet3"
)


# ------------------------------------------------------------
# Sheet4：
# Genome + taxonomy + TPM
# ------------------------------------------------------------

df_tpm <- read_excel(
  input_file,
  sheet = "Sheet4"
)


# ============================================================
# 5. 检查关键列
# ============================================================

if (!"bin_name" %in% colnames(df_gene)) {

  stop(
    " Sheet3 中没有找到 bin_name 列"
  )
}


if (!"Genome" %in% colnames(df_tpm)) {

  stop(
    " Sheet4 中没有找到 Genome 列"
  )
}


if (!"Order" %in% colnames(df_tpm)) {

  stop(
    " Sheet4 中没有找到 Order 列"
  )
}


# ============================================================
# 6. 整理 Gene 表
# ============================================================

df_gene <- df_gene %>%

  rename(
    Genome = bin_name
  )


# ------------------------------------------------------------
# Genome 之外全部视为木质素相关基因
# ------------------------------------------------------------

gene_cols <- setdiff(
  colnames(df_gene),
  "Genome"
)


cat(
  " Sheet3 中共识别到 ",
  length(gene_cols),
  " 个木质素相关基因/酶\n"
)


# ============================================================
# 7. Gene copy number 转 numeric
# ============================================================

df_gene <- df_gene %>%

  mutate(

    across(

      all_of(gene_cols),

      ~ suppressWarnings(
        as.numeric(.x)
      )

    )

  ) %>%

  mutate(

    across(

      all_of(gene_cols),

      ~ replace_na(
        .x,
        0
      )

    )

  )


# ============================================================
# 8. 自动识别 April / October TPM
# ============================================================

month04_cols <- grep(
  "^04",
  colnames(df_tpm),
  value = TRUE
)


month10_cols <- grep(
  "^10",
  colnames(df_tpm),
  value = TRUE
)


if (length(month04_cols) == 0) {

  stop(
    " 没有找到04开头的April样本"
  )
}


if (length(month10_cols) == 0) {

  stop(
    " 没有找到10开头的October样本"
  )
}


cat(
  "\nApril samples:\n",
  paste(
    month04_cols,
    collapse = ", "
  ),
  "\n"
)


cat(
  "\nOctober samples:\n",
  paste(
    month10_cols,
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 9. 计算每个 MAG 的季节平均 TPM
# ============================================================

df_tpm <- df_tpm %>%

  mutate(

    Month04 = rowMeans(

      select(
        .,
        all_of(month04_cols)
      ),

      na.rm = TRUE

    ),

    Month10 = rowMeans(

      select(
        .,
        all_of(month10_cols)
      ),

      na.rm = TRUE

    )

  )


# ------------------------------------------------------------
# NaN / Inf → 0
# ------------------------------------------------------------

df_tpm$Month04[
  !is.finite(df_tpm$Month04)
] <- 0


df_tpm$Month10[
  !is.finite(df_tpm$Month10)
] <- 0


# ============================================================
# 10. 合并 taxonomy + TPM + gene
# ============================================================

df_merge <- df_tpm %>%

  inner_join(
    df_gene,
    by = "Genome"
  ) %>%

  filter(
    !is.na(Order),
    Order != ""
  )


cat(
  " 合并后有效 MAG 数：",
  nrow(df_merge),
  "\n"
)


# ============================================================
# 11. 转成长格式
# ============================================================

df_long <- df_merge %>%

  select(

    Genome,
    Order,
    Month04,
    Month10,
    all_of(gene_cols)

  ) %>%

  pivot_longer(

    cols =
      all_of(gene_cols),

    names_to =
      "Gene",

    values_to =
      "Copy_number"

  ) %>%

  mutate(

    Copy_number =
      replace_na(
        Copy_number,
        0
      ),

    April =
      Copy_number *
      Month04,

    October =
      Copy_number *
      Month10

  ) %>%

  filter(
    Copy_number > 0
  )


# ============================================================
# 12. 原始 Order × Gene 汇总
# ============================================================

df_order_gene_original <- df_long %>%

  group_by(
    Order,
    Gene
  ) %>%

  summarise(

    April =
      sum(
        April,
        na.rm = TRUE
      ),

    October =
      sum(
        October,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  )


# ============================================================
# 13. Order 排名
#
# 根据：
# 所有基因 × April + October 总贡献
# ============================================================

order_ranking <- df_order_gene_original %>%

  group_by(Order) %>%

  summarise(

    Total =
      sum(
        April +
        October,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  ) %>%

  arrange(
    desc(Total)
  )


# ============================================================
# 14. 提取 TOP10 Orders
# ============================================================

top_orders <- order_ranking %>%

  slice_head(
    n = TOP_ORDER_N
  ) %>%

  pull(Order)


cat(
  "\n====================================\n"
)


cat(
  "TOP ",
  TOP_ORDER_N,
  " Orders:\n",
  sep = ""
)


print(
  top_orders
)


# ============================================================
# 15. Order：
# TOP10之外全部合并为 Others
# ============================================================

df_order_gene_ordergroup <- df_order_gene_original %>%

  mutate(

    Order_plot = ifelse(

      Order %in%
        top_orders,

      Order,

      "Others"

    )

  ) %>%

  group_by(
    Order_plot,
    Gene
  ) %>%

  summarise(

    April =
      sum(
        April,
        na.rm = TRUE
      ),

    October =
      sum(
        October,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  ) %>%

  rename(
    Order = Order_plot
  )


# ============================================================
# 16. Gene 排名
#
# 根据：
# 所有 Order × April + October 总贡献
#
# Order 已合并为 TOP10 + Others，
# Gene 总贡献保持不变。
# ============================================================

gene_ranking <- df_order_gene_ordergroup %>%

  group_by(Gene) %>%

  summarise(

    Total =
      sum(
        April +
        October,
        na.rm = TRUE
      ),

    April =
      sum(
        April,
        na.rm = TRUE
      ),

    October =
      sum(
        October,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  ) %>%

  filter(
    Total > 0
  ) %>%

  arrange(
    desc(Total)
  )


# ============================================================
# 17. 提取 TOP20 Genes
# ============================================================

top_genes <- gene_ranking %>%

  slice_head(
    n = TOP_GENE_N
  ) %>%

  pull(Gene)


cat(
  "\nTOP ",
  TOP_GENE_N,
  " Genes:\n",
  sep = ""
)


print(
  top_genes
)


# ============================================================
# 18. 统计其余 Gene
# ============================================================

other_genes <- gene_ranking %>%

  filter(
    !Gene %in%
      top_genes
  ) %>%

  pull(Gene)


cat(
  "\nTOP20之外被合并为 Other genes 的基因数：",
  length(other_genes),
  "\n"
)


# ============================================================
# 19. Gene：
#
# TOP20 单独保留
# 第21名以后全部合并为 Other genes
#

# ============================================================

df_order_gene <- df_order_gene_ordergroup %>%

  mutate(

    Gene_plot = ifelse(

      Gene %in%
        top_genes,

      Gene,

      "Other genes"

    )

  ) %>%

  group_by(
    Order,
    Gene_plot
  ) %>%

  summarise(

    April =
      sum(
        April,
        na.rm = TRUE
      ),

    October =
      sum(
        October,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  ) %>%

  rename(
    Gene = Gene_plot
  )


# ============================================================
# 20. 最终 Gene 顺序
#
# TOP20 +
# Other genes（如果确实存在）
# ============================================================

if (length(other_genes) > 0) {

  gene_levels <- c(
    top_genes,
    "Other genes"
  )

} else {

  gene_levels <- top_genes

}


# ============================================================
# 21. 最终 Order 顺序
#
# TOP10 +
# Others
# ============================================================

order_levels <- c(
  top_orders,
  "Others"
)


# ============================================================
# 22. 检查最终数据
# ============================================================

cat(
  "\n最终 Order 数量：",
  length(order_levels),
  "\n"
)


cat(
  "最终 Gene sector 数量：",
  length(gene_levels),
  "\n"
)


cat(
  "\n最终 Gene sectors：\n"
)


print(
  gene_levels
)


# ============================================================
# 23. 准备季节 Chord edge table
# ============================================================

prepare_chord_data <- function(
    data,
    season_col
) {


  df <- data %>%

    transmute(

      Order =
        Order,

      Gene =
        Gene,

      Value =
        .data[[season_col]]

    ) %>%

    filter(

      !is.na(Value),

      is.finite(Value),

      Value > 0

    )


  # ----------------------------------------------------------
  # 每季重新归一化到100%
  # ----------------------------------------------------------

  if (RENORMALIZE_SELECTED) {


    total_value <- sum(
      df$Value,
      na.rm = TRUE
    )


    if (total_value > 0) {


      df <- df %>%

        mutate(

          Value =
            Value /
            total_value *
            100

        )

    }

  }


  # ----------------------------------------------------------
  # 创建唯一 sector ID
  # ----------------------------------------------------------

  df <- df %>%

    mutate(

      from =
        paste0(
          "Order::",
          Order
        ),

      to =
        paste0(
          "Gene::",
          Gene
        )

    )


  return(
    df
  )

}


# ============================================================
# 24. April chord data
# ============================================================

April_chord <- prepare_chord_data(

  data =
    df_order_gene,

  season_col =
    "April"

)


# ============================================================
# 25. October chord data
# ============================================================

October_chord <- prepare_chord_data(

  data =
    df_order_gene,

  season_col =
    "October"

)


# ============================================================
# 26. 检查 Chord 数据
# ============================================================

cat(
  "\nApril edges：",
  nrow(April_chord),
  "\n"
)


cat(
  "October edges：",
  nrow(October_chord),
  "\n"
)


if (nrow(April_chord) == 0) {

  stop(
    " April 无可绘制数据"
  )
}


if (nrow(October_chord) == 0) {

  stop(
    " October 无可绘制数据"
  )
}


# ============================================================
# 27. 清理 Order 名字
# ============================================================

clean_order <- function(x) {

  gsub(
    "^o__",
    "",
    x
  )

}


# ============================================================
# 28. Sector IDs
# ============================================================

order_sector_ids <- paste0(

  "Order::",

  order_levels

)


gene_sector_ids <- paste0(

  "Gene::",

  gene_levels

)


sector_order_all <- c(

  order_sector_ids,

  gene_sector_ids

)


# ============================================================
# 29. Order 配色
# ============================================================

if (
  length(top_orders) >
    length(order_palette)
) {

  stop(
    " Order颜色不足"
  )
}


top_order_colors <- setNames(

  order_palette[
    seq_along(
      top_orders
    )
  ],

  paste0(
    "Order::",
    top_orders
  )

)


others_color <- setNames(

  OTHERS_COLOR,

  "Order::Others"

)


order_colors <- c(

  top_order_colors,

  others_color

)


# ============================================================
# 30. Gene sector 配色
#
# TOP20：浅灰蓝
# Other genes：稍深灰色
# ============================================================

top_gene_colors <- setNames(

  rep(
    GENE_COLOR,
    length(top_genes)
  ),

  paste0(
    "Gene::",
    top_genes
  )

)


# ------------------------------------------------------------
# 如果存在 Other genes
# ------------------------------------------------------------

if (length(other_genes) > 0) {

  other_gene_color <- setNames(

    OTHER_GENE_COLOR,

    "Gene::Other genes"

  )


  gene_colors <- c(

    top_gene_colors,

    other_gene_color

  )

} else {

  gene_colors <- top_gene_colors

}


# ============================================================
# 31. 全部 sector 配色
# ============================================================

grid_colors <- c(

  order_colors,

  gene_colors

)


# ============================================================
# 32. Chord 绘图函数
# ============================================================

draw_lignin_chord <- function(
    chord_data,
    title_text = ""
) {


  # ----------------------------------------------------------
  # 清除上一张图
  # ----------------------------------------------------------

  circos.clear()


  # ==========================================================
  # 当前季节实际存在的 sectors
  # ==========================================================

  current_sectors <- unique(

    c(
      chord_data$from,
      chord_data$to
    )

  )


  current_order <- sector_order_all[
    sector_order_all %in%
      current_sectors
  ]


  # ==========================================================
  # 当前 Order sectors
  # ==========================================================

  current_order_sectors <- current_order[

    grepl(
      "^Order::",
      current_order
    )

  ]


  # ==========================================================
  # 当前 Gene sectors
  # ==========================================================

  current_gene_sectors <- current_order[

    grepl(
      "^Gene::",
      current_order
    )

  ]


  # ==========================================================
  # 设置 gap
  #
  # Order：1.8°
  # Gene：1.0°
  #
  # Order与Gene之间：
  # 8°
  # ==========================================================

  gap_after <- ifelse(

    grepl(
      "^Gene::",
      current_order
    ),

    1.0,

    1.8

  )


  # ----------------------------------------------------------
  # Order 组结束
  # ----------------------------------------------------------

  if (
    length(current_order_sectors) > 0
  ) {


    last_order_position <- match(

      tail(
        current_order_sectors,
        1
      ),

      current_order

    )


    gap_after[
      last_order_position
    ] <- 8

  }


  # ----------------------------------------------------------
  # Gene 组结束
  # ----------------------------------------------------------

  gap_after[
    length(gap_after)
  ] <- 8


  # ==========================================================
  # 当前 sector 颜色
  # ==========================================================

  current_grid_colors <- grid_colors[
    current_order
  ]


  # ==========================================================
  # circos 参数
  # ==========================================================

  circos.par(

    start.degree =
      90,

    gap.after =
      gap_after,

    track.margin =
      c(
        0.003,
        0.003
      ),

    cell.padding =
      c(
        0,
        0,
        0,
        0
      ),

    points.overflow.warning =
      FALSE

  )


  # ==========================================================
  # Link 颜色
  #
  # 每条弦继承来源 Order 的颜色
  # ==========================================================

  link_colors <- sapply(

    chord_data$from,

    function(x) {


      adjustcolor(

        order_colors[
          x
        ],

        alpha.f =
          1 -
          LINK_TRANSPARENCY

      )

    }

  )


  # ==========================================================
  # 正式绘图
  # ==========================================================

  chordDiagram(

    x = chord_data %>%

      select(
        from,
        to,
        Value
      ),

    order =
      current_order,

    grid.col =
      current_grid_colors,

    col =
      link_colors,

    transparency =
      0,

    # --------------------------------------------------------
    # 无方向
    # --------------------------------------------------------

    directional =
      0,

    # --------------------------------------------------------
    # Link排序
    # --------------------------------------------------------

    link.sort =
      TRUE,

    # --------------------------------------------------------
    # 大Link置顶
    # --------------------------------------------------------

    link.largest.ontop =
      TRUE,

    # --------------------------------------------------------
    # 不画Link边界
    # --------------------------------------------------------

    link.border =
      NA,

    # --------------------------------------------------------
    # 外圈
    # --------------------------------------------------------

    annotationTrack =
      "grid",

    # --------------------------------------------------------
    # 标签轨道
    # --------------------------------------------------------

    preAllocateTracks =
      list(

        track.height =
          0.16

      )

  )


  # ==========================================================
  # 外圈标签
  # ==========================================================

  circos.trackPlotRegion(

    track.index =
      1,

    bg.border =
      NA,

    panel.fun = function(
      x,
      y
    ) {


      sector_id <-
        get.cell.meta.data(
          "sector.index"
        )


      xlim <-
        get.cell.meta.data(
          "xlim"
        )


      ylim <-
        get.cell.meta.data(
          "ylim"
        )


      # ======================================================
      # Order
      # ======================================================

      if (
        grepl(
          "^Order::",
          sector_id
        )
      ) {


        label <- gsub(
          "^Order::o__",
          "",
          sector_id
        )


        label <- gsub(
          "^Order::",
          "",
          label
        )


        label_cex <-
          0.58


        # 微生物Order使用斜体
        label_font <-
          3


      } else {


        # ====================================================
        # Gene
        # ====================================================

        label <- gsub(
          "^Gene::",
          "",
          sector_id
        )


        # ----------------------------------------------------
        # TOP20 Gene
        # ----------------------------------------------------

        if (
          label != "Other genes"
        ) {

          label_cex <-
            0.52

          label_font <-
            3

        } else {


          # --------------------------------------------------
          # Other genes
          # --------------------------------------------------

          label_cex <-
            0.55

          label_font <-
            2

        }

      }


      # ======================================================
      # 绘制文字
      # ======================================================

      circos.text(

        x =
          mean(xlim),

        y =
          ylim[1] +
          0.12,

        labels =
          label,

        facing =
          "clockwise",

        niceFacing =
          TRUE,

        adj =
          c(
            0,
            0.5
          ),

        cex =
          label_cex,

        col =
          "black",

        font =
          label_font

      )

    }

  )


  # ==========================================================
  # 标题
  # ==========================================================

  title(

    main =
      title_text,

    line =
      -1,

    cex.main =
      1.35,

    font.main =
      2,

    family =
      "serif"

  )


  # ----------------------------------------------------------
  # 清除 circos
  # ----------------------------------------------------------

  circos.clear()

}


# ============================================================
# 33. April PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Lignin_TOP20Genes_TOP10Orders_April.png"
    ),

  width =
    3600,

  height =
    3600,

  res =
    300,

  bg =
    "white"

)


par(
  mar =
    c(
      2,
      2,
      3,
      2
    )
)


draw_lignin_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


dev.off()


# ============================================================
# 34. April PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Lignin_TOP20Genes_TOP10Orders_April.pdf"
    ),

  width =
    11,

  height =
    11,

  family =
    "serif"

)


par(
  mar =
    c(
      2,
      2,
      3,
      2
    )
)


draw_lignin_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


dev.off()


# ============================================================
# 35. October PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Lignin_TOP20Genes_TOP10Orders_October.png"
    ),

  width =
    3600,

  height =
    3600,

  res =
    300,

  bg =
    "white"

)


par(
  mar =
    c(
      2,
      2,
      3,
      2
    )
)


draw_lignin_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 36. October PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Lignin_TOP20Genes_TOP10Orders_October.pdf"
    ),

  width =
    11,

  height =
    11,

  family =
    "serif"

)


par(
  mar =
    c(
      2,
      2,
      3,
      2
    )
)


draw_lignin_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 37. April + October Combined PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Lignin_TOP20Genes_TOP10Orders_combined.png"
    ),

  width =
    7000,

  height =
    3500,

  res =
    300,

  bg =
    "white"

)


par(

  mfrow =
    c(
      1,
      2
    ),

  mar =
    c(
      2,
      2,
      3,
      2
    )

)


draw_lignin_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


draw_lignin_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 38. April + October Combined PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Lignin_TOP20Genes_TOP10Orders_combined.pdf"
    ),

  width =
    22,

  height =
    11,

  family =
    "serif"

)


par(

  mfrow =
    c(
      1,
      2
    ),

  mar =
    c(
      2,
      2,
      3,
      2
    )

)


draw_lignin_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


draw_lignin_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 39. 单独输出 Order 图例
# ============================================================

legend_order_names <- c(

  clean_order(
    top_orders
  ),

  "Others"

)


legend_order_colors <- c(

  order_palette[
    seq_along(
      top_orders
    )
  ],

  OTHERS_COLOR

)


# ============================================================
# 40. Order 图例 PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Lignin_Order_legend.png"
    ),

  width =
    1800,

  height =
    2400,

  res =
    300,

  bg =
    "white"

)


par(
  mar =
    c(
      0,
      0,
      0,
      0
    )
)


plot.new()


legend(

  "center",

  legend =
    legend_order_names,

  fill =
    legend_order_colors,

  border =
    NA,

  bty =
    "n",

  title =
    "Order",

  cex =
    1.1,

  y.intersp =
    1.25

)


dev.off()


# ============================================================
# 41. Order 图例 PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Lignin_Order_legend.pdf"
    ),

  width =
    5,

  height =
    7,

  family =
    "serif"

)


par(
  mar =
    c(
      0,
      0,
      0,
      0
    )
)


plot.new()


legend(

  "center",

  legend =
    legend_order_names,

  fill =
    legend_order_colors,

  border =
    NA,

  bty =
    "n",

  title =
    "Order",

  cex =
    1,

  y.intersp =
    1.25

)


dev.off()


# ============================================================
# 42. 保存 April Chord 数据
# ============================================================

write.csv(

  April_chord,

  file =
    paste0(
      output_dir,
      "/Lignin_April_Chord_data.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 43. 保存 October Chord 数据
# ============================================================

write.csv(

  October_chord,

  file =
    paste0(
      output_dir,
      "/Lignin_October_Chord_data.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 44. 保存最终 TOP10 Order + Others
# ×
# TOP20 Gene + Other genes 数据
# ============================================================

write.csv(

  df_order_gene,

  file =
    paste0(
      output_dir,
      "/Lignin_TOP10Orders_TOP20Genes_grouped.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 45. 保存完整原始 Order × Gene 数据
# ============================================================

write.csv(

  df_order_gene_original,

  file =
    paste0(
      output_dir,
      "/Lignin_AllOrders_AllGenes_original.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 46. 保存 Order 排名
# ============================================================

write.csv(

  order_ranking,

  file =
    paste0(
      output_dir,
      "/Lignin_Order_ranking.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 47. 保存 Gene 排名
# ============================================================

write.csv(

  gene_ranking,

  file =
    paste0(
      output_dir,
      "/Lignin_Gene_ranking.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 48. 保存 TOP20 Gene 清单
# ============================================================

top_gene_table <- gene_ranking %>%

  mutate(

    Rank =
      row_number(),

    Group =
      ifelse(
        Rank <= TOP_GENE_N,
        "Top20",
        "Other genes"
      )

  )


write.csv(

  top_gene_table,

  file =
    paste0(
      output_dir,
      "/Lignin_Gene_TOP20_classification.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 49. 保存最终 sector 清单
# ============================================================

selection_table <- data.frame(

  Type =
    c(

      rep(
        "Order",
        length(order_levels)
      ),

      rep(
        "Gene",
        length(gene_levels)
      )

    ),

  Name =
    c(

      order_levels,

      gene_levels

    )

)


write.csv(

  selection_table,

  file =
    paste0(
      output_dir,
      "/Lignin_Chord_final_sectors.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 50. 保存 Order 配色
# ============================================================

order_color_table <- data.frame(

  Order =
    legend_order_names,

  Color =
    legend_order_colors

)


write.csv(

  order_color_table,

  file =
    paste0(
      output_dir,
      "/Lignin_Order_colors.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 51. 计算最终各 Order × Gene 的季节内相对贡献
#
# 计算季节内相对贡献
# ============================================================

relative_table <- df_order_gene %>%

  pivot_longer(

    cols =
      c(
        April,
        October
      ),

    names_to =
      "Season",

    values_to =
      "Potential_contribution"

  ) %>%

  group_by(
    Season,
    Gene
  ) %>%

  mutate(

    Gene_total =
      sum(
        Potential_contribution,
        na.rm = TRUE
      ),

    Relative_contribution =
      ifelse(

        Gene_total > 0,

        Potential_contribution /
          Gene_total *
          100,

        0

      )

  ) %>%

  ungroup()


write.csv(

  relative_table,

  file =
    paste0(
      output_dir,
      "/Lignin_Order_Gene_relative_contribution.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 52. 单独计算 Other genes 总贡献
#
# 用于检查聚合后是否占比过大
# ============================================================

other_gene_summary <- df_order_gene %>%

  filter(
    Gene ==
      "Other genes"
  ) %>%

  summarise(

    April =
      sum(
        April,
        na.rm = TRUE
      ),

    October =
      sum(
        October,
        na.rm = TRUE
      )

  )


write.csv(

  other_gene_summary,

  file =
    paste0(
      output_dir,
      "/Lignin_OtherGenes_summary.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)
