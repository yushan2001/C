# Cellulose analysis
#
# Modules:
#   1. Cellulase Z-score heatmap
#   2. Circular phylogenetic tree with 21 cellulose-related gene rings
#   3. Seasonal Order-cellulase chord diagrams
#
# Original paths, sheets, gene order, analysis settings, colors,
# plotting parameters, and output filenames are retained.
# ==============================================================================

# ===================== 加载包 =====================
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

# ===================== 1. 读取原始数据 =====================
file_path <- "C:/Users/余山小可爱/Desktop/纤维素原始代码数据.xlsx"
df_raw <- read_excel(file_path, sheet = "Sheet1")

# 数据结构检查
cat("===  原始数据核对 ===\n")
cat("数据维度：", nrow(df_raw), " 种纤维素酶 × ", ncol(df_raw)-1, " 个样本\n")
cat("样本列表：", paste(colnames(df_raw)[-1], collapse = ", "), "\n")
cat("纤维素酶原始顺序：", paste(df_raw$cellulase_type, collapse = ", "), "\n")
cat("\n")

# ===================== 2. 数据预处理 + Z-score标准化 =====================
# 设置行名并转换为数值矩阵
df_matrix <- df_raw %>%
  column_to_rownames("cellulase_type") %>%
  as.matrix()

# 缺失值填充为 0
df_matrix[is.na(df_matrix)] <- 0

# 按纤维素酶进行行 Z-score 标准化
# 公式：Z = (x - 行均值) / 行标准差
df_zscore <- t(apply(df_matrix, 1, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

# 转换为长格式
df_zscore_long <- as.data.frame(df_zscore) %>%
  rownames_to_column("cellulase_type") %>%
  pivot_longer(
    cols = -cellulase_type,
    names_to = "Sample",
    values_to = "Z_score"
  )

# 固定样本顺序
sample_order <- colnames(df_raw)[-1]
df_zscore_long$Sample <- factor(df_zscore_long$Sample, levels = sample_order)

# ===================== 固定纤维素酶原始顺序 =====================
# 使用原始数据行顺序
enzyme_original_order <- df_raw$cellulase_type
df_zscore_long$cellulase_type <- factor(df_zscore_long$cellulase_type, levels = enzyme_original_order)

# ===================== 3. 绘制标准化热图 =====================
heatmap_plot <- ggplot(
  df_zscore_long,
  aes(x = Sample, y = cellulase_type, fill = Z_score)
) +
  # 热图图层
  geom_tile(color = "white", linewidth = 0.2) +
  # 配色：低值为蓝色，中值为白色，高值为红色
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-3, 3),  # 限制Z值范围，避免极端值影响配色
    oob = squish,
    name = "Z-score"
  ) +
  # 坐标轴标签和标题
  labs(
    title = "Standardized Heatmap of Cellulase TPM",
    x = "Sample",
    y = "Cellulase Type"
  ) +
  # 主题设置
  theme_bw() +
  theme(
    # 标题设置
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    # 坐标轴设置
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    # 图例设置
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    # 网格线关闭
    panel.grid = element_blank()
  )

# 显示热图
print(heatmap_plot)

# ===================== 4. 保存图片 =====================
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


# ==============================================================================
# 纤维素相关基因环形进化树
# ==============================================================================

# 1. 加载依赖包
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

# ===================== 1. 参数与路径配置 =====================
file_path  <- "C:/Users/余山小可爱/Desktop/纤维素原始代码数据.xlsx"
tree_path  <- "C:/Users/余山小可爱/Desktop/gtdb_results.backbone.bac120.classify.tree"
output_dir <- "C:/Users/余山小可爱/Desktop/Output_Phylogeny_Pro/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 绘图布局与几何参数
tree_line_size <- 0.8     # 树枝线条粗细
tip_point_size <- 2.5     # 门分类 Tip 点大小
hilight_alpha  <- 0.3     # 背景高亮透明度

# 热图环形布局参数
strip_offset   <- 0.08    # 第一圈热图距离树顶端的起始距离
ring_width     <- 0.045   # 每圈热图的物理宽度
ring_gap       <- 0.006   # 圈与圈之间的空隙
border_size    <- 0.05    # 热图小方块的黑色边框粗细

# 21 个纤维素相关基因列，按内环到外环排列
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

# ===================== 2. 数据清洗与过滤 =====================
# 读取 Excel 并清除表头及 Bin 列的两端隐形空格
full_df <- read_excel(file_path, sheet = "Sheet2")
colnames(full_df) <- trimws(colnames(full_df)) 

full_df <- full_df %>%
  rename(Bin = `bin_name`) %>%  
  mutate(Bin = trimws(Bin))

# 读取进化树并保留与 Excel 匹配的样本
tree <- read.tree(tree_path)
keep_bins <- intersect(tree$tip.label, full_df$Bin)
filtered_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep_bins])

checkm_df <- full_df %>% 
  filter(Bin %in% filtered_tree$tip.label) %>%
  mutate(Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum))

# 转换数值类型，缺失值归零
for(col in gene_cols) {
  if(col %in% colnames(checkm_df)) {
    checkm_df[[col]] <- as.numeric(as.character(checkm_df[[col]]))
    checkm_df[[col]][is.na(checkm_df[[col]])] <- 0
  } else {
    warning(paste0(" 注意：在 Excel 中未找到列名：'", col, "'，请检查是否字打错了！"))
  }
}

# 获取基因拷贝数最大值
max_actual_val <- checkm_df %>% select(any_of(gene_cols)) %>% max(na.rm = TRUE)
if(is.na(max_actual_val) || max_actual_val == 0) max_actual_val <- 10

# ===================== 3. 基础树绘制与背景高亮 =====================
# 提取分类与门水平颜色映射
unique_phyla <- unique(checkm_df$Phylum)
phylum_color_map <- setNames(scales::hue_pal()(length(unique_phyla)), unique_phyla)

# 3.1 基础树
p_base <- ggtree(filtered_tree, layout = "circular", color = "black", size = tree_line_size)

# 获取树的最大半径
max_radius <- max(p_base$data$x)
tip_coords <- p_base$data %>% filter(isTip) %>% select(label, y)

current_plot <- p_base

# 3.2 计算高亮背景延伸长度
heatmap_total_width <- length(gene_cols) * (ring_width + ring_gap)
extend_length <- strip_offset + heatmap_total_width - ring_gap

# 门水平背景高亮
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

# ===================== 4. 绘制 21 层绝对数量热图环 =====================
for (i in seq_along(gene_cols)) {
  gene_name <- gene_cols[i]
  if(!gene_name %in% colnames(checkm_df)) next
  
  ring_data <- checkm_df %>% 
    select(Bin, val = !!sym(gene_name)) %>%
    inner_join(tip_coords, by = c("Bin" = "label")) %>%
    # 计算当前热图环的 X 坐标
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

# ===================== 5. 添加 Tip 点 (Phylum) =====================
# 添加门水平 Tip 点
final_plot <- current_plot + 
  new_scale_color() + 
  geom_tippoint(data = p_base$data %>% filter(isTip) %>% left_join(checkm_df, by = c("label" = "Bin")), 
                aes(x = x, y = y, color = Phylum), size = tip_point_size, alpha = 0.9, inherit.aes = FALSE) +
  scale_color_manual(values = phylum_color_map, name = "Phylum")

# ===================== 6. 导出结果与图例 =====================
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

# 移除主图图例并调整外圈边界
final_main <- final_plot + 
  theme(legend.position = "none") + 
  xlim(0, max_radius + strip_offset + heatmap_total_width + 0.4)

# 保存 PNG 和 PDF 主图
ggsave(paste0(output_dir, "Main_Tree_Cellulose_v3.png"), final_main, width = 18, height = 18, dpi = 600, device = png(type = "cairo"))
ggsave(paste0(output_dir, "Main_Tree_Cellulose_v3.pdf"), final_main, width = 18, height = 18, device = "pdf")

# 保存独立图例
ggsave(paste0(output_dir, "Legends_Cellulose_v3.pdf"), ggdraw(all_legend), width = 5, height = 10, device = "pdf")

# ============================================================
# 纤维素降解潜在贡献 Chord Diagram
#
# 三类经典 cellulases + Other cellulases
# ×
# TOP10 Orders + Others
#
# 三类经典纤维素酶：
#
# 1. Endoglucanase (EG)
# 2. Exoglucanase (ExG)
# 3. Beta-glucosidase (BGL)
#
# 其他所有纤维素相关酶/基因：
# Other cellulases
#
#
# 数据：
# Sheet2 = bin_name + cellulose-related gene/enzyme copy number
# Sheet3 = Genome + taxonomy + April/October TPM
#
#
# Potential contribution =
# MAG mean TPM × gene/enzyme copy number
#
#
# 最终显示：
#
# Order：
# TOP10 Orders + Others
#
# Cellulase：
# Endoglucanase (EG)
# Exoglucanase (ExG)
# Beta-glucosidase (BGL)
# Other cellulases
#
#
# 注意：
# 弦表示 potential contribution relationship，
# 不代表真实代谢流、物质流或微生物相互作用。
# ============================================================


# ============================================================
# 0. 加载包
# ============================================================

packages <- c(
  "readxl",
  "tidyverse",
  "circlize",
  "stringr"
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
library(stringr)


# ============================================================
# 1. 输入 / 输出路径
# ============================================================

input_file <-
  "C:/Users/余山小可爱/Desktop/纤维素.xlsx"


output_dir <-
  "C:/Users/余山小可爱/Desktop/纤维素_Chord_三类纤维素酶_EG_ExG_BGL_TOP10目"


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
# 是否将每个季节最终用于作图的全部弦
# 重新标准化为100%
# ------------------------------------------------------------

RENORMALIZE_SELECTED <- TRUE


# ------------------------------------------------------------
# Link透明度
# ------------------------------------------------------------

LINK_TRANSPARENCY <- 0.38


# ============================================================
# 3. 三类纤维素酶名称
# ============================================================

EG_NAME <-
  "Endoglucanase (EG)"


EXG_NAME <-
  "Exoglucanase (ExG)"


BGL_NAME <-
  "Beta-glucosidase (BGL)"


OTHER_NAME <-
  "Other cellulases"


# ============================================================
# 4. Order 配色
# ============================================================

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


# ============================================================
# 5. 三类纤维素酶 sector 配色
#
# 只影响外圈 sector，
# 中间弦仍然按照 Order 着色
# ============================================================

cellulase_colors_raw <- c(

  "Endoglucanase (EG)" =
    "#CFE8DD",

  "Exoglucanase (ExG)" =
    "#F5D9B8",

  "Beta-glucosidase (BGL)" =
    "#DDD5EE",

  "Other cellulases" =
    "#C5CDD2"

)


# ============================================================
# 6. 读取数据
# ============================================================

# ------------------------------------------------------------
# Sheet2
#
# bin_name +
# cellulose-related genes / enzymes
# ------------------------------------------------------------

df_gene <- read_excel(
  input_file,
  sheet = "Sheet2"
)


# ------------------------------------------------------------
# Sheet3
#
# Genome + Order + April/October TPM
# ------------------------------------------------------------

df_tpm <- read_excel(
  input_file,
  sheet = "Sheet3"
)


# ============================================================
# 7. 检查关键列
# ============================================================

if (!"bin_name" %in% colnames(df_gene)) {

  stop(
    " Sheet2 中没有找到 bin_name 列"
  )
}


if (!"Genome" %in% colnames(df_tpm)) {

  stop(
    " Sheet3 中没有找到 Genome 列"
  )
}


if (!"Order" %in% colnames(df_tpm)) {

  stop(
    " Sheet3 中没有找到 Order 列"
  )
}


# ============================================================
# 8. 整理纤维素基因表
# ============================================================

df_gene <- df_gene %>%

  rename(
    Genome = bin_name
  )


gene_cols <- setdiff(
  colnames(df_gene),
  "Genome"
)


cat(
  "\n Sheet2 中识别到 ",
  length(gene_cols),
  " 个纤维素相关酶/基因列\n"
)


# ============================================================
# 9. Copy number 转换为 numeric
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
# 10. 定义纤维素酶分类函数
#
# 分类依据：
#
# 根据列名判断：
#
# EG：
# EG
# Endoglucanase
# Endo-glucanase
# Endo_glucanase
#
# ExG：
# ExG
# Exoglucanase
# Exo-glucanase
# Exo_glucanase
#
# 注意：
# Cellobiohydrolase / CBH 不再并入 ExG，
# 它们会归入 Other cellulases。
#
# BGL：
# BGL
# Beta-glucosidase
# β-glucosidase
#
# 其余全部：
# Other cellulases
#
#
# 不自动用 GH family 推断，
# 因为某些 GH family 并不对应唯一酶功能。
# ============================================================

classify_cellulase <- function(x) {


  x2 <- tolower(x)


  x2 <- str_replace_all(
    x2,
    "β",
    "beta"
  )


  # ==========================================================
  # ExG
  #
  # 仅将明确标注为 ExG / Exoglucanase 的列归为 ExG。
  # CBH / Cellobiohydrolase 归入 Other cellulases。
  #
  # ExG 在 EG 之前判断，
  # 避免将 exoglucanase 误判为 EG。
  # ==========================================================

  if (
    str_detect(
      x2,
      "(^|[^a-z0-9])exg([^a-z0-9]|$)"
    ) |
      str_detect(
        x2,
        "exo[-_ ]?glucanase"
      )
  ) {

    return(
      EXG_NAME
    )

  }


  # ==========================================================
  # BGL
  # ==========================================================

  if (
    str_detect(
      x2,
      "(^|[^a-z0-9])bgl([^a-z0-9]|$)"
    ) |
      str_detect(
        x2,
        "beta[-_ ]?glucosidase"
      ) |
      str_detect(
        x2,
        "beta[-_ ]?d[-_ ]?glucosidase"
      )
  ) {

    return(
      BGL_NAME
    )

  }


  # ==========================================================
  # EG
  # ==========================================================

  if (
    str_detect(
      x2,
      "(^|[^a-z0-9])eg([^a-z0-9]|$)"
    ) |
      str_detect(
        x2,
        "endoglucanase"
      ) |
      str_detect(
        x2,
        "endo[-_ ]?glucanase"
      )
  ) {

    return(
      EG_NAME
    )

  }


  # ==========================================================
  # 其余
  # ==========================================================

  return(
    OTHER_NAME
  )

}


# ============================================================
# 11. 对 Sheet2 中所有酶/基因进行分类
# ============================================================

enzyme_classification <- tibble(

  Original_gene =
    gene_cols

) %>%

  mutate(

    Cellulase_class = vapply(

      Original_gene,

      classify_cellulase,

      FUN.VALUE =
        character(1)

    )

  )


cat(
  "\n========================================\n"
)


cat(
  "纤维素酶分类结果：\n"
)


print(
  enzyme_classification,
  n = Inf
)


# ============================================================
# 12. 保存分类表
#
# 保存酶分类表
# ============================================================

write.csv(

  enzyme_classification,

  file =
    paste0(
      output_dir,
      "/Cellulose_Enzyme_classification.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 13. 检查三个目标类别是否成功识别
# ============================================================

detected_classes <- unique(
  enzyme_classification$Cellulase_class
)


if (!EG_NAME %in% detected_classes) {

  warning(
    paste0(
      " 没有从 Sheet2 列名中识别到 ",
      EG_NAME,
      "。请检查 Cellulose_Enzyme_classification.csv。"
    )
  )
}


if (!EXG_NAME %in% detected_classes) {

  warning(
    paste0(
      " 没有从 Sheet2 列名中识别到 ",
      EXG_NAME,
      "。请检查 Cellulose_Enzyme_classification.csv。"
    )
  )
}


if (!BGL_NAME %in% detected_classes) {

  warning(
    paste0(
      " 没有从 Sheet2 列名中识别到 ",
      BGL_NAME,
      "。请检查 Cellulose_Enzyme_classification.csv。"
    )
  )
}


# ============================================================
# 14. 自动识别 April / October TPM 样本列
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
    " 未找到以04开头的April样本"
  )
}


if (length(month10_cols) == 0) {

  stop(
    " 未找到以10开头的October样本"
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
# 15. TPM列强制转 numeric
# ============================================================

df_tpm <- df_tpm %>%

  mutate(

    across(

      all_of(
        c(
          month04_cols,
          month10_cols
        )
      ),

      ~ suppressWarnings(
        as.numeric(.x)
      )

    )

  )


# ============================================================
# 16. 每个 MAG 计算季节平均 TPM
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


# NaN / Inf → 0
df_tpm$Month04[
  !is.finite(df_tpm$Month04)
] <- 0


df_tpm$Month10[
  !is.finite(df_tpm$Month10)
] <- 0


# ============================================================
# 17. 合并 taxonomy + TPM + cellulase
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
  "\n 合并后有效 MAG 数：",
  nrow(df_merge),
  "\n"
)


# ============================================================
# 18. 转成长格式
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
# 19. 把每个原始 Gene 映射到四个类别
# ============================================================

df_long <- df_long %>%

  left_join(

    enzyme_classification,

    by =
      c(
        "Gene" =
          "Original_gene"
      )

  )


if (any(is.na(df_long$Cellulase_class))) {

  stop(
    " 部分 Gene 无法匹配到 Cellulase_class"
  )
}


# ============================================================
# 20. 原始 Order × Gene 汇总
#
# 保留原始 Order × Gene 汇总
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
# 21. Order × Cellulase class 汇总
#
# EG / ExG / BGL / Other cellulases
# ============================================================

df_order_class_original <- df_long %>%

  group_by(
    Order,
    Cellulase_class
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
# 22. Order 排名
#
# 根据全部纤维素相关酶：
# April + October 总 potential contribution
# ============================================================

order_ranking <- df_order_class_original %>%

  group_by(Order) %>%

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

  arrange(
    desc(Total)
  )


# ============================================================
# 23. TOP10 Orders
# ============================================================

top_orders <- order_ranking %>%

  slice_head(
    n =
      TOP_ORDER_N
  ) %>%

  pull(Order)


cat(
  "\n========================================\n"
)


cat(
  "TOP ",
  TOP_ORDER_N,
  " Orders：\n",
  sep = ""
)


print(
  top_orders
)


# ============================================================
# 24. 第11名之后的 Orders → Others
# ============================================================

df_order_class <- df_order_class_original %>%

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
    Cellulase_class
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
    Order =
      Order_plot
  )


# ============================================================
# 25. 纤维素酶固定顺序
#
# 固定功能类别顺序
#

#
# EG → ExG → BGL → Other cellulases
# ============================================================

cellulase_levels_all <- c(

  EG_NAME,

  EXG_NAME,

  BGL_NAME,

  OTHER_NAME

)


# ============================================================
# 26. 只保留实际有 contribution 的类别
#
# 保留具有非零贡献的类别。
# ============================================================

class_totals <- df_order_class %>%

  group_by(
    Cellulase_class
  ) %>%

  summarise(

    Total =
      sum(
        April +
          October,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  )


cellulase_levels <- cellulase_levels_all[

  cellulase_levels_all %in%
    class_totals$Cellulase_class[
      class_totals$Total > 0
    ]

]


cat(
  "\n最终显示的 Cellulase classes：\n"
)


print(
  cellulase_levels
)


# ============================================================
# 27. 最终 Order 顺序
# ============================================================

order_levels <- c(
  top_orders,
  "Others"
)


# ============================================================
# 28. 保存 Cellulase class 总贡献
# ============================================================

cellulase_class_summary <- df_order_class %>%

  group_by(
    Cellulase_class
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

    Total =
      April +
      October,

    .groups =
      "drop"

  ) %>%

  mutate(

    Cellulase_class =
      factor(
        Cellulase_class,
        levels =
          cellulase_levels_all
      )

  ) %>%

  arrange(
    Cellulase_class
  )


write.csv(

  cellulase_class_summary,

  file =
    paste0(
      output_dir,
      "/Cellulose_Class_summary.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 29. 准备季节 Chord edge table
# ============================================================

prepare_chord_data <- function(
    data,
    season_col
) {


  df <- data %>%

    transmute(

      Order =
        Order,

      Cellulase =
        Cellulase_class,

      Value =
        .data[[season_col]]

    ) %>%

    filter(

      !is.na(Value),

      is.finite(Value),

      Value > 0

    )


  # ----------------------------------------------------------
  # 每季重新标准化为100%
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
  # Sector ID
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
          "Cellulase::",
          Cellulase
        )

    )


  return(
    df
  )

}


# ============================================================
# 30. April chord data
# ============================================================

April_chord <- prepare_chord_data(

  data =
    df_order_class,

  season_col =
    "April"

)


# ============================================================
# 31. October chord data
# ============================================================

October_chord <- prepare_chord_data(

  data =
    df_order_class,

  season_col =
    "October"

)


# ============================================================
# 32. 检查 Chord 数据
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
# 33. 清理 Order 名称
# ============================================================

clean_order <- function(x) {

  gsub(
    "^o__",
    "",
    x
  )

}


# ============================================================
# 34. 创建 Sector IDs
# ============================================================

order_sector_ids <- paste0(

  "Order::",

  order_levels

)


cellulase_sector_ids <- paste0(

  "Cellulase::",

  cellulase_levels

)


sector_order_all <- c(

  order_sector_ids,

  cellulase_sector_ids

)


# ============================================================
# 35. Order 配色
# ============================================================

if (
  length(top_orders) >
    length(order_palette)
) {

  stop(
    " Order配色数量不足"
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
# 36. Cellulase sector 配色
# ============================================================

cellulase_colors <- setNames(

  cellulase_colors_raw[
    cellulase_levels
  ],

  paste0(
    "Cellulase::",
    cellulase_levels
  )

)


# ============================================================
# 37. 所有 Sector 配色
# ============================================================

grid_colors <- c(

  order_colors,

  cellulase_colors

)


# ============================================================
# 38. Chord Diagram 绘图函数
# ============================================================

draw_cellulose_chord <- function(
    chord_data,
    title_text = ""
) {


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
  # Gap 设置
  #
  
  
  # ==========================================================

  gap_after <- ifelse(

    grepl(
      "^Cellulase::",
      current_order
    ),

    2.5,

    1.8

  )


  # ----------------------------------------------------------
  # Order组结束后大间距
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
    ] <- 10

  }


  # ----------------------------------------------------------
  # Cellulase组结束后大间距
  # ----------------------------------------------------------

  gap_after[
    length(gap_after)
  ] <- 10


  # ==========================================================
  # 当前sector配色
  # ==========================================================

  current_grid_colors <- grid_colors[
    current_order
  ]


  # ==========================================================
  # Circos参数
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
  # Link颜色
  #
  # Link 颜色由 Order 控制。
  
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
  # 绘制 Chord Diagram
  # ==========================================================

  chordDiagram(

    x =
      chord_data %>%

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
    # 无方向箭头
    # --------------------------------------------------------

    directional =
      0,

    link.sort =
      TRUE,

    link.largest.ontop =
      TRUE,

    link.border =
      NA,

    annotationTrack =
      "grid",

    preAllocateTracks =
      list(

        track.height =
          0.18

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
          0.60


        label_font <-
          3


      } else {


        # ====================================================
        # Cellulase class
        # ====================================================

        label <- gsub(
          "^Cellulase::",
          "",
          sector_id
        )


        label_cex <-
          0.65


        # 三类经典酶使用斜体
        # Other cellulases 使用普通字体
        if (
          label ==
            OTHER_NAME
        ) {

          label_font <-
            2

        } else {

          label_font <-
            3

        }

      }


      # ======================================================
      # 绘制标签
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


  circos.clear()

}


# ============================================================
# 39. April PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Cellulose_3Classes_TOP10Orders_April.png"
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


draw_cellulose_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


dev.off()


# ============================================================
# 40. April PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Cellulose_3Classes_TOP10Orders_April.pdf"
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


draw_cellulose_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


dev.off()


# ============================================================
# 41. October PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Cellulose_3Classes_TOP10Orders_October.png"
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


draw_cellulose_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 42. October PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Cellulose_3Classes_TOP10Orders_October.pdf"
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


draw_cellulose_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 43. April + October Combined PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Cellulose_3Classes_TOP10Orders_combined.png"
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


draw_cellulose_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


draw_cellulose_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 44. April + October Combined PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Cellulose_3Classes_TOP10Orders_combined.pdf"
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


draw_cellulose_chord(

  chord_data =
    April_chord,

  title_text =
    "April 2025"

)


draw_cellulose_chord(

  chord_data =
    October_chord,

  title_text =
    "October 2025"

)


dev.off()


# ============================================================
# 45. Order 图例
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
# 46. Order legend PNG
# ============================================================

png(

  filename =
    paste0(
      output_dir,
      "/Cellulose_Order_legend.png"
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
# 47. Order legend PDF
# ============================================================

pdf(

  file =
    paste0(
      output_dir,
      "/Cellulose_Order_legend.pdf"
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
# 48. 保存 April Chord data
# ============================================================

write.csv(

  April_chord,

  file =
    paste0(
      output_dir,
      "/Cellulose_April_Chord_data.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 49. 保存 October Chord data
# ============================================================

write.csv(

  October_chord,

  file =
    paste0(
      output_dir,
      "/Cellulose_October_Chord_data.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 50. 保存最终
# TOP10 Order + Others
# ×
# 3 classes + Other cellulases
# ============================================================

write.csv(

  df_order_class,

  file =
    paste0(
      output_dir,
      "/Cellulose_TOP10Orders_3Classes_grouped.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 51. 保存完整原始 Order × Gene
# ============================================================

write.csv(

  df_order_gene_original,

  file =
    paste0(
      output_dir,
      "/Cellulose_AllOrders_AllGenes_original.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 52. 保存 Order × Cellulase class 原始结果
# ============================================================

write.csv(

  df_order_class_original,

  file =
    paste0(
      output_dir,
      "/Cellulose_AllOrders_3Classes_original.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 53. 保存 Order 排名
# ============================================================

write.csv(

  order_ranking,

  file =
    paste0(
      output_dir,
      "/Cellulose_Order_ranking.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 54. 计算每个 Cellulase class 内
# 各 Order 的相对贡献
# ============================================================

relative_table <- df_order_class %>%

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
    Cellulase_class
  ) %>%

  mutate(

    Class_total =
      sum(
        Potential_contribution,
        na.rm = TRUE
      ),

    Relative_contribution =
      ifelse(

        Class_total > 0,

        Potential_contribution /
          Class_total *
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
      "/Cellulose_Order_Class_relative_contribution.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)


# ============================================================
# 55. 保存最终 sector 清单
# ============================================================

selection_table <- data.frame(

  Type =
    c(

      rep(
        "Order",
        length(order_levels)
      ),

      rep(
        "Cellulase",
        length(cellulase_levels)
      )

    ),

  Name =
    c(

      order_levels,

      cellulase_levels

    )

)


write.csv(

  selection_table,

  file =
    paste0(
      output_dir,
      "/Cellulose_Chord_final_sectors.csv"
    ),

  row.names =
    FALSE,

  fileEncoding =
    "UTF-8"

)
