# Environmental factor analysis for alpine lake samples
# PCA, Spearman correlation matrix, and NMDS/PERMANOVA
# Input: C:/Users/余山小可爱/Desktop/环境因子.xlsx (Sheet1)
# ----------------------------------------------------------------------

# Packages
# ----------------------------------------------------------------------

if (!require(readxl)) install.packages("readxl")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(vegan)) install.packages("vegan")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(GGally)) install.packages("GGally")
if (!require(ggpubr)) install.packages("ggpubr")

library(readxl)
library(tidyverse)
library(vegan)
library(ggplot2)
library(GGally)
library(ggpubr)


# ----------------------------------------------------------------------
# Input data
# ----------------------------------------------------------------------

file_path <- "C:/Users/余山小可爱/Desktop/环境因子.xlsx"

raw_data <- read_excel(
  file_path,
  sheet = "Sheet1"
)


# ----------------------------------------------------------------------
# Part I: PCA

# ----------------------------------------------------------------------


# ===================== 环境因子 PCA · 单张图 · 分季节色块 =====================


# ----------------------------------------------------------------------
# Prepare PCA data
# ----------------------------------------------------------------------

env <- raw_data


env_df <- env %>%
  column_to_rownames("指标") %>%
  t() %>%
  as.data.frame() %>%
  na.omit()


# ----------------------------------------------------------------------
# Sample groups
# 04 = 202504
# 10 = 202510
# ----------------------------------------------------------------------

env_df$Group <- ifelse(
  grepl("^04", rownames(env_df)),
  "202504",
  "202510"
)


# ----------------------------------------------------------------------
# PCA
# Z-score标准化 + rda PCA
# ----------------------------------------------------------------------

pca_env <- rda(
  scale(
    env_df[, 1:ncol(env_df)-1]
  )
)


# ----------------------------------------------------------------------
# PCA site scores
# ----------------------------------------------------------------------

site <- as.data.frame(
  scores(
    pca_env,
    display = "sites",
    scaling = 2
  )
)


colnames(site) <- c(
  "PC1",
  "PC2"
)


site$Group <- env_df$Group


# ----------------------------------------------------------------------
# Environmental variable scores
# ----------------------------------------------------------------------

arrow <- as.data.frame(
  scores(
    pca_env,
    display = "species",
    scaling = 2
  )
)


arrow$Var <- rownames(arrow)


colnames(arrow)[1:2] <- c(
  "PC1",
  "PC2"
)


# ----------------------------------------------------------------------
# Explained variance
# ----------------------------------------------------------------------

expl <- summary(pca_env)$cont$importance[2, 1:2]


lab1 <- paste0(
  "PC1 (",
  round(expl[1] * 100, 2),
  "%)"
)


lab2 <- paste0(
  "PC2 (",
  round(expl[2] * 100, 2),
  "%)"
)


# ----------------------------------------------------------------------
# PCA plot
# ----------------------------------------------------------------------

p <- ggplot() +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.3
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.3
  ) +
  
  geom_point(
    data = site,
    aes(
      PC1,
      PC2,
      fill = Group
    ),
    shape = 22,
    size = 4,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "202504" = "#FFB3DD",
      "202510" = "#ACFFFF"
    )
  ) +
  
  geom_segment(
    data = arrow,
    aes(
      0,
      0,
      xend = PC1,
      yend = PC2
    ),
    arrow = grid::arrow(
      length = grid::unit(
        0.25,
        "cm"
      )
    ),
    color = "black"
  ) +
  
  geom_text(
    data = arrow,
    aes(
      PC1,
      PC2,
      label = Var
    ),
    vjust = -0.7
  ) +
  
  labs(
    x = lab1,
    y = lab2,
    fill = "Season"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    axis.line = element_line()
  )


# ----------------------------------------------------------------------
# Save PCA plot
# ----------------------------------------------------------------------

print(p)


ggsave(
  "C:/Users/余山小可爱/Desktop/环境因子_PCA_单张分组图.png",
  p,
  width = 9,
  height = 6,
  dpi = 300
)


ggsave(
  "C:/Users/余山小可爱/Desktop/环境因子_PCA_单张分组图.pdf",
  p,
  width = 9,
  height = 6
)


# ----------------------------------------------------------------------
# Part I: PCAI

# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# 高山湖泊环境因子季节演替矩阵图
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Prepare correlation data
# ----------------------------------------------------------------------

env_data <- raw_data %>%
  column_to_rownames(var = "指标") %>%
  t() %>%
  as.data.frame()


# ----------------------------------------------------------------------
# Convert environmental variables to numeric
# ----------------------------------------------------------------------

env_data[] <- lapply(
  env_data,
  function(x)
    as.numeric(
      as.character(x)
    )
)


# ----------------------------------------------------------------------
# Season grouping
# ----------------------------------------------------------------------

env_data$Season <- ifelse(
  substr(
    rownames(env_data),
    1,
    2
  ) == "04",
  "April",
  "October"
)


env_data$Season <- factor(
  env_data$Season,
  levels = c(
    "October",
    "April"
  )
)


# ----------------------------------------------------------------------
# 修改点：
# 每个环境因子的对角线直接显示
# October vs April 箱线图
# Wilcoxon test
# 默认非配对
# ----------------------------------------------------------------------

my_diag_boxplot <- function(
    data,
    mapping,
    ...
){
  
  # ------------------------------------------------------------
  # 获取当前对角线对应的环境因子名称
  # ------------------------------------------------------------
  
  current_var <- rlang::as_name(
    mapping$x
  )
  
  
  # ------------------------------------------------------------
  #
  # x = Season
  # y = 当前环境变量
  #
  # ------------------------------------------------------------
  
  box_data <- data.frame(
    Season = data$Season,
    Value = data[[current_var]]
  )
  
  
  # ------------------------------------------------------------
  # 保持October在左、April在右
  # ------------------------------------------------------------
  
  box_data$Season <- factor(
    box_data$Season,
    levels = c(
      "October",
      "April"
    )
  )
  
  
  # ------------------------------------------------------------
  # 箱线图
  # ------------------------------------------------------------
  
  ggplot(
    box_data,
    aes(
      x = Season,
      y = Value,
      fill = Season,
      color = Season
    )
  ) +
    
    geom_boxplot(
      alpha = 0.65,
      outlier.shape = NA,
      width = 0.55
    ) +
    
    geom_jitter(
      width = 0.1,
      size = 0.8,
      alpha = 0.6
    ) +
    
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      label.x = 1.5
    ) +
    
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = 7
      ),
      axis.text.y = element_text(
        size = 7
      )
    )
}


# ----------------------------------------------------------------------
# October = #ACFFFF
# April   = #FFB3DD
# ----------------------------------------------------------------------

color_oct <- "#55C3DD"

color_apr <- "#E982AE"


# ----------------------------------------------------------------------
# Environmental variables
# Season仍然只是分组变量
# 不作为矩阵的一行/一列
# ----------------------------------------------------------------------

cols_to_plot <- setdiff(
  names(env_data),
  "Season"
)


# ----------------------------------------------------------------------
# Correlation matrix
# 1. 对角线 densityDiag
#    -> October / April箱线图
# 2. 两季颜色加深
# 其余保持：
# upper:
# Spearman相关系数
# lower:
# lm拟合
# Wilcoxon:
# wilcox.test
# ----------------------------------------------------------------------

plot_matrix_final <- ggpairs(
  
  env_data,
  
  # ------------------------------------------------------------
  #
  # ------------------------------------------------------------
  
  columns = match(
    cols_to_plot,
    names(env_data)
  ),
  
  
  mapping = aes(
    color = Season,
    fill = Season
  ),
  
  
  upper = list(
    
    # ----------------------------------------------------------
    #
    # 数值变量之间
    # Spearman相关系数
    # ----------------------------------------------------------
    
    continuous = wrap(
      "cor",
      method = "spearman",
      size = 3.5
    )
    
  ),
  
  
  lower = list(
    
    # ----------------------------------------------------------
    #
    # 数值变量之间
    # lm拟合曲线
    # ----------------------------------------------------------
    
    continuous = wrap(
      "smooth",
      alpha = 0.3,
      size = 0.7,
      method = "lm"
    )
    
  ),
  
  
  diag = list(
    
    # ----------------------------------------------------------
    #
    # densityDiag
    #
    # October / April箱线图
    # ----------------------------------------------------------
    
    continuous = my_diag_boxplot
    
  )
  
) +
  
  
# ----------------------------------------------------------------------
  
  scale_color_manual(
    values = c(
      "October" = color_oct,
      "April" = color_apr
    )
  ) +
  
  
  scale_fill_manual(
    values = c(
      "October" = color_oct,
      "April" = color_apr
    )
  ) +
  
  
# ----------------------------------------------------------------------
  
  theme_bw() +
  
  
  theme(
    
    strip.background =
      element_rect(
        fill = "white"
      ),
    
    strip.text =
      element_text(
        size = 9,
        face = "bold"
      ),
    
    axis.text =
      element_text(
        size = 7
      ),
    
    panel.grid.major =
      element_blank(),
    
    legend.position =
      "bottom"
  )


# ----------------------------------------------------------------------
# Display correlation matrix
# ----------------------------------------------------------------------

print(
  plot_matrix_final
)


# ----------------------------------------------------------------------
# Save correlation matrix
# ----------------------------------------------------------------------

ggsave(
  
  "C:/Users/余山小可爱/Desktop/最终完善版_Spearman相关性矩阵.pdf",
  
  plot = plot_matrix_final,
  
  width = 15,
  
  height = 13
)


ggsave(
  
  "C:/Users/余山小可爱/Desktop/最终完善版_Spearman相关性矩阵.png",
  
  plot = plot_matrix_final,
  
  width = 15,
  
  height = 13,
  
  dpi = 300,
  
  bg = "white"
)


# ----------------------------------------------------------------------
# Part I: PCAII

# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# NMDS input
# ----------------------------------------------------------------------

df <- raw_data


# ----------------------------------------------------------------------
# Sample columns
# ----------------------------------------------------------------------

sample_cols <- c(
  
  "04DC",
  "04JH",
  "04MH",
  "04RC",
  "04SX",
  "04TC",
  "04TJ",
  
  "10DC",
  "10JH",
  "10MH",
  "10RC",
  "10SX",
  "10TC",
  "10TJ"
)


# ----------------------------------------------------------------------
# Prepare NMDS data
# ----------------------------------------------------------------------

df_env <- df %>%
  
  filter(
    !is.na(指标) &
      指标 != ""
  ) %>%
  
  select(
    指标,
    all_of(sample_cols)
  ) %>%
  
  as.data.frame()


# ----------------------------------------------------------------------
# Environmental matrix
# 行 = samples
# 列 = environmental factors
# 和原代码逻辑完全一致
# ----------------------------------------------------------------------

rownames(df_env) <- df_env$指标


df_matrix <- df_env[
  ,
  -1,
  drop = FALSE
]


# 缺失值填充为0
df_matrix[
  is.na(df_matrix)
] <- 0


# 转置为：
# 样本 × 环境因子

df_transposed <- t(
  df_matrix
)


# ----------------------------------------------------------------------
# NMDS
# set.seed = 42
# distance = Bray
# k = 2
# trymax = 100
# ----------------------------------------------------------------------

set.seed(42)


nmds_result <- metaMDS(
  
  df_transposed,
  
  distance = "bray",
  
  k = 2,
  
  trymax = 100
)


# ----------------------------------------------------------------------
# NMDS scores and groups
# ----------------------------------------------------------------------

nmds_coords <- as.data.frame(
  
  scores(
    nmds_result,
    display = "sites"
  )
)


nmds_coords$Time <- ifelse(
  
  grepl(
    "^04",
    rownames(nmds_coords)
  ),
  
  "April 2025",
  
  "October 2025"
)


group_factor <- factor(
  nmds_coords$Time
)


# ----------------------------------------------------------------------
# PERMANOVA
# 保持你的原始代码：
# Bray-Curtis
# 999 unrestricted permutations
# ----------------------------------------------------------------------

adonis_result <- adonis2(
  
  df_transposed ~ group_factor,
  
  data = nmds_coords,
  
  method = "bray",
  
  permutations = 999
)


adonis_r2 <- round(
  
  adonis_result$R2[1],
  
  4
)


adonis_p <- round(
  
  adonis_result$`Pr(>F)`[1],
  
  4
)


# ----------------------------------------------------------------------
# Group hulls
# ----------------------------------------------------------------------

calculate_hull <- function(data){
  
  data[
    
    chull(
      data$NMDS1,
      data$NMDS2
    ),
    
  ]
}


hull_data <- nmds_coords %>%
  
  group_by(
    Time
  ) %>%
  
  do(
    calculate_hull(.)
  )


# ----------------------------------------------------------------------
# NMDS plot
# ----------------------------------------------------------------------

nmds_plot <- ggplot(
  
  nmds_coords,
  
  aes(
    x = NMDS1,
    y = NMDS2,
    color = Time,
    shape = Time
  )
  
) +
  
  
  geom_polygon(
    
    data = hull_data,
    
    aes(
      fill = Time
    ),
    
    alpha = 0.2,
    
    show.legend = FALSE
  ) +
  
  
  geom_point(
    size = 5,
    alpha = 0.8
  ) +
  
  
  geom_vline(
    
    xintercept = 0,
    
    linetype = "dashed",
    
    color = "gray50",
    
    linewidth = 0.8
  ) +
  
  
  geom_hline(
    
    yintercept = 0,
    
    linetype = "dashed",
    
    color = "gray50",
    
    linewidth = 0.8
  ) +
  
  
  
  scale_color_manual(
    
    values = c(
      "April 2025" = "#FFB3DD",
      "October 2025" = "#ACFFFF"
    )
  ) +
  
  
  scale_fill_manual(
    
    values = c(
      "April 2025" = "#FFB3DD",
      "October 2025" = "#ACFFFF"
    )
  ) +
  
  
  scale_shape_manual(
    
    values = c(
      "April 2025" = 16,
      "October 2025" = 17
    )
  ) +
  
  
  labs(
    
    title =
      "NMDS Analysis (Environmental Factor level)",
    
    x =
      "NMDS1",
    
    y =
      "NMDS2",
    
    color =
      "Time Point",
    
    shape =
      "Time Point"
  ) +
  
  
  annotate(
    
    "text",
    
    x =
      min(
        nmds_coords$NMDS1
      ) * 0.9,
    
    y =
      min(
        nmds_coords$NMDS2
      ) * 0.9,
    
    label = paste0(
      
      "Stress = ",
      
      round(
        nmds_result$stress,
        4
      ),
      
      "\n",
      
      "Adonis: R² = ",
      
      adonis_r2,
      
      ", p = ",
      
      adonis_p
    ),
    
    size = 4,
    
    fontface = "bold",
    
    hjust = 0,
    
    bbox = list(
      boxstyle = "round,pad=0.5",
      fill = "white",
      alpha = 0.9
    )
  ) +
  
  
  theme_bw() +
  
  
  theme(
    
    plot.title =
      element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      ),
    
    axis.title =
      element_text(
        size = 14,
        face = "bold"
      ),
    
    axis.text =
      element_text(
        size = 12
      ),
    
    legend.title =
      element_text(
        size = 12,
        face = "bold"
      ),
    
    legend.text =
      element_text(
        size = 10
      ),
    
    panel.grid =
      element_blank()
  )


# ----------------------------------------------------------------------
# Display NMDS
# ----------------------------------------------------------------------

print(
  nmds_plot
)


# ----------------------------------------------------------------------
# Save NMDS plot
# ----------------------------------------------------------------------

ggsave(
  
  "C:/Users/余山小可爱/Desktop/环境因子水平_NMDS图.png",
  
  plot = nmds_plot,
  
  width = 10,
  
  height = 8,
  
  dpi = 300,
  
  bg = "white"
)


ggsave(
  
  "C:/Users/余山小可爱/Desktop/环境因子水平_NMDS图.pdf",
  
  plot = nmds_plot,
  
  width = 10,
  
  height = 8,
  
  bg = "white"
)


# ----------------------------------------------------------------------
# NMDS/PERMANOVA results
# ----------------------------------------------------------------------

cat(
  "\n=== NMDS / PERMANOVA Analysis Results ===\n"
)


cat(
  "NMDS Stress：",
  round(
    nmds_result$stress,
    4
  ),
  "\n"
)


cat(
  "Adonis R²：",
  adonis_r2,
  "\n"
)


cat(
  "Adonis p-value：",
  adonis_p,
  "\n"
)


# ----------------------------------------------------------------------
