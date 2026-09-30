# Environmental RDA, hierarchical partitioning, and piecewise SEM analyses
#
# Includes:
#   1. RDA and rdacca.hp analysis for microbial community data
#   2. RDA and rdacca.hp analysis for CAZyme data
#   3. Piecewise SEM for CAZyme functional groups
#   4. Seven-variable SEM including dissolved oxygen (DO)
#
# Original input/output paths, model structures, parameters, and file names are retained.
# ======================================================================

# RDA: environmental factors and microbial community
# ======================================================================
# RDA 排序分析与层次分割（rdacca.hp）
# ======================================================================
# 1. 环境准备与包加载
required_packages <- c("readxl", "tidyverse", "vegan", "ggplot2", "ggrepel", "rdacca.hp")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(readxl)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(rdacca.hp)

# 2. 数据读取与对齐预处理
file <- "C:/Users/余山小可爱/Desktop/SEM原始代码数据.xlsx"
output_path <- "C:/Users/余山小可爱/Desktop/"

env <- read_excel(file, sheet = "Sheet1")
otu <- read_excel(file, sheet = "Sheet2")

# 整理环境因子数据
env_df <- data.frame(t(env[,-1]))
colnames(env_df) <- env$指标
rownames(env_df) <- colnames(env)[-1]

# 整理群落（物种）数据 (排除前7列非丰度信息)
otu_mat <- as.matrix(otu[,-(1:7)])

# 样本完全对齐
common <- intersect(rownames(env_df), colnames(otu_mat))
env_use <- env_df[common,]
otu_use <- otu_mat[,common]

# 对环境因子进行标准化 (Z-score 标准化：均值为 0，标准差为 1)
# 这一步能有效消除不同环境因子量纲（单位）差异对 RDA 和层次分割计算造成的严重偏差
env_scaled <- as.data.frame(scale(env_use))

# 提取季节分组
group <- ifelse(grepl("^04", common), "202504", "202510")

# 对群落数据进行 Hellinger 转化
otu_hel <- decostand(t(otu_use), method = "hellinger")

# ----------------------------------------------------------------------
# 3. 构建全因子 RDA 模型与坐标提取（使用标准化后的环境数据）
# ----------------------------------------------------------------------
rda_final <- rda(otu_hel ~ ., data = env_scaled)

# 提取作图坐标
site  <- as.data.frame(scores(rda_final, display = "sites", scaling = 2))
arrow <- as.data.frame(scores(rda_final, display = "bp", scaling = 2))

colnames(site)  <- c("RDA1","RDA2")
colnames(arrow) <- c("RDA1","RDA2")
arrow$Var <- rownames(arrow)
site$Season <- group

# ----------------------------------------------------------------------
# 4. 层次分割分析（rdacca.hp）与显著性置换检验
# ----------------------------------------------------------------------
cat("\nRunning hierarchical partitioning...\n")
# 4.1 层次分割计算（计算独立贡献率 Individual Explanation）
rda_hp_res <- rdacca.hp(otu_hel, env_scaled, method = 'RDA', type = 'adjR2', scale = FALSE)

# 4.2 置换检验获取各个因子的显著性 P 值
set.seed(123)
permu_hp_res <- permu.hp(dv = otu_hel, iv = env_scaled, method = 'RDA', type = 'adjR2', permutations = 999)

# ----------------------------------------------------------------------
# 5. 绘图与输出部分
# ----------------------------------------------------------------------
### 5.1 绘制传统的 RDA 排序图
p_rda <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  
  # 湖泊样本点，按季节填色
  geom_point(data = site, aes(RDA1, RDA2, fill = Season), 
             size = 4, shape = 21, stroke = 1) +
  
  # 环境因子箭头
  geom_segment(data = arrow, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = grid::arrow(length = grid::unit(0.2, "cm")), 
               color = "black", linewidth = 0.7) +
  geom_text_repel(data = arrow, aes(RDA1, RDA2, label = Var), size = 4, fontface = "bold") +
  
  # 配色
  scale_fill_manual(values = c("202504" = "#FFB3DD", "202510" = "#ACFFFF"), name = "Season") +
  
  # 坐标轴标签 (动态提取解释率)
  labs(
    x = paste0("RDA1 (", round(summary(rda_final)$cont$importance[2,1]*100, 1), "%)"),
    y = paste0("RDA2 (", round(summary(rda_final)$cont$importance[2,2]*100, 1), "%)")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 保存 RDA 图
ggsave(paste0(output_path, "RDA_7环境因子_仅湖泊样本点.png"), p_rda, width = 10, height = 7, dpi = 300)
ggsave(paste0(output_path, "RDA_7环境因子_仅湖泊样本点.pdf"), p_rda, width = 10, height = 7)


### 5.2 整合数据并绘制层次分割玫瑰图（兼容不同 rdacca.hp 返回格式）
# 提取独立贡献率数据
hp_df <- as.data.frame(rda_hp_res$Hier.part)

# 智能兼容新旧版 rdacca.hp 列名命名（旧版多为第4列，新版字段名可能为 I.perc(%)）
if ("I.perc(%)" %in% colnames(hp_df)) {
  hp_clean <- data.frame(variable = rownames(hp_df), percentage = hp_df[["I.perc(%)"]])
} else {
  hp_clean <- data.frame(variable = rownames(hp_df), percentage = hp_df[, ncol(hp_df)]) # 锁定最后一列百分比
}

# 自动识别新旧版 permu.hp 的返回列表名（安全提取 P 值）
if (!is.null(permu_hp_res$Result)) {
  p_df <- as.data.frame(permu_hp_res$Result)
} else if (!is.null(permu_hp_res$Signif)) {
  p_df <- as.data.frame(permu_hp_res$Signif)
} else {
  p_df <- as.data.frame(permu_hp_res) 
}

# 采用最后一列位置索引提取 P 值
p_clean <- data.frame(
  variable = rownames(p_df),
  p_value = p_df[, ncol(p_df)]
)

# 合并清洗后的完整数据
data_rose <- merge(hp_clean, p_clean, by = "variable")

# 设置显著性分组与映射色彩
data_rose$significance <- ifelse(data_rose$p_value < 0.05, "Significant", "Not Significant")
data_rose$variable <- factor(data_rose$variable, levels = rownames(hp_df)) # 保持原输入顺序

# 绘制玫瑰图
p_rose <- ggplot(data_rose, aes(x = variable, y = percentage, fill = significance)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.2) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = -0.5, size = 3.5) +
  coord_polar(start = 1.5) +
  ylim(-5, max(data_rose$percentage) + 8) + # 给中心留白和顶部标签留出空间
  
  # 显著为淡紫红色，非显著为淡蓝色
  scale_fill_manual(values = c("Significant" = "#D19FE8", "Not Significant" = "#9AD5F1")) +
  
  labs(y = "Individual Contribution to Explained Variation (I.perc % by adjR2)", x = "", fill = "Significance (permu.hp)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "gray92"),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

# 保存玫瑰图
ggsave(paste0(output_path, "RDA_层次分割贡献率玫瑰图.png"), p_rose, width = 8, height = 7, dpi = 300)
ggsave(paste0(output_path, "RDA_层次分割贡献率玫瑰图.pdf"), p_rose, width = 8, height = 7)

# ----------------------------------------------------------------------
# 6. Console summary
# ----------------------------------------------------------------------
cat("\n=================== Statistical summary ===================\n")
cat("\n[1. 各环境因子经典 ANOVA 显著性检验结果]：\n")
print(anova(rda_final, by = "term"))

cat("\n[2. 全模型 VIF 共线性检查（已标准化环境数据，建议均小于 10）]：\n")
print(vif.cca(rda_final))

cat("\n[3. 层次分割（基于 adjR2）作图数据表核对]：\n")
print(data_rose)

cat("\n========================================================\n")
cat("RDA and hierarchical partitioning analysis completed.\n")
cat("========================================================\n")
# RDA: environmental factors and CAZymes
# ======================================================================
# RDA 排序分析与层次分割（rdacca.hp）
# 数据源：Sheet5 (CAZy_Family) | 环境因子已标准化 | 统一使用 adjR2
# ======================================================================
# ----------------------------------------------------------------------
# 1. 环境准备与包加载
# ----------------------------------------------------------------------
required_packages <- c("readxl", "tidyverse", "vegan", "ggplot2", "ggrepel", "rdacca.hp")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(readxl)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(rdacca.hp)

# ----------------------------------------------------------------------
# 2. 数据读取与对齐预处理（已切换至 Sheet5）
# ----------------------------------------------------------------------
file <- "C:/Users/余山小可爱/Desktop/SEM原始代码数据.xlsx"
output_path <- "C:/Users/余山小可爱/Desktop/"

env <- read_excel(file, sheet = "Sheet1")
otu <- read_excel(file, sheet = "Sheet3") # 读取包含 CAZy_Family 的 Sheet3

# 整理环境因子数据
env_df <- data.frame(t(env[,-1]))
colnames(env_df) <- env$指标
rownames(env_df) <- colnames(env)[-1]

# 整理群落（物种）数据 (Sheet5 第一列为 CAZy_Family 非丰度信息，通过 [,-1] 排除)
otu_mat <- as.matrix(otu[,-1]) 

# 样本完全对齐
common <- intersect(rownames(env_df), colnames(otu_mat))
env_use <- env_df[common,]
otu_use <- otu_mat[,common]

# 对环境因子进行标准化 (Z-score 标准化：均值为 0，标准差为 1)
# 消除不同环境因子（如 TN 与 pH）因量纲单位差异造成的分析偏倚
env_scaled <- as.data.frame(scale(env_use))

# 提取季节分组
group <- ifelse(grepl("^04", common), "202504", "202510")

# 对群落数据进行 Hellinger 转化
otu_hel <- decostand(t(otu_use), method = "hellinger")

# ----------------------------------------------------------------------
# 3. 构建全因子 RDA 模型与坐标提取（使用标准化后的环境数据）
# ----------------------------------------------------------------------
rda_final <- rda(otu_hel ~ ., data = env_scaled)

# 提取作图坐标
site  <- as.data.frame(scores(rda_final, display = "sites", scaling = 2))
arrow <- as.data.frame(scores(rda_final, display = "bp", scaling = 2))

colnames(site)  <- c("RDA1","RDA2")
colnames(arrow) <- c("RDA1","RDA2")
arrow$Var <- rownames(arrow)
site$Season <- group

# ----------------------------------------------------------------------
# 4. 层次分割分析（rdacca.hp）与显著性置换检验
# ----------------------------------------------------------------------
cat("\nRunning hierarchical partitioning...\n")
# 4.1 层次分割计算（计算独立贡献率 Individual Explanation）
rda_hp_res <- rdacca.hp(otu_hel, env_scaled, method = 'RDA', type = 'adjR2', scale = FALSE)

# 4.2 置换检验获取各个因子的显著性 P 值
set.seed(123)
permu_hp_res <- permu.hp(dv = otu_hel, iv = env_scaled, method = 'RDA', type = 'adjR2', permutations = 999)

# ----------------------------------------------------------------------
# 5. 绘图与输出部分
# ----------------------------------------------------------------------
### 5.1 绘制传统的 RDA 排序图
p_rda <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  
  # 湖泊样本点，按季节填色
  geom_point(data = site, aes(RDA1, RDA2, fill = Season), 
             size = 4, shape = 21, stroke = 1) +
  
  # 环境因子箭头
  geom_segment(data = arrow, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = grid::arrow(length = grid::unit(0.2, "cm")), 
               color = "black", linewidth = 0.7) +
  geom_text_repel(data = arrow, aes(RDA1, RDA2, label = Var), size = 4, fontface = "bold") +
  
  # 配色
  scale_fill_manual(values = c("202504" = "#FFB3DD", "202510" = "#ACFFFF"), name = "Season") +
  
  # 坐标轴标签 (动态提取解释率)
  labs(
    x = paste0("RDA1 (", round(summary(rda_final)$cont$importance[2,1]*100, 1), "%)"),
    y = paste0("RDA2 (", round(summary(rda_final)$cont$importance[2,2]*100, 1), "%)")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 保存 RDA 图
ggsave(paste0(output_path, "RDA_7环境因子_仅湖泊样本点.png"), p_rda, width = 10, height = 7, dpi = 300)
ggsave(paste0(output_path, "RDA_7环境因子_仅湖泊样本点.pdf"), p_rda, width = 10, height = 7)


### 5.2 整合数据并绘制层次分割玫瑰图（兼容不同 rdacca.hp 返回格式）
# 提取独立贡献率数据
hp_df <- as.data.frame(rda_hp_res$Hier.part)

# 兼容新旧版因 type = "adjR2" 带来的列名变化（如 "I.perc(%)" 字段）
if ("I.perc(%)" %in% colnames(hp_df)) {
  hp_clean <- data.frame(variable = rownames(hp_df), percentage = hp_df[["I.perc(%)"]])
} else {
  hp_clean <- data.frame(variable = rownames(hp_df), percentage = hp_df[, ncol(hp_df)]) # 锁定最后一列百分比
}

# 自动识别新旧版 permu.hp 的返回列表名
if (!is.null(permu_hp_res$Result)) {
  p_df <- as.data.frame(permu_hp_res$Result)
} else if (!is.null(permu_hp_res$Signif)) {
  p_df <- as.data.frame(permu_hp_res$Signif)
} else {
  p_df <- as.data.frame(permu_hp_res) 
}

# 采用位置索引安全提取 P 值（最后 columns 一列）
p_clean <- data.frame(
  variable = rownames(p_df),
  p_value = p_df[, ncol(p_df)]
)

# 合并清洗后的完整数据
data_rose <- merge(hp_clean, p_clean, by = "variable")

# 设置显著性分组与映射色彩
data_rose$significance <- ifelse(data_rose$p_value < 0.05, "Significant", "Not Significant")
data_rose$variable <- factor(data_rose$variable, levels = rownames(hp_df)) # 保持原输入顺序

# 绘制玫瑰图
p_rose <- ggplot(data_rose, aes(x = variable, y = percentage, fill = significance)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.2) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = -0.5, size = 3.5) +
  coord_polar(start = 1.5) +
  ylim(-5, max(data_rose$percentage) + 8) + # 给中心留白和顶部标签留出空间
  
  # === 修改后的新配色：显著浅粉紫，不显著浅蓝色 ===
  scale_fill_manual(values = c(
    "Significant"     = "#D6A2E8",  # 浅粉紫色
    "Not Significant" = "#AEDDEF"   # 浅蓝色
  )) +
  
  labs(y = "Individual Contribution to Explained Variation (I.perc % by adjR2)", x = "", fill = "Significance (permu.hp)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "gray92"),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

# 保存玫瑰图
ggsave(paste0(output_path, "RDA_层次分割贡献率玫瑰图.png"), p_rose, width = 8, height = 7, dpi = 300)
ggsave(paste0(output_path, "RDA_层次分割贡献率玫瑰图.pdf"), p_rose, width = 8, height = 7)

# ----------------------------------------------------------------------
# 6. Console summary
# ----------------------------------------------------------------------
cat("\n=================== Statistical summary ===================\n")
cat("\n[1. 各环境因子经典 ANOVA 显著性检验结果]：\n")
print(anova(rda_final, by = "term"))

cat("\n[2. 全模型 VIF 共线性检查（已标准化环境数据，建议均小于 10）]：\n")
print(vif.cca(rda_final))

cat("\n[3. 层次分割（基于 adjR2）作图数据表核对]：\n")
print(data_rose)

cat("\n========================================================\n")
cat("RDA and hierarchical partitioning analysis completed.\n")
cat("========================================================\n")
# Piecewise SEM: CAZyme functional groups
# ======================================================================
# 6 环境因子模型
# 环境-群落级联与平行功能酶响应网络 SEM 模型
# + Excel report output
# ======================================================================
# ======================================================================
# 0. 安装并加载必须的 R 包
# ======================================================================
required_packages <- c(
  "readxl",
  "vegan",
  "piecewiseSEM",
  "tidyverse",
  "writexl"
)

for(pkg in required_packages){
  
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
}

library(readxl)
library(vegan)
library(piecewiseSEM)
library(tidyverse)
library(writexl)


# ======================================================================
# 1. 路径设置
# ======================================================================
file_path <- "C:/Users/余山小可爱/Desktop/SEM原始代码数据.xlsx"

output_path <- "C:/Users/余山小可爱/Desktop/SEM_6因子模型完整报告.xlsx"


if(!is.null(dev.list())){
  dev.off()
}


# ======================================================================
# 2. 读取原始数据
#
# Sheet1 = 环境因子
# Sheet2 = 微生物群落
# Sheet3 = CAZyme
# ======================================================================
df_env_raw <- read_excel(
  file_path,
  sheet = "Sheet1"
)

df_species <- read_excel(
  file_path,
  sheet = "Sheet2"
)

df_cazy <- read_excel(
  file_path,
  sheet = "Sheet3"
)


# ======================================================================
# 3. 环境因子数据转置
# ======================================================================
env_names <- as.character(
  df_env_raw[[1]]
)


df_env <- as.data.frame(
  t(
    df_env_raw[, -1]
  )
)


colnames(df_env) <- env_names


df_env$SampleID <- rownames(
  df_env
)


# 强制环境变量为 numeric
for(col in env_names){
  
  df_env[[col]] <- as.numeric(
    as.character(
      df_env[[col]]
    )
  )
}


# ======================================================================
# 4. 通用数据转置函数
#
# 自动识别：
# 04...
# 10...
#
# 输出：
# 行 = sample
# 列 = feature
# ======================================================================
extract_and_transpose <- function(df){
  
  sample_cols <- grep(
    "^(04|10)",
    colnames(df),
    value = TRUE
  )
  
  
  if(length(sample_cols) == 0){
    
    stop(
      "没有识别到以04或10开头的样本列。"
    )
  }
  
  
  sub_df <- df[
    ,
    sample_cols,
    drop = FALSE
  ]
  
  
  # 强制转换成数值
  sub_df[] <- lapply(
    sub_df,
    function(x){
      as.numeric(
        as.character(x)
      )
    }
  )
  
  
  df_t <- as.data.frame(
    t(sub_df)
  )
  
  
  rownames(df_t) <- sample_cols
  
  
  return(
    df_t
  )
}


# ======================================================================
# 5. 微生物群落降维
#
# Bray-Curtis
# PCoA
# 提取第一轴
# ======================================================================
spe_mat_t <- extract_and_transpose(
  df_species
)


spe_dist <- vegan::vegdist(
  spe_mat_t,
  method = "bray"
)


pcoa_res <- cmdscale(
  spe_dist,
  k = 2,
  eig = TRUE
)


# 检查样本顺序
if(
  !identical(
    rownames(spe_mat_t),
    rownames(df_env)
  )
){
  
  stop(
    "环境因子和微生物群落样本顺序不一致，请检查数据。"
  )
}


df_env$Species_PCoA1 <- pcoa_res$points[, 1]


# ======================================================================
# 6. 六大 CAZyme 功能组独立降维
#
# AAs
# CBMs
# CEs
# GHs
# GTs
# PLs
# ======================================================================
families <- c(
  "AAs",
  "CBMs",
  "CEs",
  "GHs",
  "GTs",
  "PLs"
)


for(fam in families){
  
  # --------------------------------------------------------------------------
  # 去掉最后的 s
  #
  # AAs  -> AA
  # CBMs -> CBM
  # CEs  -> CE
  # GHs  -> GH
  # GTs  -> GT
  # PLs  -> PL
  # --------------------------------------------------------------------------
  
  match_prefix <- sub(
    "s$",
    "",
    fam
  )
  
  
  # --------------------------------------------------------------------------
  # 匹配第一列中对应的 CAZyme family
  # --------------------------------------------------------------------------
  
  fam_rows <- grep(
    paste0(
      "^",
      match_prefix
    ),
    as.character(
      df_cazy[[1]]
    ),
    ignore.case = TRUE
  )
  
  
  cat(
    "\n",
    fam,
    "匹配到",
    length(fam_rows),
    "个CAZyme家族\n"
  )
  
  
  # ==========================================================================
  # 情况1：
  # 匹配到两个及以上家族
  #
  # PCA → PC1
  # ==========================================================================
  
  if(length(fam_rows) > 1){
    
    sub_fam_df <- df_cazy[
      fam_rows,
      ,
      drop = FALSE
    ]
    
    
    fam_mat_t <- extract_and_transpose(
      sub_fam_df
    )
    
    
    # 检查样本顺序
    if(
      !identical(
        rownames(fam_mat_t),
        rownames(df_env)
      )
    ){
      
      stop(
        paste0(
          fam,
          " 的样本顺序与环境因子数据不一致。"
        )
      )
    }
    
    
    fam_pca <- prcomp(
      fam_mat_t,
      scale. = TRUE
    )
    
    
    # Store PC1 scores
    df_env[[paste0(fam, "_PC1")]] <- fam_pca$x[, 1]
  }
  
  
  # ==========================================================================
  # 情况2：
  # 只有一个家族
  #
  # 无法做 PCA
  # → 对这一变量进行 Z-score
  # ==========================================================================
  
  else if(length(fam_rows) == 1){
    
    sub_fam_df <- df_cazy[
      fam_rows,
      ,
      drop = FALSE
    ]
    
    
    fam_mat_t <- extract_and_transpose(
      sub_fam_df
    )
    
    
    if(
      !identical(
        rownames(fam_mat_t),
        rownames(df_env)
      )
    ){
      
      stop(
        paste0(
          fam,
          " 的样本顺序与环境因子数据不一致。"
        )
      )
    }
    
    
    # Store standardized score
    df_env[[paste0(fam, "_PC1")]] <- as.numeric(
      scale(
        fam_mat_t[, 1]
      )
    )
  }
  
  
  # ==========================================================================
  # 情况3：
  # 一个家族都没找到
  # ==========================================================================
  
  else {
    
    stop(
      paste0(
        "未在 Sheet3 中匹配到任何以 ",
        match_prefix,
        " 开头的 CAZyme family。"
      )
    )
  }
}


# ======================================================================
# 7. 查看最终用于 SEM 的变量
# ======================================================================
cat(
  "\n============================================================\n"
)

cat(
  "最终 SEM 数据结构\n"
)

cat(
  "============================================================\n"
)


print(
  names(df_env)
)


cat(
  "\n样本数 = ",
  nrow(df_env),
  "\n",
  sep = ""
)


# ======================================================================
# 8. 检查 SEM 必需变量
# ======================================================================
sem_variables <- c(
  "Temperature",
  "pH",
  "Transparency",
  "TN",
  "TP",
  "DOC",
  "Species_PCoA1",
  "GHs_PC1",
  "GTs_PC1",
  "PLs_PC1",
  "CEs_PC1",
  "AAs_PC1",
  "CBMs_PC1"
)


missing_sem_variables <- setdiff(
  sem_variables,
  names(df_env)
)


if(length(missing_sem_variables) > 0){
  
  stop(
    paste0(
      "SEM数据缺少以下变量：",
      paste(
        missing_sem_variables,
        collapse = ", "
      )
    )
  )
}


# ======================================================================
# 9. 检查 NA
# ======================================================================
na_summary <- colSums(
  is.na(
    df_env[, sem_variables]
  )
)


cat(
  "\n各SEM变量NA数量：\n"
)


print(
  na_summary
)


if(any(na_summary > 0)){
  
  warning(
    "部分SEM变量存在NA，请检查输出。"
  )
}


# ======================================================================
# 10. 构建 piecewise SEM
# ======================================================================
sem_6_env <- psem(
  
  
  # ============================================================================
  # 10.1 Environment cascade
  #
  # Temperature → TN
  # ============================================================================
  
  lm(
    TN ~ Temperature,
    data = df_env
  ),
  
  
  # ============================================================================
  # 10.2 Environmental variables → microbial community
  # ============================================================================
  
  lm(
    Species_PCoA1 ~
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC,
    data = df_env
  ),
  
  
  # ============================================================================
  # 10.3 Community + environment → CAZyme groups
  # ============================================================================
  
  lm(
    GHs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC,
    data = df_env
  ),
  
  
  lm(
    GTs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC,
    data = df_env
  ),
  
  
  lm(
    PLs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC,
    data = df_env
  ),
  
  
  lm(
    CEs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC,
    data = df_env
  ),
  
  
  lm(
    AAs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC,
    data = df_env
  ),
  
  
  lm(
    CBMs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC,
    data = df_env
  ),
  
  
  # ============================================================================
  # 10.4 Residual correlations among CAZyme groups
  # ============================================================================
  
  GHs_PC1 %~~% GTs_PC1,
  GHs_PC1 %~~% PLs_PC1,
  GHs_PC1 %~~% CEs_PC1,
  GHs_PC1 %~~% AAs_PC1,
  GHs_PC1 %~~% CBMs_PC1,
  
  GTs_PC1 %~~% PLs_PC1,
  GTs_PC1 %~~% CEs_PC1,
  GTs_PC1 %~~% AAs_PC1,
  GTs_PC1 %~~% CBMs_PC1,
  
  PLs_PC1 %~~% CEs_PC1,
  PLs_PC1 %~~% AAs_PC1,
  PLs_PC1 %~~% CBMs_PC1,
  
  CEs_PC1 %~~% AAs_PC1,
  CEs_PC1 %~~% CBMs_PC1,
  
  AAs_PC1 %~~% CBMs_PC1
)


# ======================================================================
# 11. 控制台输出完整 SEM summary
# ======================================================================
cat(
  "\n============================================================\n"
)

cat(
  "PIECEWISE SEM SUMMARY\n"
)

cat(
  "============================================================\n"
)


sem_summary <- summary(
  sem_6_env
)


print(
  sem_summary
)


# ======================================================================
# ======================================================================
#
# EXCEL REPORT
#
# ======================================================================
# ======================================================================
# ======================================================================
# 12. 保存完整 summary 文本
# ======================================================================
summary_text <- capture.output(
  summary(
    sem_6_env
  )
)


summary_text_df <- data.frame(
  SEM_Report = summary_text,
  stringsAsFactors = FALSE
)


# ======================================================================
# 13. 路径系数
#
# 包含：
# Estimate
# Std.Error
# P.Value
# Std.Estimate
# ======================================================================
coefficients_df <- tryCatch(
  
  as.data.frame(
    piecewiseSEM::coefs(
      sem_6_env,
      standardize = "scale"
    )
  ),
  
  error = function(e){
    
    data.frame(
      Error = paste0(
        "Coefficients extraction failed: ",
        e$message
      )
    )
  }
)


# ======================================================================
# 14. R²
# ======================================================================
r2_df <- tryCatch(
  
  as.data.frame(
    piecewiseSEM::rsquared(
      sem_6_env
    )
  ),
  
  error = function(e){
    
    data.frame(
      Error = paste0(
        "R-squared extraction failed: ",
        e$message
      )
    )
  }
)


# ======================================================================
# 15. Fisher's C
# ======================================================================
fisher_df <- tryCatch(
  
  as.data.frame(
    piecewiseSEM::fisherC(
      sem_6_env
    )
  ),
  
  error = function(e){
    
    data.frame(
      Error = paste0(
        "Fisher C extraction failed: ",
        e$message
      )
    )
  }
)


# ======================================================================
# 16. Directed separation
# ======================================================================
dsep_df <- tryCatch(
  
  as.data.frame(
    piecewiseSEM::dSep(
      sem_6_env
    )
  ),
  
  error = function(e){
    
    data.frame(
      Error = paste0(
        "Directed separation extraction failed: ",
        e$message
      )
    )
  }
)


# ======================================================================
# 17. SEM 子模型公式
# ======================================================================
model_formula_df <- data.frame(
  
  Model = c(
    "Environment cascade",
    "Microbial community",
    "GHs",
    "GTs",
    "PLs",
    "CEs",
    "AAs",
    "CBMs"
  ),
  
  Formula = c(
    
    "TN ~ Temperature",
    
    "Species_PCoA1 ~ Temperature + pH + Transparency + TN + TP + DOC",
    
    "GHs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + TN + TP + DOC",
    
    "GTs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + TN + TP + DOC",
    
    "PLs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + TN + TP + DOC",
    
    "CEs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + TN + TP + DOC",
    
    "AAs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + TN + TP + DOC",
    
    "CBMs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + TN + TP + DOC"
  ),
  
  stringsAsFactors = FALSE
)


# ======================================================================
# 18. 残差相关列表
# ======================================================================
correlated_errors_df <- data.frame(
  
  Variable1 = c(
    "GHs_PC1",
    "GHs_PC1",
    "GHs_PC1",
    "GHs_PC1",
    "GHs_PC1",
    "GTs_PC1",
    "GTs_PC1",
    "GTs_PC1",
    "GTs_PC1",
    "PLs_PC1",
    "PLs_PC1",
    "PLs_PC1",
    "CEs_PC1",
    "CEs_PC1",
    "AAs_PC1"
  ),
  
  Variable2 = c(
    "GTs_PC1",
    "PLs_PC1",
    "CEs_PC1",
    "AAs_PC1",
    "CBMs_PC1",
    "PLs_PC1",
    "CEs_PC1",
    "AAs_PC1",
    "CBMs_PC1",
    "CEs_PC1",
    "AAs_PC1",
    "CBMs_PC1",
    "AAs_PC1",
    "CBMs_PC1",
    "CBMs_PC1"
  ),
  
  stringsAsFactors = FALSE
)


# ======================================================================
# 19. 最终 SEM input data
# ======================================================================
sem_data_df <- df_env %>%
  
  select(
    SampleID,
    Temperature,
    pH,
    Transparency,
    TN,
    TP,
    DOC,
    Species_PCoA1,
    GHs_PC1,
    GTs_PC1,
    PLs_PC1,
    CEs_PC1,
    AAs_PC1,
    CBMs_PC1
  )


# ======================================================================
# 20. CAZyme PC1基本信息
#
# 方便确认每个PC1是否正常
# ======================================================================
pc1_summary_df <- data.frame(
  
  Variable = c(
    "Species_PCoA1",
    "GHs_PC1",
    "GTs_PC1",
    "PLs_PC1",
    "CEs_PC1",
    "AAs_PC1",
    "CBMs_PC1"
  ),
  
  Mean = sapply(
    df_env[
      ,
      c(
        "Species_PCoA1",
        "GHs_PC1",
        "GTs_PC1",
        "PLs_PC1",
        "CEs_PC1",
        "AAs_PC1",
        "CBMs_PC1"
      )
    ],
    mean,
    na.rm = TRUE
  ),
  
  SD = sapply(
    df_env[
      ,
      c(
        "Species_PCoA1",
        "GHs_PC1",
        "GTs_PC1",
        "PLs_PC1",
        "CEs_PC1",
        "AAs_PC1",
        "CBMs_PC1"
      )
    ],
    sd,
    na.rm = TRUE
  )
)


# ======================================================================
# 21. Excel输出
# ======================================================================
write_xlsx(
  
  list(
    
    Summary_text =
      summary_text_df,
    
    Coefficients =
      coefficients_df,
    
    R_squared =
      r2_df,
    
    Fisher_C =
      fisher_df,
    
    Directed_separation =
      dsep_df,
    
    Model_formulas =
      model_formula_df,
    
    Correlated_errors =
      correlated_errors_df,
    
    SEM_input_data =
      sem_data_df,
    
    PC1_summary =
      pc1_summary_df
  ),
  
  path = output_path
)


# ======================================================================
# 22. Completion message
# ======================================================================
cat(
  "\n============================================================\n"
)

cat(
  " SEM analysis completed\n"
)

cat(
  " Excel report saved to:\n"
)

cat(
  output_path,
  "\n"
)

cat(
  "\nExcel Sheets：\n"
)

cat(
  "1. Summary_text\n"
)

cat(
  "2. Coefficients\n"
)

cat(
  "3. R_squared\n"
)

cat(
  "4. Fisher_C\n"
)

cat(
  "5. Directed_separation\n"
)

cat(
  "6. Model_formulas\n"
)

cat(
  "7. Correlated_errors\n"
)

cat(
  "8. SEM_input_data\n"
)

cat(
  "9. PC1_summary\n"
)

cat(
  "============================================================\n"
)
# Piecewise SEM: seven environmental variables including DO
# ======================================================================
# 7 环境因子模型（包含 DO）
# 环境-群落级联与平行功能酶响应网络 SEM 模型
#
# 环境因子：
# Temperature
# pH
# Transparency
# TN
# TP
# DOC
# DO
#
# 环境级联：
# Temperature → TN
#
# + Excel report output
# ======================================================================
# ======================================================================
# 0. 安装并加载必须的 R 包
# ======================================================================
required_packages <- c(
  "readxl",
  "vegan",
  "piecewiseSEM",
  "tidyverse",
  "writexl"
)


for(pkg in required_packages){

  if(!requireNamespace(pkg, quietly = TRUE)){

    install.packages(pkg)

  }

}


library(readxl)
library(vegan)
library(piecewiseSEM)
library(tidyverse)
library(writexl)



# ======================================================================
# 1. 路径设置
# ======================================================================
file_path <-
  "C:/Users/余山小可爱/Desktop/SEM原始代码数据.xlsx"


output_path <-
  "C:/Users/余山小可爱/Desktop/SEM_7因子_DO模型完整报告.xlsx"



if(!is.null(dev.list())){

  dev.off()

}



# ======================================================================
# 2. 读取原始数据
#
# Sheet1 = 环境因子
# Sheet2 = 微生物群落
# Sheet3 = CAZyme
# ======================================================================
df_env_raw <- read_excel(
  file_path,
  sheet = "Sheet1"
)


df_species <- read_excel(
  file_path,
  sheet = "Sheet2"
)


df_cazy <- read_excel(
  file_path,
  sheet = "Sheet3"
)



# ======================================================================
# 3. 环境因子数据转置
# ======================================================================
env_names <- as.character(
  df_env_raw[[1]]
)


df_env <- as.data.frame(
  t(
    df_env_raw[, -1]
  )
)


colnames(df_env) <- env_names


df_env$SampleID <- rownames(df_env)



# ----------------------------------------------------------------------
# 强制所有环境变量转换为 numeric
# ----------------------------------------------------------------------
for(col in env_names){

  df_env[[col]] <- as.numeric(
    as.character(
      df_env[[col]]
    )
  )

}



# ======================================================================
# 4. 通用数据转置函数
#
# 自动识别：
# 04...
# 10...
#
# 输出：
# 行 = sample
# 列 = feature
# ======================================================================
extract_and_transpose <- function(df){

  sample_cols <- grep(
    "^(04|10)",
    colnames(df),
    value = TRUE
  )


  if(length(sample_cols) == 0){

    stop(
      "没有识别到以04或10开头的样本列。"
    )

  }


  sub_df <- df[
    ,
    sample_cols,
    drop = FALSE
  ]


  # --------------------------------------------------------------------------
  # 强制转换成数值
  # --------------------------------------------------------------------------

  sub_df[] <- lapply(
    sub_df,
    function(x){

      as.numeric(
        as.character(x)
      )

    }
  )


  df_t <- as.data.frame(
    t(sub_df)
  )


  rownames(df_t) <- sample_cols


  return(df_t)

}



# ======================================================================
# 5. 微生物群落降维
#
# Bray-Curtis
# PCoA
# 提取第一轴
# ======================================================================
spe_mat_t <- extract_and_transpose(
  df_species
)



spe_dist <- vegan::vegdist(
  spe_mat_t,
  method = "bray"
)



pcoa_res <- cmdscale(
  spe_dist,
  k = 2,
  eig = TRUE
)



# ----------------------------------------------------------------------
# 检查样本顺序
# ----------------------------------------------------------------------
if(
  !identical(
    rownames(spe_mat_t),
    rownames(df_env)
  )
){

  stop(
    "环境因子和微生物群落样本顺序不一致，请检查数据。"
  )

}



df_env$Species_PCoA1 <- pcoa_res$points[, 1]



# ======================================================================
# 6. 六大 CAZyme 功能组独立降维
#
# AAs
# CBMs
# CEs
# GHs
# GTs
# PLs
# ======================================================================
families <- c(
  "AAs",
  "CBMs",
  "CEs",
  "GHs",
  "GTs",
  "PLs"
)



for(fam in families){

  # --------------------------------------------------------------------------
  # 去掉最后的 s
  #
  # AAs  -> AA
  # CBMs -> CBM
  # CEs  -> CE
  # GHs  -> GH
  # GTs  -> GT
  # PLs  -> PL
  # --------------------------------------------------------------------------

  match_prefix <- sub(
    "s$",
    "",
    fam
  )


  # --------------------------------------------------------------------------
  # 匹配第一列中对应的 CAZyme family
  # --------------------------------------------------------------------------

  fam_rows <- grep(
    paste0(
      "^",
      match_prefix
    ),
    as.character(
      df_cazy[[1]]
    ),
    ignore.case = TRUE
  )


  cat(
    "\n",
    fam,
    "匹配到",
    length(fam_rows),
    "个CAZyme家族\n"
  )


  # ==========================================================================
  # 情况1：
  # 匹配到两个及以上家族
  #
  # PCA → PC1
  # ==========================================================================

  if(length(fam_rows) > 1){

    sub_fam_df <- df_cazy[
      fam_rows,
      ,
      drop = FALSE
    ]


    fam_mat_t <- extract_and_transpose(
      sub_fam_df
    )


    # ------------------------------------------------------------------------
    # 检查样本顺序
    # ------------------------------------------------------------------------

    if(
      !identical(
        rownames(fam_mat_t),
        rownames(df_env)
      )
    ){

      stop(
        paste0(
          fam,
          " 的样本顺序与环境因子数据不一致。"
        )
      )

    }


    fam_pca <- prcomp(
      fam_mat_t,
      scale. = TRUE
    )


    # ------------------------------------------------------------------------
    # Store PC1 scores
    #
    #
    #
    # ------------------------------------------------------------------------

    df_env[[paste0(fam, "_PC1")]] <- fam_pca$x[, 1]

  }


  # ==========================================================================
  # 情况2：
  # 只有一个家族
  #
  # 无法做 PCA
  # → 对这一变量进行 Z-score
  # ==========================================================================

  else if(length(fam_rows) == 1){

    sub_fam_df <- df_cazy[
      fam_rows,
      ,
      drop = FALSE
    ]


    fam_mat_t <- extract_and_transpose(
      sub_fam_df
    )


    if(
      !identical(
        rownames(fam_mat_t),
        rownames(df_env)
      )
    ){

      stop(
        paste0(
          fam,
          " 的样本顺序与环境因子数据不一致。"
        )
      )

    }


    # ------------------------------------------------------------------------
    # Store standardized score
    # ------------------------------------------------------------------------

    df_env[[paste0(fam, "_PC1")]] <- as.numeric(
      scale(
        fam_mat_t[, 1]
      )
    )

  }


  # ==========================================================================
  # 情况3：
  # 一个家族都没找到
  # ==========================================================================

  else{

    stop(
      paste0(
        "未在 Sheet3 中匹配到任何以 ",
        match_prefix,
        " 开头的 CAZyme family。"
      )
    )

  }

}



# ======================================================================
# 7. 查看最终用于 SEM 的变量
# ======================================================================
cat(
  "\n============================================================\n"
)


cat(
  "最终 SEM 数据结构\n"
)


cat(
  "============================================================\n"
)


print(
  names(df_env)
)


cat(
  "\n样本数 = ",
  nrow(df_env),
  "\n",
  sep = ""
)



# ======================================================================
# 8. 检查 SEM 必需变量
#
# 7个环境因子：
#
# Temperature
# pH
# Transparency
# TN
# TP
# DOC
# DO
#
# 无 Chla
# ======================================================================
sem_variables <- c(
  "Temperature",
  "pH",
  "Transparency",
  "TN",
  "TP",
  "DOC",
  "DO",
  "Species_PCoA1",
  "GHs_PC1",
  "GTs_PC1",
  "PLs_PC1",
  "CEs_PC1",
  "AAs_PC1",
  "CBMs_PC1"
)



missing_sem_variables <- setdiff(
  sem_variables,
  names(df_env)
)



if(length(missing_sem_variables) > 0){

  stop(
    paste0(
      "SEM数据缺少以下变量：",
      paste(
        missing_sem_variables,
        collapse = ", "
      )
    )
  )

}



# ======================================================================
# 9. 检查 NA
# ======================================================================
na_summary <- colSums(
  is.na(
    df_env[
      ,
      sem_variables
    ]
  )
)


cat(
  "\n各SEM变量NA数量：\n"
)


print(
  na_summary
)


if(any(na_summary > 0)){

  warning(
    "部分SEM变量存在NA，请检查输出。"
  )

}



# ======================================================================
# 10. 构建 piecewise SEM
#
# 7个环境因子：
#
# Temperature
# pH
# Transparency
# TN
# TP
# DOC
# DO
#
# 环境级联：
# Temperature → TN
#
# DO作为并列环境解释变量
#
# Chla已完全删除
# ======================================================================
sem_7_env <- psem(


  # ============================================================================
  # 10.1 Environment cascade
  #
  # Temperature → TN
  # ============================================================================

  lm(
    TN ~ Temperature,
    data = df_env
  ),


  # ============================================================================
  # 10.2 Environmental variables → microbial community
  # ============================================================================

  lm(
    Species_PCoA1 ~
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC +
      DO,
    data = df_env
  ),


  # ============================================================================
  # 10.3 Community + environment → CAZyme groups
  # ============================================================================

  lm(
    GHs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC +
      DO,
    data = df_env
  ),


  lm(
    GTs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC +
      DO,
    data = df_env
  ),


  lm(
    PLs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC +
      DO,
    data = df_env
  ),


  lm(
    CEs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC +
      DO,
    data = df_env
  ),


  lm(
    AAs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC +
      DO,
    data = df_env
  ),


  lm(
    CBMs_PC1 ~
      Species_PCoA1 +
      Temperature +
      pH +
      Transparency +
      TN +
      TP +
      DOC +
      DO,
    data = df_env
  ),


  # ============================================================================
  # 10.4 Residual correlations among CAZyme groups
  # ============================================================================

  GHs_PC1 %~~% GTs_PC1,
  GHs_PC1 %~~% PLs_PC1,
  GHs_PC1 %~~% CEs_PC1,
  GHs_PC1 %~~% AAs_PC1,
  GHs_PC1 %~~% CBMs_PC1,

  GTs_PC1 %~~% PLs_PC1,
  GTs_PC1 %~~% CEs_PC1,
  GTs_PC1 %~~% AAs_PC1,
  GTs_PC1 %~~% CBMs_PC1,

  PLs_PC1 %~~% CEs_PC1,
  PLs_PC1 %~~% AAs_PC1,
  PLs_PC1 %~~% CBMs_PC1,

  CEs_PC1 %~~% AAs_PC1,
  CEs_PC1 %~~% CBMs_PC1,

  AAs_PC1 %~~% CBMs_PC1

)



# ======================================================================
# 11. 控制台输出完整 SEM summary
# ======================================================================
cat(
  "\n============================================================\n"
)


cat(
  "PIECEWISE SEM SUMMARY - 7 ENVIRONMENTAL VARIABLES\n"
)


cat(
  "============================================================\n"
)


sem_summary <- summary(
  sem_7_env
)


print(
  sem_summary
)



# ======================================================================
# ======================================================================
#
# EXCEL REPORT
#
# ======================================================================
# ======================================================================
# ======================================================================
# 12. 保存完整 summary 文本
# ======================================================================
summary_text <- capture.output(
  summary(
    sem_7_env
  )
)


summary_text_df <- data.frame(
  SEM_Report = summary_text,
  stringsAsFactors = FALSE
)



# ======================================================================
# 13. 路径系数
# ======================================================================
coefficients_df <- tryCatch(

  as.data.frame(
    piecewiseSEM::coefs(
      sem_7_env,
      standardize = "scale"
    )
  ),

  error = function(e){

    data.frame(
      Error = paste0(
        "Coefficients extraction failed: ",
        e$message
      )
    )

  }

)



# ======================================================================
# 14. R²
# ======================================================================
r2_df <- tryCatch(

  as.data.frame(
    piecewiseSEM::rsquared(
      sem_7_env
    )
  ),

  error = function(e){

    data.frame(
      Error = paste0(
        "R-squared extraction failed: ",
        e$message
      )
    )

  }

)



# ======================================================================
# 15. Fisher's C
# ======================================================================
fisher_df <- tryCatch(

  as.data.frame(
    piecewiseSEM::fisherC(
      sem_7_env
    )
  ),

  error = function(e){

    data.frame(
      Error = paste0(
        "Fisher C extraction failed: ",
        e$message
      )
    )

  }

)



# ======================================================================
# 16. Directed separation
# ======================================================================
dsep_df <- tryCatch(

  as.data.frame(
    piecewiseSEM::dSep(
      sem_7_env
    )
  ),

  error = function(e){

    data.frame(
      Error = paste0(
        "Directed separation extraction failed: ",
        e$message
      )
    )

  }

)



# ======================================================================
# 17. SEM 子模型公式
# ======================================================================
model_formula_df <- data.frame(

  Model = c(
    "Environment cascade",
    "Microbial community",
    "GHs",
    "GTs",
    "PLs",
    "CEs",
    "AAs",
    "CBMs"
  ),

  Formula = c(

    "TN ~ Temperature",

    paste0(
      "Species_PCoA1 ~ Temperature + pH + Transparency + ",
      "TN + TP + DOC + DO"
    ),

    paste0(
      "GHs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + ",
      "TN + TP + DOC + DO"
    ),

    paste0(
      "GTs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + ",
      "TN + TP + DOC + DO"
    ),

    paste0(
      "PLs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + ",
      "TN + TP + DOC + DO"
    ),

    paste0(
      "CEs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + ",
      "TN + TP + DOC + DO"
    ),

    paste0(
      "AAs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + ",
      "TN + TP + DOC + DO"
    ),

    paste0(
      "CBMs_PC1 ~ Species_PCoA1 + Temperature + pH + Transparency + ",
      "TN + TP + DOC + DO"
    )

  ),

  stringsAsFactors = FALSE

)



# ======================================================================
# 18. 残差相关列表
# ======================================================================
correlated_errors_df <- data.frame(

  Variable1 = c(
    "GHs_PC1",
    "GHs_PC1",
    "GHs_PC1",
    "GHs_PC1",
    "GHs_PC1",

    "GTs_PC1",
    "GTs_PC1",
    "GTs_PC1",
    "GTs_PC1",

    "PLs_PC1",
    "PLs_PC1",
    "PLs_PC1",

    "CEs_PC1",
    "CEs_PC1",

    "AAs_PC1"
  ),

  Variable2 = c(
    "GTs_PC1",
    "PLs_PC1",
    "CEs_PC1",
    "AAs_PC1",
    "CBMs_PC1",

    "PLs_PC1",
    "CEs_PC1",
    "AAs_PC1",
    "CBMs_PC1",

    "CEs_PC1",
    "AAs_PC1",
    "CBMs_PC1",

    "AAs_PC1",
    "CBMs_PC1",

    "CBMs_PC1"
  ),

  stringsAsFactors = FALSE

)



# ======================================================================
# 19. 最终 SEM input data
# ======================================================================
sem_data_df <- df_env %>%

  select(
    SampleID,
    Temperature,
    pH,
    Transparency,
    TN,
    TP,
    DOC,
    DO,
    Species_PCoA1,
    GHs_PC1,
    GTs_PC1,
    PLs_PC1,
    CEs_PC1,
    AAs_PC1,
    CBMs_PC1
  )



# ======================================================================
# 20. CAZyme PC1 基本信息
# ======================================================================
pc1_variables <- c(
  "Species_PCoA1",
  "GHs_PC1",
  "GTs_PC1",
  "PLs_PC1",
  "CEs_PC1",
  "AAs_PC1",
  "CBMs_PC1"
)



pc1_summary_df <- data.frame(

  Variable = pc1_variables,

  Mean = sapply(
    df_env[, pc1_variables, drop = FALSE],
    mean,
    na.rm = TRUE
  ),

  SD = sapply(
    df_env[, pc1_variables, drop = FALSE],
    sd,
    na.rm = TRUE
  ),

  stringsAsFactors = FALSE

)



# ======================================================================
# 21. 环境因子基本统计信息
#
# Chla已删除
# ======================================================================
environment_variables <- c(
  "Temperature",
  "pH",
  "Transparency",
  "TN",
  "TP",
  "DOC",
  "DO"
)



environment_summary_df <- data.frame(

  Variable =
    environment_variables,

  Mean = sapply(
    df_env[
      ,
      environment_variables,
      drop = FALSE
    ],
    mean,
    na.rm = TRUE
  ),

  SD = sapply(
    df_env[
      ,
      environment_variables,
      drop = FALSE
    ],
    sd,
    na.rm = TRUE
  ),

  Min = sapply(
    df_env[
      ,
      environment_variables,
      drop = FALSE
    ],
    min,
    na.rm = TRUE
  ),

  Max = sapply(
    df_env[
      ,
      environment_variables,
      drop = FALSE
    ],
    max,
    na.rm = TRUE
  ),

  stringsAsFactors = FALSE

)



# ======================================================================
# 22. Excel 输出
# ======================================================================
write_xlsx(

  list(

    Summary_text =
      summary_text_df,

    Coefficients =
      coefficients_df,

    R_squared =
      r2_df,

    Fisher_C =
      fisher_df,

    Directed_separation =
      dsep_df,

    Model_formulas =
      model_formula_df,

    Correlated_errors =
      correlated_errors_df,

    SEM_input_data =
      sem_data_df,

    Environment_summary =
      environment_summary_df,

    PC1_summary =
      pc1_summary_df

  ),

  path =
    output_path

)



# ======================================================================
# 23. Completion message
# ======================================================================
cat(
  "\n============================================================\n"
)


cat(
  " SEM analysis completed\n"
)


cat(
  " Model: 7 environmental variables\n"
)


cat(
  " Temperature + pH + Transparency + TN + TP + DOC + DO\n"
)


cat(
  " Chla excluded\n"
)


cat(
  " DO included\n"
)


cat(
  " Environmental cascade: Temperature → TN\n"
)


cat(
  " Excel report saved to:\n"
)


cat(
  output_path,
  "\n"
)


cat(
  "\nExcel Sheets：\n"
)


cat(
  "1. Summary_text\n"
)


cat(
  "2. Coefficients\n"
)


cat(
  "3. R_squared\n"
)


cat(
  "4. Fisher_C\n"
)


cat(
  "5. Directed_separation\n"
)


cat(
  "6. Model_formulas\n"
)


cat(
  "7. Correlated_errors\n"
)


cat(
  "8. SEM_input_data\n"
)


cat(
  "9. Environment_summary\n"
)


cat(
  "10. PC1_summary\n"
)


cat(
  "============================================================\n"
)
