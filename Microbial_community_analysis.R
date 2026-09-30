# Microbial community analysis
#
# Modules:
#   1. Alpha diversity
#   2. Genome-level NMDS and PERMANOVA
#   3. iCAMP community assembly
#   4. iCAMP 3D pie plots
#   5. MAG Z-score heatmap
#
# Original data paths, sheet names, sample order, statistical parameters,
# plotting settings, and output filenames are retained.
# ==============================================================================

# ======================================================================
# 0. 安装并加载包
# ======================================================================

packages <- c(
  "picante",
  "readxl",
  "vegan",
  "ggplot2",
  "dplyr",
  "tidyr",
  "ape",
  "plotrix",
  "scales",
  "bigmemory"
)

for(pkg in packages){

  if(!requireNamespace(pkg, quietly = TRUE)){

    install.packages(pkg)
  }
}


library(picante)
library(readxl)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ape)
library(plotrix)
library(scales)
library(bigmemory)


# 检查 iCAMP
if(!requireNamespace("iCAMP", quietly = TRUE)){

  stop(
    "当前电脑没有安装 iCAMP 包，请先安装 iCAMP 后重新运行代码。"
  )
}

library(iCAMP)


# ======================================================================
# 1. 路径
# ======================================================================

data_path <- "C:/Users/余山小可爱/Desktop/物种代码原始数据.xlsx"

tree_path <- "C:/Users/余山小可爱/Desktop/gtdb_results.backbone.bac120.classify.tree"

desktop_path <- "C:/Users/余山小可爱/Desktop"


# ======================================================================
# 2. iCAMP工作目录
#
#  核心修复
#
# 在 iCAMP_WorkSpace 中为每次分析创建独立运行目录
# ======================================================================

icamp_base_dir <- file.path(
  desktop_path,
  "iCAMP_WorkSpace"
)


if(!dir.exists(icamp_base_dir)){

  dir.create(
    icamp_base_dir,
    recursive = TRUE
  )
}


# 创建本次运行目录
work_dir <- tempfile(
  pattern = paste0(
    "run_",
    format(
      Sys.time(),
      "%Y%m%d_%H%M%S"
    ),
    "_"
  ),
  tmpdir = icamp_base_dir
)


dir.create(
  work_dir,
  recursive = TRUE
)


cat(
  "\n=============================================\n"
)

cat(
  "本次 iCAMP 工作目录：\n"
)

cat(
  work_dir,
  "\n"
)

cat(
  "=============================================\n"
)


# ======================================================================
# 3. 样本顺序
# ======================================================================

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


# ======================================================================
# 4. 读取 Sheet1
# ======================================================================

sheet1_data <- read_excel(

  data_path,

  sheet = "Sheet1"
)


# ======================================================================
# 5. 基础数据检查
# ======================================================================

if(!"Genome" %in% names(sheet1_data)){

  stop(
    "Sheet1 中没有找到 Genome 列。"
  )
}


missing_samples <- setdiff(
  sample_cols,
  names(sheet1_data)
)


if(length(missing_samples) > 0){

  stop(
    paste0(
      "Sheet1 缺少样本：",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}


# ======================================================================
# ======================================================================
#
# PART I
#
# α DIVERSITY
#
# ======================================================================
# ======================================================================


# ======================================================================
# 6. 读取系统发育树
# ======================================================================

tree <- ape::read.tree(
  tree_path
)


# ======================================================================
# 7. α多样性函数
#
# 计算参数：
#
# Richness = presence / absence
# Shannon = presence / absence
# PD = presence / absence + phylogenetic tree
# ======================================================================

calculate_alpha <- function(
    data,
    sample_cols,
    tree
){

  # ------------------------------------------------------------
  # TPM矩阵
  # ------------------------------------------------------------

  abundance_mat <- data %>%

    select(
      all_of(sample_cols)
    ) %>%

    t() %>%

    as.matrix()


  storage.mode(
    abundance_mat
  ) <- "numeric"


  colnames(
    abundance_mat
  ) <- data$Genome


  # ------------------------------------------------------------
  # Presence / absence
  # ------------------------------------------------------------

  presence_mat <- ifelse(
    abundance_mat > 0,
    1,
    0
  )


  # ------------------------------------------------------------
  # 匹配系统发育树
  # ------------------------------------------------------------

  matched_data <- match.phylo.comm(
    tree,
    presence_mat
  )


  matched_tree <- matched_data$phy

  matched_mat <- matched_data$comm


  # ------------------------------------------------------------
  # Richness
  # ------------------------------------------------------------

  richness <- specnumber(
    matched_mat
  )


  # ------------------------------------------------------------
  # Shannon
  # ------------------------------------------------------------

  shannon <- diversity(
    matched_mat,
    "shannon"
  )


  # ------------------------------------------------------------
  # PD
  # ------------------------------------------------------------

  pd_result <- pd(
    matched_mat,
    matched_tree,
    include.root = FALSE
  )


  pd_value <- pd_result$PD


  # ------------------------------------------------------------
  # 输出
  # ------------------------------------------------------------

  alpha_res <- data.frame(

    Sample = rownames(
      matched_mat
    ),

    Time = ifelse(

      grepl(
        "^04",
        rownames(matched_mat)
      ),

      "202504",

      "202510"
    ),

    Richness = richness,

    Shannon = shannon,

    Phylogenetic = pd_value
  )


  return(
    alpha_res
  )
}


# ======================================================================
# 8. 计算 α多样性
# ======================================================================

alpha_result <- calculate_alpha(

  sheet1_data,

  sample_cols,

  tree
)


# ======================================================================
# 9. 长格式
# ======================================================================

alpha_long <- alpha_result %>%

  pivot_longer(

    cols = c(
      Richness,
      Shannon,
      Phylogenetic
    ),

    names_to = "Index",

    values_to = "Value"
  )


# ======================================================================
# 10. Wilcoxon
#
# 非配对 Wilcoxon 检验
# ======================================================================

p_values <- alpha_long %>%

  group_by(
    Index
  ) %>%

  summarise(

    p_value = wilcox.test(
      Value ~ Time
    )$p.value,

    .groups = "drop"
  ) %>%

  mutate(

    label = ifelse(

      p_value < 0.001,

      "p < 0.001",

      paste0(
        "p = ",
        signif(
          p_value,
          2
        )
      )
    )
  )


# ======================================================================
# 11. P值位置
# ======================================================================

y_pos_data <- alpha_long %>%

  group_by(
    Index
  ) %>%

  summarise(

    y_pos = max(
      Value,
      na.rm = TRUE
    ) * 1.15,

    .groups = "drop"
  )


p_values <- left_join(

  p_values,

  y_pos_data,

  by = "Index"
)


# ======================================================================
# 12. α多样性绘图
# ======================================================================

p_alpha <- ggplot(

  alpha_long,

  aes(
    x = Time,
    y = Value,
    fill = Time
  )

) +

  geom_boxplot(

    alpha = 0.7,

    width = 0.5
  ) +

  geom_jitter(

    width = 0.1,

    size = 1.8,

    alpha = 0.8,

    color = "black"
  ) +

  facet_wrap(

    ~Index,

    scales = "free_y",

    ncol = 3
  ) +

  scale_fill_manual(

    values = c(

      "202504" = "#FFB3DD",

      "202510" = "#ACFFFF"
    )
  ) +

  labs(

    title = "Comparison of α-diversity between 202504 and 202510",

    x = "",

    y = "α-diversity index value",

    fill = ""
  ) +

  geom_text(

    data = p_values,

    aes(

      x = 1.5,

      y = y_pos,

      label = label
    ),

    inherit.aes = FALSE,

    size = 4.2,

    fontface = "bold"
  ) +

  theme_minimal() +

  theme(

    plot.title = element_text(

      hjust = 0.5,

      size = 14,

      face = "bold"
    ),

    axis.text = element_text(
      size = 11
    ),

    axis.title.y = element_text(
      size = 12
    ),

    strip.text = element_text(

      size = 12,

      face = "bold"
    ),

    legend.position = "none",

    panel.grid.minor = element_blank()
  )


print(
  p_alpha
)


# ======================================================================
# 13. 保存 α多样性
# ======================================================================

ggsave(

  file.path(
    desktop_path,
    "α多样性差异图_修正版.png"
  ),

  plot = p_alpha,

  width = 12,

  height = 5,

  dpi = 300,

  bg = "white"
)


ggsave(

  file.path(
    desktop_path,
    "α多样性差异图_修正版.pdf"
  ),

  plot = p_alpha,

  width = 12,

  height = 5,

  bg = "white"
)


write.csv(

  alpha_result,

  file.path(
    desktop_path,
    "α多样性结果_修正版.csv"
  ),

  row.names = FALSE
)


write.csv(

  p_values,

  file.path(
    desktop_path,
    "α多样性_Wilcoxon_P值_修正版.csv"
  ),

  row.names = FALSE
)

# ======================================================================
# ======================================================================
#
# PART II
#
# GENOME NMDS + PERMANOVA
#
# ======================================================================
# ======================================================================


# ======================================================================
# 14. Genome矩阵
# ======================================================================

df_genome <- sheet1_data %>%

  filter(

    !is.na(Genome) &

      Genome != ""
  ) %>%

  select(

    Genome,

    all_of(sample_cols)
  ) %>%

  as.data.frame()


rownames(
  df_genome
) <- df_genome$Genome


df_matrix <- df_genome[
  ,
  -1,
  drop = FALSE
]


df_matrix[] <- lapply(

  df_matrix,

  function(x){

    as.numeric(
      as.character(x)
    )
  }
)


df_matrix[
  is.na(df_matrix)
] <- 0


df_transposed <- t(
  as.matrix(
    df_matrix
  )
)


# ======================================================================
# 15. NMDS
#
# 参数：
# Bray
# k = 2
# trymax = 100
# ======================================================================

set.seed(
  42
)


nmds_result <- metaMDS(

  df_transposed,

  distance = "bray",

  k = 2,

  trymax = 100
)


# ======================================================================
# 16. NMDS坐标
# ======================================================================

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


# ======================================================================
# 17. PERMANOVA
#
# 参数：
# Bray
# 999 unrestricted permutations
# ======================================================================

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


# ======================================================================
# 18. Convex hull
# ======================================================================

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


# ======================================================================
# 19. NMDS绘图
# ======================================================================

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

    title = "NMDS Analysis (Genome level)",

    x = "NMDS1",

    y = "NMDS2",

    color = "Time Point",

    shape = "Time Point"
  ) +

  annotate(

    "text",

    x = min(
      nmds_coords$NMDS1
    ) * 0.9,

    y = min(
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

    plot.title = element_text(

      hjust = 0.5,

      size = 16,

      face = "bold"
    ),

    axis.title = element_text(

      size = 14,

      face = "bold"
    ),

    axis.text = element_text(
      size = 12
    ),

    legend.title = element_text(

      size = 12,

      face = "bold"
    ),

    legend.text = element_text(
      size = 10
    ),

    panel.grid = element_blank()
  )


print(
  nmds_plot
)


# ======================================================================
# 20. 保存 NMDS
# ======================================================================

ggsave(

  file.path(
    desktop_path,
    "Genome水平_NMDS图.png"
  ),

  plot = nmds_plot,

  width = 10,

  height = 8,

  dpi = 300,

  bg = "white"
)


ggsave(

  file.path(
    desktop_path,
    "Genome水平_NMDS图.pdf"
  ),

  plot = nmds_plot,

  width = 10,

  height = 8,

  bg = "white"
)

# ======================================================================
# ======================================================================
#
# PART III
#
# iCAMP
#
# 分析流程：
# Sheet3环境因子 + dniche + phylogenetic binning + ps.bin
# + icamp.big + icamp.bins + icamp.boot
#
# 输出接口：final_table
# final_table 用于后续 3D Pie 绘图
# ======================================================================
# ======================================================================


# ======================================================================
# 21. iCAMP参数
# ======================================================================

rand_time <- 1000
nworker <- 4
ds_value <- 0.2
bin_size_limit <- 5
sig_index <- "SES.RC"


# ======================================================================
# 22. iCAMP community matrix
# ======================================================================

raw_data <- sheet1_data

comm_df <- raw_data[
  ,
  sample_cols
]

comm_df[] <- lapply(
  comm_df,
  function(x){
    as.numeric(
      as.character(x)
    )
  }
)

comm_df[is.na(comm_df)] <- 0

comm <- as.matrix(
  t(
    as.data.frame(
      comm_df
    )
  )
)

colnames(comm) <- as.character(raw_data$Genome)
storage.mode(comm) <- "numeric"

# 删除所有样本中总丰度为0的MAG
comm <- comm[
  ,
  colSums(comm, na.rm = TRUE) > 0,
  drop = FALSE
]

if(any(rowSums(comm) <= 0)){
  stop("存在总丰度为0的样本，请检查 Sheet1。")
}


# ======================================================================
# 23. 读取 Sheet3 环境因子
#
# Sheet3格式：
# 指标  04DC 04JH 04MH 04RC 04SX 04TC 04TJ
#       10DC 10JH 10MH 10RC 10SX 10TC 10TJ
# ======================================================================

env_raw <- read_excel(
  data_path,
  sheet = "Sheet3"
)

if(!"指标" %in% names(env_raw)){
  stop("Sheet3 中没有找到“指标”列。")
}

missing_env_samples <- setdiff(
  sample_cols,
  names(env_raw)
)

if(length(missing_env_samples) > 0){
  stop(
    paste0(
      "Sheet3 缺少样本：",
      paste(missing_env_samples, collapse = ", ")
    )
  )
}

env_values <- env_raw[
  ,
  sample_cols
]

env_values[] <- lapply(
  env_values,
  function(x){
    as.numeric(
      as.character(x)
    )
  }
)

env <- as.data.frame(
  t(
    as.matrix(env_values)
  )
)

colnames(env) <- as.character(env_raw$指标)
rownames(env) <- sample_cols

if(anyNA(env)){
  stop("Sheet3 环境数据存在 NA，请检查原始数据。")
}

# 删除零方差环境变量
env_sd <- apply(
  env,
  2,
  sd,
  na.rm = TRUE
)

zero_var_env <- names(
  env_sd[
    is.na(env_sd) | env_sd == 0
  ]
)

if(length(zero_var_env) > 0){
  cat("\n以下环境变量无变异，将删除：\n")
  print(zero_var_env)

  env <- env[
    ,
    !colnames(env) %in% zero_var_env,
    drop = FALSE
  ]
}

# Z-score标准化
env_scaled <- as.data.frame(
  scale(env)
)

rownames(env_scaled) <- rownames(env)


# ======================================================================
# 24. Season分组
# ======================================================================

treat <- data.frame(
  Season = c(
    rep("April", 7),
    rep("October", 7)
  ),
  row.names = sample_cols,
  stringsAsFactors = FALSE
)

# 统一样本顺序
comm <- comm[
  sample_cols,
  ,
  drop = FALSE
]

env_scaled <- env_scaled[
  sample_cols,
  ,
  drop = FALSE
]

treat <- treat[
  sample_cols,
  ,
  drop = FALSE
]


# ======================================================================
# 25. 系统发育树匹配
# ======================================================================

tree_icamp <- ape::read.tree(
  tree_path
)

common_ids <- intersect(
  colnames(comm),
  tree_icamp$tip.label
)

cat("\n=============================================\n")
cat("iCAMP MAG匹配结果\n")
cat("群落矩阵MAG数：", ncol(comm), "\n")
cat("系统发育树tip数：", length(tree_icamp$tip.label), "\n")
cat("共同MAG数：", length(common_ids), "\n")

if(length(common_ids) == 0){
  stop("群落矩阵 Genome ID 与系统发育树 tip 完全无法匹配。")
}

comm <- comm[
  ,
  common_ids,
  drop = FALSE
]

tree_icamp <- ape::keep.tip(
  tree_icamp,
  common_ids
)

# iCAMP bin汇总需要的分类对象
clas <- data.frame(
  TaxonID = colnames(comm),
  row.names = colnames(comm),
  stringsAsFactors = FALSE
)


# ======================================================================
# 26. 检查本次独立工作目录
# ======================================================================

old_pd_files <- c(
  "pd.bin",
  "pd.desc",
  "pd.taxon.name.csv"
)

existing_pd_files <- old_pd_files[
  file.exists(
    file.path(
      work_dir,
      old_pd_files
    )
  )
]

if(length(existing_pd_files) > 0){
  stop(
    paste0(
      "新的 iCAMP 工作目录中意外存在旧 pd 文件：",
      paste(existing_pd_files, collapse = ", ")
    )
  )
}


# ======================================================================
# 27. 生成系统发育距离大矩阵
# ======================================================================

cat("\n运行 pdist.big()...\n")

pd.big <- pdist.big(
  tree = tree_icamp,
  wd = work_dir,
  nworker = nworker,
  memory.G = 10
)

# ======================================================================
# 28. 环境生态位差异 dniche
# ======================================================================

cat("\n运行 dniche()...\n")

niche.dif <- dniche(
  env = env_scaled,
  comm = comm,
  method = "niche.value",
  nworker = nworker,
  out.dist = FALSE,
  bigmemo = TRUE,
  nd.wd = work_dir
)

# ======================================================================
# 29. 如有需要，对系统发育树定根
# ======================================================================

if(!ape::is.rooted(tree_icamp)){

  cat("\n系统发育树未定根，运行 midpoint.root.big()...\n")

  rt <- midpoint.root.big(
    tree = tree_icamp,
    pd.desc = pd.big$pd.file,
    pd.spname = pd.big$tip.label,
    pd.wd = pd.big$pd.wd,
    nworker = nworker
  )

  tree_icamp <- rt$tree

}


# ======================================================================
# 30. 系统发育分箱
# ======================================================================

cat("\n运行 taxa.binphy.big()...\n")

phylobin <- taxa.binphy.big(
  tree = tree_icamp,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  ds = ds_value,
  bin.size.limit = bin_size_limit,
  nworker = nworker
)

sp.bin <- phylobin$sp.bin[
  ,
  3,
  drop = FALSE
]

# ======================================================================
# 31. 系统发育信号检测 ps.bin
# ======================================================================

row_total <- rowSums(comm)

sp.ra <- colMeans(
  sweep(
    comm,
    1,
    row_total,
    "/"
  )
)

spname.use <- colnames(comm)

cat("\n运行 ps.bin()...\n")

binps <- ps.bin(
  sp.bin = sp.bin,
  sp.ra = sp.ra,
  spname.use = spname.use,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  nd.list = niche.dif$nd,
  nd.spname = niche.dif$names,
  ndbig.wd = niche.dif$nd.wd,
  cor.method = "pearson",
  r.cut = 0.1,
  p.cut = 0.05,
  min.spn = 5
)

write.csv(
  data.frame(
    ds = ds_value,
    n.min = bin_size_limit,
    binps$Index
  ),
  file.path(
    desktop_path,
    "iCAMP_Phylogenetic_Signal_Summary.csv"
  ),
  row.names = FALSE
)

if(!is.null(binps$detail)){

  detail.df <- cbind(
    BinID = rownames(binps$detail),
    binps$detail
  )

  write.csv(
    detail.df,
    file.path(
      desktop_path,
      "iCAMP_Phylogenetic_Signal_Detail.csv"
    ),
    row.names = FALSE
  )
}

# ======================================================================
# 32. iCAMP核心计算
#
# 参数：
# rand = 1000
# nworker = 4
# ds = 0.2
# bin.size.limit = 5
# sig.index = SES.RC
# ======================================================================

cat("\n运行 icamp.big()...\n")

icres <- icamp.big(
  comm = comm,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  rand = rand_time,
  tree = tree_icamp,
  prefix = "Laojun_MAG",
  ds = ds_value,
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "bMPD",
  sig.index = sig_index,
  bin.size.limit = bin_size_limit,
  nworker = nworker,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = FALSE,
  ignore.zero = TRUE,
  output.wd = work_dir,
  correct.special = TRUE,
  unit.sum = rowSums(comm),
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL
)

# ======================================================================
# 33. 检查 SES.RC 输出
# ======================================================================

if(is.null(icres$bNRIiRCa)){
  stop(
    paste0(
      "当前 sig.index = ",
      sig_index,
      "，但 icres 中没有生成 bNRIiRCa。"
    )
  )
}

write.csv(
  icres$bNRIiRCa,
  file.path(
    desktop_path,
    "iCAMP_Pairwise_Result.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 34. icamp.bins：正式计算组水平五种组装过程贡献
# ======================================================================

cat("\n运行 icamp.bins()...\n")

icbin <- icamp.bins(
  icamp.detail = icres$detail,
  treat = treat,
  clas = clas,
  silent = FALSE,
  boot = TRUE,
  rand.time = rand_time,
  between.group = TRUE
)

save(
  icbin,
  file = file.path(
    work_dir,
    "Laojun_MAG_iCAMP_Summary.rda"
  )
)

if(!is.null(icbin$Pt)){
  write.csv(
    icbin$Pt,
    file.path(
      desktop_path,
      "iCAMP_ProcessImportance_EachGroup.csv"
    ),
    row.names = FALSE
  )
}

if(!is.null(icbin$Ptk)){
  write.csv(
    icbin$Ptk,
    file.path(
      desktop_path,
      "iCAMP_ProcessImportance_EachBin_EachGroup.csv"
    ),
    row.names = FALSE
  )
}

if(!is.null(icbin$Ptuv)){
  write.csv(
    icbin$Ptuv,
    file.path(
      desktop_path,
      "iCAMP_ProcessImportance_EachTurnover.csv"
    ),
    row.names = FALSE
  )
}

if(!is.null(icbin$BPtk)){
  write.csv(
    icbin$BPtk,
    file.path(
      desktop_path,
      "iCAMP_BinContributeToProcess_EachGroup.csv"
    ),
    row.names = FALSE
  )
}

# ======================================================================
# 35. Bootstrap：April vs October
#
# SES.RC 对应：icres$bNRIiRCa
# ======================================================================

cat("\n运行 icamp.boot()...\n")

icboot <- icamp.boot(
  icamp.result = icres$bNRIiRCa,
  treat = treat,
  rand.time = rand_time,
  compare = TRUE,
  silent = FALSE,
  between.group = TRUE,
  ST.estimation = TRUE
)

save(
  icboot,
  file = file.path(
    work_dir,
    "Laojun_MAG_iCAMP_Bootstrap.rda"
  )
)

if(!is.null(icboot$summary)){
  write.csv(
    icboot$summary,
    file.path(
      desktop_path,
      "iCAMP_Bootstrap_Summary.csv"
    ),
    row.names = FALSE
  )
}

if(!is.null(icboot$compare)){
  write.csv(
    icboot$compare,
    file.path(
      desktop_path,
      "iCAMP_April_October_Comparison.csv"
    ),
    row.names = FALSE
  )
}

# ======================================================================
# 36. 生成第一版绘图继续使用的 final_table
#
# 使用 icamp.bins() 的组水平结果构建 final_table：
# Season / HeS / HoS / DL / HD / DR
# ======================================================================

if(is.null(icbin$Pt)){
  stop("icamp.bins() 没有生成 icbin$Pt，无法构建 final_table。")
}

if(!"Group" %in% names(icbin$Pt)){
  stop("icbin$Pt 中没有 Group 列，无法提取 April / October。")
}

required_process_cols <- c(
  "HeS",
  "HoS",
  "DL",
  "HD",
  "DR"
)

missing_process_cols <- setdiff(
  required_process_cols,
  names(icbin$Pt)
)

if(length(missing_process_cols) > 0){
  stop(
    paste0(
      "icbin$Pt 缺少过程列：",
      paste(missing_process_cols, collapse = ", ")
    )
  )
}

final_table <- icbin$Pt[
  icbin$Pt$Group %in% c("April", "October"),
  c(
    "Group",
    "HeS",
    "HoS",
    "DL",
    "HD",
    "DR"
  ),
  drop = FALSE
]

final_table <- final_table[
  match(
    c("April", "October"),
    final_table$Group
  ),
  ,
  drop = FALSE
]

colnames(final_table)[1] <- "Season"
rownames(final_table) <- NULL

if(nrow(final_table) != 2 || anyNA(final_table$Season)){
  stop("无法从 icbin$Pt 正确提取 April 和 October 两组结果。")
}

# ======================================================================
# 将五种组装过程转换为 numeric
# ======================================================================

for(col in required_process_cols){
  final_table[[col]] <- as.numeric(as.character(final_table[[col]]))
}

if(anyNA(final_table[, required_process_cols, drop = FALSE])){
  stop(
    paste0(
      "HeS / HoS / DL / HD / DR 转换为 numeric 后出现 NA。",
      "请运行 str(icbin$Pt) 检查原始结果。"
    )
  )
}

# 检查每个季节五种过程之和
process_sum <- rowSums(
  final_table[, required_process_cols, drop = FALSE]
)

cat("\n=== 五种组装过程之和 ===\n")
print(
  data.frame(
    Season = final_table$Season,
    Sum = process_sum
  )
)

# 保存比例值
write.csv(
  final_table,
  file.path(
    desktop_path,
    "iCAMP_Final_Result_Table.csv"
  ),
  row.names = FALSE
)

# 输出百分比版本
final_table_percent <- final_table

final_table_percent[, required_process_cols] <-
  final_table_percent[, required_process_cols, drop = FALSE] * 100

write.csv(
  final_table_percent,
  file.path(
    desktop_path,
    "iCAMP_Final_Result_Table_percent.csv"
  ),
  row.names = FALSE
)

# 保存参数
parameter_table <- data.frame(
  Parameter = c(
    "rand",
    "nworker",
    "ds",
    "bin.size.limit",
    "sig.index",
    "MAG_number",
    "Sample_number",
    "Environmental_variable_number"
  ),
  Value = c(
    rand_time,
    nworker,
    ds_value,
    bin_size_limit,
    sig_index,
    ncol(comm),
    nrow(comm),
    ncol(env_scaled)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  parameter_table,
  file.path(
    desktop_path,
    "iCAMP_Parameters.csv"
  ),
  row.names = FALSE
)

# 保存本次工作目录
writeLines(
  work_dir,
  con = file.path(
    desktop_path,
    "iCAMP_本次运行目录.txt"
  )
)

cat("\n=== iCAMP Final Result ===\n")
print(final_table)

cat("\n=== iCAMP Final Result (%) ===\n")
print(final_table_percent)

# ======================================================================
# ======================================================================
#
# PART IV
#
# iCAMP 3D PIE
#
# ======================================================================
# ======================================================================


# ======================================================================
# 34. 数据整理
# ======================================================================

season_summary <- final_table %>%

  mutate(

    across(

      c(
        HeS,
        HoS,
        DL,
        HD,
        DR
      ),

      ~as.numeric(
        as.character(.)
      )
    )
  ) %>%

  group_by(
    Season
  ) %>%

  summarise(

    across(

      c(
        HeS,
        HoS,
        DL,
        HD,
        DR
      ),

      sum,

      na.rm = TRUE
    ),

    .groups = "drop"
  )


# 按 April -> October 排序
season_summary <- season_summary[
  match(c("April", "October"), season_summary$Season),
  ,
  drop = FALSE
]

if(anyNA(season_summary$Season)){
  stop("season_summary 中没有同时找到 April 和 October。")
}


# ======================================================================
# 35. 图形参数
#
# 绘图参数
# ======================================================================

labels_vector <- c(
  "HeS",
  "HoS",
  "DL",
  "HD",
  "DR"
)


my_colors <- c(

  "#A9D1DF",

  "#E9E7C3",

  "#F3B2C3",

  "#8DE489",

  "#D3D4EE"
)


legend_labels <- c(

  "HeS (Heterogeneous Selection)",

  "HoS (Homogeneous Selection)",

  "DL (Dispersal Limitation)",

  "HD (Homogenizing Dispersal)",

  "Drift and Others"
)


# ======================================================================
# 36. April
# ======================================================================

season1_name <- as.character(
  season_summary$Season[1]
)


data_s1 <- as.numeric(
  season_summary[
    1,
    -1
  ]
)


pct_s1 <- round(

  data_s1 /
    sum(data_s1) *
    100,

  1
)


labels_s1 <- paste(

  labels_vector,

  "\n",

  pct_s1,

  "%",

  sep = ""
)


# ======================================================================
# 37. October
# ======================================================================

season2_name <- as.character(
  season_summary$Season[2]
)


data_s2 <- as.numeric(
  season_summary[
    2,
    -1
  ]
)


pct_s2 <- round(

  data_s2 /
    sum(data_s2) *
    100,

  1
)


labels_s2 <- paste(

  labels_vector,

  "\n",

  pct_s2,

  "%",

  sep = ""
)


# ======================================================================
# 38. 3D Pie函数
#
# 绘图参数
# ======================================================================

draw_3d_pies_with_legend <- function(){


  layout(

    matrix(
      c(
        1,
        2,
        3
      ),
      nrow = 1
    ),

    widths = c(
      4,
      4,
      3
    )
  )


  # April

  par(
    mar = c(
      3,
      2,
      5,
      2
    )
  )


  pie3D(

    data_s1,

    labels = labels_s1,

    explode = 0.08,

    col = my_colors,

    main = paste(
      "Season:",
      season1_name
    ),

    labelcex = 0.9,

    theta = pi / 3,

    height = 0.1
  )


  # October

  par(
    mar = c(
      3,
      2,
      5,
      2
    )
  )


  pie3D(

    data_s2,

    labels = labels_s2,

    explode = 0.08,

    col = my_colors,

    main = paste(
      "Season:",
      season2_name
    ),

    labelcex = 0.9,

    theta = pi / 3,

    height = 0.1
  )


  # Legend

  plot.new()


  par(
    mar = c(
      1,
      0,
      1,
      1
    )
  )


  legend(

    x = "center",

    legend = legend_labels,

    fill = my_colors,

    bty = "n",

    cex = 1.0,

    y.intersp = 2.0,

    title = "Processes"
  )


  layout(1)
}


# ======================================================================
# 39. 保存 iCAMP饼图
# ======================================================================

pdf_file <- file.path(

  desktop_path,

  "Season_3D_Pie_with_Legend.pdf"
)


png_file <- file.path(

  desktop_path,

  "Season_3D_Pie_with_Legend.png"
)


pdf(

  file = pdf_file,

  width = 12,

  height = 5.5
)


draw_3d_pies_with_legend()


dev.off()


png(

  file = png_file,

  width = 12,

  height = 5.5,

  units = "in",

  res = 300
)


draw_3d_pies_with_legend()


dev.off()

# ======================================================================
# ======================================================================
#
# PART V
#
# MAG HEATMAP
#
# Sheet2
#
# ======================================================================
# ======================================================================


# ======================================================================
# 40. 读取 Sheet2
# ======================================================================

heatmap_data <- read_excel(

  data_path,

  sheet = "Sheet2"
)


# ======================================================================
# 41. 检查 MAGs
# ======================================================================

if(!"MAGs" %in% names(heatmap_data)){

  stop(
    "Sheet2 中没有找到 MAGs 列。"
  )
}


# ======================================================================
# 42. 保存原始值
# ======================================================================

df_raw <- heatmap_data %>%

  select(

    MAGs,

    all_of(sample_cols)
  )


# ======================================================================
# 43. Z-score
#
# 参数：
# 按行标准化
# ======================================================================

df_mat <- heatmap_data %>%

  select(
    all_of(sample_cols)
  )


df_mat[] <- lapply(

  df_mat,

  function(x){

    as.numeric(
      as.character(x)
    )
  }
)


df_scaled <- t(

  scale(
    t(
      as.matrix(df_mat)
    )
  )
)


df_final <- cbind(

  MAGs = heatmap_data$MAGs,

  as.data.frame(
    df_scaled
  )
)


# ======================================================================
# 44. 长格式
# ======================================================================

df_long <- df_final %>%

  pivot_longer(

    cols = -MAGs,

    names_to = "Sample",

    values_to = "Value"
  ) %>%

  mutate(

    Sample = factor(

      Sample,

      levels = sample_cols
    ),

    MAGs = factor(

      MAGs,

      levels = unique(
        MAGs
      )
    )
  )


# ======================================================================
# 45. 判断原始值是否为0
# ======================================================================

df_long <- df_long %>%

  left_join(

    df_raw %>%

      pivot_longer(

        cols = -MAGs,

        names_to = "Sample",

        values_to = "Raw"
      ),

    by = c(
      "MAGs",
      "Sample"
    )
  ) %>%

  mutate(

    IsZero =
      Raw == 0
  )


# ======================================================================
# 46. Heatmap
#
# 参数保持原始：
#
# low = #99AADD
# mid = white
# high = #E6A6B2
# limits = -4 ~ 4
# zero = grey
# ======================================================================

p_heatmap <- ggplot(

  df_long,

  aes(
    x = Sample,
    y = MAGs
  )

) +

  geom_tile(

    aes(
      fill = Value
    ),

    color = NA
  ) +

  geom_tile(

    data = filter(
      df_long,
      IsZero
    ),

    fill = "#D3D3D3",

    color = NA
  ) +

  scale_fill_gradient2(

    low = "#99AADD",

    mid = "#FFFFFF",

    high = "#E6A6B2",

    midpoint = 0,

    limits = c(
      -4,
      4
    ),

    breaks = c(
      -4,
      -2,
      0,
      2,
      4
    )
  ) +

  labs(

    title = "Normalized Heatmap (Z-score)",

    x = "",

    y = ""
  ) +

  theme_bw() +

  theme(

    plot.title = element_text(

      hjust = 0.5,

      size = 16,

      face = "bold"
    ),

    axis.text.x = element_text(

      angle = 45,

      hjust = 1,

      size = 11
    ),

    axis.text.y = element_text(
      size = 10
    ),

    panel.grid = element_blank()
  ) +

  scale_y_discrete(

    limits = rev(
      unique(
        df_long$MAGs
      )
    )
  )


print(
  p_heatmap
)


# ======================================================================
# 47. 保存热图
# ======================================================================

ggsave(

  file.path(
    desktop_path,
    "MAGs_标准化热图_0值灰色.png"
  ),

  plot = p_heatmap,

  width = 16,

  height = 18,

  dpi = 300,

  bg = "white"
)


ggsave(

  file.path(
    desktop_path,
    "MAGs_标准化热图_0值灰色.pdf"
  ),

  plot = p_heatmap,

  width = 16,

  height = 18,

  bg = "white"
)
