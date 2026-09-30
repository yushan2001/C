# CAZyme-MAG co-occurrence network analysis
# Spearman correlation networks, Louvain modules, Zi-Pi roles,
# and network topology statistics
# ======================================================================

# Packages
# ----------------------------------------------------------------------
required_packages <- c(
  "readxl", "dplyr", "igraph", "Hmisc", "tibble", "tidyr",
  "ggplot2", "openxlsx"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
file_path <- "C:/Users/余山小可爱/Desktop/网络原始代码数据.xlsx"
out_dir   <- "C:/Users/余山小可爱/Desktop"

# 相关性筛选阈值
r_cutoff <- 0.6
p_cutoff <- 0.05

# Multiple-testing adjustment
p_adjust_method <- "none"

# Louvain edge weighting
# TRUE: use correlation coefficient R as the edge weight
use_edge_weight_for_louvain <- TRUE

# Random seed for Louvain reproducibility
louvain_seed <- 123

# Zi-Pi thresholds
zi_cutoff <- 2.5
pi_cutoff <- 0.62

roles_levels <- c("Network hubs", "Module hubs", "Connectors", "Peripherals")


# Excel output cache
excel_sheets_list <- list()

# ----------------------------------------------------------------------
# Input data
# ----------------------------------------------------------------------
mags        <- read_excel(file_path, sheet = 1)
cazyme_data <- read_excel(file_path, sheet = 2)

mags$MAGs        <- as.character(mags$MAGs)
cazyme_data$MAGs <- as.character(cazyme_data$MAGs)

# Sample columns beginning with 04 or 10
sample_cols <- grep("^(04|10)", colnames(mags), value = TRUE)

if (length(sample_cols) == 0) {
  stop("没有识别到以 04 或 10 开头的样本列，请检查 mags 表格列名。")
}

# CAZyme family columns
# Exclude non-CAZyme taxonomy/identifier columns
cazyme_names <- grep("^(GH|GT|PL|CE|AA|CBM)", colnames(cazyme_data), value = TRUE)

if (length(cazyme_names) == 0) {
  stop("没有识别到 CAZyme 家族列，请检查 sheet 10 的列名是否以 GH/GT/PL/CE/AA/CBM 开头。")
}

# Convert abundance and copy-number columns to numeric
mags <- mags %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.))))

cazyme_data <- cazyme_data %>%
  mutate(across(all_of(cazyme_names), ~ suppressWarnings(as.numeric(.))))

# MAG identifiers
mag_names <- mags$MAGs

# ----------------------------------------------------------------------
# MAG taxonomy
# ----------------------------------------------------------------------
mag_tax_map <- mags %>%
  select(MAGs, Phylum, Order) %>%
  distinct() %>%
  rename(ID = MAGs)

# ----------------------------------------------------------------------
# MAG abundance in long format
# ----------------------------------------------------------------------
mag_long <- mags %>%
  select(MAGs, all_of(sample_cols)) %>%
  pivot_longer(
    cols = -MAGs,
    names_to = "Sample",
    values_to = "Abundance"
  )

# ----------------------------------------------------------------------
# CAZyme family abundance matrix
# CAZyme abundance = Σ MAG abundance × CAZyme copy number
# ----------------------------------------------------------------------
cazyme_mat_raw <- cazyme_data %>%
  select(MAGs, all_of(cazyme_names)) %>%
  pivot_longer(
    cols = -MAGs,
    names_to = "Gene",
    values_to = "Copy_Num"
  ) %>%
  inner_join(mag_long, by = "MAGs") %>%
  mutate(
    Copy_Num   = ifelse(is.na(Copy_Num), 0, Copy_Num),
    Abundance  = ifelse(is.na(Abundance), 0, Abundance),
    Gene_Abund = Abundance * Copy_Num
  ) %>%
  group_by(Sample, Gene) %>%
  summarise(
    Total_Abund = sum(Gene_Abund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Total_Abund,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample") %>%
  as.data.frame()

# ----------------------------------------------------------------------
# MAG abundance matrix
# ----------------------------------------------------------------------
mag_mat <- mags %>%
  select(MAGs, all_of(sample_cols)) %>%
  column_to_rownames("MAGs") %>%
  t() %>%
  as.data.frame()

# Align sample order
cazyme_mat <- as.data.frame(
  matrix(
    0,
    nrow = nrow(mag_mat),
    ncol = ncol(cazyme_mat_raw),
    dimnames = list(rownames(mag_mat), colnames(cazyme_mat_raw))
  )
)

common_samples <- intersect(rownames(mag_mat), rownames(cazyme_mat_raw))
cazyme_mat[common_samples, ] <- cazyme_mat_raw[common_samples, , drop = FALSE]

# ----------------------------------------------------------------------
# Combined MAG and CAZyme abundance matrix
# ----------------------------------------------------------------------
combined_cazyme <- cbind(mag_mat, cazyme_mat)

df_04_cazyme <- combined_cazyme[grepl("^04", rownames(combined_cazyme)), , drop = FALSE]
df_10_cazyme <- combined_cazyme[grepl("^10", rownames(combined_cazyme)), , drop = FALSE]

# ======================================================================
# MAG-CAZyme positive-correlation network
# ======================================================================
make_cazyme_network <- function(
    mat,
    time_point,
    tax_map,
    mag_names,
    cazyme_names,
    out_dir,
    r_cutoff = 0.6,
    p_cutoff = 0.05,
    p_adjust_method = "none"
) {
  
  cat(paste0("\n开始构建 [", time_point, "] MAG-CAZyme 正相关网络...\n"))
  
  # Remove zero-variance columns
  keep_var <- sapply(mat, function(x) {
    x <- as.numeric(x)
    v <- var(x, na.rm = TRUE)
    is.finite(v) && v > 0
  })
  
  mat_clean <- mat[, keep_var, drop = FALSE]
  
  all_g <- colnames(mat_clean)
  m_in  <- intersect(mag_names, all_g)
  l_in  <- intersect(cazyme_names, all_g)
  
  if (length(m_in) == 0 || length(l_in) == 0) {
    cat(">>> [", time_point, "] MAG 或 CAZyme 节点不足，跳过。\n")
    return(NULL)
  }
  
  # Spearman correlations
  corr <- rcorr(as.matrix(mat_clean), type = "spearman")
  r_mat <- corr$r
  p_mat <- corr$P
  
  # Retain MAG-CAZyme correlations only
  pair_grid <- expand.grid(
    Source = m_in,
    Target = l_in,
    stringsAsFactors = FALSE
  )
  
  pair_grid$R <- mapply(
    function(a, b) r_mat[a, b],
    pair_grid$Source,
    pair_grid$Target
  )
  
  pair_grid$P <- mapply(
    function(a, b) p_mat[a, b],
    pair_grid$Source,
    pair_grid$Target
  )
  
  pair_grid$P_adj <- p.adjust(pair_grid$P, method = p_adjust_method)
  
  edges_df <- pair_grid %>%
    filter(
      !is.na(R),
      !is.na(P_adj),
      R >= r_cutoff,
      P_adj < p_cutoff
    ) %>%
    select(Source, Target, R, P, P_adj)
  
  if (nrow(edges_df) == 0) {
    cat(">>> [", time_point, "] 无显著正相关边。\n")
    return(NULL)
  }
  
  # Node table
  node_ids <- sort(unique(c(edges_df$Source, edges_df$Target)))
  
  nodes_df <- tibble(ID = node_ids) %>%
    mutate(
      Group = ifelse(ID %in% mag_names, "Microbe", "Gene")
    ) %>%
    left_join(tax_map, by = "ID") %>%
    mutate(
      Type = ifelse(Group == "Gene", "CAZyme_Gene", Phylum),
      Type = ifelse(is.na(Type), "Unclassified", Type),
      Order = ifelse(Group == "Gene", ID, Order),
      Order = ifelse(is.na(Order), "Unclassified", Order)
    ) %>%
    select(ID, Type, Order, Group)
  
  # Export node and edge tables
  nodes_raw_path <- file.path(out_dir, paste0("CAZyme网络_节点_raw_", time_point, ".csv"))
  edges_path     <- file.path(out_dir, paste0("CAZyme网络_边_", time_point, ".csv"))
  
  write.csv(nodes_df, nodes_raw_path, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(edges_df, edges_path, row.names = FALSE, fileEncoding = "UTF-8")
  
  cat(">>> [", time_point, "] 网络导出完成！节点：", nrow(nodes_df),
      " 边：", nrow(edges_df), "\n")
  
  return(list(nodes = nodes_df, edges = edges_df))
}

# ======================================================================
# Zi-Pi roles and network topology
# ======================================================================
calculate_zi_pi_cazyme <- function(
    time_point,
    out_dir,
    zi_cutoff = 2.5,
    pi_cutoff = 0.62,
    use_edge_weight_for_louvain = TRUE,
    louvain_seed = 123
) {
  
  cat(paste0("\n===== CAZyme 共现网络 [", time_point, "] Zi-Pi 与拓扑指标计算 =====\n"))
  
  nodes_raw_path  <- file.path(out_dir, paste0("CAZyme网络_节点_raw_", time_point, ".csv"))
  edges_path      <- file.path(out_dir, paste0("CAZyme网络_边_", time_point, ".csv"))
  nodes_zipi_path <- file.path(out_dir, paste0("CAZyme网络_节点_ZiPi_", time_point, ".csv"))
  
  if (!file.exists(nodes_raw_path) || !file.exists(edges_path)) {
    cat(">>> [", time_point, "] 节点或边文件不存在，跳过。\n")
    return(NULL)
  }
  
  nodes_df <- read.csv(nodes_raw_path, stringsAsFactors = FALSE, check.names = FALSE)
  edges_df <- read.csv(edges_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  nodes_df <- nodes_df %>%
    drop_na(ID) %>%
    distinct(ID, .keep_all = TRUE) %>%
    select(ID, Type, Order, Group)
  
  edges_df <- edges_df %>%
    drop_na(Source, Target) %>%
    filter(Source %in% nodes_df$ID, Target %in% nodes_df$ID)
  
  if (nrow(edges_df) == 0) {
    cat(">>> [", time_point, "] 有效边为 0，跳过。\n")
    return(NULL)
  }
  
  # Build graph with R, P, and adjusted P values
  g <- graph_from_data_frame(
    d = edges_df[, c("Source", "Target", "R", "P", "P_adj")],
    directed = FALSE,
    vertices = nodes_df
  )
  
  # Correlation coefficient R as edge weight
  E(g)$weight <- as.numeric(E(g)$R)
  E(g)$weight[!is.finite(E(g)$weight)] <- 1
  
  # Louvain random seed
  set.seed(louvain_seed)
  
  if (use_edge_weight_for_louvain) {
    comm <- cluster_louvain(g, weights = E(g)$weight)
  } else {
    comm <- cluster_louvain(g, weights = NA)
  }
  
  V(g)$module <- membership(comm)
  
  all_nodes <- V(g)$name
  g_deg     <- degree(g)
  g_str     <- strength(g, weights = E(g)$weight)
  
  zi <- numeric(length(all_nodes))
  pi <- numeric(length(all_nodes))
  names(zi) <- all_nodes
  names(pi) <- all_nodes
  
# ----------------------------------------------------------------------
  # Zi and Pi
# ----------------------------------------------------------------------
  for (nn in all_nodes) {
    
    v_nn <- V(g)[name == nn]
    m_i  <- V(g)[v_nn]$module
    
    inn_names <- neighbors(g, v_nn)$name
    
    # Within-module degree for node nn
    k_in <- sum(V(g)[name %in% inn_names]$module == m_i)
    
    # Within-module degree distribution
    same_m <- V(g)[V(g)$module == m_i]$name
    
    inner_k <- sapply(same_m, function(x) {
      v_x <- V(g)[name == x]
      xn  <- neighbors(g, v_x)$name
      sum(V(g)[name %in% xn]$module == m_i)
    })
    
    mk <- mean(inner_k, na.rm = TRUE)
    sk <- sd(inner_k, na.rm = TRUE)
    
    zi[nn] <- ifelse(sk == 0 || is.na(sk), 0, (k_in - mk) / sk)
    
    # Participation coefficient across modules
    kt <- length(inn_names)
    
    if (kt == 0) {
      pi[nn] <- 0
    } else {
      neighbor_modules <- V(g)[name %in% inn_names]$module
      tbl <- table(neighbor_modules)
      pi[nn] <- 1 - sum((tbl / kt)^2)
    }
  }
  
  zipi <- tibble(
    ID        = all_nodes,
    Module_ID = as.integer(V(g)$module),
    Zi        = as.numeric(zi[all_nodes]),
    Pi        = as.numeric(pi[all_nodes]),
    Degree    = as.numeric(g_deg[all_nodes]),
    Strength  = as.numeric(g_str[all_nodes])
  ) %>%
    mutate(
      Roles = case_when(
        Zi >= zi_cutoff & Pi >= pi_cutoff ~ "Network hubs",
        Zi >= zi_cutoff & Pi <  pi_cutoff ~ "Module hubs",
        Zi <  zi_cutoff & Pi >= pi_cutoff ~ "Connectors",
        TRUE                              ~ "Peripherals"
      ),
      Roles = factor(Roles, levels = roles_levels)
    )
  
  res <- nodes_df %>%
    left_join(zipi, by = "ID") %>%
    mutate(
      Group = ifelse(Type == "CAZyme_Gene", "Gene", "Microbe")
    )
  
  # Export Zi-Pi node results
  # Plot original Pi and Zi coordinates without jitter
  write.csv(res, nodes_zipi_path, row.names = FALSE, fileEncoding = "UTF-8")
  
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
  p_base <- ggplot(res, aes(x = Pi, y = Zi)) +
    geom_hline(
      yintercept = zi_cutoff,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.6
    ) +
    geom_vline(
      xintercept = pi_cutoff,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.6
    ) +
    annotate("text", x = 0.31, y = 3.2, label = "Module hubs", size = 3.5, fontface = "bold") +
    annotate("text", x = 0.85, y = 3.2, label = "Network hubs", size = 3.5, fontface = "bold") +
    annotate("text", x = 0.31, y = 0.5, label = "Peripherals", size = 3.5, fontface = "bold") +
    annotate("text", x = 0.85, y = 0.5, label = "Connectors", size = 3.5, fontface = "bold") +
    geom_point(aes(color = Type, shape = Group), size = 3, alpha = 0.8) +
    scale_shape_manual(values = c("Microbe" = 16, "Gene" = 15)) +
    coord_cartesian(xlim = c(0, 0.9), clip = "off") +
    labs(
      title = paste0("CAZyme Co-occurrence Network Zi-Pi ", time_point),
      x = "Among-module connectivity (Pi)",
      y = "Within-module connectivity (Zi)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )
  
  ggsave(
    file.path(out_dir, paste0("CAZyme_ZiPi_", time_point, ".png")),
    p_base,
    width = 7.5,
    height = 6.5,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, paste0("CAZyme_ZiPi_", time_point, ".pdf")),
    p_base,
    width = 7.5,
    height = 6.5
  )
  
# ----------------------------------------------------------------------
  # Global network metrics
# ----------------------------------------------------------------------
  network_metrics <- data.frame(
    Property = c(
      "Total Nodes (节点总数)",
      "Total Edges (边总数)",
      "Network Density (网络密度)",
      "Average Degree (平均度)",
      "Average Strength (平均加权度)",
      "Average Path Length (平均路径长度)",
      "Network Diameter (网络直径)",
      "Average Clustering Coefficient (平均聚类系数)",
      "Modularity (模块化指数)",
      "Number of Modules (模块数量)",
      "Number of Components (连通分量数量)",
      "Largest Component Size (最大连通分量节点数)"
    ),
    Value = c(
      vcount(g),
      ecount(g),
      edge_density(g, loops = FALSE),
      mean(degree(g)),
      mean(strength(g, weights = E(g)$weight)),
      mean_distance(g, directed = FALSE, unconnected = TRUE, weights = NA),
      diameter(g, directed = FALSE, unconnected = TRUE, weights = NA),
      transitivity(g, type = "localaverage", isolates = "zero"),
      modularity(g, membership(comm), weights = if (use_edge_weight_for_louvain) E(g)$weight else NA),
      length(unique(V(g)$module)),
      components(g)$no,
      max(components(g)$csize)
    ),
    stringsAsFactors = FALSE
  )
  
# ----------------------------------------------------------------------
  # Node-role summary
# ----------------------------------------------------------------------
  role_count <- res %>%
    mutate(Roles = factor(as.character(Roles), levels = roles_levels)) %>%
    group_by(Roles, .drop = FALSE) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(Roles = as.character(Roles))
  
  role_percent <- role_count %>%
    mutate(
      Percentage_Value = round(Count / sum(Count) * 100, 4),
      Percentage = paste0(round(Count / sum(Count) * 100, 2), "%")
    )
  
  microbe_role <- res %>%
    filter(Group == "Microbe") %>%
    mutate(Roles = factor(as.character(Roles), levels = roles_levels)) %>%
    group_by(Type, Roles, .drop = FALSE) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = Roles,
      values_from = Count,
      values_fill = 0
    )
  
  gene_role <- res %>%
    filter(Group == "Gene") %>%
    mutate(Roles = factor(as.character(Roles), levels = roles_levels)) %>%
    group_by(ID, Roles, .drop = FALSE) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = Roles,
      values_from = Count,
      values_fill = 0
    )
  
# ----------------------------------------------------------------------
  # Zi-Pi and degree summaries
# ----------------------------------------------------------------------
  zipi_summary <- res %>%
    summarise(
      Mean_Zi = mean(Zi, na.rm = TRUE),
      SD_Zi   = sd(Zi, na.rm = TRUE),
      Max_Zi  = max(Zi, na.rm = TRUE),
      Min_Zi  = min(Zi, na.rm = TRUE),
      Mean_Pi = mean(Pi, na.rm = TRUE),
      SD_Pi   = sd(Pi, na.rm = TRUE),
      Max_Pi  = max(Pi, na.rm = TRUE),
      Min_Pi  = min(Pi, na.rm = TRUE)
    )
  
  degree_summary <- res %>%
    group_by(Group) %>%
    summarise(
      Mean_Degree   = mean(Degree, na.rm = TRUE),
      SD_Degree     = sd(Degree, na.rm = TRUE),
      Max_Degree    = max(Degree, na.rm = TRUE),
      Min_Degree    = min(Degree, na.rm = TRUE),
      Mean_Strength = mean(Strength, na.rm = TRUE),
      SD_Strength   = sd(Strength, na.rm = TRUE),
      Max_Strength  = max(Strength, na.rm = TRUE),
      Min_Strength  = min(Strength, na.rm = TRUE),
      .groups = "drop"
    )
  
# ----------------------------------------------------------------------
  # Module size
# ----------------------------------------------------------------------
  module_size <- res %>%
    group_by(Module_ID) %>%
    summarise(
      Total_Nodes   = n(),
      Microbe_Count = sum(Group == "Microbe"),
      Gene_Count    = sum(Group == "Gene"),
      .groups = "drop"
    ) %>%
    arrange(desc(Total_Nodes))
  
# ----------------------------------------------------------------------
  # Top nodes
# ----------------------------------------------------------------------
  top_degree <- res %>%
    arrange(desc(Degree)) %>%
    slice_head(n = 10) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
  top_strength <- res %>%
    arrange(desc(Strength)) %>%
    slice_head(n = 10) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
  top_zi <- res %>%
    arrange(desc(Zi)) %>%
    slice_head(n = 10) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
  top_pi <- res %>%
    arrange(desc(Pi)) %>%
    slice_head(n = 10) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
# ----------------------------------------------------------------------
  # Key-role nodes
# ----------------------------------------------------------------------
  module_hub_list <- res %>%
    filter(Roles == "Module hubs") %>%
    arrange(desc(Zi), desc(Degree)) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
  connector_list <- res %>%
    filter(Roles == "Connectors") %>%
    arrange(desc(Pi), desc(Degree)) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
  network_hub_list <- res %>%
    filter(Roles == "Network hubs") %>%
    arrange(desc(Zi), desc(Pi)) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
  key_role_nodes <- res %>%
    filter(Roles != "Peripherals") %>%
    arrange(Roles, desc(Zi), desc(Pi)) %>%
    select(ID, Type, Group, Module_ID, Degree, Strength, Zi, Pi, Roles)
  
# ----------------------------------------------------------------------
  # Coordinate overlap check
  # N > 1 indicates overlapping rounded Zi-Pi coordinates
# ----------------------------------------------------------------------
  coordinate_overlap <- res %>%
    mutate(
      Zi_round = round(Zi, 3),
      Pi_round = round(Pi, 3)
    ) %>%
    group_by(Pi_round, Zi_round, Roles) %>%
    summarise(
      N = n(),
      IDs = paste(ID, collapse = "; "),
      Groups = paste(Group, collapse = "; "),
      Types = paste(Type, collapse = "; "),
      .groups = "drop"
    ) %>%
    filter(N > 1) %>%
    arrange(desc(N), desc(Zi_round), desc(Pi_round))
  
# ----------------------------------------------------------------------
  # Excel output
  # Four sheets are retained for each season:
  # 1) Summary: network and node-role statistics
  # 2) RoleDetail: microbial phyla and CAZyme node roles
  # 3) KeyNodes: non-peripheral nodes
  # 4) TopNodes: top Degree / Strength / Zi / Pi
# ----------------------------------------------------------------------

  # Network, role, Zi-Pi, degree, and strength summaries
  summary_network <- network_metrics %>%
    transmute(
      Section = "Network metrics",
      Metric = Property,
      Group = NA_character_,
      Value = as.character(Value)
    )

  summary_roles <- role_percent %>%
    transmute(
      Section = "Role summary",
      Metric = Roles,
      Group = NA_character_,
      Value = paste0(Count, " (", Percentage, ")")
    )

  summary_zipi <- zipi_summary %>%
    pivot_longer(
      cols = everything(),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    transmute(
      Section = "Zi-Pi summary",
      Metric = Metric,
      Group = NA_character_,
      Value = as.character(Value)
    )

  summary_degree <- degree_summary %>%
    pivot_longer(
      cols = -Group,
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    transmute(
      Section = "Degree/Strength summary",
      Metric = Metric,
      Group = as.character(Group),
      Value = as.character(Value)
    )

  summary_all <- bind_rows(
    summary_network,
    summary_roles,
    summary_zipi,
    summary_degree
  )

  # Role details
  # Microbes summarized by phylum; CAZymes summarized by family ID
  role_detail_microbe <- microbe_role %>%
    mutate(
      Category = "Microbe_Phylum",
      Name = Type
    ) %>%
    select(Category, Name, all_of(roles_levels))

  role_detail_gene <- gene_role %>%
    mutate(
      Category = "CAZyme_Gene",
      Name = ID
    ) %>%
    select(Category, Name, all_of(roles_levels))

  role_detail <- bind_rows(
    role_detail_microbe,
    role_detail_gene
  )

  # Key nodes: Module hubs, Connectors, and Network hubs
  key_nodes_sheet <- key_role_nodes

  # Combine four Top-10 node rankings
  top_nodes_sheet <- bind_rows(
    top_degree %>% mutate(Ranking_Type = "Top Degree", .before = 1),
    top_strength %>% mutate(Ranking_Type = "Top Strength", .before = 1),
    top_zi %>% mutate(Ranking_Type = "Top Zi", .before = 1),
    top_pi %>% mutate(Ranking_Type = "Top Pi", .before = 1)
  )

  # Store tables for Excel export
  excel_sheets_list[[paste0(time_point, "_Summary")]]    <<- summary_all
  excel_sheets_list[[paste0(time_point, "_RoleDetail")]] <<- role_detail
  excel_sheets_list[[paste0(time_point, "_KeyNodes")]]   <<- key_nodes_sheet
  excel_sheets_list[[paste0(time_point, "_TopNodes")]]   <<- top_nodes_sheet

  cat(">>> [", time_point, "] Zi-Pi 和拓扑指标计算完成！\n")
  cat("    Module hubs: ", nrow(module_hub_list), "\n")
  cat("    Connectors: ", nrow(connector_list), "\n")
  cat("    Network hubs: ", nrow(network_hub_list), "\n")
  
  return(res)
}

# ======================================================================
# Build networks for both seasons
# ======================================================================
make_cazyme_network(
  mat = df_04_cazyme,
  time_point = "202504",
  tax_map = mag_tax_map,
  mag_names = mag_names,
  cazyme_names = cazyme_names,
  out_dir = out_dir,
  r_cutoff = r_cutoff,
  p_cutoff = p_cutoff,
  p_adjust_method = p_adjust_method
)

make_cazyme_network(
  mat = df_10_cazyme,
  time_point = "202510",
  tax_map = mag_tax_map,
  mag_names = mag_names,
  cazyme_names = cazyme_names,
  out_dir = out_dir,
  r_cutoff = r_cutoff,
  p_cutoff = p_cutoff,
  p_adjust_method = p_adjust_method
)

# ======================================================================
# Run Zi-Pi analysis for both seasons
# ======================================================================
res_202504 <- calculate_zi_pi_cazyme(
  time_point = "202504",
  out_dir = out_dir,
  zi_cutoff = zi_cutoff,
  pi_cutoff = pi_cutoff,
  use_edge_weight_for_louvain = use_edge_weight_for_louvain,
  louvain_seed = louvain_seed
)

res_202510 <- calculate_zi_pi_cazyme(
  time_point = "202510",
  out_dir = out_dir,
  zi_cutoff = zi_cutoff,
  pi_cutoff = pi_cutoff,
  use_edge_weight_for_louvain = use_edge_weight_for_louvain,
  louvain_seed = louvain_seed
)

# ======================================================================
# ======================================================================
output_excel_path <- file.path(out_dir, "CAZyme_ZiPi_Statistics_Simplified.xlsx")

cat("\n正在写入 Excel 文件...\n")

wb <- createWorkbook()

for (sheet_name in names(excel_sheets_list)) {
  
  # Excel worksheet names are limited to 31 characters
  safe_sheet_name <- substr(sheet_name, 1, 31)
  
  addWorksheet(wb, sheetName = safe_sheet_name)
  
  writeData(
    wb,
    sheet = safe_sheet_name,
    x = excel_sheets_list[[sheet_name]],
    withFilter = TRUE
  )
  
  if (ncol(excel_sheets_list[[sheet_name]]) > 0) {
    setColWidths(
      wb,
      sheet = safe_sheet_name,
      cols = 1:ncol(excel_sheets_list[[sheet_name]]),
      widths = "auto"
    )
  }
}

saveWorkbook(wb, output_excel_path, overwrite = TRUE)
