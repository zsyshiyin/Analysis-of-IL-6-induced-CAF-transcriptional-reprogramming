# ============================================================================
# 转录重编程基因集核心转录因子分析（规范化完整版）- 4组热图版（排列顺序修正）
# 主图：DoRothEA ≥ 4 | 补充图：DoRothEA ≥ 1 或文献已知
# 条形图：文献已知TF名称加粗 | 热图：显示全部4组样品（按GROUP_ORDER排列）
# 新增：额外输出主图和补充图TF的活性矩阵 (Table_S5, Table_S6)
# ============================================================================

# 随机种子
set.seed(42)
rm(list = ls())
set.seed(42)

# ============================================================================
# 工作路径与输出目录
# ============================================================================
work_dir <- "D:/zsy/SX/Fomal-final/32-CIPS-Tumor-TF"
setwd(work_dir)
output_dir <- "TF_analysis_results_DoRothEA_standardized"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("========================================\n")
cat("转录重编程基因集核心转录因子分析\n")
cat("========================================\n\n")

# ============================================================================
# ★★★ 全部可视化参数 ★★★
# ============================================================================

FIG_DPI <- 300

# ---- 阈值优化图 ----
fig0_width <- 9; fig0_height <- 6
fig0_title <- "Correlation Threshold Optimization"
fig0_title_size <- 14; fig0_title_face <- "plain"
fig0_axis_title_size <- 12; fig0_axis_text_size <- 10
fig0_avg_line_color <- "#2C7BB6"; fig0_median_line_color <- "#E4945A"
fig0_target_line_color <- "red"; fig0_selected_line_color <- "#C60036"
fig0_line_size <- 1.2; fig0_point_size <- 2

# ---- 主图条形图 (Figure 1: DoRothEA ≥ 4) ----
fig1_width <- 12; fig1_height <- 6
fig1_title <- "Core Transcription Factors (DoRothEA ≥ 4)"
fig1_title_size <- 20; fig1_title_face <- "plain"
fig1_subtitle_size <- 10; fig1_subtitle_color <- "#7F8C8D"
fig1_axis_title_size <- 24; fig1_axis_text_size <- 24
fig1_bar_width <- 0.8; fig1_bar_alpha <- 0.9
fig1_color_validated <- "#C60036"; fig1_color_inferred <- "#FDE0DD"
fig1_show_values <- TRUE; fig1_value_size <- 8
fig1_value_color <- "#2C3E50"

# ---- 主图表达热图 (Figure 2: 全部4组, 按GROUP_ORDER排列) ----
fig2_width <- 5; fig2_height <- 2
fig2_title <- "Core Transcription Factors - Expression"
fig2_color_low <- "#2C7BB6"; fig2_color_mid <- "white"; fig2_color_high <- "#D7191C"
fig2_fontsize_row <- 12; fig2_fontsize_col <- 12
fig2_cluster_rows <- FALSE; fig2_cluster_cols <- FALSE
fig2_show_rownames <- TRUE; fig2_show_colnames <- TRUE
fig2_scale_range <- 2
fig2_annotation_colors <- list(
  Group = c("DMEM" = "#2AA7DE", "Con" = "#25377F", "I6" = "#C60036", "R2" = "#E4945A")
)

# ---- 主图活性热图 (Figure 3: 全部4组, 按GROUP_ORDER排列) ----
fig3_width <- 5; fig3_height <- 2
fig3_title <- "Core Transcription Factors - Activity"
fig3_color_low <- "#2C7BB6"; fig3_color_mid <- "white"; fig3_color_high <- "#D7191C"
fig3_fontsize_row <- 12; fig3_fontsize_col <- 12
fig3_cluster_rows <- FALSE; fig3_cluster_cols <- FALSE
fig3_show_rownames <- TRUE; fig3_show_colnames <- TRUE
fig3_scale_range <- 2

# ---- 补充条形图 (Figure S1: DoRothEA ≥ 1 或文献已知) ----
figS1_width <- 9; figS1_height <- 6
figS1_title <- "Core Transcription Factors (All Candidates)"
figS1_title_size <- 16; figS1_title_face <- "plain"
figS1_subtitle_size <- 10; figS1_subtitle_color <- "#7F8C8D"
figS1_axis_title_size <- 18; figS1_axis_text_size <- 16
figS1_bar_width <- 0.8; figS1_bar_alpha <- 0.9
figS1_color_validated <- "#C60036"; figS1_color_inferred <- "#FDE0DD"
figS1_show_values <- TRUE; figS1_value_size <- 5
figS1_value_color <- "#2C3E50"

# ---- 补充表达热图 (Figure S2, 全部4组, 按GROUP_ORDER排列) ----
figS2_width <- 8; figS2_height <- 6
figS2_title <- "Core Transcription Factors - Expression (All Candidates)"
figS2_color_low <- "#2C7BB6"; figS2_color_mid <- "white"; figS2_color_high <- "#D7191C"
figS2_fontsize_row <- 11; figS2_fontsize_col <- 10
figS2_cluster_rows <- FALSE; figS2_cluster_cols <- FALSE
figS2_show_rownames <- TRUE; figS2_show_colnames <- TRUE
figS2_scale_range <- 2

# ---- 补充活性热图 (Figure S3, 全部4组, 按GROUP_ORDER排列) ----
figS3_width <- 6; figS3_height <- 2
figS3_title <- "Core Transcription Factors - Activity (All Candidates)"
figS3_color_low <- "#2C7BB6"; figS3_color_mid <- "white"; figS3_color_high <- "#D7191C"
figS3_fontsize_row <- 12; figS3_fontsize_col <- 12
figS3_cluster_rows <- FALSE; figS3_cluster_cols <- FALSE
figS3_show_rownames <- TRUE; figS3_show_colnames <- TRUE
figS3_scale_range <- 2

# ============================================================================
# 筛选参数
# ============================================================================
ENABLE_ACTIVITY_FILTER <- TRUE
MIN_TARGET_GENES <- 5
ACTIVITY_FC_THRESHOLD <- 1.1
AUTO_SELECT_THRESHOLD <- TRUE
DESIRED_TARGETS_PER_TF <- 12
THRESHOLD_MIN <- 0.80; THRESHOLD_MAX <- 0.95; THRESHOLD_STEP <- 0.01
MIN_DOROTHEA_FOR_MAIN <- 4

# ============================================================================
# 数据文件
# ============================================================================
EXPRESSION_FILE <- "expression_matrix.csv"
GENE_SET_FILE <- "gene_sets.csv"

GROUP_MAPPING <- list("DMEM" = "DMEM", "Con" = "Con", "I6" = "I6", "R2" = "R2")
GROUP_ORDER <- c("DMEM", "Con", "I6", "R2")
GROUP_COLORS <- c("DMEM" = "#2AA7DE", "Con" = "#25377F", "I6" = "#C60036", "R2" = "#E4945A")
EXPERIMENTAL_GROUPS <- c("I6")
CONTROL_GROUP <- "Con"

MAX_TFS_TO_ANALYZE <- 200
CORRELATION_PVALUE_THRESHOLD <- 0.05

# ============================================================================
# CAF已知转录因子
# ============================================================================
CAF_KNOWN_TFS <- c(
  "TP53", "NFKB1", "RELA", "STAT3", "JUN", "FOS",
  "SNAI1", "TWIST1", "HIF1A", "GATA6", "SMAD4", "E2F3",
  "RBPJ", "NFE2L2", "FOXO1", "FOXO3",
  "MYC", "CEBPB", "RUNX1", "RUNX2", "ETS1", "TCF7L2", "LEF1"
)
CAF_KNOWN_TFS <- unique(CAF_KNOWN_TFS)

# ============================================================================
# 加载R包
# ============================================================================
cat("加载R包...\n")
required_packages <- c("ggplot2", "pheatmap", "dplyr", "tidyr", "RColorBrewer")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

if (!requireNamespace("dorothea", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("dorothea")
}
library(dorothea)
data("dorothea_hs", package = "dorothea")
cat("✓ 所有包加载完成\n\n")

# ============================================================================
# 1. DoRothEA regulons
# ============================================================================
cat("1. 构建DoRothEA regulons...\n")
all_tfs <- unique(dorothea_hs$tf)
dorothea_regulons <- dorothea_hs %>%
  filter(confidence %in% c('A', 'B', 'C')) %>%
  dplyr::select(tf, target, confidence, mor)

tf_target_list_dorothea <- dorothea_regulons %>%
  group_by(tf) %>%
  summarise(targets = list(unique(target)), .groups = 'drop') %>%
  filter(lengths(targets) >= 1)
tf_target_list_dorothea <- setNames(tf_target_list_dorothea$targets, tf_target_list_dorothea$tf)

# ============================================================================
# 2. 读取数据
# ============================================================================
cat("2. 读取表达数据...\n")
expr_data <- read.csv(EXPRESSION_FILE, stringsAsFactors = FALSE, check.names = FALSE)
rownames(expr_data) <- expr_data[, 1]; expr_data <- expr_data[, -1]
expression_matrix <- as.matrix(expr_data)
mode(expression_matrix) <- "numeric"
expression_matrix[is.na(expression_matrix)] <- 0

sample_names <- colnames(expression_matrix)
sample_groups <- data.frame(SampleID = sample_names, Group = NA, stringsAsFactors = FALSE)
for (prefix in names(GROUP_MAPPING)) {
  idx <- grep(paste0("^", prefix), sample_names)
  if (length(idx) > 0) sample_groups$Group[idx] <- GROUP_MAPPING[[prefix]]
}
sample_groups <- sample_groups[!is.na(sample_groups$Group), ]
sample_groups$Group <- factor(sample_groups$Group, levels = GROUP_ORDER)
expression_matrix <- expression_matrix[, sample_groups$SampleID, drop = FALSE]

# 预先构建按GROUP_ORDER排列的样品顺序（用于热图）
sample_groups_ordered <- sample_groups %>%
  mutate(Group = factor(Group, levels = GROUP_ORDER)) %>%
  arrange(Group)

# ============================================================================
# 3. 读取基因集
# ============================================================================
if (file.exists(GENE_SET_FILE)) {
  gene_set <- unique(read.csv(GENE_SET_FILE, stringsAsFactors = FALSE)[, 1])
} else {
  ctrl <- sample_groups$SampleID[sample_groups$Group == CONTROL_GROUP]
  trt <- sample_groups$SampleID[sample_groups$Group == EXPERIMENTAL_GROUPS[1]]
  logFC <- rowMeans(expression_matrix[, trt]) - rowMeans(expression_matrix[, ctrl])
  gene_set <- rownames(expression_matrix)[abs(logFC) >= 1]
}
gene_set <- intersect(gene_set, rownames(expression_matrix))
cat(sprintf("  基因集: %d 基因\n", length(gene_set)))

# ============================================================================
# 4. 筛选高变异TF
# ============================================================================
cat("3. 筛选TF...\n")
gene_vars <- apply(expression_matrix, 1, sd, na.rm = TRUE)
available_tfs <- intersect(all_tfs, rownames(expression_matrix))
if (length(available_tfs) > MAX_TFS_TO_ANALYZE) {
  tf_vars <- sort(gene_vars[available_tfs], decreasing = TRUE)
  available_tfs <- names(tf_vars)[1:MAX_TFS_TO_ANALYZE]
}

# ============================================================================
# 5. Spearman相关性
# ============================================================================
cat("4. Spearman相关性...\n")
gene_set_expr <- expression_matrix[gene_set, , drop = FALSE]
tf_expr <- expression_matrix[available_tfs, , drop = FALSE]
cor_matrix <- cor(t(tf_expr), t(gene_set_expr), method = "spearman", use = "pairwise.complete.obs")
n_samples <- ncol(expression_matrix)
t_stat <- cor_matrix * sqrt((n_samples - 2) / (1 - cor_matrix^2))
p_matrix <- 2 * pt(abs(t_stat), df = n_samples - 2, lower.tail = FALSE)

# ============================================================================
# 6. 自动阈值优化
# ============================================================================
cat("5. 阈值优化...\n")
if (AUTO_SELECT_THRESHOLD) {
  test_thresholds <- seq(THRESHOLD_MIN, THRESHOLD_MAX, by = THRESHOLD_STEP)
  best_threshold <- THRESHOLD_MIN; best_deviation <- Inf
  threshold_results <- data.frame()
  
  for (thr in test_thresholds) {
    sig <- abs(cor_matrix) >= thr & p_matrix < CORRELATION_PVALUE_THRESHOLD
    target_counts <- rowSums(sig, na.rm = TRUE)
    avg_targets <- mean(target_counts[target_counts > 0])
    med_targets <- median(target_counts[target_counts > 0])
    
    threshold_results <- rbind(threshold_results, data.frame(
      Threshold = thr, TFs = sum(target_counts > 0),
      Avg = avg_targets, Median = med_targets
    ))
    
    deviation <- abs(avg_targets - DESIRED_TARGETS_PER_TF)
    if (deviation < best_deviation) { best_deviation <- deviation; best_threshold <- thr }
  }
  correlation_threshold <- best_threshold
} else {
  correlation_threshold <- 0.85
}

# ============================================================================
# 7. 提取显著关系
# ============================================================================
cat("6. 提取显著关系...\n")
significant <- abs(cor_matrix) >= correlation_threshold & p_matrix < CORRELATION_PVALUE_THRESHOLD
tf_target_genes <- list(); tf_target_cors <- list(); valid_tfs <- c()

for (i in 1:nrow(significant)) {
  tf <- rownames(significant)[i]
  targets <- colnames(significant)[significant[i, ]]
  if (length(targets) >= MIN_TARGET_GENES) {
    tf_target_genes[[tf]] <- targets
    tf_target_cors[[tf]] <- cor_matrix[i, significant[i, ]]
    valid_tfs <- c(valid_tfs, tf)
  }
}

# ============================================================================
# 8. DoRothEA验证 + 综合评分 + 活性
# ============================================================================
cat("7. DoRothEA验证与评分...\n")
dorothea_target_counts <- sapply(valid_tfs, function(tf) {
  if (tf %in% names(tf_target_list_dorothea)) {
    return(length(intersect(tf_target_list_dorothea[[tf]], gene_set)))
  }
  return(0)
})
dorothea_target_genes <- sapply(valid_tfs, function(tf) {
  if (tf %in% names(tf_target_list_dorothea)) {
    return(paste(intersect(tf_target_list_dorothea[[tf]], gene_set), collapse = ";"))
  }
  return("")
})

tf_scores <- data.frame(
  TF = valid_tfs,
  Target_Count = sapply(valid_tfs, function(tf) length(tf_target_genes[[tf]])),
  Mean_Abs_Cor = sapply(valid_tfs, function(tf) mean(abs(tf_target_cors[[tf]]))),
  Dorothea_Targets = dorothea_target_counts,
  Dorothea_Target_Genes = dorothea_target_genes,
  Correlation_Target_Genes = sapply(valid_tfs, function(tf) {
    paste(tf_target_genes[[tf]], collapse = ";")
  }),
  stringsAsFactors = FALSE
)
tf_scores$Score <- tf_scores$Target_Count * tf_scores$Mean_Abs_Cor
tf_scores$Has_Dorothea <- tf_scores$Dorothea_Targets > 0
tf_scores$Literature_Known <- tf_scores$TF %in% CAF_KNOWN_TFS

# TF活性
tf_activity <- matrix(0, nrow = length(valid_tfs), ncol = ncol(expression_matrix))
rownames(tf_activity) <- valid_tfs
colnames(tf_activity) <- colnames(expression_matrix)
for (tf in valid_tfs) {
  targets <- tf_target_genes[[tf]]
  if (length(targets) > 0) {
    tf_activity[tf, ] <- colMeans(expression_matrix[targets, , drop = FALSE], na.rm = TRUE)
  }
}

tf_activity_by_group <- data.frame(TF = valid_tfs, stringsAsFactors = FALSE)
for (grp in GROUP_ORDER) {
  grp_samples <- sample_groups$SampleID[sample_groups$Group == grp]
  tf_activity_by_group[[grp]] <- rowMeans(tf_activity[valid_tfs, grp_samples, drop = FALSE], na.rm = TRUE)
}
for (grp in EXPERIMENTAL_GROUPS) {
  tf_activity_by_group[[paste0(grp, "_FC")]] <- 
    tf_activity_by_group[[grp]] / (tf_activity_by_group[[CONTROL_GROUP]] + 0.01)
}
tf_activity_by_group$Max_Activity_FC <- apply(
  tf_activity_by_group[, paste0(EXPERIMENTAL_GROUPS, "_FC"), drop = FALSE], 1, max, na.rm = TRUE
)

# ============================================================================
# 9. 筛选
# ============================================================================
cat("8. 筛选...\n")
tf_scores$Activity_Pass <- tf_scores$TF %in% 
  tf_activity_by_group$TF[tf_activity_by_group$Max_Activity_FC >= ACTIVITY_FC_THRESHOLD]

final_core_tfs <- tf_scores
if (ENABLE_ACTIVITY_FILTER) final_core_tfs <- final_core_tfs %>% filter(Activity_Pass)
final_core_tfs <- final_core_tfs %>%
  left_join(tf_activity_by_group[, c("TF", "Max_Activity_FC")], by = "TF")

# ============================================================================
# 10. 展示列表
# ============================================================================
cat("9. 构建展示列表...\n")

all_display <- final_core_tfs %>%
  filter(Has_Dorothea | Literature_Known) %>%
  mutate(
    Has_Dorothea_Int = as.integer(Has_Dorothea),
    Is_Literature_Int = as.integer(Literature_Known)
  ) %>%
  arrange(desc(Has_Dorothea_Int), desc(Dorothea_Targets), 
          desc(Is_Literature_Int), desc(Score), TF) %>%
  mutate(Source = case_when(
    Has_Dorothea & Literature_Known ~ "DoRothEA + Literature",
    Has_Dorothea ~ "DoRothEA-validated",
    TRUE ~ "Literature-reported"
  ))

main_display <- all_display %>% filter(Dorothea_Targets >= MIN_DOROTHEA_FOR_MAIN)
if (nrow(main_display) < 5) {
  main_display <- all_display %>% filter(Dorothea_Targets >= 2)
}

cat(sprintf("  主图TF (DoRothEA ≥ %d): %d 个\n", MIN_DOROTHEA_FOR_MAIN, nrow(main_display)))
cat(sprintf("  补充图TF (全部候选): %d 个\n", nrow(all_display)))

# ============================================================================
# 11. 保存表格
# ============================================================================
write.csv(tf_scores, file.path(output_dir, "Table_S1_All_TF_scores.csv"), row.names = FALSE)
write.csv(final_core_tfs, file.path(output_dir, "Table_S2_Final_candidate_TFs.csv"), row.names = FALSE)
write.csv(main_display, file.path(output_dir, "Table_S3_Main_Figure_TFs.csv"), row.names = FALSE)
write.csv(all_display, file.path(output_dir, "Table_S4_Supplementary_Figure_TFs.csv"), row.names = FALSE)

# 额外输出主图TF活性矩阵
tf_activity_main <- as.data.frame(tf_activity[main_display$TF, , drop = FALSE])
tf_activity_main <- cbind(TF = rownames(tf_activity_main), tf_activity_main)
write.csv(tf_activity_main, file.path(output_dir, "Table_S5_Main_TF_activity_matrix.csv"), row.names = FALSE)

# 额外输出补充图TF活性矩阵
tf_activity_supp <- as.data.frame(tf_activity[all_display$TF, , drop = FALSE])
tf_activity_supp <- cbind(TF = rownames(tf_activity_supp), tf_activity_supp)
write.csv(tf_activity_supp, file.path(output_dir, "Table_S6_Supplementary_TF_activity_matrix.csv"), row.names = FALSE)

# ============================================================================
# 辅助函数：绘制条形图（文献已知TF名称加粗）
# ============================================================================
draw_barplot <- function(display_df, fig_params, output_prefix) {
  display_tfs <- display_df$TF
  n_tfs <- length(display_tfs)
  
  plot_data <- display_df %>%
    mutate(
      TF = factor(TF, levels = rev(TF)),
      Unvalidated = Target_Count - Dorothea_Targets
    ) %>%
    dplyr::select(TF, Dorothea_Targets, Unvalidated, Literature_Known, Target_Count) %>%
    pivot_longer(cols = c(Dorothea_Targets, Unvalidated), 
                 names_to = "Type", values_to = "Count") %>%
    mutate(Type = factor(Type, levels = c("Unvalidated", "Dorothea_Targets"),
                         labels = c("Correlation-only", "DoRothEA-validated")))
  
  label_validated <- plot_data %>%
    group_by(TF) %>%
    arrange(TF, desc(Type)) %>%
    mutate(cumsum = cumsum(Count)) %>%
    filter(Type == "DoRothEA-validated", Count > 0)
  
  label_total <- plot_data %>%
    group_by(TF) %>%
    summarise(Total = sum(Count), .groups = 'drop') %>%
    mutate(TF = factor(TF, levels = levels(plot_data$TF)))
  
  subtitle <- paste0("Dark: DoRothEA-validated | Light: Correlation-only | Bold: CAF literature-known TF")
  
  p <- ggplot(plot_data, aes(x = TF, y = Count, fill = Type)) +
    geom_bar(stat = "identity", width = fig_params$bar_width, alpha = fig_params$bar_alpha) +
    coord_flip() +
    scale_fill_manual(values = c("DoRothEA-validated" = fig_params$color_validated,
                                 "Correlation-only" = fig_params$color_inferred),
                      name = "Target Type") +
    labs(title = fig_params$title, subtitle = subtitle, x = NULL, y = "Number of Target Genes") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = fig_params$title_face, size = fig_params$title_size),
          plot.subtitle = element_text(hjust = 0.5, color = fig_params$subtitle_color, size = fig_params$subtitle_size),
          axis.title = element_text(size = fig_params$axis_title_size),
          axis.text.y = element_text(size = fig_params$axis_text_size,
                                     face = ifelse(rev(display_df$Literature_Known), "bold", "plain")),
          axis.text.x = element_text(size = fig_params$axis_text_size),
          legend.position = "right")
  
  if (fig_params$show_values && nrow(label_validated) > 0) {
    p <- p + geom_text(data = label_validated,
                       aes(x = TF, y = cumsum, label = Count),
                       inherit.aes = FALSE, hjust = 1.1, color = "white",
                       size = fig_params$value_size)
  }
  if (fig_params$show_values && nrow(label_total) > 0) {
    p <- p + geom_text(data = label_total,
                       aes(x = TF, y = Total, label = Total),
                       inherit.aes = FALSE, hjust = 1.1, color = fig_params$value_color,
                       size = fig_params$value_size)
  }
  
  y_max <- ifelse(nrow(label_total) > 0, max(label_total$Total) * 1.08, 10)
  p <- p + scale_y_continuous(limits = c(0, y_max), expand = c(0, 0))
  
  height_adj <- max(fig_params$height, n_tfs * 0.5)
  
  ggsave(file.path(output_dir, paste0(output_prefix, "_barplot.pdf")), p, 
         width = fig_params$width, height = height_adj, dpi = FIG_DPI)
  ggsave(file.path(output_dir, paste0(output_prefix, "_barplot.png")), p, 
         width = fig_params$width, height = height_adj, dpi = FIG_DPI, bg = "white")
  
  cat(sprintf("  ✓ %s: 条形图 (%d TFs)\n", output_prefix, n_tfs))
  return(p)
}

# ============================================================================
# 辅助函数：绘制热图（全部4组，按GROUP_ORDER排列，无文献标注）
# ============================================================================
draw_heatmap <- function(display_df, matrix_to_plot, fig_params, output_prefix, plot_type) {
  display_tfs <- display_df$TF
  
  # 按 GROUP_ORDER 排列的样品顺序提取表达矩阵
  plot_data <- matrix_to_plot[display_tfs, sample_groups_ordered$SampleID, drop = FALSE]
  plot_scaled <- t(scale(t(plot_data)))
  plot_scaled[is.na(plot_scaled)] <- 0
  plot_scaled[plot_scaled > fig_params$scale_range] <- fig_params$scale_range
  plot_scaled[plot_scaled < -fig_params$scale_range] <- -fig_params$scale_range
  
  # 使用按 GROUP_ORDER 排列的分组注释
  anno_col <- data.frame(
    Group = sample_groups_ordered$Group,
    row.names = sample_groups_ordered$SampleID
  )
  
  height_adj <- max(fig_params$height, length(display_tfs) * ifelse(plot_type == "activity", 0.3, 0.2))
  
  pdf(file.path(output_dir, paste0(output_prefix, "_", plot_type, "_heatmap.pdf")), 
      width = fig_params$width, height = height_adj)
  pheatmap(plot_scaled,
           color = colorRampPalette(c(fig_params$color_low, fig_params$color_mid, fig_params$color_high))(100),
           annotation_col = anno_col,
           annotation_colors = fig2_annotation_colors,
           cluster_rows = fig_params$cluster_rows, cluster_cols = fig_params$cluster_cols,
           show_rownames = fig_params$show_rownames, show_colnames = fig_params$show_colnames,
           fontsize_row = fig_params$fontsize_row, fontsize_col = fig_params$fontsize_col,
           main = fig_params$title)
  dev.off()
  
  png(file.path(output_dir, paste0(output_prefix, "_", plot_type, "_heatmap.png")), 
      width = fig_params$width, height = height_adj, units = "in", res = FIG_DPI, bg = "white")
  pheatmap(plot_scaled,
           color = colorRampPalette(c(fig_params$color_low, fig_params$color_mid, fig_params$color_high))(100),
           annotation_col = anno_col,
           annotation_colors = fig2_annotation_colors,
           cluster_rows = fig_params$cluster_rows, cluster_cols = fig_params$cluster_cols,
           show_rownames = fig_params$show_rownames, show_colnames = fig_params$show_colnames,
           fontsize_row = fig_params$fontsize_row, fontsize_col = fig_params$fontsize_col,
           main = fig_params$title)
  dev.off()
  
  cat(sprintf("  ✓ %s: %s热图 (%d TFs)\n", output_prefix, plot_type, length(display_tfs)))
}

# ============================================================================
# ★★★ 生成所有图片 ★★★
# ============================================================================
cat("\n10. 生成图片...\n\n")

# ---- 阈值优化图 ----
if (exists("threshold_results") && nrow(threshold_results) > 0) {
  p0 <- ggplot(threshold_results, aes(x = Threshold)) +
    geom_line(aes(y = Avg, color = "Average"), size = fig0_line_size) +
    geom_line(aes(y = Median, color = "Median"), size = fig0_line_size, linetype = "dashed") +
    geom_hline(yintercept = DESIRED_TARGETS_PER_TF, linetype = "dashed", 
               color = fig0_target_line_color, alpha = 0.7) +
    geom_vline(xintercept = best_threshold, linetype = "dotted", 
               color = fig0_selected_line_color, size = 1) +
    scale_color_manual(values = c("Average" = fig0_avg_line_color, "Median" = fig0_median_line_color)) +
    labs(title = fig0_title, x = "Correlation Threshold (|r|)", y = "Target Genes per TF", color = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = fig0_title_face, size = fig0_title_size),
          axis.title = element_text(size = fig0_axis_title_size),
          axis.text = element_text(size = fig0_axis_text_size),
          legend.position = "bottom")
  ggsave(file.path(output_dir, "Figure0_Threshold_Optimization.pdf"), p0, 
         width = fig0_width, height = fig0_height, dpi = FIG_DPI)
  ggsave(file.path(output_dir, "Figure0_Threshold_Optimization.png"), p0, 
         width = fig0_width, height = fig0_height, dpi = FIG_DPI, bg = "white")
  cat("✓ Figure 0: 阈值优化图\n")
}

# ---- 主图 ----
if (nrow(main_display) > 0) {
  main_bar_params <- list(
    width = fig1_width, height = fig1_height,
    title = fig1_title, title_size = fig1_title_size, title_face = fig1_title_face,
    subtitle_size = fig1_subtitle_size, subtitle_color = fig1_subtitle_color,
    axis_title_size = fig1_axis_title_size, axis_text_size = fig1_axis_text_size,
    bar_width = fig1_bar_width, bar_alpha = fig1_bar_alpha,
    color_validated = fig1_color_validated, color_inferred = fig1_color_inferred,
    show_values = fig1_show_values, value_size = fig1_value_size, value_color = fig1_value_color
  )
  draw_barplot(main_display, main_bar_params, "Figure1_Main")
  
  main_expr_params <- list(
    width = fig2_width, height = fig2_height,
    title = fig2_title,
    color_low = fig2_color_low, color_mid = fig2_color_mid, color_high = fig2_color_high,
    fontsize_row = fig2_fontsize_row, fontsize_col = fig2_fontsize_col,
    cluster_rows = fig2_cluster_rows, cluster_cols = fig2_cluster_cols,
    show_rownames = fig2_show_rownames, show_colnames = fig2_show_colnames,
    scale_range = fig2_scale_range
  )
  draw_heatmap(main_display, expression_matrix, main_expr_params, "Figure2_Main", "expression")
  
  main_act_params <- list(
    width = fig3_width, height = fig3_height,
    title = fig3_title,
    color_low = fig3_color_low, color_mid = fig3_color_mid, color_high = fig3_color_high,
    fontsize_row = fig3_fontsize_row, fontsize_col = fig3_fontsize_col,
    cluster_rows = fig3_cluster_rows, cluster_cols = fig3_cluster_cols,
    show_rownames = fig3_show_rownames, show_colnames = fig3_show_colnames,
    scale_range = fig3_scale_range
  )
  draw_heatmap(main_display, tf_activity, main_act_params, "Figure3_Main", "activity")
}

# ---- 补充图 ----
if (nrow(all_display) > 0) {
  supp_bar_params <- list(
    width = figS1_width, height = figS1_height,
    title = figS1_title, title_size = figS1_title_size, title_face = figS1_title_face,
    subtitle_size = figS1_subtitle_size, subtitle_color = figS1_subtitle_color,
    axis_title_size = figS1_axis_title_size, axis_text_size = figS1_axis_text_size,
    bar_width = figS1_bar_width, bar_alpha = figS1_bar_alpha,
    color_validated = figS1_color_validated, color_inferred = figS1_color_inferred,
    show_values = figS1_show_values, value_size = figS1_value_size, value_color = figS1_value_color
  )
  draw_barplot(all_display, supp_bar_params, "FigureS1_Supplementary")
  
  supp_expr_params <- list(
    width = figS2_width, height = figS2_height,
    title = figS2_title,
    color_low = figS2_color_low, color_mid = figS2_color_mid, color_high = figS2_color_high,
    fontsize_row = figS2_fontsize_row, fontsize_col = figS2_fontsize_col,
    cluster_rows = figS2_cluster_rows, cluster_cols = figS2_cluster_cols,
    show_rownames = figS2_show_rownames, show_colnames = figS2_show_colnames,
    scale_range = figS2_scale_range
  )
  draw_heatmap(all_display, expression_matrix, supp_expr_params, "FigureS2_Supplementary", "expression")
  
  supp_act_params <- list(
    width = figS3_width, height = figS3_height,
    title = figS3_title,
    color_low = figS3_color_low, color_mid = figS3_color_mid, color_high = figS3_color_high,
    fontsize_row = figS3_fontsize_row, fontsize_col = figS3_fontsize_col,
    cluster_rows = figS3_cluster_rows, cluster_cols = figS3_cluster_cols,
    show_rownames = figS3_show_rownames, show_colnames = figS3_show_colnames,
    scale_range = figS3_scale_range
  )
  draw_heatmap(all_display, tf_activity, supp_act_params, "FigureS3_Supplementary", "activity")
}

# ============================================================================
# 12. 分析报告
# ============================================================================
sink(file.path(output_dir, "Analysis_Report.txt"))

cat("================================================================================\n")
cat("        转录重编程基因集核心转录因子分析报告\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("随机种子: 42\n\n")

cat("一、方法学\n")
cat("--------------------------------------------------------------------------------\n")
cat("  1. TF来源: DoRothEA全基因组数据库 (", length(all_tfs), "个TF)\n", sep="")
cat("  2. 靶基因推断: Spearman秩相关, |r| > ", correlation_threshold, ", P < ", CORRELATION_PVALUE_THRESHOLD, "\n", sep="")
cat("  3. 最少靶基因数: ≥", MIN_TARGET_GENES, "\n", sep="")
cat("  4. 靶基因验证: DoRothEA regulons (A/B/C级)\n")
cat("  5. TF活性: 靶基因平均表达\n")
cat("  6. 主图: DoRothEA靶基因 ≥ ", MIN_DOROTHEA_FOR_MAIN, " | 补充图: 全部候选\n\n", sep="")

cat("二、结果统计\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("  最终通过筛选: %d 个TF\n", nrow(final_core_tfs)))
cat(sprintf("  DoRothEA验证: %d (%.1f%%)\n", sum(final_core_tfs$Has_Dorothea), mean(final_core_tfs$Has_Dorothea)*100))
cat(sprintf("  文献已知: %d\n", sum(final_core_tfs$Literature_Known)))
cat(sprintf("  主图展示 (≥%d): %d\n", MIN_DOROTHEA_FOR_MAIN, nrow(main_display)))
cat(sprintf("  补充图展示: %d\n\n", nrow(all_display)))

cat("三、主图TF\n")
cat("--------------------------------------------------------------------------------\n")
for (i in 1:nrow(main_display)) {
  cat(sprintf("  %2d. %s: DoRothEA=%d, 总靶基因=%d\n",
              i, main_display$TF[i], main_display$Dorothea_Targets[i], main_display$Target_Count[i]))
}

cat("\n四、输出文件\n")
cat("--------------------------------------------------------------------------------\n")
cat("  Figure0 - 阈值优化图\n")
cat("  Figure1_Main_barplot - 主图条形图 (Bold = CAF文献已知)\n")
cat("  Figure2_Main_expression_heatmap - 主图表达热图 (DMEM/Con/I6/R2)\n")
cat("  Figure3_Main_activity_heatmap - 主图活性热图 (DMEM/Con/I6/R2)\n")
cat("  FigureS1_Supplementary_barplot - 补充条形图\n")
cat("  FigureS2_Supplementary_expression_heatmap - 补充表达热图\n")
cat("  FigureS3_Supplementary_activity_heatmap - 补充活性热图\n")
cat("  Table_S1 - 所有TF评分\n")
cat("  Table_S2 - 最终候选TF\n")
cat("  Table_S3 - 主图展示TF\n")
cat("  Table_S4 - 补充图展示TF\n")
cat("  Table_S5 - 主图TF活性矩阵\n")
cat("  Table_S6 - 补充图TF活性矩阵\n")

cat("\n五、引用\n")
cat("--------------------------------------------------------------------------------\n")
cat("  DoRothEA: Garcia-Alonso L, et al. Genome Research, 2019\n")
cat("  VIPER: Alvarez MJ, et al. Nature Genetics, 2016\n")

sink()

cat("\n========================================\n")
cat("分析完成！\n")
cat("输出目录:", output_dir, "\n")
cat("========================================\n")