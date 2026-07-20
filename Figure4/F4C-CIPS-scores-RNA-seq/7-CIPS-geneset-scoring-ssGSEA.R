# ============================================================================
# 基因集表达评分分析（ssGSEA）- Student's t-test版本
# 功能：使用GSEA官方ssGSEA算法计算表达评分，箱线图展示各组差异
# 统计方法：双尾独立样本Student's t检验（假设方差齐性）
# ============================================================================

# 清空环境变量
rm(list = ls())

# ========================== 工作路径设定 ====================================
work_dir <- "D:/zsy/SX/Fomal-final/7-CIPS-heatmap"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")

# ========================== 创建输出目录 ====================================
output_dir <- "gene_set_scoring_ssGSEA"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# ============================================================================
# ====================== 用户自定义参数区域 ==================================
# ============================================================================

# ----------------------------- 输入文件 -------------------------------------
expression_file <- "gene_expression_matrix.csv"
gene_list_file <- "gene_pathway_stats.csv"
gene_id_col <- "GeneID"
list_gene_col <- "GeneID"
gene_set_name <- "Target_GeneSet"

# ----------------------------- 分组映射 -------------------------------------
group_mapping <- list(
  "Con" = "Con",
  "I6"  = "I6",
  "R2"  = "R2"
)

group_order <- c("Con", "I6", "R2")

group_colors <- c(
  "Con" = "#25377F",
  "I6"  = "#C60036",
  "R2"  = "#E4945A"
)

# ----------------------------- ssGSEA参数 ----------------------------------
ssgsea_alpha <- 0.25           # GSEA权重参数（0.25为默认值）
perform_zscore <- TRUE         # 是否对评分进行Z-score标准化

# ----------------------------- 箱线图参数 -----------------------------------
boxplot_width <- 6
boxplot_height <- 5.5
boxplot_title <- paste0("Gene Set Enrichment Score (", gene_set_name, ")")
boxplot_title_size <- 20
boxplot_x_label <- "Group"
boxplot_y_label <- "ssGSEA Enrichment Score"
boxplot_axis_title_size <- 22
boxplot_axis_text_size <- 20
boxplot_show_points <- TRUE
boxplot_point_size <- 3
boxplot_point_alpha <- 0.6
boxplot_jitter_width <- 0.3
boxplot_box_width <- 0.5
boxplot_box_alpha <- 0.7

# ----------------------------- Y轴范围参数 -----------------------------------
y_axis_limits <- NULL

# ----------------------------- p值标注参数 -----------------------------------
show_pvalue_labels <- TRUE
pval_label_size <- 6
pval_label_color <- "black"
pval_line_color <- "black"
pval_line_size <- 0.5
pval_bracket_height <- 0.05
pval_text_y_offset <- 0.08
pval_step_increase <- 0.17

# ----------------------------- 统计检验参数 ---------------------------------
stat_test_method <- "t.test"       # "t.test" 或 "wilcox.test"
p_adjust_method <- "BH"            # "BH", "bonferroni", "holm", "none"
significance_level <- 0.05
show_p_value_type <- "adjusted"    # "adjusted" 或 "raw"

# ----------------------------- 随机种子 -------------------------------------
random_seed <- 42
set.seed(random_seed)

# ----------------------------- 输出文件 -------------------------------------
output_score_file <- file.path(output_dir, paste0("ssgsea_scores_", gene_set_name, ".csv"))
output_boxplot <- file.path(output_dir, paste0("ssgsea_boxplot_", gene_set_name, ".pdf"))
output_summary <- file.path(output_dir, "ssgsea_analysis_summary.txt")
output_session_info <- file.path(output_dir, "session_info.txt")

# ============================================================================
# ====================== 加载必要的包 ========================================
# ============================================================================

cat("\n========== 加载必要的包 ==========\n")

required_packages <- c("dplyr", "ggplot2", "reshape2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
  cat(sprintf("  ✓ %s (v%s)\n", pkg, as.character(packageVersion(pkg))))
}

cat("✓ 所有包已加载\n")

# ============================================================================
# ====================== ssGSEA算法实现 =====================================
# ============================================================================

cat("\n========== 定义ssGSEA算法 ==========\n")

ssgsea_official <- function(expression_matrix, gene_set, alpha = 0.25, verbose = TRUE) {
  # =========================================================================
  # 基于GSEA官方算法的ssGSEA实现
  # 参考文献: 
  #   Subramanian et al., PNAS 2005
  #   Barbie et al., Nature 2009
  # =========================================================================
  
  n_genes <- nrow(expression_matrix)
  n_samples <- ncol(expression_matrix)
  sample_names <- colnames(expression_matrix)
  
  gs_idx <- which(rownames(expression_matrix) %in% gene_set)
  n_gs <- length(gs_idx)
  
  if (verbose) {
    cat(sprintf("总基因数: %d\n", n_genes))
    cat(sprintf("基因集大小: %d\n", n_gs))
    cat(sprintf("样品数: %d\n", n_samples))
    cat(sprintf("alpha参数: %.2f\n", alpha))
  }
  
  if (n_gs == 0) stop("没有基因集基因在表达矩阵中")
  
  scores <- numeric(n_samples)
  names(scores) <- sample_names
  
  for (s in 1:n_samples) {
    expr <- expression_matrix[, s]
    order_idx <- order(expr, decreasing = TRUE)
    ranked_expr <- expr[order_idx]
    ranked_genes <- rownames(expression_matrix)[order_idx]
    
    gs_positions <- which(ranked_genes %in% gene_set)
    
    if (length(gs_positions) == 0) {
      scores[s] <- 0
      next
    }
    
    N <- n_genes
    Nh <- length(gs_positions)
    
    es <- 0
    max_es <- 0
    min_es <- 0
    
    for (i in 1:N) {
      if (i %in% gs_positions) {
        weight <- abs(ranked_expr[i]) ^ alpha
        es <- es + weight
      } else {
        es <- es - (1 / (N - Nh))
      }
      if (es > max_es) max_es <- es
      if (es < min_es) min_es <- es
    }
    
    scores[s] <- if (abs(max_es) >= abs(min_es)) max_es else min_es
    
    if (verbose && s <= 3) {
      cat(sprintf("\n样品 %s: ES = %.6f\n", sample_names[s], scores[s]))
    }
  }
  
  if (max(scores, na.rm = TRUE) > min(scores, na.rm = TRUE)) {
    scores_norm <- 2 * (scores - min(scores, na.rm = TRUE)) / 
      (max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE)) - 1
    if (verbose) cat("\n分数归一化到[-1, 1]区间\n")
  } else {
    scores_norm <- scores
  }
  
  names(scores_norm) <- sample_names
  return(scores_norm)
}

cat("✓ ssGSEA算法定义完成\n")

# ========================== 读取表达矩阵 ====================================
cat("\n========== 读取表达矩阵 ==========\n")

if (!file.exists(expression_file)) {
  stop("错误：未找到表达矩阵文件: ", expression_file)
}

expr_data <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("表达矩阵维度: %d 行, %d 列\n", nrow(expr_data), ncol(expr_data)))

if (gene_id_col %in% colnames(expr_data)) {
  rownames(expr_data) <- expr_data[[gene_id_col]]
  expr_data <- expr_data[, !colnames(expr_data) %in% gene_id_col]
} else {
  rownames(expr_data) <- expr_data[, 1]
  expr_data <- expr_data[, -1]
}

expression_matrix <- as.matrix(expr_data)
mode(expression_matrix) <- "numeric"

if (any(is.na(expression_matrix))) {
  cat(sprintf("警告：发现 %d 个缺失值，已替换为0\n", sum(is.na(expression_matrix))))
  expression_matrix[is.na(expression_matrix)] <- 0
}

cat(sprintf("表达矩阵维度: %d 基因 × %d 样品\n", nrow(expression_matrix), ncol(expression_matrix)))

# ========================== 读取基因列表 ====================================
cat("\n========== 读取基因列表 ==========\n")

if (!file.exists(gene_list_file)) {
  stop("错误：未找到基因列表文件: ", gene_list_file)
}

gene_data <- read.csv(gene_list_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("基因列表维度: %d 行, %d 列\n", nrow(gene_data), ncol(gene_data)))

if (list_gene_col %in% colnames(gene_data)) {
  gene_set <- gene_data[[list_gene_col]]
} else {
  gene_set <- gene_data[, 1]
}

gene_set <- toupper(trimws(gene_set))
gene_set <- unique(gene_set[!is.na(gene_set) & gene_set != ""])

cat(sprintf("提供基因集大小: %d\n", length(gene_set)))

# ========================== 筛选匹配基因 ====================================
cat("\n========== 筛选匹配基因 ==========\n")

rownames(expression_matrix) <- toupper(rownames(expression_matrix))
available_genes <- intersect(gene_set, rownames(expression_matrix))
missing_genes <- setdiff(gene_set, rownames(expression_matrix))

cat(sprintf("匹配的基因数: %d / %d (%.1f%%)\n", 
            length(available_genes), length(gene_set),
            ifelse(length(gene_set) > 0, length(available_genes)/length(gene_set)*100, 0)))

if (length(available_genes) == 0) stop("错误：没有匹配的基因！")

if (length(missing_genes) > 0) {
  cat(sprintf("注意：%d 个基因未在表达矩阵中找到\n", length(missing_genes)))
}

# ========================== 运行ssGSEA分析 ==================================
cat("\n========== 运行ssGSEA分析 ==========\n")

gsva_start_time <- Sys.time()

ssgsea_scores <- ssgsea_official(
  expression_matrix = expression_matrix,
  gene_set = available_genes,
  alpha = ssgsea_alpha,
  verbose = TRUE
)

gsva_end_time <- Sys.time()
cat(sprintf("\n✓ ssGSEA分析完成，耗时: %.2f 秒\n", 
            difftime(gsva_end_time, gsva_start_time, units = "secs")))

cat(sprintf("原始ssGSEA评分范围: [%.4f, %.4f]\n", 
            min(ssgsea_scores, na.rm = TRUE), max(ssgsea_scores, na.rm = TRUE)))

# ========================== 提取评分数据 ====================================
cat("\n========== 提取评分数据 ==========\n")

scores <- ssgsea_scores

if (perform_zscore && sd(scores, na.rm = TRUE) > 0) {
  cat("执行Z-score标准化...\n")
  scores <- as.numeric(scale(scores))
  names(scores) <- names(ssgsea_scores)
  cat(sprintf("标准化后评分范围: [%.4f, %.4f]\n", 
              min(scores, na.rm = TRUE), max(scores, na.rm = TRUE)))
}

# ========================== 创建评分数据框 ====================================
cat("\n========== 创建评分数据框 ==========\n")

score_df <- data.frame(
  Sample = names(scores),
  Score = as.numeric(scores),
  stringsAsFactors = FALSE
)

score_df$Group <- NA
for (prefix in names(group_mapping)) {
  matching_idx <- grep(paste0("^", prefix), score_df$Sample)
  if (length(matching_idx) > 0) {
    score_df$Group[matching_idx] <- group_mapping[[prefix]]
  }
}

score_df <- score_df[!is.na(score_df$Group), ]
score_df$Group <- factor(score_df$Group, levels = group_order)

cat("分组信息:\n")
print(table(score_df$Group))

write.csv(score_df, output_score_file, row.names = FALSE)
cat(sprintf("✓ 评分结果已保存: %s\n", output_score_file))

# ========================== 评分统计 ====================================
cat("\n========== 计算评分统计 ==========\n")

score_stats <- score_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    Median = median(Score, na.rm = TRUE),
    IQR = IQR(Score, na.rm = TRUE),
    .groups = "drop"
  )
print(score_stats)

# ========================== 统计检验（Student's t-test，假设方差齐性）===========
cat("\n========== 计算统计检验（双尾Student's t-test，含多重比较校正）==========\n")

groups <- levels(score_df$Group)
pairwise_pvalues <- data.frame()

for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    group1_data <- score_df$Score[score_df$Group == groups[i]]
    group2_data <- score_df$Score[score_df$Group == groups[j]]
    
    if (length(group1_data) >= 2 && length(group2_data) >= 2) {
      if (stat_test_method == "t.test") {
        # 使用 var.equal = TRUE（Student's t-test，假设方差齐性）
        test_result <- t.test(group1_data, group2_data, var.equal = TRUE)
        test_name <- "Student's t-test"
      } else if (stat_test_method == "wilcox.test") {
        test_result <- wilcox.test(group1_data, group2_data, exact = FALSE)
        test_name <- "Mann-Whitney U test"
      }
      
      pairwise_pvalues <- rbind(pairwise_pvalues, data.frame(
        Group1 = groups[i],
        Group2 = groups[j],
        Comparison = paste(groups[i], "vs", groups[j]),
        Raw_P_value = test_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 多重比较校正
if (nrow(pairwise_pvalues) > 0) {
  if (p_adjust_method != "none" && nrow(pairwise_pvalues) > 1) {
    pairwise_pvalues$Adjusted_P_value <- p.adjust(pairwise_pvalues$Raw_P_value, 
                                                  method = p_adjust_method)
    cat(sprintf("使用 %s 方法进行多重比较校正\n", p_adjust_method))
  } else {
    pairwise_pvalues$Adjusted_P_value <- pairwise_pvalues$Raw_P_value
    cat("未进行多重比较校正\n")
  }
  
  if (show_p_value_type == "adjusted") {
    pairwise_pvalues$Display_P_value <- pairwise_pvalues$Adjusted_P_value
  } else {
    pairwise_pvalues$Display_P_value <- pairwise_pvalues$Raw_P_value
  }
  
  pairwise_pvalues$P_label <- sapply(pairwise_pvalues$Display_P_value, function(p) {
    if (p < 0.0001) return("p < 0.0001")
    if (p < 0.001) return(sprintf("p = %.4f", p))
    if (p < 0.01) return(sprintf("p = %.4f", p))
    return(sprintf("p = %.3f", p))
  })
  
  pairwise_pvalues$Significant <- pairwise_pvalues$Display_P_value < significance_level
  
  pairwise_pvalues$Significance <- ifelse(pairwise_pvalues$Display_P_value < 0.001, "***",
                                          ifelse(pairwise_pvalues$Display_P_value < 0.01, "**",
                                                 ifelse(pairwise_pvalues$Display_P_value < 0.05, "*", "ns")))
}

cat(sprintf("\n两两比较结果（%s，显示%s值）:\n", 
            test_name,
            ifelse(show_p_value_type == "adjusted", "校正后p", "原始p")))
print(pairwise_pvalues[, c("Comparison", "P_label", "Significance")])

# ========================== 计算Y轴范围 ====================================
if (is.null(y_axis_limits)) {
  y_max <- max(score_df$Score, na.rm = TRUE)
  y_min <- min(score_df$Score, na.rm = TRUE)
  y_range <- y_max - y_min
  y_range <- ifelse(y_range == 0, 1, y_range)
  
  n_comparisons <- nrow(pairwise_pvalues)
  y_max_expanded <- y_max + y_range * (0.2 + n_comparisons * pval_step_increase)
  y_axis_limits <- c(y_min, y_max_expanded)
  cat(sprintf("自动计算Y轴范围: [%.2f, %.2f]\n", y_axis_limits[1], y_axis_limits[2]))
}

# ========================== 绘制箱线图 ====================================
cat("\n========== 绘制箱线图 ==========\n")

p <- ggplot(score_df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(width = boxplot_box_width, 
               alpha = boxplot_box_alpha,
               outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", 
               width = boxplot_box_width * 0.5,
               size = 0.3) +
  scale_fill_manual(values = group_colors, name = "Group") +
  labs(title = boxplot_title, 
       x = boxplot_x_label, 
       y = boxplot_y_label) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = boxplot_title_size, face = "bold"),
    axis.title = element_text(size = boxplot_axis_title_size),
    axis.text = element_text(size = boxplot_axis_text_size, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = ifelse(length(group_order) > 3, "right", "none")
  ) +
  coord_cartesian(ylim = y_axis_limits, clip = "off")

if (boxplot_show_points) {
  p <- p + 
    geom_jitter(width = boxplot_jitter_width, 
                size = boxplot_point_size, 
                alpha = boxplot_point_alpha, 
                aes(color = Group),
                shape = 16) +
    scale_color_manual(values = group_colors, guide = "none")
}

# ========================== 添加p值标注 ====================================
if (show_pvalue_labels && nrow(pairwise_pvalues) > 0) {
  
  x_positions <- 1:length(group_order)
  names(x_positions) <- group_order
  
  y_max_data <- max(score_df$Score, na.rm = TRUE)
  y_range_data <- diff(range(score_df$Score, na.rm = TRUE))
  y_range_data <- ifelse(y_range_data == 0, 1, y_range_data)
  
  for (i in 1:nrow(pairwise_pvalues)) {
    comp <- pairwise_pvalues[i, ]
    
    x1 <- x_positions[comp$Group1]
    x2 <- x_positions[comp$Group2]
    
    y_position <- y_max_data + y_range_data * (pval_bracket_height + i * pval_step_increase)
    
    line_color <- ifelse(comp$Significant, pval_line_color, "grey60")
    text_color <- ifelse(comp$Significant, pval_label_color, "grey60")
    
    p <- p +
      annotate("segment", x = x1, xend = x1, 
               y = y_position - y_range_data * pval_bracket_height * 0.5, 
               yend = y_position,
               color = line_color, size = pval_line_size) +
      annotate("segment", x = x2, xend = x2, 
               y = y_position - y_range_data * pval_bracket_height * 0.5, 
               yend = y_position,
               color = line_color, size = pval_line_size) +
      annotate("segment", x = x1, xend = x2, 
               y = y_position, yend = y_position,
               color = line_color, size = pval_line_size) +
      annotate("text", x = (x1 + x2) / 2, 
               y = y_position + y_range_data * pval_text_y_offset,
               label = comp$P_label, 
               size = pval_label_size, 
               color = text_color,
               fontface = ifelse(comp$Significant, "bold", "plain"))
  }
}

# ========================== 保存图片 ====================================
ggsave(output_boxplot, p, width = boxplot_width, height = boxplot_height, dpi = 300)
cat(sprintf("✓ 箱线图已保存: %s\n", output_boxplot))

ggsave(file.path(output_dir, paste0("ssgsea_boxplot_", gene_set_name, ".png")), 
       p, width = boxplot_width, height = boxplot_height, dpi = 300)
cat(sprintf("✓ PNG预览图已保存\n"))

# ========================== 保存会话信息 ====================================
sink(output_session_info)
cat("================ R会话信息 ================\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("随机种子:", random_seed, "\n")
cat("ssGSEA算法: GSEA官方实现\n")
cat("alpha参数:", ssgsea_alpha, "\n")
cat("统计检验: 双尾独立样本Student's t-test (var.equal = TRUE)\n\n")
sessionInfo()
sink()

cat(sprintf("✓ 会话信息已保存: %s\n", output_session_info))

# ========================== 生成分析报告 ====================================
cat("\n========== 生成分析报告 ==========\n")

sink(output_summary)

cat("================================================================================\n")
cat("             基因集富集评分分析（ssGSEA）报告\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("随机种子:", random_seed, "\n\n")

cat("一、输入数据概况\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("表达矩阵维度: %d 基因 × %d 样品\n", nrow(expression_matrix), ncol(expression_matrix)))
cat(sprintf("基因集名称: %s\n", gene_set_name))
cat(sprintf("提供基因数: %d\n", length(gene_set)))
cat(sprintf("匹配基因数: %d (%.1f%%)\n", 
            length(available_genes), 
            ifelse(length(gene_set) > 0, length(available_genes)/length(gene_set)*100, 0)))

cat("\n二、ssGSEA分析参数\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("算法: 基于GSEA官方算法的ssGSEA (Subramanian et al., PNAS 2005)\n"))
cat(sprintf("alpha参数: %.2f\n", ssgsea_alpha))
cat(sprintf("评分归一化: [-1, 1]区间\n"))
cat(sprintf("Z-score标准化: %s\n", ifelse(perform_zscore, "是", "否")))

cat("\n三、分组信息\n")
cat("--------------------------------------------------------------------------------\n")
print(table(score_df$Group))

cat("\n四、ssGSEA评分统计\n")
cat("--------------------------------------------------------------------------------\n")
print(as.data.frame(score_stats))

cat("\n五、统计检验结果\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("检验方法: 双尾独立样本Student's t-test\n"))
cat(sprintf("多重比较校正方法: %s\n", p_adjust_method))
cat(sprintf("显著性阈值: p < %s\n", significance_level))
cat(sprintf("显示值类型: %s\n\n", show_p_value_type))

if (nrow(pairwise_pvalues) > 0) {
  display_cols <- c("Comparison", "Raw_P_value", "Adjusted_P_value", "P_label", 
                    "Significance", "Significant")
  display_cols <- intersect(display_cols, colnames(pairwise_pvalues))
  print(pairwise_pvalues[, display_cols, drop = FALSE])
}

cat("\n六、引用信息\n")
cat("--------------------------------------------------------------------------------\n")
cat("ssGSEA算法请引用:\n")
cat("1. Subramanian A, et al. (2005) Gene set enrichment analysis. PNAS, 102:15545-15550.\n")
cat("2. Barbie DA, et al. (2009) Systematic RNA interference reveals that oncogenic\n")
cat("   KRAS-driven cancers require TBK1. Nature, 462:108-112.\n")

cat("\n================================================================================\n")
cat("                            报告结束\n")
cat("================================================================================\n")

sink()

cat(sprintf("\n✓ 分析报告已保存: %s\n", output_summary))

cat("\n========================================\n")
cat("ssGSEA分析完成！\n")
cat("输出目录:", output_dir, "\n")
cat("========================================\n")

if (interactive()) {
  print(p)
}