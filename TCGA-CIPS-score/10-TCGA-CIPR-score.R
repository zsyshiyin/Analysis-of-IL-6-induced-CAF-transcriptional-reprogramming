# =============================================================================
# 生信分析脚本：ssGSEA基因集富集分析
# 功能：计算TCGA数据库中各样品的基因集表达评分
# 输入：mRNA.txt（表达矩阵，第一列基因名，第一行样品名）
#       gene_set.txt（基因列表）
# 输出：ssGSEA评分矩阵
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/10-TCGA-CIPR-score"  # 请修改此路径
setwd(work_dir)

expression_file <- "mRNA.txt"
gene_set_file <- "gene_set.txt"
group_file <- "group_info.csv"

group_sample_col <- "Sample"
group_group_col <- "Group"

output_file <- "ssGSEA_scores_diagnostic.csv"
output_plot <- "ssGSEA_distribution_diagnostic.pdf"

# ====================== 加载包 ======================
cat("\n========== 加载包 ==========\n")

library(ggplot2)
library(dplyr)

# ====================== 数据读取 ======================
cat("\n========== 读取表达矩阵 ==========\n")

first_line <- readLines(expression_file, n = 1)
if (grepl("\t", first_line)) {
  expr_data_raw <- read.table(expression_file, header = TRUE, 
                              stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
} else {
  expr_data_raw <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)
}

cat(sprintf("原始数据维度: %d 行, %d 列\n", nrow(expr_data_raw), ncol(expr_data_raw)))

# 提取基因名和表达矩阵
gene_names <- expr_data_raw[, 1]
expr_data <- expr_data_raw[, -1, drop = FALSE]
rownames(expr_data) <- gene_names

# 转换为数值矩阵
expr_matrix <- as.matrix(expr_data)
mode(expr_matrix) <- "numeric"

cat(sprintf("表达矩阵维度: %d 基因, %d 样品\n", nrow(expr_matrix), ncol(expr_matrix)))

# ====================== 读取基因集 ======================
cat("\n========== 读取基因集 ==========\n")

geneset_genes <- readLines(gene_set_file)
geneset_genes <- geneset_genes[!grepl("^#", geneset_genes) & geneset_genes != ""]
geneset_genes <- toupper(trimws(geneset_genes))

cat(sprintf("基因集大小: %d\n", length(geneset_genes)))

# 匹配基因
rownames(expr_matrix) <- toupper(rownames(expr_matrix))
available_genes <- intersect(geneset_genes, rownames(expr_matrix))
cat(sprintf("匹配到的基因数: %d\n", length(available_genes)))

if (length(available_genes) == 0) stop("没有匹配的基因！")

# ====================== 读取分组信息 ======================
cat("\n========== 读取分组信息 ==========\n")

has_group <- FALSE
group_info <- NULL

if (file.exists(group_file)) {
  group_info <- read.csv(group_file, stringsAsFactors = FALSE)
  colnames(group_info) <- toupper(colnames(group_info))
  sample_col <- toupper(group_sample_col)
  group_col <- toupper(group_group_col)
  
  if (sample_col %in% colnames(group_info) && group_col %in% colnames(group_info)) {
    group_info[[sample_col]] <- toupper(gsub("[[:space:]]", "", group_info[[sample_col]]))
    has_group <- TRUE
    cat("分组信息加载成功\n")
    print(table(group_info[[group_col]]))
  }
}

# ====================== 改进的ssGSEA实现 ======================
cat("\n========== 执行改进版ssGSEA ==========\n")

ssgsea_improved <- function(expr_matrix, gene_set, alpha = 0.25) {
  
  n_genes <- nrow(expr_matrix)
  n_samples <- ncol(expr_matrix)
  n_gs <- length(gene_set)
  
  scores <- numeric(n_samples)
  names(scores) <- colnames(expr_matrix)
  
  cat(sprintf("基因集大小: %d, 总基因数: %d\n", n_gs, n_genes))
  
  for (s in 1:n_samples) {
    expr <- expr_matrix[, s]
    expr[is.na(expr)] <- 0
    
    # 按表达值降序排序
    order_idx <- order(expr, decreasing = TRUE)
    ranked_genes <- rownames(expr_matrix)[order_idx]
    ranked_expr <- expr[order_idx]
    
    # 基因集中基因的位置
    gs_positions <- which(ranked_genes %in% gene_set)
    
    if (length(gs_positions) == 0) {
      scores[s] <- 0
      next
    }
    
    # 计算加权富集分数
    running_sum <- 0
    max_es <- 0
    min_es <- 0
    
    # 计算基因集中基因的权重（基于表达值）
    gs_weights <- abs(ranked_expr[gs_positions]) ^ alpha
    sum_gs_weights <- sum(gs_weights)
    
    # 如果所有权重为0，使用等权重
    if (sum_gs_weights == 0) {
      gs_weights <- rep(1, length(gs_positions))
      sum_gs_weights <- length(gs_positions)
    }
    
    # 计算非基因集基因的步长
    step_not_gs <- 1 / (n_genes - n_gs)
    
    for (i in 1:n_genes) {
      if (i %in% gs_positions) {
        # 基因集中的基因：加上加权贡献
        idx <- which(gs_positions == i)
        running_sum <- running_sum + (gs_weights[idx] / sum_gs_weights)
      } else {
        # 非基因集中的基因：减去步长
        running_sum <- running_sum - step_not_gs
      }
      
      if (running_sum > max_es) max_es <- running_sum
      if (running_sum < min_es) min_es <- running_sum
    }
    
    # 取最大绝对偏差
    scores[s] <- ifelse(abs(max_es) >= abs(min_es), max_es, min_es)
    
    if (s %% 50 == 0) cat(sprintf("  处理 %d/%d\n", s, n_samples))
  }
  
  return(scores)
}

# 尝试不同的alpha值
cat("\n尝试alpha=0.25...\n")
scores_025 <- ssgsea_improved(expr_matrix, available_genes, alpha = 0.25)

cat("\n尝试alpha=0.5...\n")
scores_05 <- ssgsea_improved(expr_matrix, available_genes, alpha = 0.5)

cat("\n尝试alpha=1.0...\n")
scores_10 <- ssgsea_improved(expr_matrix, available_genes, alpha = 1.0)

# ====================== 结果诊断 ======================
cat("\n========== 结果诊断 ==========\n")

diagnose_scores <- function(scores, name) {
  cat(sprintf("\n--- %s ---\n", name))
  cat("范围:", range(scores), "\n")
  cat("均值:", mean(scores), "\n")
  cat("中位数:", median(scores), "\n")
  cat("标准差:", sd(scores), "\n")
  cat("唯一值数量:", length(unique(scores)), "\n")
  
  # 检查双峰分布
  hist_data <- hist(scores, breaks = 30, plot = FALSE)
  peaks <- which(diff(sign(diff(hist_data$counts))) == -2) + 1
  cat("检测到的峰值数量:", length(peaks), "\n")
  
  return(list(scores = scores, peaks = length(peaks)))
}

diag_025 <- diagnose_scores(scores_025, "alpha=0.25")
diag_05 <- diagnose_scores(scores_05, "alpha=0.5")
diag_10 <- diagnose_scores(scores_10, "alpha=1.0")

# 选择最佳的alpha值（基于标准差最大，即区分度最好）
stds <- c(sd(scores_025), sd(scores_05), sd(scores_10))
best_alpha <- c(0.25, 0.5, 1.0)[which.max(stds)]
scores_final <- list(scores_025, scores_05, scores_10)[[which.max(stds)]]

cat(sprintf("\n选择alpha=%.2f (标准差最大: %.4f)\n", best_alpha, max(stds)))

# ====================== 结果整理 ======================
cat("\n========== 结果整理 ==========\n")

# 创建结果数据框
scores_df <- data.frame(
  Sample = colnames(expr_matrix),
  ssGSEA_Score = scores_final,
  stringsAsFactors = FALSE
)

# 添加分组信息
if (has_group) {
  sample_col <- toupper(group_sample_col)
  group_col <- toupper(group_group_col)
  
  scores_df$Sample_upper <- toupper(gsub("[[:space:]]", "", scores_df$Sample))
  group_info$Sample_upper <- toupper(gsub("[[:space:]]", "", group_info[[sample_col]]))
  
  scores_df <- scores_df %>%
    left_join(group_info[, c("Sample_upper", group_col)], by = "Sample_upper") %>%
    select(-Sample_upper)
  
  colnames(scores_df)[ncol(scores_df)] <- "Group"
}

# 按评分排序
scores_df <- scores_df[order(scores_df$ssGSEA_Score, decreasing = TRUE), ]

# ====================== 可视化 ======================
cat("\n========== 绘制诊断图 ==========\n")

pdf(output_plot, width = 12, height = 10)

# 设置多图布局
par(mfrow = c(2, 2))

# 1. 不同alpha值的分布比较
hist(scores_025, breaks = 30, main = "alpha=0.25", xlab = "ssGSEA Score", col = "lightblue")
hist(scores_05, breaks = 30, main = "alpha=0.5", xlab = "ssGSEA Score", col = "lightgreen")
hist(scores_10, breaks = 30, main = "alpha=1.0", xlab = "ssGSEA Score", col = "lightcoral")

# 2. 最终选择的alpha值的分布
hist(scores_final, breaks = 30, main = paste("Selected alpha =", best_alpha), 
     xlab = "ssGSEA Score", col = "steelblue")

dev.off()

# 使用ggplot2绘制更详细的图
pdf("ssGSEA_detailed.pdf", width = 10, height = 8)

p <- ggplot(scores_df, aes(x = ssGSEA_Score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, 
                 fill = "steelblue", color = "black", alpha = 0.7) +
  geom_density(color = "darkred", size = 1) +
  labs(title = paste("ssGSEA Score Distribution (alpha =", best_alpha, ")"),
       subtitle = paste("Samples:", nrow(scores_df), "| Unique values:", 
                        length(unique(scores_final))),
       x = "ssGSEA Score", y = "Density") +
  theme_minimal()

# 如果有分组信息，添加分组颜色
if (has_group && "Group" %in% colnames(scores_df)) {
  p <- ggplot(scores_df, aes(x = ssGSEA_Score, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, 
                   alpha = 0.7, position = "identity") +
    geom_density(aes(color = Group), size = 1, alpha = 0.5) +
    labs(title = paste("ssGSEA Score Distribution by Group (alpha =", best_alpha, ")"),
         x = "ssGSEA Score", y = "Density") +
    theme_minimal()
}

print(p)
dev.off()

cat("诊断图已保存\n")

# ====================== 保存结果 ======================
cat("\n========== 保存结果 ==========\n")

write.csv(scores_df, file = output_file, row.names = FALSE)
cat("结果已保存到:", output_file, "\n")

# ====================== 最终评估 ======================
cat("\n========== 最终评估 ==========\n")

unique_count <- length(unique(scores_final))
peak_count <- diag_10$peaks  # 使用最后一个诊断的峰值数

if (peak_count >= 2) {
  cat("⚠️ 检测到多峰分布\n")
  if (has_group) {
    cat("可能原因: 样品存在真实的生物学分组\n")
    cat("建议: 按分组分别分析或使用分组信息进行后续分析\n")
  } else {
    cat("可能原因: 算法偏差或数据质量问题\n")
    cat("建议: 尝试其他基因集富集方法如GSVA\n")
  }
} else if (unique_count < 10) {
  cat("⚠️ 评分区分度太低\n")
  cat("建议: 考虑扩大基因集或使用其他方法\n")
} else {
  cat("✓ 评分分布正常\n")
}

cat("\n==================== 完成 ====================\n")