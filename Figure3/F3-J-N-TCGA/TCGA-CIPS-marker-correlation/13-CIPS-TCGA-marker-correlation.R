# =============================================================================
# 生信分析脚本：基因表达与ssGSEA评分的相关性分析（Spearman相关）
# 功能：分析指定基因的表达水平与ssGSEA评分的相关性，绘制散点图
# 输入：mRNA.txt（TCGA转录组数据）
#       ssGSEA_scores_diagnostic.csv（ssGSEA评分文件）
# 输出：相关性散点图、相关性统计结果
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/13-CIPS-TCGA-CAF-marker-correlation"  # 请修改此路径
setwd(work_dir)

# 2. 创建输出文件夹
output_dir <- "gene_ssGSEA_correlation_analysis_ACTA2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 3. 输入文件名设置
expression_file <- "mRNA.txt"                    # TCGA转录组数据
score_file <- "ssGSEA_scores_diagnostic.csv"     # ssGSEA评分文件

# 4. 列名设置
score_id_col <- "Sample"                           # 样品ID列
score_value_col <- "ssGSEA_Score"                  # 评分列

# 5. 要分析的基因列表（可设置单个或多个基因）
# 方式1：直接指定基因名（推荐）
#genes_of_interest <- c("TP53", "EGFR", "KRAS", "BRAF", "PIK3CA", "PTEN", "CDH1", "VIM")

# 方式2：如果只想分析单个基因，使用：
genes_of_interest <- c("ACTA2")
# 方式3：如果要分析所有基因，设置为：genes_of_interest <- "all"

# 6. 相关性分析方法
cor_method <- "spearman"          # 相关性方法: "spearman", "pearson", "kendall"
p_adjust_method <- "none"           # p值校正方法: "BH", "bonferroni", "holm", "none"

# 7. 图片设置（完全可自定义）
plot_width <- 7                    # 单图宽度（英寸）
plot_height <- 8                   # 单图高度（英寸）
multi_plot_width <- 12             # 多图布局宽度（英寸）
multi_plot_height <- 12            # 多图布局高度（英寸）
point_pch <- 16                    # 点形状（16=实心圆）
point_cex <- 1.5                   # 点大小
point_color <- "#CA0E12"           # 点的颜色
line_color <- "black"              # 回归线颜色
line_lwd <- 4                      # 回归线宽度
title_cex <- 1.5                   # 标题字体大小
label_cex <- 1.8                   # 轴标签字体大小
axis_cex <- 1.5                    # 轴刻度字体大小
mgp_values <- c(2.5, 0.8, 0)      # 坐标轴标题位置

# 8. 输出文件设置
output_correlation_details <- file.path(output_dir, "gene_ssGSEA_correlation_details.csv")
output_summary <- file.path(output_dir, "analysis_summary.txt")
output_single_gene_dir <- file.path(output_dir, "single_gene_plots")
output_combined_plot <- file.path(output_dir, "combined_gene_correlation.pdf")

# 创建单基因图文件夹
if (!dir.exists(output_single_gene_dir)) {
  dir.create(output_single_gene_dir, recursive = TRUE)
  cat("创建单基因图文件夹:", output_single_gene_dir, "\n")
}

# ====================== 加载必要的包 ======================
cat("\n========== 加载必要的R包 ==========\n")

packages <- c("dplyr", "tidyr", "stringr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取表达矩阵 ======================
cat("\n========== 读取表达矩阵 ==========\n")

# 检测文件分隔符
first_line <- readLines(expression_file, n = 1)
if (grepl("\t", first_line)) {
  expr_data <- read.table(expression_file, 
                          header = TRUE, 
                          stringsAsFactors = FALSE, 
                          check.names = FALSE,
                          sep = "\t",
                          row.names = 1)
  cat("检测到制表符分隔符，第一列为基因名\n")
} else {
  expr_data <- read.csv(expression_file, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE,
                        row.names = 1)
  cat("检测到逗号分隔符，第一列为基因名\n")
}

cat(sprintf("表达矩阵维度: %d 个基因, %d 个样品\n", nrow(expr_data), ncol(expr_data)))

# ====================== 读取ssGSEA评分数据 ======================
cat("\n========== 读取ssGSEA评分数据 ==========\n")

score_data <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("评分数据维度: %d 行, %d 列\n", nrow(score_data), ncol(score_data)))

# ====================== ID标准化和匹配 ======================
cat("\n========== ID标准化和匹配 ==========\n")

expr_samples <- colnames(expr_data)
score_samples <- score_data[[score_id_col]]

# 样品名标准化（转为大写并去除空格）
expr_samples_clean <- toupper(gsub("[[:space:]]", "", expr_samples))
score_samples_clean <- toupper(gsub("[[:space:]]", "", score_samples))

# 尝试多种匹配方法
common_samples <- intersect(expr_samples_clean, score_samples_clean)

if (length(common_samples) == 0) {
  # 尝试截取TCGA ID前12位
  expr_samples_trim <- substr(expr_samples_clean, 1, 12)
  score_samples_trim <- substr(score_samples_clean, 1, 12)
  common_samples <- intersect(expr_samples_trim, score_samples_trim)
  cat("使用TCGA ID前12位进行匹配\n")
}

if (length(common_samples) == 0) {
  stop("错误：无法匹配样品ID！请检查样品ID格式。")
}

# 创建匹配的数据框
matched_data <- data.frame(
  Expression_Sample = expr_samples[match(common_samples, expr_samples_clean)],
  Score_Sample = score_samples[match(common_samples, score_samples_clean)],
  Match_ID = common_samples,
  stringsAsFactors = FALSE
)

cat(sprintf("\n成功匹配的样品数: %d\n", nrow(matched_data)))

# ====================== 对齐数据 ======================
cat("\n========== 对齐数据 ==========\n")

# 提取表达矩阵中匹配的样品
expr_aligned <- expr_data[, matched_data$Expression_Sample, drop = FALSE]

# 提取评分
score_aligned <- score_data[match(matched_data$Score_Sample, score_data[[score_id_col]]), score_value_col]
score_aligned <- as.numeric(score_aligned)

# 移除NA
valid_idx <- !is.na(score_aligned)
expr_aligned <- expr_aligned[, valid_idx, drop = FALSE]
score_aligned <- score_aligned[valid_idx]
matched_data <- matched_data[valid_idx, ]

cat(sprintf("最终分析样品数: %d\n", ncol(expr_aligned)))

# ====================== 处理要分析的基因 ======================
cat("\n========== 处理要分析的基因 ==========\n")

# 基因名标准化
rownames(expr_aligned) <- toupper(rownames(expr_aligned))
genes_of_interest <- toupper(genes_of_interest)

# 确定实际要分析的基因
if (length(genes_of_interest) == 1 && genes_of_interest[1] == "ALL") {
  # 分析所有基因
  genes_to_analyze <- rownames(expr_aligned)
  cat("分析模式：分析所有基因\n")
} else {
  # 分析指定基因
  genes_to_analyze <- intersect(genes_of_interest, rownames(expr_aligned))
  cat("分析模式：分析指定基因\n")
}

# 检查基因是否存在
missing_genes <- setdiff(genes_of_interest, genes_to_analyze)
if (length(missing_genes) > 0) {
  cat(sprintf("警告：以下 %d 个基因在表达矩阵中不存在: %s\n", 
              length(missing_genes), paste(missing_genes, collapse = ", ")))
}

if (length(genes_to_analyze) == 0) {
  stop("错误：没有找到任何要分析的基因！")
}

cat(sprintf("成功找到 %d 个基因进行分析\n", length(genes_to_analyze)))
cat("基因列表:\n")
print(genes_to_analyze)

# 提取要分析的基因的表达数据
expr_selected <- expr_aligned[genes_to_analyze, , drop = FALSE]

# ====================== 创建分析数据框 ======================
cat("\n========== 创建分析数据框 ==========\n")

# 转置表达数据
expr_t <- as.data.frame(t(expr_selected))
expr_t$ssGSEA_Score <- score_aligned

cat(sprintf("分析数据框维度: %d 行, %d 列\n", nrow(expr_t), ncol(expr_t)))

# ====================== 计算相关性 ======================
cat("\n========== 计算相关性 ==========\n")

# 存储相关性结果
cor_results <- data.frame(
  Gene = genes_to_analyze,
  Correlation = NA,
  P_value = NA,
  Lower_CI = NA,
  Upper_CI = NA,
  Mean_Expression = NA,
  SD_Expression = NA,
  stringsAsFactors = FALSE
)

# 计算每个基因的相关性
for (i in seq_along(genes_to_analyze)) {
  gene <- genes_to_analyze[i]
  
  # 提取数据
  gene_expr <- expr_t[[gene]]
  score <- expr_t$ssGSEA_Score
  
  # 移除缺失值
  valid_idx <- complete.cases(gene_expr, score)
  
  if (sum(valid_idx) >= 3) {
    # 计算相关性（对于Spearman，设置conf.int = TRUE以计算置信区间）
    # 注意：Spearman的置信区间是通过bootstrap计算的，可能较慢
    if (cor_method == "spearman") {
      test <- cor.test(gene_expr[valid_idx], score[valid_idx], 
                       method = cor_method, 
                       exact = FALSE,
                       conf.int = TRUE)
    } else {
      test <- cor.test(gene_expr[valid_idx], score[valid_idx], 
                       method = cor_method,
                       conf.int = TRUE)
    }
    
    cor_results$Correlation[i] <- test$estimate
    cor_results$P_value[i] <- test$p.value
    
    # 检查是否有置信区间
    if (!is.null(test$conf.int)) {
      cor_results$Lower_CI[i] <- test$conf.int[1]
      cor_results$Upper_CI[i] <- test$conf.int[2]
    } else {
      cor_results$Lower_CI[i] <- NA
      cor_results$Upper_CI[i] <- NA
    }
    
    cor_results$Mean_Expression[i] <- mean(gene_expr[valid_idx], na.rm = TRUE)
    cor_results$SD_Expression[i] <- sd(gene_expr[valid_idx], na.rm = TRUE)
  } else {
    cat(sprintf("警告：基因 %s 的有效数据点不足3个，跳过\n", gene))
  }
}

# 移除无效结果
cor_results <- cor_results[!is.na(cor_results$Correlation), ]

# 添加显著性标记
cor_results$Significant <- cor_results$P_value < 0.05
cor_results$Signif_mark <- ifelse(cor_results$P_value < 0.001, "***",
                                  ifelse(cor_results$P_value < 0.01, "**",
                                         ifelse(cor_results$P_value < 0.05, "*", "ns")))

# p值校正
if (p_adjust_method != "none" && nrow(cor_results) > 1) {
  cor_results$Adjusted_P_value <- p.adjust(cor_results$P_value, method = p_adjust_method)
  cor_results$Adjusted_Significant <- cor_results$Adjusted_P_value < 0.05
} else {
  cor_results$Adjusted_P_value <- cor_results$P_value
  cor_results$Adjusted_Significant <- cor_results$Significant
}

# 按相关性绝对值排序
cor_results <- cor_results[order(-abs(cor_results$Correlation)), ]

# ====================== 保存相关性结果 ======================
cat("\n========== 保存相关性结果 ==========\n")

write.csv(cor_results, file = output_correlation_details, row.names = FALSE)
cat("相关性详细结果已保存到:", output_correlation_details, "\n")

# 显示结果摘要
cat("\n相关性结果摘要:\n")
cat(sprintf("总基因数: %d\n", nrow(cor_results)))
cat(sprintf("显著相关基因数 (p<0.05): %d (%.1f%%)\n", 
            sum(cor_results$Significant, na.rm = TRUE),
            100 * sum(cor_results$Significant, na.rm = TRUE) / nrow(cor_results)))
cat(sprintf("校正后显著相关基因数: %d (%.1f%%)\n",
            sum(cor_results$Adjusted_Significant, na.rm = TRUE),
            100 * sum(cor_results$Adjusted_Significant, na.rm = TRUE) / nrow(cor_results)))

# 显示相关性最强的基因
cat("\n相关性最强的5个基因:\n")
top5 <- head(cor_results, 5)
for (i in 1:nrow(top5)) {
  if (top5$P_value[i] < 0.001) {
    p_disp <- "p<0.001"
  } else {
    p_disp <- sprintf("p=%.3f", top5$P_value[i])
  }
  cat(sprintf("  %s: r = %.3f%s, %s\n", 
              top5$Gene[i], top5$Correlation[i], 
              top5$Signif_mark[i], p_disp))
}

# ====================== 绘制单个基因散点图 ======================
cat("\n========== 绘制单个基因散点图 ==========\n")

# 绘图函数
plot_gene_correlation <- function(data, gene_name, cor_info) {
  # 获取该基因的相关性信息
  gene_cor <- cor_info[cor_info$Gene == gene_name, ]
  
  if (nrow(gene_cor) == 0 || is.na(gene_cor$Correlation)) {
    cat(sprintf("警告：无法找到基因 %s 的相关性信息\n", gene_name))
    return(NULL)
  }
  
  # 准备数据
  x <- data$ssGSEA_Score
  y <- data[[gene_name]]
  
  # 移除缺失值
  valid_idx <- complete.cases(x, y)
  x <- x[valid_idx]
  y <- y[valid_idx]
  
  # 拟合线性模型
  z <- lm(y ~ x)
  
  # 格式化p值
  if (gene_cor$P_value < 0.001) {
    p_display <- "p < 0.001"
  } else {
    p_display <- sprintf("p = %.3f", gene_cor$P_value)
  }
  
  # 创建标题（只显示原始p值，不显示校正p值）
  main_title <- sprintf("%s\nSpearman's r = %.3f%s\n%s", 
                        gene_name,
                        gene_cor$Correlation, 
                        gene_cor$Signif_mark,
                        p_display)
  
  # 绘制图形
  plot(x, y, 
       type = "p", 
       pch = point_pch, 
       cex = point_cex, 
       col = point_color,
       main = main_title,
       xlab = "ssGSEA Score", 
       ylab = paste0(gene_name, " Expression"),
       cex.main = title_cex * 0.8,
       cex.lab = label_cex * 0.8,
       cex.axis = axis_cex * 0.8,
       mgp = mgp_values)
  
  # 添加回归线
  lines(x[order(x)], fitted(z)[order(x)], col = line_color, lwd = line_lwd)
  
  # 添加置信区间
  pred <- predict(z, interval = "confidence")
  lines(x[order(x)], pred[order(x), 2], col = "gray", lty = 2, lwd = 1)
  lines(x[order(x)], pred[order(x), 3], col = "gray", lty = 2, lwd = 1)
}

# 绘制单个基因图（每个基因单独保存）
if (nrow(cor_results) > 0) {
  cat("\n正在绘制单基因散点图...\n")
  for (i in 1:nrow(cor_results)) {
    gene <- cor_results$Gene[i]
    
    # 生成输出文件名
    output_file <- file.path(output_single_gene_dir, 
                             sprintf("%s_correlation.pdf", gene))
    
    # 绘制图形
    pdf(output_file, width = plot_width, height = plot_height)
    plot_gene_correlation(expr_t, gene, cor_results)
    dev.off()
    
    cat(sprintf("  已保存: %s\n", basename(output_file)))
  }
  cat(sprintf("\n共保存 %d 个单基因散点图\n", nrow(cor_results)))
}

# ====================== 绘制组合散点图 ======================
cat("\n========== 绘制组合散点图 ==========\n")

if (nrow(cor_results) > 0) {
  # 计算需要的页面数（每页最多4个图）
  n_genes <- nrow(cor_results)
  genes_per_page <- 4
  n_pages <- ceiling(n_genes / genes_per_page)
  
  pdf(output_combined_plot, width = multi_plot_width, height = multi_plot_height)
  
  for (page in 1:n_pages) {
    # 设置2x2布局
    par(mfrow = c(2, 2))
    
    start_idx <- (page - 1) * genes_per_page + 1
    end_idx <- min(page * genes_per_page, n_genes)
    
    for (i in start_idx:end_idx) {
      gene <- cor_results$Gene[i]
      plot_gene_correlation(expr_t, gene, cor_results)
    }
    
    # 如果最后一页不足4个图，填充空白
    if (page == n_pages && (end_idx - start_idx + 1) < 4) {
      remaining <- 4 - (end_idx - start_idx + 1)
      for (j in 1:remaining) {
        plot.new()
      }
    }
  }
  
  dev.off()
  cat("组合散点图已保存到:", output_combined_plot, "\n")
}

# ====================== 生成分析报告 ======================
cat("\n========== 生成分析报告 ==========\n")

sink(output_summary)

cat("========================================\n")
cat("   基因表达与ssGSEA评分相关性分析报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、数据概况\n")
cat("------------\n")
cat(sprintf("表达矩阵基因数: %d\n", nrow(expr_data)))
cat(sprintf("表达矩阵样品数: %d\n", ncol(expr_data)))
cat(sprintf("ssGSEA评分样品数: %d\n", nrow(score_data)))
cat(sprintf("成功匹配样品数: %d\n", ncol(expr_aligned)))
cat(sprintf("要分析的基因数: %d\n", length(genes_of_interest)))
cat(sprintf("成功分析的基因数: %d\n", nrow(cor_results)))

cat("\n二、相关性分析方法\n")
cat("------------------\n")
cat(sprintf("相关性方法: %s\n", cor_method))
cat(sprintf("p值校正方法: %s\n", ifelse(p_adjust_method == "none", "未校正", p_adjust_method)))

cat("\n三、相关性结果摘要\n")
cat("------------------\n")
cat(sprintf("总基因数: %d\n", nrow(cor_results)))
cat(sprintf("显著相关基因数 (p<0.05): %d (%.1f%%)\n", 
            sum(cor_results$Significant, na.rm = TRUE),
            100 * sum(cor_results$Significant, na.rm = TRUE) / nrow(cor_results)))
cat(sprintf("校正后显著相关基因数: %d (%.1f%%)\n",
            sum(cor_results$Adjusted_Significant, na.rm = TRUE),
            100 * sum(cor_results$Adjusted_Significant, na.rm = TRUE) / nrow(cor_results)))

cat("\n四、正相关最强的5个基因\n")
cat("------------------------\n")
top5_positive <- head(cor_results[cor_results$Correlation > 0, ], 5)
if (nrow(top5_positive) > 0) {
  for (i in 1:nrow(top5_positive)) {
    if (top5_positive$P_value[i] < 0.001) {
      p_disp <- "p<0.001"
    } else {
      p_disp <- sprintf("p=%.3f", top5_positive$P_value[i])
    }
    cat(sprintf("  %s: r = %.3f%s, %s\n", 
                top5_positive$Gene[i], top5_positive$Correlation[i], 
                top5_positive$Signif_mark[i], p_disp))
  }
} else {
  cat("  无正相关基因\n")
}

cat("\n五、负相关最强的5个基因\n")
cat("------------------------\n")
top5_negative <- head(cor_results[cor_results$Correlation < 0, ], 5)
if (nrow(top5_negative) > 0) {
  for (i in 1:nrow(top5_negative)) {
    if (top5_negative$P_value[i] < 0.001) {
      p_disp <- "p<0.001"
    } else {
      p_disp <- sprintf("p=%.3f", top5_negative$P_value[i])
    }
    cat(sprintf("  %s: r = %.3f%s, %s\n", 
                top5_negative$Gene[i], top5_negative$Correlation[i], 
                top5_negative$Signif_mark[i], p_disp))
  }
} else {
  cat("  无负相关基因\n")
}

cat("\n六、输出文件列表\n")
cat("----------------\n")
cat("1.", basename(output_correlation_details), "- 相关性详细结果\n")
cat("2.", basename(output_combined_plot), "- 组合散点图\n")
cat("3.", basename(output_single_gene_dir), "/ - 单基因散点图文件夹\n")
cat("   包含", nrow(cor_results), "个单基因散点图\n")

cat("\n========================================\n")
cat("                完成\n")
cat("========================================\n")

sink()

cat("\n分析报告已保存到:", output_summary, "\n")
cat("\n=================== 分析完成 ===================\n")