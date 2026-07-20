# =============================================================================
# 置换检验分析：验证抑制剂对CIPS基因集的特异性抑制
# 功能：比较CIPS基因集的logFC是否显著大于随机基因集
# 输入：CIPS基因列表 + 所有基因的logFC
# 输出：置换检验直方图、统计结果
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/3.6-Rh2-Permutation test"  # 请修改此路径
setwd(work_dir)

# 2. 创建输出文件夹
output_dir <- "Permutation_test_results-1"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 3. 输入文件设置
cips_gene_file <- "CIPS_gene_list.txt"        # CIPS基因列表（每行一个基因名）
logfc_file <- "all_genes_logFC.csv"           # 所有基因的logFC（Gene, logFC两列）

# 4. 置换检验参数
n_permutations <- 10000                        # 置换次数
test_statistic <- "median"                     # 统计量类型: "median", "mean"
random_seed <- 1234                            # 随机种子

# 5. 图片设置（可自定义）
plot_width <- 8                                # 图片宽度（英寸）
plot_height <- 6                               # 图片高度（英寸）

# 颜色设置
hist_color <- "steelblue"                      # 直方图颜色
hist_alpha <- 0.7                              # 直方图透明度
line_color <- "red"                            # 红竖线颜色
line_lwd <- 1.5                                # 红竖线宽度
density_color <- "darkblue"                    # 密度曲线颜色
density_lwd <- 1                               # 密度曲线宽度

# X轴范围设置
x_axis_min <- -0.7                             # X轴最小值
x_axis_max <- 0.1                              # X轴最大值
x_axis_breaks <- seq(-0.7, 0.1, by = 0.1)     # X轴刻度间隔

# 字体大小设置
title_size <- 24                               # 主标题字体大小
subtitle_size <- 10                            # 副标题字体大小
axis_title_size <- 24                          # 坐标轴标题字体大小
axis_text_size <- 24                           # 坐标轴刻度字体大小

# 6. 输出文件设置
output_histogram <- file.path(output_dir, "Permutation_test_histogram.pdf")
output_results <- file.path(output_dir, "Permutation_test_results.csv")
output_summary <- file.path(output_dir, "analysis_summary.txt")

# ====================== 加载必要包 ======================
cat("\n========== 加载必要包 ==========\n")

packages <- c("dplyr", "ggplot2", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取数据 ======================
cat("\n========== 读取数据 ==========\n")

# 读取CIPS基因列表
if (file.exists(cips_gene_file)) {
  cips_genes <- readLines(cips_gene_file, warn = FALSE)
  cips_genes <- toupper(trimws(cips_genes))
  cips_genes <- cips_genes[cips_genes != ""]
  cips_genes <- cips_genes[!grepl("^#", cips_genes)]
  cat(sprintf("CIPS基因数: %d\n", length(cips_genes)))
} else {
  stop("错误：未找到CIPS基因列表文件！")
}

# 读取所有基因的logFC
if (file.exists(logfc_file)) {
  logfc_data <- read.csv(logfc_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # 检查列名
  if (!"Gene" %in% colnames(logfc_data) || !"logFC" %in% colnames(logfc_data)) {
    gene_col <- grep("gene|Gene|GENE", colnames(logfc_data), ignore.case = TRUE)
    logfc_col <- grep("logFC|log_fc|fc|log2FoldChange", colnames(logfc_data), ignore.case = TRUE)
    
    if (length(gene_col) > 0 && length(logfc_col) > 0) {
      colnames(logfc_data)[gene_col[1]] <- "Gene"
      colnames(logfc_data)[logfc_col[1]] <- "logFC"
    } else {
      stop("错误：无法识别基因名和logFC列！")
    }
  }
  
  # 基因名转大写
  logfc_data$Gene <- toupper(trimws(logfc_data$Gene))
  
  # 确保logFC列为数值型
  if (!is.numeric(logfc_data$logFC)) {
    cat("警告：logFC列不是数值型，正在转换...\n")
    logfc_data$logFC <- as.character(logfc_data$logFC)
    logfc_data$logFC[logfc_data$logFC %in% c("", "NA", "NULL")] <- NA
    logfc_data$logFC <- as.numeric(logfc_data$logFC)
  }
  
  # 移除NA和无限值
  logfc_data <- logfc_data[!is.na(logfc_data$logFC), ]
  logfc_data <- logfc_data[is.finite(logfc_data$logFC), ]
  
  cat(sprintf("logFC数据: %d 个基因\n", nrow(logfc_data)))
  
} else {
  stop("错误：未找到logFC文件！")
}

# ====================== 数据验证 ======================
cat("\n========== 数据验证 ==========\n")

cips_in_data <- intersect(cips_genes, logfc_data$Gene)
missing_genes <- setdiff(cips_genes, cips_in_data)

cat(sprintf("CIPS基因在logFC数据中的数量: %d/%d\n", length(cips_in_data), length(cips_genes)))
if (length(missing_genes) > 0) {
  cat(sprintf("警告：%d个CIPS基因不在logFC文件中\n", length(missing_genes)))
}

cips_genes_use <- cips_in_data
if (length(cips_genes_use) == 0) {
  stop("错误：没有找到任何CIPS基因！")
}

# 提取CIPS基因的logFC
cips_logfc <- logfc_data$logFC[match(cips_genes_use, logfc_data$Gene)]
cips_logfc <- as.numeric(cips_logfc)
cips_logfc <- cips_logfc[is.finite(cips_logfc)]

background_logfc <- logfc_data$logFC
background_logfc <- as.numeric(background_logfc)
background_logfc <- background_logfc[is.finite(background_logfc)]

cat(sprintf("\nCIPS基因logFC统计:\n"))
cat(sprintf("  有效基因数: %d\n", length(cips_logfc)))
cat(sprintf("  中位数: %.4f\n", median(cips_logfc)))
cat(sprintf("  均值: %.4f\n", mean(cips_logfc)))

# ====================== 置换检验 ======================
cat("\n========== 执行置换检验 ==========\n")

set.seed(random_seed)

if (test_statistic == "median") {
  actual_stat <- median(cips_logfc)
  stat_name <- "Median"
} else {
  actual_stat <- mean(cips_logfc)
  stat_name <- "Mean"
}

cat(sprintf("CIPS实际%s: %.6f\n", stat_name, actual_stat))

permuted_stats <- numeric(n_permutations)
n_cips_genes <- length(cips_genes_use)

for (i in 1:n_permutations) {
  random_indices <- sample(length(background_logfc), n_cips_genes, replace = FALSE)
  random_logfc <- background_logfc[random_indices]
  
  if (test_statistic == "median") {
    permuted_stats[i] <- median(random_logfc)
  } else {
    permuted_stats[i] <- mean(random_logfc)
  }
  
  if (i %% 1000 == 0) {
    cat(sprintf("  已完成: %d/%d\n", i, n_permutations))
  }
}

# 计算P值和Z-score
p_value <- sum(permuted_stats <= actual_stat) / n_permutations
z_score <- (actual_stat - mean(permuted_stats)) / sd(permuted_stats)

cat(sprintf("P值: %.6f\n", p_value))
cat(sprintf("Z-score: %.4f\n", z_score))

# 格式化P值显示
if (p_value < 0.0001) {
  p_display <- "P < 0.0001"
} else if (p_value < 0.001) {
  p_display <- "P < 0.001"
} else if (p_value < 0.01) {
  p_display <- sprintf("P = %.3f", p_value)
} else {
  p_display <- sprintf("P = %.4f", p_value)
}

# ====================== 绘制直方图 ======================
cat("\n========== 绘制直方图 ==========\n")

plot_df <- data.frame(Permuted_Stat = permuted_stats)

# 创建标题
if (actual_stat < 0) {
  direction <- "downregulation"
} else {
  direction <- "upregulation"
}

main_title <- sprintf("Permutation Test: Specificity of CIPS %s", direction)

# 副标题包含所有统计信息
subtitle <- sprintf("%s log2FC = %.4f | %s | Z-score = %.2f | Based on %d permutations | CIPS size = %d", 
                    stat_name, actual_stat, p_display, z_score, n_permutations, n_cips_genes)

# 创建图形
p <- ggplot(plot_df, aes(x = Permuted_Stat)) +
  # 直方图
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 50, 
                 fill = hist_color, 
                 alpha = hist_alpha,
                 color = "black",
                 linewidth = 0.3) +
  # 密度曲线
  geom_density(color = density_color, 
               linewidth = density_lwd, 
               alpha = 0.5) +
  # 红竖线（只保留线，不添加标注）
  geom_vline(xintercept = actual_stat, 
             color = line_color, 
             linewidth = line_lwd, 
             linetype = "solid") +
  # 设置X轴范围
  coord_cartesian(xlim = c(x_axis_min, x_axis_max)) +
  # 设置X轴刻度和标签（强制使用小数格式）
  scale_x_continuous(
    breaks = x_axis_breaks,
    labels = scales::label_number(accuracy = 0.1)
  ) +
  # 标签和标题
  labs(title = main_title,
       subtitle = subtitle,
       x = paste(stat_name, "log2 Fold Change of random gene sets"),
       y = "Density") +
  # 主题设置
  theme_minimal() +
  theme(plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = subtitle_size, hjust = 0.5),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5))

# 保存图片
ggsave(output_histogram, plot = p, width = plot_width, height = plot_height, dpi = 300)
cat(sprintf("直方图已保存: %s\n", output_histogram))

# ====================== 保存结果 ======================
cat("\n========== 保存结果 ==========\n")

results_df <- data.frame(
  Metric = c(
    "CIPS_gene_count",
    "Background_gene_count",
    "Permutation_iterations",
    "CIPS_actual_statistic",
    "Permutation_mean",
    "Permutation_sd",
    "Permutation_median",
    "Permutation_2.5_percentile",
    "Permutation_97.5_percentile",
    "P_value_one_tailed",
    "Z_score",
    "Random_seed"
  ),
  Value = c(
    n_cips_genes,
    length(background_logfc),
    n_permutations,
    sprintf("%.6f", actual_stat),
    sprintf("%.6f", mean(permuted_stats)),
    sprintf("%.6f", sd(permuted_stats)),
    sprintf("%.6f", median(permuted_stats)),
    sprintf("%.6f", quantile(permuted_stats, 0.025)),
    sprintf("%.6f", quantile(permuted_stats, 0.975)),
    sprintf("%.6f", p_value),
    sprintf("%.4f", z_score),
    random_seed
  ),
  stringsAsFactors = FALSE
)

write.csv(results_df, file = output_results, row.names = FALSE)
cat(sprintf("结果已保存: %s\n", output_results))

# ====================== 生成分析报告 ======================
cat("\n========== 生成分析报告 ==========\n")

sink(output_summary)

cat("========================================\n")
cat("   置换检验分析报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、输入数据概况\n")
cat("----------------\n")
cat(sprintf("背景基因总数: %d\n", length(background_logfc)))
cat(sprintf("CIPS基因数: %d\n", length(cips_genes_use)))

cat("\n二、CIPS基因logFC分布\n")
cat("---------------------\n")
cat(sprintf("  中位数: %.4f\n", median(cips_logfc)))
cat(sprintf("  均值: %.4f\n", mean(cips_logfc)))
cat(sprintf("  最小值: %.4f\n", min(cips_logfc)))
cat(sprintf("  最大值: %.4f\n", max(cips_logfc)))

cat("\n三、置换检验结果\n")
cat("----------------\n")
cat(sprintf("置换次数: %d\n", n_permutations))
cat(sprintf("CIPS实际%s: %.6f\n", stat_name, actual_stat))
cat(sprintf("随机分布均值: %.6f\n", mean(permuted_stats)))
cat(sprintf("随机分布标准差: %.6f\n", sd(permuted_stats)))
cat(sprintf("P值: %.6f\n", p_value))
cat(sprintf("Z-score: %.4f\n", z_score))

if (p_value < 0.001) {
  cat("\n结论: *** 非常显著 ***\n")
}

cat("\n========================================\n")
cat("                完成\n")
cat("========================================\n")

sink()

cat("\n分析报告已保存: ", output_summary, "\n")
cat("\n=================== 分析完成 ====================\n")