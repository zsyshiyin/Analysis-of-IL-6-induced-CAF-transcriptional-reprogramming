# ============================================================================
# 基因集表达评分分析（GSVA/ssGSEA）
# 功能：根据提供的基因列表，计算每个样品的表达评分，并用箱线图展示各组差异
# ============================================================================

# 清空环境变量
rm(list = ls())

# ========================== 工作路径设定 ====================================
work_dir <- "D:/zsy/SX/Fomal-final/27-CIPS-Tumor-heatmap"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")

# ========================== 创建输出目录 ====================================
output_dir <- "gene_set_scoring"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# ============================================================================
# ====================== 用户自定义参数区域 ==================================
# ============================================================================

expression_file <- "Tumor-mRNA.csv"
gene_list_file <- "gene_pathway_stats-xin.csv"
gene_id_col <- "GeneID"
list_gene_col <- "GeneID"

group_mapping <- list(
  "DMEM" = "DMEM", 
  "Con" = "Con",
  "I6" = "I6",
  "R2" = "R2"
)

group_order <- c("DMEM", "Con", "I6", "R2")

group_colors <- c(
  "DMEM" = "#2AA7DE", 
  "Con" = "#25377F",
  "I6" = "#C60036",
  "R2" = "#E4945A"
)

scoring_method <- "mean"
normalize_score <- TRUE
log2_transform <- TRUE
log_pseudocount <- 1

# ----------------------------- 箱线图参数 -----------------------------------
boxplot_width <- 8
boxplot_height <- 5.5
boxplot_title <- "Gene Set Expression Score"
boxplot_title_size <- 24
boxplot_x_label <- "Group"
boxplot_y_label <- "Expression Score"
boxplot_axis_title_size <- 22
boxplot_axis_text_size <- 20
boxplot_show_points <- TRUE
boxplot_point_size <- 3
boxplot_point_alpha <- 0.6
boxplot_jitter_width <- 0.3
boxplot_box_width <- 0.5
boxplot_box_alpha <- 0.7

# ----------------------------- Y轴范围参数 -----------------------------------
# 设置为NULL表示自动，设置数值如c(-2, 3)表示固定范围
y_axis_limits <- c(-0.5, 3.5)   # 自定义Y轴范围，为p值标注留出空间

# ----------------------------- p值标注参数 -----------------------------------
show_pvalue_labels <- TRUE           # 是否显示p值标注
pval_label_size <- 6                 # p值标签字体大小
pval_label_color <- "black"          # p值标签颜色
pval_line_color <- "black"           # 标注线颜色
pval_line_size <- 0.5                # 标注线粗细
pval_bracket_height <- 0.05          # 括号高度（相对于数据范围的比例）
pval_text_y_offset <- 0.15           # 文字相对于括号的偏移量
pval_step_increase <- 0.30            # 多个比较时的高度递增

# ----------------------------- 统计检验参数 ---------------------------------
stat_test_method <- "t.test"    # "t.test" 或 "wilcox.test"

# ----------------------------- 输出文件 -------------------------------------
output_score_file <- file.path(output_dir, "gene_set_scores.csv")
output_boxplot <- file.path(output_dir, "gene_set_score_boxplot.pdf")
output_summary <- file.path(output_dir, "analysis_summary.txt")

# ============================================================================
# ====================== 加载必要的包 ========================================
# ============================================================================

cat("\n========== 加载必要的包 ==========\n")

library(dplyr)
library(ggplot2)
library(reshape2)

cat("✓ 所有包已加载\n")

# ========================== 读取表达矩阵 ====================================
cat("\n========== 读取表达矩阵 ==========\n")

if (!file.exists(expression_file)) {
  stop("错误：未找到表达矩阵文件")
}

expr_data <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("表达矩阵维度: %d 行, %d 列\n", nrow(expr_data), ncol(expr_data)))

# 设置行名
if (gene_id_col %in% colnames(expr_data)) {
  rownames(expr_data) <- expr_data[[gene_id_col]]
  expr_data <- expr_data[, !colnames(expr_data) %in% gene_id_col]
} else {
  rownames(expr_data) <- expr_data[, 1]
  expr_data <- expr_data[, -1]
}

# 转换为数值矩阵
expression_matrix <- as.matrix(expr_data)
mode(expression_matrix) <- "numeric"
expression_matrix[is.na(expression_matrix)] <- 0

cat(sprintf("表达矩阵维度: %d 基因 × %d 样品\n", nrow(expression_matrix), ncol(expression_matrix)))

# ========================== 读取基因列表 ====================================
cat("\n========== 读取基因列表 ==========\n")

if (!file.exists(gene_list_file)) {
  stop("错误：未找到基因列表文件")
}

gene_data <- read.csv(gene_list_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("基因列表维度: %d 行, %d 列\n", nrow(gene_data), ncol(gene_data)))

# 提取基因列表
if (list_gene_col %in% colnames(gene_data)) {
  gene_set <- gene_data[[list_gene_col]]
} else {
  gene_set <- gene_data[, 1]
}

# 标准化基因名
gene_set <- toupper(trimws(gene_set))
rownames(expression_matrix) <- toupper(rownames(expression_matrix))

# 筛选匹配基因
available_genes <- intersect(gene_set, rownames(expression_matrix))
cat(sprintf("基因集大小: %d\n", length(gene_set)))
cat(sprintf("匹配的基因数: %d\n", length(available_genes)))

if (length(available_genes) == 0) {
  stop("错误：没有匹配的基因！")
}

# ========================== 数据预处理 ====================================
cat("\n========== 数据预处理 ==========\n")

# 提取匹配基因的表达数据
expr_subset <- expression_matrix[available_genes, , drop = FALSE]

# Log2转换
if (log2_transform) {
  cat("执行log2转换...\n")
  if (min(expr_subset, na.rm = TRUE) >= 0) {
    expr_subset <- log2(expr_subset + log_pseudocount)
    cat("  ✓ log2转换完成\n")
  }
}

# ========================== 计算表达评分 ====================================
cat("\n========== 计算表达评分 ==========\n")

if (scoring_method == "mean") {
  scores <- colMeans(expr_subset, na.rm = TRUE)
} else if (scoring_method == "median") {
  scores <- apply(expr_subset, 2, median, na.rm = TRUE)
}

# 标准化评分
if (normalize_score && sd(scores) > 0) {
  scores <- scale(scores)
  scores <- as.numeric(scores)
}

cat(sprintf("评分范围: [%.4f, %.4f]\n", min(scores), max(scores)))

# ========================== 创建评分数据框 ====================================
score_df <- data.frame(
  Sample = colnames(expr_subset),
  Score = as.numeric(scores),
  stringsAsFactors = FALSE
)

# 添加分组信息
score_df$Group <- NA
for (prefix in names(group_mapping)) {
  matching_idx <- grep(paste0("^", prefix), score_df$Sample)
  if (length(matching_idx) > 0) {
    score_df$Group[matching_idx] <- group_mapping[[prefix]]
  }
}

# 移除未识别的样品
score_df <- score_df[!is.na(score_df$Group), ]
score_df$Group <- factor(score_df$Group, levels = group_order)

cat("分组信息:\n")
print(table(score_df$Group))

# 保存评分结果
write.csv(score_df, output_score_file, row.names = FALSE)
cat(sprintf("✓ 评分结果已保存: %s\n", output_score_file))

# ========================== 评分统计 ====================================
score_stats <- score_df %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
print(score_stats)

# ========================== 计算p值 ====================================
cat("\n========== 计算统计检验p值 ==========\n")

# 计算两两比较的p值
pairwise_pvalues <- data.frame()

groups <- levels(score_df$Group)
for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    group1_data <- score_df$Score[score_df$Group == groups[i]]
    group2_data <- score_df$Score[score_df$Group == groups[j]]
    
    if (length(group1_data) >= 2 && length(group2_data) >= 2) {
      if (stat_test_method == "t.test") {
        test_result <- t.test(group1_data, group2_data)
      } else {
        test_result <- wilcox.test(group1_data, group2_data)
      }
      
      # 格式化p值显示（星号标记）
      p_val <- test_result$p.value
      if (p_val < 0.0001) {
        p_label <- "****"
      } else if (p_val < 0.001) {
        p_label <- "***"
      } else if (p_val < 0.01) {
        p_label <- "**"
      } else if (p_val < 0.05) {
        p_label <- "*"
      } else {
        p_label <- "ns"
      }
      
      pairwise_pvalues <- rbind(pairwise_pvalues, data.frame(
        Group1 = groups[i],
        Group2 = groups[j],
        P_value = test_result$p.value,
        P_label = p_label,
        Significant = test_result$p.value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n两两比较p值:\n")
print(pairwise_pvalues)

# ========================== 计算Y轴范围 ====================================
if (is.null(y_axis_limits)) {
  # 自动计算范围，为p值标注留出空间
  y_max <- max(score_df$Score, na.rm = TRUE)
  y_min <- min(score_df$Score, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # 根据比较数量增加顶部空间
  n_comparisons <- nrow(pairwise_pvalues)
  y_max_expanded <- y_max + y_range * (0.2 + n_comparisons * pval_step_increase)
  y_axis_limits <- c(y_min, y_max_expanded)
} else {
  # 使用用户自定义的范围
  cat(sprintf("使用自定义Y轴范围: [%.2f, %.2f]\n", y_axis_limits[1], y_axis_limits[2]))
}

# ========================== 绘制箱线图 ====================================
cat("\n========== 绘制箱线图 ==========\n")

# 创建箱线图
p <- ggplot(score_df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(width = boxplot_box_width, alpha = boxplot_box_alpha) +
  scale_fill_manual(values = group_colors, name = "Group") +
  labs(title = boxplot_title, x = boxplot_x_label, y = boxplot_y_label) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = boxplot_title_size, face = "bold"),
        axis.title = element_text(size = boxplot_axis_title_size),
        axis.text = element_text(size = boxplot_axis_text_size)) +
  coord_cartesian(ylim = y_axis_limits)

# 添加散点
if (boxplot_show_points) {
  p <- p + geom_jitter(width = boxplot_jitter_width, size = boxplot_point_size, 
                       alpha = boxplot_point_alpha, aes(color = Group)) +
    scale_color_manual(values = group_colors, guide = "none")
}

# ========================== 添加p值标注 ====================================
if (show_pvalue_labels && nrow(pairwise_pvalues) > 0) {
  
  # 获取各组在x轴上的位置
  x_positions <- 1:length(group_order)
  names(x_positions) <- group_order
  
  # 计算Y轴基础位置
  y_max_data <- max(score_df$Score, na.rm = TRUE)
  y_range_data <- diff(range(score_df$Score, na.rm = TRUE))
  
  # 为每个比较分配Y位置
  for (i in 1:nrow(pairwise_pvalues)) {
    comp <- pairwise_pvalues[i, ]
    
    # 获取x坐标
    x1 <- x_positions[comp$Group1]
    x2 <- x_positions[comp$Group2]
    
    # 计算Y位置（根据比较顺序递增）
    y_position <- y_max_data + y_range_data * (pval_bracket_height + i * pval_step_increase)
    
    # 添加括号线
    p <- p +
      # 左侧竖线
      annotate("segment", x = x1, xend = x1, 
               y = y_position - y_range_data * pval_bracket_height * 0.5, 
               yend = y_position,
               color = pval_line_color, size = pval_line_size) +
      # 右侧竖线
      annotate("segment", x = x2, xend = x2, 
               y = y_position - y_range_data * pval_bracket_height * 0.5, 
               yend = y_position,
               color = pval_line_color, size = pval_line_size) +
      # 横线
      annotate("segment", x = x1, xend = x2, 
               y = y_position, yend = y_position,
               color = pval_line_color, size = pval_line_size) +
      # p值标签
      annotate("text", x = (x1 + x2) / 2, 
               y = y_position + y_range_data * pval_text_y_offset,
               label = comp$P_label, 
               size = pval_label_size, color = pval_label_color)
  }
}

# ========================== 保存图片 ====================================
ggsave(output_boxplot, p, width = boxplot_width, height = boxplot_height, dpi = 300)
cat(sprintf("✓ 箱线图已保存: %s\n", output_boxplot))

# ========================== 生成报告 ====================================
cat("\n========== 生成分析报告 ==========\n")

sink(output_summary)

cat("================================================================================\n")
cat("                    基因集表达评分分析报告\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、输入数据概况\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("表达矩阵基因数: %d\n", nrow(expression_matrix)))
cat(sprintf("表达矩阵样品数: %d\n", ncol(expression_matrix)))
cat(sprintf("基因集大小: %d\n", length(gene_set)))
cat(sprintf("匹配的基因数: %d\n", length(available_genes)))

cat("\n二、评分方法\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("评分方法: %s\n", scoring_method))
cat(sprintf("Log2转换: %s\n", ifelse(log2_transform, "是", "否")))
cat(sprintf("Z-score标准化: %s\n", ifelse(normalize_score, "是", "否")))

cat("\n三、分组信息\n")
cat("--------------------------------------------------------------------------------\n")
print(table(score_df$Group))

cat("\n四、评分统计\n")
cat("--------------------------------------------------------------------------------\n")
print(score_stats)

cat("\n五、统计检验结果\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("检验方法: %s\n", stat_test_method))
for (i in 1:nrow(pairwise_pvalues)) {
  cat(sprintf("  %s vs %s: %s\n", 
              pairwise_pvalues$Group1[i], pairwise_pvalues$Group2[i],
              pairwise_pvalues$P_label[i]))
}

cat("\n六、图形参数\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("Y轴范围: [%.2f, %.2f]\n", y_axis_limits[1], y_axis_limits[2]))
cat(sprintf("图片尺寸: %.1f × %.1f 英寸\n", boxplot_width, boxplot_height))

cat("\n七、输出文件\n")
cat("--------------------------------------------------------------------------------\n")
cat("1.", basename(output_score_file), "- 表达评分结果\n")
cat("2.", basename(output_boxplot), "- 箱线图\n")

sink()

cat("\n========================================\n")
cat("分析完成！\n")
cat("输出目录:", output_dir, "\n")
cat("========================================\n")