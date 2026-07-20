# ============================================================
# ChIP-Atlas 富集分析结果可视化（自动读取TF列表版）
# 功能：自动读取工作目录下的TF列表文件，箱线图展示富集分数分布
# ============================================================

# -------------------- 1. 环境准备 --------------------
library(ggplot2)
library(dplyr)
library(readr)

# ============================================================
# ★★★ 用户自定义参数区域 ★★★
# ============================================================

# ---------- 路径设置 ----------
work_dir <- "D:/zsy/SX/GEO-ChIP-seq"
input_file <- file.path(work_dir, "ChIP_Atlas_enrichment_results.csv")
output_dir <- file.path(work_dir, "figures-xin")

# ---------- TF列表文件（放在工作目录下，自动读取）----------
# 支持的文件命名格式：
#   - key_tfs.csv / Key_TFs.csv / key_tf_list.csv
#   - de_tfs.csv / DE_TFs.csv / de_tf_list.csv
# 如果找不到文件，则使用下方的手动列表作为备选

# 备选手动列表（仅在自动读取失败时使用）
key_tfs_backup <- c("STAT3", "NFKB1", "RELA")
de_tfs_backup <- c("E2F1", "SNAI2", "SOX9", "SMAD3", "JUN", "FOXM1", "TCF4", "TP53")

# ---------- 筛选参数 ----------
fe_threshold <- 2  # 富集倍数阈值

# ---------- 图片尺寸 ----------
boxplot_width <- 4
boxplot_height <- 5

# ---------- 字体大小 ----------
boxplot_title_size <- 12
boxplot_axis_title_size <- 18
boxplot_axis_text_size <- 18
pval_label_size <- 10

# ---------- 箱线图参数 ----------
boxplot_title <- "CIPS Gene Enrichment Across Transcription Factor Targets"
boxplot_x_label <- ""
boxplot_y_label <- "Fold Enrichment"
boxplot_show_points <- TRUE
boxplot_point_size <- 2.5
boxplot_point_alpha <- 0.6
boxplot_jitter_width <- 0.2
boxplot_box_width <- 0.5
boxplot_box_alpha <- 0.7

# ---------- 配色方案（Key TFs按顺序分配颜色）----------
key_tf_colors <- c("#C60036", "#E4945A", "#00A087", "#4DBBD5", "#7E6148", "#6A51A3")
color_de_tfs <- "#25377F"

# ---------- Y轴范围参数 ----------
y_axis_limits <- NULL   # NULL表示自动计算

# ---------- p值标注参数 ----------
show_pvalue_labels <- TRUE
pval_label_color <- "black"
pval_line_color <- "black"
pval_line_size <- 0.5
pval_bracket_height <- 0.01
pval_text_y_offset <- 0.05
pval_step_increase <- 0.15
stat_test_method <- "t.test"  # "t.test" 或 "wilcox.test"

# ---------- 输出文件 ----------
output_boxplot <- file.path(output_dir, "ChIP_Atlas_Boxplot_Key_vs_DE.pdf")
output_stats <- file.path(output_dir, "ChIP_Atlas_Statistics_Key_vs_DE.txt")

# ============================================================
# 以下代码无需修改
# ============================================================

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

cat("========================================\n")
cat("ChIP-Atlas 富集分析 - 箱线图绘制\n")
cat("========================================\n")

# -------------------- 2. 自动读取TF列表文件 --------------------
cat("\n========== 自动读取TF列表文件 ==========\n")

# 定义可能的关键词文件名模式
key_tf_patterns <- c("key_tfs.csv", "Key_TFs.csv", "key_tf_list.csv", 
                     "Key_TF.csv", "key.csv", "*key*.csv")
de_tf_patterns <- c("de_tfs.csv", "DE_TFs.csv", "de_tf_list.csv", 
                    "DE_TF.csv", "de.csv", "*DE*.csv", "*de*.csv")

# 函数：查找并读取TF列表文件
find_and_read_tf_file <- function(patterns, default_list, list_name) {
  all_files <- list.files(work_dir, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  
  # 按模式匹配
  for (pattern in patterns) {
    matched <- grep(pattern, all_files, value = TRUE, ignore.case = TRUE)
    if (length(matched) > 0) {
      file_to_read <- matched[1]
      cat(sprintf("  - 找到%s文件: %s\n", list_name, file_to_read))
      
      # 读取文件
      tf_data <- read_csv(file.path(work_dir, file_to_read), show_col_types = FALSE)
      
      # 尝试识别TF列
      if ("TF" %in% colnames(tf_data)) {
        tf_list <- unique(tf_data$TF)
      } else if ("GeneSymbol" %in% colnames(tf_data)) {
        tf_list <- unique(tf_data$GeneSymbol)
      } else if ("Symbol" %in% colnames(tf_data)) {
        tf_list <- unique(tf_data$Symbol)
      } else {
        tf_list <- unique(tf_data[[1]])  # 使用第一列
      }
      
      tf_list <- tf_list[!is.na(tf_list) & tf_list != ""]
      cat(sprintf("    读取到 %d 个TF\n", length(tf_list)))
      return(tf_list)
    }
  }
  
  # 没找到文件，使用备选列表
  cat(sprintf("  - 未找到%s文件，使用备选手动列表 (%d 个TF)\n", list_name, length(default_list)))
  return(default_list)
}

# 读取Key TFs
key_tfs <- find_and_read_tf_file(key_tf_patterns, key_tfs_backup, "Key TFs")
cat(sprintf("\nKey TFs (%d个): %s\n", length(key_tfs), paste(key_tfs, collapse = ", ")))

# 读取DE-TFs
de_tfs <- find_and_read_tf_file(de_tf_patterns, de_tfs_backup, "DE-TFs")
cat(sprintf("DE-TFs (%d个): %s\n", length(de_tfs), paste(de_tfs, collapse = ", ")))

# -------------------- 3. 读取ChIP-Atlas数据 --------------------
cat("\n========== 读取ChIP-Atlas数据 ==========\n")

if (!file.exists(input_file)) {
  stop("输入文件不存在: ", input_file)
}

enrichment_data <- read.csv(input_file, stringsAsFactors = FALSE)
cat("  - 共读取", nrow(enrichment_data), "条记录\n")

# 检查必需的列
if (!"Feature" %in% colnames(enrichment_data)) {
  stop("ChIP-Atlas数据中缺少 'Feature' 列")
}
if (!"Fold_Enrichment" %in% colnames(enrichment_data)) {
  stop("ChIP-Atlas数据中缺少 'Fold_Enrichment' 列")
}

# 筛选
enrichment_data <- enrichment_data[enrichment_data$Fold_Enrichment > fe_threshold, ]
cat("  - 筛选后剩余", nrow(enrichment_data), "条记录\n")

# -------------------- 4. 数据分类（只保留 Key TFs 和 DE-TFs）--------------------
enrichment_data$category <- NA
enrichment_data$category[enrichment_data$Feature %in% key_tfs] <- "Key TFs"
enrichment_data$category[enrichment_data$Feature %in% de_tfs] <- "DE-TFs"

# 移除 Others
enrichment_data <- enrichment_data[!is.na(enrichment_data$category), ]

# 为箱线图创建分组变量
enrichment_data$tf_group <- enrichment_data$Feature
enrichment_data$tf_group[enrichment_data$Feature %in% de_tfs] <- "DE-TFs"

# 设置因子顺序
group_order <- c(key_tfs, "DE-TFs")
enrichment_data$tf_group <- factor(enrichment_data$tf_group, levels = group_order)

cat("\n各组样本量:\n")
print(table(enrichment_data$tf_group))

# -------------------- 5. 定义配色 --------------------
# 为Key TFs分配颜色（循环使用颜色列表）
key_colors <- setNames(
  key_tf_colors[1:length(key_tfs)],
  key_tfs
)
group_colors <- c(key_colors, "DE-TFs" = color_de_tfs)

cat("\n配色方案:\n")
print(group_colors)

# -------------------- 6. 计算统计检验p值 --------------------
cat("\n========== 计算统计检验p值 ==========\n")

# 需要比较的组对：每个 Key TF 与 DE-TFs 比较
pairwise_pvalues <- data.frame()

for (key_tf in key_tfs) {
  g1 <- key_tf
  g2 <- "DE-TFs"
  
  # 检查两组是否都有数据
  if (g1 %in% levels(enrichment_data$tf_group) && 
      g2 %in% levels(enrichment_data$tf_group)) {
    
    d1 <- enrichment_data$Fold_Enrichment[enrichment_data$tf_group == g1]
    d2 <- enrichment_data$Fold_Enrichment[enrichment_data$tf_group == g2]
    
    if (length(d1) >= 2 && length(d2) >= 2) {
      if (stat_test_method == "t.test") {
        test_result <- t.test(d1, d2)
      } else {
        test_result <- wilcox.test(d1, d2)
      }
      
      # 生成p值标签（星号表示）
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
        Group1 = g1,
        Group2 = g2,
        P_value = p_val,
        P_label = p_label,
        Significant = p_val < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(pairwise_pvalues) > 0) {
  cat("\n两两比较p值 (各Key TF vs DE-TFs):\n")
  print(pairwise_pvalues)
} else {
  cat("\n警告：没有足够的数据进行统计检验\n")
}

# -------------------- 7. 计算Y轴范围 --------------------
y_max <- max(enrichment_data$Fold_Enrichment, na.rm = TRUE)
y_min <- min(enrichment_data$Fold_Enrichment, na.rm = TRUE)
y_range <- y_max - y_min

if (is.null(y_axis_limits)) {
  # 自动计算范围，为p值标注留出空间
  n_comparisons <- nrow(pairwise_pvalues)
  if (n_comparisons > 0) {
    y_max_expanded <- y_max + y_range * (0.15 + n_comparisons * pval_step_increase)
  } else {
    y_max_expanded <- y_max * 1.1
  }
  y_axis_limits <- c(y_min, y_max_expanded)
  cat(sprintf("\n自动计算Y轴范围: [%.2f, %.2f]\n", y_axis_limits[1], y_axis_limits[2]))
} else {
  cat(sprintf("\n使用自定义Y轴范围: [%.2f, %.2f]\n", y_axis_limits[1], y_axis_limits[2]))
}

# -------------------- 8. 绘制箱线图 --------------------
cat("\n========== 绘制箱线图 ==========\n")

# 基础箱线图
p <- ggplot(enrichment_data, aes(x = tf_group, y = Fold_Enrichment, fill = tf_group)) +
  geom_boxplot(width = boxplot_box_width, alpha = boxplot_box_alpha, outlier.shape = NA) +
  scale_fill_manual(values = group_colors, name = "TF Group") +
  labs(title = boxplot_title, x = boxplot_x_label, y = boxplot_y_label) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = boxplot_title_size, face = "bold"),
    axis.title = element_text(size = boxplot_axis_title_size),
    axis.text = element_text(size = boxplot_axis_text_size),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  coord_cartesian(ylim = y_axis_limits)

# 添加散点
if (boxplot_show_points) {
  p <- p + geom_jitter(width = boxplot_jitter_width, size = boxplot_point_size,
                       alpha = boxplot_point_alpha, shape = 21, 
                       color = "grey30", aes(fill = tf_group)) +
    scale_fill_manual(values = group_colors, guide = "none")
}

# -------------------- 9. 添加p值标注 --------------------
if (show_pvalue_labels && nrow(pairwise_pvalues) > 0) {
  
  # 获取各组在x轴上的位置
  x_positions <- 1:length(group_order)
  names(x_positions) <- group_order
  
  # 为每个比较分配Y位置
  for (i in 1:nrow(pairwise_pvalues)) {
    comp <- pairwise_pvalues[i, ]
    
    # 获取x坐标
    x1 <- x_positions[comp$Group1]
    x2 <- x_positions[comp$Group2]
    
    # 计算Y位置
    y_position <- y_max + y_range * (pval_bracket_height + i * pval_step_increase)
    
    # 添加标注
    p <- p +
      # 左侧竖线
      annotate("segment", x = x1, xend = x1,
               y = y_position - y_range * pval_bracket_height * 0.3,
               yend = y_position,
               color = pval_line_color, linewidth = pval_line_size) +
      # 右侧竖线
      annotate("segment", x = x2, xend = x2,
               y = y_position - y_range * pval_bracket_height * 0.3,
               yend = y_position,
               color = pval_line_color, linewidth = pval_line_size) +
      # 横线
      annotate("segment", x = x1, xend = x2,
               y = y_position, yend = y_position,
               color = pval_line_color, linewidth = pval_line_size) +
      # p值标签
      annotate("text", x = (x1 + x2) / 2,
               y = y_position + y_range * pval_text_y_offset,
               label = comp$P_label,
               size = pval_label_size, color = pval_label_color)
  }
}

# -------------------- 10. 保存图片 --------------------
ggsave(output_boxplot, p, width = boxplot_width, height = boxplot_height, dpi = 300)
cat(sprintf("\n✓ 箱线图已保存: %s\n", output_boxplot))

# 显示图片
print(p)

# -------------------- 11. 输出统计报告 --------------------
sink(output_stats)

cat("================================================================================\n")
cat("                ChIP-Atlas 富集分析统计报告\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、数据概况\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("总记录数: %d\n", nrow(enrichment_data)))
cat(sprintf("富集倍数阈值: > %.1f\n", fe_threshold))

cat("\n二、TF列表来源\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("Key TFs (%d个): %s\n", length(key_tfs), paste(key_tfs, collapse = ", ")))
cat(sprintf("DE-TFs (%d个): %s\n", length(de_tfs), paste(de_tfs, collapse = ", ")))

cat("\n三、各组样本量\n")
cat("--------------------------------------------------------------------------------\n")
print(table(enrichment_data$tf_group))

cat("\n四、各组富集分数统计\n")
cat("--------------------------------------------------------------------------------\n")
stats_summary <- enrichment_data %>%
  group_by(tf_group) %>%
  summarise(
    Median = median(Fold_Enrichment),
    Mean = mean(Fold_Enrichment),
    SD = sd(Fold_Enrichment),
    Min = min(Fold_Enrichment),
    Max = max(Fold_Enrichment),
    n = n(),
    .groups = "drop"
  )
print(stats_summary)

if (nrow(pairwise_pvalues) > 0) {
  cat("\n五、统计检验结果 (各Key TF vs DE-TFs)\n")
  cat("--------------------------------------------------------------------------------\n")
  cat(sprintf("检验方法: %s\n", stat_test_method))
  cat("\np值标注说明: **** = p < 0.0001, *** = p < 0.001, ** = p < 0.01, * = p < 0.05, ns = 不显著\n\n")
  for (i in 1:nrow(pairwise_pvalues)) {
    cat(sprintf("  %s vs %s: %s (p = %.2e)\n",
                pairwise_pvalues$Group1[i], pairwise_pvalues$Group2[i],
                pairwise_pvalues$P_label[i], pairwise_pvalues$P_value[i]))
  }
} else {
  cat("\n五、统计检验结果\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("无足够数据进行统计检验\n")
}

cat("\n\n六、图形参数\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("Y轴范围: [%.2f, %.2f]\n", y_axis_limits[1], y_axis_limits[2]))
cat(sprintf("图片尺寸: %.1f × %.1f 英寸\n", boxplot_width, boxplot_height))

sink()

cat(sprintf("\n✓ 统计报告已保存: %s\n", output_stats))

cat("\n========================================\n")
cat("分析完成！\n")
cat("========================================\n")