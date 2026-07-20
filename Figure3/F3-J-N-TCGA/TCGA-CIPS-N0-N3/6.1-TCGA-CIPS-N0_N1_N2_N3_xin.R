# =============================================================================
# 生信分析脚本：ssGSEA评分与N分级关联分析（N0, N1, N2, N3四组分析）
# 功能：比较N0、N1、N2、N3四组的评分差异，绘制箱线图并进行多重比较
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/13-TCGA-clinical"  # 请修改此路径
setwd(work_dir)

# 2. 创建输出文件夹
output_dir <- "N_stage_analysis_N0_N1_N2_N3"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 3. 输入文件名设置
clinical_file <- "TCGA_clinical.csv"        # 临床数据文件
score_file <- "ssGSEA_scores_diagnostic.csv" # ssGSEA评分文件

# 4. 列名设置（可自定义）
clinical_id_col <- "Id"           # 样品ID列
clinical_N <- "N"                  # N分期列

# ssGSEA评分文件列名
score_id_col <- "Sample"           # 样品ID列
score_value_col <- "ssGSEA_Score"  # 评分列

# 5. N分期映射设置（将原始分期映射到标准N0、N1、N2、N3）
# 请根据实际数据中的N分期值进行修改
n_mapping_rules <- list(
  "N0" = c("N0", "N0 (i-)", "N0 (i+)", "N0 (mol+)"),
  "N1" = c("N1", "N1a", "N1b", "N1c"),
  "N2" = c("N2", "N2a", "N2b"),
  "N3" = c("N3", "N3a", "N3b", "N3c")
)

# 6. 组别顺序和标签设置
group_order <- c("N0", "N1", "N2", "N3")                     # 组别显示顺序
group_labels <- c("N0", "N1", "N2", "N3")                    # 组别显示标签
group_colors <- c("#25377F", "#2AA7DE", "#2E86AB", "#C60036") # N0, N1, N2, N3的颜色

# 7. 箱线图设置（完全可自定义）
boxplot_width <- 10
boxplot_height <- 10
boxplot_title <- "ssGSEA Score by N Stage (N0, N1, N2, N3)"
x_label <- "N Stage"
y_label <- "CIPS Expression ssGSEA Score"
legend_title <- "N Stage"

# 点设置
show_points <- TRUE                                  # 是否显示散点
point_jitter <- 0.2                                   # 点抖动宽度
point_size <- 1.5                                     # 点大小
point_alpha <- 0.3                                    # 点透明度

# 均值点设置
show_mean_label <- TRUE                               # 是否显示均值标签
mean_label_digits <- 3                                # 均值小数位数
mean_label_vjust <- -0.2                             # 均值标签垂直位置
mean_label_hjust <- 0.5                                # 均值标签水平位置
mean_label_size <- 8                                 # 均值标签大小
mean_label_color <- "black"                            # 均值标签颜色
mean_label_fontface <- "plain"                         # 均值标签字体

# ----------------------------- Y轴范围参数（仿第二个代码）-----------------------------------
# 设置为NULL表示自动计算并自动为p值标注留出空间
# 设置数值如c(-2, 3)表示固定范围
y_axis_limits <- NULL          # 修改这里：NULL=自动，c(ymin, ymax)=固定范围
y_axis_breaks <- NULL                                  # Y轴刻度
y_axis_expand <- c(0.3, 1.0)                          # Y轴扩展比例（增加顶部空间用于显示p值）

# ----------------------------- p值标注参数（仿第二个代码）-----------------------------------
pval_label_size <- 9                                   # p值标签字体大小
pval_label_color <- "black"                            # p值标签颜色
pval_line_color <- "black"                             # 标注线颜色
pval_line_size <- 0.7                                  # 标注线粗细
pval_bracket_height <- 0.05                            # 括号高度（相对于数据范围的比例）
pval_text_y_offset <- 0.1                             # 文字相对于括号的偏移量
pval_step_increase <- 0.17                             # 多个比较时的高度递增

# 主题设置
plot_title_hjust <- 0.5                                # 主标题水平对齐
plot_title_size <- 20                                   # 主标题大小
plot_title_face <- "bold"                               # 主标题字体
plot_title_margin <- 10                                 # 主标题底部边距

subtitle_hjust <- 0.5                                   # 副标题水平对齐
subtitle_size <- 18                                     # 副标题大小
subtitle_color <- "darkred"                             # 副标题颜色

x_axis_title_size <- 32                                 # X轴标题大小
x_axis_text_size <- 32                                 # X轴刻度大小
x_axis_text_angle <- 0                                  # X轴刻度角度
x_axis_text_hjust <- 0.5                                # X轴刻度水平对齐
x_axis_text_vjust <- 0.5                                # X轴刻度垂直对齐

y_axis_title_size <- 32                                 # Y轴标题大小
y_axis_title_margin <- 5                                # Y轴标题右边距
y_axis_text_size <- 32                                  # Y轴刻度大小
y_axis_text_angle <- 0                                  # Y轴刻度角度
y_axis_text_hjust <- 0.5                                # Y轴刻度水平对齐
y_axis_text_vjust <- 0.5                                # Y轴刻度垂直对齐

# 多重比较p值标注设置
p_value_method <- "t.test"                         # 检验方法: "wilcox.test" (Mann-Whitney U检验) 或 "t.test"
p_adjust_method <- "none"                                 # p值校正方法: "none", "BH", "bonferroni", "holm"等
show_significance_only <- FALSE                         # 是否只显示有显著差异的比较（p < 0.05）
p_value_digits <- 3                                     # p值显示的小数位数

# 8. 堆叠柱状图设置
barplot_width <- 10
barplot_height <- 6
barplot_colors <- c("#2E8B57", "#DC143C")              # 低表达（绿色），高表达（红色）
barplot_title <- "Proportion of High/Low Expression by N Stage"
barplot_xlab <- "N Stage"
barplot_ylab <- "Number of Samples"
barplot_legend_title <- "Expression Level"

# 9. 输出文件设置
output_boxplot <- file.path(output_dir, "N_stage_boxplot_N0_N1_N2_N3.pdf")
output_barplot <- file.path(output_dir, "N_stage_barplot_stacked_N0_N1_N2_N3.pdf")
output_stats <- file.path(output_dir, "N_stage_statistics.csv")
output_data <- file.path(output_dir, "merged_N_stage_data.csv")
output_mapping <- file.path(output_dir, "N_stage_mapping.csv")
output_summary <- file.path(output_dir, "analysis_summary.txt")
output_pairwise_stats <- file.path(output_dir, "pairwise_comparison_stats.csv")

# ====================== 加载必要的包 ======================
cat("\n========== 加载必要的R包 ==========\n")

packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "ggpubr", "rstatix")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取数据 ======================
cat("\n========== 读取数据 ==========\n")

clinical_data <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("临床数据维度: %d 行, %d 列\n", nrow(clinical_data), ncol(clinical_data)))

score_data <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("评分数据维度: %d 行, %d 列\n", nrow(score_data), ncol(score_data)))

# ====================== ID标准化和匹配 ======================
cat("\n========== ID标准化和匹配 ==========\n")

clinical_data$PatientID <- substr(clinical_data[[clinical_id_col]], 1, 12)
score_data$PatientID <- substr(score_data[[score_id_col]], 1, 12)

common_patients <- intersect(clinical_data$PatientID, score_data$PatientID)
cat(sprintf("共同患者数量: %d\n", length(common_patients)))

if (length(common_patients) == 0) {
  stop("错误：没有匹配的患者ID！请检查ID格式。")
}

merged_data <- clinical_data %>%
  filter(PatientID %in% common_patients) %>%
  left_join(score_data[, c("PatientID", score_value_col)], by = "PatientID") %>%
  distinct(PatientID, .keep_all = TRUE)

cat(sprintf("合并后数据维度: %d 行, %d 列\n", nrow(merged_data), ncol(merged_data)))

# ====================== 数据预处理 ======================
cat("\n========== 数据预处理 ==========\n")

colnames(merged_data)[colnames(merged_data) == score_value_col] <- "ssGSEA_Score"

# ====================== N分期处理（映射到N0、N1、N2、N3）======================
cat("\n========== N分期处理（映射到N0、N1、N2、N3）==========\n")

# 显示原始N分期分布
cat("原始N分期分布:\n")
original_n_table <- table(merged_data[[clinical_N]], useNA = "ifany")
print(original_n_table)

# 创建N分组（映射到标准分期）
merged_data$N_group <- NA_character_

for (std_n in names(n_mapping_rules)) {
  merged_data$N_group[merged_data[[clinical_N]] %in% n_mapping_rules[[std_n]]] <- std_n
}

# 设置组别顺序和标签
merged_data$N_group <- factor(merged_data$N_group, levels = group_order)
if (!is.null(group_labels) && length(group_labels) == length(levels(merged_data$N_group))) {
  levels(merged_data$N_group) <- group_labels
}

cat("\nN分组分布:\n")
n_group_table <- table(merged_data$N_group, useNA = "ifany")
print(n_group_table)

# 统计每组样本数
n_group_counts <- merged_data %>%
  filter(!is.na(N_group)) %>%
  group_by(N_group) %>%
  summarise(n = n()) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

cat("\n每组样本数:\n")
print(n_group_counts)

# 保存N分期映射关系
n_mapping_df <- merged_data %>%
  select(original_N = all_of(clinical_N), N_group) %>%
  distinct() %>%
  arrange(original_N)

write.csv(n_mapping_df, file = output_mapping, row.names = FALSE)
cat(sprintf("\n完整映射关系已保存到: %s\n", output_mapping))

# 检查是否有足够的样本进行分析
valid_data <- merged_data %>% filter(!is.na(N_group))
if (nrow(valid_data) == 0) {
  stop("错误：没有有效的N分期数据！")
}

if (n_distinct(valid_data$N_group) < 2) {
  stop("错误：N分组后只有一个组，无法进行比较！")
}

# ====================== 评分分组（用于堆叠柱状图）======================
cat("\n========== 评分分组 ==========\n")

# 根据评分中位数将患者分为高低两组
score_cutoff <- median(merged_data$ssGSEA_Score, na.rm = TRUE)
merged_data$Score_group <- ifelse(merged_data$ssGSEA_Score >= score_cutoff, "High", "Low")
merged_data$Score_group <- factor(merged_data$Score_group, levels = c("Low", "High"))
cat(sprintf("使用中位数分组 (cutoff = %.4f)\n", score_cutoff))

cat("评分分组统计:\n")
print(table(merged_data$Score_group))

# ====================== 准备堆叠柱状图数据 ======================
cat("\n========== 准备堆叠柱状图数据 ==========\n")

barplot_data <- merged_data %>%
  filter(!is.na(N_group)) %>%
  group_by(N_group, Score_group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(N_group, Score_group, fill = list(count = 0))

barplot_summary <- barplot_data %>%
  group_by(N_group) %>%
  mutate(
    total = sum(count),
    proportion = count / total,
    percentage = paste0(round(100 * proportion, 1), "%")
  )

cat("\n每组高/低表达样本数:\n")
print(barplot_summary)

# ====================== 绘制堆叠柱状图 ======================
cat("\n========== 绘制堆叠柱状图 ==========\n")

p_bar <- ggplot(barplot_data, aes(x = N_group, y = count, fill = Score_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = barplot_colors, 
                    labels = c("Low Expression", "High Expression")) +
  labs(title = barplot_title,
       x = barplot_xlab,
       y = barplot_ylab,
       fill = barplot_legend_title) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = "top")

# 添加数量标签
p_bar <- p_bar + 
  geom_text(data = barplot_data %>% filter(count > 0),
            aes(label = count, y = count),
            position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold")

# 保存堆叠柱状图
pdf(output_barplot, width = barplot_width, height = barplot_height)
print(p_bar)
dev.off()
cat("堆叠柱状图已保存到:", output_barplot, "\n")

# ====================== 多重比较和箱线图绘制 ======================
cat("\n========== 多重比较和箱线图绘制 ==========\n")

boxplot_data <- merged_data %>% filter(!is.na(N_group))

if (nrow(boxplot_data) > 0) {
  # 计算每组统计量
  stats_data <- boxplot_data %>%
    group_by(N_group) %>%
    summarise(
      n = n(),
      median = median(ssGSEA_Score, na.rm = TRUE),
      q1 = quantile(ssGSEA_Score, 0.25, na.rm = TRUE),
      q3 = quantile(ssGSEA_Score, 0.75, na.rm = TRUE),
      mean = mean(ssGSEA_Score, na.rm = TRUE),
      sd = sd(ssGSEA_Score, na.rm = TRUE)
    )
  
  cat("\n每组统计量:\n")
  print(stats_data)
  
  # Kruskal-Wallis检验（非参数ANOVA）
  kruskal_test <- kruskal.test(ssGSEA_Score ~ N_group, data = boxplot_data)
  kw_p_value <- kruskal_test$p.value
  
  cat("\nKruskal-Wallis检验结果:\n")
  cat(sprintf("  Kruskal-Wallis检验 p值: %.3e\n", kw_p_value))
  
  # 计算所有组间的两两比较
  # 获取所有组别组合
  groups <- levels(boxplot_data$N_group)
  group_pairs <- combn(groups, 2, simplify = FALSE)
  
  # 存储两两比较结果
  pairwise_results <- data.frame(
    Group1 = character(),
    Group2 = character(),
    P_value = numeric(),
    Adjusted_P_value = numeric(),
    Significant = logical(),
    P_label = character(),  # 新增：格式化后的p值标签
    stringsAsFactors = FALSE
  )
  
  # 计算每对组别的p值
  for (pair in group_pairs) {
    group1_data <- boxplot_data$ssGSEA_Score[boxplot_data$N_group == pair[1]]
    group2_data <- boxplot_data$ssGSEA_Score[boxplot_data$N_group == pair[2]]
    
    if (length(group1_data) > 0 && length(group2_data) > 0) {
      # 根据选择的方法进行检验
      if (p_value_method == "t.test") {
        test_result <- t.test(group1_data, group2_data)
      } else {
        test_result <- wilcox.test(group1_data, group2_data)
      }
      
      p_raw <- test_result$p.value
      
      # 格式化p值标签
      if (p_raw < 0.001) {
        p_label <- "p < 0.001"
      } else {
        p_label <- sprintf("p = %.3f", p_raw)
      }
      
      pairwise_results <- rbind(pairwise_results, data.frame(
        Group1 = pair[1],
        Group2 = pair[2],
        P_value = p_raw,
        Adjusted_P_value = NA,
        Significant = FALSE,
        P_label = p_label,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 对p值进行多重校正
  if (nrow(pairwise_results) > 0 && p_adjust_method != "none") {
    pairwise_results$Adjusted_P_value <- p.adjust(pairwise_results$P_value, method = p_adjust_method)
    pairwise_results$Significant <- pairwise_results$Adjusted_P_value < 0.05
    
    # 更新标签
    for (i in 1:nrow(pairwise_results)) {
      if (pairwise_results$Adjusted_P_value[i] < 0.001) {
        pairwise_results$P_label[i] <- "p < 0.001"
      } else {
        pairwise_results$P_label[i] <- sprintf("p = %.3f", pairwise_results$Adjusted_P_value[i])
      }
    }
  } else if (nrow(pairwise_results) > 0) {
    pairwise_results$Adjusted_P_value <- pairwise_results$P_value
    pairwise_results$Significant <- pairwise_results$P_value < 0.05
  }
  
  # 保存两两比较结果
  write.csv(pairwise_results, file = output_pairwise_stats, row.names = FALSE)
  cat(sprintf("\n两两比较结果已保存到: %s\n", output_pairwise_stats))
  
  cat("\n两两比较结果（", ifelse(p_adjust_method == "none", "未校正", paste0(p_adjust_method, "校正")), "）:\n", sep = "")
  print(pairwise_results)
  
  # ====================== 计算Y轴范围（仿第二个代码）======================
  cat("\n========== 计算Y轴范围 ==========\n")
  
  y_max_data <- max(boxplot_data$ssGSEA_Score, na.rm = TRUE)
  y_min_data <- min(boxplot_data$ssGSEA_Score, na.rm = TRUE)
  y_range_data <- y_max_data - y_min_data
  
  # 筛选需要显示的比较
  if (show_significance_only) {
    comparisons_to_show <- pairwise_results[pairwise_results$Significant, ]
  } else {
    comparisons_to_show <- pairwise_results
  }
  
  n_comparisons <- nrow(comparisons_to_show)
  
  if (is.null(y_axis_limits)) {
    # 自动计算范围，为p值标注留出空间
    y_max_expanded <- y_max_data + y_range_data * (0.15 + n_comparisons * pval_step_increase)
    y_axis_limits <- c(y_min_data, y_max_expanded)
    cat(sprintf("自动计算Y轴范围: [%.3f, %.3f]\n", y_axis_limits[1], y_axis_limits[2]))
  } else {
    # 使用用户自定义的范围
    cat(sprintf("使用自定义Y轴范围: [%.3f, %.3f]\n", y_axis_limits[1], y_axis_limits[2]))
  }
  
  # ====================== 创建基础箱线图 ======================
  p <- ggplot(boxplot_data, aes(x = N_group, y = ssGSEA_Score, fill = N_group))
  
  # 箱线图
  p <- p + geom_boxplot(width = 0.8, alpha = 0.7, outlier.shape = NA)
  
  # 散点
  if (show_points) {
    p <- p + geom_jitter(width = point_jitter, size = point_size, 
                         alpha = point_alpha, aes(color = N_group))
  }
  
  # 均值点
  p <- p + stat_summary(fun = mean, geom = "point", 
                        shape = 18, size = 4, color = "black")
  
  # 均值标签
  if (show_mean_label) {
    p <- p + stat_summary(
      fun = mean, 
      geom = "text", 
      aes(label = sprintf(paste0("%.", mean_label_digits, "f"), ..y..)),
      vjust = mean_label_vjust,
      hjust = mean_label_hjust,
      size = mean_label_size,
      color = mean_label_color,
      fontface = mean_label_fontface
    )
  }
  
  # 颜色设置
  n_groups <- length(unique(boxplot_data$N_group))
  p <- p + scale_fill_manual(values = group_colors[1:n_groups]) +
    scale_color_manual(values = group_colors[1:n_groups])
  
  # 设置Y轴范围（使用coord_cartesian，仿第二个代码）
  p <- p + coord_cartesian(ylim = y_axis_limits)
  
  # ====================== 添加p值标注（仿第二个代码）======================
  if (n_comparisons > 0) {
    
    # 获取各组在x轴上的位置
    x_positions <- 1:length(group_order)
    names(x_positions) <- group_order
    
    # 为每个比较添加标注
    for (i in 1:n_comparisons) {
      comp <- comparisons_to_show[i, ]
      
      # 获取x坐标
      x1 <- x_positions[comp$Group1]
      x2 <- x_positions[comp$Group2]
      
      # 计算Y位置（根据比较顺序递增）
      y_position <- y_max_data + y_range_data * (pval_bracket_height + i * pval_step_increase)
      
      # 确保Y位置在图表范围内
      if (y_position > y_axis_limits[2]) {
        y_position <- y_axis_limits[2] - y_range_data * 0.05
      }
      
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
  
  # 添加全局p值作为副标题
  if (kw_p_value < 0.001) {
    subtitle <- sprintf("Kruskal-Wallis test p < 0.001")
  } else {
    subtitle <- sprintf("Kruskal-Wallis test p = %.3f", kw_p_value)
  }
  
  # 主题设置
  p <- p + labs(
    title = boxplot_title,
    subtitle = subtitle,
    x = x_label,
    y = y_label,
    fill = legend_title,
    color = legend_title
  ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = plot_title_hjust,
        size = plot_title_size,
        face = plot_title_face,
        margin = margin(b = plot_title_margin)
      ),
      plot.subtitle = element_text(
        hjust = subtitle_hjust,
        size = subtitle_size,
        color = subtitle_color
      ),
      axis.title.x = element_text(
        size = x_axis_title_size
      ),
      axis.text.x = element_text(
        size = x_axis_text_size,
        angle = x_axis_text_angle,
        hjust = x_axis_text_hjust,
        vjust = x_axis_text_vjust
      ),
      axis.title.y = element_text(
        size = y_axis_title_size,
        margin = margin(r = y_axis_title_margin)
      ),
      axis.text.y = element_text(
        size = y_axis_text_size,
        angle = y_axis_text_angle,
        hjust = y_axis_text_hjust,
        vjust = y_axis_text_vjust
      ),
      legend.position = "right"
    )
  
  # 保存箱线图
  pdf(output_boxplot, width = boxplot_width, height = boxplot_height)
  print(p)
  dev.off()
  cat("箱线图已保存到:", output_boxplot, "\n")
  
  # 保存统计结果
  stats_df <- data.frame(
    Test = "Kruskal-Wallis test",
    P_value = kw_p_value,
    Groups_compared = paste(groups, collapse = " vs "),
    Pairwise_comparison_method = ifelse(p_value_method == "t.test", "t-test", "Wilcoxon test"),
    P_adjust_method = p_adjust_method
  )
  write.csv(stats_df, file = output_stats, row.names = FALSE)
}

# ====================== 保存合并后的数据 ======================
cat("\n========== 保存合并数据 ==========\n")

write.csv(merged_data, file = output_data, row.names = FALSE)
cat("合并后的数据已保存到:", output_data, "\n")

# ====================== 生成分析总结 ======================
cat("\n========== 生成分析总结 ==========\n")

sink(output_summary)

cat("========================================\n")
cat("   ssGSEA评分与N分级分析报告\n")
cat("      （N0, N1, N2, N3四组分析）\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time()), "\n\n")

cat("一、数据概况\n")
cat("------------\n")
cat(sprintf("总样本数: %d\n", nrow(merged_data)))
cat(sprintf("有效N分期样本数: %d\n", sum(!is.na(merged_data$N_group))))

cat("\n二、N分期分布\n")
cat("--------------\n")
cat("原始N分期分布:\n")
print(original_n_table)
cat("\nN分组分布:\n")
print(n_group_table)

cat("\n三、每组高/低表达分布\n")
cat("----------------------\n")
print(barplot_summary)

cat("\n四、箱线图统计量\n")
cat("----------------\n")
print(stats_data)

cat("\n五、统计检验结果\n")
cat("----------------\n")
cat(sprintf("Kruskal-Wallis检验 p值: %.3e\n", kw_p_value))
if (kw_p_value < 0.05) {
  cat("  注：组间存在显著差异\n")
} else {
  cat("  注：组间无显著差异\n")
}

cat(sprintf("\n两两比较方法: %s\n", ifelse(p_value_method == "t.test", "t-test", "Wilcoxon rank-sum test")))
cat(sprintf("多重比较校正方法: %s\n", ifelse(p_adjust_method == "none", "未校正", p_adjust_method)))

cat("\n六、两两比较结果\n")
cat("----------------\n")
print(pairwise_results[, c("Group1", "Group2", "P_label", "Significant")])

cat("\n七、Y轴范围设置\n")
cat("----------------\n")
if (is.null(y_axis_limits)) {
  cat("Y轴范围: 自动计算\n")
} else {
  cat(sprintf("Y轴范围: [%.3f, %.3f] (自定义)\n", y_axis_limits[1], y_axis_limits[2]))
}

cat("\n八、输出文件列表\n")
cat("----------------\n")
cat("1.", basename(output_boxplot), "- N分期箱线图 (N0, N1, N2, N3，含多重比较p值)\n")
cat("2.", basename(output_barplot), "- 堆叠柱状图\n")
cat("3.", basename(output_stats), "- 统计结果\n")
cat("4.", basename(output_data), "- 合并后的数据\n")
cat("5.", basename(output_mapping), "- N分期映射关系\n")
cat("6.", basename(output_pairwise_stats), "- 两两比较统计结果\n")

cat("\n========================================\n")
cat("                完成\n")
cat("========================================\n")

sink()

cat("分析总结已保存到:", output_summary, "\n")
cat("\n=================== 分析完成 ===================\n")