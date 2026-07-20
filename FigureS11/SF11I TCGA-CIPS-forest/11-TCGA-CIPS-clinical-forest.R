# ============================================================================
# COX回归森林图绘制脚本
# 功能：从COX回归结果绘制森林图
# ============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
work_dir <- "D:/zsy/SX/Fomal-final/11-TCGA-clinical/survival_analysis-TNM"
setwd(work_dir)

# 创建输出文件夹
output_dir <- "forest_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 输入文件
cox_results_file <- file.path(output_dir, "COX_regression_results.csv")  # COX回归结果文件

# 如果没有该文件，尝试其他位置
if (!file.exists(cox_results_file)) {
  cox_results_file <- "COX_regression_results.csv"
}

# 输出文件
output_forest <- file.path(output_dir, "forest_plot.pdf")
output_forest_png <- file.path(output_dir, "forest_plot.png")
output_forest_ggplot <- file.path(output_dir, "forest_plot_ggplot.pdf")
output_simple <- file.path(output_dir, "forest_plot_significant.pdf")
output_table <- file.path(output_dir, "forest_data_summary.csv")

# ====================== 加载必要的包 ======================
cat("\n========== 加载必要的R包 ==========\n")

packages <- c("forestplot", "ggplot2", "dplyr", "tidyr", "gridExtra")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取COX结果 ======================
cat("\n========== 读取COX回归结果 ==========\n")

if (!file.exists(cox_results_file)) {
  stop("COX结果文件不存在: ", cox_results_file, "\n请先运行COX分析生成结果文件")
}

cox_data <- read.csv(cox_results_file, stringsAsFactors = FALSE)
cat("COX结果数据维度:", dim(cox_data), "\n")
cat("列名:\n")
print(colnames(cox_data))

# 显示前几行数据
cat("\n数据前几行:\n")
print(head(cox_data))

# ====================== 数据准备 ======================
cat("\n========== 准备绘图数据 ==========\n")

# 选择要绘制的模型（通常选择多因素模型）
if ("Type" %in% colnames(cox_data)) {
  # 优先使用多因素模型
  if ("Multivariate" %in% cox_data$Type) {
    plot_data <- cox_data[cox_data$Type == "Multivariate", ]
    cat("使用多因素模型结果\n")
  } else {
    plot_data <- cox_data[cox_data$Type == "Univariate", ]
    cat("使用单因素模型结果\n")
  }
} else if ("Model" %in% colnames(cox_data)) {
  plot_data <- cox_data[cox_data$Model == "Stage Model", ]
  cat("使用Stage模型结果\n")
} else {
  plot_data <- cox_data
  cat("使用全部结果\n")
}

# 确保有数据
if (nrow(plot_data) == 0) {
  stop("没有找到绘图数据，请检查COX结果文件")
}

# 清理变量名，使其更易读
plot_data$Variable_clean <- plot_data$Variable

# 变量名映射（根据您的数据调整）
name_mapping <- list(
  "ssGSEA_Score" = "ssGSEA Score (per 0.1 increase)",
  "age" = "Age (per year)",
  "genderMale" = "Gender (Male vs Female)",
  "genderFemale" = "Gender (Female vs Male)",
  "Stage_combinedII" = "Stage II vs I",
  "Stage_combinedIII" = "Stage III vs I",
  "Stage_combinedIV" = "Stage IV vs I",
  "stage_groupAdvanced" = "Stage (Advanced vs Early)",
  "N_groupN\\+" = "Lymph Node (N+ vs N0)",
  "M_simpleM1" = "Metastasis (M1 vs M0)",
  "grade_groupHigh Grade" = "Grade (High vs Low)",
  "T_groupT3-4" = "Tumor (T3-4 vs T1-2)"
)

# 应用变量名映射
for (orig in names(name_mapping)) {
  plot_data$Variable_clean <- gsub(orig, name_mapping[[orig]], plot_data$Variable_clean)
}

# 创建用于森林图的数据框
forest_data <- data.frame(
  Variable = plot_data$Variable_clean,
  HR = as.numeric(plot_data$HR),
  Lower = as.numeric(plot_data$HR_lower),
  Upper = as.numeric(plot_data$HR_upper),
  P = as.numeric(plot_data$P_value),
  stringsAsFactors = FALSE
)

# 移除无效行（HR为NA或无限值）
forest_data <- forest_data[!is.na(forest_data$HR) & is.finite(forest_data$HR), ]

# 按HR排序（可选）
forest_data <- forest_data[order(forest_data$HR, decreasing = TRUE), ]

# 添加显著性标记
forest_data$Significant <- ifelse(forest_data$P < 0.05, 
                                  ifelse(forest_data$P < 0.01, 
                                         ifelse(forest_data$P < 0.001, "***", "**"), "*"), 
                                  "ns")
forest_data$P_text <- ifelse(forest_data$P < 0.001, "<0.001", 
                             ifelse(forest_data$P < 0.01, sprintf("%.3f", forest_data$P),
                                    sprintf("%.3f", forest_data$P)))

cat("绘图数据预览:\n")
print(forest_data)

# ====================== 方法1：使用forestplot包 ======================
cat("\n========== 绘制森林图（方法1：forestplot包）==========\n")

# 准备表格数据
tabletext <- cbind(
  c("Variable", forest_data$Variable),
  c("HR (95% CI)", 
    paste0(sprintf("%.2f", forest_data$HR), " (", 
           sprintf("%.2f", forest_data$Lower), "-", 
           sprintf("%.2f", forest_data$Upper), ")")),
  c("P-value", forest_data$P_text)
)

# 添加显著性标记到变量名
for (i in 1:nrow(forest_data)) {
  if (forest_data$Significant[i] != "ns") {
    forest_data$Variable[i] <- paste0(forest_data$Variable[i], forest_data$Significant[i])
  }
}

tabletext <- cbind(
  c("Variable", forest_data$Variable),
  c("HR (95% CI)", 
    paste0(sprintf("%.2f", forest_data$HR), " (", 
           sprintf("%.2f", forest_data$Lower), "-", 
           sprintf("%.2f", forest_data$Upper), ")")),
  c("P-value", forest_data$P_text)
)

# 创建森林图
tryCatch({
  pdf(output_forest, width = 10, height = max(6, nrow(forest_data) * 0.4))
  
  forestplot(labeltext = tabletext,
             mean = c(NA, forest_data$HR),
             lower = c(NA, forest_data$Lower),
             upper = c(NA, forest_data$Upper),
             graph.pos = 2,
             graphwidth = unit(0.4, "npc"),
             col = fpColors(box = "#C60036", lines = "black", zero = "gray50"),
             xlab = "Hazard Ratio (95% CI)",
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab = gpar(cex = 0.9)),
             lwd.zero = 1,
             lwd.ci = 1.5,
             boxsize = 0.2,
             zero = 1,
             colgap = unit(5, "mm"),
             lineheight = "auto",
             title = "Multivariate Cox Regression Analysis")
  
  dev.off()
  cat("森林图已保存到:", output_forest, "\n")
}, error = function(e) {
  cat("forestplot包绘制失败:", e$message, "\n")
})

# ====================== 方法2：使用ggplot2 ======================
cat("\n========== 绘制森林图（方法2：ggplot2）==========\n")

# 设置因子顺序
forest_data$Variable <- factor(forest_data$Variable, 
                               levels = rev(forest_data$Variable))

# 创建ggplot森林图
p <- ggplot(forest_data, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.2, color = "gray30", linewidth = 0.8) +
  geom_point(aes(color = Significant), size = 3) +
  scale_color_manual(values = c("***" = "#C60036", "**" = "#C60036", 
                                "*" = "#C60036", "ns" = "gray60"),
                     labels = c("***" = "p < 0.001", "**" = "p < 0.01", 
                                "*" = "p < 0.05", "ns" = "Not significant")) +
  scale_x_log10() +
  labs(title = "Multivariate Cox Regression Analysis",
       x = "Hazard Ratio (95% CI)",
       y = "",
       color = "Significance") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "top")

# 添加HR值标签
p <- p + geom_text(aes(label = sprintf("%.2f", HR), x = Upper + 0.05), 
                   hjust = 0, size = 3)

# 保存PNG格式
png(output_forest_png, width = 10, height = max(6, nrow(forest_data) * 0.4), 
    units = "in", res = 300)
print(p)
dev.off()
cat("森林图（PNG）已保存到:", output_forest_png, "\n")

# 保存PDF
pdf(output_forest_ggplot, width = 10, height = max(6, nrow(forest_data) * 0.4))
print(p)
dev.off()
cat("森林图（ggplot PDF）已保存到:", output_forest_ggplot, "\n")

# ====================== 方法3：简化版森林图 ======================
cat("\n========== 绘制简化版森林图 ==========\n")

# 创建简化版数据（只显示显著变量）
sig_data <- forest_data[forest_data$Significant != "ns", ]

if (nrow(sig_data) > 0) {
  sig_data$Variable <- factor(sig_data$Variable, levels = rev(sig_data$Variable))
  
  p_simple <- ggplot(sig_data, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "gray30") +
    geom_point(color = "#C60036", size = 3) +
    scale_x_log10() +
    labs(title = "Significant Factors in Multivariate Cox Regression",
         x = "Hazard Ratio (95% CI)",
         y = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 10))
  
  pdf(output_simple, width = 8, height = max(5, nrow(sig_data) * 0.4))
  print(p_simple)
  dev.off()
  cat("简化版森林图（仅显著变量）已保存到:", output_simple, "\n")
} else {
  cat("没有显著变量，跳过简化版森林图\n")
}

# ====================== 生成数据汇总表 ======================
cat("\n========== 生成数据汇总表 ==========\n")

summary_table <- forest_data %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, Lower, Upper),
    P_value = ifelse(P < 0.001, "<0.001", sprintf("%.3f", P))
  ) %>%
  select(Variable, HR_CI, P_value, Significant)

write.csv(summary_table, file = output_table, row.names = FALSE)
cat("数据汇总表已保存到:", output_table, "\n")

# ====================== 打印结果摘要 ======================
cat("\n========== 结果摘要 ==========\n")
cat("\n森林图数据:\n")
print(summary_table)

cat("\n森林图文件已生成:\n")
cat("  1. ", output_forest, " - forestplot包版本\n")
cat("  2. ", output_forest_ggplot, " - ggplot2版本\n")
cat("  3. ", output_forest_png, " - PNG格式\n")
if (nrow(sig_data) > 0) {
  cat("  4. ", output_simple, " - 仅显著变量版本\n")
}
cat("  5. ", output_table, " - 数据汇总表\n")

cat("\n=================== 森林图绘制完成 ===================\n")