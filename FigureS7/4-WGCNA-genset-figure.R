# ============================================================================
# WGCNA分析 - 模块3：模块与功能基因集重叠分析（完整修复版）
# ============================================================================

# 清空环境变量
rm(list = ls())

# ========================== 工作路径设定 ====================================
work_dir <- "D:/zsy/SX/Fomal-final/27.5-CAF-tumor-WGCNA"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")

# ========================== 创建输出目录 ====================================
output_dir <- "WGCNA_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 创建子目录
gsea_dir <- file.path(output_dir, "gene_set_enrichment")
if (!dir.exists(gsea_dir)) {
  dir.create(gsea_dir, recursive = TRUE)
}

# ========================== 参数设置 =========================================
# 可视化参数
plot_width <- 12
plot_height <- 10

# 热图参数（可自定义）
heatmap_colors <- c("#2C7BB6", "white", "#D7191C")  # 热图颜色
heatmap_font_row <- 0.8    # 行标签字体大小
heatmap_font_col <- 1.0    # 列标签字体大小
heatmap_margins <- c(8, 8) # 边距

# 输出文件
output_barplot <- file.path(gsea_dir, "module_enrichment_fold_enrichment.pdf")
output_bubble <- file.path(gsea_dir, "module_enrichment_bubble.pdf")
output_summary <- file.path(gsea_dir, "module_enrichment_summary.csv")

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

library(ggplot2)
library(dplyr)
library(reshape2)

# 检查并加载pheatmap包（用于带聚类的热图）
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
if (has_pheatmap) {
  library(pheatmap)
  cat("✓ pheatmap包已加载\n")
} else {
  cat("⚠ pheatmap包未安装，将使用基础heatmap（无聚类版本）\n")
  cat("  如需带聚类的热图，请运行: install.packages('pheatmap')\n")
}

cat("✓ 所有包已加载\n")

# ========================== 读取富集结果 ====================================
cat("\n========== 读取富集结果 ==========\n")

enrichment_file <- file.path(gsea_dir, "module_gene_set_enrichment.csv")

if (!file.exists(enrichment_file)) {
  stop("错误：未找到富集结果文件，请先运行模块3")
}

enrichment_results <- read.csv(enrichment_file, stringsAsFactors = FALSE)
cat(sprintf("富集结果维度: %d 行, %d 列\n", nrow(enrichment_results), ncol(enrichment_results)))

# ========================== 按模块汇总统计 ====================================
cat("\n========== 按模块汇总统计 ==========\n")

# 计算每个模块的汇总指标（基于富集倍数和显著性）
module_summary <- enrichment_results %>%
  group_by(Module) %>%
  summarise(
    # 平均富集倍数
    Mean_Fold_Enrichment = mean(Fold_Enrichment, na.rm = TRUE),
    # 最大富集倍数
    Max_Fold_Enrichment = max(Fold_Enrichment, na.rm = TRUE),
    # 平均显著性（-log10 FDR）
    Mean_Significance = mean(-log10(FDR + 1e-10), na.rm = TRUE),
    # 最大显著性
    Max_Significance = max(-log10(FDR + 1e-10), na.rm = TRUE),
    # 显著富集的基因集数量
    Significant_GS_Count = sum(FDR < 0.05, na.rm = TRUE),
    # 总基因集数
    Total_GS_Count = n(),
    # 有富集（FE>1）的基因集数量
    Enriched_GS_Count = sum(Fold_Enrichment > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_Fold_Enrichment))

# 添加排名
module_summary$Rank <- 1:nrow(module_summary)

# 限制富集倍数显示范围（避免极端值影响）
module_summary$Max_Fold_Enrichment <- pmin(module_summary$Max_Fold_Enrichment, 10)
module_summary$Mean_Fold_Enrichment <- pmin(module_summary$Mean_Fold_Enrichment, 5)

# 保存汇总结果
write.csv(module_summary, output_summary, row.names = FALSE)
cat(sprintf("✓ 模块汇总已保存: %s\n", output_summary))

# 显示结果
cat("\n模块汇总统计（Top 10，按最大富集倍数排序）:\n")
print(head(module_summary[, c("Module", "Max_Fold_Enrichment", "Mean_Fold_Enrichment", 
                              "Mean_Significance", "Significant_GS_Count")], 10))

# ========================== 图1：最大富集倍数柱状图 ====================================
cat("\n========== 绘制最大富集倍数柱状图 ==========\n")

# 准备数据
plot_data <- module_summary %>%
  mutate(
    Module = factor(Module, levels = rev(Module)),
    # 标记是否有显著富集
    Has_Significant = Significant_GS_Count > 0
  )

p1 <- ggplot(plot_data, aes(x = Module, y = Max_Fold_Enrichment, fill = Has_Significant)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1fx", Max_Fold_Enrichment)), 
            hjust = -0.2, size = 3) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#D7191C", "FALSE" = "#2C7BB6"),
                    labels = c("TRUE" = "Has Significant Enrichment", 
                               "FALSE" = "No Significant Enrichment"),
                    name = "Significance") +
  labs(title = "Module Enrichment Summary",
       subtitle = "Maximum Fold Enrichment across all gene sets (Dashed line: FE = 1)",
       x = "Module", y = "Maximum Fold Enrichment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "top")

ggsave(output_barplot, p1, width = plot_width, height = plot_height, dpi = 300)
cat("✓ 最大富集倍数柱状图已保存:", output_barplot, "\n")

# ========================== 图2：气泡图（富集倍数 vs 显著性） ====================================
cat("\n========== 绘制气泡图 ==========\n")

# 准备气泡图数据
bubble_data <- module_summary %>%
  mutate(
    Module = factor(Module, levels = rev(Module)),
    # 限制气泡大小范围
    Bubble_Size = pmin(Significant_GS_Count, 20),
    # 标记是否显著
    Has_Significant = Significant_GS_Count > 0
  )

p2 <- ggplot(bubble_data, aes(x = Mean_Fold_Enrichment, y = Module)) +
  geom_point(aes(size = Bubble_Size, color = Mean_Significance), alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_color_gradient(low = "lightblue", high = "red", 
                       name = "Mean -log10(FDR)",
                       limits = c(0, max(bubble_data$Mean_Significance, na.rm = TRUE))) +
  scale_size_continuous(range = c(3, 12), 
                        name = "Significant Gene Sets Count",
                        breaks = pretty(bubble_data$Bubble_Size)) +
  labs(title = "Module Enrichment Summary",
       subtitle = "X-axis: Mean Fold Enrichment | Color: Mean Significance | Size: Significant GS Count",
       x = "Mean Fold Enrichment", y = "Module") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right") +
  xlim(0, max(bubble_data$Mean_Fold_Enrichment, na.rm = TRUE) * 1.1)

ggsave(output_bubble, p2, width = plot_width, height = plot_height, dpi = 300)
cat("✓ 气泡图已保存:", output_bubble, "\n")

# ========================== 图3：简化版排名图 ====================================
cat("\n========== 绘制简化版排名图 ==========\n")

# 创建综合得分（富集倍数 × 显著性）
module_summary$Combined_Score <- module_summary$Max_Fold_Enrichment * 
  module_summary$Max_Significance

# 按综合得分排序
top_modules <- module_summary %>%
  arrange(desc(Combined_Score)) %>%
  head(15)

p3 <- ggplot(top_modules, aes(x = reorder(Module, Combined_Score), y = Combined_Score)) +
  geom_bar(stat = "identity", aes(fill = Max_Significance), width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "red", name = "Max -log10(FDR)") +
  geom_text(aes(label = sprintf("FE=%.1fx", Max_Fold_Enrichment)), 
            hjust = -0.2, size = 3) +
  labs(title = "Top Modules by Combined Score",
       subtitle = "Combined Score = Max Fold Enrichment × Max Significance\nValues show Max Fold Enrichment",
       x = "Module", y = "Combined Score") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right")

rank_file <- file.path(gsea_dir, "module_enrichment_ranking.pdf")
ggsave(rank_file, p3, width = 10, height = 8, dpi = 300)
cat("✓ 简化版排名图已保存:", rank_file, "\n")

# ========================== 图4：热图（修复聚类问题） ====================================
cat("\n========== 绘制热图 ==========\n")

# 准备热图数据
heatmap_data <- module_summary %>%
  select(Module, Max_Fold_Enrichment, Mean_Significance, Significant_GS_Count) %>%
  mutate(
    Max_Fold_Enrichment_norm = pmin(Max_Fold_Enrichment, 5) / 5,
    Mean_Significance_norm = Mean_Significance / max(Mean_Significance, na.rm = TRUE),
    Significant_GS_Count_norm = Significant_GS_Count / max(Significant_GS_Count, na.rm = TRUE)
  )

heatmap_matrix <- as.matrix(heatmap_data[, c("Max_Fold_Enrichment_norm", 
                                             "Mean_Significance_norm", 
                                             "Significant_GS_Count_norm")])
rownames(heatmap_matrix) <- heatmap_data$Module
colnames(heatmap_matrix) <- c("Fold Enrichment", "Significance", "Sig GS Count")

# 检查数据是否有问题（修复NA处理）
cat("\n热图数据诊断:\n")
cat("  数据维度:", paste(dim(heatmap_matrix), collapse = " x "), "\n")
cat("  检查缺失值:", any(is.na(heatmap_matrix)), "\n")
cat("  检查无限值:", any(is.infinite(heatmap_matrix)), "\n")

# 安全地检查每个模块的数据变异情况（处理NA）
row_variances <- apply(heatmap_matrix, 1, function(x) {
  if(all(!is.na(x))) {
    return(var(x))
  } else {
    return(NA)
  }
})

cat("  行方差范围:", paste(range(row_variances, na.rm = TRUE), collapse = " ~ "), "\n")

# 检查是否存在方差为0的行（排除NA）
if(any(!is.na(row_variances) & row_variances == 0)) {
  zero_var_modules <- names(row_variances[!is.na(row_variances) & row_variances == 0])
  cat("  ⚠ 警告：以下模块所有指标值完全相同（方差为0）:\n")
  cat("    ", paste(zero_var_modules, collapse = ", "), "\n")
} else {
  cat("  ✓ 所有模块数据均有变异\n")
}

# 动态调整高度
heatmap_height <- max(6, nrow(heatmap_matrix) * 0.3)

# 方案1：基础heatmap（无聚类版本，确保不会出错）
heatmap_file <- file.path(gsea_dir, "module_enrichment_heatmap.pdf")
pdf(heatmap_file, width = 8, height = heatmap_height)

heatmap(heatmap_matrix,
        main = "Module Enrichment Summary Heatmap",
        col = colorRampPalette(heatmap_colors)(100),
        scale = "none",
        Rowv = NA,      # 禁用行聚类（避免NA/NaN/Inf错误）
        Colv = NA,      # 禁用列聚类
        cexRow = heatmap_font_row,
        cexCol = heatmap_font_col,
        margins = heatmap_margins)

dev.off()
cat("✓ 热图已保存（无聚类版本）:", heatmap_file, "\n")

# 方案2：带聚类的热图（使用pheatmap包，如果可用）
if (has_pheatmap) {
  cat("\n尝试绘制带聚类的热图...\n")
  
  heatmap_cluster_file <- file.path(gsea_dir, "module_enrichment_heatmap_clustered.pdf")
  
  # 为pheatmap准备数据（使用原始值，不归一化）
  pheatmap_matrix <- as.matrix(heatmap_data[, c("Max_Fold_Enrichment", 
                                                "Mean_Significance", 
                                                "Significant_GS_Count")])
  rownames(pheatmap_matrix) <- heatmap_data$Module
  colnames(pheatmap_matrix) <- c("Max Fold Enrichment", "Mean -log10(FDR)", "Significant GS Count")
  
  # 限制极端值
  pheatmap_matrix[, "Max Fold Enrichment"] <- pmin(pheatmap_matrix[, "Max Fold Enrichment"], 5)
  
  # 创建自定义颜色
  my_colors <- colorRampPalette(heatmap_colors)(100)
  
  # 检查是否存在方差为0的行并添加微小噪声（安全处理NA）
  has_zero_var <- any(!is.na(row_variances) & row_variances == 0)
  if(has_zero_var) {
    cat("  检测到方差为0的行，添加微小噪声以允许聚类...\n")
    # 为所有行添加微小噪声（避免NA）
    noise_matrix <- matrix(rnorm(nrow(pheatmap_matrix) * ncol(pheatmap_matrix), 0, 1e-6),
                           nrow = nrow(pheatmap_matrix), ncol = ncol(pheatmap_matrix))
    pheatmap_matrix <- pheatmap_matrix + noise_matrix
  }
  
  # 绘制带聚类的热图
  pheatmap(pheatmap_matrix,
           main = "Module Enrichment Summary Heatmap (Clustered)",
           color = my_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           display_numbers = FALSE,
           fontsize_row = 8,
           fontsize_col = 10,
           filename = heatmap_cluster_file,
           width = 8,
           height = heatmap_height)
  
  cat("✓ 带聚类的热图已保存:", heatmap_cluster_file, "\n")
  
} else {
  cat("\n⚠ 跳过聚类热图（pheatmap包未安装）\n")
  cat("  如需聚类热图，请运行: install.packages('pheatmap')\n")
}

# ========================== 可选：图5 - 标准化后的热图数据表 ====================================
cat("\n========== 保存热图数据 ==========\n")

# 保存热图使用的数据
heatmap_data_file <- file.path(gsea_dir, "heatmap_source_data.csv")
heatmap_export <- data.frame(
  Module = rownames(heatmap_matrix),
  Fold_Enrichment_norm = heatmap_matrix[, "Fold Enrichment"],
  Significance_norm = heatmap_matrix[, "Significance"],
  Sig_GS_Count_norm = heatmap_matrix[, "Sig GS Count"],
  # 同时保留原始值
  Fold_Enrichment_raw = module_summary$Max_Fold_Enrichment[match(rownames(heatmap_matrix), module_summary$Module)],
  Significance_raw = module_summary$Mean_Significance[match(rownames(heatmap_matrix), module_summary$Module)],
  Sig_GS_Count_raw = module_summary$Significant_GS_Count[match(rownames(heatmap_matrix), module_summary$Module)]
)
write.csv(heatmap_export, heatmap_data_file, row.names = FALSE)
cat("✓ 热图源数据已保存:", heatmap_data_file, "\n")

# ========================== 生成报告 ====================================
cat("\n========== 生成分析报告 ==========\n")

sink(file.path(gsea_dir, "module_summary_report.txt"))

cat("================================================================================\n")
cat("                模块-基因集富集汇总报告（富集倍数版）\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、Top 10 模块（按最大富集倍数排序）\n")
cat("--------------------------------------------------------------------------------\n")
top10 <- head(module_summary, 10)
for (i in 1:nrow(top10)) {
  cat(sprintf("  %d. %s: 最大FE=%.2fx, 平均FE=%.2fx, 平均显著性=%.2f, 显著基因集数=%d/%d\n", 
              i, top10$Module[i], 
              top10$Max_Fold_Enrichment[i],
              top10$Mean_Fold_Enrichment[i],
              top10$Mean_Significance[i],
              top10$Significant_GS_Count[i],
              top10$Total_GS_Count[i]))
}

cat("\n二、富集倍数解读\n")
cat("--------------------------------------------------------------------------------\n")
cat("  FE > 2: 强富集\n")
cat("  1.5 < FE < 2: 中等富集\n")
cat("  1 < FE < 1.5: 弱富集\n")
cat("  FE = 1: 随机水平\n")
cat("  FE < 1: 耗竭\n")

cat("\n三、输出文件\n")
cat("--------------------------------------------------------------------------------\n")
cat("1.", basename(output_summary), "- 模块汇总统计表\n")
cat("2.", basename(output_barplot), "- 最大富集倍数柱状图\n")
cat("3.", basename(output_bubble), "- 富集倍数 vs 显著性气泡图\n")
cat("4.", basename(rank_file), "- 综合得分排名图\n")
cat("5.", basename(heatmap_file), "- 热图（无聚类）\n")
if (has_pheatmap) {
  cat("6.", "module_enrichment_heatmap_clustered.pdf", "- 热图（带聚类）\n")
  cat("7.", basename(heatmap_data_file), "- 热图源数据\n")
} else {
  cat("6.", basename(heatmap_data_file), "- 热图源数据\n")
}

sink()

cat("\n========================================\n")
cat("模块汇总可视化完成！\n")
cat("输出目录:", gsea_dir, "\n")
cat("========================================\n")