# ============================================================================
# WGCNA分析 - 模块4：模块特征基因差异分析
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
me_dir <- file.path(output_dir, "module_eigengene_analysis")
if (!dir.exists(me_dir)) {
  dir.create(me_dir, recursive = TRUE)
}

# ========================== 参数设置 =========================================
# 可视化参数
boxplot_width <- 8
boxplot_height <- 6

# 输出文件
output_rdata <- file.path(me_dir, "module_eigengene_analysis.RData")
output_boxplot <- file.path(me_dir, "module_eigengene_boxplot.pdf")
output_heatmap <- file.path(me_dir, "module_eigengene_heatmap.pdf")
output_stats <- file.path(me_dir, "module_eigengene_statistics.csv")

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)

# ========================== 加载数据 ====================================
cat("\n========== 加载数据 ==========\n")

load(file.path(output_dir, "02_network_modules.RData"))
load(file.path(output_dir, "01_data_prepared.RData"))
load(file.path(output_dir, "03_module_enrichment.RData"))

cat(sprintf("Top模块: %s\n", paste(top_modules, collapse = ", ")))

# ========================== 提取Top模块的特征基因 ====================================
cat("\n========== 提取Top模块的特征基因 ==========\n")

# 只保留Top模块的特征基因
top_MEs <- MEs[, colnames(MEs) %in% paste0("ME", top_modules), drop = FALSE]

# 添加样品分组信息
me_data <- data.frame(
  Sample = rownames(top_MEs),
  Group = sample_groups$Group[match(rownames(top_MEs), sample_groups$SampleID)],
  top_MEs,
  stringsAsFactors = FALSE
)

cat("特征基因数据:\n")
print(head(me_data))

# ========================== 计算各组统计量 ====================================
cat("\n========== 计算各组统计量 ==========\n")

me_stats <- data.frame()

for (mod in colnames(top_MEs)) {
  for (grp in unique(sample_groups$Group)) {
    grp_values <- me_data[me_data$Group == grp, mod]
    
    me_stats <- rbind(me_stats, data.frame(
      Module = mod,
      Group = grp,
      Mean = mean(grp_values, na.rm = TRUE),
      SD = sd(grp_values, na.rm = TRUE),
      Median = median(grp_values, na.rm = TRUE),
      Min = min(grp_values, na.rm = TRUE),
      Max = max(grp_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(me_stats, output_stats, row.names = FALSE)
cat(sprintf("特征基因统计已保存: %s\n", output_stats))

# ========================== 统计检验 ====================================
cat("\n========== 统计检验 ==========\n")

stat_results <- data.frame()

for (mod in colnames(top_MEs)) {
  # 提取数据
  control_vals <- me_data[me_data$Group == "Control", mod]
  treatment_groups <- setdiff(unique(me_data$Group), "Control")
  
  for (trt in treatment_groups) {
    treatment_vals <- me_data[me_data$Group == trt, mod]
    
    # Wilcoxon检验
    if (length(control_vals) >= 2 && length(treatment_vals) >= 2) {
      test_result <- wilcox.test(control_vals, treatment_vals)
      p_val <- test_result$p.value
      fc <- mean(treatment_vals, na.rm = TRUE) / mean(control_vals, na.rm = TRUE)
      
      stat_results <- rbind(stat_results, data.frame(
        Module = mod,
        Comparison = paste(trt, "vs Control"),
        Fold_Change = fc,
        P_value = p_val,
        Significant = p_val < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
}

stat_results$Signif_mark <- ifelse(stat_results$P_value < 0.001, "***",
                                   ifelse(stat_results$P_value < 0.01, "**",
                                          ifelse(stat_results$P_value < 0.05, "*", "ns")))

cat("\n统计检验结果:\n")
print(stat_results)

# ========================== 图1：特征基因箱线图 ====================================
cat("\n========== 绘制特征基因箱线图 ==========\n")

# 转换数据格式
plot_data <- melt(me_data, id.vars = c("Sample", "Group"), 
                  variable.name = "Module", value.name = "Eigengene")

# 只保留Top模块
plot_data <- plot_data[plot_data$Module %in% colnames(top_MEs), ]

# 设置因子顺序
plot_data$Group <- factor(plot_data$Group, levels = unique(sample_groups$Group))

p1 <- ggplot(plot_data, aes(x = Group, y = Eigengene, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  facet_wrap(~Module, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Control" = "#2E86AB", 
                               "DMEM" = "#1F8A4C",
                               "I6" = "#C60036", 
                               "R2" = "#E4945A")) +
  labs(title = "Module Eigengene Expression Across Groups",
       x = "Group", y = "Module Eigengene") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none")

ggsave(output_boxplot, p1, width = boxplot_width, height = boxplot_height, dpi = 300)
cat("✓ 特征基因箱线图已保存\n")

# ========================== 图2：特征基因热图 ====================================
cat("\n========== 绘制特征基因热图 ==========\n")

# 准备热图数据
heatmap_data <- top_MEs
colnames(heatmap_data) <- gsub("ME", "", colnames(heatmap_data))

# Z-score标准化
heatmap_scaled <- t(scale(t(heatmap_data)))
heatmap_scaled[heatmap_scaled > 2] <- 2
heatmap_scaled[heatmap_scaled < -2] <- -2

# 创建注释
annotation_col <- data.frame(
  Group = sample_groups$Group,
  row.names = sample_groups$SampleID
)

annotation_colors <- list(
  Group = c("Control" = "#2E86AB", "DMEM" = "#1F8A4C", 
            "I6" = "#C60036", "R2" = "#E4945A")
)

pdf(output_heatmap, width = 8, height = 6)
pheatmap(t(heatmap_scaled),
         color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Module Eigengene Heatmap",
         fontsize_row = 10,
         fontsize_col = 8)
dev.off()
cat("✓ 特征基因热图已保存\n")

# ========================== 保存结果 ====================================
cat("\n========== 保存结果 ==========\n")

save(me_data, me_stats, stat_results, top_MEs, file = output_rdata)

cat(sprintf("特征基因分析结果已保存: %s\n", output_rdata))

cat("\n模块4完成！\n")