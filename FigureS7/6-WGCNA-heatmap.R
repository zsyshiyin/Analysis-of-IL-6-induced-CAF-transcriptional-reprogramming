# ============================================================================
# WGCNA分析 - 模块6：模块-分组关系热图（经典WGCNA热图）
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
module_trait_dir <- file.path(output_dir, "module_trait_analysis")
if (!dir.exists(module_trait_dir)) {
  dir.create(module_trait_dir, recursive = TRUE)
}

# ========================== 参数设置 =========================================
# 可视化参数
heatmap_width <- 10
heatmap_height <- 8

# 输出文件
output_rdata <- file.path(module_trait_dir, "module_trait_analysis.RData")
output_heatmap <- file.path(module_trait_dir, "module_trait_heatmap.pdf")
output_correlation <- file.path(module_trait_dir, "module_trait_correlation.csv")
output_summary <- file.path(module_trait_dir, "module_trait_summary.txt")

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

library(WGCNA)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(reshape2)

enableWGCNAThreads()

cat("✓ 所有包已加载\n")

# ========================== 加载数据 ====================================
cat("\n========== 加载数据 ==========\n")

load(file.path(output_dir, "02_network_modules.RData"))
load(file.path(output_dir, "01_data_prepared.RData"))

cat(sprintf("模块数: %d\n", ncol(MEs)))
cat(sprintf("样品数: %d\n", nrow(MEs)))

# ========================== 查看实际分组 ====================================
cat("\n========== 实际分组信息 ==========\n")

# 获取实际存在的分组
all_groups <- unique(sample_groups$Group)
cat("实际分组:", paste(all_groups, collapse = ", "), "\n")
cat("分组数量:", length(all_groups), "\n")

# ========================== 创建分组设计矩阵（只包含实际分组） ====================================
cat("\n========== 创建设计矩阵 ==========\n")

# 只创建实际分组的二值变量
trait_matrix <- data.frame(row.names = sample_groups$SampleID)

for (grp in all_groups) {
  trait_matrix[[grp]] <- as.numeric(sample_groups$Group == grp)
}

# 可选：添加一个数值型变量（处理组 vs 对照组）- 只在有Control组时添加
#if ("Control" %in% all_groups) {
#trait_matrix$Treatment <- as.numeric(sample_groups$Group != "Control")
#cat("添加了Treatment变量（处理组 vs 对照组）\n")
#}

cat(sprintf("设计矩阵维度: %d 样品 × %d 变量\n", nrow(trait_matrix), ncol(trait_matrix)))
cat("设计矩阵列名:", paste(colnames(trait_matrix), collapse = ", "), "\n")

# ========================== 计算模块与分组的相关性 ====================================
cat("\n========== 计算模块与分组的相关性 ==========\n")

# 计算相关性
module_trait_cor <- cor(MEs, trait_matrix, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(MEs))

# 整理结果
module_trait_results <- data.frame()

for (mod in rownames(module_trait_cor)) {
  for (trait in colnames(module_trait_cor)) {
    module_trait_results <- rbind(module_trait_results, data.frame(
      Module = gsub("ME", "", mod),
      Trait = trait,
      Correlation = module_trait_cor[mod, trait],
      PValue = module_trait_pvalue[mod, trait],
      stringsAsFactors = FALSE
    ))
  }
}

# 添加显著性标记
module_trait_results$Significant <- module_trait_results$PValue < 0.05
module_trait_results$Signif_mark <- ifelse(module_trait_results$PValue < 0.001, "***",
                                           ifelse(module_trait_results$PValue < 0.01, "**",
                                                  ifelse(module_trait_results$PValue < 0.05, "*", "")))

# 保存结果
write.csv(module_trait_results, output_correlation, row.names = FALSE)
cat(sprintf("✓ 模块-分组相关性已保存: %s\n", output_correlation))

# 显示结果
cat("\n模块-分组相关性（Top 10）:\n")
print(head(module_trait_results[order(-abs(module_trait_results$Correlation)), ], 10))

# ========================== 图1：模块-分组关系热图（经典WGCNA热图） ====================================
cat("\n========== 绘制模块-分组关系热图 ==========\n")

# 准备热图数据
heatmap_data <- module_trait_cor

# 简化模块名称（去掉ME前缀）
rownames(heatmap_data) <- gsub("ME", "", rownames(heatmap_data))

# 创建显著性标签矩阵
signif_labels <- matrix("", nrow = nrow(heatmap_data), ncol = ncol(heatmap_data))
for (i in 1:nrow(module_trait_pvalue)) {
  for (j in 1:ncol(module_trait_pvalue)) {
    if (module_trait_pvalue[i, j] < 0.05) {
      if (module_trait_pvalue[i, j] < 0.001) {
        signif_labels[i, j] <- "***"
      } else if (module_trait_pvalue[i, j] < 0.01) {
        signif_labels[i, j] <- "**"
      } else {
        signif_labels[i, j] <- "*"
      }
    }
  }
}

# 动态调整热图尺寸
n_modules <- nrow(heatmap_data)
n_traits <- ncol(heatmap_data)
heatmap_height <- max(6, n_modules * 0.3)
heatmap_width <- max(8, n_traits * 1.5)

# 绘制热图
pdf(output_heatmap, width = heatmap_width, height = heatmap_height)

pheatmap(heatmap_data,
         color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         display_numbers = signif_labels,
         number_color = "black",
         number_size = 8,
         main = "Module-Trait Relationships",
         fontsize_row = 10,
         fontsize_col = 12)

dev.off()
cat("✓ 模块-分组关系热图已保存:", output_heatmap, "\n")

# ========================== 图2：模块特征基因箱线图 ====================================
cat("\n========== 绘制模块特征基因箱线图 ==========\n")

# 准备数据
me_plot_data <- data.frame(
  Sample = rownames(MEs),
  Group = sample_groups$Group,
  MEs,
  stringsAsFactors = FALSE
)

# 转换为长格式
me_plot_long <- melt(me_plot_data, 
                     id.vars = c("Sample", "Group"),
                     variable.name = "Module",
                     value.name = "Eigengene")

# 只保留非灰色模块
me_plot_long <- me_plot_long[!grepl("grey", me_plot_long$Module), ]

if (nrow(me_plot_long) > 0) {
  # 选择Top模块（基于与Treatment的相关性）
  if ("Treatment" %in% colnames(module_trait_cor)) {
    top_modules <- module_trait_results %>%
      filter(Trait == "Treatment") %>%
      arrange(desc(abs(Correlation))) %>%
      head(9) %>%
      pull(Module)
  } else if (length(all_groups) > 1) {
    # 如果没有Treatment，使用第一个非Control组
    other_group <- all_groups[all_groups != "Control"][1]
    top_modules <- module_trait_results %>%
      filter(Trait == other_group) %>%
      arrange(desc(abs(Correlation))) %>%
      head(9) %>%
      pull(Module)
  } else {
    # 如果只有一个组，按模块大小选择
    module_sizes <- table(moduleColors)
    top_modules <- names(sort(module_sizes, decreasing = TRUE))[1:9]
    top_modules <- top_modules[top_modules != "grey"]
  }
  
  # 过滤数据
  me_plot_filtered <- me_plot_long[me_plot_long$Module %in% paste0("ME", top_modules), ]
  me_plot_filtered$Module <- gsub("ME", "", me_plot_filtered$Module)
  
  # 绘制箱线图
  n_modules_plot <- length(unique(me_plot_filtered$Module))
  n_rows <- ceiling(n_modules_plot / 3)
  plot_height <- max(6, n_rows * 3)
  
  p <- ggplot(me_plot_filtered, aes(x = Group, y = Eigengene, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    facet_wrap(~Module, scales = "free_y", nrow = n_rows) +
    scale_fill_manual(values = setNames(
      c("#2E86AB", "#1F8A4C", "#C60036", "#E4945A")[1:length(all_groups)],
      all_groups)) +
    labs(title = "Module Eigengene Expression by Group",
         x = "Group", y = "Module Eigengene") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 10, face = "bold"),
          legend.position = "none")
  
  boxplot_file <- file.path(module_trait_dir, "module_eigengene_boxplot.pdf")
  ggsave(boxplot_file, p, width = 12, height = plot_height, dpi = 300)
  cat("✓ 模块特征基因箱线图已保存:", boxplot_file, "\n")
}

# ========================== 图3：Top模块的特征基因表达热图 ====================================
cat("\n========== 绘制Top模块的特征基因表达热图 ==========\n")

# 选择与分组显著相关的模块
if ("Treatment" %in% colnames(module_trait_cor)) {
  significant_modules <- module_trait_results %>%
    filter(Significant == TRUE, Trait == "Treatment") %>%
    arrange(desc(abs(Correlation))) %>%
    head(10) %>%
    pull(Module)
} else if (length(all_groups) > 1) {
  other_group <- all_groups[all_groups != "Control"][1]
  significant_modules <- module_trait_results %>%
    filter(Significant == TRUE, Trait == other_group) %>%
    arrange(desc(abs(Correlation))) %>%
    head(10) %>%
    pull(Module)
} else {
  significant_modules <- c()
}

if (length(significant_modules) == 0) {
  # 如果没有显著模块，选择相关性最高的
  if ("Treatment" %in% colnames(module_trait_cor)) {
    significant_modules <- module_trait_results %>%
      filter(Trait == "Treatment") %>%
      arrange(desc(abs(Correlation))) %>%
      head(5) %>%
      pull(Module)
  } else {
    significant_modules <- module_trait_results %>%
      arrange(desc(abs(Correlation))) %>%
      head(5) %>%
      pull(Module)
  }
}

cat("Top模块:", paste(significant_modules, collapse = ", "), "\n")

if (length(significant_modules) > 0) {
  # 提取这些模块的特征基因
  top_MEs <- MEs[, paste0("ME", significant_modules), drop = FALSE]
  
  # Z-score标准化
  top_MEs_scaled <- t(scale(t(top_MEs)))
  top_MEs_scaled[is.na(top_MEs_scaled)] <- 0
  top_MEs_scaled[top_MEs_scaled > 2] <- 2
  top_MEs_scaled[top_MEs_scaled < -2] <- -2
  
  # 创建注释
  annotation_col <- data.frame(
    Group = sample_groups$Group,
    row.names = sample_groups$SampleID
  )
  
  # 只设置实际存在的分组颜色
  group_colors_used <- c()
  for (grp in all_groups) {
    if (grp == "Control") group_colors_used[grp] <- "#2E86AB"
    else if (grp == "DMEM") group_colors_used[grp] <- "#1F8A4C"
    else if (grp == "I6") group_colors_used[grp] <- "#C60036"
    else if (grp == "R2") group_colors_used[grp] <- "#E4945A"
    else group_colors_used[grp] <- "gray"
  }
  
  annotation_colors <- list(Group = group_colors_used)
  
  # 绘制热图
  heatmap_file <- file.path(module_trait_dir, "top_modules_eigengene_heatmap.pdf")
  pdf(heatmap_file, width = 10, height = max(6, length(significant_modules) * 0.4))
  pheatmap(t(top_MEs_scaled),
           color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           main = "Top Module Eigengenes",
           fontsize_row = 10,
           fontsize_col = 8)
  dev.off()
  cat("✓ Top模块特征基因热图已保存:", heatmap_file, "\n")
}

# ========================== 生成报告 ====================================
cat("\n========== 生成分析报告 ==========\n")

sink(output_summary)

cat("================================================================================\n")
cat("                   模块-分组关系分析报告\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、数据概况\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("实际分组: %s\n", paste(all_groups, collapse = ", ")))
cat(sprintf("分组数量: %d\n", length(all_groups)))
cat(sprintf("模块数: %d\n", ncol(MEs)))
cat(sprintf("样品数: %d\n", nrow(MEs)))

cat("\n二、模块-分组相关性（Top 10）\n")
cat("--------------------------------------------------------------------------------\n")
top_cors <- module_trait_results %>%
  arrange(desc(abs(Correlation))) %>%
  head(10)
for (i in 1:nrow(top_cors)) {
  cat(sprintf("  %s vs %s: r = %.3f, p = %.2e %s\n", 
              top_cors$Module[i], top_cors$Trait[i], 
              top_cors$Correlation[i], top_cors$PValue[i], 
              top_cors$Signif_mark[i]))
}

cat("\n三、输出文件\n")
cat("--------------------------------------------------------------------------------\n")
cat("1.", basename(output_correlation), "- 模块-分组相关性表格\n")
cat("2.", basename(output_heatmap), "- 模块-分组关系热图\n")
if (exists("boxplot_file")) {
  cat("3.", basename(boxplot_file), "- 模块特征基因箱线图\n")
}
if (exists("heatmap_file")) {
  cat("4.", basename(heatmap_file), "- Top模块特征基因热图\n")
}

sink()

cat("\n========================================\n")
cat("模块-分组关系分析完成！\n")
cat("输出目录:", module_trait_dir, "\n")
cat("========================================\n")