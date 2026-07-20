# ============================================================================
# WGCNA分析 - 模块2：网络构建和模块识别（小样本优化版）
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

# ========================== 参数设置 =========================================
# 网络构建参数（针对小样本优化）
power <- NULL                    # 软阈值幂次（NULL表示自动选择）
min_module_size <- 5             # 最小模块大小（小样本适当降低）
merge_cut_height <- 0.3          # 模块合并阈值
deep_split <- 2                  # 模块检测深度

# 输出文件
output_rdata <- file.path(output_dir, "02_network_modules.RData")
output_power_plot <- file.path(output_dir, "power_selection.pdf")
output_dendrogram <- file.path(output_dir, "gene_dendrogram.pdf")
output_module_colors <- file.path(output_dir, "module_colors.csv")

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

library(WGCNA)
library(dplyr)

enableWGCNAThreads()

# ========================== 加载预处理数据 ====================================
cat("\n========== 加载预处理数据 ==========\n")

load(file.path(output_dir, "01_data_prepared.RData"))
cat(sprintf("加载数据: %d 基因 × %d 样品\n", nrow(expr_matrix_norm), ncol(expr_matrix_norm)))

# ========================== 选择软阈值 ====================================
cat("\n========== 选择软阈值 ==========\n")

# 计算软阈值
powers <- c(1:20)
sft <- pickSoftThreshold(
  t(expr_matrix_norm),
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

# 绘制软阈值选择图
pdf(output_power_plot, width = 12, height = 6)
par(mfrow = c(1, 2))

# 左图：拟合指数
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.7, col = "red")  # 小样本降低阈值到0.7

# 右图：平均连接度
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", 
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

dev.off()
cat("✓ 软阈值选择图已保存\n")

# 选择最佳软阈值（第一个R^2 > 0.7的power）
if (is.null(power)) {
  valid_powers <- which(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] > 0.7)
  if (length(valid_powers) > 0) {
    power <- sft$fitIndices[valid_powers[1], 1]
  } else {
    power <- 4  # 小样本默认值
    cat("警告：未找到R^2 > 0.7的软阈值，使用默认值 power = 4\n")
  }
}
cat(sprintf("选择的软阈值: power = %d\n", power))

# ========================== 构建网络 ====================================
cat("\n========== 构建网络 ==========\n")

# 计算邻接矩阵
adjacency <- adjacency(t(expr_matrix_norm), power = power, type = "signed")

# 转换为TOM矩阵
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# 基因聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 绘制基因聚类树
pdf(output_dendrogram, width = 12, height = 8)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering dendrogram", 
     labels = FALSE, hang = 0.04)
dev.off()
cat("✓ 基因聚类树已保存\n")

# ========================== 模块识别 ====================================
cat("\n========== 模块识别 ==========\n")

# 动态剪枝识别模块
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deep_split,
  pamRespectsDendro = FALSE,
  minClusterSize = min_module_size
)

# 获取模块颜色
moduleColors <- labels2colors(dynamicMods)
cat("模块统计:\n")
print(table(moduleColors))

# 计算模块特征基因
MEs <- moduleEigengenes(t(expr_matrix_norm), colors = moduleColors)$eigengenes

# 合并相似模块
merge <- mergeCloseModules(t(expr_matrix_norm), moduleColors, cutHeight = merge_cut_height)
moduleColors <- merge$colors
MEs <- merge$newMEs

cat("\n合并后模块统计:\n")
module_counts <- table(moduleColors)
print(module_counts)

# 保存模块颜色
module_colors_df <- data.frame(
  Gene = rownames(expr_matrix_norm),
  Module = moduleColors,
  stringsAsFactors = FALSE
)
write.csv(module_colors_df, output_module_colors, row.names = FALSE)

# ========================== 保存结果 ====================================
cat("\n========== 保存结果 ==========\n")

save(adjacency, TOM, dissTOM, geneTree, dynamicMods, moduleColors, MEs,
     power, file = output_rdata)

cat(sprintf("网络和模块结果已保存: %s\n", output_rdata))
cat(sprintf("  模块数: %d\n", length(unique(moduleColors[moduleColors != "grey"]))))
cat(sprintf("  灰色模块（未分配）基因数: %d\n", sum(moduleColors == "grey")))

cat("\n模块2完成！\n")