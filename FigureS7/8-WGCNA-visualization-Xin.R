# ============================================================================
# WGCNA可视化 - 完整版（可自定义所有参数）
# 包括：样本聚类热图、模块划分、模块-表型关系、MM-GS筛选、富集气泡图
# ============================================================================

# 清空环境变量
rm(list = ls())

# ========================== 工作路径设定 ====================================
work_dir <- "D:/zsy/SX/Fomal-final/3-CAF-WGCNA"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")

# ========================== 创建输出目录 ====================================
output_dir <- "WGCNA_visualization_Xin_Xin"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# ============================================================================
# ====================== 用户自定义参数区域 ==================================
# ============================================================================

# ----------------------------- 全局参数 -------------------------------------
global_dpi <- 300

# ----------------------------- 图1：样本聚类树参数 ---------------------------
sample_tree_width <- 14
sample_tree_height <- 10
sample_tree_main <- "Sample dendrogram and trait heatmap"
sample_tree_cex_dendroLabels <- 0.8
sample_tree_cex_colorLabels <- 1.2
sample_tree_cex_rowText <- 0.8

# ----------------------------- 图2：基因聚类树参数 ---------------------------
gene_dendro_width <- 10
gene_dendro_height <- 8
gene_dendro_main <- "Gene dendrogram and module colors"
gene_dendro_hang <- 0.03
gene_dendro_guideHang <- 0.05

# ----------------------------- 图3：模块-性状热图参数 -------------------------
module_trait_width <- 6
module_trait_height <- 8
module_trait_main <- "Module-trait relationships"
module_trait_cex_text <- 1.2
module_trait_zlim <- c(-1, 1)
module_trait_colors <- c("blue", "white", "red")  # 低、中、高颜色

# p值显示格式
module_trait_show_pvalue <- TRUE  # 是否显示p值
module_trait_cor_digits <- 2      # 相关性小数位数
module_trait_p_digits <- 3        # p值小数位数

# ----------------------------- 图4：模块富集气泡图参数 -------------------------
bubble_width <- 8
bubble_height <- 6
bubble_title <- "Module Enrichment Summary"
bubble_title_size <- 22
bubble_x_label <- "Maximum Fold Enrichment"
bubble_y_label <- "Module"
bubble_axis_title_size <- 20
bubble_axis_text_size <- 18
bubble_point_alpha <- 0.8
bubble_color_low <- "#25377F"
bubble_color_high <- "#C60036"
bubble_size_min <- 3
bubble_size_max <- 12

# ----------------------------- 图5：模块特征基因树和热图参数 -------------------
eigengene_width <- 10
eigengene_height <- 12
eigengene_marDendro <- c(0, 4, 1, 2)
eigengene_marHeatmap <- c(3, 4, 1, 2)
eigengene_cex_lab <- 0.8
eigengene_xLabelsAngle <- 90

# ----------------------------- 图6：TOM plot参数 -------------------------------
tom_plot_width <- 12
tom_plot_height <- 10
tom_plot_main <- "Network heatmap plot"
tom_plot_nSelect <- 400                    # 显示的基因数
tom_plot_power <- 7                        # TOM幂次
tom_plot_colors <- c("red", "orange", "lemonchiffon")  # 热图颜色
tom_plot_seed <- 10                        # 随机种子

# ============================================================================
# ====================== 加载必要的包 ========================================
# ============================================================================

cat("\n========== 加载必要的包 ==========\n")

library(WGCNA)
library(ggplot2)
library(dplyr)
library(gplots)

options(allowWGCNAThreads = FALSE)

cat("✓ 所有包已加载\n")

# ========================== 加载WGCNA结果 ====================================
cat("\n========== 加载WGCNA结果 ==========\n")

load(file.path("WGCNA_analysis", "02_network_modules.RData"))
load(file.path("WGCNA_analysis", "01_data_prepared.RData"))

# 加载富集结果
enrichment_file <- file.path("WGCNA_analysis", "gene_set_enrichment", "module_gene_set_enrichment.csv")
if (file.exists(enrichment_file)) {
  enrichment_results <- read.csv(enrichment_file, stringsAsFactors = FALSE)
  cat("富集结果加载成功\n")
} else {
  cat("未找到富集结果文件\n")
  enrichment_results <- NULL
}

cat("数据加载成功\n")
cat(sprintf("模块数: %d\n", length(unique(moduleColors))))
cat(sprintf("基因数: %d\n", nrow(expr_matrix_norm)))
cat(sprintf("样品数: %d\n", ncol(expr_matrix_norm)))

# ========================== 获取分组信息 ====================================
cat("\n========== 获取分组信息 ==========\n")

sample_names <- colnames(expr_matrix_norm)
group_mapping <- list("Con" = "Control", "DMEM" = "DMEM", "I6" = "I6", "R2" = "R2")

sample_groups <- data.frame(SampleID = sample_names, Group = NA, stringsAsFactors = FALSE)
for (prefix in names(group_mapping)) {
  matching_samples <- grep(paste0("^", prefix), sample_names, value = TRUE)
  if (length(matching_samples) > 0) {
    sample_groups$Group[sample_groups$SampleID %in% matching_samples] <- group_mapping[[prefix]]
  }
}
sample_groups <- sample_groups[!is.na(sample_groups$Group), ]
groups <- unique(sample_groups$Group)

cat("分组信息:\n")
print(table(sample_groups$Group))

# ========================== 创建设计矩阵（性状矩阵） ====================================
datTraits <- data.frame(row.names = sample_groups$SampleID)
for (grp in groups) {
  datTraits[[grp]] <- as.numeric(sample_groups$Group == grp)
}

# 确保样品顺序一致
expr_samples <- colnames(expr_matrix_norm)
traitSamples <- rownames(datTraits)
traitRows <- match(expr_samples, traitSamples)
datTraits <- datTraits[traitRows, , drop = FALSE]

cat("性状矩阵维度:", dim(datTraits), "\n")

# ========================== 计算模块特征基因 ====================================
MEs <- moduleEigengenes(t(expr_matrix_norm), colors = moduleColors)$eigengenes

# ========================== 图1：样本聚类树 + 性状热图 ====================================
cat("\n========== 图1：样本聚类树与性状热图 ==========\n")

sampleTree2 <- hclust(dist(t(expr_matrix_norm)), method = "average")

# 创建性状颜色矩阵
traitColors <- data.frame(matrix(0, nrow = nrow(datTraits), ncol = ncol(datTraits)))
for (i in 1:ncol(datTraits)) {
  traitColors[, i] <- as.numeric(datTraits[, i])
}
colnames(traitColors) <- colnames(datTraits)

pdf(file.path(output_dir, "1_Sample_dendrogram_and_trait_heatmap.pdf"), 
    width = sample_tree_width, height = sample_tree_height)

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = sample_tree_main,
                    cex.colorLabels = sample_tree_cex_colorLabels,
                    cex.dendroLabels = sample_tree_cex_dendroLabels,
                    cex.rowText = sample_tree_cex_rowText)

dev.off()
cat("✓ 图1已保存\n")

# ========================== 图2：基因聚类树与模块颜色 ====================================
cat("\n========== 图2：基因聚类树与模块颜色 ==========\n")

pdf(file.path(output_dir, "2_Gene_dendrogram_and_module_colors.pdf"), 
    width = gene_dendro_width, height = gene_dendro_height)

plotDendroAndColors(geneTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = gene_dendro_hang,
                    addGuide = TRUE,
                    guideHang = gene_dendro_guideHang,
                    main = gene_dendro_main)

dev.off()
cat("✓ 图2已保存\n")

# ========================== 图3：模块-性状关系热图 ====================================
cat("\n========== 图3：模块-性状关系热图 ==========\n")

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

# 创建显示标签
if (module_trait_show_pvalue) {
  textMatrix <- paste(signif(moduleTraitCor, module_trait_cor_digits), "\n(",
                      signif(moduleTraitPvalue, module_trait_p_digits), ")", sep = "")
} else {
  textMatrix <- paste(signif(moduleTraitCor, module_trait_cor_digits), sep = "")
}
dim(textMatrix) <- dim(moduleTraitCor)

# 简化模块名称
rownames(moduleTraitCor) <- gsub("ME", "", rownames(moduleTraitCor))

# ========== 自定义热图颜色 ==========
# 方式1：使用预定义颜色名称
# module_trait_colors <- c("blue", "white", "red")  # 蓝-白-红渐变

# 方式2：使用十六进制颜色代码
module_trait_colors <- c("#4A60B5", "white", "#E04A6E")  # 深蓝-白-深红

# 方式3：使用自定义颜色向量
# module_trait_colors <- c("purple", "black", "yellow")  # 紫-黑-黄

# 根据用户设置的颜色创建渐变
if (length(module_trait_colors) == 3) {
  # 三色渐变
  heatmap_colors <- colorRampPalette(module_trait_colors)(50)
} else if (length(module_trait_colors) == 2) {
  # 两色渐变
  heatmap_colors <- colorRampPalette(module_trait_colors)(50)
} else {
  # 默认使用红-白-蓝渐变
  heatmap_colors <- blueWhiteRed(50)
}

pdf(file.path(output_dir, "3_Module_trait_relationships.pdf"), 
    width = module_trait_width, height = module_trait_height)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = heatmap_colors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = module_trait_cex_text,
               zlim = module_trait_zlim,
               main = module_trait_main)

dev.off()
cat("✓ 图3已保存\n")

# ========================== 图4：模块富集气泡图 ====================================
cat("\n========== 图4：模块富集气泡图 ==========\n")

if (!is.null(enrichment_results) && nrow(enrichment_results) > 0) {
  module_enrichment_summary <- enrichment_results %>%
    group_by(Module) %>%
    summarise(
      Max_Fold_Enrichment = max(Fold_Enrichment, na.rm = TRUE),
      Min_FDR = min(FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      neg_log10_FDR = -log10(Min_FDR + 1e-10),
      neg_log10_FDR = pmin(neg_log10_FDR, 10)
    )
  
  # 确保所有模块都显示
  all_modules <- unique(moduleColors[moduleColors != "grey"])
  module_enrichment_summary <- module_enrichment_summary %>%
    right_join(data.frame(Module = all_modules), by = "Module") %>%
    mutate(
      Max_Fold_Enrichment = ifelse(is.na(Max_Fold_Enrichment), 1, Max_Fold_Enrichment),
      neg_log10_FDR = ifelse(is.na(neg_log10_FDR), 0, neg_log10_FDR)
    )
  
  # 修正：颜色代表FDR（-log10(FDR)），大小代表富集分数（Fold Enrichment）
  p4 <- ggplot(module_enrichment_summary, aes(x = Max_Fold_Enrichment, y = Module)) +
    geom_point(aes(size = Max_Fold_Enrichment, color = neg_log10_FDR), alpha = bubble_point_alpha) +
    scale_color_gradient(low = bubble_color_low, high = bubble_color_high, name = "-log10(FDR)") +
    scale_size_continuous(range = c(bubble_size_min, bubble_size_max), name = "Fold Enrichment") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#C60036", alpha = 0.5) +
    labs(title = bubble_title, x = bubble_x_label, y = bubble_y_label) +
    theme_bw() +
    theme(plot.title = element_text(size = bubble_title_size, face = "bold", hjust = 0.5),
          axis.title = element_text(size = bubble_axis_title_size),
          axis.text = element_text(size = bubble_axis_text_size))
  
  ggsave(file.path(output_dir, "4_Module_enrichment_bubble.pdf"), 
         p4, width = bubble_width, height = bubble_height, dpi = global_dpi)
  cat("✓ 图4已保存\n")
} else {
  cat("未找到富集结果，跳过图4\n")
}

# ========================== 图5：模块特征基因树和热图 ====================================
cat("\n========== 图5：模块特征基因树和热图 ==========\n")

pdf(file.path(output_dir, "5_Eigengene_dendrogram_and_heatmap.pdf"), 
    width = eigengene_width, height = eigengene_height)

plotEigengeneNetworks(MEs, "Eigengene dendrogram and heatmap",
                      marDendro = eigengene_marDendro,
                      marHeatmap = eigengene_marHeatmap,
                      cex.lab = eigengene_cex_lab,
                      xLabelsAngle = eigengene_xLabelsAngle)

dev.off()
cat("✓ 图5已保存\n")

# ========================== 图6：TOM plot ====================================
cat("\n========== 图6：TOM plot ==========\n")

# 确保dissTOM存在
if (!exists("dissTOM")) {
  cat("计算dissTOM矩阵...\n")
  if (!exists("TOM")) {
    cat("计算TOM矩阵（这可能需要几分钟）...\n")
    adjacency <- adjacency(t(expr_matrix_norm), power = 6, type = "signed")
    TOM <- TOMsimilarity(adjacency, TOMType = "signed")
  }
  dissTOM <- 1 - TOM
}

# 随机选择基因
nSelect <- min(tom_plot_nSelect, nrow(dissTOM))
set.seed(tom_plot_seed)
select <- sample(nrow(dissTOM), size = nSelect)
selectTOM <- dissTOM[select, select]
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]

# 对TOM进行幂次转换
plotDiss <- selectTOM^tom_plot_power
diag(plotDiss) <- NA

pdf(file.path(output_dir, "6_Network_heatmap_plot.pdf"), 
    width = tom_plot_width, height = tom_plot_height)

TOMplot(plotDiss, selectTree, selectColors,
        col = colorpanel(250, tom_plot_colors[1], tom_plot_colors[2], tom_plot_colors[3]),
        main = tom_plot_main)

dev.off()
cat("✓ 图6已保存\n")

# ========================== 生成报告 ====================================
cat("\n========== 生成可视化报告 ==========\n")

sink(file.path(output_dir, "visualization_report.txt"))

cat("================================================================================\n")
cat("                    WGCNA可视化分析报告\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、输出文件列表\n")
cat("--------------------------------------------------------------------------------\n")
cat("1. 1_Sample_dendrogram_and_trait_heatmap.pdf - 样本聚类树 + 性状热图\n")
cat("2. 2_Gene_dendrogram_and_module_colors.pdf - 基因聚类树 + 模块颜色\n")
cat("3. 3_Module_trait_relationships.pdf - 模块-性状关系热图\n")
cat("4. 4_Module_enrichment_bubble.pdf - 模块富集气泡图\n")
cat("5. 5_Eigengene_dendrogram_and_heatmap.pdf - 模块特征基因树和热图\n")
cat("6. 6_Network_heatmap_plot.pdf - 基因网络热图（TOM plot）\n")

cat("\n二、模块统计\n")
cat("--------------------------------------------------------------------------------\n")
module_counts <- table(moduleColors)
for (mod in names(module_counts)) {
  cat(sprintf("模块 %s: %d 个基因\n", mod, module_counts[mod]))
}

cat("\n三、当前使用的参数\n")
cat("--------------------------------------------------------------------------------\n")
cat("图1尺寸:", sample_tree_width, "x", sample_tree_height, "\n")
cat("图2尺寸:", gene_dendro_width, "x", gene_dendro_height, "\n")
cat("图3尺寸:", module_trait_width, "x", module_trait_height, "\n")
cat("图4尺寸:", bubble_width, "x", bubble_height, "\n")
cat("图5尺寸:", eigengene_width, "x", eigengene_height, "\n")
cat("图6尺寸:", tom_plot_width, "x", tom_plot_height, "\n")
cat("TOM plot显示基因数:", nSelect, "\n")

sink()

cat("\n========================================\n")
cat("WGCNA可视化完成！\n")
cat("输出目录:", output_dir, "\n")
cat("========================================\n")