# ============================================================
# GSE199515 人类 TNBC 成纤维细胞 scRNA-seq 分析
# 目的：验证 CIPS 基因集在 CAF 中的特异性高表达
# 修改：
#   1. 基因集评分改用 ssGSEA（修复版）
#   2. 添加多重比较校正（BH方法）
#   3. 设置随机种子
#   4. 仅显示P值数值，无星号标注
#   5. P值位置自定义参数移至代码最前端
# ============================================================

# -------------------- 0. 设置随机种子 --------------------
set.seed(42)
cat("随机种子已设置: 42\n")

# ============================================================================
# ====================== 用户自定义参数区域 ==================================
# ============================================================================

# -------------------- 工作路径 --------------------
work_dir <- "D:/zsy/SX/GEO-scRNA-seq"
cips_file <- "D:/zsy/SX/GEO-scRNA-seq/gene_pathway_stats.csv"

# -------------------- CAF cluster 手动指定 --------------------
caf_clusters <- c("0", "1", "2")  # 请根据实际CAF Score打印结果修改

# -------------------- IL6/ACTA2 双阳性阈值 --------------------
il6_acta2_threshold <- 0.5

# -------------------- 主图A：UMAP参数 --------------------
umap_width <- 10
umap_height <- 6
umap_pt_size <- 0.6

# 三群体配色
color_double_positive <- "#C60036"    # IL6+ ACTA2+ 双阳性CAF
color_other_caf <- "#25377F"          # 其他CAF
color_other_cells <- "grey85"         # 其他细胞

# 三群体名称
name_double_positive <- "IL6+ ACTA2+ CAF"
name_other_caf <- "Other CAF"
name_other_cells <- "non-cells"

# -------------------- 主图B：箱线图参数 --------------------
boxplot_width <- 6
boxplot_height <- 7

# P值标注参数（自定义位置）
pval_text_size <- 6          # P值字体大小
pval_text_color <- "black"     # P值文字颜色
pval_offset <- 0.05            # 第一个P值的Y轴偏移比例（相对于数据范围）
pval_spacing <- 0.15           # 多个P值之间的间距比例
pval_line_color <- "grey40"    # 括号线颜色
pval_line_size <- 0.5          # 括号线粗细

# 自定义每个比较的精确Y位置（设为NULL则自动计算）
# 格式：list("组1 vs 组2" = y位置, ...)
custom_pval_positions <- NULL
# 使用示例（取消注释并修改数值）：
# custom_pval_positions <- list(
#   "IL6+ ACTA2+ CAF vs Other CAF" = 1.8,
#   "IL6+ ACTA2+ CAF vs Other cells" = 2.0,
#   "Other CAF vs Other cells" = 2.2
# )

# -------------------- 附属图：CAF亚型比较图参数 --------------------
subtype_box_width <- 9
subtype_box_height <- 6
subtype_pval_text_size <- 4
subtype_pval_offset <- 0.06
subtype_pval_spacing <- 0.10

# -------------------- 配色方案 --------------------
caf_other_colors <- c("CAF" = "#E64B35", "Other" = "#4DBBD5")
sample_colors <- c("TNBC1" = "#3B9AB2", "TNBC2" = "#E1AF00", "TNBC3" = "#F21A00")
cips_heatmap_colors <- c("grey90", "lightblue", "orange", "red", "darkred")

# -------------------- 统计检验参数 --------------------
p_adjust_method <- "BH"        # 多重比较校正方法
significance_level <- 0.05     # 显著性阈值

# ============================================================================
# ====================== 环境准备 ============================================
# ============================================================================

library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(pheatmap)

setwd(work_dir)

# -------------------- 读取 CIPS 基因列表 --------------------
cips_genes_human <- read.csv(cips_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
cips_genes_human <- unique(cips_genes_human[!is.na(cips_genes_human) & cips_genes_human != ""])

cat("共读取", length(cips_genes_human), "个 CIPS 基因\n")
cat("前10个基因:", paste(head(cips_genes_human, 10), collapse = ", "), "\n")

# ============================================================================
# ====================== ssGSEA 算法实现 =====================================
# ============================================================================

cat("\n========== 定义 ssGSEA 算法 ==========\n")

ssgsea_score <- function(expression_matrix, gene_set, alpha = 0.25, verbose = TRUE) {
  # 基于GSEA官方算法的ssGSEA实现
  # 参考文献: Subramanian et al., PNAS 2005; Barbie et al., Nature 2009
  
  expression_matrix <- as.matrix(expression_matrix)
  mode(expression_matrix) <- "numeric"
  
  n_genes <- nrow(expression_matrix)
  n_cells <- ncol(expression_matrix)
  cell_names <- colnames(expression_matrix)
  
  available_genes <- intersect(gene_set, rownames(expression_matrix))
  n_gs <- length(available_genes)
  
  if (verbose) {
    cat(sprintf("  总基因数: %d\n", n_genes))
    cat(sprintf("  基因集命中数: %d\n", n_gs))
    cat(sprintf("  细胞数: %d\n", n_cells))
    cat(sprintf("  alpha: %.2f\n", alpha))
  }
  
  if (n_gs == 0) stop("无基因集基因在表达矩阵中")
  
  scores <- numeric(n_cells)
  names(scores) <- cell_names
  
  for (s in 1:n_cells) {
    expr <- expression_matrix[, s]
    order_idx <- order(expr, decreasing = TRUE)
    ranked_expr <- expr[order_idx]
    ranked_genes <- rownames(expression_matrix)[order_idx]
    
    gs_positions <- which(ranked_genes %in% available_genes)
    
    if (length(gs_positions) == 0) {
      scores[s] <- 0
      next
    }
    
    N <- n_genes
    Nh <- length(gs_positions)
    
    es <- 0
    max_es <- 0
    min_es <- 0
    
    for (i in 1:N) {
      if (i %in% gs_positions) {
        weight <- abs(ranked_expr[i]) ^ alpha
        es <- es + weight
      } else {
        es <- es - (1 / (N - Nh))
      }
      if (es > max_es) max_es <- es
      if (es < min_es) min_es <- es
    }
    
    scores[s] <- if (abs(max_es) >= abs(min_es)) max_es else min_es
    
    if (verbose && s %% 500 == 0) {
      cat(sprintf("  已处理 %d/%d 个细胞...\n", s, n_cells))
    }
  }
  
  if (max(scores, na.rm = TRUE) > min(scores, na.rm = TRUE)) {
    scores_norm <- 2 * (scores - min(scores, na.rm = TRUE)) / 
      (max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE)) - 1
    if (verbose) cat("  分数归一化到[-1, 1]区间\n")
  } else {
    scores_norm <- scores
  }
  
  names(scores_norm) <- cell_names
  return(scores_norm)
}

cat("✓ ssGSEA 算法定义完成\n")

# ============================================================================
# ====================== 多重比较校正函数 =====================================
# ============================================================================

calculate_adjusted_pairwise <- function(data_df, score_col, group_col, 
                                        method = "BH", sig_level = 0.05) {
  groups <- sort(unique(data_df[[group_col]]))
  n_groups <- length(groups)
  
  results <- data.frame()
  
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      g1 <- groups[i]
      g2 <- groups[j]
      
      s1 <- data_df[[score_col]][data_df[[group_col]] == g1]
      s2 <- data_df[[score_col]][data_df[[group_col]] == g2]
      
      if (length(s1) >= 3 && length(s2) >= 3) {
        test <- wilcox.test(s1, s2)
        results <- rbind(results, data.frame(
          Group1 = g1,
          Group2 = g2,
          Comparison = paste(g1, "vs", g2),
          Raw_P = test$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(results) > 0) {
    results$Adjusted_P <- p.adjust(results$Raw_P, method = method)
    results$Significant <- results$Adjusted_P < sig_level
    
    results$P_label <- sapply(results$Adjusted_P, function(p) {
      if (p < 0.0001) return("P < 0.0001")
      if (p < 0.001) return(sprintf("P = %.4f", p))
      if (p < 0.01) return(sprintf("P = %.4f", p))
      return(sprintf("P = %.3f", p))
    })
  }
  
  return(results)
}

# ============================================================================
# ====================== 添加P值标注的通用函数 ================================
# ============================================================================

add_pvalue_annotations <- function(plot_obj, pairwise_results, group_positions, 
                                   data_y_max, data_y_range,
                                   text_size = 4.5, text_color = "black",
                                   offset = 0.06, spacing = 0.10,
                                   line_color = "grey40", line_size = 0.3,
                                   custom_positions = NULL) {
  if (nrow(pairwise_results) == 0) return(plot_obj)
  
  for (k in 1:nrow(pairwise_results)) {
    comp <- pairwise_results[k, ]
    
    x1 <- group_positions[comp$Group1]
    x2 <- group_positions[comp$Group2]
    
    # 计算Y位置
    if (!is.null(custom_positions)) {
      comp_key <- paste(comp$Group1, "vs", comp$Group2)
      comp_key_rev <- paste(comp$Group2, "vs", comp$Group1)
      
      if (comp_key %in% names(custom_positions)) {
        y_pos <- custom_positions[[comp_key]]
      } else if (comp_key_rev %in% names(custom_positions)) {
        y_pos <- custom_positions[[comp_key_rev]]
      } else {
        y_pos <- data_y_max + data_y_range * (offset + k * spacing)
      }
    } else {
      y_pos <- data_y_max + data_y_range * (offset + k * spacing)
    }
    
    # 添加括号线
    plot_obj <- plot_obj +
      annotate("segment", x = x1, xend = x1, 
               y = y_pos - data_y_range * 0.02, yend = y_pos,
               size = line_size, color = line_color) +
      annotate("segment", x = x2, xend = x2, 
               y = y_pos - data_y_range * 0.02, yend = y_pos,
               size = line_size, color = line_color) +
      annotate("segment", x = x1, xend = x2, 
               y = y_pos, yend = y_pos,
               size = line_size, color = line_color)
    
    # 添加P值文字（仅数值）
    plot_obj <- plot_obj +
      annotate("text", x = (x1 + x2) / 2, 
               y = y_pos + data_y_range * 0.02,
               label = comp$P_label, 
               size = text_size, color = text_color)
  }
  
  return(plot_obj)
}

# ============================================================================
# ====================== 数据读取 ============================================
# ============================================================================

sample_names <- c("TNBC1", "TNBC2", "TNBC3")
seurat_list <- list()

for (i in 1:length(sample_names)) {
  sample <- sample_names[i]
  sample_path <- file.path(work_dir, sample)
  
  cat("\n读取样本:", sample, "\n")
  cat("路径:", sample_path, "\n")
  
  if (!dir.exists(sample_path)) stop("找不到文件夹: ", sample_path)
  
  counts <- Read10X(sample_path)
  
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample,
    min.cells = 3,
    min.features = 200
  )
  
  seurat_obj$orig_sample <- sample
  seurat_list[[sample]] <- seurat_obj
  
  cat("样本", sample, "读取完成，细胞数:", ncol(seurat_obj), "\n")
}

# ============================================================================
# ====================== 合并与质控 ==========================================
# ============================================================================

seurat_merged <- merge(seurat_list[["TNBC1"]], 
                       y = seurat_list[c("TNBC2", "TNBC3")],
                       add.cell.ids = sample_names)

seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")

p1 <- VlnPlot(seurat_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)
ggsave("01_QC_before_filtering.pdf", p1, width = 12, height = 5)

seurat_merged <- subset(
  seurat_merged,
  subset = nFeature_RNA > 200 & 
    nFeature_RNA < 6000 & 
    percent.mt < 20
)

cat("\n质控后剩余细胞数:", ncol(seurat_merged), "\n")

# ============================================================================
# ====================== 标准化与降维 ========================================
# ============================================================================

seurat_merged <- NormalizeData(seurat_merged, 
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)

seurat_merged <- FindVariableFeatures(seurat_merged,
                                      selection.method = "vst",
                                      nfeatures = 2000)

seurat_merged <- ScaleData(seurat_merged, 
                           vars.to.regress = c("percent.mt", "nFeature_RNA"))

seurat_merged <- RunPCA(seurat_merged, npcs = 50)

p2 <- ElbowPlot(seurat_merged, ndims = 50)
ggsave("02_ElbowPlot.pdf", p2, width = 6, height = 4)

n_pcs <- 20

seurat_merged <- RunUMAP(seurat_merged, dims = 1:n_pcs, seed.use = 42)
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:n_pcs)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5, random.seed = 42)

# ============================================================================
# ====================== 细胞类型鉴定 ========================================
# ============================================================================

caf_markers <- c("ACTA2", "FAP", "PDGFRA", "COL1A1", "COL1A2", "PDPN", "S100A4", "VIM")
other_markers <- c("PECAM1", "VWF", "PTPRC", "CD3D", "CD68", "EPCAM", "KRT8", "KRT18")

p3 <- FeaturePlot(seurat_merged, 
                  features = c("ACTA2", "FAP", "COL1A1", "PECAM1", "PTPRC", "EPCAM"),
                  ncol = 3, cols = c("grey90", "red"))
ggsave("03_Marker_FeaturePlot.pdf", p3, width = 14, height = 10)

p4 <- DotPlot(seurat_merged, features = c(caf_markers, other_markers),
              group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
ggsave("04_Marker_DotPlot.pdf", p4, width = 14, height = 7)

cluster_avg <- AverageExpression(seurat_merged, features = caf_markers, 
                                 group.by = "seurat_clusters")
cluster_avg_matrix <- as.matrix(cluster_avg$RNA)
caf_score <- colMeans(cluster_avg_matrix)

cat("\n========== 各 Cluster 的 CAF Score ==========\n")
for (i in 1:length(caf_score)) {
  cat("Cluster", names(caf_score)[i], ":", round(caf_score[i], 3), "\n")
}

cat("\n鉴定为 CAF 的 cluster:", paste(caf_clusters, collapse = ", "), "\n")

seurat_merged$cell_type <- "Other"
seurat_merged$cell_type[seurat_merged$seurat_clusters %in% caf_clusters] <- "CAF"

cat("\n========== 细胞类型分布 ==========\n")
print(table(seurat_merged$cell_type))

# ============================================================================
# ====================== 计算 CIPS ssGSEA 评分 ===============================
# ============================================================================

cat("\n========== 计算 CIPS ssGSEA 评分 ==========\n")

seurat_merged <- JoinLayers(seurat_merged)
expr_data <- GetAssayData(seurat_merged, layer = "data")

cips_present <- intersect(cips_genes_human, rownames(expr_data))
cat("CIPS 基因命中率:", length(cips_present), "/", length(cips_genes_human), 
    "(", round(length(cips_present)/length(cips_genes_human)*100, 1), "%)\n")

cips_scores <- ssgsea_score(expr_data, cips_present, alpha = 0.25, verbose = TRUE)
seurat_merged$CIPS_Score <- cips_scores[colnames(seurat_merged)]

cat(sprintf("ssGSEA 评分范围: [%.4f, %.4f]\n", 
            min(seurat_merged$CIPS_Score), max(seurat_merged$CIPS_Score)))

# ============================================================================
# ====================== 可视化 ==============================================
# ============================================================================

# UMAP 按样本
p5 <- DimPlot(seurat_merged, reduction = "umap", group.by = "orig_sample",
              cols = sample_colors, pt.size = 0.5) +
  labs(title = "Samples (TNBC1-3)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", aspect.ratio = 1)
ggsave("05_UMAP_by_sample.pdf", p5, width = 8, height = 6)

# UMAP 按 cluster
p6 <- DimPlot(seurat_merged, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 5, pt.size = 0.5) +
  labs(title = "Clusters") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", aspect.ratio = 1)
ggsave("06_UMAP_by_cluster.pdf", p6, width = 8, height = 6)

# UMAP 按细胞类型
p7 <- DimPlot(seurat_merged, reduction = "umap", group.by = "cell_type",
              cols = caf_other_colors, pt.size = 0.5, order = c("CAF", "Other")) +
  labs(title = "Cell Type Annotation") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", aspect.ratio = 1)
ggsave("07_UMAP_by_celltype.pdf", p7, width = 8, height = 6)

# UMAP CIPS 评分
p8 <- FeaturePlot(seurat_merged, features = "CIPS_Score", pt.size = 0.5, order = TRUE) +
  scale_color_gradientn(colors = cips_heatmap_colors) +
  labs(title = "CIPS ssGSEA Score") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", aspect.ratio = 1)
ggsave("08_UMAP_CIPS_score.pdf", p8, width = 8, height = 6)

# ============================================================================
# ====================== CAF vs Other 箱线图 =================================
# ============================================================================

cat("\n========== CAF vs Other 统计检验 ==========\n")

caf_scores <- seurat_merged$CIPS_Score[seurat_merged$cell_type == "CAF"]
other_scores <- seurat_merged$CIPS_Score[seurat_merged$cell_type == "Other"]

wilcox_caf_vs_other <- wilcox.test(caf_scores, other_scores)
p_raw <- wilcox_caf_vs_other$p.value

cat(sprintf("CAF vs Other Wilcoxon 检验: P = %.4e\n", p_raw))

p_label <- ifelse(p_raw < 0.0001, "P < 0.0001",
                  ifelse(p_raw < 0.001, sprintf("P = %.4f", p_raw),
                         ifelse(p_raw < 0.01, sprintf("P = %.4f", p_raw),
                                sprintf("P = %.3f", p_raw))))

p9 <- VlnPlot(seurat_merged, features = "CIPS_Score", group.by = "cell_type",
              pt.size = 0, cols = caf_other_colors) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5) +
  annotate("text", x = 1.5, y = max(seurat_merged$CIPS_Score) * 1.1,
           label = p_label, size = 5) +
  labs(title = "CIPS ssGSEA Score: CAF vs Other Cells",
       y = "CIPS ssGSEA Score", x = "") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12))
ggsave("09_CIPS_CAF_vs_Other_boxplot.pdf", p9, width = 6, height = 6)

# ============================================================================
# ====================== 样本间比较 ==========================================
# ============================================================================

cat("\n========== 样本间 CAF CIPS 评分比较 ==========\n")

caf_cells <- subset(seurat_merged, cell_type == "CAF")
n_caf <- ncol(caf_cells)
cat("CAF 细胞数量:", n_caf, "\n")

if (n_caf > 0) {
  caf_df <- data.frame(Sample = caf_cells$orig_sample, CIPS_Score = caf_cells$CIPS_Score)
  
  sample_pairwise <- calculate_adjusted_pairwise(caf_df, "CIPS_Score", "Sample", 
                                                 method = p_adjust_method)
  
  cat("\n样本间两两比较结果（BH校正后）:\n")
  print(sample_pairwise[, c("Comparison", "Raw_P", "Adjusted_P", "P_label")])
  
  p10 <- ggplot(caf_df, aes(x = Sample, y = CIPS_Score, fill = Sample)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
    scale_fill_manual(values = sample_colors) +
    labs(title = "CIPS ssGSEA Score Across TNBC Samples (CAF only)",
         subtitle = paste("Wilcoxon test with", p_adjust_method, "correction"),
         y = "CIPS ssGSEA Score", x = "") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "none", axis.text.x = element_text(size = 12, face = "bold"))
  
  ggsave("10_CIPS_across_samples.pdf", p10, width = 8, height = 6)
}

# ============================================================================
# ====================== 热图 ================================================
# ============================================================================

cat("\n========== 生成 CIPS 基因热图 ==========\n")

cips_expr <- GetAssayData(seurat_merged, layer = "data")
cips_present_filtered <- intersect(cips_present, rownames(cips_expr))
cips_expr <- cips_expr[cips_present_filtered, ]

var_genes <- apply(cips_expr, 1, var, na.rm = TRUE)
top50_genes <- names(sort(var_genes, decreasing = TRUE))[1:min(50, length(cips_present_filtered))]

clusters <- as.character(seurat_merged$seurat_clusters)
cluster_list <- sort(unique(clusters))

avg_matrix <- matrix(0, nrow = length(top50_genes), ncol = length(cluster_list))
rownames(avg_matrix) <- top50_genes
colnames(avg_matrix) <- cluster_list

for (i in 1:length(cluster_list)) {
  cl <- cluster_list[i]
  cells_in_cl <- which(clusters == cl)
  if (length(cells_in_cl) > 0) {
    genes_available <- intersect(top50_genes, rownames(cips_expr))
    if (length(genes_available) > 0) {
      sub_matrix <- cips_expr[genes_available, cells_in_cl, drop = FALSE]
      avg_matrix[genes_available, i] <- rowMeans(sub_matrix, na.rm = TRUE)
    }
  }
}

avg_matrix <- avg_matrix[rowSums(avg_matrix, na.rm = TRUE) > 0, ]

anno_col <- data.frame(
  Cell_Type = ifelse(colnames(avg_matrix) %in% caf_clusters, "CAF", "Other"),
  row.names = colnames(avg_matrix)
)
anno_colors <- list(Cell_Type = c(CAF = "#E64B35", Other = "#4DBBD5"))

pheatmap(avg_matrix,
         scale = "row",
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize_row = 7,
         fontsize_col = 10,
         main = "Top Variable CIPS Genes (ssGSEA method)",
         filename = "11_CIPS_heatmap.pdf",
         width = 10,
         height = 10)

# ============================================================================
# ====================== IL6+ ACTA2+ 双阳性分析 ==============================
# ============================================================================

cat("\n========== IL6+ ACTA2+ 双阳性 CAF 亚型分析 ==========\n")

has_il6 <- "IL6" %in% rownames(seurat_merged)
has_acta2 <- "ACTA2" %in% rownames(seurat_merged)

if (has_il6 && has_acta2) {
  
  expr_il6 <- GetAssayData(seurat_merged, layer = "data")["IL6", ]
  expr_acta2 <- GetAssayData(seurat_merged, layer = "data")["ACTA2", ]
  
  seurat_merged$caf_subtype <- "Other"
  seurat_merged$caf_subtype[expr_il6 > il6_acta2_threshold & expr_acta2 > il6_acta2_threshold] <- "IL6+_ACTA2+"
  seurat_merged$caf_subtype[expr_il6 > il6_acta2_threshold & expr_acta2 <= il6_acta2_threshold] <- "IL6+_ACTA2-"
  seurat_merged$caf_subtype[expr_il6 <= il6_acta2_threshold & expr_acta2 > il6_acta2_threshold] <- "IL6-_ACTA2+"
  
  cat("\n细胞亚型统计:\n")
  print(table(seurat_merged$caf_subtype))
  
  seurat_merged$highlight <- "Other"
  seurat_merged$highlight[seurat_merged$caf_subtype == "IL6+_ACTA2+"] <- "IL6+ ACTA2+ CAF"
  
  highlight_colors <- c("IL6+ ACTA2+ CAF" = "#E64B35", "Other" = "grey90")
  
  p12 <- DimPlot(seurat_merged, reduction = "umap", group.by = "highlight",
                 cols = highlight_colors, pt.size = 0.5,
                 order = c("IL6+ ACTA2+ CAF", "Other")) +
    labs(title = "IL6+ ACTA2+ Double-Positive CAFs") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right", aspect.ratio = 1)
  ggsave("12_UMAP_IL6_ACTA2_double_positive.pdf", p12, width = 8, height = 6)
  
  p13 <- FeaturePlot(seurat_merged, features = c("IL6", "ACTA2"),
                     ncol = 2, pt.size = 0.5, order = TRUE,
                     cols = c("grey90", "red")) +
    plot_annotation(title = "Expression of IL6 and ACTA2")
  ggsave("13_FeaturePlot_IL6_ACTA2.pdf", p13, width = 12, height = 5)
  
  # CAF 亚型 CIPS 评分比较
  caf_cells <- subset(seurat_merged, cell_type == "CAF")
  
  compare_df <- data.frame(
    Subtype = caf_cells$caf_subtype,
    CIPS_Score = caf_cells$CIPS_Score
  )
  
  subtype_counts <- table(compare_df$Subtype)
  valid_subtypes <- names(subtype_counts[subtype_counts >= 10])
  compare_df_filtered <- compare_df[compare_df$Subtype %in% valid_subtypes, ]
  
  if (length(valid_subtypes) >= 2) {
    
    subtype_pairwise <- calculate_adjusted_pairwise(
      compare_df_filtered, "CIPS_Score", "Subtype", method = p_adjust_method
    )
    
    cat("\nCAF 亚型两两比较结果（BH校正后）:\n")
    print(subtype_pairwise[, c("Comparison", "Raw_P", "Adjusted_P", "P_label")])
    
    subtype_colors <- c(
      "IL6+_ACTA2+" = "#E64B35",
      "IL6+_ACTA2-" = "#4DBBD5",
      "IL6-_ACTA2+" = "#00A087",
      "Other" = "grey70"
    )
    
    p14 <- ggplot(compare_df_filtered, aes(x = Subtype, y = CIPS_Score, fill = Subtype)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.15, fill = "white", alpha = 0.5) +
      scale_fill_manual(values = subtype_colors) +
      labs(title = "CIPS ssGSEA Score Across CAF Subtypes",
           subtitle = paste("Wilcoxon test with", p_adjust_method, "correction"),
           y = "CIPS ssGSEA Score", x = "") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            legend.position = "none",
            axis.text.x = element_text(size = 11, angle = 15, hjust = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))
    
    # 添加P值标注
    if (nrow(subtype_pairwise) > 0) {
      subtype_positions <- setNames(1:length(valid_subtypes), valid_subtypes)
      y_max_sub <- max(compare_df_filtered$CIPS_Score, na.rm = TRUE)
      y_range_sub <- diff(range(compare_df_filtered$CIPS_Score, na.rm = TRUE))
      if (y_range_sub == 0) y_range_sub <- 1
      
      p14 <- add_pvalue_annotations(
        p14, subtype_pairwise, subtype_positions,
        y_max_sub, y_range_sub,
        text_size = subtype_pval_text_size,
        offset = subtype_pval_offset,
        spacing = subtype_pval_spacing
      )
    }
    
    ggsave("14_CIPS_score_by_CAF_subtype.pdf", p14, 
           width = subtype_box_width, height = subtype_box_height)
  }
}

# ============================================================================
# ★★★ 正文主图 ★★★
# ============================================================================

cat("\n========================================\n")
cat("开始生成正文主图（三群体分析）\n")
cat("========================================\n")

# 定义三群体
seurat_merged$group_three <- name_other_cells
seurat_merged$group_three[
  seurat_merged$cell_type == "CAF" & seurat_merged$caf_subtype == "IL6+_ACTA2+"
] <- name_double_positive
seurat_merged$group_three[
  seurat_merged$cell_type == "CAF" & seurat_merged$caf_subtype != "IL6+_ACTA2+"
] <- name_other_caf

cat("\n三群体统计:\n")
print(table(seurat_merged$group_three))

# -------------------- 主图 A：UMAP --------------------
colors_three <- c(color_double_positive, color_other_caf, color_other_cells)
names(colors_three) <- c(name_double_positive, name_other_caf, name_other_cells)

seurat_merged$group_three_factor <- factor(
  seurat_merged$group_three,
  levels = c(name_other_cells, name_other_caf, name_double_positive)
)

p_main_A <- DimPlot(seurat_merged, reduction = "umap", group.by = "group_three_factor",
                    pt.size = umap_pt_size, order = TRUE) +
  scale_color_manual(values = colors_three) +
  labs(tag = "A") +
  theme_minimal(base_size = 16) +
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    legend.position = c(1.2, 0.85),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("Main_Figure_A_ThreeGroup_UMAP.pdf", p_main_A, 
       width = umap_width, height = umap_height)

# -------------------- 主图 B：箱线图 --------------------
plot_df <- data.frame(
  Group = seurat_merged$group_three,
  CIPS_Score = seurat_merged$CIPS_Score
)
plot_df$Group <- factor(plot_df$Group,
                        levels = c(name_other_cells, name_other_caf, name_double_positive))

three_group_results <- calculate_adjusted_pairwise(
  plot_df, "CIPS_Score", "Group", method = p_adjust_method
)

cat("\n三群体 CIPS 评分比较（BH校正后）:\n")
print(three_group_results[, c("Comparison", "Raw_P", "Adjusted_P", "P_label")])

colors_box <- c(color_double_positive, color_other_caf, color_other_cells)
names(colors_box) <- c(name_double_positive, name_other_caf, name_other_cells)

p_main_B <- ggplot(plot_df, aes(x = Group, y = CIPS_Score, fill = Group)) +
  geom_violin(alpha = 0.75, trim = FALSE, linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 0.6,
               outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = colors_box) +
  labs(tag = "B", y = "CIPS ssGSEA Score", x = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(size = 14, angle = 15, hjust = 1),
    axis.title.y = element_text(size = 16),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))

# 添加P值标注（仅数值，无星号）
if (nrow(three_group_results) > 0) {
  group_positions <- c(name_double_positive = 3, name_other_caf = 2, name_other_cells = 1)
  y_max_main <- max(plot_df$CIPS_Score, na.rm = TRUE)
  y_range_main <- diff(range(plot_df$CIPS_Score, na.rm = TRUE))
  if (y_range_main == 0) y_range_main <- 1
  
  p_main_B <- add_pvalue_annotations(
    p_main_B, three_group_results, group_positions,
    y_max_main, y_range_main,
    text_size = pval_text_size,
    text_color = pval_text_color,
    offset = pval_offset,
    spacing = pval_spacing,
    line_color = pval_line_color,
    line_size = pval_line_size,
    custom_positions = custom_pval_positions
  )
}

ggsave("Main_Figure_B_CIPS_ThreeGroup_Boxplot.pdf", p_main_B,
       width = boxplot_width, height = boxplot_height)

# ============================================================================
# ====================== 保存统计报告 ========================================
# ============================================================================

sink("Statistics_ThreeGroup.txt")

cat("========================================\n")
cat("三群体 CIPS ssGSEA 评分统计报告\n")
cat("方法: Wilcoxon rank-sum test + BH 校正\n")
cat("========================================\n\n")

cat("1. 三群体细胞数量\n")
print(table(seurat_merged$group_three))

cat("\n2. 各群体 CIPS ssGSEA 评分统计\n")
for (grp in names(table(seurat_merged$group_three))) {
  scores <- seurat_merged$CIPS_Score[seurat_merged$group_three == grp]
  if (length(scores) > 0) {
    cat("\n", grp, ":\n", sep = "")
    cat("  均值:", round(mean(scores), 4), "\n")
    cat("  标准差:", round(sd(scores), 4), "\n")
    cat("  中位数:", round(median(scores), 4), "\n")
    cat("  细胞数:", length(scores), "\n")
  }
}

cat("\n3. 两两比较结果（BH 校正）\n")
print(three_group_results[, c("Comparison", "Raw_P", "Adjusted_P", "P_label")])

cat("\n4. 各样本中三群体分布\n")
print(table(seurat_merged$group_three, seurat_merged$orig_sample))

sink()

# -------------------- 保存最终对象 --------------------
saveRDS(seurat_merged, file = "GSE199515_seurat_final.rds")

cat("\n统计报告已保存: Statistics_ThreeGroup.txt\n")
cat("最终 Seurat 对象已保存为: GSE199515_seurat_final.rds\n")
cat("\n========================================\n")
cat("★ 所有分析完成！★\n")
cat("========================================\n")