# ============================================================
# 空间转录组分析 - 修复内存问题版
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(future)

# 增加内存限制
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB

# 路径设置
work_dir <- "D:/zsy/SX/GEO-mIHC"
data_dir <- file.path(work_dir, "filtered_feature_bc_matrix")
cips_file <- "D:/zsy/SX/GEO-mIHC/gene_pathway_stats.csv"
setwd(work_dir)

cat("========================================\n")
cat("空间转录组分析（内存优化版）\n")
cat("========================================\n")

# -------------------- 读取数据（如果已有对象可跳过）--------------------
if (!exists("bc_spatial")) {
  cat("\n正在读取数据...\n")
  spatial_image <- Read10X_Image(
    image.dir = file.path(data_dir, "spatial"),
    image.name = "tissue_lowres_image.png"
  )
  expression_matrix <- Read10X(data.dir = data_dir)
  bc_spatial <- CreateSeuratObject(counts = expression_matrix, assay = "Spatial")
  bc_spatial[["breast_cancer_slice"]] <- spatial_image
}

# -------------------- 质控 --------------------
cat("\n质控...\n")
bc_spatial[["percent.mt"]] <- PercentageFeatureSet(bc_spatial, pattern = "^MT-")
bc_spatial <- subset(bc_spatial, subset = nFeature_Spatial > 200 & percent.mt < 30)
cat("质控后 Spots:", ncol(bc_spatial), "\n")

# -------------------- 标准化（两种方法任选）--------------------

# 方法一：尝试 SCTransform（如果内存足够）
cat("\n尝试 SCTransform...\n")
bc_spatial <- tryCatch({
  SCTransform(bc_spatial, assay = "Spatial", verbose = FALSE)
}, error = function(e) {
  cat("SCTransform 失败，改用 LogNormalize...\n")
  bc_spatial <- NormalizeData(bc_spatial, assay = "Spatial")
  bc_spatial <- FindVariableFeatures(bc_spatial, assay = "Spatial", nfeatures = 2000)
  bc_spatial <- ScaleData(bc_spatial, assay = "Spatial")
  return(bc_spatial)
})

# -------------------- 降维与聚类 --------------------
cat("降维与聚类...\n")

# 确定使用哪个 assay
assay_use <- ifelse("SCT" %in% names(bc_spatial@assays), "SCT", "Spatial")

bc_spatial <- RunPCA(bc_spatial, assay = assay_use, verbose = FALSE)
bc_spatial <- FindNeighbors(bc_spatial, reduction = "pca", dims = 1:30)
bc_spatial <- FindClusters(bc_spatial, verbose = FALSE, resolution = 0.5)
bc_spatial <- RunUMAP(bc_spatial, reduction = "pca", dims = 1:30)
cat("聚类数:", length(unique(bc_spatial$seurat_clusters)), "\n")

# -------------------- 计算评分 --------------------
cat("\n计算评分...\n")

# CIPS 评分
if (file.exists(cips_file)) {
  cips_genes <- read.csv(cips_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  cips_genes <- unique(cips_genes[!is.na(cips_genes) & cips_genes != ""])
  cips_present <- intersect(cips_genes, rownames(bc_spatial))
  cat("  CIPS 命中:", length(cips_present), "/", length(cips_genes), "\n")
  
  if (length(cips_present) > 0) {
    bc_spatial <- AddModuleScore(bc_spatial, features = list(CIPS = cips_present), 
                                 name = "CIPS_Score", ctrl = 100)
  }
}

# IL6+ACTA2 评分
if ("IL6" %in% rownames(bc_spatial) && "ACTA2" %in% rownames(bc_spatial)) {
  bc_spatial <- AddModuleScore(bc_spatial, features = list(IL6_ACTA2 = c("IL6", "ACTA2")),
                               name = "IL6_ACTA2_Score", ctrl = 2)
  cat("  IL6+ACTA2 评分已计算\n")
}

# 肿瘤和基质评分
tumor_markers <- c("EPCAM", "KRT8", "KRT18", "KRT19")
caf_markers <- c("ACTA2", "COL1A1", "FAP", "PDPN")

tumor_present <- intersect(tumor_markers, rownames(bc_spatial))
caf_present <- intersect(caf_markers, rownames(bc_spatial))

if (length(tumor_present) > 0) {
  bc_spatial <- AddModuleScore(bc_spatial, features = list(Tumor = tumor_present),
                               name = "Tumor_Score", ctrl = min(length(tumor_present), 5))
}

if (length(caf_present) > 0) {
  bc_spatial <- AddModuleScore(bc_spatial, features = list(Stroma = caf_present),
                               name = "Stroma_Score", ctrl = min(length(caf_present), 5))
}

# -------------------- 保存 --------------------
saveRDS(bc_spatial, file = file.path(work_dir, "bc_spatial_processed.rds"))
cat("\n对象已保存: bc_spatial_processed.rds\n")

# -------------------- 测试可视化 --------------------
cat("\n生成测试图...\n")

p_he <- SpatialDimPlot(bc_spatial, alpha = 0) + 
  labs(title = "H&E Stained Tissue Section") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

ggsave("Test_Spatial_HE.pdf", p_he, width = 6, height = 6)

if ("CIPS_Score1" %in% colnames(bc_spatial@meta.data)) {
  p_cips <- SpatialFeaturePlot(bc_spatial, features = "CIPS_Score1",
                               alpha = c(0.1, 1)) +
    scale_fill_gradientn(colors = c("grey90", "lightblue", "darkred")) +
    labs(title = "CIPS Score") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("Test_Spatial_CIPS.pdf", p_cips, width = 6, height = 6)
}

cat("\n完成！\n")