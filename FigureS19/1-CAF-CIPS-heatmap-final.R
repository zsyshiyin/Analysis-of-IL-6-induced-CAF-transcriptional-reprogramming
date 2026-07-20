# ============================================================================
# 基因表达热图可视化脚本（支持指定基因列表 - 4组样品版）
# ============================================================================

cat("\n========== 检查并安装必要的包 ==========\n")

required_packages <- c("pheatmap", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
  cat(sprintf("✓ %s 已加载\n", pkg))
}

# ========================== 用户自定义参数 ============================
work_path <- "D:/zsy/SX/Fomal-final/27-CIPS-Tumor-heatmap"
setwd(work_path)

expression_file <- "Tumor-mRNA-final.csv"
gene_list_file <- "gene_pathway_stats.csv"

sample_groups_config <- list(
  "DMEM" = "DMEM", "Con" = "Con", "I6" = "I6", "R2" = "R2"
)
group_order <- c("DMEM", "Con", "I6", "R2")
group_colors_custom <- c(
  "DMEM" = "#2AA7DE", "Con" = "#25377F", "I6" = "#C60036", "R2" = "#E4945A"
)

heatmap_color_low <- "#2C7BB6"
heatmap_color_mid <- "white"
heatmap_color_high <- "#D7191C"

show_rownames <- TRUE
fontsize_row <- 16
cluster_cols <- FALSE
cluster_rows <- TRUE

# 无聚类热图
show_no_cluster_heatmap <- TRUE

output_dir <- "heatmap_results_4groups_final2"

# ========================== 创建输出目录 ====================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ========================== 读取表达矩阵 ====================================
cat("\n========== 读取表达矩阵 ==========\n")

expression_data <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("原始数据维度: %d 基因 × %d 列\n", nrow(expression_data), ncol(expression_data)))

rownames(expression_data) <- expression_data[, 1]
expression_data <- expression_data[, -1]
expression_matrix <- as.matrix(expression_data)
mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- toupper(rownames(expression_matrix))

cat(sprintf("表达矩阵维度: %d 基因 × %d 样本\n", nrow(expression_matrix), ncol(expression_matrix)))

# ========================== 准备基因列表 ====================================
cat("\n========== 准备基因列表 ==========\n")

if (!is.null(gene_list_file) && file.exists(gene_list_file)) {
  cat("从文件读取基因列表:", gene_list_file, "\n")
  if (grepl("\\.csv$", gene_list_file)) {
    gene_list_data <- read.csv(gene_list_file, stringsAsFactors = FALSE)
    gene_col <- intersect(colnames(gene_list_data), c("GeneID", "Gene", "gene", "Symbol", "symbol"))
    if (length(gene_col) == 0) gene_col <- colnames(gene_list_data)[1]
    genes_to_plot <- gene_list_data[[gene_col[1]]]
  } else {
    genes_to_plot <- readLines(gene_list_file)
    genes_to_plot <- genes_to_plot[!grepl("^#", genes_to_plot) & genes_to_plot != ""]
  }
  cat(sprintf("从文件读取 %d 个基因\n", length(genes_to_plot)))
}

genes_to_plot_upper <- toupper(genes_to_plot)
cat(sprintf("指定基因总数: %d\n", length(genes_to_plot_upper)))

# ========================== 筛选基因 ====================================
cat("\n========== 筛选基因 ==========\n")

available_genes <- intersect(genes_to_plot_upper, rownames(expression_matrix))
missing_genes <- setdiff(genes_to_plot_upper, rownames(expression_matrix))

cat(sprintf("在表达矩阵中找到的基因: %d / %d\n", length(available_genes), length(genes_to_plot_upper)))

if (length(missing_genes) > 0) {
  cat("\n不在表达矩阵中的基因:\n")
  print(head(missing_genes, 20))
  if (length(missing_genes) > 20) cat(sprintf("... 共 %d 个基因缺失\n", length(missing_genes)))
}

if (length(available_genes) == 0) stop("错误：没有找到任何指定基因！")

expression_matrix <- expression_matrix[available_genes, , drop = FALSE]
cat(sprintf("提取基因数: %d\n", nrow(expression_matrix)))

# ========================== 识别样本分组 ====================================
cat("\n========== 识别样本分组 ==========\n")

sample_names <- colnames(expression_matrix)
sample_groups <- data.frame(SampleID = sample_names, Group = NA, stringsAsFactors = FALSE)

for (prefix in names(sample_groups_config)) {
  matching_samples <- grep(paste0("^", prefix), sample_names, value = TRUE)
  if (length(matching_samples) > 0) {
    sample_groups$Group[sample_groups$SampleID %in% matching_samples] <- sample_groups_config[[prefix]]
  }
}

sample_groups <- sample_groups[!is.na(sample_groups$Group), ]
expression_matrix <- expression_matrix[, sample_groups$SampleID, drop = FALSE]

sample_groups$Group <- factor(sample_groups$Group, levels = group_order)
sample_groups <- sample_groups[order(sample_groups$Group), ]
expression_matrix <- expression_matrix[, sample_groups$SampleID, drop = FALSE]

cat("分组信息:\n")
print(table(sample_groups$Group))

# ========================== 数据预处理 ====================================
cat("\n========== 数据预处理 ==========\n")

expression_matrix[is.na(expression_matrix)] <- 0
expression_matrix[is.infinite(expression_matrix)] <- 0

if (min(expression_matrix) >= 0) {
  expression_matrix <- log2(expression_matrix + 1)
  cat("已应用log2转换\n")
}

heatmap_data <- t(scale(t(expression_matrix)))
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data[heatmap_data > 3] <- 3
heatmap_data[heatmap_data < -3] <- -3

cat(sprintf("最终数据: %d 基因 × %d 样本\n", nrow(heatmap_data), ncol(heatmap_data)))

# ========================== 创建注释 ====================================
annotation_col <- data.frame(Group = sample_groups$Group, row.names = sample_groups$SampleID)
existing_groups <- unique(sample_groups$Group)
group_colors <- group_colors_custom[existing_groups]
annotation_colors <- list(Group = group_colors)

# ========================== 生成热图 ====================================
cat("\n========== 生成热图 ==========\n")

heatmap_colors <- colorRampPalette(c(heatmap_color_low, heatmap_color_mid, heatmap_color_high))(100)
n_genes <- nrow(heatmap_data)
n_samples <- ncol(heatmap_data)

height <- if (n_genes > 100) max(8, n_genes * 0.08) else if (n_genes > 50) max(8, n_genes * 0.1) else max(4, n_genes * 0.2)
width <- max(6, n_samples * 0.15)
fontsize_row_auto <- if (n_genes > 100) 5 else if (n_genes > 50) 6 else fontsize_row

title_main <- sprintf("Gene Expression Heatmap (%d genes)", n_genes)

# ========== 1. 聚类热图 ==========
png(file.path(output_dir, "gene_expression_heatmap_clustered.png"), 
    width = width, height = height, units = "in", res = 150)
pheatmap(heatmap_data, color = heatmap_colors, annotation_col = annotation_col,
         annotation_colors = annotation_colors, show_rownames = show_rownames,
         cluster_rows = TRUE, cluster_cols = cluster_cols,
         fontsize_row = fontsize_row_auto, fontsize_col = 10, 
         main = paste(title_main, "(Clustered)"), border_color = NA)
dev.off()
cat("✓ 聚类热图 PNG 已保存\n")

pdf(file.path(output_dir, "gene_expression_heatmap_clustered.pdf"), 
    width = width, height = height)
pheatmap(heatmap_data, color = heatmap_colors, annotation_col = annotation_col,
         annotation_colors = annotation_colors, show_rownames = show_rownames,
         cluster_rows = TRUE, cluster_cols = cluster_cols,
         fontsize_row = fontsize_row_auto, fontsize_col = 10, 
         main = paste(title_main, "(Clustered)"), border_color = NA)
dev.off()
cat("✓ 聚类热图 PDF 已保存\n")

# ========== 2. 无聚类热图 ==========
if (show_no_cluster_heatmap) {
  cat("\n生成无聚类热图...\n")
  
  png(file.path(output_dir, "gene_expression_heatmap_no_cluster.png"), 
      width = width, height = height, units = "in", res = 150)
  pheatmap(heatmap_data, color = heatmap_colors, annotation_col = annotation_col,
           annotation_colors = annotation_colors, show_rownames = show_rownames,
           cluster_rows = FALSE, cluster_cols = FALSE,
           fontsize_row = fontsize_row_auto, fontsize_col = 10, 
           main = paste(title_main, "(Sample Order: DMEM, Con, I6, R2)"), 
           border_color = NA)
  dev.off()
  cat("✓ 无聚类热图 PNG 已保存\n")
  
  pdf(file.path(output_dir, "gene_expression_heatmap_no_cluster.pdf"), 
      width = width, height = height)
  pheatmap(heatmap_data, color = heatmap_colors, annotation_col = annotation_col,
           annotation_colors = annotation_colors, show_rownames = show_rownames,
           cluster_rows = FALSE, cluster_cols = FALSE,
           fontsize_row = fontsize_row_auto, fontsize_col = 10, 
           main = paste(title_main, "(Sample Order: DMEM, Con, I6, R2)"), 
           border_color = NA)
  dev.off()
  cat("✓ 无聚类热图 PDF 已保存\n")
}

# ========================== 保存数据 ====================================
cat("\n========== 保存处理后的数据 ==========\n")

heatmap_data_df <- as.data.frame(heatmap_data)
heatmap_data_df <- cbind(GeneID = rownames(heatmap_data_df), heatmap_data_df)
write.csv(heatmap_data_df, file.path(output_dir, "normalized_data.csv"), row.names = FALSE)

expr_data_df <- as.data.frame(expression_matrix)
expr_data_df <- cbind(GeneID = rownames(expr_data_df), expr_data_df)
write.csv(expr_data_df, file.path(output_dir, "expression_data_filtered.csv"), row.names = FALSE)

write.csv(sample_groups, file.path(output_dir, "sample_groups.csv"), row.names = FALSE)

gene_list_used <- data.frame(GeneID = rownames(expression_matrix), stringsAsFactors = FALSE)
write.csv(gene_list_used, file.path(output_dir, "genes_used.csv"), row.names = FALSE)

# ========================== 输出报告 ====================================
cat("\n========================================\n")
cat("热图可视化完成！\n")
cat("\n输出文件:\n")
cat("  ├── gene_expression_heatmap_clustered.png\n")
cat("  ├── gene_expression_heatmap_clustered.pdf\n")
if (show_no_cluster_heatmap) {
  cat("  ├── gene_expression_heatmap_no_cluster.png\n")
  cat("  ├── gene_expression_heatmap_no_cluster.pdf\n")
}
cat("  ├── normalized_data.csv\n")
cat("  ├── expression_data_filtered.csv\n")
cat("  ├── sample_groups.csv\n")
cat("  └── genes_used.csv\n")

cat(sprintf("\n指定基因: %d | 实际使用: %d | 缺失: %d\n", 
            length(genes_to_plot_upper), nrow(heatmap_data), length(missing_genes)))
cat("========================================\n")