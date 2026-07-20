# ============================================================================
# 基因表达热图可视化脚本（支持指定基因列表）
# 功能：对指定的基因列表进行热图可视化，显示基因名称
# ============================================================================

# ========================== 检查并安装必要的包 ====================================
cat("\n========== 检查并安装必要的包 ==========\n")

required_packages <- c("pheatmap", "ggplot2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
  cat(sprintf("✓ %s 已加载\n", pkg))
}

# ========================== 用户自定义参数设置区域 ============================

# 工作路径
work_path <- "D:/zsy/SX/Fomal-final/7-CIPS-heatmap"  # 请修改为实际工作路径
setwd(work_path)

# 输入文件参数
expression_file <- "gene_expression_matrix.csv"  # 处理后的表达矩阵文件

# ===== 基因列表设置（两种方式任选其一）=====
# 方式1：直接指定基因列表（推荐）
#genes_to_plot <- c("CDH2", "VIM", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2",
#"MMP2", "MMP9", "COL1A1", "COL1A2", "COL3A1", "ACTA2", "TGFB1",
#"CDH1", "EPCAM", "DSP", "OCLN", "KRT19", "KRT18")

# 方式2：从文件读取基因列表（如果指定了文件路径，将覆盖方式1）
gene_list_file <- "gene_pathway_stats.csv"  # 设置为文件路径，如 "gene_list.txt" 或 "gene_list.csv"
# gene_list_file <- "my_gene_list.txt"

# 样品分组信息
sample_groups_config <- list(
  "Con" = "Con",
  "I6" = "IL-6",
  "R2" = "(20S)G-Rh2"
)

# 组别颜色
group_colors_custom <- c(
  "Con" = "#25377F",
  "IL-6" = "#C60036",
  "(20S)G-Rh2" = "#E4945A"
)

# 热图颜色
heatmap_color_low <- "#2C7BB6"
heatmap_color_mid <- "white"
heatmap_color_high <- "#D7191C"

# 是否显示基因名称（如果基因数>200，建议设为FALSE，但这里我们指定了列表，通常不会太多）
show_rownames <- TRUE

# 基因名显示大小（如果基因名太长可以调小）
fontsize_row <- 8

# 输出目录
output_dir <- "heatmap_results_zhufigure2"

# ========== 前30个基因热图的自定义参数 ==========
top30_heatmap_params <- list(
  enabled = TRUE,                    # 是否生成前30个基因的热图
  logFC_column = "logFC",            # logFC列名（根据您的文件列名设置）
  use_abs_logFC = TRUE,              # 是否使用logFC的绝对值（TRUE: 上下调都考虑, FALSE: 只考虑上调）
  top_n = 30,                        # 选择前N个基因
  width = 6,                         # 热图宽度（英寸）
  height = 6,                       # 热图高度（英寸）
  fontsize_row = 9,                  # 基因名字体大小
  fontsize_col = 10,                 # 样本名字体大小
  show_rownames = TRUE,              # 是否显示基因名
  cluster_rows = TRUE,               # 是否对行聚类
  cluster_cols = TRUE                # 是否对列聚类
)

# ========================== 创建输出目录 ====================================
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ========================== 读取表达矩阵 ====================================
cat("\n========== 读取表达矩阵 ==========\n")

expression_data <- read.csv(expression_file, 
                            stringsAsFactors = FALSE,
                            check.names = FALSE)

cat(sprintf("原始数据维度: %d 基因 × %d 样本\n", nrow(expression_data), ncol(expression_data)))

# 将第一列设为行名
rownames(expression_data) <- expression_data[, 1]
expression_data <- expression_data[, -1]

# 转换为数值矩阵
expression_matrix <- as.matrix(expression_data)
mode(expression_matrix) <- "numeric"

# ========================== 读取基因统计信息（logFC等） ====================================
cat("\n========== 读取基因统计信息 ==========\n")

gene_stats <- NULL
if (file.exists(gene_list_file)) {
  cat("从文件读取基因统计信息:", gene_list_file, "\n")
  gene_stats <- read.csv(gene_list_file, stringsAsFactors = FALSE)
  cat(sprintf("基因统计文件维度: %d 基因 × %d 列\n", nrow(gene_stats), ncol(gene_stats)))
  
  # 检查logFC列是否存在
  if (!top30_heatmap_params$logFC_column %in% colnames(gene_stats)) {
    stop(sprintf("错误: 文件中未找到列 '%s'，请检查列名！", top30_heatmap_params$logFC_column))
  }
  
  cat(sprintf("✓ 找到logFC列: %s\n", top30_heatmap_params$logFC_column))
  logFC_values <- gene_stats[[top30_heatmap_params$logFC_column]]
  cat(sprintf("  logFC范围: [%.3f, %.3f]\n", min(logFC_values, na.rm=TRUE), max(logFC_values, na.rm=TRUE)))
  
} else {
  stop("错误: 未找到基因统计文件，请确保文件存在！")
}

# ========================== 读取或使用指定的基因列表 ====================================
cat("\n========== 准备基因列表 ==========\n")

# 优先从文件读取
if (!is.null(gene_list_file) && file.exists(gene_list_file)) {
  cat("从文件读取基因列表:", gene_list_file, "\n")
  
  # 判断文件格式
  if (grepl("\\.csv$", gene_list_file)) {
    gene_list_data <- read.csv(gene_list_file, stringsAsFactors = FALSE)
    # 尝试找到基因列
    gene_col <- NULL
    for (col in c("GeneID", "Gene", "gene", "Symbol", "symbol")) {
      if (col %in% colnames(gene_list_data)) {
        gene_col <- col
        break
      }
    }
    if (is.null(gene_col)) {
      gene_col <- colnames(gene_list_data)[1]  # 使用第一列
    }
    genes_to_plot <- gene_list_data[[gene_col]]
  } else {
    # txt文件，每行一个基因
    genes_to_plot <- readLines(gene_list_file)
    genes_to_plot <- genes_to_plot[!grepl("^#", genes_to_plot) & genes_to_plot != ""]
  }
  
  cat(sprintf("从文件读取 %d 个基因\n", length(genes_to_plot)))
}

cat(sprintf("指定基因列表: %d 个基因\n", length(genes_to_plot)))
cat("前10个基因:\n")
print(head(genes_to_plot, 10))

# ========================== 筛选基因 ====================================
cat("\n========== 筛选基因 ==========\n")

# 转换为大写以便匹配
genes_to_plot_upper <- toupper(genes_to_plot)
rownames_upper <- toupper(rownames(expression_matrix))

# 找出在表达矩阵中存在的基因
available_genes <- intersect(genes_to_plot_upper, rownames_upper)
missing_genes <- setdiff(genes_to_plot_upper, rownames_upper)

cat(sprintf("指定基因总数: %d\n", length(genes_to_plot)))
cat(sprintf("在表达矩阵中找到的基因: %d\n", length(available_genes)))

if (length(missing_genes) > 0) {
  cat("\n警告：以下基因不在表达矩阵中:\n")
  print(head(missing_genes, 20))
  if (length(missing_genes) > 20) {
    cat(sprintf("... 共 %d 个基因缺失\n", length(missing_genes)))
  }
}

if (length(available_genes) == 0) {
  stop("错误：没有找到任何指定基因！")
}

# 提取指定基因的表达矩阵
expression_matrix <- expression_matrix[tolower(rownames(expression_matrix)) %in% tolower(available_genes), , drop = FALSE]

# 按原顺序排序（如果可能）
original_order <- genes_to_plot_upper[genes_to_plot_upper %in% available_genes]
current_order <- rownames(expression_matrix)
new_order <- original_order[original_order %in% current_order]
expression_matrix <- expression_matrix[new_order, , drop = FALSE]

cat(sprintf("\n最终提取 %d 个基因\n", nrow(expression_matrix)))

# ========================== 识别样本分组 ====================================
cat("\n========== 识别样本分组 ==========\n")

sample_names <- colnames(expression_matrix)
sample_groups <- data.frame(
  SampleID = sample_names,
  Group = NA,
  stringsAsFactors = FALSE
)

for (prefix in names(sample_groups_config)) {
  matching_samples <- grep(paste0("^", prefix), sample_names, value = TRUE)
  if (length(matching_samples) > 0) {
    sample_groups$Group[sample_groups$SampleID %in% matching_samples] <- 
      sample_groups_config[[prefix]]
  }
}

# 移除未识别的样本
sample_groups <- sample_groups[!is.na(sample_groups$Group), ]
expression_matrix <- expression_matrix[, sample_groups$SampleID, drop = FALSE]

cat("分组信息:\n")
print(table(sample_groups$Group))

# ========================== 数据清理和标准化 ====================================
cat("\n========== 数据预处理 ==========\n")

# 1. 处理缺失值和无限值
cat("清理数据...\n")
expression_matrix[is.na(expression_matrix)] <- 0
expression_matrix[is.infinite(expression_matrix)] <- 0

# 2. Log2转换
cat("Log2转换...\n")
if (min(expression_matrix) >= 0) {
  expression_matrix <- log2(expression_matrix + 1)
  cat("  已应用log2转换\n")
}

# 3. Z-score标准化（对基因进行标准化）
cat("Z-score标准化...\n")
heatmap_data <- t(scale(t(expression_matrix)))

# 处理标准化后的NA（当标准差为0时会出现）
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data[is.infinite(heatmap_data)] <- 0

# 4. 限制范围（避免极端值影响颜色显示）
heatmap_data[heatmap_data > 3] <- 3
heatmap_data[heatmap_data < -3] <- -3

cat("数据预处理完成\n")
cat(sprintf("最终数据维度: %d 基因 × %d 样本\n", nrow(heatmap_data), ncol(heatmap_data)))

# ========================== 根据logFC排序基因 ====================================
cat("\n========== 根据logFC排序基因 ==========\n")

# 创建logFC数据框
gene_logFC <- data.frame(
  GeneID = rownames(expression_matrix),
  stringsAsFactors = FALSE
)

# 匹配logFC值（不区分大小写）
gene_logFC$GeneID_upper <- toupper(gene_logFC$GeneID)
gene_stats$GeneID_upper <- toupper(gene_stats[[1]])

# 添加logFC值
gene_logFC$logFC <- NA
for (i in 1:nrow(gene_logFC)) {
  match_idx <- which(gene_stats$GeneID_upper == gene_logFC$GeneID_upper[i])
  if (length(match_idx) > 0) {
    gene_logFC$logFC[i] <- gene_stats[match_idx[1], top30_heatmap_params$logFC_column]
  }
}

# 检查是否有缺失的logFC
missing_logFC <- sum(is.na(gene_logFC$logFC))
if (missing_logFC > 0) {
  cat(sprintf("警告: %d 个基因没有找到logFC值\n", missing_logFC))
}

# 计算绝对值logFC
if (top30_heatmap_params$use_abs_logFC) {
  gene_logFC$sort_value <- abs(gene_logFC$logFC)
  cat("使用logFC绝对值进行排序\n")
} else {
  gene_logFC$sort_value <- gene_logFC$logFC
  cat("使用logFC原始值进行排序\n")
}

# 按logFC绝对值排序（从大到小）
gene_logFC_sorted <- gene_logFC[order(-gene_logFC$sort_value), ]

cat("\nlogFC绝对值最大的前10个基因:\n")
print(head(gene_logFC_sorted[, c("GeneID", "logFC", "sort_value")], 10))

# ========================== 创建注释 ====================================
annotation_col <- data.frame(
  Group = sample_groups$Group,
  row.names = sample_groups$SampleID
)

# 设置分组颜色
existing_groups <- unique(sample_groups$Group)
group_colors <- group_colors_custom[existing_groups]
annotation_colors <- list(Group = group_colors)

# ========================== 生成所有基因的热图 ====================================
cat("\n========== 生成所有基因的热图 ==========\n")

# 创建颜色
heatmap_colors <- colorRampPalette(c(heatmap_color_low, heatmap_color_mid, heatmap_color_high))(100)

# 设置标题
title_text <- sprintf("Gene Expression Heatmap (%d genes)", nrow(heatmap_data))

# 确定图片高度（根据基因数动态调整）
n_genes <- nrow(heatmap_data)
if (n_genes > 100) {
  height <- max(11, n_genes * 0.1)
  fontsize_row <- 7
  cat(sprintf("基因数较多(%d)，调整高度为 %.1f英寸，字体大小为 %d\n", n_genes, height, fontsize_row))
} else if (n_genes > 50) {
  height <- max(8, n_genes * 0.1)
  fontsize_row <- 6
} else {
  height <- max(6, n_genes * 0.15)
  fontsize_row <- 8
}

# 设置宽度（根据样本数动态调整）
n_samples <- ncol(heatmap_data)
width <- max(6, n_samples * 0.2)

cat(sprintf("图片尺寸: %.1f × %.1f 英寸\n", width, height))

# 生成PNG
png_file <- file.path(output_dir, "gene_expression_heatmap.png")
cat("正在生成PNG热图...\n")

tryCatch({
  png(png_file, width = width, height = height, units = "in", res = 150)
  pheatmap(heatmap_data,
           color = heatmap_colors,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = show_rownames,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = title_text,
           fontsize_row = fontsize_row,
           fontsize_col = 10,
           border_color = NA)
  dev.off()
  cat(sprintf("✓ PNG热图已保存: %s\n", png_file))
}, error = function(e) {
  cat(sprintf("✗ PNG生成失败: %s\n", e$message))
  dev.off()
})

# 生成PDF
pdf_file <- file.path(output_dir, "gene_expression_heatmap.pdf")
cat("正在生成PDF热图...\n")

tryCatch({
  pdf(pdf_file, width = width, height = height)
  pheatmap(heatmap_data,
           color = heatmap_colors,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = show_rownames,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = title_text,
           fontsize_row = fontsize_row,
           fontsize_col = 10,
           border_color = NA)
  dev.off()
  cat(sprintf("✓ PDF热图已保存: %s\n", pdf_file))
}, error = function(e) {
  cat(sprintf("✗ PDF生成失败: %s\n", e$message))
  dev.off()
})

# ========================== 生成无聚类的热图（备选）============================
cat("\n========== 生成无聚类热图（备选）==========\n")

# 生成一个按原始顺序排列、无聚类的热图
png_file_no_cluster <- file.path(output_dir, "gene_expression_heatmap_no_cluster.png")

tryCatch({
  png(png_file_no_cluster, width = width, height = height, units = "in", res = 150)
  pheatmap(heatmap_data,
           color = heatmap_colors,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = show_rownames,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = paste(title_text, "(No clustering)"),
           fontsize_row = fontsize_row,
           fontsize_col = 10,
           border_color = NA)
  dev.off()
  cat(sprintf("✓ 无聚类PNG热图已保存: %s\n", png_file_no_cluster))
}, error = function(e) {
  cat(sprintf("✗ 无聚类热图生成失败: %s\n", e$message))
  dev.off()
})

# ==================== 生成logFC最大的前30个基因的热图 ====================
if (top30_heatmap_params$enabled && nrow(heatmap_data) >= top30_heatmap_params$top_n) {
  cat("\n========== 生成logFC最大的前30个基因的热图 ==========\n")
  
  top_n <- top30_heatmap_params$top_n
  
  # 过滤掉没有logFC值的基因
  gene_logFC_filtered <- gene_logFC_sorted[!is.na(gene_logFC_sorted$logFC), ]
  
  if (nrow(gene_logFC_filtered) < top_n) {
    cat(sprintf("警告: 只有 %d 个基因有logFC值，少于 %d\n", nrow(gene_logFC_filtered), top_n))
    top_n <- nrow(gene_logFC_filtered)
  }
  
  # 选择logFC绝对值最大的前N个基因
  top_genes <- head(gene_logFC_filtered$GeneID, top_n)
  
  # 提取这些基因的表达数据
  heatmap_data_top30 <- heatmap_data[top_genes, , drop = FALSE]
  
  # 按logFC绝对值排序（从大到小）
  sort_order <- match(top_genes, gene_logFC_filtered$GeneID)
  heatmap_data_top30 <- heatmap_data_top30[order(sort_order), , drop = FALSE]
  
  # 获取对应的logFC值
  top_genes_logFC <- gene_logFC_filtered[match(rownames(heatmap_data_top30), gene_logFC_filtered$GeneID), ]
  
  cat(sprintf("提取 %d 个logFC绝对值最大的基因\n", nrow(heatmap_data_top30)))
  cat("\n选择的基因及logFC值:\n")
  print(top_genes_logFC[, c("GeneID", "logFC", "sort_value")])
  
  cat(sprintf("\n使用自定义参数: 宽度=%.1f英寸, 高度=%.1f英寸, 基因字体大小=%d\n", 
              top30_heatmap_params$width, top30_heatmap_params$height, 
              top30_heatmap_params$fontsize_row))
  
  # 生成标题
  if (top30_heatmap_params$use_abs_logFC) {
    title_suffix <- "by |logFC| (Most differentially expressed)"
  } else {
    title_suffix <- "by logFC (Highest up-regulated)"
  }
  
  # 生成前30个基因的PNG热图
  png_file_top30 <- file.path(output_dir, "top30_logFC_genes_heatmap.png")
  cat("正在生成logFC最大前30个基因的PNG热图...\n")
  
  tryCatch({
    png(png_file_top30, 
        width = top30_heatmap_params$width, 
        height = top30_heatmap_params$height, 
        units = "in", 
        res = 150)
    pheatmap(heatmap_data_top30,
             color = heatmap_colors,
             annotation_col = annotation_col,
             annotation_colors = annotation_colors,
             show_rownames = top30_heatmap_params$show_rownames,
             cluster_rows = top30_heatmap_params$cluster_rows,
             cluster_cols = top30_heatmap_params$cluster_cols,
             main = sprintf("Top %d Genes by |logFC|\n%s", nrow(heatmap_data_top30), title_suffix),
             fontsize_row = top30_heatmap_params$fontsize_row,
             fontsize_col = top30_heatmap_params$fontsize_col,
             border_color = NA)
    dev.off()
    cat(sprintf("✓ logFC最大前30个基因PNG热图已保存: %s\n", png_file_top30))
  }, error = function(e) {
    cat(sprintf("✗ PNG生成失败: %s\n", e$message))
    dev.off()
  })
  
  # 生成前30个基因的PDF热图
  pdf_file_top30 <- file.path(output_dir, "top30_logFC_genes_heatmap.pdf")
  cat("正在生成logFC最大前30个基因的PDF热图...\n")
  
  tryCatch({
    pdf(pdf_file_top30, 
        width = top30_heatmap_params$width, 
        height = top30_heatmap_params$height)
    pheatmap(heatmap_data_top30,
             color = heatmap_colors,
             annotation_col = annotation_col,
             annotation_colors = annotation_colors,
             show_rownames = top30_heatmap_params$show_rownames,
             cluster_rows = top30_heatmap_params$cluster_rows,
             cluster_cols = top30_heatmap_params$cluster_cols,
             main = sprintf("Top %d Genes by |logFC|\n%s", nrow(heatmap_data_top30), title_suffix),
             fontsize_row = top30_heatmap_params$fontsize_row,
             fontsize_col = top30_heatmap_params$fontsize_col,
             border_color = NA)
    dev.off()
    cat(sprintf("✓ logFC最大前30个基因PDF热图已保存: %s\n", pdf_file_top30))
  }, error = function(e) {
    cat(sprintf("✗ PDF生成失败: %s\n", e$message))
    dev.off()
  })
  
  # 可选：生成无聚类的热图
  if (top30_heatmap_params$cluster_rows || top30_heatmap_params$cluster_cols) {
    png_file_top30_no_cluster <- file.path(output_dir, "top30_logFC_genes_heatmap_no_cluster.png")
    cat("正在生成logFC最大前30个基因无聚类PNG热图...\n")
    
    tryCatch({
      png(png_file_top30_no_cluster, 
          width = top30_heatmap_params$width, 
          height = top30_heatmap_params$height, 
          units = "in", 
          res = 150)
      pheatmap(heatmap_data_top30,
               color = heatmap_colors,
               annotation_col = annotation_col,
               annotation_colors = annotation_colors,
               show_rownames = top30_heatmap_params$show_rownames,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               main = sprintf("Top %d Genes by |logFC| (No clustering)\n%s", nrow(heatmap_data_top30), title_suffix),
               fontsize_row = top30_heatmap_params$fontsize_row,
               fontsize_col = top30_heatmap_params$fontsize_col,
               border_color = NA)
      dev.off()
      cat(sprintf("✓ logFC最大前30个基因无聚类PNG热图已保存: %s\n", png_file_top30_no_cluster))
    }, error = function(e) {
      cat(sprintf("✗ 无聚类热图生成失败: %s\n", e$message))
      dev.off()
    })
  }
  
} else if (top30_heatmap_params$enabled && nrow(heatmap_data) < top30_heatmap_params$top_n) {
  cat("\n========== 警告：前30个基因热图 ==========\n")
  cat(sprintf("警告：总基因数(%d)少于%d，跳过前30个基因热图生成\n", nrow(heatmap_data), top30_heatmap_params$top_n))
}

# ========================== 保存数据 ====================================
cat("\n========== 保存处理后的数据 ==========\n")

# 保存标准化数据
heatmap_data_df <- as.data.frame(heatmap_data)
heatmap_data_df <- cbind(GeneID = rownames(heatmap_data_df), heatmap_data_df)
write.csv(heatmap_data_df, file.path(output_dir, "normalized_data.csv"), row.names = FALSE)

# 保存原始表达数据
expr_data_df <- as.data.frame(expression_matrix)
expr_data_df <- cbind(GeneID = rownames(expr_data_df), expr_data_df)
write.csv(expr_data_df, file.path(output_dir, "expression_data_filtered.csv"), row.names = FALSE)

# 保存样本分组
write.csv(sample_groups, file.path(output_dir, "sample_groups.csv"), row.names = FALSE)

# 保存基因列表（实际使用的基因）
gene_list_used <- data.frame(
  GeneID = rownames(expression_matrix),
  Order = 1:nrow(expression_matrix),
  stringsAsFactors = FALSE
)
write.csv(gene_list_used, file.path(output_dir, "genes_used.csv"), row.names = FALSE)

# 保存logFC排序结果
write.csv(gene_logFC_sorted, file.path(output_dir, "genes_with_logFC_sorted.csv"), row.names = FALSE)

# 保存前30个logFC最大的基因列表
if (top30_heatmap_params$enabled && nrow(heatmap_data) >= top30_heatmap_params$top_n) {
  top30_genes_df <- top_genes_logFC
  top30_genes_df$Rank <- 1:nrow(top30_genes_df)
  write.csv(top30_genes_df, file.path(output_dir, "top30_logFC_genes_list.csv"), row.names = FALSE)
}

# ========================== 生成报告 ====================================
cat("\n========================================\n")
cat("热图可视化完成！\n")
cat("\n输出目录:", output_dir, "\n")
cat("\n输出文件:\n")
cat("  ├── gene_expression_heatmap.png\n")
cat("  ├── gene_expression_heatmap.pdf\n")
cat("  ├── gene_expression_heatmap_no_cluster.png\n")
if (top30_heatmap_params$enabled && nrow(heatmap_data) >= top30_heatmap_params$top_n) {
  cat("  ├── top30_logFC_genes_heatmap.png\n")
  cat("  ├── top30_logFC_genes_heatmap.pdf\n")
  if (top30_heatmap_params$cluster_rows || top30_heatmap_params$cluster_cols) {
    cat("  ├── top30_logFC_genes_heatmap_no_cluster.png\n")
  }
  cat("  ├── top30_logFC_genes_list.csv\n")
}
cat("  ├── normalized_data.csv\n")
cat("  ├── expression_data_filtered.csv\n")
cat("  ├── sample_groups.csv\n")
cat("  ├── genes_used.csv\n")
cat("  └── genes_with_logFC_sorted.csv\n")

cat("\n数据统计:\n")
cat(sprintf("  指定基因数: %d\n", length(genes_to_plot)))
cat(sprintf("  实际使用基因数: %d\n", nrow(expression_matrix)))
cat(sprintf("  样本数: %d\n", ncol(expression_matrix)))
cat("  分组分布:\n")
for(grp in names(table(sample_groups$Group))) {
  cat(sprintf("    %s: %d\n", grp, table(sample_groups$Group)[grp]))
}

cat("\n基因名称显示:\n")
cat(sprintf("  是否显示基因名: %s\n", ifelse(show_rownames, "是", "否")))
if (show_rownames) {
  cat(sprintf("  基因名字体大小: %d\n", fontsize_row))
}

if (top30_heatmap_params$enabled && nrow(heatmap_data) >= top30_heatmap_params$top_n) {
  cat("\nlogFC最大前30个基因热图自定义参数:\n")
  cat(sprintf("  logFC列名: %s\n", top30_heatmap_params$logFC_column))
  cat(sprintf("  使用绝对值: %s\n", ifelse(top30_heatmap_params$use_abs_logFC, "是", "否")))
  cat(sprintf("  选择基因数: %d\n", top30_heatmap_params$top_n))
  cat(sprintf("  热图宽度: %.1f 英寸\n", top30_heatmap_params$width))
  cat(sprintf("  热图高度: %.1f 英寸\n", top30_heatmap_params$height))
  cat(sprintf("  基因字体大小: %d\n", top30_heatmap_params$fontsize_row))
  cat(sprintf("  样本字体大小: %d\n", top30_heatmap_params$fontsize_col))
  cat(sprintf("  是否显示基因名: %s\n", ifelse(top30_heatmap_params$show_rownames, "是", "否")))
  cat(sprintf("  是否行聚类: %s\n", ifelse(top30_heatmap_params$cluster_rows, "是", "否")))
  cat(sprintf("  是否列聚类: %s\n", ifelse(top30_heatmap_params$cluster_cols, "是", "否")))
}

cat("\n========================================\n")