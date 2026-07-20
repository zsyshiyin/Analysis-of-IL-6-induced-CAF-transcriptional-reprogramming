# =============================================================================
# 生信分析脚本：差异蛋白与差异基因重合分析及表达相关性
# 功能：分析差异蛋白与差异基因的重合情况，绘制Venn图，并分析重合分子的表达相关性
# 输入：I6_vs_Con_DEGs.csv（差异蛋白数据）
#       I6-vs-Con.known.DEG.csv（差异基因数据）
# 输出：Venn图、重合分子列表、相关性散点图、统计分析结果
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/17-CAF-Sec-mRNA-venn-correlation"  # 请修改此路径
setwd(work_dir)

# 2. 创建输出文件夹
output_dir <- "protein_gene_overlap_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 3. 输入文件名设置
protein_file <- "I6_vs_Con_DEGs.csv"           # 差异蛋白数据文件
gene_file <- "I6-vs-Con.known.DEG.csv"         # 差异基因数据文件

# 4. 列名设置（可自定义）
# 差异蛋白文件列名
protein_id_col <- "GeneID"                  # 蛋白ID列名（请根据实际修改）
protein_name_col <- "GeneID"              # 蛋白名称列名（请根据实际修改）
protein_logFC_col <- "logFC"           # log2FC列名
protein_pvalue_col <- "PValue"                  # p值列名
protein_adj_pvalue_col <- "FDR"                # 校正后p值列名

# 差异基因文件列名
gene_id_col <- "GeneID"                        # 基因ID列名（请根据实际修改）
gene_name_col <- "Gene_name"                    # 基因名称列名（请根据实际修改）
gene_logFC_col <- "logFC"              # log2FC列名
gene_pvalue_col <- "Pvalue"                     # p值列名
gene_adj_pvalue_col <- "FDR"                   # 校正后p值列名

# 5. 匹配方式设置
# 用于匹配的列：可以选择使用ID匹配或名称匹配
use_id_for_match <- TRUE                        # TRUE: 使用ID匹配, FALSE: 使用名称匹配
protein_match_col <- protein_id_col             # 蛋白匹配列（如果use_id_for_match=TRUE）
gene_match_col <- gene_id_col                   # 基因匹配列（如果use_id_for_match=TRUE）

# 如果不使用ID匹配，则使用名称匹配
if (!use_id_for_match) {
  protein_match_col <- protein_name_col
  gene_match_col <- gene_name_col
}

# 6. 显著性阈值设置
protein_pvalue_threshold <- 0.05                # 蛋白显著p值阈值
protein_logFC_threshold <- 1                    # 蛋白显著log2FC阈值（绝对值）
gene_pvalue_threshold <- 0.05                   # 基因显著p值阈值
gene_logFC_threshold <- 1                       # 基因显著log2FC阈值（绝对值）

# 7. Venn图设置
venn_colors <- c("#E41A1C", "#F39B7F")           # 蛋白组和转录组的颜色
venn_alpha <- 0.5                                # 透明度
venn_title <- "Overlap of Differentially Expressed\nProteins and Genes"
venn_font_size <- 1.2                            # 字体大小

# 8. 相关性分析设置
cor_method <- "spearman"                         # 相关性方法: "spearman", "pearson"
cor_p_adjust_method <- "BH"                      # p值校正方法

# 9. 散点图设置（完全可自定义）
scatter_width <- 6                               # 单图宽度（英寸）
scatter_height <- 6                              # 单图高度（英寸）
point_pch <- 16                                  # 点形状
point_cex <- 1.5                                 # 点大小
point_color <- "#CA0E12"                         # 点的颜色
line_color <- "black"                            # 回归线颜色
line_lwd <- 3                                    # 回归线宽度
title_cex <- 1.2                                 # 标题字体大小
label_cex <- 1.2                                 # 轴标签字体大小
axis_cex <- 1.0                                  # 轴刻度字体大小

# 10. 输出文件设置
output_venn <- file.path(output_dir, "protein_gene_overlap_venn.pdf")
output_overlap_list <- file.path(output_dir, "overlap_molecules_list.csv")
output_correlation_results <- file.path(output_dir, "correlation_analysis_results.csv")
output_scatter_dir <- file.path(output_dir, "correlation_scatter_plots")
output_summary <- file.path(output_dir, "analysis_summary.txt")

# 创建散点图文件夹
if (!dir.exists(output_scatter_dir)) {
  dir.create(output_scatter_dir, recursive = TRUE)
  cat("创建散点图文件夹:", output_scatter_dir, "\n")
}

# ====================== 加载必要的包 ======================
cat("\n========== 加载必要的R包 ==========\n")

packages <- c("dplyr", "tidyr", "ggplot2", "VennDiagram", "ggpubr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取差异蛋白数据 ======================
cat("\n========== 读取差异蛋白数据 ==========\n")

protein_data <- read.csv(protein_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("差异蛋白数据维度: %d 行, %d 列\n", nrow(protein_data), ncol(protein_data)))
cat("列名:\n")
print(colnames(protein_data))

# ====================== 读取差异基因数据 ======================
cat("\n========== 读取差异基因数据 ==========\n")

gene_data <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("差异基因数据维度: %d 行, %d 列\n", nrow(gene_data), ncol(gene_data)))
cat("列名:\n")
print(colnames(gene_data))

# ====================== 数据预处理 ======================
cat("\n========== 数据预处理 ==========\n")

# 确保匹配列存在
if (!(protein_match_col %in% colnames(protein_data))) {
  stop(sprintf("错误：蛋白数据中找不到列 '%s'", protein_match_col))
}
if (!(gene_match_col %in% colnames(gene_data))) {
  stop(sprintf("错误：基因数据中找不到列 '%s'", gene_match_col))
}

# 标准化匹配列（转为大写，去除空格）
protein_data$Match_ID <- toupper(gsub("[[:space:]]", "", protein_data[[protein_match_col]]))
gene_data$Match_ID <- toupper(gsub("[[:space:]]", "", gene_data[[gene_match_col]]))

# 筛选显著差异的蛋白
if (protein_pvalue_col %in% colnames(protein_data)) {
  protein_significant <- protein_data %>%
    filter(!!sym(protein_pvalue_col) < protein_pvalue_threshold) %>%
    filter(abs(!!sym(protein_logFC_col)) >= protein_logFC_threshold)
} else {
  # 如果没有p值列，则使用所有数据
  protein_significant <- protein_data
  cat("警告：未找到p值列，使用所有蛋白数据\n")
}

# 筛选显著差异的基因
if (gene_pvalue_col %in% colnames(gene_data)) {
  gene_significant <- gene_data %>%
    filter(!!sym(gene_pvalue_col) < gene_pvalue_threshold) %>%
    filter(abs(!!sym(gene_logFC_col)) >= gene_logFC_threshold)
} else {
  # 如果没有p值列，则使用所有数据
  gene_significant <- gene_data
  cat("警告：未找到p值列，使用所有基因数据\n")
}

cat(sprintf("\n显著差异蛋白数: %d\n", nrow(protein_significant)))
cat(sprintf("显著差异基因数: %d\n", nrow(gene_significant)))

# ====================== 提取分子名称列表 ======================
cat("\n========== 提取分子名称列表 ==========\n")

protein_names <- unique(protein_significant$Match_ID)
gene_names <- unique(gene_significant$Match_ID)

cat(sprintf("蛋白唯一分子数: %d\n", length(protein_names)))
cat(sprintf("基因唯一分子数: %d\n", length(gene_names)))

# ====================== 分析重合情况 ======================
cat("\n========== 分析重合情况 ==========\n")

# 计算重合分子
overlap_molecules <- intersect(protein_names, gene_names)
protein_only <- setdiff(protein_names, gene_names)
gene_only <- setdiff(gene_names, protein_names)

cat(sprintf("\n重合分子数: %d\n", length(overlap_molecules)))
cat(sprintf("蛋白特有分子数: %d\n", length(protein_only)))
cat(sprintf("基因特有分子数: %d\n", length(gene_only)))

if (length(overlap_molecules) > 0) {
  cat("\n重合分子列表:\n")
  print(overlap_molecules)
}

# ====================== 创建重合分子数据框 ======================
cat("\n========== 创建重合分子数据框 ==========\n")

if (length(overlap_molecules) > 0) {
  # 提取重合分子的详细信息
  overlap_data <- data.frame(
    Match_ID = overlap_molecules,
    stringsAsFactors = FALSE
  )
  
  # 添加蛋白信息
  protein_overlap <- protein_significant %>%
    filter(Match_ID %in% overlap_molecules) %>%
    select(Match_ID, all_of(c(protein_id_col, protein_name_col, 
                              protein_logFC_col, protein_pvalue_col, 
                              protein_adj_pvalue_col)))
  
  # 添加基因信息
  gene_overlap <- gene_significant %>%
    filter(Match_ID %in% overlap_molecules) %>%
    select(Match_ID, all_of(c(gene_id_col, gene_name_col,
                              gene_logFC_col, gene_pvalue_col,
                              gene_adj_pvalue_col)))
  
  # 合并数据
  overlap_data <- overlap_data %>%
    left_join(protein_overlap, by = "Match_ID") %>%
    left_join(gene_overlap, by = "Match_ID")
  
  # 重命名列
  colnames(overlap_data)[colnames(overlap_data) == protein_logFC_col] <- "Protein_Log2FC"
  colnames(overlap_data)[colnames(overlap_data) == protein_pvalue_col] <- "Protein_Pvalue"
  colnames(overlap_data)[colnames(overlap_data) == protein_adj_pvalue_col] <- "Protein_AdjPvalue"
  colnames(overlap_data)[colnames(overlap_data) == gene_logFC_col] <- "Gene_Log2FC"
  colnames(overlap_data)[colnames(overlap_data) == gene_pvalue_col] <- "Gene_Pvalue"
  colnames(overlap_data)[colnames(overlap_data) == gene_adj_pvalue_col] <- "Gene_AdjPvalue"
  
  # 保存重合分子列表
  write.csv(overlap_data, file = output_overlap_list, row.names = FALSE)
  cat("\n重合分子列表已保存到:", output_overlap_list, "\n")
  
  # 显示重合分子统计
  cat("\n重合分子表达变化统计:\n")
  cat("蛋白表达变化:\n")
  print(summary(overlap_data$Protein_Log2FC))
  cat("\n基因表达变化:\n")
  print(summary(overlap_data$Gene_Log2FC))
}

# ====================== 绘制Venn图 ======================
cat("\n========== 绘制Venn图 ==========\n")

# 创建列表用于Venn图
venn_list <- list(
  Proteins = protein_names,
  Genes = gene_names
)

# 绘制Venn图
pdf(output_venn, width = 8, height = 8)

# 使用VennDiagram包绘制
grid.newpage()
venn_plot <- draw.pairwise.venn(
  area1 = length(protein_names),
  area2 = length(gene_names),
  cross.area = length(overlap_molecules),
  category = c("Proteins", "Genes"),
  fill = venn_colors,
  alpha = venn_alpha,
  lty = "blank",
  cex = venn_font_size,
  cat.cex = venn_font_size,
  cat.pos = c(0, 0),
  cat.dist = c(0.05, 0.05),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
  ext.text = FALSE,
  scaled = FALSE
)

# 添加标题
grid.text(venn_title, x = 0.5, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))

dev.off()
cat("Venn图已保存到:", output_venn, "\n")

# ====================== 相关性分析 ======================
cat("\n========== 相关性分析 ==========\n")

if (length(overlap_molecules) > 0) {
  # 准备用于相关性分析的数据
  # 注意：这里假设蛋白和基因的表达值已经包含在原始数据中
  # 如果没有表达值，需要单独读取表达矩阵文件
  
  # 检查是否有表达值列
  if (("Protein_Log2FC" %in% colnames(overlap_data)) && 
      ("Gene_Log2FC" %in% colnames(overlap_data))) {
    
    # 使用log2FC进行相关性分析
    x <- overlap_data$Protein_Log2FC
    y <- overlap_data$Gene_Log2FC
    
    # 移除缺失值
    valid_idx <- complete.cases(x, y)
    x_valid <- x[valid_idx]
    y_valid <- y[valid_idx]
    valid_molecules <- overlap_data$Match_ID[valid_idx]
    
    cat(sprintf("\n用于相关性分析的分子数: %d\n", length(valid_molecules)))
    
    if (length(valid_molecules) >= 3) {
      # 计算相关性
      cor_test <- cor.test(x_valid, y_valid, method = cor_method)
      
      cat("\n相关性分析结果:\n")
      cat(sprintf("  相关性方法: %s\n", cor_method))
      cat(sprintf("  相关系数: %.4f\n", cor_test$estimate))
      cat(sprintf("  p值: %.4e\n", cor_test$p.value))
      
      # 创建相关性结果数据框
      correlation_result <- data.frame(
        Analysis = "Protein_Log2FC vs Gene_Log2FC",
        Correlation_Method = cor_method,
        Correlation_Coefficient = cor_test$estimate,
        P_value = cor_test$p.value,
        Sample_Size = length(valid_molecules),
        stringsAsFactors = FALSE
      )
      
      # 计算每个分子的表达一致性
      overlap_data$Expression_Consistency <- ifelse(
        overlap_data$Protein_Log2FC * overlap_data$Gene_Log2FC > 0,
        "Same Direction",
        "Opposite Direction"
      )
      
      # 添加显著性信息
      correlation_result$Significance <- ifelse(cor_test$p.value < 0.001, "***",
                                                ifelse(cor_test$p.value < 0.01, "**",
                                                       ifelse(cor_test$p.value < 0.05, "*", "ns")))
      
      # 保存相关性结果
      write.csv(correlation_result, file = output_correlation_results, row.names = FALSE)
      cat("\n相关性分析结果已保存到:", output_correlation_results, "\n")
      
      # ====================== 绘制相关性散点图 ======================
      cat("\n========== 绘制相关性散点图 ==========\n")
      
      # 绘制整体相关性散点图
      scatter_plot_file <- file.path(output_scatter_dir, "overall_correlation_scatter.pdf")
      
      pdf(scatter_plot_file, width = scatter_width, height = scatter_height)
      
      # 准备绘图数据
      plot_data <- data.frame(
        Protein_Log2FC = x_valid,
        Gene_Log2FC = y_valid,
        Molecule = valid_molecules
      )
      
      # 拟合线性模型
      z <- lm(Gene_Log2FC ~ Protein_Log2FC, data = plot_data)
      
      # 格式化标题
      if (cor_test$p.value < 0.001) {
        p_display <- "p < 0.001"
      } else {
        p_display <- sprintf("p = %.3f", cor_test$p.value)
      }
      
      main_title <- sprintf("Protein vs Gene Expression Correlation\nSpearman's r = %.3f%s, %s", 
                            cor_test$estimate, 
                            correlation_result$Significance,
                            p_display)
      
      # 绘制散点图
      plot(plot_data$Protein_Log2FC, plot_data$Gene_Log2FC,
           type = "p",
           pch = point_pch,
           cex = point_cex,
           col = point_color,
           main = main_title,
           xlab = "Protein Log2FC",
           ylab = "Gene Log2FC",
           cex.main = title_cex,
           cex.lab = label_cex,
           cex.axis = axis_cex)
      
      # 添加回归线
      lines(plot_data$Protein_Log2FC[order(plot_data$Protein_Log2FC)], 
            fitted(z)[order(plot_data$Protein_Log2FC)], 
            col = line_color, 
            lwd = line_lwd)
      
      # 添加参考线
      abline(h = 0, v = 0, col = "gray", lty = 2, lwd = 1)
      
      dev.off()
      cat("整体相关性散点图已保存到:", scatter_plot_file, "\n")
      
      # 为每个重合分子绘制单独的散点图（需要原始表达数据）
      # 注意：这部分需要每个分子的原始表达值，如果数据可用，可以添加
      
    } else {
      cat("\n警告：有效数据点不足3个，无法进行相关性分析\n")
    }
    
  } else {
    cat("\n警告：未找到表达值列（Protein_Log2FC和Gene_Log2FC），无法进行相关性分析\n")
    cat("如果需要分析表达相关性，请确保数据中包含表达值列\n")
  }
  
} else {
  cat("\n没有重合分子，跳过相关性分析\n")
}

# ====================== 生成分析报告 ======================
cat("\n========== 生成分析报告 ==========\n")

sink(output_summary)

cat("========================================\n")
cat("   差异蛋白与差异基因重合分析报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、数据概况\n")
cat("------------\n")
cat(sprintf("差异蛋白总数: %d\n", nrow(protein_data)))
cat(sprintf("显著差异蛋白数 (p<%.2f, |log2FC|>=%.1f): %d\n", 
            protein_pvalue_threshold, protein_logFC_threshold, nrow(protein_significant)))
cat(sprintf("差异基因总数: %d\n", nrow(gene_data)))
cat(sprintf("显著差异基因数 (p<%.2f, |log2FC|>=%.1f): %d\n", 
            gene_pvalue_threshold, gene_logFC_threshold, nrow(gene_significant)))

cat("\n二、重合分析结果\n")
cat("----------------\n")
cat(sprintf("蛋白唯一分子数: %d\n", length(protein_names)))
cat(sprintf("基因唯一分子数: %d\n", length(gene_names)))
cat(sprintf("重合分子数: %d\n", length(overlap_molecules)))
cat(sprintf("重合率 (占蛋白): %.1f%%\n", 100 * length(overlap_molecules) / length(protein_names)))
cat(sprintf("重合率 (占基因): %.1f%%\n", 100 * length(overlap_molecules) / length(gene_names)))

if (length(overlap_molecules) > 0) {
  cat("\n三、重合分子列表\n")
  cat("----------------\n")
  for (i in 1:min(20, length(overlap_molecules))) {
    cat(sprintf("  %s\n", overlap_molecules[i]))
  }
  if (length(overlap_molecules) > 20) {
    cat(sprintf("  ... 还有 %d 个分子\n", length(overlap_molecules) - 20))
  }
  
  cat("\n四、相关性分析结果\n")
  cat("------------------\n")
  if (exists("cor_test")) {
    cat(sprintf("相关性方法: %s\n", cor_method))
    cat(sprintf("相关系数: %.4f\n", cor_test$estimate))
    cat(sprintf("p值: %.4e\n", cor_test$p.value))
    cat(sprintf("样本量: %d\n", length(valid_molecules)))
    
    # 表达一致性统计
    if ("Expression_Consistency" %in% colnames(overlap_data)) {
      consistency_table <- table(overlap_data$Expression_Consistency, useNA = "ifany")
      cat("\n表达方向一致性:\n")
      print(consistency_table)
    }
  } else {
    cat("未进行相关性分析\n")
  }
}

cat("\n五、输出文件列表\n")
cat("----------------\n")
cat("1.", basename(output_venn), "- Venn图\n")
cat("2.", basename(output_overlap_list), "- 重合分子列表\n")
if (file.exists(output_correlation_results)) {
  cat("3.", basename(output_correlation_results), "- 相关性分析结果\n")
}
cat("4.", basename(output_scatter_dir), "/ - 相关性散点图文件夹\n")
cat("5.", basename(output_summary), "- 分析报告\n")

cat("\n========================================\n")
cat("                完成\n")
cat("========================================\n")

sink()

cat("\n分析报告已保存到:", output_summary, "\n")
cat("\n=================== 分析完成 ===================\n")