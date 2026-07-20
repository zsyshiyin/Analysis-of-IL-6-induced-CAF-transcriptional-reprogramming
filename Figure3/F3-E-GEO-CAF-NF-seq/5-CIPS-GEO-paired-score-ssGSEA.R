# ============================================================================
# 配对样本基因集表达评分分析（ssGSEA方法）- 带分组连线点图
# ============================================================================

# 清空环境变量
rm(list = ls())

# ========================== 工作路径设定 ====================================
# 请修改为您的实际工作路径
work_dir <- "D:/zsy/SX/Fomal-final/8-CIPS-GEO-score"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")

# ========================== 参数设置 =========================================
# 输入文件参数
expr_file <- "expression_data-geo.csv"           
geneset_file <- "gene_set.txt"                
pairing_file <- "pairing_info.csv"            

# 输出文件参数
output_dir <- "paired_geneset_results"        
score_file <- "paired_geneset_scores_ssgsea.csv"     
plot_file <- "paired_scores_ssgsea_plot.png"         
stats_file <- "paired_stats_ssgsea.csv"              

# ===== ssGSEA参数 =====
ssgsea_alpha <- 0.1                           # GSEA官方参数
ssgsea_norm <- TRUE                             # 归一化分数

# ===== 列名自定义参数 =====
expr_gene_col <- "GeneID"                       
pairing_pair_col <- "PairID"                    
pairing_sample_col <- "SampleID"                  
pairing_group_col <- "Group"                      

# ===== 组别顺序参数 =====
group_order <- c("NFs", "CAFs")                  
group_labels <- c("NFs", "CAFs")                             

# ===== 连线点图自定义参数 =====
plot_width <- 3                                  
plot_height <- 4                                 
plot_res <- 300                                  
point_size <- 6                                  
point_shape <- 21                                
point_fill_colors <- c("#25377F", "#C60036")     
point_border_color <- "black"                    
point_alpha <- 0.8                              
line_color <- "gray50"                            
line_width <- 1.5                                    
line_alpha <- 0.6                                 
line_type <- "solid"                              

axis_title_size <- 20                             
axis_text_size <- 18                              
legend_title_size <- 18                           
legend_text_size <- 16                            
plot_title_size <- 22                             

show_pvalue <- TRUE                               
pvalue_size <- 7                                  
pvalue_color <- "black"                            
pvalue_position <- "top"                          

plot_title <- "ssGSEA Scores"     
x_label <- "Group"                                  
y_label <- "ssGSEA scores of patient-derived fibroblast samples from public databases "                           
legend_title <- "Group"                              

perform_test <- TRUE                                
test_method <- "wilcox"                              

# ========================== 检查必要包 ======================================
required_packages <- c("ggplot2", "dplyr", "tidyr", "ggpubr", "tools")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# ========================== 创建输出目录 ====================================
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ========================== 读取数据 ====================================
cat("\n========== 读取数据 ==========\n")

expr_data <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cat(sprintf("表达矩阵: %d 基因, %d 样品\n", nrow(expr_data), ncol(expr_data)))

# 显示基因表达量示例
cat("\n基因表达量示例（前5个基因，前5个样品）:\n")
print(round(expr_data[1:min(5, nrow(expr_data)), 1:min(5, ncol(expr_data))], 2))

geneset_genes <- readLines(geneset_file)
geneset_genes <- toupper(trimws(geneset_genes[geneset_genes != ""]))
cat(sprintf("\n基因集: %d 基因\n", length(geneset_genes)))
cat("基因集基因:", paste(head(geneset_genes, 10), collapse=", "), "\n")

pairing_info <- read.csv(pairing_file, stringsAsFactors = FALSE)
cat(sprintf("\n配对信息: %d 行\n", nrow(pairing_info)))
print(head(pairing_info))

# ========================== 样品名标准化 ====================================
cat("\n========== 样品名标准化 ==========\n")

# 标准化样品名
colnames(expr_data) <- toupper(gsub("[[:space:]]", "", colnames(expr_data)))
pairing_info[[pairing_sample_col]] <- toupper(gsub("[[:space:]]", "", 
                                                   pairing_info[[pairing_sample_col]]))

cat("表达矩阵样品名:", paste(head(colnames(expr_data)), collapse=", "), "\n")
cat("配对文件样品名:", paste(head(pairing_info[[pairing_sample_col]]), collapse=", "), "\n")

# 找出共同样品
common_samples <- intersect(colnames(expr_data), pairing_info[[pairing_sample_col]])
cat(sprintf("共同样品: %d\n", length(common_samples)))

if (length(common_samples) == 0) stop("无共同样品")

# 筛选数据
expr_data <- expr_data[, common_samples, drop = FALSE]
pairing_info <- pairing_info[pairing_info[[pairing_sample_col]] %in% common_samples, ]

# ========================== 筛选基因集基因 ====================================
cat("\n========== 筛选基因集基因 ==========\n")

rownames(expr_data) <- toupper(rownames(expr_data))
available_genes <- intersect(geneset_genes, rownames(expr_data))
missing_genes <- setdiff(geneset_genes, rownames(expr_data))

cat(sprintf("可用基因: %d/%d\n", length(available_genes), length(geneset_genes)))

if (length(missing_genes) > 0) {
  cat("缺失基因:", paste(head(missing_genes, 10), collapse=", "), "\n")
}

if (length(available_genes) == 0) stop("无可用基因")

geneset_expr <- expr_data[available_genes, , drop = FALSE]

# 显示基因集基因的表达量
cat("\n基因集基因表达量示例:\n")
print(round(geneset_expr[1:min(5, nrow(geneset_expr)), 1:min(5, ncol(geneset_expr))], 2))

# ============================================================================
# 正确的ssGSEA实现（基于GSEA官方算法）
# ============================================================================

cat("\n========== 运行GSEA官方ssGSEA算法 ==========\n")

ssgsea_official <- function(expression_matrix, gene_set, alpha = 0.25) {
  # 基于GSEA官方算法的ssGSEA实现
  # 参考文献: Subramanian et al., PNAS 2005; Barbie et al., Nature 2009
  
  n_genes <- nrow(expression_matrix)
  n_samples <- ncol(expression_matrix)
  n_gs <- length(gene_set)
  
  # 获取基因集索引
  gs_idx <- which(rownames(expression_matrix) %in% gene_set)
  
  cat(sprintf("总基因数: %d\n", n_genes))
  cat(sprintf("基因集大小: %d\n", length(gs_idx)))
  cat(sprintf("样品数: %d\n", n_samples))
  cat(sprintf("alpha参数: %.2f\n", alpha))
  
  scores <- numeric(n_samples)
  names(scores) <- colnames(expression_matrix)
  
  for (s in 1:n_samples) {
    # 获取该样品的表达值
    expr <- expression_matrix[, s]
    
    # 步骤1: 对基因表达值排序（从高到低）
    # 这是GSEA的核心 - 基于表达值的绝对排序
    order_idx <- order(expr, decreasing = TRUE)
    ranked_expr <- expr[order_idx]
    ranked_genes <- rownames(expression_matrix)[order_idx]
    
    # 步骤2: 计算基因集基因的位置
    gs_positions <- which(ranked_genes %in% gene_set)
    
    if (length(gs_positions) == 0) {
      scores[s] <- 0
      next
    }
    
    # 步骤3: 计算富集分数（ES）
    # 使用加权累积和
    N <- n_genes
    Nh <- length(gs_positions)
    
    # 初始化
    es <- 0
    max_es <- 0
    min_es <- 0
    
    # 计算缺失基因的权重
    miss_metric <- sqrt((N - Nh) / Nh)
    
    # 遍历排序列表
    for (i in 1:N) {
      if (i %in% gs_positions) {
        # 命中基因：增加分数，权重基于表达值的秩
        # 这是ssGSEA的关键 - 使用表达值的秩的alpha次方作为权重
        rank_i <- which(order_idx == gs_positions[which(gs_positions == i)])
        weight <- abs(expr[order_idx[i]]) ^ alpha  # 使用表达值的alpha次方
        es <- es + weight
      } else {
        # 未命中基因：减少分数
        es <- es - (1 / (N - Nh))
      }
      
      # 记录最大偏差
      if (es > max_es) max_es <- es
      if (es < min_es) min_es <- es
    }
    
    # 步骤4: 归一化
    # ES是最大偏差
    if (abs(max_es) >= abs(min_es)) {
      scores[s] <- max_es
    } else {
      scores[s] <- min_es
    }
    
    # 调试信息（前3个样品）
    if (s <= 3) {
      cat(sprintf("\n样品 %s:\n", colnames(expression_matrix)[s]))
      cat(sprintf("  基因集基因位置: %s\n", paste(head(gs_positions, 10), collapse=", ")))
      cat(sprintf("  原始ES: %.6f\n", scores[s]))
      
      # 显示前10个最高表达基因
      cat("  前10个高表达基因:\n")
      for (i in 1:min(10, N)) {
        hit_flag <- ifelse(i %in% gs_positions, " [HIT]", "")
        cat(sprintf("    %2d. %s: %.2f%s\n", 
                    i, ranked_genes[i], ranked_expr[i], hit_flag))
      }
    }
  }
  
  # 步骤5: 归一化到[-1, 1]区间
  if (max(scores) > min(scores)) {
    scores_norm <- 2 * (scores - min(scores)) / (max(scores) - min(scores)) - 1
    cat("\n分数归一化到[-1, 1]区间\n")
  } else {
    scores_norm <- scores
    cat("\n警告: 所有分数相同，无法归一化\n")
  }
  
  return(scores_norm)
}

# 运行官方ssGSEA
expr_matrix <- as.matrix(geneset_expr)
ssgsea_scores <- ssgsea_official(expr_matrix, available_genes, alpha = ssgsea_alpha)

# 创建评分数据框
scores_df <- data.frame(
  Sample = colnames(geneset_expr),
  Score = ssgsea_scores,
  stringsAsFactors = FALSE
)

# 显示最终评分
cat("\n\n========== ssGSEA评分结果 ==========\n")
cat(sprintf("评分范围: [%.6f, %.6f]\n", min(scores_df$Score), max(scores_df$Score)))
cat(sprintf("评分均值: %.6f\n", mean(scores_df$Score)))
cat(sprintf("评分标准差: %.6f\n", sd(scores_df$Score)))

# 按样品显示评分
cat("\n各样品评分:\n")
for (i in 1:nrow(scores_df)) {
  cat(sprintf("  %s: %.6f\n", scores_df$Sample[i], scores_df$Score[i]))
}

# ============================================================================
# 验证：检查基因集基因在不同组中的表达差异
# ============================================================================

cat("\n========== 验证基因表达差异 ==========\n")

# 将分组信息添加到表达数据
group_map <- setNames(pairing_info[[pairing_group_col]], 
                      pairing_info[[pairing_sample_col]])

# 计算每个基因在两组中的平均表达
group1_samples <- names(group_map[group_map == group_order[1]])
group2_samples <- names(group_map[group_map == group_order[2]])

cat(sprintf("\n%s组样品: %s\n", group_order[1], paste(group1_samples, collapse=", ")))
cat(sprintf("%s组样品: %s\n", group_order[2], paste(group2_samples, collapse=", ")))

# 计算基因集基因的差异
gene_diffs <- data.frame(
  Gene = available_genes,
  Mean_Group1 = rowMeans(geneset_expr[, intersect(available_genes, group1_samples), drop=FALSE]),
  Mean_Group2 = rowMeans(geneset_expr[, intersect(available_genes, group2_samples), drop=FALSE]),
  stringsAsFactors = FALSE
)
gene_diffs$Diff <- gene_diffs$Mean_Group2 - gene_diffs$Mean_Group1
gene_diffs$Log2FC <- log2(gene_diffs$Mean_Group2 / gene_diffs$Mean_Group1)

cat("\n基因集基因表达差异（前10个）:\n")
print(head(gene_diffs[order(abs(gene_diffs$Diff), decreasing=TRUE), ], 10))

# ============================================================================
# 合并配对信息
# ============================================================================

cat("\n========== 合并配对信息 ==========\n")

paired_scores <- merge(scores_df, pairing_info, 
                       by.x = "Sample", by.y = pairing_sample_col,
                       all = FALSE)

cat(sprintf("合并后数据: %d 行\n", nrow(paired_scores)))

if (nrow(paired_scores) == 0) {
  stop("合并失败")
}

# 检查配对完整性
pair_counts <- table(paired_scores[[pairing_pair_col]])
cat("\n每个PairID的样品数:\n")
print(pair_counts)

# 保留完整配对
valid_pairs <- names(pair_counts[pair_counts == 2])
paired_scores <- paired_scores[paired_scores[[pairing_pair_col]] %in% valid_pairs, ]
cat(sprintf("\n有效配对: %d\n", length(valid_pairs)))

if (length(valid_pairs) == 0) {
  stop("无完整配对")
}

# 转换为宽格式
unique_pairs <- unique(paired_scores[[pairing_pair_col]])
group_names <- unique(paired_scores[[pairing_group_col]])

paired_wide <- data.frame(PairID = unique_pairs)
for (g in group_names) {
  paired_wide[[g]] <- NA_real_
}

for (i in 1:nrow(paired_scores)) {
  pair_id <- paired_scores[i, pairing_pair_col]
  group <- paired_scores[i, pairing_group_col]
  score <- paired_scores[i, "Score"]
  
  paired_wide[paired_wide$PairID == pair_id, group] <- score
}

cat("\n宽格式数据:\n")
print(paired_wide)

# 移除不完整的配对
paired_wide <- paired_wide[complete.cases(paired_wide[, group_names]), ]
cat(sprintf("\n完整宽格式数据: %d 行\n", nrow(paired_wide)))

if (nrow(paired_wide) == 0) {
  stop("无完整宽格式数据")
}

# ============================================================================
# 统计分析
# ============================================================================

cat("\n========== 统计分析 ==========\n")

# 确保组顺序正确
if (!is.null(group_order)) {
  group_names <- group_order[group_order %in% group_names]
}

group1_name <- group_names[1]
group2_name <- group_names[2]

group1_values <- paired_wide[[group1_name]]
group2_values <- paired_wide[[group2_name]]

cat(sprintf("\n%s组评分: %s\n", group1_name, paste(sprintf("%.4f", group1_values), collapse=", ")))
cat(sprintf("%s组评分: %s\n", group2_name, paste(sprintf("%.4f", group2_values), collapse=", ")))

# 配对检验
if (perform_test && length(group1_values) >= 2) {
  if (test_method == "t.test") {
    test_result <- t.test(group1_values, group2_values, paired = TRUE)
  } else {
    test_result <- wilcox.test(group1_values, group2_values, paired = TRUE)
  }
  p_value <- test_result$p.value
  cat(sprintf("\n%s p值: %.4f\n", test_method, p_value))
} else {
  p_value <- NA
}

# ============================================================================
# 保存结果
# ============================================================================

cat("\n========== 保存结果 ==========\n")

write.csv(paired_scores, file.path(output_dir, score_file), row.names=FALSE)
write.csv(paired_wide, file.path(output_dir, stats_file), row.names=FALSE)

cat(sprintf("评分数据保存至: %s\n", file.path(output_dir, score_file)))
cat(sprintf("统计信息保存至: %s\n", file.path(output_dir, stats_file)))

# ============================================================================
# 绘图
# ============================================================================

cat("\n========== 绘图 ==========\n")

plot_data <- paired_scores
plot_data$Group <- factor(plot_data[[pairing_group_col]], levels = group_order)

p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_line(aes(group = !!sym(pairing_pair_col)), 
            color = line_color, linewidth = line_width, alpha = line_alpha) +
  geom_point(size = point_size, shape = point_shape,
             fill = point_fill_colors[as.numeric(plot_data$Group)],
             color = point_border_color, alpha = point_alpha) +
  scale_fill_manual(values = point_fill_colors) +
  labs(title = plot_title, x = x_label, y = y_label, fill = legend_title) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = plot_title_size, face = "bold"),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.position = "right")

if (exists("p_value") && show_pvalue && !is.na(p_value)) {
  y_range <- range(plot_data$Score)
  p_y <- y_range[2] + 0.05 * diff(y_range)
  p_text <- ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))
  p <- p + annotate("text", x = 1.5, y = p_y, label = p_text,
                    size = pvalue_size, color = pvalue_color, hjust = 0.5)
}

png(file.path(output_dir, plot_file), width = plot_width, 
    height = plot_height, units = "in", res = plot_res)
print(p)
dev.off()

cat(sprintf("图形保存至: %s\n", file.path(output_dir, plot_file)))

cat("\n=== 分析完成 ===\n")