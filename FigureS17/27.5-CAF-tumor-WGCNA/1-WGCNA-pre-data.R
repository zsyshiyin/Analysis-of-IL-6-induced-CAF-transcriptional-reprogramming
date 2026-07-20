# ============================================================================
# WGCNA分析 - 模块1：数据准备和预处理（自动适应样品组数）
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
# 输入文件
caf_expression_file <- "Tumor-mRNA.csv"  # CAF转录组表达矩阵

# 分组识别规则（根据列名前缀）
# 请根据您的实际数据修改，只保留实际存在的分组
# 注意：不要有空元素！
group_mapping <- list(
  "DMEM" = "DMEM",    # 如果没有DMEM组，请注释掉
  "Con" = "Control",
  "I6" = "I6", 
  "R2" = "R2"        # 如果没有R2组，请注释掉
)

# 过滤参数
min_expr <- 1                    # 最小表达量
min_sample_frac <- 0.5           # 至少多少比例的样品表达
min_var_percentile <- 0.25       # 保留变异度前25%的基因

# 输出文件
output_rdata <- file.path(output_dir, "01_data_prepared.RData")
output_group_info <- file.path(output_dir, "sample_groups.csv")

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

required_packages <- c("WGCNA", "dplyr", "tidyr", "ggplot2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "WGCNA") {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("WGCNA")
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
  }
  cat(sprintf("✓ %s\n", pkg))
}

# 允许WGCNA使用多线程
enableWGCNAThreads()

# ========================== 读取CAF表达数据 ====================================
cat("\n========== 读取CAF表达数据 ==========\n")

if (!file.exists(caf_expression_file)) {
  stop("错误：未找到CAF表达矩阵文件")
}

# 读取数据
caf_data <- read.csv(caf_expression_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("原始数据维度: %d 基因 × %d 列\n", nrow(caf_data), ncol(caf_data)))
cat("列名:", paste(head(colnames(caf_data), 10), collapse = ", "), "\n")

# 设置行名
rownames(caf_data) <- caf_data[, 1]
caf_data <- caf_data[, -1]

# 转换为数值矩阵
expr_matrix <- as.matrix(caf_data)
mode(expr_matrix) <- "numeric"
expr_matrix[is.na(expr_matrix)] <- 0

cat(sprintf("表达矩阵维度: %d 基因 × %d 样品\n", nrow(expr_matrix), ncol(expr_matrix)))

# ========================== 自动识别样品分组 ====================================
cat("\n========== 自动识别样品分组 ==========\n")

sample_names <- colnames(expr_matrix)
sample_groups <- data.frame(
  SampleID = sample_names,
  Group = NA,
  stringsAsFactors = FALSE
)

# 根据列名前缀识别分组
for (prefix in names(group_mapping)) {
  matching_samples <- grep(paste0("^", prefix), sample_names, value = TRUE)
  if (length(matching_samples) > 0) {
    sample_groups$Group[sample_groups$SampleID %in% matching_samples] <- group_mapping[[prefix]]
  }
}

# 检查是否有未识别的样品
unassigned <- is.na(sample_groups$Group)
if (any(unassigned)) {
  cat("\n警告：以下样品未能识别分组:\n")
  cat(paste(sample_groups$SampleID[unassigned], collapse = ", "), "\n")
  cat("请检查group_mapping设置\n")
  
  # 移除未识别的样品
  sample_groups <- sample_groups[!unassigned, ]
  expr_matrix <- expr_matrix[, sample_groups$SampleID, drop = FALSE]
}

# 获取分组信息
groups <- unique(sample_groups$Group)
n_groups <- length(groups)
cat(sprintf("\n识别到 %d 个分组:\n", n_groups))
for (grp in groups) {
  n_samples <- sum(sample_groups$Group == grp)
  cat(sprintf("  %s: %d 个样品\n", grp, n_samples))
}

# 保存样品分组信息
write.csv(sample_groups, output_group_info, row.names = FALSE)
cat(sprintf("样品分组信息已保存: %s\n", output_group_info))

# ========================== 数据过滤 ====================================
cat("\n========== 数据过滤 ==========\n")

# 1. 过滤低表达基因
cat("过滤低表达基因...\n")
keep_by_expr <- rowMeans(expr_matrix) >= min_expr
cat(sprintf("  过滤前基因数: %d\n", nrow(expr_matrix)))
expr_matrix <- expr_matrix[keep_by_expr, ]
cat(sprintf("  过滤后基因数: %d\n", nrow(expr_matrix)))

# 2. 过滤缺失值过多的基因
cat("\n过滤缺失值过多的基因...\n")
na_frac <- apply(expr_matrix, 1, function(x) sum(is.na(x)) / length(x))
keep_by_na <- na_frac < (1 - min_sample_frac)
expr_matrix <- expr_matrix[keep_by_na, ]
cat(sprintf("  过滤后基因数: %d\n", nrow(expr_matrix)))

# 3. 保留变异度高的基因
cat("\n保留变异度高的基因...\n")
gene_var <- apply(expr_matrix, 1, var, na.rm = TRUE)
var_threshold <- quantile(gene_var, min_var_percentile, na.rm = TRUE)
keep_by_var <- gene_var >= var_threshold
expr_matrix <- expr_matrix[keep_by_var, ]
cat(sprintf("  保留变异度前 %.0f%% 的基因\n", (1 - min_var_percentile) * 100))
cat(sprintf("  过滤后基因数: %d\n", nrow(expr_matrix)))

# ========================== 数据标准化 ====================================
cat("\n========== 数据标准化 ==========\n")

# 对每个基因进行Z-score标准化
expr_matrix_norm <- t(scale(t(expr_matrix)))
expr_matrix_norm[is.na(expr_matrix_norm)] <- 0

cat("数据标准化完成\n")

# ========================== 保存预处理结果 ====================================
cat("\n========== 保存预处理结果 ==========\n")

save(expr_matrix, expr_matrix_norm, sample_groups, groups, n_groups, 
     file = output_rdata)

cat(sprintf("预处理数据已保存: %s\n", output_rdata))
cat(sprintf("  基因数: %d\n", nrow(expr_matrix)))
cat(sprintf("  样品数: %d\n", ncol(expr_matrix)))
cat(sprintf("  分组数: %d\n", n_groups))

cat("\n模块1完成！\n")