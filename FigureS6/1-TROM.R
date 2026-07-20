# ============================================================================
# TROM方法分析三个样品组之间的转录谱相似度
# ============================================================================

# 清空环境变量
rm(list = ls())

# ========================== 工作路径设定 ====================================
# 请修改为您的实际工作路径
work_dir <- "D:/zsy/SX/Fomal-Xin/1-TROM"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")


# ========================== 参数设置（所有可自定义参数集中在此）=============
# 输入文件参数
expr_file <- "gene_expression.csv"        # 基因表达矩阵文件
group_file <- "sample_groups.csv"         # 样本分组信息文件
custom_gene_file <- NULL                   # 自定义基因集文件（NULL表示不使用）

# 输出文件参数
output_dir <- "TROM_results"               # 输出文件夹
trom_score_file <- "trom_scores.csv"       # TROM得分矩阵输出文件
heatmap_file <- "trom_heatmap.png"         # 热图输出文件

# TROM分析参数
z_threshold <- NULL                         # Z-score阈值（NULL表示自动选择）
use_orthologs <- FALSE                      # 是否仅使用直系同源基因
force_reinstall_trom <- FALSE                # 是否强制重新安装TROM
auto_select_z <- TRUE                        # 是否自动选择Z-score阈值

# 数据标准化参数
perform_zscore <- TRUE                       # 是否手动进行Z-score标准化
zscore_method <- "row"                        # "row"按基因标准化，"col"按样品标准化

# 归一化参数
normalize_similarity <- TRUE                  # 是否归一化相似度矩阵到[0,1]区间
similarity_min <- NULL                        # 自定义归一化最小值（NULL则自动）
similarity_max <- NULL                        # 自定义归一化最大值（NULL则自动）

# 调试参数
verbose_trom <- TRUE                          # 显示TROM详细过程
validate_matrix <- TRUE                        # 是否验证矩阵合理性

# 热图自定义参数 - 基础
heatmap_width <- 10                          # 热图宽度（英寸）
heatmap_height <- 8                          # 热图高度（英寸）
heatmap_res <- 300                           # 热图分辨率（DPI）
heatmap_color_low <- "#25377F"                   # 低值颜色
heatmap_color_mid <- "white"                  # 中值颜色
heatmap_color_high <- "#C60036"                    # 高值颜色
heatmap_fontsize_row <- 12                     # 行标签字体大小
heatmap_fontsize_col <- 12                     # 列标签字体大小
heatmap_main_title <- "TROM Similarity Between Sample Groups"  # 主标题
heatmap_xlab <- "Sample Groups"                # X轴标签
heatmap_ylab <- "Sample Groups"                # Y轴标签
heatmap_key_label <- "Similarity Score"        # 图例标签

# ===== 新增：热图数值显示参数 =====
# 样品级别热图数值显示
sample_display_numbers <- FALSE                # 样品级别热图是否显示数值（通常不显示，因为格子太多）
sample_number_color <- "white"                  # 样品级别热图数值颜色
sample_number_fontsize <- 16                    # 样品级别热图数值大小

# 组间平均热图数值显示
group_display_numbers <- TRUE                   # 组间平均热图是否显示数值
group_number_color <- "white"                    # 组间平均热图数值颜色
group_number_fontsize <- 16                     # 组间平均热图数值大小
group_number_format <- "%.3f"                    # 组间平均热图数值格式

# 高级数值颜色控制
use_conditional_colors <- FALSE                  # 是否使用条件颜色（根据数值改变文字颜色）
conditional_thresholds <- c(0.7, 0.4)            # 条件颜色阈值（高/中/低）
conditional_colors <- c("white", "black", "darkblue")  # 对应的颜色

# 列名自定义参数（与输入文件匹配）
expr_gene_col <- "GeneID"                      # 表达矩阵中基因ID的列名
group_sample_col <- "SampleID"                  # 分组文件中样品名的列名
group_group_col <- "Group"                       # 分组文件中组别的列名

# 数据预处理参数
remove_non_numeric <- TRUE                      # 是否自动移除非数值列
convert_to_numeric <- TRUE                       # 是否尝试将字符型转换为数值型
na_threshold <- 0.2                              # 缺失值阈值，超过此比例的行将被移除
decimal_separator <- "."                          # 小数点分隔符
thousand_separator <- ""                          # 千位分隔符
remove_units <- TRUE                              # 是否移除单位符号

# ========================== 检查必要包 ======================================
cat("\n========== 检查必要包 ==========\n")

required_packages <- c("pheatmap", "RColorBrewer", "openxlsx")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat(sprintf("安装CRAN包: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/", dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    cat(sprintf("✓ %s 已安装\n", pkg))
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("安装BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# ========================== 检查TROM包 ======================================
cat("\n========== 检查TROM包 ==========\n")

trom_available <- require("TROM", character.only = TRUE, quietly = TRUE)
use_alternative <- FALSE

if (trom_available && !force_reinstall_trom) {
  cat("✓ TROM包已安装，直接使用\n")
  use_alternative <- FALSE
} else {
  cat("TROM包未安装或强制重新安装，尝试安装...\n")
  # 尝试安装TROM
  trom_installed <- FALSE
  tryCatch({
    if (!require("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org/")
    }
    remotes::install_github("Vivianstats/TROM", force = TRUE, upgrade = "never", quiet = TRUE)
    if (require("TROM", character.only = TRUE, quietly = TRUE)) {
      trom_installed <- TRUE
      cat("✓ TROM安装成功！\n")
      use_alternative <- FALSE
    }
  }, error = function(e) {
    cat("TROM安装失败，将使用相关性方法\n")
    use_alternative <- TRUE
  })
  
  if (!trom_installed) {
    use_alternative <- TRUE
  }
}

# ========================== 创建输出目录 ====================================
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# ========================== 读取表达矩阵 ====================================
cat("\n========== 读取表达矩阵 ==========\n")

if (!file.exists(expr_file)) {
  stop("表达矩阵文件不存在: ", expr_file)
}

# 以字符型读取所有数据，避免自动转换
expr_data_raw <- read.csv(expr_file, 
                          check.names = FALSE, 
                          stringsAsFactors = FALSE,
                          colClasses = "character")
cat(sprintf("原始表达矩阵维度: %d 行, %d 列\n", nrow(expr_data_raw), ncol(expr_data_raw)))

# ========================== 处理表达矩阵 ====================================
cat("\n========== 处理表达矩阵 ==========\n")

# 处理第一列（基因ID）
first_col_values <- expr_data_raw[[1]]
cat(sprintf("第一列前几个值: %s\n", paste(head(first_col_values), collapse = ", ")))

# 假设第一列是基因ID，将其设为行名
if (is.character(first_col_values)) {
  # 清理基因ID
  clean_ids <- trimws(first_col_values)
  clean_ids <- gsub("[^[:alnum:]._-]", "_", clean_ids)
  clean_ids <- make.unique(clean_ids)
  rownames(expr_data_raw) <- clean_ids
  expr_data <- expr_data_raw[, -1, drop = FALSE]
  cat("✓ 已将第一列设为行名（基因ID）\n")
} else {
  expr_data <- expr_data_raw
  cat("第一列不是字符型，保持原样\n")
}

# ===== 强制数值转换 =====
cat("\n强制转换所有列为数值型...\n")

convert_to_numeric_value <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  
  if (is.character(x)) {
    if (remove_units) {
      x <- gsub("\\s*(kb|MB|GB|bp|Kb|Mb|Gb|KB|K|M|G|个|条|个基因|个转录本)\\s*", "", x, ignore.case = TRUE)
    }
    
    if (thousand_separator != "") {
      x <- gsub(thousand_separator, "", x)
    }
    
    if (decimal_separator != ".") {
      x <- gsub(decimal_separator, ".", x)
    }
    
    x <- gsub("[^-0-9.eE]", "", x)
    x[x == ""] <- NA
  }
  
  result <- suppressWarnings(as.numeric(x))
  return(result)
}

original_colnames <- colnames(expr_data)
expr_data_numeric <- as.data.frame(matrix(NA, 
                                          nrow = nrow(expr_data), 
                                          ncol = ncol(expr_data)))
colnames(expr_data_numeric) <- original_colnames
rownames(expr_data_numeric) <- rownames(expr_data)

for (i in 1:ncol(expr_data)) {
  col_name <- original_colnames[i]
  raw_values <- expr_data[[i]]
  converted <- convert_to_numeric_value(raw_values)
  expr_data_numeric[[i]] <- converted
}

expr_data <- expr_data_numeric

# 检查是否还有非数值列
numeric_cols <- sapply(expr_data, is.numeric)
cat(sprintf("\n转换后数值型列: %d/%d\n", sum(numeric_cols), ncol(expr_data)))

if (sum(numeric_cols) == 0) {
  debug_file <- file.path(output_dir, "debug_original_data.csv")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(expr_data_raw, file = debug_file, row.names = FALSE)
  stop("无法将数据转换为数值型，原始数据已保存到: ", debug_file)
}

# ===== 处理缺失值 =====
cat("\n处理缺失值...\n")

na_count <- sum(is.na(expr_data))
if (na_count > 0) {
  cat(sprintf("发现 %d 个缺失值\n", na_count))
  
  na_per_row <- rowSums(is.na(expr_data)) / ncol(expr_data)
  rows_to_remove <- which(na_per_row > na_threshold)
  
  if (length(rows_to_remove) > 0) {
    cat(sprintf("移除缺失值过多的基因: %d 个\n", length(rows_to_remove)))
    expr_data <- expr_data[-rows_to_remove, ]
  }
  
  if (any(is.na(expr_data))) {
    row_means <- rowMeans(expr_data, na.rm = TRUE)
    for (i in 1:nrow(expr_data)) {
      na_idx <- is.na(expr_data[i, ])
      if (any(na_idx)) {
        expr_data[i, na_idx] <- row_means[i]
      }
    }
    cat("已用行均值填充剩余缺失值\n")
  }
}

# ========================== 读取分组信息 ============================
cat("\n========== 读取分组信息 ==========\n")

if (!file.exists(group_file)) {
  stop("分组文件不存在: ", group_file)
}

# 读取分组文件
group_info <- read.csv(group_file, 
                       stringsAsFactors = FALSE,
                       check.names = FALSE,
                       na.strings = c("", "NA", "NULL", "null", "NaN", "nan"))

cat("分组文件原始内容（前10行）:\n")
print(head(group_info, 10))

# 检查列名
cat("\n分组文件列名:", paste(colnames(group_info), collapse = ", "), "\n")

# 检查指定的列名是否存在
if (!group_sample_col %in% colnames(group_info)) {
  cat(sprintf("警告：指定的样品名列名 '%s' 不存在\n", group_sample_col))
  # 尝试自动识别
  if (ncol(group_info) >= 2) {
    group_sample_col <- colnames(group_info)[1]
    group_group_col <- colnames(group_info)[2]
    cat(sprintf("自动使用第一列作为样品名: %s，第二列作为组别: %s\n", 
                group_sample_col, group_group_col))
  } else {
    stop("无法识别分组文件的列")
  }
}

# 清理数据
group_info[[group_sample_col]] <- trimws(as.character(group_info[[group_sample_col]]))
group_info[[group_group_col]] <- trimws(as.character(group_info[[group_group_col]]))

# 移除空值行
group_info <- group_info[!is.na(group_info[[group_sample_col]]) & 
                           group_info[[group_sample_col]] != "" &
                           !is.na(group_info[[group_group_col]]) & 
                           group_info[[group_group_col]] != "", ]

cat(sprintf("\n清理后分组信息: %d 个样品\n", nrow(group_info)))

# ========================== 匹配样品名 ====================================
cat("\n========== 样品名匹配 ==========\n")

expr_samples <- colnames(expr_data)
group_samples <- group_info[[group_sample_col]]

cat("表达矩阵样品名（前10个）:", paste(head(expr_samples, 10), collapse = ", "), "\n")
cat("分组文件样品名（前10个）:", paste(head(group_samples, 10), collapse = ", "), "\n")

# 找出共同样品
common_samples <- intersect(expr_samples, group_samples)
cat(sprintf("\n共同样品数量: %d\n", length(common_samples)))

if (length(common_samples) == 0) {
  # 尝试模糊匹配
  cat("\n尝试模糊匹配...\n")
  
  # 标准化样品名
  expr_samples_norm <- toupper(gsub("[[:space:]]", "", expr_samples))
  group_samples_norm <- toupper(gsub("[[:space:]]", "", group_samples))
  
  # 尝试匹配
  matches <- match(group_samples_norm, expr_samples_norm)
  matched_idx <- which(!is.na(matches))
  
  if (length(matched_idx) > 0) {
    cat(sprintf("找到 %d 个模糊匹配\n", length(matched_idx)))
    # 使用标准化后的名称
    colnames(expr_data) <- expr_samples_norm
    group_samples <- group_samples_norm
    common_samples <- intersect(colnames(expr_data), group_samples)
  } else {
    stop("无法匹配任何样品名，请检查分组文件")
  }
}

# 只保留共同样品
expr_data <- expr_data[, common_samples, drop = FALSE]
group_info <- group_info[group_samples %in% common_samples, , drop = FALSE]

# 按表达矩阵的顺序重新排列分组信息
group_info <- group_info[match(colnames(expr_data), group_info[[group_sample_col]]), ]

groups <- group_info[[group_group_col]]
unique_groups <- unique(groups)
cat("检测到组别:", paste(unique_groups, collapse = ", "), "\n")

# 检查每个组的样品数
for (g in unique_groups) {
  n_samples <- sum(groups == g)
  cat(sprintf("  组 %s: %d 个样品\n", g, n_samples))
}

# ========================== 数据标准化 ====================================
cat("\n========== 数据标准化 ==========\n")

if (perform_zscore) {
  cat("对数据进行Z-score标准化...\n")
  expr_data_scaled <- expr_data
  
  if (zscore_method == "row") {
    # 按基因（行）标准化
    for (i in 1:nrow(expr_data)) {
      row_mean <- mean(as.numeric(expr_data[i, ]), na.rm = TRUE)
      row_sd <- sd(as.numeric(expr_data[i, ]), na.rm = TRUE)
      if (row_sd > 0) {
        expr_data_scaled[i, ] <- (expr_data[i, ] - row_mean) / row_sd
      }
    }
    cat("  已完成按基因Z-score标准化\n")
  } else {
    # 按样品（列）标准化
    expr_data_scaled <- scale(expr_data)
    cat("  已完成按样品Z-score标准化\n")
  }
  
  cat(sprintf("标准化后数据范围: [%.3f, %.3f]\n", 
              min(expr_data_scaled, na.rm = TRUE), 
              max(expr_data_scaled, na.rm = TRUE)))
  
  expr_data_for_analysis <- expr_data_scaled
} else {
  expr_data_for_analysis <- expr_data
}

# ========================== 准备TROM输入数据 ================================
cat("\n========== 准备TROM输入数据 ==========\n")

expr_with_id <- cbind(rownames(expr_data_for_analysis), expr_data_for_analysis)
colnames(expr_with_id)[1] <- "GeneID"

# ========================== 自动选择Z-score阈值 ==============================
cat("\n========== 自动选择Z-score阈值 ==========\n")

if (!use_alternative && auto_select_z && is.null(z_threshold)) {
  cat("使用choose.z()函数自动选择最优阈值...\n")
  
  tryCatch({
    # 使用choose.z函数自动选择阈值
    # mode = FALSE 表示使用mode + sd方法（官方推荐）
    z_threshold <- choose.z(expr_with_id, mode = FALSE)
    cat(sprintf("✓ 自动选择的Z-score阈值: %.3f\n", z_threshold))
    
    # 显示阈值选择的详细信息
    cat("\n阈值选择原理：基于mode + sd方法，选择能产生最稀疏但稳定的转录组对应图的阈值\n")
    
  }, error = function(e) {
    cat("自动选择阈值失败:", e$message, "\n")
    cat("使用默认阈值 1.5\n")
    z_threshold <- 1.5
  })
  
} else if (!use_alternative && !is.null(z_threshold)) {
  cat(sprintf("使用手动设置的阈值: %.3f\n", z_threshold))
} else if (!use_alternative) {
  cat("使用默认阈值 1.5\n")
  z_threshold <- 1.5
}

# 如果TROM不可用，提示将使用相关性
if (use_alternative) {
  cat("TROM包不可用，将使用Pearson相关性方法\n")
}

# ========================== 检查Z-score阈值的合理性 ==========================
if (!use_alternative && exists("z_threshold")) {
  cat("\n========== 验证阈值合理性 ==========\n")
  
  # 计算每个样品的关联基因数量
  z_scores <- scale(expr_data_for_analysis)
  genes_per_sample <- colSums(abs(z_scores) > z_threshold, na.rm = TRUE)
  
  cat(sprintf("阈值 %.3f 下:\n", z_threshold))
  cat(sprintf("  每个样品平均关联基因数: %.1f\n", mean(genes_per_sample)))
  cat(sprintf("  关联基因数范围: %d - %d\n", min(genes_per_sample), max(genes_per_sample)))
  
  if (mean(genes_per_sample) < 5) {
    cat("⚠️ 警告：关联基因数太少（<5），结果可能不稳健\n")
    cat("  建议：考虑使用更低的阈值或相关性方法\n")
  } else if (mean(genes_per_sample) > 100) {
    cat("⚠️ 警告：关联基因数太多（>100），可能包含噪声\n")
    cat("  建议：考虑使用更高的阈值\n")
  } else {
    cat("✓ 阈值选择合理\n")
  }
}

# ========================== 相似度分析 ====================================
cat("\n========== 运行相似度分析 ==========\n")

similarity_matrix <- NULL
method_name <- NULL

if (!use_alternative) {
  cat(sprintf("使用TROM方法分析，Z-score阈值: %.3f\n", z_threshold))
  
  tryCatch({
    similarity_matrix <- ws.trom(sp_gene_expr = expr_with_id, z_thre = z_threshold)
    method_name <- "TROM"
    cat("✓ TROM分析成功！\n")
    
  }, error = function(e) {
    cat("TROM分析出错:", e$message, "\n")
    cat("切换到替代方法...\n")
    use_alternative <- TRUE
  })
}

if (use_alternative || is.null(similarity_matrix)) {
  cat("使用Pearson相关性作为替代方法\n")
  
  expr_matrix <- as.matrix(expr_data_for_analysis)
  similarity_matrix <- cor(expr_matrix, method = "pearson", use = "pairwise.complete.obs")
  method_name <- "Correlation"
  cat("✓ 相关性分析成功！\n")
}

# ========================== 相似度矩阵验证和归一化 ============================
cat("\n========== 相似度矩阵验证和归一化 ==========\n")

# 验证函数
validate_similarity_matrix <- function(mat, method_name) {
  cat(sprintf("\n--- %s 相似度矩阵验证 ---\n", method_name))
  
  # 维度
  cat(sprintf("矩阵维度: %d x %d\n", nrow(mat), ncol(mat)))
  
  # 数值范围
  mat_range <- range(mat, na.rm = TRUE)
  cat(sprintf("数值范围: [%.3f, %.3f]\n", mat_range[1], mat_range[2]))
  
  # 对角线
  diag_values <- diag(mat)
  cat(sprintf("对角线范围: [%.3f, %.3f]\n", min(diag_values), max(diag_values)))
  
  # 非对角线值
  non_diag <- mat
  diag(non_diag) <- NA
  non_diag_values <- as.vector(non_diag)
  non_diag_values <- non_diag_values[!is.na(non_diag_values)]
  
  if (length(non_diag_values) > 0) {
    cat(sprintf("非对角线: 均值=%.3f, 中位数=%.3f, 标准差=%.3f\n", 
                mean(non_diag_values), median(non_diag_values), sd(non_diag_values)))
  }
  
  # 合理性判断
  cat("\n合理性判断:\n")
  if (mat_range[2] > 10) {
    cat("  ⚠️ 数值过大（>10），需要归一化\n")
    return(FALSE)
  } else if (mat_range[2] <= 1 && mat_range[1] >= -1) {
    cat("  ✓ 数值范围合理（在[-1,1]区间内）\n")
    return(TRUE)
  } else if (mat_range[2] <= 1 && mat_range[1] >= 0) {
    cat("  ✓ 数值范围合理（在[0,1]区间内）\n")
    return(TRUE)
  } else {
    cat("  ⚠️ 数值范围异常\n")
    return(FALSE)
  }
}

# 归一化函数
normalize_matrix <- function(mat, min_val = NULL, max_val = NULL) {
  if (is.null(min_val)) {
    min_val <- min(mat, na.rm = TRUE)
  }
  if (is.null(max_val)) {
    max_val <- max(mat, na.rm = TRUE)
  }
  
  if (max_val > min_val) {
    mat_norm <- (mat - min_val) / (max_val - min_val)
    cat(sprintf("归一化范围: [%.3f, %.3f] -> [0, 1]\n", min_val, max_val))
    return(mat_norm)
  } else {
    cat("警告：最大值等于最小值，无法归一化\n")
    return(mat)
  }
}

# 验证相似度矩阵
is_reasonable <- validate_similarity_matrix(similarity_matrix, method_name)

# 如果需要归一化
if (!is_reasonable && normalize_similarity) {
  cat("\n对相似度矩阵进行归一化处理...\n")
  similarity_matrix <- normalize_matrix(similarity_matrix, similarity_min, similarity_max)
  validate_similarity_matrix(similarity_matrix, paste0(method_name, "_normalized"))
}

# ========================== 计算组间平均相似度 ============================
cat("\n========== 计算组间平均相似度 ==========\n")

group_labels <- groups
names(group_labels) <- colnames(expr_data)

# 获取有效样品
valid_samples <- intersect(colnames(similarity_matrix), names(group_labels))

if (length(valid_samples) == 0) {
  stop("相似度矩阵样品名与分组信息不匹配！")
}

group_labels_valid <- group_labels[valid_samples]

# 按组别排序
group_order <- order(group_labels_valid)
samples_ordered <- valid_samples[group_order]
group_labels_ordered <- group_labels_valid[samples_ordered]

# 重新排序相似度矩阵
similarity_ordered <- similarity_matrix[samples_ordered, samples_ordered]

# 计算组间平均
n_groups <- length(unique_groups)
group_mean_matrix <- matrix(NA, nrow = n_groups, ncol = n_groups,
                            dimnames = list(unique_groups, unique_groups))

for (i in 1:n_groups) {
  for (j in 1:n_groups) {
    group_i_samples <- names(group_labels_valid)[group_labels_valid == unique_groups[i]]
    group_j_samples <- names(group_labels_valid)[group_labels_valid == unique_groups[j]]
    
    if (length(group_i_samples) > 0 && length(group_j_samples) > 0) {
      sub_matrix <- similarity_matrix[group_i_samples, group_j_samples, drop = FALSE]
      group_mean_matrix[i, j] <- mean(as.matrix(sub_matrix), na.rm = TRUE)
    }
  }
}

cat("\n组间平均相似度矩阵:\n")
print(round(group_mean_matrix, 4))

# 检查组间平均值的合理性
if (any(group_mean_matrix > 1, na.rm = TRUE) || any(group_mean_matrix < -1, na.rm = TRUE)) {
  cat("\n⚠️ 警告：组间平均相似度超出[-1,1]范围！\n")
  
  if (normalize_similarity) {
    cat("对组间平均矩阵进行归一化...\n")
    group_mean_matrix <- normalize_matrix(group_mean_matrix)
    cat("\n归一化后的组间平均相似度矩阵:\n")
    print(round(group_mean_matrix, 4))
  }
}

# ========================== 保存结果 ====================================
cat("\n========== 保存结果 ==========\n")

# 保存原始相似度矩阵
write.csv(similarity_matrix, 
          file = file.path(output_dir, paste0("full_", method_name, "_", trom_score_file)),
          quote = FALSE, row.names = TRUE)

# 保存组间平均相似度矩阵
write.csv(group_mean_matrix, 
          file = file.path(output_dir, paste0("group_mean_", method_name, "_", trom_score_file)),
          quote = FALSE, row.names = TRUE)

# 保存分析参数
params_file <- file.path(output_dir, "analysis_parameters.txt")
sink(params_file)
cat("TROM分析参数\n")
cat("============\n\n")
cat("日期:", date(), "\n")
cat("方法:", method_name, "\n")
if (!use_alternative) {
  cat("Z-score阈值:", z_threshold, "\n")
  cat("阈值选择方式:", ifelse(auto_select_z && is.null(z_threshold), "自动", "手动"), "\n")
}
cat("数据标准化:", ifelse(perform_zscore, paste0("是 (", zscore_method, ")"), "否"), "\n")
cat("归一化:", ifelse(normalize_similarity, "是", "否"), "\n")
cat("\n组别信息:\n")
for (g in unique_groups) {
  cat(sprintf("  %s: %d个样品\n", g, sum(groups == g)))
}
sink()

cat(sprintf("结果已保存到: %s\n", output_dir))

# ========================== 生成热图 ====================================
cat("\n========== 生成热图 ==========\n")

if (nrow(similarity_ordered) > 0) {
  
  # 定义颜色函数
  color_palette <- colorRampPalette(c(heatmap_color_low, heatmap_color_mid, heatmap_color_high))(100)
  
  # ===== 样品级别热图 =====
  cat("\n--- 生成样品级别热图 ---\n")
  
  sample_heatmap_file <- gsub(".png", paste0("_", method_name, ".png"), heatmap_file)
  png(file.path(output_dir, sample_heatmap_file), 
      width = heatmap_width, 
      height = heatmap_height, 
      units = "in", 
      res = heatmap_res)
  
  # 创建注释
  annotation_col <- data.frame(Group = group_labels_ordered)
  rownames(annotation_col) <- colnames(similarity_ordered)
  
  # 定义组别颜色
  if (n_groups <= 8) {
    group_colors <- brewer.pal(n_groups, "Set1")
  } else {
    group_colors <- rainbow(n_groups)
  }
  names(group_colors) <- unique_groups
  ann_colors <- list(Group = group_colors)
  
  # 绘制样品级别热图
  pheatmap(similarity_ordered,
           color = color_palette,
           annotation_col = annotation_col,
           annotation_row = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = heatmap_fontsize_row,
           fontsize_col = heatmap_fontsize_col,
           main = paste(method_name, "Similarity Between Samples"),
           xlab = heatmap_xlab,
           ylab = heatmap_ylab,
           legend_labels = heatmap_key_label,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = NA,
           # 数值显示相关参数
           display_numbers = sample_display_numbers,
           number_color = sample_number_color,
           fontsize_number = sample_number_fontsize)
  
  dev.off()
  cat(sprintf("样品级别热图已保存: %s\n", file.path(output_dir, sample_heatmap_file)))
  cat(sprintf("  数值显示: %s, 颜色: %s, 大小: %d\n", 
              ifelse(sample_display_numbers, "是", "否"),
              sample_number_color, sample_number_fontsize))
  
  # ===== 组间平均热图 =====
  if (all(!is.na(group_mean_matrix))) {
    cat("\n--- 生成组间平均热图 ---\n")
    
    group_heatmap_file <- gsub(".png", paste0("_", method_name, "_group_mean.png"), heatmap_file)
    png(file.path(output_dir, group_heatmap_file), 
        width = 6, 
        height = 5, 
        units = "in", 
        res = heatmap_res)
    
    # 根据数值范围动态设置显示格式
    if (max(group_mean_matrix, na.rm = TRUE) < 0.01) {
      current_number_format <- "%.4f"
    } else {
      current_number_format <- group_number_format
    }
    
    # 如果启用条件颜色，创建颜色矩阵
    if (use_conditional_colors && nrow(group_mean_matrix) <= 10) {
      # 创建颜色矩阵
      number_color_matrix <- matrix(group_number_color, 
                                    nrow = nrow(group_mean_matrix), 
                                    ncol = ncol(group_mean_matrix))
      
      # 根据阈值设置颜色
      for (i in 1:nrow(group_mean_matrix)) {
        for (j in 1:ncol(group_mean_matrix)) {
          if (!is.na(group_mean_matrix[i, j])) {
            val <- group_mean_matrix[i, j]
            if (val >= conditional_thresholds[1]) {
              number_color_matrix[i, j] <- conditional_colors[1]  # 高值颜色
            } else if (val >= conditional_thresholds[2]) {
              number_color_matrix[i, j] <- conditional_colors[2]  # 中值颜色
            } else {
              number_color_matrix[i, j] <- conditional_colors[3]  # 低值颜色
            }
          }
        }
      }
      
      # 绘制组间平均热图（带条件颜色）
      pheatmap(group_mean_matrix,
               color = color_palette,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize_row = 12,
               fontsize_col = 12,
               main = paste("Average", method_name, "Score Between Groups"),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               display_numbers = group_display_numbers,
               number_color = number_color_matrix,
               fontsize_number = group_number_fontsize,
               number_format = current_number_format,
               border_color = "grey60")
      
      cat("  使用条件颜色: 高值(>=", conditional_thresholds[1], ")=", conditional_colors[1], 
          ", 中值(>=", conditional_thresholds[2], ")=", conditional_colors[2], 
          ", 低值=", conditional_colors[3], "\n", sep="")
      
    } else {
      # 绘制组间平均热图（单一颜色）
      pheatmap(group_mean_matrix,
               color = color_palette,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize_row = 12,
               fontsize_col = 12,
               main = paste("Average", method_name, "Score Between Groups"),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               display_numbers = group_display_numbers,
               number_color = group_number_color,
               fontsize_number = group_number_fontsize,
               number_format = current_number_format,
               border_color = "grey60")
    }
    
    dev.off()
    cat(sprintf("组间平均热图已保存: %s\n", file.path(output_dir, group_heatmap_file)))
    cat(sprintf("  数值显示: %s, 颜色: %s, 大小: %d, 格式: %s\n", 
                ifelse(group_display_numbers, "是", "否"),
                ifelse(use_conditional_colors, "条件颜色", group_number_color),
                group_number_fontsize, current_number_format))
    
    # ===== 额外：生成带数值矩阵的CSV文件（方便查看）=====
    numbers_file <- file.path(output_dir, paste0("group_mean_values_", method_name, ".csv"))
    write.csv(round(group_mean_matrix, 4), file = numbers_file, quote = FALSE)
    cat(sprintf("数值矩阵已保存到: %s\n", numbers_file))
  }
}

# ========================== 输出统计摘要 ====================================
cat("\n========== 分析完成 ==========\n")
cat("输出文件保存在:", output_dir, "\n")
cat("使用的方法:", method_name, "\n")
if (!use_alternative) {
  cat("Z-score阈值:", z_threshold, "\n")
}

# 输出组间相似度总结
if (n_groups >= 2 && all(!is.na(group_mean_matrix))) {
  cat("\n组间相似度总结:\n")
  for (i in 1:(n_groups-1)) {
    for (j in (i+1):n_groups) {
      cat(sprintf("  %s vs %s: 平均%s得分 = %.4f\n", 
                  unique_groups[i], unique_groups[j], method_name,
                  group_mean_matrix[i, j]))
    }
  }
}

# 输出组内相似度
if (all(!is.na(diag(group_mean_matrix)))) {
  cat("\n组内相似度:\n")
  for (i in 1:n_groups) {
    cat(sprintf("  %s组内: 平均%s得分 = %.4f\n", 
                unique_groups[i], method_name,
                group_mean_matrix[i, i]))
  }
}

cat("\n=== 脚本运行完成 ===\n")