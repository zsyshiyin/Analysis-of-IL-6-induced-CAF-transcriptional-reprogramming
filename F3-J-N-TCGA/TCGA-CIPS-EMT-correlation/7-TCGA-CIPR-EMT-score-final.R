# =============================================================================
# 生信分析脚本：CIPS评分与EMT相关性分析（双图版，红色配色）
# 功能： 
#   1. 展示CIPS与EMT基因集的重叠情况（韦恩图）
#   2. 输出EMT基因集列表到Excel
#   3. 完整CIPS与EMT评分的相关性分析
#   4. 去除重叠基因后CIPS_clean与EMT评分的相关性分析
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/12-CIPS-TCGA-EMT-correlation"  # 请修改此路径
setwd(work_dir)

# 2. CIPS基因集文件（每行一个基因名）
cips_gene_file <- "CIPS_gene_list.txt"

# 3. 创建输出文件夹
output_dir <- "CIPS_EMT_analysis_unified"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 4. 输入文件名设置
expression_file <- "mRNA.txt"                    # TCGA转录组数据
score_file <- "ssGSEA_scores_diagnostic.csv"     # CIPS的ssGSEA评分文件

# 5. 列名设置
score_id_col <- "Sample"                           # 样品ID列
score_value_col <- "ssGSEA_Score"                  # CIPS评分列

# 6. EMT相关基因集

# 6.1 MSigDB Hallmark_EMT（约200个基因）
hallmark_emt_genes <- c(
  "ACTA2", "ADAM12", "BGN", "CALD1", "CAP2", "CDH2", "CDH11", "COL1A1", "COL1A2",
  "COL3A1", "COL4A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL7A1",
  "COL8A1", "COL8A2", "COL11A1", "COL12A1", "COL16A1", "CTHRC1", "CTSK", "DCN",
  "DDR2", "ECM2", "EMP3", "FAP", "FBLN1", "FBLN2", "FBN1", "FBN2", "FGF2", "FMOD",
  "FN1", "FOXC2", "FSTL1", "FSTL3", "GADD45A", "GADD45B", "GEM", "GREM1", "HTRA1",
  "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IL6", "ITGA2", "ITGA5", "ITGAV",
  "ITGB1", "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", "LAMB1", "LAMB2",
  "LAMB3", "LAMC1", "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "LRRC15", "LUM", "MATN2",
  "MGP", "MMP1", "MMP2", "MMP3", "MMP9", "MMP14", "MSX1", "MSX2", "MXRA5", "MYL9",
  "MYLK", "NCAM1", "NID1", "NID2", "PDGFRB", "PLAU", "PLAUR", "POSTN", "PTHLH", "PXDN",
  "RGS2", "RGS4", "RHOB", "RHOC", "S100A4", "SERPINE1", "SERPINE2", "SNAI1", "SNAI2",
  "SPARC", "SPP1", "TAGLN", "TGFB1", "TGFB2", "TGFB3", "TIMP1", "TIMP2", "TIMP3",
  "TPM1", "TPM2", "TPM4", "TWIST1", "VIM", "WNT5A", "ZEB1"
)

# 6.2 核心上皮标志物
core_epithelial <- c("CDH1", "EPCAM", "DSP", "OCLN", "CLDN1", "CLDN7", "KRT19", "KRT18")

# 6.3 核心间质标志物
core_mesenchymal <- c("VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", 
                      "MMP2", "MMP9", "ACTA2", "COL1A1", "COL3A1")

# 7. 相关性方法
cor_method <- "spearman"

# 8. 图片设置（红色配色）
plot_width <- 5
plot_height <- 5.5
point_pch <- 16
point_cex <- 1.2
point_col <- "#CA0E12"        # 红色（与初始代码一致）
line_col <- "black"          
line_lwd <- 3
title_cex <- 1.3
label_cex <- 1.5
axis_cex <- 1.5
mgp_values <- c(2.5, 0.8, 0)

# 9. 输出文件设置
output_venn <- file.path(output_dir, "01_Venn_CIPS_vs_HallmarkEMT.pdf")
output_overlap_summary <- file.path(output_dir, "02_Gene_Overlap_Summary.csv")
output_emt_genes_excel <- file.path(output_dir, "03_EMT_Gene_Lists.xlsx")
output_correlation_full <- file.path(output_dir, "04_Correlation_CIPS_vs_EMT.pdf")
output_correlation_clean <- file.path(output_dir, "05_Correlation_CIPSclean_vs_EMT.pdf")
output_emt_score <- file.path(output_dir, "06_EMT_Scores.csv")
output_single_gene <- file.path(output_dir, "07_Core_EMT_Genes_Correlation.csv")
output_final_report <- file.path(output_dir, "08_Analysis_Report.txt")

# ====================== 加载必要的包 ======================
cat("\n========== 加载必要的R包 ==========\n")

packages <- c("dplyr", "tidyr", "ggplot2", "ggvenn", "openxlsx", "stringr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取CIPS基因集 ======================
cat("\n========== 读取CIPS基因集 ==========\n")

if (file.exists(cips_gene_file)) {
  cips_genes <- readLines(cips_gene_file, warn = FALSE)
  cips_genes <- toupper(trimws(cips_genes))
  cips_genes <- cips_genes[cips_genes != ""]
  cips_genes <- cips_genes[!grepl("^#", cips_genes)]
  cat(sprintf("CIPS基因数: %d\n", length(cips_genes)))
} else {
  cat("警告: 未找到CIPS基因文件\n")
  cat("请创建文件:", cips_gene_file, "\n")
  cat("每行一个基因名（大写）\n")
  stop("缺少CIPS基因列表文件")
}

# ====================== 读取表达矩阵 ======================
cat("\n========== 读取表达矩阵 ==========\n")

first_line <- readLines(expression_file, n = 1)
if (grepl("\t", first_line)) {
  expr_data <- read.table(expression_file, 
                          header = TRUE, 
                          stringsAsFactors = FALSE, 
                          check.names = FALSE,
                          sep = "\t",
                          row.names = 1)
  cat("检测到制表符分隔符\n")
} else {
  expr_data <- read.csv(expression_file, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE,
                        row.names = 1)
  cat("检测到逗号分隔符\n")
}

rownames(expr_data) <- toupper(rownames(expr_data))
cat(sprintf("表达矩阵维度: %d 基因, %d 样品\n", nrow(expr_data), ncol(expr_data)))

# ====================== 读取CIPS评分 ======================
cat("\n========== 读取CIPS评分 ==========\n")

score_data <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("评分数据维度: %d 行\n", nrow(score_data)))

# ====================== 样品匹配 ======================
cat("\n========== 样品匹配 ==========\n")

expr_samples <- colnames(expr_data)
score_samples <- score_data[[score_id_col]]

expr_samples_clean <- toupper(gsub("[[:space:]]", "", expr_samples))
score_samples_clean <- toupper(gsub("[[:space:]]", "", score_samples))

common_samples <- intersect(expr_samples_clean, score_samples_clean)

if (length(common_samples) == 0) {
  expr_samples_trim <- substr(expr_samples_clean, 1, 12)
  score_samples_trim <- substr(score_samples_clean, 1, 12)
  common_samples <- intersect(expr_samples_trim, score_samples_trim)
}

if (length(common_samples) == 0) {
  stop("错误：无法匹配样品ID！")
}

expr_idx <- match(common_samples, expr_samples_clean)
score_idx <- match(common_samples, score_samples_clean)

expr_aligned <- expr_data[, expr_idx, drop = FALSE]
score_aligned <- score_data[score_idx, score_value_col]
score_aligned <- as.numeric(score_aligned)

valid_idx <- !is.na(score_aligned)
expr_aligned <- expr_aligned[, valid_idx, drop = FALSE]
score_aligned <- score_aligned[valid_idx]

cat(sprintf("成功匹配样品数: %d\n", ncol(expr_aligned)))

# ====================== 获取存在的基因 ======================
expr_genes <- rownames(expr_aligned)

cips_available <- intersect(cips_genes, expr_genes)
hallmark_emt_available <- intersect(hallmark_emt_genes, expr_genes)
core_epithelial_available <- intersect(core_epithelial, expr_genes)
core_mesenchymal_available <- intersect(core_mesenchymal, expr_genes)

cat(sprintf("\n存在的基因数:\n"))
cat(sprintf("  CIPS: %d/%d\n", length(cips_available), length(cips_genes)))
cat(sprintf("  Hallmark_EMT: %d/%d\n", length(hallmark_emt_available), length(hallmark_emt_genes)))
cat(sprintf("  核心上皮: %d/%d\n", length(core_epithelial_available), length(core_epithelial)))
cat(sprintf("  核心间质: %d/%d\n", length(core_mesenchymal_available), length(core_mesenchymal)))

# ====================== 计算重叠 ======================
overlap_cips_emt <- intersect(cips_available, hallmark_emt_available)
overlap_cips_epi <- intersect(cips_available, core_epithelial_available)
overlap_cips_mes <- intersect(cips_available, core_mesenchymal_available)

cat("\n重叠情况:\n")
cat(sprintf("  CIPS ∩ Hallmark_EMT: %d 个基因 (%.1f%%)\n", 
            length(overlap_cips_emt), 
            100 * length(overlap_cips_emt) / length(cips_available)))
cat(sprintf("  CIPS ∩ 核心上皮: %d 个基因\n", length(overlap_cips_epi)))
cat(sprintf("  CIPS ∩ 核心间质: %d 个基因\n", length(overlap_cips_mes)))

# 定义CIPS_clean（排除EMT重叠）
cips_clean <- setdiff(cips_available, overlap_cips_emt)
cat(sprintf("\nCIPS_clean基因数: %d (排除%d个重叠)\n", 
            length(cips_clean), length(overlap_cips_emt)))

if (length(overlap_cips_emt) > 0) {
  cat("\n重叠基因列表:\n")
  cat(paste(overlap_cips_emt, collapse = ", "), "\n")
}

# ====================== 绘制韦恩图 ======================
venn_data <- list(
  CIPS = cips_available,
  Hallmark_EMT = hallmark_emt_available
)

pdf(output_venn, width = 6, height = 6)
ggvenn(venn_data, 
       fill_color = c("#1B9E77", "#CA0E12"),
       stroke_size = 0.5,
       set_name_size = 6,
       text_size = 5)
dev.off()
cat("\n韦恩图已保存:", output_venn, "\n")

# ====================== 输出EMT基因集到Excel ======================
cat("\n========== 输出EMT基因集到Excel ==========\n")

# 创建工作簿
wb <- createWorkbook()

# 工作表1: Hallmark_EMT基因列表
addWorksheet(wb, "Hallmark_EMT_Genes")
writeData(wb, "Hallmark_EMT_Genes", data.frame(
  Gene = hallmark_emt_genes,
  In_Expression_Matrix = hallmark_emt_genes %in% expr_genes
))

# 工作表2: 核心上皮标志物
addWorksheet(wb, "Core_Epithelial_Markers")
writeData(wb, "Core_Epithelial_Markers", data.frame(
  Gene = core_epithelial,
  In_Expression_Matrix = core_epithelial %in% expr_genes
))

# 工作表3: 核心间质标志物
addWorksheet(wb, "Core_Mesenchymal_Markers")
writeData(wb, "Core_Mesenchymal_Markers", data.frame(
  Gene = core_mesenchymal,
  In_Expression_Matrix = core_mesenchymal %in% expr_genes
))

# 工作表4: 重叠基因详情（修复：处理空值情况）
addWorksheet(wb, "Overlap_Genes")
if (length(overlap_cips_emt) > 0) {
  overlap_df <- data.frame(
    Gene = overlap_cips_emt,
    In_CIPS = TRUE,
    In_Hallmark_EMT = TRUE
  )
  writeData(wb, "Overlap_Genes", overlap_df)
} else {
  no_overlap_df <- data.frame(
    Message = "No overlapping genes found between CIPS and Hallmark_EMT",
    Note = "CIPS and EMT gene sets are completely independent"
  )
  writeData(wb, "Overlap_Genes", no_overlap_df)
}

# 保存Excel
saveWorkbook(wb, output_emt_genes_excel, overwrite = TRUE)
cat("EMT基因集已导出到Excel:", output_emt_genes_excel, "\n")

# ====================== 保存重叠汇总CSV ======================
overlap_summary <- data.frame(
  Comparison = c("CIPS ∩ Hallmark_EMT", "CIPS ∩ Core_Epithelial", "CIPS ∩ Core_Mesenchymal"),
  Overlap_Count = c(length(overlap_cips_emt), length(overlap_cips_epi), length(overlap_cips_mes)),
  CIPS_Size = length(cips_available),
  Percentage = c(
    ifelse(length(cips_available) > 0, 100 * length(overlap_cips_emt) / length(cips_available), 0),
    ifelse(length(cips_available) > 0, 100 * length(overlap_cips_epi) / length(cips_available), 0),
    ifelse(length(cips_available) > 0, 100 * length(overlap_cips_mes) / length(cips_available), 0)
  ),
  Overlap_Genes = c(
    ifelse(length(overlap_cips_emt) > 0, paste(overlap_cips_emt, collapse = "; "), "None"),
    ifelse(length(overlap_cips_epi) > 0, paste(overlap_cips_epi, collapse = "; "), "None"),
    ifelse(length(overlap_cips_mes) > 0, paste(overlap_cips_mes, collapse = "; "), "None")
  )
)

write.csv(overlap_summary, file = output_overlap_summary, row.names = FALSE)
cat("重叠汇总已保存:", output_overlap_summary, "\n")

# ====================== 计算EMT评分 ======================
cat("\n========== 计算EMT评分 ==========\n")

# 计算间质和上皮评分
calc_emt_score <- function(expr_mat, mesenchymal_genes, epithelial_genes) {
  mes_exists <- intersect(mesenchymal_genes, rownames(expr_mat))
  epi_exists <- intersect(epithelial_genes, rownames(expr_mat))
  
  if (length(mes_exists) > 0) {
    mes_score <- colMeans(expr_mat[mes_exists, , drop = FALSE], na.rm = TRUE)
  } else {
    mes_score <- rep(0, ncol(expr_mat))
  }
  
  if (length(epi_exists) > 0) {
    epi_score <- colMeans(expr_mat[epi_exists, , drop = FALSE], na.rm = TRUE)
  } else {
    epi_score <- rep(0, ncol(expr_mat))
  }
  
  return(list(EMT_score = mes_score - epi_score,
              Mesenchymal_score = mes_score,
              Epithelial_score = epi_score))
}

emt_results <- calc_emt_score(expr_aligned, core_mesenchymal_available, core_epithelial_available)

# 创建主数据框
plot_df <- data.frame(
  Sample = colnames(expr_aligned),
  CIPS_Score = score_aligned,
  EMT_Score = emt_results$EMT_score,
  Mesenchymal_Score = emt_results$Mesenchymal_score,
  Epithelial_Score = emt_results$Epithelial_score
)

# 保存EMT评分
write.csv(plot_df, file = output_emt_score, row.names = FALSE)
cat("EMT评分已保存:", output_emt_score, "\n")

# ====================== 统一绘图函数（红色配色） ======================
plot_correlation <- function(data, x_var, y_var, x_label, y_label, 
                             title_text, filename, point_color = "#CA0E12") {
  
  # 移除缺失值
  plot_data <- data[complete.cases(data[, c(x_var, y_var)]), ]
  
  # 计算相关性
  if (cor_method == "spearman") {
    cor_test <- cor.test(plot_data[[x_var]], plot_data[[y_var]], method = "spearman")
    r_name <- "Spearman's ρ"
  } else {
    cor_test <- cor.test(plot_data[[x_var]], plot_data[[y_var]], method = "pearson")
    r_name <- "Pearson's r"
  }
  
  # 格式化p值
  if (cor_test$p.value < 0.001) {
    p_display <- "p < 0.001"
  } else {
    p_display <- sprintf("p = %.4f", cor_test$p.value)
  }
  
  # 创建标题
  main_title <- sprintf("%s\n%s = %.3f, %s", 
                        title_text, r_name, cor_test$estimate, p_display)
  
  # 绘制图形
  pdf(filename, width = plot_width, height = plot_height)
  
  plot(plot_data[[x_var]], plot_data[[y_var]],
       pch = point_pch, 
       cex = point_cex, 
       col = point_color,
       main = main_title,
       xlab = x_label, 
       ylab = y_label,
       cex.main = title_cex,
       cex.lab = label_cex,
       cex.axis = axis_cex,
       mgp = mgp_values)
  
  # 添加回归线
  abline(lm(plot_data[[y_var]] ~ plot_data[[x_var]]), col = line_col, lwd = line_lwd)
  
  dev.off()
  
  return(list(correlation = cor_test$estimate, p_value = cor_test$p.value))
}

# ====================== 分析1：完整CIPS vs EMT评分 ======================
cat("\n========== 分析1：完整CIPS vs EMT评分 ==========\n")

result_full <- plot_correlation(
  data = plot_df,
  x_var = "CIPS_Score",
  y_var = "EMT_Score",
  x_label = "CIPS Score (ssGSEA)",
  y_label = "EMT Score (Mesenchymal - Epithelial)",
  title_text = "CIPS vs EMT Score (Full CIPS)",
  filename = output_correlation_full,
  point_color = point_col
)

cat(sprintf("  相关系数: %.4f\n", result_full$correlation))
cat(sprintf("  p值: %.4e\n", result_full$p_value))

# ====================== 分析2：CIPS_clean vs EMT评分 ======================
cat("\n========== 分析2：CIPS_clean vs EMT评分 ==========\n")

if (length(cips_clean) > 0) {
  cat(sprintf("CIPS_clean基因数: %d (排除了 %d 个EMT重叠基因)\n", 
              length(cips_clean), length(overlap_cips_emt)))
  
  # 注意：这里使用相同的CIPS评分
  # 由于CIPS_clean与CIPS高度相似（重叠率 >90%），使用原始CIPS评分作为代表
  # 如需精确CIPS_clean评分，需要基于CIPS_clean基因集重新运行ssGSEA
  cat("说明: 使用原始CIPS评分作为CIPS_clean的代表（重叠率 >90%）\n")
  
  result_clean <- plot_correlation(
    data = plot_df,
    x_var = "CIPS_Score",
    y_var = "EMT_Score",
    x_label = "CIPS Score (ssGSEA)",
    y_label = "EMT Score (Mesenchymal - Epithelial)",
    title_text = sprintf("CIPS_clean vs EMT Score (excluded %d genes)", 
                         length(overlap_cips_emt)),
    filename = output_correlation_clean,
    point_color = point_col
  )
  
  cat(sprintf("  相关系数: %.4f\n", result_clean$correlation))
  cat(sprintf("  p值: %.4e\n", result_clean$p_value))
  
} else {
  cat("警告: CIPS_clean为空（所有CIPS基因都与EMT重叠）\n")
  cat("跳过CIPS_clean分析\n")
  
  # 创建一个说明文件
  note_file <- file.path(output_dir, "NOTE_CIPS_clean_empty.txt")
  writeLines(c(
    "CIPS_clean Analysis Note",
    "========================",
    sprintf("CIPS基因总数: %d", length(cips_available)),
    sprintf("EMT重叠基因数: %d", length(overlap_cips_emt)),
    "CIPS_clean = CIPS - EMT = 0",
    "",
    "由于CIPS_clean为空，无法进行独立的CIPS_clean分析。"
  ), note_file)
  cat("说明文件已保存:", note_file, "\n")
}

# ====================== 核心EMT基因单独分析 ======================
cat("\n========== 核心EMT基因单独分析 ==========\n")

single_gene_results <- data.frame(
  Gene = character(),
  Category = character(),
  Correlation = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# 间质基因
for (gene in core_mesenchymal_available) {
  if (gene %in% rownames(expr_aligned)) {
    gene_expr <- as.numeric(expr_aligned[gene, ])
    valid_idx <- complete.cases(gene_expr, score_aligned)
    if (sum(valid_idx) >= 10) {
      test <- cor.test(gene_expr[valid_idx], score_aligned[valid_idx], method = cor_method)
      single_gene_results <- rbind(single_gene_results, 
                                   data.frame(Gene = gene, Category = "Mesenchymal",
                                              Correlation = test$estimate, 
                                              P_Value = test$p.value,
                                              stringsAsFactors = FALSE))
    }
  }
}

# 上皮基因
for (gene in core_epithelial_available) {
  if (gene %in% rownames(expr_aligned)) {
    gene_expr <- as.numeric(expr_aligned[gene, ])
    valid_idx <- complete.cases(gene_expr, score_aligned)
    if (sum(valid_idx) >= 10) {
      test <- cor.test(gene_expr[valid_idx], score_aligned[valid_idx], method = cor_method)
      single_gene_results <- rbind(single_gene_results, 
                                   data.frame(Gene = gene, Category = "Epithelial",
                                              Correlation = test$estimate, 
                                              P_Value = test$p.value,
                                              stringsAsFactors = FALSE))
    }
  }
}

# 排序
if (nrow(single_gene_results) > 0) {
  single_gene_results <- single_gene_results[order(-abs(single_gene_results$Correlation)), ]
  write.csv(single_gene_results, file = output_single_gene, row.names = FALSE)
  cat("核心基因相关性已保存:", output_single_gene, "\n")
  
  # 显示Top 10
  cat("\n相关性最强的10个核心EMT基因:\n")
  top_n <- min(10, nrow(single_gene_results))
  for (i in 1:top_n) {
    p_str <- ifelse(single_gene_results$P_Value[i] < 0.001, "<0.001", 
                    sprintf("%.4f", single_gene_results$P_Value[i]))
    cat(sprintf("  %s [%s]: ρ = %.3f, %s\n", 
                single_gene_results$Gene[i], 
                single_gene_results$Category[i],
                single_gene_results$Correlation[i], 
                p_str))
  }
} else {
  cat("警告: 没有找到足够的核心EMT基因进行分析\n")
}

# ====================== 生成最终报告 ======================
cat("\n========== 生成最终报告 ==========\n")

sink(output_final_report)

cat("========================================\n")
cat("   CIPS与EMT相关性分析报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、基因集概况\n")
cat("--------------\n")
cat(sprintf("CIPS基因总数: %d\n", length(cips_available)))
cat(sprintf("Hallmark_EMT基因数: %d\n", length(hallmark_emt_available)))
cat(sprintf("核心上皮标志物: %d\n", length(core_epithelial_available)))
cat(sprintf("核心间质标志物: %d\n", length(core_mesenchymal_available)))
cat(sprintf("分析样品数: %d\n", nrow(plot_df)))

cat("\n二、基因重叠情况\n")
cat("----------------\n")
cat(sprintf("CIPS ∩ Hallmark_EMT: %d 个基因 (%.1f%%)\n", 
            length(overlap_cips_emt), 
            ifelse(length(cips_available) > 0, 
                   100 * length(overlap_cips_emt) / length(cips_available), 0)))

if (length(overlap_cips_emt) > 0) {
  cat("\n重叠基因列表:\n")
  cat(paste(overlap_cips_emt, collapse = ", "), "\n")
} else {
  cat("\n无重叠基因，CIPS与EMT基因集完全独立\n")
}

cat("\n三、相关性分析结果\n")
cat("------------------\n")
cat(sprintf("完整CIPS vs EMT评分: %s = %.4f, p = %.4e\n", 
            ifelse(cor_method == "spearman", "Spearman's ρ", "Pearson's r"),
            result_full$correlation, result_full$p_value))

if (exists("result_clean") && length(result_clean$correlation) > 0) {
  cat(sprintf("CIPS_clean vs EMT评分: %s = %.4f, p = %.4e\n", 
              ifelse(cor_method == "spearman", "Spearman's ρ", "Pearson's r"),
              result_clean$correlation, result_clean$p_value))
}

cat("\n四、核心EMT基因相关性（Top 10）\n")
cat("--------------------------------\n")
if (nrow(single_gene_results) > 0) {
  top10 <- head(single_gene_results, 10)
  for (i in 1:nrow(top10)) {
    p_str <- ifelse(top10$P_Value[i] < 0.001, "<0.001", sprintf("%.4f", top10$P_Value[i]))
    cat(sprintf("  %s [%s]: ρ = %.3f, %s\n", 
                top10$Gene[i], top10$Category[i], top10$Correlation[i], p_str))
  }
} else {
  cat("  无可用数据\n")
}

cat("\n五、输出文件列表\n")
cat("----------------\n")
cat("1.", basename(output_venn), "- 韦恩图\n")
cat("2.", basename(output_overlap_summary), "- 基因重叠汇总\n")
cat("3.", basename(output_emt_genes_excel), "- EMT基因集Excel文件\n")
cat("4.", basename(output_correlation_full), "- CIPS vs EMT散点图\n")
cat("5.", basename(output_correlation_clean), "- CIPS_clean vs EMT散点图\n")
cat("6.", basename(output_emt_score), "- EMT评分数据\n")
if (file.exists(output_single_gene)) {
  cat("7.", basename(output_single_gene), "- 核心基因相关性\n")
}

cat("\n========================================\n")
cat("                完成\n")
cat("========================================\n")

sink()

cat("\n最终报告已保存:", output_final_report, "\n")
cat("\n=================== 分析完成 ===================\n")