# ==================== GMT 靶向通路富集分析（独立运行版） ====================
# 用途：使用自定义 GMT 文件，分析分泌组与转录组在特定通路上的功能关系
# 特点：针对分泌蛋白数量少的情况，放宽筛选阈值，聚焦靶向通路

rm(list = ls())
gc()

# -------------------- 1. 环境准备与路径设置 --------------------
work_dir <- "D:/zsy/SX/Fomal-final/16.7-CIPS-Sec-gmt"  # 请修改为实际路径
setwd(work_dir)

output_dir <- file.path(work_dir, "GMT_Targeted_Pathway_Analysis-Xin")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 加载所需包
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  library(forcats)
  library(ggrepel)
  library(igraph)
  library(ggraph)
})

# 解决函数冲突
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# 安全打印函数（避免 na.print 错误）
safe_cat_df <- function(df, title = NULL, max_rows = 10) {
  if (!is.null(title)) cat("\n", title, "\n", sep = "")
  if (nrow(df) == 0) {
    cat("  (无数据)\n")
    return()
  }
  n_show <- min(max_rows, nrow(df))
  for (i in 1:n_show) {
    row_str <- paste(names(df), "=", sapply(df[i, ], as.character), collapse = ", ")
    cat("  ", i, ". ", row_str, "\n", sep = "")
  }
  if (nrow(df) > max_rows) {
    cat("  ... (共", nrow(df), "行，仅显示前", max_rows, "行)\n")
  }
}

# -------------------- 2. 自定义可视化参数 --------------------
plot_config <- list(
  dotplot_width = 6,
  dotplot_height = 7,
  barplot_width = 10,
  barplot_height = 7,
  dpi = 300,
  base_size = 16,
  title_size = 12,
  color_trans = "#E64B35",
  color_sec = "#4DBBD5",
  color_shared = "#7E6148"
)

# -------------------- 3. 输入文件读取（自动适配列名） --------------------
cat("\n==================== 读取输入文件 ====================\n")

# 转录组差异基因
deg_trans <- read_csv("transcriptome_deg.csv", show_col_types = FALSE)
cat("转录组文件列名:", paste(colnames(deg_trans), collapse = ", "), "\n")

if (!"GeneSymbol" %in% colnames(deg_trans)) {
  stop("转录组文件必须包含 'GeneSymbol' 列")
}

trans_genes <- unique(deg_trans$GeneSymbol)
cat("转录组差异基因数:", length(trans_genes), "\n")

# 自动检测 log2FC 和 P 值列
fc_col <- intersect(colnames(deg_trans), c("log2FC", "logFC", "Log2FC", "LogFC", "fold_change", "FoldChange"))
pval_col <- intersect(colnames(deg_trans), c("P_value", "P.Value", "pvalue", "p_value", "Pval", "adj.P.Val", "FDR", "padj"))

if (length(fc_col) > 0) {
  fc_col <- fc_col[1]
  cat("使用列 '", fc_col, "' 作为倍数变化列\n", sep = "")
} else {
  warning("未找到 log2FC 相关列，将创建虚拟列")
  deg_trans$log2FC <- 0
  fc_col <- "log2FC"
}

if (length(pval_col) > 0) {
  pval_col <- pval_col[1]
  cat("使用列 '", pval_col, "' 作为 P 值列\n", sep = "")
} else {
  warning("未找到 P 值相关列，将创建虚拟列")
  deg_trans$P_value <- 1
  pval_col <- "P_value"
}

# 提取关键信息
trans_fc <- deg_trans %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  select(GeneSymbol, log2FC = !!sym(fc_col), P_value = !!sym(pval_col))

# 差异分泌蛋白
sec_proteins <- read_csv("secretome_proteins.csv", show_col_types = FALSE)
cat("分泌蛋白文件列名:", paste(colnames(sec_proteins), collapse = ", "), "\n")

if (!"GeneSymbol" %in% colnames(sec_proteins)) {
  stop("分泌蛋白文件必须包含 'GeneSymbol' 列")
}
sec_genes <- unique(sec_proteins$GeneSymbol)
cat("差异分泌蛋白数:", length(sec_genes), "\n")

# 背景基因
bg_genes_df <- read_csv("background_genes.csv", show_col_types = FALSE)
if (!"GeneSymbol" %in% colnames(bg_genes_df)) {
  stop("背景基因文件必须包含 'GeneSymbol' 列")
}
bg_genes <- unique(bg_genes_df$GeneSymbol)
cat("背景基因总数:", length(bg_genes), "\n")

# -------------------- 4. 自动查找所有 GMT 文件 --------------------
cat("\n==================== 查找 GMT 文件 ====================\n")

gmt_files <- list.files(
  path = work_dir,
  pattern = "\\.gmt$",
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = TRUE
)

if (length(gmt_files) == 0) {
  stop("工作目录及其子文件夹中未找到任何 .gmt 文件")
}

cat("找到", length(gmt_files), "个 GMT 文件:\n")
for (i in seq_along(gmt_files)) {
  cat("  ", i, ": ", basename(gmt_files[i]), "\n", sep = "")
}

# -------------------- 5. GMT 文件读取函数 --------------------
read_gmt <- function(gmt_path) {
  if (!file.exists(gmt_path)) {
    warning("GMT 文件不存在: ", gmt_path)
    return(NULL)
  }
  
  gmt_lines <- readLines(gmt_path, warn = FALSE)
  gmt_list <- list()
  
  for (line in gmt_lines) {
    if (nchar(line) == 0) next
    
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 3) next
    
    pathway_name <- fields[1]
    genes <- fields[-(1:2)]
    genes <- genes[genes != "" & !is.na(genes)]
    
    if (length(genes) > 0) {
      gmt_list[[pathway_name]] <- genes
    }
  }
  
  return(gmt_list)
}

# -------------------- 6. 超几何检验函数 --------------------
hypergeometric_test <- function(gene_set, pathway_genes, background_size) {
  overlap <- intersect(gene_set, pathway_genes)
  k <- length(overlap)
  
  if (k == 0) {
    return(data.frame(
      Overlap_Count = 0,
      Overlap_Genes = "",
      P_value = 1,
      Enrichment_Fold = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  m <- length(pathway_genes)
  n <- background_size - m
  q <- length(gene_set)
  
  p_value <- phyper(k - 1, m, n, q, lower.tail = FALSE)
  expected <- q * (m / background_size)
  enrichment_fold <- k / expected
  
  return(data.frame(
    Overlap_Count = k,
    Overlap_Genes = paste(overlap, collapse = ";"),
    P_value = p_value,
    Enrichment_Fold = enrichment_fold,
    stringsAsFactors = FALSE
  ))
}

# -------------------- 7. 批量处理所有 GMT 文件 --------------------
cat("\n==================== 批量处理 GMT 文件 ====================\n")

all_enrichment_results <- list()
file_summary <- data.frame(
  File = character(),
  Total_Pathways = integer(),
  Pathways_With_Trans_Overlap = integer(),
  Pathways_With_Sec_Overlap = integer(),
  Pathways_With_Both_Overlap = integer(),
  stringsAsFactors = FALSE
)

for (f in seq_along(gmt_files)) {
  gmt_file <- gmt_files[f]
  file_name <- basename(gmt_file)
  file_base <- tools::file_path_sans_ext(file_name)
  
  cat("\n处理文件 ", f, "/", length(gmt_files), ": ", file_name, "\n", sep = "")
  
  custom_pathways <- read_gmt(gmt_file)
  
  if (is.null(custom_pathways) || length(custom_pathways) == 0) {
    cat("  警告：文件为空或格式错误，跳过\n")
    next
  }
  
  cat("  通路数:", length(custom_pathways), "\n")
  
  # 靶向富集分析
  targeted_results <- list()
  
  for (i in seq_along(custom_pathways)) {
    pathway_name <- names(custom_pathways)[i]
    pathway_genes <- custom_pathways[[pathway_name]]
    
    trans_test <- hypergeometric_test(trans_genes, pathway_genes, length(bg_genes))
    sec_test <- hypergeometric_test(sec_genes, pathway_genes, length(bg_genes))
    
    targeted_results[[i]] <- data.frame(
      File = file_name,
      Pathway = pathway_name,
      Total_Genes_in_Pathway = length(pathway_genes),
      Trans_Overlap_Count = trans_test$Overlap_Count,
      Trans_Overlap_Genes = trans_test$Overlap_Genes,
      Trans_P_value = trans_test$P_value,
      Trans_Enrichment_Fold = trans_test$Enrichment_Fold,
      Sec_Overlap_Count = sec_test$Overlap_Count,
      Sec_Overlap_Genes = sec_test$Overlap_Genes,
      Sec_P_value = sec_test$P_value,
      Sec_Enrichment_Fold = sec_test$Enrichment_Fold,
      stringsAsFactors = FALSE
    )
  }
  
  file_results <- bind_rows(targeted_results) %>%
    mutate(
      Has_Trans_Overlap = Trans_Overlap_Count > 0,
      Has_Sec_Overlap = Sec_Overlap_Count > 0,
      Both_Have_Overlap = Has_Trans_Overlap & Has_Sec_Overlap,
      Trans_Significant = Trans_P_value < 0.1,
      Sec_Significant = Sec_P_value < 0.1
    ) %>%
    arrange(desc(Both_Have_Overlap), Trans_P_value, Sec_P_value)
  
  all_enrichment_results[[f]] <- file_results
  
  file_summary <- file_summary %>%
    add_row(
      File = file_name,
      Total_Pathways = nrow(file_results),
      Pathways_With_Trans_Overlap = sum(file_results$Has_Trans_Overlap),
      Pathways_With_Sec_Overlap = sum(file_results$Has_Sec_Overlap),
      Pathways_With_Both_Overlap = sum(file_results$Both_Have_Overlap)
    )
  
  cat("  转录组有交集:", sum(file_results$Has_Trans_Overlap), 
      ", 分泌组有交集:", sum(file_results$Has_Sec_Overlap),
      ", 两层均有:", sum(file_results$Both_Have_Overlap), "\n")
  
  # 为每个 GMT 文件创建子文件夹
  file_output_dir <- file.path(output_dir, file_base)
  if (!dir.exists(file_output_dir)) dir.create(file_output_dir, recursive = TRUE)
  
  write_csv(file_results, file.path(file_output_dir, 
                                    paste0(file_base, "_enrichment_results.csv")))
}

# 合并所有结果
all_results_combined <- bind_rows(all_enrichment_results)

# 保存合并结果（简化文件名）
write_csv(all_results_combined, file.path(output_dir, "00_All_GMT_combined.csv"))
write_csv(file_summary, file.path(output_dir, "00_GMT_summary.csv"))

cat("\n==================== 批量处理完成 ====================\n")
for (i in 1:nrow(file_summary)) {
  cat(sprintf("  %s: 总=%d, 转录交集=%d, 分泌交集=%d, 两层均有=%d\n",
              file_summary$File[i],
              file_summary$Total_Pathways[i],
              file_summary$Pathways_With_Trans_Overlap[i],
              file_summary$Pathways_With_Sec_Overlap[i],
              file_summary$Pathways_With_Both_Overlap[i]))
}

# -------------------- 8. 筛选重点通路 --------------------
cat("\n==================== 筛选重点通路 ====================\n")

priority_pathways <- all_results_combined %>%
  filter(Both_Have_Overlap | Trans_Significant) %>%
  arrange(File, desc(Both_Have_Overlap), Trans_P_value)

cat("重点通路总数:", nrow(priority_pathways), "\n")

if (nrow(priority_pathways) > 0) {
  write_csv(priority_pathways, file.path(output_dir, "01_Priority_pathways.csv"))
}

# -------------------- 9. 可视化：两层组学均显著富集通路气泡图（P < 0.05） --------------------
cat("\n==================== 生成两层组学均显著富集通路气泡图 ====================\n")

# 筛选条件：分泌组显著（P < 0.05）且 转录组显著（P < 0.05）
sec_priority <- all_results_combined %>%
  filter(Has_Sec_Overlap, Has_Trans_Overlap,
         Sec_P_value < 0.05, Trans_P_value < 0.05) %>%
  group_by(Pathway) %>%
  summarise(
    Sec_Enrichment_Fold = max(Sec_Enrichment_Fold, na.rm = TRUE),
    Sec_Overlap_Count = max(Sec_Overlap_Count, na.rm = TRUE),
    Sec_P_value = min(Sec_P_value, na.rm = TRUE),
    Trans_Enrichment_Fold = max(Trans_Enrichment_Fold, na.rm = TRUE),
    Trans_Overlap_Count = max(Trans_Overlap_Count, na.rm = TRUE),
    Trans_P_value = min(Trans_P_value, na.rm = TRUE),
    Files = paste(unique(File), collapse = "; "),
    .groups = 'drop'
  ) %>%
  # 按分泌组富集倍数降序排列
  arrange(desc(Sec_Enrichment_Fold)) %>%
  mutate(
    Trans_LogP = -log10(Trans_P_value + 1e-300),
    Sec_LogP = -log10(Sec_P_value + 1e-300),
    Trans_LogP = ifelse(is.infinite(Trans_LogP) | Trans_LogP > 10, 10, Trans_LogP),
    Sec_LogP = ifelse(is.infinite(Sec_LogP) | Sec_LogP > 10, 10, Sec_LogP),
    Significance = case_when(
      Sec_P_value < 0.01 & Trans_P_value < 0.01 ~ "**",
      Sec_P_value < 0.05 & Trans_P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

if (nrow(sec_priority) > 0) {
  
  n_pathways <- nrow(sec_priority)
  cat("共找到", n_pathways, "个两层组学均显著富集的通路（P < 0.05）\n")
  
  # 准备绘图数据
  trans_data <- sec_priority %>%
    select(Pathway, 
           Overlap_Count = Trans_Overlap_Count, 
           LogP = Trans_LogP, 
           Enrichment_Fold = Trans_Enrichment_Fold) %>%
    mutate(Source = "Transcriptome")
  
  sec_data <- sec_priority %>%
    select(Pathway, 
           Overlap_Count = Sec_Overlap_Count, 
           LogP = Sec_LogP, 
           Enrichment_Fold = Sec_Enrichment_Fold) %>%
    mutate(Source = "Secretome")
  
  plot_data <- bind_rows(trans_data, sec_data) %>%
    mutate(
      Pathway = factor(Pathway, levels = rev(sec_priority$Pathway)),
      Source = factor(Source, levels = c("Transcriptome", "Secretome"))
    )
  
  # 获取实际最大交集数
  actual_max_overlap <- max(plot_data$Overlap_Count, na.rm = TRUE)
  
  # 手动生成图例断点（确保从 1 开始）
  if (actual_max_overlap <= 5) {
    my_breaks <- 1:actual_max_overlap
  } else if (actual_max_overlap <= 10) {
    my_breaks <- c(1, 2, 5, actual_max_overlap)
  } else if (actual_max_overlap <= 20) {
    my_breaks <- c(1, 5, 10, 15, actual_max_overlap)
  } else {
    my_breaks <- c(1, 10, 20, 30, actual_max_overlap)
  }
  my_breaks <- sort(unique(my_breaks))
  
  # 根据通路数量动态调整图片高度
  plot_height <- max(6, n_pathways * 0.45)
  
  # 绘制气泡图
  dot_plot <- ggplot(plot_data, 
                     aes(x = Source, y = Pathway, size = Overlap_Count, color = LogP)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(
      name = "Overlap\nGene Count",
      range = c(3, 12),
      breaks = my_breaks,
      limits = c(1, actual_max_overlap)
    ) +
    scale_color_gradientn(
      name = expression(-log[10](P)),
      colors = c("#FDE0DD", "#F39B7F", "#E64B35", "#A50026"),
      limits = c(0, max(plot_data$LogP, na.rm = TRUE))
    ) +
    facet_grid(. ~ Source, scales = "free_x", space = "free_x") +
    labs(
      title = "Pathways Significantly Enriched in Both Omics Layers",
      subtitle = paste0(n_pathways, " pathways with P < 0.05 in both secretome and transcriptome"),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      strip.text = element_text(face = "plain", size = 12, color = "white"),
      strip.background = element_rect(fill = "#2C3E50"),
      axis.text.y = element_text(size = max(8, 11 - n_pathways * 0.03), 
                                 face = "plain", color = "#2C3E50"),
      axis.text.x = element_text(size = 12, color = "#34495E"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line.x = element_line(color = "#BDC3C7", linewidth = 0.3),
      plot.title = element_text(face = "plain", hjust = 0.5, size = plot_config$title_size),
      plot.subtitle = element_text(hjust = 0.5, color = "#7F8C8D", size = 10),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(face = "plain", size = 12),
      panel.spacing = unit(1, "lines")
    )
  
  ggsave(file.path(output_dir, "02_Both_Omics_Significant_Pathways_dotplot.pdf"),
         dot_plot, width = 8, height = plot_height, dpi = plot_config$dpi, limitsize = FALSE)
  ggsave(file.path(output_dir, "02_Both_Omics_Significant_Pathways_dotplot.png"),
         dot_plot, width = 8, height = plot_height, dpi = plot_config$dpi, bg = "white", limitsize = FALSE)
  
  cat("气泡图已保存（高度:", round(plot_height, 1), "英寸）。\n")
  
  # 保存数据表
  sec_priority_output <- sec_priority %>%
    select(Pathway, 
           Sec_Enrichment_Fold, Sec_Overlap_Count, Sec_P_value,
           Trans_Enrichment_Fold, Trans_Overlap_Count, Trans_P_value,
           Significance, Files) %>%
    arrange(desc(Sec_Enrichment_Fold))
  
  write_csv(sec_priority_output, file.path(output_dir, "02_Both_Omics_Significant_Pathways_data.csv"))
  
  cat("\n--- 两层组学均显著富集通路（P < 0.05，共", n_pathways, "个）---\n")
  for (i in 1:nrow(sec_priority_output)) {
    sig_mark <- sec_priority_output$Significance[i]
    cat(sprintf("  %d. %s %s: 分泌(倍数=%.2f, 基因=%d, P=%.4f) | 转录(倍数=%.2f, 基因=%d, P=%.4f)\n",
                i,
                sec_priority_output$Pathway[i],
                sig_mark,
                sec_priority_output$Sec_Enrichment_Fold[i],
                sec_priority_output$Sec_Overlap_Count[i],
                sec_priority_output$Sec_P_value[i],
                sec_priority_output$Trans_Enrichment_Fold[i],
                sec_priority_output$Trans_Overlap_Count[i],
                sec_priority_output$Trans_P_value[i]))
  }
  
} else {
  cat("警告：没有找到两层组学均显著富集的通路（P < 0.05）。\n")
  
  # 输出当前数据分布情况供参考
  cat("\n数据分布统计:\n")
  
  cat("\n1. 分泌组有交集且显著的统计:\n")
  sec_stats <- all_results_combined %>%
    filter(Has_Sec_Overlap) %>%
    summarise(
      Total = n(),
      P_lt_05 = sum(Sec_P_value < 0.05, na.rm = TRUE),
      P_lt_01 = sum(Sec_P_value < 0.1, na.rm = TRUE),
      Min_P = min(Sec_P_value, na.rm = TRUE)
    )
  cat("  总通路数:", sec_stats$Total, "\n")
  cat("  P < 0.05:", sec_stats$P_lt_05, "\n")
  cat("  P < 0.1:", sec_stats$P_lt_01, "\n")
  cat("  最小 P 值:", sec_stats$Min_P, "\n")
  
  cat("\n2. 转录组有交集且显著的统计:\n")
  trans_stats <- all_results_combined %>%
    filter(Has_Trans_Overlap) %>%
    summarise(
      Total = n(),
      P_lt_05 = sum(Trans_P_value < 0.05, na.rm = TRUE),
      P_lt_01 = sum(Trans_P_value < 0.1, na.rm = TRUE),
      Min_P = min(Trans_P_value, na.rm = TRUE)
    )
  cat("  总通路数:", trans_stats$Total, "\n")
  cat("  P < 0.05:", trans_stats$P_lt_05, "\n")
  cat("  P < 0.1:", trans_stats$P_lt_01, "\n")
  cat("  最小 P 值:", trans_stats$Min_P, "\n")
  
  cat("\n3. 两层均有交集的统计:\n")
  both_stats <- all_results_combined %>%
    filter(Has_Sec_Overlap, Has_Trans_Overlap) %>%
    summarise(
      Total = n(),
      Both_P_lt_05 = sum(Sec_P_value < 0.05 & Trans_P_value < 0.05, na.rm = TRUE),
      Sec_P_lt_05 = sum(Sec_P_value < 0.05, na.rm = TRUE),
      Trans_P_lt_05 = sum(Trans_P_value < 0.05, na.rm = TRUE)
    )
  cat("  总通路数:", both_stats$Total, "\n")
  cat("  两层均显著(P<0.05):", both_stats$Both_P_lt_05, "\n")
  cat("  仅分泌组显著:", both_stats$Sec_P_lt_05 - both_stats$Both_P_lt_05, "\n")
  cat("  仅转录组显著:", both_stats$Trans_P_lt_05 - both_stats$Both_P_lt_05, "\n")
}
# -------------------- 10. 补充条形图：两层均显著通路的富集倍数对比 --------------------
if (nrow(sec_priority) > 0) {
  
  n_pathways <- nrow(sec_priority)
  plot_height <- max(6, n_pathways * 0.45)
  
  # 准备对比数据
  bar_data <- sec_priority %>%
    select(Pathway, Sec_Enrichment_Fold, Trans_Enrichment_Fold) %>%
    pivot_longer(cols = c(Sec_Enrichment_Fold, Trans_Enrichment_Fold),
                 names_to = "Omics",
                 values_to = "Enrichment_Fold") %>%
    mutate(
      Omics = ifelse(grepl("Sec", Omics), "Secretome", "Transcriptome"),
      Pathway = factor(Pathway, levels = rev(sec_priority$Pathway))
    )
  
  bar_plot <- ggplot(bar_data, 
                     aes(x = Enrichment_Fold, y = Pathway, fill = Omics)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#7F8C8D", alpha = 0.6) +
    scale_fill_manual(
      name = "Omics Layer",
      values = c("Secretome" = plot_config$color_sec, 
                 "Transcriptome" = plot_config$color_trans)
    ) +
    labs(
      title = "Enrichment Fold Comparison: Both Omics Significant Pathways",
      subtitle = paste0(n_pathways, " pathways with P < 0.05 in both layers"),
      x = "Enrichment Fold", y = NULL
    ) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      axis.text.y = element_text(size = max(8, 11 - n_pathways * 0.05), 
                                 face = "bold", color = "#2C3E50"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = plot_config$title_size),
      plot.subtitle = element_text(hjust = 0.5, color = "#7F8C8D", size = 9),
      legend.position = "top"
    )
  
  ggsave(file.path(output_dir, "03_Both_Omics_Enrichment_Comparison_barplot.pdf"),
         bar_plot, width = 10, height = plot_height, dpi = plot_config$dpi, limitsize = FALSE)
  ggsave(file.path(output_dir, "03_Both_Omics_Enrichment_Comparison_barplot.png"),
         bar_plot, width = 10, height = plot_height, dpi = plot_config$dpi, bg = "white", limitsize = FALSE)
  
  cat("富集倍数对比条形图已保存。\n")
}

# -------------------- 11. 详细报告 --------------------
cat("\n==================== 生成详细报告 ====================\n")

sec_overlap_detailed <- all_results_combined %>%
  filter(Has_Sec_Overlap) %>%
  select(File, Pathway, Total_Genes_in_Pathway,
         Sec_Overlap_Count, Sec_Overlap_Genes,
         Trans_Overlap_Count, Trans_P_value, Both_Have_Overlap) %>%
  arrange(File, desc(Both_Have_Overlap), Trans_P_value)

write_csv(sec_overlap_detailed, file.path(output_dir, "04_Secreted_pathway_mapping.csv"))

# 基因-通路映射
gene_pathway_mapping <- sec_overlap_detailed %>%
  filter(Sec_Overlap_Count > 0) %>%
  select(File, Pathway, Sec_Overlap_Genes) %>%
  separate_rows(Sec_Overlap_Genes, sep = ";") %>%
  rename(GeneSymbol = Sec_Overlap_Genes) %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  distinct()

write_csv(gene_pathway_mapping, file.path(output_dir, "05_Gene_to_Pathway_mapping.csv"))

# -------------------- 12. 分析摘要 --------------------
sink(file.path(output_dir, "06_Analysis_summary.txt"))
cat("========== GMT 批量靶向通路富集分析摘要 ==========\n\n")
cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("--- 输入数据 ---\n")
cat("转录组差异基因数:", length(trans_genes), "\n")
cat("差异分泌蛋白数:", length(sec_genes), "\n")
cat("背景基因总数:", length(bg_genes), "\n")
cat("GMT 文件数:", length(gmt_files), "\n\n")

cat("--- GMT 文件处理结果 ---\n")
for (i in 1:nrow(file_summary)) {
  cat(sprintf("  %s: 总通路=%d, 转录交集=%d, 分泌交集=%d, 两层均有=%d\n",
              file_summary$File[i],
              file_summary$Total_Pathways[i],
              file_summary$Pathways_With_Trans_Overlap[i],
              file_summary$Pathways_With_Sec_Overlap[i],
              file_summary$Pathways_With_Both_Overlap[i]))
}

cat("\n--- 重点通路（两层均有交集）---\n")
priority_both <- all_results_combined %>%
  filter(Both_Have_Overlap) %>%
  arrange(Trans_P_value) %>%
  select(File, Pathway, Trans_Overlap_Count, Sec_Overlap_Count, Trans_P_value)

if (nrow(priority_both) > 0) {
  for (i in 1:min(20, nrow(priority_both))) {
    cat(sprintf("  %d. [%s] %s: 转录组=%d, 分泌组=%d, P=%.3f\n",
                i,
                priority_both$File[i],
                priority_both$Pathway[i],
                priority_both$Trans_Overlap_Count[i],
                priority_both$Sec_Overlap_Count[i],
                priority_both$Trans_P_value[i]))
  }
} else {
  cat("  无两层组学均有交集的通路。\n")
}

sink()

cat("\n==================== 分析完成 ====================\n")
cat("所有结果已保存至:", output_dir, "\n")