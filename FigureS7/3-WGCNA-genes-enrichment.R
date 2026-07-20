# ============================================================================
# WGCNA分析 - 模块3：模块与功能基因集重叠分析（修复版）
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

# 创建子目录
gsea_dir <- file.path(output_dir, "gene_set_enrichment")
if (!dir.exists(gsea_dir)) {
  dir.create(gsea_dir, recursive = TRUE)
}

# ========================== 参数设置 =========================================
# GMT文件模式
gmt_pattern <- "\\.gmt$"
gmt_recursive <- FALSE

# 统计阈值
pvalue_threshold <- 0.05
top_n_modules <- 2

# 输出文件
output_rdata <- file.path(output_dir, "03_module_enrichment.RData")
output_enrichment <- file.path(gsea_dir, "module_gene_set_enrichment.csv")
output_top_modules <- file.path(gsea_dir, "top_enriched_modules.csv")

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

library(dplyr)
library(tidyr)

# ========================== 加载WGCNA结果 ====================================
cat("\n========== 加载WGCNA结果 ==========\n")

load(file.path(output_dir, "02_network_modules.RData"))
load(file.path(output_dir, "01_data_prepared.RData"))

# 获取每个模块的基因列表
module_genes <- list()
for (mod in unique(moduleColors)) {
  module_genes[[mod]] <- rownames(expr_matrix_norm)[moduleColors == mod]
}

cat("模块基因数统计:\n")
for (mod in names(module_genes)) {
  cat(sprintf("  %s: %d 个基因\n", mod, length(module_genes[[mod]])))
}

# ========================== 查找并读取GMT文件 ====================================
cat("\n========== 查找并读取GMT文件 ==========\n")

# 获取所有GMT文件
gmt_files <- list.files(path = ".", 
                        pattern = gmt_pattern, 
                        recursive = gmt_recursive, 
                        ignore.case = TRUE,
                        full.names = TRUE)

if (length(gmt_files) == 0) {
  cat("未找到GMT文件，使用示例基因集\n")
  # 创建示例基因集
  all_gene_sets <- list(
    "Epithelial_Mesenchymal_Transition" = c("SNAI1", "SNAI2", "ZEB1", "VIM", "CDH2", "FN1", "MMP2", "MMP9"),
    "Proliferation" = c("MKI67", "PCNA", "CCNB1", "CCND1", "CDK1", "CDK2"),
    "Apoptosis" = c("BAX", "BCL2", "CASP3", "CASP9", "P53")
  )
  gene_set_info <- data.frame(
    GeneSet = names(all_gene_sets),
    Description = names(all_gene_sets),
    Source = "Example",
    Size = sapply(all_gene_sets, length),
    stringsAsFactors = FALSE
  )
} else {
  cat(sprintf("找到 %d 个GMT文件:\n", length(gmt_files)))
  for (i in 1:length(gmt_files)) {
    cat(sprintf("  %d. %s\n", i, basename(gmt_files[i])))
  }
  
  # 读取所有GMT文件
  all_gene_sets <- list()
  gene_set_info <- data.frame()
  
  for (gmt_file in gmt_files) {
    file_name <- basename(gmt_file)
    cat(sprintf("\n读取文件: %s\n", file_name))
    
    gmt_lines <- readLines(gmt_file, warn = FALSE)
    
    for (line in gmt_lines) {
      if (nchar(trimws(line)) == 0) next
      
      parts <- strsplit(line, "\t")[[1]]
      
      if (length(parts) >= 3) {
        gene_set_name <- parts[1]
        gene_set_desc <- parts[2]
        genes <- toupper(parts[3:length(parts)])
        genes <- genes[genes != "" & !is.na(genes)]
        
        if (length(genes) > 0) {
          all_gene_sets[[gene_set_name]] <- genes
          gene_set_info <- rbind(gene_set_info, data.frame(
            GeneSet = gene_set_name,
            Description = gene_set_desc,
            Source = file_name,
            Size = length(genes),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  cat(sprintf("\n总共读取 %d 个基因集\n", length(all_gene_sets)))
}

# ========================== 超几何检验 ====================================
cat("\n========== 计算模块与基因集的重叠显著性 ==========\n")

# 背景基因数
background_size <- nrow(expr_matrix_norm)

# 存储富集结果
enrichment_results <- data.frame()

for (mod in names(module_genes)) {
  module_genes_list <- module_genes[[mod]]
  module_size <- length(module_genes_list)
  
  for (gs_name in names(all_gene_sets)) {
    gs_genes <- all_gene_sets[[gs_name]]
    gs_size <- length(gs_genes)
    
    # 计算重叠
    overlap <- intersect(module_genes_list, gs_genes)
    overlap_size <- length(overlap)
    
    if (overlap_size > 0) {
      # 超几何检验
      p_value <- phyper(overlap_size - 1, gs_size, background_size - gs_size, 
                        module_size, lower.tail = FALSE)
      
      # 富集倍数
      expected <- (module_size / background_size) * gs_size
      fold_enrichment <- overlap_size / expected
      
      enrichment_results <- rbind(enrichment_results, data.frame(
        Module = mod,
        Module_Size = module_size,
        GeneSet = gs_name,
        GeneSet_Size = gs_size,
        GeneSet_Source = ifelse(exists("gene_set_info"), 
                                gene_set_info$Source[gene_set_info$GeneSet == gs_name], 
                                "Example"),
        Overlap_Size = overlap_size,
        Overlap_Genes = paste(overlap, collapse = "; "),
        Expected = expected,
        Fold_Enrichment = fold_enrichment,
        PValue = p_value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 计算FDR
if (nrow(enrichment_results) > 0) {
  enrichment_results$FDR <- p.adjust(enrichment_results$PValue, method = "BH")
  enrichment_results <- enrichment_results %>% arrange(PValue, desc(Fold_Enrichment))
  enrichment_results$Significant <- enrichment_results$FDR < pvalue_threshold
  enrichment_results$Signif_mark <- ifelse(enrichment_results$FDR < 0.001, "***",
                                           ifelse(enrichment_results$FDR < 0.01, "**",
                                                  ifelse(enrichment_results$FDR < 0.05, "*", "ns")))
}

# 保存完整结果
write.csv(enrichment_results, output_enrichment, row.names = FALSE)
cat(sprintf("富集分析结果已保存: %s (%d 行)\n", output_enrichment, nrow(enrichment_results)))

# ========================== 选择Top显著模块 ====================================
cat("\n========== 选择Top显著模块 ==========\n")

if (nrow(enrichment_results) > 0) {
  # 按模块汇总显著富集的基因集数量
  module_summary <- enrichment_results %>%
    filter(Significant == TRUE) %>%
    group_by(Module) %>%
    summarise(
      Significant_GS_Count = n(),
      Max_Fold_Enrichment = max(Fold_Enrichment, na.rm = TRUE),
      Min_FDR = min(FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Significant_GS_Count), Min_FDR)
  
  cat("\n模块富集摘要:\n")
  print(module_summary)
  
  # 选择Top模块
  if (nrow(module_summary) > 0) {
    top_modules <- head(module_summary$Module, top_n_modules)
  } else {
    # 如果没有显著富集，选择模块大小最大的
    module_sizes <- sapply(module_genes, length)
    top_modules <- names(sort(module_sizes, decreasing = TRUE))[1:min(top_n_modules, length(module_sizes))]
  }
} else {
  # 如果没有富集结果，选择模块大小最大的
  module_sizes <- sapply(module_genes, length)
  top_modules <- names(sort(module_sizes, decreasing = TRUE))[1:min(top_n_modules, length(module_sizes))]
}

cat(sprintf("\n选择的Top %d 模块: %s\n", top_n_modules, paste(top_modules, collapse = ", ")))

# 保存Top模块信息
write.csv(module_summary, output_top_modules, row.names = FALSE)

# ========================== 保存结果 ====================================
cat("\n========== 保存结果 ==========\n")

save(enrichment_results, module_summary, top_modules, all_gene_sets,
     file = output_rdata)

cat(sprintf("富集分析结果已保存: %s\n", output_rdata))

cat("\n模块3完成！\n")
cat(sprintf("选择的Top模块: %s\n", paste(top_modules, collapse = ", ")))