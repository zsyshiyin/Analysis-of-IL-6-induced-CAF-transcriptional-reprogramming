# ============================================================================
# 差异基因与GSEA通路匹配分析（自动识别GMT文件）+ 双基因集Venn图
# ============================================================================

# 清空环境变量
rm(list = ls())

# ========================== 工作路径设定 ====================================
# 请修改为您的实际工作路径
work_dir <- "D:/zsy/SX/Fomal-final/5-CIPS-shaixuan"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")

# ========================== 参数设置 =========================================
# 输入文件参数
deg_file <- "I6-vs-Con.known.up.DEG.csv"                # 差异基因列表文件

# GMT文件参数
gmt_pattern <- "\\.gmt$"                     # GMT文件匹配模式（正则表达式）
gmt_recursive <- FALSE                        # 是否在子目录中搜索GMT文件
gmt_ignore_case <- TRUE                        # 是否忽略大小写

# 输出文件参数
output_dir <- "pathway_mapping_results_Venn"      # 输出文件夹
output_file <- "gene_pathway_mapping.csv"    # 主输出文件
summary_file <- "pathway_summary.csv"        # 通路汇总文件
gene_stats_file <- "gene_pathway_stats.csv"  # 基因统计文件
unmatched_file <- "unmatched_genes.csv"      # 未匹配基因文件

# 列名自定义参数（差异基因文件）
deg_gene_col <- "GeneID"                      # 基因ID列名
deg_logFC_col <- "logFC"                       # logFC列名
deg_pval_col <- "PValue"                       # P值列名
deg_fdr_col <- "FDR"                           # FDR列名

# 筛选参数
pval_threshold <- 0.05                         # P值阈值
fdr_threshold <- 0.05                          # FDR阈值
min_logFC <- 0                                  # 最小logFC绝对值（0表示不筛选）
use_adjusted_pvalue <- TRUE                     # 是否使用调整后P值（FDR）进行筛选

# ===== 新增：数据处理参数 =====
remove_duplicates <- TRUE                       # 是否移除重复的基因-通路对
include_all_genes <- FALSE                      # 是否包含未匹配到通路的基因
include_pathway_source <- TRUE                   # 是否包含通路来源信息
include_pathway_description <- TRUE              # 是否包含通路描述
verbose <- TRUE                                  # 是否显示详细过程
save_detailed_stats <- TRUE                      # 是否保存详细统计信息

# 基因ID处理参数
gene_id_type <- "symbol"                         # 基因ID类型（symbol/ensembl/entrez）
remove_version <- TRUE                           # 是否移除Ensembl ID的版本号（如 ENSG000001.4 -> ENSG000001）
standardize_gene_ids <- TRUE                      # 是否标准化基因ID（转换为大写）

# ========================== Venn图参数（新增） ================================
# 是否生成Venn图
generate_venn <- TRUE                            # 是否生成Venn图

# Venn图数据集名称
venn_name1 <- "CAFs_vs_NFs_UP_DEGs"              # 数据集1名称（差异基因）
venn_name2 <- "Invasion-related gene set genes"    # 数据集2名称（GMT所有基因）

# Venn图外观参数
venn_width <- 6                                  # 图片宽度（英寸）
venn_height <- 3                                 # 图片高度（英寸）
venn_title <- "CIPS screening"        # 标题
venn_title_size <- 16                            # 标题字体大小（磅）
venn_cat_cex <- 1.2                              # 分类标签字体大小
venn_num_cex <- 2                              # 数字字体大小
venn_color1 <- "#E41A1C"                         # 数据集1颜色（红色）
venn_color2 <- "#377EB8"                         # 数据集2颜色（蓝色）
venn_alpha <- 0.3                                # 透明度
venn_lty <- "blank"                              # 边框线型
venn_fontface <- "plain"                          # 字体样式
venn_cat_dist <- 0.05                            # 分类标签距离
venn_margin <- 0.05                              # 边距

# ========================== 检查必要包 ======================================
cat("\n========== 检查必要包 ==========\n")

required_packages <- c("dplyr", "tidyr", "stringr", "tools", "VennDiagram", "ggplot2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat(sprintf("安装CRAN包: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  } else {
    cat(sprintf("✓ %s 已安装\n", pkg))
  }
}

# ========================== 创建输出目录 ====================================
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# ========================== 读取差异基因数据 ================================
cat("\n========== 读取差异基因数据 ==========\n")

if (!file.exists(deg_file)) {
  stop("差异基因文件不存在: ", deg_file)
}

# 读取差异基因数据
deg_data <- read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("原始差异基因数据维度: %d 行, %d 列\n", nrow(deg_data), ncol(deg_data)))
cat("列名:", paste(colnames(deg_data), collapse = ", "), "\n")

# 检查必要的列是否存在
required_cols <- c(deg_gene_col, deg_logFC_col, deg_pval_col, deg_fdr_col)
missing_cols <- required_cols[!required_cols %in% colnames(deg_data)]

if (length(missing_cols) > 0) {
  cat("警告: 以下指定的列名不存在:", paste(missing_cols, collapse = ", "), "\n")
  cat("可用的列名:", paste(colnames(deg_data), collapse = ", "), "\n")
  
  # 尝试自动匹配
  cat("\n尝试自动匹配列名...\n")
  
  # 基因列匹配
  gene_col_candidates <- grep("gene|symbol|id|name", colnames(deg_data), ignore.case = TRUE, value = TRUE)
  if (length(gene_col_candidates) > 0) {
    deg_gene_col <- gene_col_candidates[1]
    cat(sprintf("  自动匹配基因列为: %s\n", deg_gene_col))
  }
  
  # logFC列匹配
  logFC_candidates <- grep("logfc|log2fc|foldchange|fc", colnames(deg_data), ignore.case = TRUE, value = TRUE)
  if (length(logFC_candidates) > 0) {
    deg_logFC_col <- logFC_candidates[1]
    cat(sprintf("  自动匹配logFC列为: %s\n", deg_logFC_col))
  }
  
  # P值列匹配
  pval_candidates <- grep("pvalue|pval|p.val|p_value|p value", colnames(deg_data), ignore.case = TRUE, value = TRUE)
  if (length(pval_candidates) > 0) {
    deg_pval_col <- pval_candidates[1]
    cat(sprintf("  自动匹配P值列为: %s\n", deg_pval_col))
  }
  
  # FDR列匹配
  fdr_candidates <- grep("fdr|padj|adj.p|qvalue|q.val", colnames(deg_data), ignore.case = TRUE, value = TRUE)
  if (length(fdr_candidates) > 0) {
    deg_fdr_col <- fdr_candidates[1]
    cat(sprintf("  自动匹配FDR列为: %s\n", deg_fdr_col))
  }
}

# 提取需要的列
deg_selected <- deg_data[, c(deg_gene_col, deg_logFC_col, deg_pval_col, deg_fdr_col)]
colnames(deg_selected) <- c("GeneID", "logFC", "PValue", "FDR")

# 处理基因ID
cat("\n处理基因ID...\n")

# 移除NA值
deg_selected <- deg_selected[!is.na(deg_selected$GeneID) & deg_selected$GeneID != "", ]

# 标准化基因ID（转换为大写）
if (standardize_gene_ids) {
  deg_selected$GeneID <- toupper(trimws(as.character(deg_selected$GeneID)))
  cat("✓ 基因ID已转换为大写\n")
}

# 移除Ensembl ID版本号
if (remove_version && gene_id_type == "ensembl") {
  deg_selected$GeneID <- gsub("\\.[0-9]+$", "", deg_selected$GeneID)
  cat("✓ 已移除Ensembl ID版本号\n")
}

# 筛选差异基因
cat("\n应用筛选阈值:\n")
cat(sprintf("  P值阈值: %.3f\n", pval_threshold))
cat(sprintf("  FDR阈值: %.3f\n", fdr_threshold))
cat(sprintf("  |logFC|最小值: %.2f\n", min_logFC))

if (use_adjusted_pvalue) {
  deg_filtered <- deg_selected %>%
    filter(FDR <= fdr_threshold) %>%
    filter(abs(logFC) >= min_logFC)
} else {
  deg_filtered <- deg_selected %>%
    filter(PValue <= pval_threshold) %>%
    filter(abs(logFC) >= min_logFC)
}

cat(sprintf("\n筛选后差异基因数量: %d\n", nrow(deg_filtered)))

# 显示上下调基因数量
up_genes <- sum(deg_filtered$logFC > 0)
down_genes <- sum(deg_filtered$logFC < 0)
cat(sprintf("  上调基因: %d\n", up_genes))
cat(sprintf("  下调基因: %d\n", down_genes))

# ========================== 查找GMT文件 ====================================
cat("\n========== 查找GMT文件 ==========\n")

# 获取当前工作目录下的所有GMT文件
gmt_files <- list.files(path = ".", 
                        pattern = gmt_pattern, 
                        recursive = gmt_recursive, 
                        ignore.case = gmt_ignore_case,
                        full.names = TRUE)

if (length(gmt_files) == 0) {
  stop("未找到任何GMT文件！请确保工作目录下存在.gmt文件")
}

cat(sprintf("找到 %d 个GMT文件:\n", length(gmt_files)))
for (i in 1:length(gmt_files)) {
  cat(sprintf("  %d. %s\n", i, basename(gmt_files[i])))
}

# ========================== 读取GMT文件 ====================================
cat("\n========== 读取GMT文件 ==========\n")

# 初始化通路数据框
all_pathways <- data.frame(
  PathwayID = character(),
  PathwayName = character(),
  Source = character(),
  Genes = character(),
  stringsAsFactors = FALSE
)

# 用于存储所有GMT基因的集合
all_gmt_genes <- c()

# 读取每个GMT文件
for (gmt_file in gmt_files) {
  file_name <- basename(gmt_file)
  cat(sprintf("\n读取文件: %s\n", file_name))
  
  # 读取GMT文件
  gmt_lines <- readLines(gmt_file, warn = FALSE)
  cat(sprintf("  包含 %d 个通路\n", length(gmt_lines)))
  
  # 解析每一行
  file_pathways <- data.frame()
  
  for (line in gmt_lines) {
    # 跳过空行
    if (nchar(trimws(line)) == 0) next
    
    # 分割行（GMT格式：通路名\t描述\t基因1\t基因2\t...）
    parts <- strsplit(line, "\t")[[1]]
    
    if (length(parts) >= 3) {
      pathway_id <- parts[1]
      pathway_desc <- parts[2]
      genes <- parts[3:length(parts)]
      
      # 移除空基因
      genes <- genes[genes != "" & !is.na(genes)]
      
      if (length(genes) > 0) {
        # 标准化基因ID
        if (standardize_gene_ids) {
          genes <- toupper(trimws(genes))
        }
        
        # 添加到所有GMT基因集合
        all_gmt_genes <- c(all_gmt_genes, genes)
        
        file_pathways <- rbind(file_pathways, data.frame(
          PathwayID = pathway_id,
          PathwayName = pathway_desc,
          Source = file_name,
          Genes = I(list(genes)),
          GeneCount = length(genes),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  cat(sprintf("  成功解析 %d 个通路\n", nrow(file_pathways)))
  all_pathways <- rbind(all_pathways, file_pathways)
}

# 去重GMT基因
all_gmt_genes <- unique(all_gmt_genes)
cat(sprintf("\n总计: %d 个通路，来自 %d 个GMT文件\n", 
            nrow(all_pathways), length(unique(all_pathways$Source))))
cat(sprintf("GMT所有基因集去重后基因总数: %d\n", length(all_gmt_genes)))

# ========================== Venn图分析（使用差异基因 vs GMT所有基因） ================================
if (generate_venn) {
  cat("\n========== Venn图分析 ==========\n")
  
  # 提取上调基因列表作为第一个数据集
  genes_set1 <- deg_filtered$GeneID[deg_filtered$logFC > 0]
  cat(sprintf("数据集1（%s）: %d 个基因\n", venn_name1, length(genes_set1)))
  
  # 第二个数据集：GMT所有基因
  genes_set2 <- all_gmt_genes
  cat(sprintf("数据集2（%s）: %d 个基因\n", venn_name2, length(genes_set2)))
  
  # 计算交集和差集
  overlap_genes <- intersect(genes_set1, genes_set2)
  only_set1 <- setdiff(genes_set1, genes_set2)
  only_set2 <- setdiff(genes_set2, genes_set1)
  
  cat(sprintf("\nVenn图统计:\n"))
  cat(sprintf("  数据集1（%s）: %d 个基因\n", venn_name1, length(genes_set1)))
  cat(sprintf("  数据集2（%s）: %d 个基因\n", venn_name2, length(genes_set2)))
  cat(sprintf("  重合基因数: %d\n", length(overlap_genes)))
  cat(sprintf("  数据集1特有: %d\n", length(only_set1)))
  cat(sprintf("  数据集2特有: %d\n", length(only_set2)))
  
  # 保存Venn图结果
  if (length(overlap_genes) > 0) {
    overlap_df <- data.frame(GeneID = overlap_genes, stringsAsFactors = FALSE)
    write.csv(overlap_df, file.path(output_dir, "venn_overlap_genes.csv"), row.names = FALSE)
    cat(sprintf("✓ 重合基因已保存: %s\n", file.path(output_dir, "venn_overlap_genes.csv")))
  }
  
  # 保存数据集1特有基因
  if (length(only_set1) > 0) {
    only1_df <- data.frame(GeneID = only_set1, stringsAsFactors = FALSE)
    write.csv(only1_df, file.path(output_dir, "venn_unique_set1_genes.csv"), row.names = FALSE)
    cat(sprintf("✓ 数据集1特有基因已保存: %s\n", file.path(output_dir, "venn_unique_set1_genes.csv")))
  }
  
  # 保存数据集2特有基因
  if (length(only_set2) > 0) {
    only2_df <- data.frame(GeneID = only_set2, stringsAsFactors = FALSE)
    write.csv(only2_df, file.path(output_dir, "venn_unique_set2_genes.csv"), row.names = FALSE)
    cat(sprintf("✓ 数据集2特有基因已保存: %s\n", file.path(output_dir, "venn_unique_set2_genes.csv")))
  }
  
  # 绘制Venn图
  cat("\n绘制Venn图...\n")
  
  venn_output_pdf <- file.path(output_dir, "venn_diagram.pdf")
  venn_output_png <- file.path(output_dir, "venn_diagram.png")
  
  # 绘制PDF
  pdf(venn_output_pdf, width = venn_width, height = venn_height)
  grid.newpage()
  
  venn_plot <- draw.pairwise.venn(
    area1 = length(genes_set1),
    area2 = length(genes_set2),
    cross.area = length(overlap_genes),
    category = c(venn_name1, venn_name2),
    fill = c(venn_color1, venn_color2),
    alpha = venn_alpha,
    lty = venn_lty,
    cex = venn_num_cex,
    cat.cex = venn_cat_cex,
    cat.pos = c(0, 0),
    cat.dist = c(venn_cat_dist, venn_cat_dist),
    cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
    ext.text = FALSE,
    scaled = FALSE,
    fontface = venn_fontface,
    fontfamily = "sans"
  )
  
  # 添加标题
  grid.text(venn_title, x = 0.5, y = 0.98, 
            gp = gpar(fontsize = venn_title_size, fontface = "bold"))
  
  dev.off()
  cat(sprintf("✓ Venn图已保存: %s\n", venn_output_pdf))
  
  # 绘制PNG
  png(venn_output_png, width = venn_width * 100, height = venn_height * 100, 
      res = 100, units = "px")
  grid.newpage()
  
  draw.pairwise.venn(
    area1 = length(genes_set1),
    area2 = length(genes_set2),
    cross.area = length(overlap_genes),
    category = c(venn_name1, venn_name2),
    fill = c(venn_color1, venn_color2),
    alpha = venn_alpha,
    lty = venn_lty,
    cex = venn_num_cex,
    cat.cex = venn_cat_cex,
    cat.pos = c(0, 0),
    cat.dist = c(venn_cat_dist, venn_cat_dist),
    cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
    ext.text = FALSE,
    scaled = FALSE,
    fontface = venn_fontface,
    fontfamily = "sans"
  )
  
  grid.text(venn_title, x = 0.5, y = 0.98, 
            gp = gpar(fontsize = venn_title_size, fontface = "bold"))
  
  dev.off()
  cat(sprintf("✓ Venn图(PNG)已保存: %s\n", venn_output_png))
}

# ========================== 匹配基因与通路 ====================================
cat("\n========== 匹配基因与通路 ==========\n")

# 创建基因-通路映射表
gene_pathway_mapping <- data.frame()

# 创建进度条
total_pathways <- nrow(all_pathways)
pb_width <- 50

for (i in 1:total_pathways) {
  # 显示进度
  if (verbose && i %% max(1, floor(total_pathways/20)) == 0) {
    progress <- floor(i / total_pathways * pb_width)
    cat(sprintf("\r[%-*s] %d/%d 通路", pb_width, 
                paste(rep("=", progress), collapse = ""), 
                i, total_pathways))
    flush.console()
  }
  
  pathway <- all_pathways[i, ]
  pathway_genes <- unlist(pathway$Genes)
  
  # 找出该通路中存在于差异基因列表中的基因
  matched_genes <- intersect(pathway_genes, deg_filtered$GeneID)
  
  if (length(matched_genes) > 0) {
    # 获取这些基因的差异表达信息
    gene_info <- deg_filtered[deg_filtered$GeneID %in% matched_genes, ]
    
    # 创建映射记录
    for (gene in matched_genes) {
      gene_data <- gene_info[gene_info$GeneID == gene, ]
      
      # 确保gene_data只有一行
      if (nrow(gene_data) > 1) {
        gene_data <- gene_data[1, ]
      }
      
      gene_pathway_mapping <- rbind(gene_pathway_mapping, data.frame(
        GeneID = gene,
        logFC = gene_data$logFC,
        PValue = gene_data$PValue,
        FDR = gene_data$FDR,
        Regulation = ifelse(gene_data$logFC > 0, "Up", "Down"),
        PathwayID = pathway$PathwayID,
        PathwayName = pathway$PathwayName,
        PathwaySource = pathway$Source,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n")  # 换行

# 移除重复的基因-通路对
if (remove_duplicates && nrow(gene_pathway_mapping) > 0) {
  original_count <- nrow(gene_pathway_mapping)
  gene_pathway_mapping <- gene_pathway_mapping %>%
    distinct(GeneID, PathwayID, .keep_all = TRUE)
  cat(sprintf("移除重复项: %d -> %d\n", original_count, nrow(gene_pathway_mapping)))
}

# 找出未匹配的基因
unmatched_genes <- setdiff(deg_filtered$GeneID, unique(gene_pathway_mapping$GeneID))

cat(sprintf("\n匹配结果:\n"))
cat(sprintf("  匹配上的基因数: %d\n", length(unique(gene_pathway_mapping$GeneID))))
cat(sprintf("  未匹配的基因数: %d\n", length(unmatched_genes)))
cat(sprintf("  匹配上的通路数: %d\n", length(unique(gene_pathway_mapping$PathwayID))))
cat(sprintf("  总映射关系数: %d\n", nrow(gene_pathway_mapping)))

# ========================== 生成统计信息 ====================================
cat("\n========== 生成统计信息 ==========\n")

if (nrow(gene_pathway_mapping) > 0) {
  
  # 1. 通路汇总统计
  pathway_summary <- gene_pathway_mapping %>%
    group_by(PathwayID, PathwayName, PathwaySource) %>%
    summarise(
      GeneCount = n_distinct(GeneID),
      UpGenes = sum(Regulation == "Up"),
      DownGenes = sum(Regulation == "Down"),
      AvgLogFC = mean(logFC, na.rm = TRUE),
      MinFDR = min(FDR, na.rm = TRUE),
      GeneList = paste(sort(unique(GeneID)), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(GeneCount))
  
  cat(sprintf("通路汇总: %d 个通路包含匹配基因\n", nrow(pathway_summary)))
  
  # 2. 基因统计（每个基因出现在多少个通路中）
  gene_stats <- gene_pathway_mapping %>%
    group_by(GeneID, logFC, PValue, FDR, Regulation) %>%
    summarise(
      PathwayCount = n_distinct(PathwayID),
      PathwayList = paste(sort(unique(PathwayName)), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(PathwayCount))
  
  cat(sprintf("基因统计: %d 个基因匹配到通路\n", nrow(gene_stats)))
  
  # 3. 按来源统计
  source_stats <- gene_pathway_mapping %>%
    group_by(PathwaySource) %>%
    summarise(
      PathwayCount = n_distinct(PathwayID),
      GeneCount = n_distinct(GeneID),
      MappingCount = n(),
      .groups = "drop"
    )
  
  cat("\n各GMT文件匹配统计:\n")
  print(source_stats)
  
} else {
  cat("没有找到任何匹配的基因-通路对！\n")
  cat("可能的原因:\n")
  cat("  1. 基因ID格式不匹配（如symbol vs ensembl）\n")
  cat("  2. 筛选阈值太严格\n")
  cat("  3. 通路文件中不包含您的差异基因\n")
}

# ========================== 保存结果 ====================================
cat("\n========== 保存结果 ==========\n")

if (nrow(gene_pathway_mapping) > 0) {
  
  # 保存主映射表
  write.csv(gene_pathway_mapping, 
            file = file.path(output_dir, output_file),
            row.names = FALSE, quote = FALSE)
  cat(sprintf("主映射表已保存: %s (%d行)\n", 
              file.path(output_dir, output_file), nrow(gene_pathway_mapping)))
  
  # 保存通路汇总
  write.csv(pathway_summary, 
            file = file.path(output_dir, summary_file),
            row.names = FALSE, quote = FALSE)
  cat(sprintf("通路汇总已保存: %s (%d行)\n", 
              file.path(output_dir, summary_file), nrow(pathway_summary)))
  
  # 保存基因统计
  write.csv(gene_stats, 
            file = file.path(output_dir, gene_stats_file),
            row.names = FALSE, quote = FALSE)
  cat(sprintf("基因统计已保存: %s (%d行)\n", 
              file.path(output_dir, gene_stats_file), nrow(gene_stats)))
  
  # 保存未匹配基因
  if (length(unmatched_genes) > 0) {
    unmatched_df <- deg_filtered[deg_filtered$GeneID %in% unmatched_genes, ]
    write.csv(unmatched_df, 
              file = file.path(output_dir, unmatched_file),
              row.names = FALSE, quote = FALSE)
    cat(sprintf("未匹配基因已保存: %s (%d行)\n", 
                file.path(output_dir, unmatched_file), nrow(unmatched_df)))
  }
  
  # 保存来源统计
  source_stats_file <- file.path(output_dir, "source_statistics.csv")
  write.csv(source_stats, file = source_stats_file, row.names = FALSE, quote = FALSE)
  cat(sprintf("来源统计已保存: %s\n", source_stats_file))
  
} else if (include_all_genes) {
  # 如果没有匹配，但用户想保存所有基因
  write.csv(deg_filtered, 
            file = file.path(output_dir, unmatched_file),
            row.names = FALSE, quote = FALSE)
  cat(sprintf("未找到匹配，差异基因已保存到: %s\n", 
              file.path(output_dir, unmatched_file)))
}

# 保存分析参数
params_file <- file.path(output_dir, "mapping_parameters.txt")
sink(params_file)
cat("差异基因-通路匹配分析参数\n")
cat("========================\n\n")
cat("分析时间:", date(), "\n")
cat("工作目录:", getwd(), "\n\n")
cat("输入文件:\n")
cat("  差异基因文件:", deg_file, "\n")
cat("  GMT文件数:", length(gmt_files), "\n\n")
cat("筛选参数:\n")
cat("  P值阈值:", pval_threshold, "\n")
cat("  FDR阈值:", fdr_threshold, "\n")
cat("  最小|logFC|:", min_logFC, "\n")
cat("  使用FDR筛选:", use_adjusted_pvalue, "\n")
cat("  移除重复项:", remove_duplicates, "\n\n")
cat("基因处理:\n")
cat("  基因ID类型:", gene_id_type, "\n")
cat("  标准化为大写:", standardize_gene_ids, "\n")
cat("  移除Ensembl版本:", remove_version, "\n\n")
cat("统计信息:\n")
cat("  原始差异基因数:", nrow(deg_selected), "\n")
cat("  筛选后差异基因数:", nrow(deg_filtered), "\n")
cat("  上调基因数:", up_genes, "\n")
cat("  下调基因数:", down_genes, "\n")
cat("  通路总数:", nrow(all_pathways), "\n")
if (nrow(gene_pathway_mapping) > 0) {
  cat("  匹配基因数:", length(unique(gene_pathway_mapping$GeneID)), "\n")
  cat("  匹配通路数:", length(unique(gene_pathway_mapping$PathwayID)), "\n")
  cat("  总映射数:", nrow(gene_pathway_mapping), "\n")
}
cat("  未匹配基因数:", length(unmatched_genes), "\n")

if (generate_venn) {
  cat("\nVenn图统计:\n")
  cat("  数据集1（", venn_name1, "）: ", length(genes_set1), " 个基因\n")
  cat("  数据集2（", venn_name2, "）: ", length(genes_set2), " 个基因\n")
  cat("  重合基因数: ", length(overlap_genes), "\n")
}
sink()
cat(sprintf("分析参数已保存: %s\n", params_file))

# ========================== 输出简要报告 ====================================
cat("\n========== 分析完成 ==========\n")
cat("输出文件保存在:", output_dir, "\n")
cat("\n生成的文件:\n")
if (nrow(gene_pathway_mapping) > 0) {
  cat(sprintf("  1. %s - 基因-通路映射表\n", output_file))
  cat(sprintf("  2. %s - 通路汇总统计\n", summary_file))
  cat(sprintf("  3. %s - 基因统计（每个基因出现的通路）\n", gene_stats_file))
  if (length(unmatched_genes) > 0) {
    cat(sprintf("  4. %s - 未匹配基因列表\n", unmatched_file))
  }
  cat(sprintf("  5. source_statistics.csv - 各GMT文件统计\n"))
} else {
  cat(sprintf("  1. %s - 未匹配基因列表\n", unmatched_file))
}
cat(sprintf("  6. mapping_parameters.txt - 分析参数\n"))

if (generate_venn) {
  cat(sprintf("  7. venn_diagram.pdf - Venn图\n"))
  cat(sprintf("  8. venn_diagram.png - Venn图(PNG)\n"))
  if (length(overlap_genes) > 0) {
    cat(sprintf("  9. venn_overlap_genes.csv - 重合基因列表\n"))
    cat(sprintf("  10. venn_unique_set1_genes.csv - 数据集1特有基因列表\n"))
    cat(sprintf("  11. venn_unique_set2_genes.csv - 数据集2特有基因列表\n"))
  }
}

# 显示前10个匹配的通路
if (nrow(gene_pathway_mapping) > 0) {
  cat("\n包含最多匹配基因的前10个通路:\n")
  top_pathways <- pathway_summary %>%
    head(10) %>%
    select(PathwayName, GeneCount, UpGenes, DownGenes)
  print(top_pathways)
  
  cat("\n匹配最多通路的前10个基因:\n")
  top_genes <- gene_stats %>%
    head(10) %>%
    select(GeneID, logFC, Regulation, PathwayCount)
  print(top_genes)
}

cat("\n=== 脚本运行完成 ===\n")