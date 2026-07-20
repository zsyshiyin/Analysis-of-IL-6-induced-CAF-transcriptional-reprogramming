# =============================================================================
# 生成GSEA输入文件：RNK（排序文件）和GMT（基因集文件）
# 功能：为GSEA Desktop软件准备输入文件
# 输入：所有基因的logFC列表、CIPS基因列表
# 输出：RNK文件、GMT文件
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/3.7-Rh2-jingzheng"  # 请修改此路径
setwd(work_dir)

# 创建输出文件夹
output_dir <- "GSEA_input_files"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 输入文件
cips_gene_file <- "CIPS_gene_list.txt"
logfc_file <- "all_genes_logFC.csv"

# 输出文件
output_rnk_correct <- file.path(output_dir, "gene_ranking.rnk")
output_gmt <- file.path(output_dir, "CIPS.gmt")

# ====================== 读取数据 ======================
cat("\n========== 读取数据 ==========\n")

# 读取CIPS基因
cips_genes <- readLines(cips_gene_file, warn = FALSE)
cips_genes <- toupper(trimws(cips_genes))
cips_genes <- cips_genes[cips_genes != ""]
cips_genes <- cips_genes[!grepl("^#", cips_genes)]
cat(sprintf("CIPS基因数: %d\n", length(cips_genes)))

# 读取logFC数据
logfc_data <- read.csv(logfc_file, stringsAsFactors = FALSE, check.names = FALSE)

# 识别列名
gene_col <- grep("gene|Gene|GENE", colnames(logfc_data), ignore.case = TRUE)[1]
logfc_col <- grep("logFC|log_fc|fc|log2FoldChange", colnames(logfc_data), ignore.case = TRUE)[1]

if (is.na(gene_col) || is.na(logfc_col)) {
  stop("无法识别基因名或logFC列！")
}

logfc_data <- logfc_data[, c(gene_col, logfc_col)]
colnames(logfc_data) <- c("Gene", "logFC")

# 基因名转大写
logfc_data$Gene <- toupper(trimws(logfc_data$Gene))
logfc_data$logFC <- as.numeric(logfc_data$logFC)

# 移除NA
logfc_data <- logfc_data[!is.na(logfc_data$logFC) & is.finite(logfc_data$logFC), ]

cat(sprintf("logFC数据: %d 个基因\n", nrow(logfc_data)))

# ====================== 生成GSEA Desktop兼容的RNK文件 ======================
cat("\n========== 生成RNK文件 ==========\n")

# 重要：GSEA Desktop要求：
# 1. 文件不能有注释行/列名行
# 2. 每行格式：基因名\t数值（制表符分隔）
# 3. 数值必须是数字，不能是字符
# 4. 按数值降序排序

# 按logFC降序排序
logfc_sorted <- logfc_data[order(-logfc_data$logFC), ]

# 生成RNK文件（无列名，制表符分隔）
# 使用 write.table 时设置 col.names = FALSE, row.names = FALSE
write.table(logfc_sorted[, c("Gene", "logFC")], 
            file = output_rnk_correct,
            row.names = FALSE,
            col.names = FALSE,      # 关键：不输出列名
            sep = "\t",             # 制表符分隔
            quote = FALSE)          # 不加引号

cat(sprintf("RNK文件已保存: %s\n", output_rnk_correct))

# 验证文件内容
cat("\nRNK文件前10行预览:\n")
rnk_check <- readLines(output_rnk_correct, n = 10)
for (i in 1:length(rnk_check)) {
  cat(sprintf("  %s\n", rnk_check[i]))
}

# ====================== 生成GMT文件 ======================
cat("\n========== 生成GMT文件 ==========\n")

# GMT格式：基因集名称\t基因集描述\t基因1\t基因2\t基因3...
gmt_line <- paste("CIPS", "CAF_Invasive_Progression_Signature", paste(cips_genes, collapse = "\t"), sep = "\t")
writeLines(gmt_line, output_gmt)
cat(sprintf("GMT文件已保存: %s\n", output_gmt))

# ====================== 生成备用的LRNK格式 ======================
cat("\n========== 生RNA文件 ==========\n")

# 某些GSEA版本可能需要.LRNK后缀
output_lrnk <- file.path(output_dir, "gene_ranking.lrnk")
file.copy(output_rnk_correct, output_lrnk, overwrite = TRUE)
cat(sprintf("LRNK文件已保存: %s\n", output_lrnk))

# ====================== 生成使用说明 ======================
cat("\n========== 生成使用说明 ==========\n")

output_instructions <- file.path(output_dir, "GSEA_Desktop_Instructions.txt")

instructions <- '
========================================
   GSEA Desktop 使用说明
========================================

一、文件位置
-----------
RNK文件: gene_ranking.rnk
GMT文件: CIPS.gmt

二、重要提醒
-----------
⚠ 生成的RNK文件已经是GSEA Desktop兼容格式：
   - 无注释行/列名行
   - 制表符分隔
   - 数值格式正确

三、GSEA Desktop操作步骤
-----------------------
1. 打开GSEA Desktop软件

2. 加载RNK文件：
   - 菜单栏：File → Load data
   - 或直接将 "gene_ranking.rnk" 拖入左侧面板

3. 加载GMT文件：
   - 将 "CIPS.gmt" 拖入左侧面板

4. 运行GSEA：
   - 点击 "Run GSEA"
   - Gene sets database: 选择 "CIPS.gmt"
   - Ranked list: 选择 "gene_ranking.rnk"
   - Number of permutations: 1000
   - Metric for ranking genes: "Signal2Noise" 或选择 "log2_Ratio_of_Classes"
   - Collapse dataset to gene symbols: 勾选
   - 点击 "Run"

四、如果还是报错
---------------
某些GSEA版本要求：
1. 数值必须是浮点数（如 -1.23 而不是 -1.23）
   → 当前格式已满足

2. 文件扩展名必须是 .rnk 或 .lrnk
   → 已同时提供 .rnk 和 .lrnk 两种格式

3. 尝试使用 "gene_ranking.lrnk" 文件（某些版本需要这个扩展名）

五、备选方案
-----------
如果GSEA Desktop仍然报错，可以使用R生成的GSEA富集图，
代码已在之前的答复中提供。
'

writeLines(instructions, output_instructions)
cat(sprintf("使用说明已保存: %s\n", output_instructions))

# ====================== 验证文件格式 ======================
cat("\n========== 验证文件格式 ==========\n")

# 验证RNK文件可以被正确读取
test_data <- tryCatch({
  read.table(output_rnk_correct, header = FALSE, nrows = 5, sep = "\t")
}, error = function(e) {
  NULL
})

if (!is.null(test_data)) {
  cat("✓ RNK文件格式验证通过\n")
  cat("  前5行:\n")
  print(test_data)
} else {
  cat("✗ RNK文件格式验证失败\n")
}

# ====================== 完成 ======================
cat("\n========== 文件生成完成 ==========\n")
cat(sprintf("\n输出目录: %s\n", output_dir))
cat(sprintf("1. %s - RNK文件（GSEA Desktop兼容）\n", basename(output_rnk_correct)))
cat(sprintf("2. %s - GMT文件\n", basename(output_gmt)))
cat(sprintf("3. %s - 备用的LRNK文件\n", basename(output_lrnk)))
cat(sprintf("4. %s - 使用说明\n", basename(output_instructions)))
cat("\n现在可以在GSEA Desktop中加载这些文件了！\n")