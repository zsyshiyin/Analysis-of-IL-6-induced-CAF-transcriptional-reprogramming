library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

# 设置工作目录
setwd("D:\\zsy\\SX\\Fomal SX\\ssGSEA\\IL6-final")

# 1. 读取所有.tsv通路文件
files <- grep("\\.tsv$", dir(), value = TRUE)
if(length(files) == 0) stop("未找到任何.tsv通路文件")

# 读取函数（带错误处理）
read_tsv_safe <- function(f) {
  tryCatch({
    df <- read.delim(f, check.names = FALSE)
    df$.filename <- f  # 保留原始文件名
    df
  }, error = function(e) {
    warning(paste("文件读取失败:", f, "错误:", e$message))
    NULL
  })
}

data <- lapply(files, read_tsv_safe)
data <- data[!sapply(data, is.null)] # 移除失败的文件

# 2. 合并数据并提取纯净通路名称
dataSet <- ldply(data, data.frame)
dataSet$pathway <- gsub("\\.tsv$", "", basename(dataSet$.filename))

# 3. 读取NES-p.txt统计文件 - 并获取通路顺序
nes_data <- NULL
if(file.exists("NES-p-v.txt")) {
  nes_data <- tryCatch({
    df <- read.delim("NES-p-v.txt", check.names = FALSE)
    # 规范列名（兼容不同大小写和空格）
    colnames(df) <- c("GeneSet", "ES", "NES", "Nominal_p", "FDR_q", "FWER_p")
    df
  }, error = function(e) {
    warning(paste("NES-p.txt读取失败:", e$message))
    NULL
  })
} else {
  stop("找不到NES-p-v.txt文件")
}

# 4. 关键修改：严格按照NES-p-v.txt文件的顺序排列通路
# 直接从NES-p-v.txt获取通路顺序
nes_pathway_order <- as.character(nes_data$GeneSet)
cat("\n=== NES-p-v.txt中的通路顺序 ===\n")
for(i in seq_along(nes_pathway_order)) {
  cat(sprintf("%2d: %s\n", i, nes_pathway_order[i]))
}

# 只保留在数据中存在的通路，并保持NES-p-v.txt的顺序
pathway_names <- nes_pathway_order[nes_pathway_order %in% unique(dataSet$pathway)]

# 检查是否有缺失的通路
missing_in_data <- nes_pathway_order[!nes_pathway_order %in% unique(dataSet$pathway)]
if(length(missing_in_data) > 0) {
  warning(paste("以下通路在.tsv文件中不存在:", 
                paste(missing_in_data, collapse = ", ")))
}

# 检查是否有额外的通路在数据中存在但不在NES-p-v.txt中
extra_in_data <- unique(dataSet$pathway)[!unique(dataSet$pathway) %in% nes_pathway_order]
if(length(extra_in_data) > 0) {
  warning(paste("以下通路在.tsv文件中存在但不在NES-p-v.txt中:", 
                paste(extra_in_data, collapse = ", ")))
}

cat("\n=== 最终显示顺序（按NES-p-v.txt顺序）===\n")
for(i in seq_along(pathway_names)) {
  cat(sprintf("%2d: %s\n", i, pathway_names[i]))
}

# 5. 精确匹配通路名称与统计值
match_stats <- function(pathway_names, nes_df) {
  result <- data.frame(
    Pathway = pathway_names,
    ES = NA, NES = NA, Nominal_p = NA, FDR_q = NA, FWER_p = NA,
    stringsAsFactors = FALSE
  )
  
  if(!is.null(nes_df)) {
    # 去除文件名和GeneSet中的干扰字符后匹配
    clean_name <- function(x) tolower(gsub("[^[:alnum:]]", "", x))
    
    for(i in seq_along(pathway_names)) {
      # 先尝试直接匹配
      match_idx <- which(nes_df$GeneSet == pathway_names[i])
      
      # 如果直接匹配失败，尝试模糊匹配
      if(length(match_idx) == 0) {
        clean_path <- clean_name(pathway_names[i])
        clean_nes <- clean_name(nes_df$GeneSet)
        match_idx <- which(clean_nes == clean_path)
      }
      
      if(length(match_idx) > 0) {
        result[i, c("ES", "NES", "Nominal_p", "FDR_q", "FWER_p")] <- 
          nes_df[match_idx[1], c("ES", "NES", "Nominal_p", "FDR_q", "FWER_p")]
      }
    }
  }
  return(result)
}

stats_data <- match_stats(pathway_names, nes_data)

# 6. 创建高区分度颜色 - 严格按照pathway_names的顺序
gseaCol <- colorRampPalette(c("#E4945A", "#CA0E12", "#25377F"))(length(pathway_names))
names(gseaCol) <- pathway_names  # 直接按NES-p-v.txt顺序命名

# 7. 绘制GSEA图 - 按照pathway_names的顺序绘制
pGsea <- ggplot(dataSet, aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES, 
                             color = factor(pathway, levels = pathway_names))) +
  geom_line(size = 1.2) +
  scale_color_manual(values = gseaCol) +
  labs(x = "", y = "Enrichment Score") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.margin = unit(c(1, 1, 0.5, 0.5), "line")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")

# 8. 绘制基因位置图 - 按照pathway_names的顺序绘制
pGene <- ggplot(dataSet, aes(x = RANK.IN.GENE.LIST, 
                             y = factor(pathway, levels = pathway_names), 
                             fill = factor(pathway, levels = pathway_names))) +
  geom_tile(height = 0.7) +
  scale_fill_manual(values = gseaCol) +
  labs(x = "Gene Rank (high → low expression)", y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 10, vjust = -0.5),
    legend.position = "none",
    plot.margin = unit(c(0.5, 1, 1, 0.5), "line")
  )

# 9. 创建图例表格函数
create_legend_table <- function(df, colors) {
  # 添加颜色列
  df <- cbind(Color = "", df)
  
  # 创建基本表格
  tab <- tableGrob(df, 
                   theme = ttheme_minimal(
                     core = list(
                       fg_params = list(hjust = 0, x = 0.05, fontsize = 8),
                       padding = unit(c(2, 2), "mm")
                     ),
                     colhead = list(
                       fg_params = list(hjust = 0, x = 0.05, fontsize = 8, fontface = "bold"),
                       padding = unit(c(2, 2), "mm")
                     )
                   ),
                   rows = NULL)
  
  # 添加颜色块
  for (i in seq_along(colors)) {
    tab <- gtable::gtable_add_grob(tab,
                                   rectGrob(gp = gpar(fill = colors[i], col = "black")),
                                   t = i + 1,
                                   l = 1,
                                   b = i + 1,
                                   r = 1)
  }
  
  # 设置Pathway列宽度
  col_widths <- tab$widths
  col_widths[2] <- unit(13.5, "cm")
  tab$widths <- col_widths
  
  return(tab)
}

# 10. 准备图例数据 - 严格按照pathway_names的顺序
legend_df <- data.frame(
  Pathway = pathway_names,
  stringsAsFactors = FALSE
)

if(!all(is.na(stats_data$ES))) legend_df$ES <- format(stats_data$ES, digits = 3)
if(!all(is.na(stats_data$NES))) legend_df$NES <- format(stats_data$NES, digits = 3)
if(!all(is.na(stats_data$Nominal_p))) legend_df$`Nominal p` <- format(stats_data$Nominal_p, digits = 3)
if(!all(is.na(stats_data$FDR_q))) legend_df$`FDR q` <- format(stats_data$FDR_q, digits = 3)
if(!all(is.na(stats_data$FWER_p))) legend_df$`FWER p` <- format(stats_data$FWER_p, digits = 3)

# 11. 创建图例表格
legend_table <- create_legend_table(legend_df, gseaCol)

# 12. 组合图形
gGsea <- ggplot_gtable(ggplot_build(pGsea))
gGene <- ggplot_gtable(ggplot_build(pGene))
maxWidth <- grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths <- as.list(maxWidth)
gGene$widths <- as.list(maxWidth)

# 13. 调整布局
final_plot <- grid.arrange(
  arrangeGrob(gGsea, gGene, nrow = 2, heights = c(0.8, 0.2)),
  legend_table,
  ncol = 2,
  widths = c(0.3, 0.3),
  padding = unit(0.2, "cm")
)

# 14. 输出PDF
pdf('multipleGSEA.pdf', width = 18, height = 6)
grid.draw(final_plot)
dev.off()

# 15. 打印匹配结果
cat("\n=== 通路统计值匹配结果（按NES-p-v.txt顺序）===\n")
print(data.frame(Pathway = pathway_names, 
                 NES = stats_data$NES,
                 `Nominal p` = stats_data$Nominal_p))