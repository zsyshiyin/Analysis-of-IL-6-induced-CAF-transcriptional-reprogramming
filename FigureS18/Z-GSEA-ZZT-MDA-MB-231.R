# GSEA比较分析柱形图 - 横向渐变柱条版本
# 功能：绘制横向柱形图比较GSEA分析的NES值
# 柱条颜色为左右渐变，粗细表示q值大小，透明度表示p值

# 加载必要的包
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(gtable)) install.packages("gtable")
if (!require(stringr)) install.packages("stringr")  # 新增：用于文本换行处理
library(ggplot2)
library(dplyr)
library(gridExtra)
library(gtable)
library(stringr)  # 新增：用于文本换行处理

# ==================== 参数设置区域 ====================
# 设置工作路径 - 所有文件的读取和保存都将基于此路径
setwd("D:\\zsy\\SX\\MRC-5-Z-GSEA\\MCF-7-qinxiqianyi")  # 请替换为您的实际工作路径

# 设置您的数据文件路径（可以使用相对路径或绝对路径）
# 如果数据文件在工作目录下，可以直接使用文件名
data_file_path <- "GSEAzzt-MCF7.csv"  # 修改为相对路径

# 渐变颜色设置 - 根据NES正负值和比较组设置不同的渐变配色
gradient_colors <- list(
  # CAF-CM_vs_NF-CM 的配色
  "CAF-CM_vs_NF-CM" = list(
    positive = c(left = "#25377F", right = "#CA0E12"),  # 红色渐变（正值）
    negative = c(left = "#CA0E12", right = "#25377F")   # 蓝色渐变（负值）
  ),
  # CAF(+G-Rh2)-CM_vs_CAF-CM 的配色
  "CAF(+G-Rh2)-CM_vs_CAF-CM" = list(
    positive = c(left = "#CA0E12", right = "#E4945A"),  # 青色渐变（正值）
    negative = c(left = "#E4945A", right = "#CA0E12")   # 橙色渐变（负值）
  )
)

# NES范围设置
nes_limits <- c(-2.5, 2.5)

# 字体大小设置
font_sizes <- list(
  axis_title = 18,
  axis_text = 14,
  strip_text = 20,
  pathway_text = 25  # 通路名称增大2个字号（从20增大到24）
)

# 图形尺寸设置
plot_width <- 10
plot_height <- 24
plot_dpi <- 300

x_axis_physical_width <- 6.8  # 直接设置横坐标轴的物理宽度，增大此值使横坐标轴变宽
y_label_width <- 3.5        # 左侧Y轴标签区域固定宽度（一般不需要修改

# 输出文件设置（使用相对路径，将保存在工作目录下）
main_output_file <- "GSEA_comparison_gradient_bar_plot_main.png"
main_output_pdf <- "GSEA_comparison_gradient_bar_plot_main.pdf"
legend_output_file <- "GSEA_comparison_gradient_bar_plot_legends.png"
legend_output_pdf <- "GSEA_comparison_gradient_bar_plot_legends.pdf"

# ==================== 文本换行处理函数 ====================

#' 将长文本按单词数量换行，每行不超过3个单词
#' @param text 原始文本
#' @param max_words 每行最大单词数，默认为3
#' @return 换行后的文本
wrap_pathway_text <- function(text, max_words = 3) {
  # 如果文本为空，直接返回
  if (is.na(text) || text == "") return(text)
  
  # 按空格分割单词
  words <- unlist(strsplit(as.character(text), "\\s+"))
  
  # 如果单词数量小于等于最大单词数，直接返回原文本
  if (length(words) <= max_words) {
    return(text)
  }
  
  # 将单词分组，每组最多max_words个单词
  n_groups <- ceiling(length(words) / max_words)
  wrapped_lines <- character(n_groups)
  
  for (i in 1:n_groups) {
    start_idx <- (i - 1) * max_words + 1
    end_idx <- min(i * max_words, length(words))
    wrapped_lines[i] <- paste(words[start_idx:end_idx], collapse = " ")
  }
  
  # 用换行符连接各组
  return(paste(wrapped_lines, collapse = "\n"))
}

#' 批量处理通路名称文本换行
#' @param df 包含Pathway列的数据框
#' @return 处理后的数据框
process_pathway_text_wrap <- function(df) {
  df <- df %>%
    mutate(
      # 为每条通路创建换行后的文本
      Pathway_wrapped = sapply(Pathway, wrap_pathway_text, max_words = 3)
    )
  return(df)
}

# ==================== 数据处理函数 ====================

#' 准备GSEA数据用于绘图
prepare_gsea_data <- function(df) {
  # 首先处理通路名称换行
  df <- process_pathway_text_wrap(df)
  
  # 限制NES范围
  df$NES <- pmin(pmax(df$NES, nes_limits[1]), nes_limits[2])
  
  # 创建柱条宽度映射 - q值越小，柱条越粗
  # 调整宽度映射，确保q值<0.8都能显示明显的条带，宽度范围0.1-0.75
  df$bar_width <- 0.1 + (1 - pmin(df$FDR, 0.8)/0.8) * 0.65  # 宽度范围：0.1-0.75
  
  # 创建柱条透明度映射 - p值越小，柱条越不透明
  # 修正透明度映射：p值=0时alpha=1.0（完全不透明），p值=0.7时alpha=0.1
  df$bar_alpha <- 1.0 - (pmin(df$NOM_pVal, 0.7)/0.7) * 0.9  # 透明度范围：0.1-1.0
  
  # 为渐变创建位置信息
  df$group_id <- as.numeric(factor(df$Pathway)) * 100 + as.numeric(factor(df$Comparison))
  
  # 检查透明度映射
  cat("透明度映射检查:\n")
  cat("p值范围:", range(df$NOM_pVal), "\n")
  cat("alpha值范围:", range(df$bar_alpha), "\n")
  cat("p=0时的alpha:", df$bar_alpha[df$NOM_pVal == min(df$NOM_pVal)][1], "\n")
  cat("p=0.2时的alpha:", 1.0 - (0.2/0.2) * 0.6, "\n")
  
  return(df)
}

# ==================== 渐变柱条绘制函数 ====================

#' 创建渐变柱条数据
create_gradient_data <- function(df) {
  gradient_data <- data.frame()
  
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    comparison <- row$Comparison
    nes_direction <- ifelse(row$NES >= 0, "positive", "negative")
    
    # 获取该比较组和NES方向的渐变颜色
    colors <- gradient_colors[[comparison]][[nes_direction]]
    
    # 创建渐变数据点
    n_points <- 20  # 渐变点数
    for (j in 1:n_points) {
      # 计算渐变位置和颜色
      progress <- (j - 1) / (n_points - 1)
      
      # 线性插值计算颜色
      r1 <- col2rgb(colors["left"])[1]
      g1 <- col2rgb(colors["left"])[2]
      b1 <- col2rgb(colors["left"])[3]
      
      r2 <- col2rgb(colors["right"])[1]
      g2 <- col2rgb(colors["right"])[2]
      b2 <- col2rgb(colors["right"])[3]
      
      r <- round(r1 + (r2 - r1) * progress)
      g <- round(g1 + (g2 - g1) * progress)
      b <- round(b1 + (b2 - b1) * progress)
      
      color <- rgb(r, g, b, maxColorValue = 255)
      
      # 计算柱条位置
      if (row$NES >= 0) {
        x_start <- (j - 1) * (row$NES / n_points)
        x_end <- j * (row$NES / n_points)
      } else {
        x_start <- row$NES + (j - 1) * (-row$NES / n_points)
        x_end <- row$NES + j * (-row$NES / n_points)
      }
      
      gradient_data <- rbind(gradient_data, data.frame(
        Pathway = row$Pathway,
        Pathway_wrapped = row$Pathway_wrapped,  # 添加换行后的通路名称
        Category = row$Category,
        Comparison = row$Comparison,
        NES = row$NES,
        FDR = row$FDR,
        NOM_pVal = row$NOM_pVal,
        bar_width = row$bar_width,
        bar_alpha = row$bar_alpha,
        group_id = row$group_id,
        x_start = x_start,
        x_end = x_end,
        color = color,
        gradient_pos = j,
        nes_direction = nes_direction
      ))
    }
  }
  
  return(gradient_data)
}

# ==================== 修改后的绘图函数 ====================

#' 绘制GSEA比较柱形图（无图例）
#' 修改：确保Category标签按照输入文件中的顺序排列
plot_gsea_comparison_no_legend <- function(df) {
  
  # 准备数据
  df_prepared <- prepare_gsea_data(df)
  
  # 创建渐变数据
  gradient_data <- create_gradient_data(df_prepared)
  
  # 修正：确保通路顺序正确，避免重复
  # 按Category和NES排序，使用原始Pathway进行排序
  pathway_order <- df_prepared %>%
    arrange(Category, NES) %>%
    distinct(Pathway, .keep_all = TRUE) %>%
    pull(Pathway)
  
  # 创建通路到数值的映射（使用原始Pathway作为标识）
  pathway_levels <- unique(pathway_order)
  pathway_to_num <- setNames(seq_along(pathway_levels), pathway_levels)
  
  # 为渐变数据添加数值化的y坐标
  # 注意：这里使用原始Pathway进行匹配
  gradient_data$pathway_num <- pathway_to_num[gradient_data$Pathway]
  
  # 为y轴标签创建数据框，使用换行后的通路名称
  y_labels_df <- df_prepared %>%
    distinct(Pathway, .keep_all = TRUE) %>%
    arrange(match(Pathway, pathway_levels)) %>%
    select(Pathway, Pathway_wrapped)
  
  y_labels <- setNames(y_labels_df$Pathway_wrapped, y_labels_df$Pathway)
  
  # ========== 修改点：获取Category的原始顺序 ==========
  # 获取Category在输入文件中出现的原始顺序
  category_order <- unique(df_prepared$Category)
  
  # 将Category转换为因子，并按照原始顺序设置水平
  df_prepared$Category <- factor(df_prepared$Category, levels = category_order)
  gradient_data$Category <- factor(gradient_data$Category, levels = category_order)
  
  # 创建图形 - 完全移除所有图例
  p <- ggplot() +
    # 添加NES=0参考线
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    # 添加渐变柱条 - 使用geom_rect
    geom_rect(
      data = gradient_data,
      aes(
        xmin = x_start,
        xmax = x_end,
        ymin = pathway_num - bar_width/2,
        ymax = pathway_num + bar_width/2,
        fill = color,
        alpha = bar_alpha
      )
    ) +
    # 颜色标度 - 使用identity
    scale_fill_identity() +
    # 透明度标度 - 不显示图例
    scale_alpha_continuous(guide = "none", range = c(0.3, 1.0)) +
    # 坐标轴设置
    scale_x_continuous(
      limits = nes_limits, 
      breaks = seq(nes_limits[1], nes_limits[2], 1.0)
    ) +
    scale_y_continuous(
      breaks = 1:length(pathway_levels),
      labels = y_labels[pathway_levels]  # 使用换行后的通路名称
    ) +
    # 主题设置 - 白色背景和表格线，完全移除图例
    theme_bw() +
    theme(
      # 调整y轴文本：增大字体，左对齐，支持多行文本
      axis.text.y = element_text(
        size = font_sizes$pathway_text, 
        color = "black",
        hjust = 1,  # 右对齐
        vjust = 0.5,  # 垂直居中
        lineheight = 0.9,  # 调整行间距
        margin = margin(r = 5, l = 5)  # 增加右边距
      ),
      axis.text.x = element_text(size = font_sizes$axis_text, color = "black"),
      axis.title.x = element_text(size = font_sizes$axis_title, face = "bold", 
                                  margin = margin(t = 15)),
      axis.title.y = element_blank(),
      legend.position = "none",  # 完全移除所有图例
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      strip.text = element_text(size = font_sizes$strip_text, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      # 增加整体边距，特别是左边距以容纳换行后的长文本
      plot.margin = margin(1, 1, 1, 2, "cm")
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)"
    ) +
    # 按分类分面，并按照原始顺序排列分类
    facet_grid(Category ~ ., scales = "free_y", space = "free_y")
  
  return(p)
}

# ==================== 修改后的图例函数 ====================

#' 创建比较组图例 - 渐变版本
#' 修改为使用gradient_colors中实际定义的比较组名称
create_comparison_legend <- function() {
  # 创建渐变数据
  gradient_legend_data <- data.frame()
  
  # 使用gradient_colors中实际定义的比较组名称
  comparisons <- names(gradient_colors)
  
  for (comp_idx in 1:length(comparisons)) {
    comparison <- comparisons[comp_idx]
    colors <- gradient_colors[[comparison]]
    
    # 创建10个渐变片段 - 保持10个片段以保证渐变程度不变
    for (j in 1:10) {
      progress <- (j - 1) / 9
      
      # 计算渐变颜色
      r1 <- col2rgb(colors["left"])[1]
      g1 <- col2rgb(colors["left"])[2]
      b1 <- col2rgb(colors["left"])[3]
      
      r2 <- col2rgb(colors["right"])[1]
      g2 <- col2rgb(colors["right"])[2]
      b2 <- col2rgb(colors["right"])[3]
      
      r <- round(r1 + (r2 - r1) * progress)
      g <- round(g1 + (g2 - g1) * progress)
      b <- round(b1 + (b2 - b1) * progress)
      
      color <- rgb(r, g, b, maxColorValue = 255)
      
      gradient_legend_data <- rbind(gradient_legend_data, data.frame(
        Comparison = comparison,
        xmin = 0.8 + (j - 1) * 0.01,
        xmax = 0.8 + j * 0.01,
        ymin = 2.2 - comp_idx * 0.6,
        ymax = 2.4 - comp_idx * 0.6,
        color = color
      ))
    }
  }
  
  # 创建标签数据
  label_data <- data.frame(
    Comparison = comparisons,
    x = 0.92,
    y = seq(1.76, by = -0.6, length.out = length(comparisons))
  )
  
  p <- ggplot() +
    # 添加渐变矩形
    geom_rect(
      data = gradient_legend_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color)
    ) +
    scale_fill_identity() +
    # 添加标签 - 位于横线正右侧
    geom_text(
      data = label_data,
      aes(x = x, y = y, label = Comparison),
      hjust = 0, size = 3.5, color = "black", fontface = "bold"
    ) +
    xlim(0.5, 2.0) +
    ylim(0.5, 2.5) +
    labs(title = "比较组") +
    theme_void() +
    theme(
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 5))
    )
  
  return(p)
}

#' 创建FDR图例 - 渐变粗细竖线
create_fdr_legend <- function() {
  # 创建渐变粗细数据
  fdr_data <- data.frame()
  n_segments <- 20
  
  for (i in 1:n_segments) {
    progress <- (i - 1) / (n_segments - 1)
    # 从细到粗的宽度变化
    width <- 0.01 + progress * 0.05  # 宽度范围：0.05-0.35
    
    fdr_data <- rbind(fdr_data, data.frame(
      segment = i,
      xmin = 0.97 - width/2,  # 中心对齐，宽度渐变
      xmax = 0.97 + width/2,
      ymin = 0.8 + (i - 1) * (2.2 / n_segments),
      ymax = 0.8 + i * (2.2 / n_segments)
    ))
  }
  
  p <- ggplot() +
    # 添加渐变粗细竖线
    geom_rect(
      data = fdr_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "black"
    ) +
    # 添加标签
    annotate("text", x = 1.01, y = 3.0, label = "0", size = 3, hjust = 0) +
    annotate("text", x = 1.01, y = 0.8, label = "0.25", size = 3, hjust = 0) +
    xlim(0.5, 1.8) +
    ylim(0.5, 3.2) +
    labs(title = "FDR q-val") +
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 3))
    )
  
  return(p)
}

#' 创建p-value图例 - 渐变透明度竖线
create_pval_legend <- function() {
  # 创建渐变透明度数据
  pval_data <- data.frame()
  n_segments <- 20
  
  for (i in 1:n_segments) {
    progress <- (i - 1) / (n_segments - 1)
    alpha <- 0.3 + progress * 0.7  # 从浅到深
    
    pval_data <- rbind(pval_data, data.frame(
      segment = i,
      xmin = 0.95,
      xmax = 1.00,
      ymin = 0.8 + (i - 1) * (2.2 / n_segments),
      ymax = 0.8 + i * (2.2 / n_segments),
      alpha = alpha
    ))
  }
  
  p <- ggplot() +
    # 添加渐变透明度竖线
    geom_rect(
      data = pval_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha),
      fill = "black"
    ) +
    scale_alpha_identity() +
    # 添加标签
    annotate("text", x = 1.01, y = 3.0, label = "0", size = 3, hjust = 0) +
    annotate("text", x = 1.01, y = 0.8, label = ">0.05", size = 3, hjust = 0) +
    xlim(0.5, 1.8) +
    ylim(0.5, 3.2) +
    labs(title = "NOM p-val") +
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 3))
    )
  
  return(p)
}

# ==================== 主程序 ====================

# 显示当前工作路径
cat("当前工作路径:", getwd(), "\n")

# 加载数据
cat("正在加载数据...\n")
cat("数据文件路径:", data_file_path, "\n")
tryCatch({
  gsea_data <- read.csv(data_file_path, stringsAsFactors = FALSE)
  cat("数据加载成功!\n")
}, error = function(e) {
  stop("数据加载失败: ", e$message, "\n请检查文件路径是否正确")
})

# 检查数据格式
required_columns <- c("Pathway", "Category", "Comparison", "NES", "FDR", "NOM_pVal")
missing_columns <- setdiff(required_columns, colnames(gsea_data))
if (length(missing_columns) > 0) {
  stop("数据缺少必要的列: ", paste(missing_columns, collapse = ", "))
}

# 检查数据
cat("数据预览:\n")
print(head(gsea_data))
cat("\n数据摘要:\n")
cat("通路数量:", length(unique(gsea_data$Pathway)), "\n")
cat("分类数量:", length(unique(gsea_data$Category)), "\n")
cat("比较组:", paste(unique(gsea_data$Comparison), collapse = ", "), "\n")
cat("NES范围:", paste(round(range(gsea_data$NES), 2), collapse = " 到 "), "\n")
cat("FDR范围:", paste(round(range(gsea_data$FDR), 4), collapse = " 到 "), "\n")
cat("NOM p-value范围:", paste(round(range(gsea_data$NOM_pVal), 4), collapse = " 到 "), "\n")

# 预览处理后的通路名称（示例）
cat("\n通路名称换行处理示例（前5条）:\n")
processed_sample <- process_pathway_text_wrap(gsea_data)
for (i in 1:min(5, nrow(processed_sample))) {
  cat("原始:", processed_sample$Pathway[i], "\n")
  cat("换行后:\n", processed_sample$Pathway_wrapped[i], "\n\n")
}

# 绘制主图形（无图例）
cat("正在生成主图形...\n")
main_plot <- plot_gsea_comparison_no_legend(gsea_data)

# 显示主图形
print(main_plot)

# 保存主图形 - 横坐标轴绝对宽度控制
main_output_file <- "GSEA_comparison_gradient_bar_plot_main.png"
cat("保存主图形到:", main_output_file, "\n")

# ========== 修改这里的ggsave ==========
total_plot_width <- x_axis_physical_width + y_label_width
ggsave(main_output_file, main_plot, 
       width = total_plot_width, 
       height = plot_height, 
       dpi = plot_dpi)

# 同时保存PDF版本
ggsave("GSEA_comparison_gradient_bar_plot_main.pdf", 
       main_plot, 
       width = total_plot_width, 
       height = plot_height)

cat("横坐标轴物理宽度设置为:", x_axis_physical_width, "英寸\n")
cat("整体图片宽度:", total_plot_width, "英寸\n")

# 创建图例
cat("正在生成图例...\n")
comparison_legend <- create_comparison_legend()
fdr_legend <- create_fdr_legend()
pval_legend <- create_pval_legend()

# 组合所有图例
all_legends <- grid.arrange(
  comparison_legend,
  fdr_legend,
  pval_legend,
  ncol = 1,
  heights = c(0.3, 0.35, 0.35)
)

# 保存组合图例
cat("保存图例到:", file.path(getwd(), legend_output_file), "\n")
ggsave(legend_output_file, all_legends, width = 4, height = 6, dpi = plot_dpi)

# 同时保存PDF版本
ggsave(legend_output_pdf, all_legends, width = 4, height = 6)
cat("保存图例PDF到:", file.path(getwd(), legend_output_pdf), "\n")

cat("\n=== 图形生成完成! ===\n")
cat("所有输出文件已保存到工作目录:\n", getwd(), "\n")
cat("主图:", main_output_file, "\n")
cat("主图PDF:", main_output_pdf, "\n")
cat("图例:", legend_output_file, "\n")
cat("图例PDF:", legend_output_pdf, "\n")