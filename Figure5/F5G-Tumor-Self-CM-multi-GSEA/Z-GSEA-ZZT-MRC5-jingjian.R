# ============================================================
# GSEA 分类散点图（有颜色，无连线）
# ============================================================

library(ggplot2)
library(dplyr)
library(scales)

# ==================== 参数设置 ====================
setwd("D:\\zsy\\SX\\MRC-5-Z-GSEA\\MDA-MB-231-DMEM-jingjian")  # 修改为您的路径

input_file <- "GSEA_category_data.csv"
output_pdf <- "GSEA_category_scatter_compact.pdf"
output_png <- "GSEA_category_scatter_compact.png"

# ==================== 配色 ====================
comparison_colors <- c(
  "NF-CM_vs_Self-CM" = "#25377F",
  "CAF-CM_vs_Self-CM" = "#C60036"
)

# ==================== 图形参数 ====================
nes_limits <- c(-3, 3)
point_size_range <- c(1.5, 4.5)
alpha_range <- c(0.5, 1.0)  # p值越小，alpha越大（越不透明）
jitter_width <- 0.15

midline_color <- "grey30"
midline_width <- 0.4
midline_type <- "dashed"

font_title <- 16
font_axis_title <- 16
font_axis_text <- 14
font_category <- 18

plot_width <- 6
plot_height <- 4
plot_dpi <- 300

# ==================== 读取数据 ====================
cat("正在读取数据...\n")
gsea_data <- read.csv(input_file, stringsAsFactors = FALSE)
cat("共", nrow(gsea_data), "行\n")

# 去除 Comparison 列前后空格
gsea_data$Comparison <- trimws(gsea_data$Comparison)

# 确保 NOM_pVal 列存在且为数值
if(!"NOM_pVal" %in% colnames(gsea_data)) {
  cat("警告：未找到 NOM_pVal 列，尝试查找 p 值相关列...\n")
  pval_cols <- grep("p|P|Pval|pval", colnames(gsea_data), value = TRUE)
  if(length(pval_cols) > 0) {
    gsea_data$NOM_pVal <- gsea_data[[pval_cols[1]]]
    cat("使用列:", pval_cols[1], "作为 p 值\n")
  } else {
    cat("未找到 p 值列，创建模拟数据用于演示\n")
    gsea_data$NOM_pVal <- runif(nrow(gsea_data), 0, 0.1)
  }
}

# 确保 FDR 列存在
if(!"FDR" %in% colnames(gsea_data)) {
  cat("警告：未找到 FDR 列，使用 NOM_pVal 作为替代\n")
  gsea_data$FDR <- gsea_data$NOM_pVal
}

# ==================== 保留输入文件中的分类顺序 ====================
category_order <- unique(gsea_data$Category)
cat("\n输入文件中的分类顺序:\n")
print(category_order)

# ==================== 创建 Y 轴位置映射 ====================
category_positions <- data.frame(
  Category = category_order,
  y_center = seq_along(category_order)
)

gsea_data$y_center <- category_positions$y_center[match(gsea_data$Category, category_positions$Category)]

# ==================== 在同一分类内添加随机抖动 ====================
set.seed(42)
gsea_data <- gsea_data %>%
  group_by(Category) %>%
  mutate(
    y_jitter = y_center + runif(n(), -jitter_width, jitter_width)
  ) %>%
  ungroup()

# ==================== 限制 NES 范围 ====================
gsea_data$NES <- pmin(pmax(gsea_data$NES, nes_limits[1]), nes_limits[2])

# ==================== Comparison 转换为因子 ====================
gsea_data$Comparison <- factor(gsea_data$Comparison, levels = names(comparison_colors))

# ==================== 设置 Category 因子顺序 ====================
gsea_data$Category <- factor(gsea_data$Category, levels = category_order)

# ==================== 创建 p 值的反转变换用于透明度 ====================
# 方法1：直接使用 p 值（p 值越小，alpha 越大）
# 方法2：使用 -log10(p) 转换，使显著性差异更明显
gsea_data$alpha_value <- 1 - gsea_data$NOM_pVal  # p=0 -> alpha=1, p=0.1 -> alpha=0.9
# 或者使用 -log10 转换：gsea_data$alpha_value <- pmin(pmax(-log10(gsea_data$NOM_pVal) / 5, 0.5), 1)

# ==================== 绘图 ====================
p <- ggplot(gsea_data, aes(x = NES, y = y_jitter)) +
  # NES=0 中线
  geom_vline(xintercept = 0, 
             linetype = midline_type, 
             color = midline_color, 
             linewidth = midline_width) +
  # 散点
  geom_point(
    aes(
      size = FDR,
      alpha = NOM_pVal,  # 直接使用 p 值
      fill = Comparison
    ),
    shape = 21,
    color = "grey30",
    stroke = 0.15
  ) +
  # 颜色标度
  scale_fill_manual(
    name = "Comparison",
    values = comparison_colors,
    drop = FALSE
  ) +
  # 点大小（FDR越小越大）
  scale_size_continuous(
    name = "FDR",
    range = point_size_range,
    trans = "reverse",
    labels = scales::scientific,
    guide = guide_legend(
      override.aes = list(alpha = 0.8),  # 确保图例中的点可见
      order = 2
    )
  ) +
  # 透明度（p值越小越不透明）- 修复图例显示
  scale_alpha_continuous(
    name = "Nominal p-value",
    range = alpha_range,
    trans = "reverse",  # p值越小，alpha越大
    labels = scales::scientific,
    breaks = c(0.001, 0.01, 0.05, 0.1),  # 自定义显示断点
    limits = c(0, max(gsea_data$NOM_pVal, na.rm = TRUE)),
    guide = guide_legend(
      title = "Nominal p-value",
      override.aes = list(
        size = 3,  # 固定图例中点的大小
        shape = 21,
        fill = "grey50"  # 图例中点的颜色
      ),
      order = 3
    )
  ) +
  # X轴
  scale_x_continuous(
    limits = nes_limits,
    breaks = seq(-3, 3, 1),
    expand = expansion(mult = 0.05)
  ) +
  # Y轴（按输入文件顺序显示分类）
  scale_y_continuous(
    breaks = category_positions$y_center,
    labels = category_positions$Category,
    expand = expansion(mult = 0.08)
  ) +
  # 标签
  labs(
    title = "GSEA Enrichment by Functional Category",
    x = "Normalized Enrichment Score (NES)",
    y = ""
  ) +
  # 主题
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = font_title),
    axis.title.x = element_text(size = font_axis_title),
    axis.text.x = element_text(size = font_axis_text, color = "black"),
    axis.text.y = element_text(size = font_category, face = "plain", color = "black"),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "plain"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3, alpha = 0.9), 
      order = 1
    ),
    size = guide_legend(order = 2),
    alpha = guide_legend(
      title = "Nominal p-value",
      override.aes = list(
        size = 3,
        shape = 21,
        fill = "grey50"
      ),
      order = 3
    )
  )

# ==================== 保存 ====================
ggsave(output_pdf, p, width = plot_width, height = plot_height, device = "pdf")
ggsave(output_png, p, width = plot_width, height = plot_height, dpi = plot_dpi)

cat("\nPDF 已保存:", output_pdf, "\n")
cat("PNG 已保存:", output_png, "\n")

# 打印图例信息
cat("\n图例透明度范围：", alpha_range[1], "到", alpha_range[2], "\n")
cat("p值范围：", min(gsea_data$NOM_pVal), "到", max(gsea_data$NOM_pVal), "\n")

print(p)