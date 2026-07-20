# ============================================================
# 空间转录组分析：通用CAF评分 vs CIPS评分 vs IL6+ACTA2+评分
# 功能：计算三种评分，并在同一区域进行比较
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(tidyr)

# ============================================================================
# ★★★ 全部可视化参数自定义区域 ★★★
# ============================================================================

# ---------- 路径设置 ----------
work_dir <- "D:/zsy/SX/GEO-mIHC"
output_dir <- file.path(work_dir, "figures_CAF_comparison-final")

# ========================================================================
# 图1：通用CAF评分空间分布
# ========================================================================
fig1_width <- 7
fig1_height <- 7
fig1_pt_size <- 1.5
fig1_alpha_min <- 0.1
fig1_alpha_max <- 1.0
fig1_title <- "General CAF Signature"
fig1_subtitle <- "Spatial Distribution"
fig1_title_size <- 14
fig1_title_face <- "bold"
fig1_subtitle_size <- 11
fig1_legend_title <- "CAF\nScore"
fig1_legend_text_size <- 10
fig1_legend_title_size <- 10
fig1_legend_height <- 0.6
fig1_color_low <- "grey90"
fig1_color_mid <- "lightgreen"
fig1_color_high <- "darkgreen"

# ========================================================================
# 图2：CIPS vs 通用CAF并排
# ========================================================================
fig2_width <- 14
fig2_height <- 7
fig2_pt_size <- 1.5
fig2_alpha_min <- 0.1
fig2_alpha_max <- 1.0
fig2_title_left <- "CIPS Signature"
fig2_title_right <- "General CAF Signature"
fig2_title_size <- 14
fig2_title_face <- "bold"
fig2_subtitle_size <- 11
fig2_legend_title_left <- "CIPS\nScore"
fig2_legend_title_right <- "CAF\nScore"
fig2_legend_text_size <- 10
fig2_legend_title_size <- 10
fig2_legend_height <- 0.5
fig2_cips_low <- "#25377F"
fig2_cips_mid <- "#2AA7DE"
fig2_cips_high <- "#C60036"
fig2_caf_low <- "grey90"
fig2_caf_mid <- "lightgreen"
fig2_caf_high <- "darkgreen"
fig2_tag_levels <- "A"

# ========================================================================
# 图3：IL6+ACTA2 vs 通用CAF并排
# ========================================================================
fig3_width <- 14
fig3_height <- 7
fig3_pt_size <- 1.5
fig3_alpha_min <- 0.1
fig3_alpha_max <- 1.0
fig3_title_left <- "IL6+ACTA2+ Signature"
fig3_title_right <- "General CAF Signature"
fig3_title_size <- 14
fig3_title_face <- "bold"
fig3_subtitle_size <- 11
fig3_legend_title_left <- "IL6+ACTA2\nScore"
fig3_legend_title_right <- "CAF\nScore"
fig3_legend_text_size <- 10
fig3_legend_title_size <- 10
fig3_legend_height <- 0.5
fig3_il6_low <- "grey90"
fig3_il6_mid <- "orange"
fig3_il6_high <- "#E64B35"
fig3_caf_low <- "grey90"
fig3_caf_mid <- "lightgreen"
fig3_caf_high <- "darkgreen"
fig3_tag_levels <- "A"

# ========================================================================
# 图4：三图并排（CIPS + IL6+ACTA2 + 通用CAF）
# ========================================================================
fig4_width <- 19
fig4_height <- 7
fig4_pt_size <- 1.2
fig4_alpha_min <- 0.1
fig4_alpha_max <- 1.0
fig4_title_left <- "CIPS"
fig4_title_mid <- "IL6+ACTA2+"
fig4_title_right <- "General CAF"
fig4_title_size <- 14
fig4_title_face <- "bold"
fig4_subtitle_size <- 11
fig4_legend_height <- 0.4
fig4_legend_text_size <- 10
fig4_legend_title_size <- 10
fig4_cips_low <- "grey90"
fig4_cips_mid <- "lightblue"
fig4_cips_high <- "darkred"
fig4_il6_low <- "grey90"
fig4_il6_mid <- "orange"
fig4_il6_high <- "#E64B35"
fig4_caf_low <- "grey90"
fig4_caf_mid <- "lightgreen"
fig4_caf_high <- "darkgreen"
fig4_tag_levels <- "A"

# ========================================================================
# 图5：各区域CIPS vs 通用CAF箱线图
# ========================================================================
fig5_width <- 7
fig5_height <- 5
fig5_title <- "CIPS vs General CAF Scores by Region"
fig5_title_size <- 18
fig5_title_face <- "plain"
fig5_subtitle_size <- 8
fig5_subtitle_color <- "#7F8C8D"
fig5_x_label <- ""
fig5_y_label <- "Module Score"
fig5_axis_title_size <- 16
fig5_axis_text_size <- 16
fig5_axis_text_angle <- 0
fig5_axis_text_hjust <- 0.5
fig5_axis_text_face <- "plain"
fig5_legend_title <- "Signature"
fig5_legend_text_size <- 16
fig5_legend_title_size <- 16
fig5_box_width <- 0.5
fig5_box_alpha <- 0.7
fig5_outlier_size <- 0.5
fig5_color_cips <- "#C60036"
fig5_color_caf <- "darkgreen"

# ========================================================================
# 图6：各区域IL6+ACTA2 vs 通用CAF箱线图
# ========================================================================
fig6_width <- 7
fig6_height <- 5
fig6_title <- "IL6+ACTA2+ vs General CAF Scores by Region"
fig6_title_size <- 18
fig6_title_face <- "plain"
fig6_subtitle_size <- 8
fig6_subtitle_color <- "#7F8C8D"
fig6_x_label <- ""
fig6_y_label <- "Module Score"
fig6_axis_title_size <- 16
fig6_axis_text_size <- 16
fig6_axis_text_angle <- 0
fig6_axis_text_hjust <- 0.5
fig6_axis_text_face <- "plain"
fig6_legend_title <- "Signature"
fig6_legend_text_size <- 16
fig6_legend_title_size <- 16
fig6_box_width <- 0.6
fig6_box_alpha <- 0.7
fig6_outlier_size <- 0.5
fig6_color_il6 <- "#E64B35"
fig6_color_caf <- "darkgreen"

# ========================================================================
# 图7：区域分类图
# ========================================================================
fig7_width <- 7
fig7_height <- 7
fig7_pt_size <- 1.5
fig7_title <- "Regional Classification"
fig7_title_size <- 14
fig7_title_face <- "bold"
fig7_subtitle_size <- 11
fig7_legend_text_size <- 10
fig7_legend_title_size <- 10
fig7_color_tumor <- "#C0651C"
fig7_color_interface <- "#D7191C"
fig7_color_stroma <- "#2D5A4A"
fig7_color_mixed <- "grey70"

# ========================================================================
# 图8：三区域 CIPS/CAF 比值箱线图
# ========================================================================
fig8_width <- 5
fig8_height <- 4
fig8_title <- "CIPS / CAF Ratio by Region"
fig8_title_size <- 16
fig8_title_face <- "bold"
fig8_subtitle_size <- 6
fig8_subtitle_color <- "#7F8C8D"
fig8_y_label <- "CIPS / CAF Ratio"
fig8_x_label <- ""
fig8_axis_text_size <- 20
fig8_axis_text_face <- "plain"
fig8_axis_title_size <- 20
fig8_legend_position <- "none"
fig8_violin_alpha <- 0.7
fig8_box_width <- 0.15
fig8_outlier_size <- 0.5
fig8_color_tumor <- "#B55A1A"
fig8_color_interface <- "#D7191C"
fig8_color_stroma <- "#2D5A4A"

# ========================================================================
# 图9：交界区 vs 基质区 CIPS/CAF 比值
# ========================================================================
fig9_width <- 5
fig9_height <- 5.5
fig9_title <- "CIPS / CAF Ratio"
fig9_title_size <- 14
fig9_title_face <- "bold"
fig9_subtitle_size <- 9
fig9_subtitle_color <- "#7F8C8D"
fig9_y_label <- "CIPS / CAF Ratio"
fig9_x_label <- ""
fig9_axis_text_size <- 12
fig9_axis_text_face <- "bold"
fig9_axis_title_size <- 12
fig9_legend_position <- "none"
fig9_violin_alpha <- 0.7
fig9_box_width <- 0.15
fig9_outlier_size <- 0.5

# ========================================================================
# 图10：CIPS/CAF 比值空间分布
# ========================================================================
fig10_width <- 7
fig10_height <- 6
fig10_pt_size <- 1.5
fig10_alpha_min <- 0.1
fig10_alpha_max <- 1.0
fig10_title <- "CIPS / CAF Ratio"
fig10_title_size <- 14
fig10_title_face <- "bold"
fig10_subtitle <- "Spatial distribution of CAF pro-invasive activity"
fig10_subtitle_size <- 10
fig10_legend_text_size <- 10
fig10_legend_title_size <- 10
fig10_color_low <- "grey90"
fig10_color_mid_low <- "lightblue"
fig10_color_mid_high <- "orange"
fig10_color_high <- "darkred"
fig10_legend_title <- "CIPS/CAF\nRatio"

# ========================================================================
# 图11：三区域 IL6+ACTA2/CAF 比值箱线图
# ========================================================================
fig11_width <- 5
fig11_height <- 4
fig11_title <- "IL6+ACTA2 / CAF Ratio by Region"
fig11_title_size <- 16
fig11_title_face <- "plain"
fig11_subtitle_size <- 6
fig11_subtitle_color <- "#7F8C8D"
fig11_y_label <- "IL6+ACTA2 / CAF Ratio"
fig11_x_label <- ""
fig11_axis_text_size <- 20
fig11_axis_text_face <- "plain"
fig11_axis_title_size <- 20
fig11_legend_position <- "none"
fig11_violin_alpha <- 0.7
fig11_box_width <- 0.15
fig11_outlier_size <- 0.5
fig11_color_tumor <- "#B55A1A"
fig11_color_interface <- "#D7191C"
fig11_color_stroma <- "#2D5A4A"
fig11_ratio_max <- 5

# ========================================================================
# 图12：交界区 vs 基质区 IL6+ACTA2/CAF 比值
# ========================================================================
fig12_width <- 5
fig12_height <- 5.5
fig12_title <- "IL6+ACTA2 / CAF Ratio"
fig12_title_size <- 14
fig12_title_face <- "bold"
fig12_subtitle_size <- 9
fig12_subtitle_color <- "#7F8C8D"
fig12_y_label <- "IL6+ACTA2 / CAF Ratio"
fig12_x_label <- ""
fig12_axis_text_size <- 12
fig12_axis_text_face <- "bold"
fig12_axis_title_size <- 12
fig12_legend_position <- "none"
fig12_violin_alpha <- 0.7
fig12_box_width <- 0.15
fig12_outlier_size <- 0.5

# ========================================================================
# 图13：IL6+ACTA2/CAF 比值空间分布
# ========================================================================
fig13_width <- 7
fig13_height <- 6
fig13_pt_size <- 1.5
fig13_alpha_min <- 0.1
fig13_alpha_max <- 1.0
fig13_title <- "IL6+ACTA2 / CAF Ratio"
fig13_title_size <- 14
fig13_title_face <- "bold"
fig13_subtitle <- "Spatial distribution of CAF inflammatory activation"
fig13_subtitle_size <- 10
fig13_legend_text_size <- 10
fig13_legend_title_size <- 10
fig13_color_low <- "grey90"
fig13_color_mid_low <- "lightblue"
fig13_color_mid_high <- "orange"
fig13_color_high <- "darkred"
fig13_legend_title <- "IL6+ACTA2/CAF\nRatio"

# ========================================================================
# 通用主题参数
# ========================================================================
theme_spatial_base_size <- 12
theme_spatial_legend_height <- 0.6

# ========================================================================
# 区域定义阈值
# ========================================================================
tumor_threshold_low <- 0.3
tumor_threshold_high <- 0.5
stroma_threshold_low <- 0.3
stroma_threshold_high <- 0.5

# ============================================================================
# 以下代码无需修改
# ============================================================================

# 辅助函数：格式化p值
format_p_value <- function(p_val) {
  if (is.na(p_val) || is.null(p_val)) return("P = NA")
  if (p_val < 0.0001) return("P < 0.0001")
  return(paste0("P = ", sprintf("%.4f", p_val)))
}

# 辅助函数：计算所有两两比较p值并生成副标题
build_pairwise_subtitle <- function(data, group_col, value_col, test_method = "t.test") {
  groups <- unique(data[[group_col]])
  n_groups <- length(groups)
  subtitle_parts <- c()
  
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      g1 <- groups[i]
      g2 <- groups[j]
      v1 <- data[[value_col]][data[[group_col]] == g1]
      v2 <- data[[value_col]][data[[group_col]] == g2]
      if (length(v1) >= 3 && length(v2) >= 3) {
        if (test_method == "t.test") {
          p_val <- t.test(v1, v2)$p.value
        } else {
          p_val <- wilcox.test(v1, v2)$p.value
        }
        subtitle_parts <- c(subtitle_parts, paste0(g1, " vs ", g2, ": ", format_p_value(p_val)))
      }
    }
  }
  paste(subtitle_parts, collapse = " | ")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

cat("========================================\n")
cat("CAF评分对比分析\n")
cat("========================================\n")

# -------------------- 1. 读取数据 --------------------
bc_spatial <- readRDS(file.path(work_dir, "bc_spatial_processed.rds"))

caf_markers <- c("ACTA2", "FAP", "PDGFRA", "COL1A1", "COL1A2", "PDPN", "S100A4", "VIM", "TNC", "POSTN")
caf_markers_present <- intersect(caf_markers, rownames(bc_spatial))
cat("\n通用CAF标志基因命中率:", length(caf_markers_present), "/", length(caf_markers), "\n")

# -------------------- 2. 计算通用CAF评分 --------------------
if (!"CAF_Score1" %in% colnames(bc_spatial@meta.data)) {
  cat("\n计算通用CAF评分...\n")
  bc_spatial <- AddModuleScore(
    bc_spatial,
    features = list(CAF = caf_markers_present),
    name = "CAF_Score",
    ctrl = 100
  )
}

# -------------------- 3. 检查评分 --------------------
has_cips <- "CIPS_Score1" %in% colnames(bc_spatial@meta.data)
has_il6_acta2 <- "IL6_ACTA2_Score1" %in% colnames(bc_spatial@meta.data)
has_caf <- "CAF_Score1" %in% colnames(bc_spatial@meta.data)
has_tumor <- "Tumor_Score1" %in% colnames(bc_spatial@meta.data)
has_stroma <- "Stroma_Score1" %in% colnames(bc_spatial@meta.data)

cat("\n评分检查:\n")
cat("  CIPS评分:", ifelse(has_cips, "存在", "不存在"), "\n")
cat("  IL6+ACTA2评分:", ifelse(has_il6_acta2, "存在", "不存在"), "\n")
cat("  通用CAF评分:", ifelse(has_caf, "存在", "不存在"), "\n")
cat("  肿瘤评分:", ifelse(has_tumor, "存在", "不存在"), "\n")
cat("  基质评分:", ifelse(has_stroma, "存在", "不存在"), "\n")

# -------------------- 4. 定义区域 --------------------
if (has_tumor && has_stroma) {
  tumor_scores <- bc_spatial$Tumor_Score1
  stroma_scores <- bc_spatial$Stroma_Score1
  
  tumor_scores_norm <- (tumor_scores - min(tumor_scores)) / (max(tumor_scores) - min(tumor_scores))
  stroma_scores_norm <- (stroma_scores - min(stroma_scores)) / (max(stroma_scores) - min(stroma_scores))
  
  bc_spatial$region <- "Mixed/Other"
  bc_spatial$region[tumor_scores_norm > tumor_threshold_high & stroma_scores_norm < stroma_threshold_low] <- "Tumor Core"
  bc_spatial$region[stroma_scores_norm > stroma_threshold_high & tumor_scores_norm < tumor_threshold_low] <- "Stroma Rich"
  bc_spatial$region[tumor_scores_norm >= tumor_threshold_low & tumor_scores_norm <= tumor_threshold_high &
                      stroma_scores_norm >= stroma_threshold_low & stroma_scores_norm <= stroma_threshold_high] <- "Interface"
  
  bc_spatial$region <- factor(bc_spatial$region, levels = c("Interface", "Tumor Core", "Stroma Rich", "Mixed/Other"))
  
  cat("\n区域分布:\n")
  print(table(bc_spatial$region))
}

# -------------------- 5. 通用主题 --------------------
theme_spatial <- theme_minimal(base_size = theme_spatial_base_size) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.height = unit(theme_spatial_legend_height, "cm"),
    legend.text = element_text(size = fig1_legend_text_size),
    legend.title = element_text(size = fig1_legend_title_size),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 1
  )

theme_box <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "#7F8C8D"),
    legend.position = "right",
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

# ============================================================
# 图1：通用CAF评分空间分布
# ============================================================
if (has_caf) {
  cat("\n生成图1: 通用CAF评分空间分布...\n")
  p1 <- SpatialFeaturePlot(bc_spatial, features = "CAF_Score1",
                           alpha = c(fig1_alpha_min, fig1_alpha_max),
                           pt.size.factor = fig1_pt_size) +
    scale_fill_gradientn(colors = c(fig1_color_low, fig1_color_mid, fig1_color_high), name = fig1_legend_title) +
    labs(title = fig1_title, subtitle = fig1_subtitle) +
    theme_spatial +
    theme(plot.title = element_text(size = fig1_title_size, face = fig1_title_face),
          plot.subtitle = element_text(size = fig1_subtitle_size))
  ggsave(file.path(output_dir, "Fig1_General_CAF_Spatial.pdf"), p1, width = fig1_width, height = fig1_height)
  cat("  已保存: Fig1_General_CAF_Spatial.pdf\n")
}

# ============================================================
# 图2：CIPS vs 通用CAF并排
# ============================================================
if (has_cips && has_caf) {
  cat("\n生成图2: CIPS vs 通用CAF并排比较...\n")
  p2_left <- SpatialFeaturePlot(bc_spatial, features = "CIPS_Score1",
                                alpha = c(fig2_alpha_min, fig2_alpha_max),
                                pt.size.factor = fig2_pt_size) +
    scale_fill_gradientn(colors = c(fig2_cips_low, fig2_cips_mid, fig2_cips_high), name = fig2_legend_title_left) +
    labs(title = fig2_title_left) +
    theme_spatial +
    theme(plot.title = element_text(size = fig2_title_size, face = fig2_title_face),
          legend.key.height = unit(fig2_legend_height, "cm"))
  
  p2_right <- SpatialFeaturePlot(bc_spatial, features = "CAF_Score1",
                                 alpha = c(fig2_alpha_min, fig2_alpha_max),
                                 pt.size.factor = fig2_pt_size) +
    scale_fill_gradientn(colors = c(fig2_caf_low, fig2_caf_mid, fig2_caf_high), name = fig2_legend_title_right) +
    labs(title = fig2_title_right) +
    theme_spatial +
    theme(plot.title = element_text(size = fig2_title_size, face = fig2_title_face),
          legend.key.height = unit(fig2_legend_height, "cm"))
  
  p2 <- p2_left | p2_right
  p2 <- p2 + plot_annotation(tag_levels = fig2_tag_levels)
  ggsave(file.path(output_dir, "Fig2_CIPS_vs_CAF_Spatial.pdf"), p2, width = fig2_width, height = fig2_height)
  cat("  已保存: Fig2_CIPS_vs_CAF_Spatial.pdf\n")
}

# ============================================================
# 图3：IL6+ACTA2 vs 通用CAF并排
# ============================================================
if (has_il6_acta2 && has_caf) {
  cat("\n生成图3: IL6+ACTA2 vs 通用CAF并排比较...\n")
  p3_left <- SpatialFeaturePlot(bc_spatial, features = "IL6_ACTA2_Score1",
                                alpha = c(fig3_alpha_min, fig3_alpha_max),
                                pt.size.factor = fig3_pt_size) +
    scale_fill_gradientn(colors = c(fig3_il6_low, fig3_il6_mid, fig3_il6_high), name = fig3_legend_title_left) +
    labs(title = fig3_title_left) +
    theme_spatial +
    theme(plot.title = element_text(size = fig3_title_size, face = fig3_title_face),
          legend.key.height = unit(fig3_legend_height, "cm"))
  
  p3_right <- SpatialFeaturePlot(bc_spatial, features = "CAF_Score1",
                                 alpha = c(fig3_alpha_min, fig3_alpha_max),
                                 pt.size.factor = fig3_pt_size) +
    scale_fill_gradientn(colors = c(fig3_caf_low, fig3_caf_mid, fig3_caf_high), name = fig3_legend_title_right) +
    labs(title = fig3_title_right) +
    theme_spatial +
    theme(plot.title = element_text(size = fig3_title_size, face = fig3_title_face),
          legend.key.height = unit(fig3_legend_height, "cm"))
  
  p3 <- p3_left | p3_right
  p3 <- p3 + plot_annotation(tag_levels = fig3_tag_levels)
  ggsave(file.path(output_dir, "Fig3_IL6ACTA2_vs_CAF_Spatial.pdf"), p3, width = fig3_width, height = fig3_height)
  cat("  已保存: Fig3_IL6ACTA2_vs_CAF_Spatial.pdf\n")
}

# ============================================================
# 图4：三图并排
# ============================================================
if (has_cips && has_il6_acta2 && has_caf) {
  cat("\n生成图4: 三图并排比较...\n")
  p4_left <- SpatialFeaturePlot(bc_spatial, features = "CIPS_Score1",
                                alpha = c(fig4_alpha_min, fig4_alpha_max),
                                pt.size.factor = fig4_pt_size) +
    scale_fill_gradientn(colors = c(fig4_cips_low, fig4_cips_mid, fig4_cips_high), name = "CIPS") +
    labs(title = fig4_title_left) +
    theme_spatial +
    theme(plot.title = element_text(size = fig4_title_size, face = fig4_title_face),
          legend.key.height = unit(fig4_legend_height, "cm"))
  
  p4_mid <- SpatialFeaturePlot(bc_spatial, features = "IL6_ACTA2_Score1",
                               alpha = c(fig4_alpha_min, fig4_alpha_max),
                               pt.size.factor = fig4_pt_size) +
    scale_fill_gradientn(colors = c(fig4_il6_low, fig4_il6_mid, fig4_il6_high), name = "IL6+ACTA2") +
    labs(title = fig4_title_mid) +
    theme_spatial +
    theme(plot.title = element_text(size = fig4_title_size, face = fig4_title_face),
          legend.key.height = unit(fig4_legend_height, "cm"))
  
  p4_right <- SpatialFeaturePlot(bc_spatial, features = "CAF_Score1",
                                 alpha = c(fig4_alpha_min, fig4_alpha_max),
                                 pt.size.factor = fig4_pt_size) +
    scale_fill_gradientn(colors = c(fig4_caf_low, fig4_caf_mid, fig4_caf_high), name = "CAF") +
    labs(title = fig4_title_right) +
    theme_spatial +
    theme(plot.title = element_text(size = fig4_title_size, face = fig4_title_face),
          legend.key.height = unit(fig4_legend_height, "cm"))
  
  p4 <- p4_left | p4_mid | p4_right
  p4 <- p4 + plot_annotation(tag_levels = fig4_tag_levels)
  ggsave(file.path(output_dir, "Fig4_ThreePanel_Comparison.pdf"), p4, width = fig4_width, height = fig4_height)
  cat("  已保存: Fig4_ThreePanel_Comparison.pdf\n")
}

# ============================================================
# 图5：各区域CIPS vs 通用CAF箱线图
# ============================================================
if (has_cips && has_caf && has_tumor && has_stroma) {
  cat("\n生成图5: 各区域CIPS与通用CAF评分比较...\n")
  
  plot_df <- data.frame(
    Region = bc_spatial$region,
    CIPS = bc_spatial$CIPS_Score1,
    CAF = bc_spatial$CAF_Score1
  )
  plot_df <- plot_df[plot_df$Region %in% c("Tumor Core", "Interface", "Stroma Rich"), ]
  plot_df$Region <- factor(plot_df$Region, levels = c("Tumor Core", "Interface", "Stroma Rich"))
  
  plot_df_long <- pivot_longer(plot_df, cols = c("CIPS", "CAF"), names_to = "Signature", values_to = "Score")
  plot_df_long$Signature <- factor(plot_df_long$Signature, levels = c("CIPS", "CAF"))
  
  sig_colors <- c("CIPS" = fig5_color_cips, "CAF" = fig5_color_caf)
  
  # 计算所有两两比较p值
  p_subtitle_parts <- c()
  for (sig in c("CIPS", "CAF")) {
    sub_df <- plot_df_long[plot_df_long$Signature == sig, ]
    pairwise_sub <- build_pairwise_subtitle(sub_df, "Region", "Score")
    p_subtitle_parts <- c(p_subtitle_parts, paste0(sig, ": ", pairwise_sub))
  }
  p_subtitle_full <- paste(p_subtitle_parts, collapse = " | ")
  
  p5 <- ggplot(plot_df_long, aes(x = Region, y = Score, fill = Signature)) +
    geom_boxplot(width = fig5_box_width, alpha = fig5_box_alpha, outlier.size = fig5_outlier_size) +
    scale_fill_manual(values = sig_colors, name = fig5_legend_title) +
    labs(title = fig5_title, subtitle = p_subtitle_full, x = fig5_x_label, y = fig5_y_label) +
    theme_box +
    theme(plot.title = element_text(size = fig5_title_size, face = fig5_title_face),
          plot.subtitle = element_text(size = fig5_subtitle_size, color = fig5_subtitle_color),
          axis.title = element_text(size = fig5_axis_title_size),
          axis.text.x = element_text(size = fig5_axis_text_size, angle = fig5_axis_text_angle, 
                                     hjust = fig5_axis_text_hjust, face = fig5_axis_text_face),
          legend.text = element_text(size = fig5_legend_text_size),
          legend.title = element_text(size = fig5_legend_title_size))
  
  ggsave(file.path(output_dir, "Fig5_CIPS_vs_CAF_by_Region_Boxplot.pdf"), p5, width = fig5_width, height = fig5_height)
  cat("  已保存: Fig5_CIPS_vs_CAF_by_Region_Boxplot.pdf\n")
}

# ============================================================
# 图6：各区域IL6+ACTA2 vs 通用CAF箱线图
# ============================================================
if (has_il6_acta2 && has_caf && has_tumor && has_stroma) {
  cat("\n生成图6: 各区域IL6+ACTA2与通用CAF评分比较...\n")
  
  plot_df <- data.frame(
    Region = bc_spatial$region,
    IL6_ACTA2 = bc_spatial$IL6_ACTA2_Score1,
    CAF = bc_spatial$CAF_Score1
  )
  plot_df <- plot_df[plot_df$Region %in% c("Tumor Core", "Interface", "Stroma Rich"), ]
  plot_df$Region <- factor(plot_df$Region, levels = c("Tumor Core", "Interface", "Stroma Rich"))
  
  plot_df_long <- pivot_longer(plot_df, cols = c("IL6_ACTA2", "CAF"), names_to = "Signature", values_to = "Score")
  plot_df_long$Signature <- factor(plot_df_long$Signature, levels = c("IL6_ACTA2", "CAF"))
  
  sig_colors <- c("IL6_ACTA2" = fig6_color_il6, "CAF" = fig6_color_caf)
  
  p_subtitle_parts <- c()
  for (sig in c("IL6_ACTA2", "CAF")) {
    sub_df <- plot_df_long[plot_df_long$Signature == sig, ]
    pairwise_sub <- build_pairwise_subtitle(sub_df, "Region", "Score")
    p_subtitle_parts <- c(p_subtitle_parts, paste0(sig, ": ", pairwise_sub))
  }
  p_subtitle_full <- paste(p_subtitle_parts, collapse = " | ")
  
  p6 <- ggplot(plot_df_long, aes(x = Region, y = Score, fill = Signature)) +
    geom_boxplot(width = fig6_box_width, alpha = fig6_box_alpha, outlier.size = fig6_outlier_size) +
    scale_fill_manual(values = sig_colors, name = fig6_legend_title) +
    labs(title = fig6_title, subtitle = p_subtitle_full, x = fig6_x_label, y = fig6_y_label) +
    theme_box +
    theme(plot.title = element_text(size = fig6_title_size, face = fig6_title_face),
          plot.subtitle = element_text(size = fig6_subtitle_size, color = fig6_subtitle_color),
          axis.title = element_text(size = fig6_axis_title_size),
          axis.text.x = element_text(size = fig6_axis_text_size, angle = fig6_axis_text_angle, 
                                     hjust = fig6_axis_text_hjust, face = fig6_axis_text_face),
          legend.text = element_text(size = fig6_legend_text_size),
          legend.title = element_text(size = fig6_legend_title_size))
  
  ggsave(file.path(output_dir, "Fig6_IL6ACTA2_vs_CAF_by_Region_Boxplot.pdf"), p6, width = fig6_width, height = fig6_height)
  cat("  已保存: Fig6_IL6ACTA2_vs_CAF_by_Region_Boxplot.pdf\n")
}

# ============================================================
# 图7：区域分类图
# ============================================================
if (has_tumor && has_stroma) {
  cat("\n生成图7: 区域分类图...\n")
  region_colors <- c("Tumor Core" = fig7_color_tumor, "Stroma Rich" = fig7_color_stroma,
                     "Interface" = fig7_color_interface, "Mixed/Other" = fig7_color_mixed)
  p7 <- SpatialDimPlot(bc_spatial, group.by = "region", cols = region_colors, pt.size.factor = fig7_pt_size) +
    labs(title = fig7_title) +
    theme_spatial +
    theme(plot.title = element_text(size = fig7_title_size, face = fig7_title_face))
  ggsave(file.path(output_dir, "Fig7_Regional_Classification.pdf"), p7, width = fig7_width, height = fig7_height)
  cat("  已保存: Fig7_Regional_Classification.pdf\n")
}

# ============================================================
# CIPS/CAF 比值分析
# ============================================================
cat("\n========== 生成 CIPS/CAF 比值分析 ==========\n")

if (has_cips && has_caf && has_tumor && has_stroma) {
  
  bc_spatial$CIPS_CAF_ratio <- bc_spatial$CIPS_Score1 / bc_spatial$CAF_Score1
  
  plot_df_cips <- data.frame(
    Region = bc_spatial$region,
    CIPS = bc_spatial$CIPS_Score1,
    CAF = bc_spatial$CAF_Score1,
    Ratio = bc_spatial$CIPS_CAF_ratio
  )
  plot_df_cips <- plot_df_cips[plot_df_cips$Region %in% c("Tumor Core", "Interface", "Stroma Rich"), ]
  plot_df_cips$Region <- factor(plot_df_cips$Region, levels = c("Tumor Core", "Interface", "Stroma Rich"))
  plot_df_cips_clean <- plot_df_cips[is.finite(plot_df_cips$Ratio) & plot_df_cips$Ratio > 0 & plot_df_cips$Ratio < 2, ]
  
  # 图8：三区域 CIPS/CAF 比值
  ratio_cips_colors <- c("Tumor Core" = fig8_color_tumor, "Interface" = fig8_color_interface, "Stroma Rich" = fig8_color_stroma)
  pval_sub8 <- build_pairwise_subtitle(plot_df_cips_clean, "Region", "Ratio")
  
  p8 <- ggplot(plot_df_cips_clean, aes(x = Region, y = Ratio, fill = Region)) +
    geom_violin(alpha = fig8_violin_alpha, trim = FALSE) +
    geom_boxplot(width = fig8_box_width, fill = "white", alpha = 0.5, outlier.size = fig8_outlier_size) +
    scale_fill_manual(values = ratio_cips_colors) +
    labs(title = fig8_title, subtitle = pval_sub8, y = fig8_y_label, x = fig8_x_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = fig8_legend_position,
          plot.title = element_text(hjust = 0.5, face = fig8_title_face, size = fig8_title_size),
          plot.subtitle = element_text(hjust = 0.5, size = fig8_subtitle_size, color = fig8_subtitle_color),
          axis.text.x = element_text(size = fig8_axis_text_size, face = fig8_axis_text_face),
          axis.title = element_text(size = fig8_axis_title_size))
  ggsave(file.path(output_dir, "Fig8_CIPS_CAF_Ratio_All_Regions.pdf"), p8, width = fig8_width, height = fig8_height)
  cat("  已保存: Fig8_CIPS_CAF_Ratio_All_Regions.pdf\n")
  
  # 图9：交界区 vs 基质区
  plot_df_cips_two <- plot_df_cips_clean[plot_df_cips_clean$Region %in% c("Interface", "Stroma Rich"), ]
  plot_df_cips_two$Region <- factor(plot_df_cips_two$Region, levels = c("Interface", "Stroma Rich"))
  
  pval_sub9 <- build_pairwise_subtitle(plot_df_cips_two, "Region", "Ratio")
  
  p9 <- ggplot(plot_df_cips_two, aes(x = Region, y = Ratio, fill = Region)) +
    geom_violin(alpha = fig9_violin_alpha, trim = FALSE) +
    geom_boxplot(width = fig9_box_width, fill = "white", alpha = 0.5, outlier.size = fig9_outlier_size) +
    scale_fill_manual(values = ratio_cips_colors) +
    labs(title = fig9_title, subtitle = pval_sub9, y = fig9_y_label, x = fig9_x_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = fig9_legend_position,
          plot.title = element_text(hjust = 0.5, face = fig9_title_face, size = fig9_title_size),
          plot.subtitle = element_text(hjust = 0.5, size = fig9_subtitle_size, color = fig9_subtitle_color),
          axis.text.x = element_text(size = fig9_axis_text_size, face = fig9_axis_text_face),
          axis.title = element_text(size = fig9_axis_title_size))
  ggsave(file.path(output_dir, "Fig9_CIPS_CAF_Ratio_Interface_vs_Stroma.pdf"), p9, width = fig9_width, height = fig9_height)
  cat("  已保存: Fig9_CIPS_CAF_Ratio_Interface_vs_Stroma.pdf\n")
  
  # 图10：CIPS/CAF 比值空间分布
  p10 <- SpatialFeaturePlot(bc_spatial, features = "CIPS_CAF_ratio",
                            alpha = c(fig10_alpha_min, fig10_alpha_max), pt.size.factor = fig10_pt_size) +
    scale_fill_gradientn(colors = c(fig10_color_low, fig10_color_mid_low, fig10_color_mid_high, fig10_color_high),
                         name = fig10_legend_title) +
    labs(title = fig10_title, subtitle = fig10_subtitle) +
    theme_spatial +
    theme(plot.title = element_text(size = fig10_title_size, face = fig10_title_face),
          plot.subtitle = element_text(size = fig10_subtitle_size))
  ggsave(file.path(output_dir, "Fig10_CIPS_CAF_Ratio_Spatial.pdf"), p10, width = fig10_width, height = fig10_height)
  cat("  已保存: Fig10_CIPS_CAF_Ratio_Spatial.pdf\n")
}

# ============================================================
# IL6+ACTA2/CAF 比值分析
# ============================================================
cat("\n========== 生成 IL6+ACTA2 / CAF 比值分析 ==========\n")

if (has_il6_acta2 && has_caf && has_tumor && has_stroma) {
  
  bc_spatial$IL6_ACTA2_CAF_ratio <- bc_spatial$IL6_ACTA2_Score1 / bc_spatial$CAF_Score1
  
  plot_df_il6 <- data.frame(
    Region = bc_spatial$region,
    Ratio = bc_spatial$IL6_ACTA2_CAF_ratio
  )
  plot_df_il6 <- plot_df_il6[plot_df_il6$Region %in% c("Tumor Core", "Interface", "Stroma Rich"), ]
  plot_df_il6$Region <- factor(plot_df_il6$Region, levels = c("Tumor Core", "Interface", "Stroma Rich"))
  plot_df_il6_clean <- plot_df_il6[is.finite(plot_df_il6$Ratio) & plot_df_il6$Ratio > 0 & plot_df_il6$Ratio < fig11_ratio_max, ]
  
  # 图11：三区域 IL6+ACTA2/CAF 比值
  ratio_il6_colors <- c("Tumor Core" = fig11_color_tumor, "Interface" = fig11_color_interface, "Stroma Rich" = fig11_color_stroma)
  pval_sub11 <- build_pairwise_subtitle(plot_df_il6_clean, "Region", "Ratio")
  
  p11 <- ggplot(plot_df_il6_clean, aes(x = Region, y = Ratio, fill = Region)) +
    geom_violin(alpha = fig11_violin_alpha, trim = FALSE) +
    geom_boxplot(width = fig11_box_width, fill = "white", alpha = 0.5, outlier.size = fig11_outlier_size) +
    scale_fill_manual(values = ratio_il6_colors) +
    labs(title = fig11_title, subtitle = pval_sub11, y = fig11_y_label, x = fig11_x_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = fig11_legend_position,
          plot.title = element_text(hjust = 0.5, face = fig11_title_face, size = fig11_title_size),
          plot.subtitle = element_text(hjust = 0.5, size = fig11_subtitle_size, color = fig11_subtitle_color),
          axis.text.x = element_text(size = fig11_axis_text_size, face = fig11_axis_text_face),
          axis.title = element_text(size = fig11_axis_title_size))
  ggsave(file.path(output_dir, "Fig11_IL6_ACTA2_CAF_Ratio_All_Regions.pdf"), p11, width = fig11_width, height = fig11_height)
  cat("  已保存: Fig11_IL6_ACTA2_CAF_Ratio_All_Regions.pdf\n")
  
  # 图12：交界区 vs 基质区
  plot_df_il6_two <- plot_df_il6_clean[plot_df_il6_clean$Region %in% c("Interface", "Stroma Rich"), ]
  plot_df_il6_two$Region <- factor(plot_df_il6_two$Region, levels = c("Interface", "Stroma Rich"))
  
  pval_sub12 <- build_pairwise_subtitle(plot_df_il6_two, "Region", "Ratio")
  
  p12 <- ggplot(plot_df_il6_two, aes(x = Region, y = Ratio, fill = Region)) +
    geom_violin(alpha = fig12_violin_alpha, trim = FALSE) +
    geom_boxplot(width = fig12_box_width, fill = "white", alpha = 0.5, outlier.size = fig12_outlier_size) +
    scale_fill_manual(values = ratio_il6_colors) +
    labs(title = fig12_title, subtitle = pval_sub12, y = fig12_y_label, x = fig12_x_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = fig12_legend_position,
          plot.title = element_text(hjust = 0.5, face = fig12_title_face, size = fig12_title_size),
          plot.subtitle = element_text(hjust = 0.5, size = fig12_subtitle_size, color = fig12_subtitle_color),
          axis.text.x = element_text(size = fig12_axis_text_size, face = fig12_axis_text_face),
          axis.title = element_text(size = fig12_axis_title_size))
  ggsave(file.path(output_dir, "Fig12_IL6_ACTA2_CAF_Ratio_Interface_vs_Stroma.pdf"), p12, width = fig12_width, height = fig12_height)
  cat("  已保存: Fig12_IL6_ACTA2_CAF_Ratio_Interface_vs_Stroma.pdf\n")
  
  # 图13：IL6+ACTA2/CAF 比值空间分布
  p13 <- SpatialFeaturePlot(bc_spatial, features = "IL6_ACTA2_CAF_ratio",
                            alpha = c(fig13_alpha_min, fig13_alpha_max), pt.size.factor = fig13_pt_size) +
    scale_fill_gradientn(colors = c(fig13_color_low, fig13_color_mid_low, fig13_color_mid_high, fig13_color_high),
                         name = fig13_legend_title) +
    labs(title = fig13_title, subtitle = fig13_subtitle) +
    theme_spatial +
    theme(plot.title = element_text(size = fig13_title_size, face = fig13_title_face),
          plot.subtitle = element_text(size = fig13_subtitle_size))
  ggsave(file.path(output_dir, "Fig13_IL6_ACTA2_CAF_Ratio_Spatial.pdf"), p13, width = fig13_width, height = fig13_height)
  cat("  已保存: Fig13_IL6_ACTA2_CAF_Ratio_Spatial.pdf\n")
}

# ============================================================
# 输出统计信息
# ============================================================
sink(file.path(output_dir, "CAF_Comparison_Statistics.txt"))
cat("========================================\n")
cat("CAF评分对比统计报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("1. 通用CAF标志基因:\n")
cat("   - 原始基因数:", length(caf_markers), "\n")
cat("   - 命中基因数:", length(caf_markers_present), "\n")
cat("   - 命中率:", round(length(caf_markers_present)/length(caf_markers)*100, 1), "%\n\n")
if (has_tumor && has_stroma) {
  cat("2. 区域分布:\n")
  print(table(bc_spatial$region))
}
sink()

cat("\n统计报告已保存: CAF_Comparison_Statistics.txt\n")
cat("\n========================================\n")
cat("所有图片生成完成！\n")
cat("输出目录:", output_dir, "\n")
cat("========================================\n")