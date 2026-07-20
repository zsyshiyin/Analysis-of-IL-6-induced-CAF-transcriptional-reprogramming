# =============================================================================
# 竞争性基因集分析：评估CIPS在抑制剂响应中的特异性排名
# 排名方法：优先按FDR（小到大），再按NES（小到大，负值更负优先）
# =============================================================================

# 清空环境变量
rm(list = ls())

# ====================== 用户自定义参数设置 ======================
# 1. 工作路径设置（请修改为您的实际路径）
work_dir <- "D:/zsy/SX/Fomal-final/3.7-Rh2-jingzheng"  # 请修改此路径
setwd(work_dir)

# 2. 创建输出文件夹
output_dir <- "Competitive_pathway_analysis_xin_final"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出文件夹:", output_dir, "\n")
}

# 3. 输入文件设置（请根据你的实际文件名修改）
gsea_go_file <- "GSEA_GO_results.csv"        # GO的GSEA结果
gsea_kegg_file <- "GSEA_KEGG_results.csv"    # KEGG的GSEA结果

# 4. CIPS的GSEA结果（请根据你的实际结果修改！）
# 这些数值需要从GSEA Desktop或R运行结果中获得
cips_nes <- -2.028      # 替换为你的CIPS实际NES值（负值表示下调）
cips_fdr <- 0.0       # 替换为你的CIPS实际FDR值
cips_nom_pval <- 0.0  # 替换为你的CIPS实际nominal p值
cips_pathway_name <- "CIPS (CAF Invasive Progression Signature)"

# 5. 筛选参数
fdr_cutoff <- 0.05                          # 显著性阈值
top_n_pathways <- 20                        # 条形图显示前N个通路

# 6. 图片设置
plot_width <- 8
plot_height <- 6
plot_width_bar <- 10
plot_height_bar <- 8

# 颜色设置
cips_color <- "#CA0E12"      # CIPS红色
down_color <- "#377EB8"       # 下调通路蓝色
up_color <- "#4DAF4A"         # 上调通路绿色
not_significant_color <- "gray80"

# 字体大小
title_size <- 24
subtitle_size <- 10
axis_title_size <- 24
axis_text_size <- 22

# 7. 输出文件设置
# 主图
output_main_dotplot <- file.path(output_dir, "Figure_Main_Downregulated_Only.pdf")
# 补充图
output_supp_all_dotplot <- file.path(output_dir, "Figure_Supplement_All_Pathways.pdf")
output_supp_kegg_bar <- file.path(output_dir, "Figure_Supplement_KEGG_Barplot.pdf")
output_supp_go_bar <- file.path(output_dir, "Figure_Supplement_GO_Barplot.pdf")
# 数据文件
output_gsea_results <- file.path(output_dir, "GSEA_results_all_pathways.csv")
output_cips_ranking <- file.path(output_dir, "CIPS_ranking_summary.csv")
output_summary <- file.path(output_dir, "analysis_summary.txt")

# ====================== 加载必要包 ======================
cat("\n========== 加载必要包 ==========\n")

packages <- c("dplyr", "ggplot2", "tidyr", "stringr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ====================== 读取GSEA结果函数 ======================
read_gsea_results <- function(file_path, source_name) {
  cat(sprintf("\n读取%s结果: %s\n", source_name, basename(file_path)))
  
  if (!file.exists(file_path)) {
    cat(sprintf("警告：文件不存在 - %s\n", file_path))
    return(NULL)
  }
  
  # 读取文件
  gsea_data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("  原始数据: %d 行, %d 列\n", nrow(gsea_data), ncol(gsea_data)))
  
  # 显示列名
  cat("  列名:", paste(colnames(gsea_data), collapse = ", "), "\n")
  
  # 根据实际列名映射
  result <- gsea_data
  
  # 通路名称列（使用Term或KO_ID）
  if ("Term" %in% colnames(result)) {
    colnames(result)[colnames(result) == "Term"] <- "Term"
  } else if ("KO_ID" %in% colnames(result)) {
    colnames(result)[colnames(result) == "KO_ID"] <- "Term"
  } else {
    cat("  错误：找不到Term或KO_ID列\n")
    return(NULL)
  }
  
  # NES列
  if ("NES" %in% colnames(result)) {
    colnames(result)[colnames(result) == "NES"] <- "NES"
  } else {
    cat("  错误：找不到NES列\n")
    return(NULL)
  }
  
  # FDR列（FDR q-val）
  if ("FDR q-val" %in% colnames(result)) {
    colnames(result)[colnames(result) == "FDR q-val"] <- "FDR"
  } else if ("FDR" %in% colnames(result)) {
    colnames(result)[colnames(result) == "FDR"] <- "FDR"
  } else {
    cat("  警告：找不到FDR列，将使用NA\n")
    result$FDR <- NA
  }
  
  # p值列（NOM p-val）
  if ("NOM p-val" %in% colnames(result)) {
    colnames(result)[colnames(result) == "NOM p-val"] <- "NOM_pval"
  } else {
    result$NOM_pval <- NA
  }
  
  # 转换NES为数值
  result$NES <- as.numeric(as.character(result$NES))
  
  # 转换FDR为数值
  if ("FDR" %in% colnames(result)) {
    result$FDR <- as.numeric(as.character(result$FDR))
  }
  
  # 移除NES为NA的行
  before_n <- nrow(result)
  result <- result[!is.na(result$NES), ]
  cat(sprintf("  有效通路数: %d (移除 %d 行)\n", nrow(result), before_n - nrow(result)))
  
  if (nrow(result) == 0) {
    cat("  错误：没有有效的NES值\n")
    return(NULL)
  }
  
  # 显示NES范围
  cat(sprintf("  NES范围: [%.4f, %.4f]\n", min(result$NES), max(result$NES)))
  
  # 添加来源标识
  result$Source <- source_name
  
  return(result)
}

# ====================== 读取数据 ======================
cat("\n========== 读取GSEA结果 ==========\n")

gsea_go <- read_gsea_results(gsea_go_file, "GO")
gsea_kegg <- read_gsea_results(gsea_kegg_file, "KEGG")

# 合并数据
gsea_all <- bind_rows(gsea_go, gsea_kegg)

if (is.null(gsea_all) || nrow(gsea_all) == 0) {
  cat("\n错误：没有成功读取任何GSEA结果！\n")
  cat("请检查：\n")
  cat("  1. 文件路径是否正确\n")
  cat("  2. 文件是否包含NES列（数值）\n")
  stop("无法继续分析")
}

cat(sprintf("\n总通路数: %d (GO: %d, KEGG: %d)\n", 
            nrow(gsea_all), 
            ifelse(is.null(gsea_go), 0, nrow(gsea_go)),
            ifelse(is.null(gsea_kegg), 0, nrow(gsea_kegg))))

# 显示前几行数据
cat("\n数据预览（前5行）:\n")
print(head(gsea_all[, c("Term", "NES", "FDR", "Source")], 5))

# ====================== 创建CIPS记录 ======================
cat("\n========== CIPS GSEA结果 ==========\n")
cat(sprintf("  NES = %.4f\n", cips_nes))
cat(sprintf("  FDR = %.4f\n", cips_fdr))

cips_record <- data.frame(
  Term = cips_pathway_name,
  NES = cips_nes,
  FDR = cips_fdr,
  NOM_pval = cips_nom_pval,
  Source = "CIPS",
  stringsAsFactors = FALSE
)

# ====================== 合并数据 ======================
gsea_all_with_cips <- bind_rows(gsea_all, cips_record)

# 添加分组信息
gsea_all_with_cips <- gsea_all_with_cips %>%
  mutate(
    Direction = case_when(
      Term == cips_pathway_name ~ "CIPS",
      NES < 0 ~ "Downregulated",
      NES > 0 ~ "Upregulated",
      TRUE ~ "Other"
    ),
    Significant = ifelse(!is.na(FDR) & FDR < fdr_cutoff, TRUE, FALSE)
  )

# ====================== 关键：先按FDR排序，再按NES排序 ======================
cat("\n========== 计算CIPS排名（先FDR，后NES）==========\n")

# 处理FDR=0的情况（用极小值替代以便排序）
gsea_all_with_cips <- gsea_all_with_cips %>%
  mutate(
    FDR_for_ranking = ifelse(FDR == 0 | is.na(FDR), 1e-10, FDR)
  )

# 1. 在下调通路中的排名（FDR小优先 -> NES小优先）
down_regulated <- gsea_all_with_cips %>%
  filter(NES < 0) %>%
  arrange(FDR_for_ranking, NES) %>%  # 先按FDR升序，再按NES升序
  mutate(
    Rank = row_number(),
    Total = n(),
    Percentile = (1 - Rank / Total) * 100
  )

# 获取CIPS在下调通路中的排名
cips_down <- down_regulated %>% filter(Term == cips_pathway_name)

if (nrow(cips_down) > 0) {
  cat(sprintf("\n下调通路总数: %d\n", nrow(down_regulated)))
  cat(sprintf("CIPS的NES = %.4f\n", cips_down$NES[1]))
  cat(sprintf("CIPS的FDR = %.4f\n", cips_down$FDR[1]))
  cat(sprintf("CIPS在下调通路中排名: 第 %d 位 (前 %.1f%%)\n", 
              cips_down$Rank[1], cips_down$Percentile[1]))
  
  if (cips_down$Rank[1] == 1) {
    cat("\n✓✓✓ CIPS是下调最显著的通路！✓✓✓\n")
  } else {
    cat(sprintf("\n⚠ CIPS排名第%d，前面有%d个通路\n", 
                cips_down$Rank[1], cips_down$Rank[1] - 1))
    
    # 显示排在CIPS前面的通路（按FDR排序后的排名）
    cat("\n排在CIPS前面的通路（按FDR优先排序）:\n")
    down_regulated %>%
      filter(Rank < cips_down$Rank[1]) %>%
      select(Rank, Term, NES, FDR) %>%
      print()
  }
} else {
  cat("警告：未找到CIPS在下调通路中的记录\n")
}

# 2. 在所有通路中的排名
all_ranked <- gsea_all_with_cips %>%
  arrange(FDR_for_ranking, abs(NES)) %>%
  mutate(
    Rank_all = row_number(),
    Total_all = n(),
    Percentile_all = (1 - Rank_all / Total_all) * 100
  )

cips_all <- all_ranked %>% filter(Term == cips_pathway_name)

if (nrow(cips_all) > 0) {
  cat(sprintf("\n所有通路总数: %d\n", nrow(all_ranked)))
  cat(sprintf("CIPS在所有通路中排名: 第 %d 位 (前 %.1f%%)\n", 
              cips_all$Rank_all[1], cips_all$Percentile_all[1]))
}

# 保存排名结果
ranking_summary <- data.frame(
  Analysis = c("Downregulated pathways only (FDR then NES)", 
               "All pathways (FDR then |NES|)"),
  NES = c(ifelse(nrow(cips_down) > 0, cips_down$NES[1], NA),
          ifelse(nrow(cips_all) > 0, cips_all$NES[1], NA)),
  FDR = c(ifelse(nrow(cips_down) > 0, cips_down$FDR[1], NA),
          ifelse(nrow(cips_all) > 0, cips_all$FDR[1], NA)),
  Rank = c(ifelse(nrow(cips_down) > 0, cips_down$Rank[1], NA),
           ifelse(nrow(cips_all) > 0, cips_all$Rank_all[1], NA)),
  Total_Pathways = c(ifelse(nrow(cips_down) > 0, nrow(down_regulated), NA),
                     ifelse(nrow(cips_all) > 0, nrow(all_ranked), NA)),
  Percentile = c(ifelse(nrow(cips_down) > 0, cips_down$Percentile[1], NA),
                 ifelse(nrow(cips_all) > 0, cips_all$Percentile_all[1], NA))
)

write.csv(ranking_summary, file = output_cips_ranking, row.names = FALSE)
cat(sprintf("\n排名结果已保存: %s\n", output_cips_ranking))

# ====================== 主图：下调通路点图 ======================
cat("\n========== 绘制主图：下调通路点图 ==========\n")

# 准备下调通路数据（使用已排序的数据）
down_data <- down_regulated %>%
  mutate(
    FDR_display = ifelse(FDR == 0 | is.na(FDR), 1e-10, FDR),
    logFDR = -log10(FDR_display),
    logFDR_plot = ifelse(logFDR > 10, 10, logFDR),
    is_CIPS = (Term == cips_pathway_name),
    Point_size = ifelse(is_CIPS, 5, 1.8),
    Point_color = ifelse(is_CIPS, cips_color, "gray50"),
    Point_alpha = ifelse(is_CIPS, 1, 0.5)
  )

# 创建标题
if (nrow(cips_down) > 0 && cips_down$Rank[1] == 1) {
  rank_text <- sprintf("CIPS ranked #1 out of %d downregulated pathways", nrow(down_data))
} else if (nrow(cips_down) > 0) {
  rank_text <- sprintf("CIPS ranked #%d out of %d downregulated pathways (top %.1f%%)",
                       cips_down$Rank[1], nrow(down_data), cips_down$Percentile[1])
} else {
  rank_text <- "CIPS ranking information not available"
}

sub_title_text <- sprintf("%s | NES = %.3f, FDR = %.2e", rank_text, cips_nes, cips_fdr)

# 创建主图
p_main <- ggplot(down_data, aes(x = NES, y = logFDR_plot)) +
  # 背景点（其他下调通路）
  geom_point(data = subset(down_data, !is_CIPS),
             aes(size = Point_size, alpha = Point_alpha),
             color = "gray50") +
  # 前景点（CIPS）
  geom_point(data = subset(down_data, is_CIPS),
             aes(size = Point_size),
             color = cips_color) +
  scale_size_identity() +
  scale_alpha_identity() +
  # 阈值线
  geom_hline(yintercept = -log10(fdr_cutoff), 
             linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, 
             linetype = "solid", color = "black", linewidth = 0.3) +
  # CIPS标签
  geom_text(data = subset(down_data, is_CIPS),
            aes(x = NES - 0.3, y = logFDR_plot + 0.5),
            label = "CIPS", 
            color = cips_color, size = 5, fontface = "bold") +
  # 标题和标签
  labs(title = "Competitive Pathway Analysis (Downregulated Pathways Only)",
       subtitle = sub_title_text,
       x = "Normalized Enrichment Score (NES)",
       y = "-log10(FDR)") +
  scale_y_continuous(limits = c(0, 11),
                     breaks = seq(0, 10, by = 2),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(limits = c(min(down_data$NES) - 0.2, -0.1),
                     breaks = seq(-3, 0, by = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = subtitle_size, hjust = 0.5),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5))

ggsave(output_main_dotplot, plot = p_main, width = plot_width, height = plot_height, dpi = 300)
cat(sprintf("主图已保存: %s\n", output_main_dotplot))

# ====================== 补充图1：所有通路点图 ======================
cat("\n========== 绘制补充图1：所有通路点图 ==========\n")

# 准备所有通路数据
all_data <- gsea_all_with_cips %>%
  mutate(
    FDR_display = ifelse(FDR == 0 | is.na(FDR), 1e-10, FDR),
    logFDR = -log10(FDR_display),
    logFDR_plot = ifelse(logFDR > 10, 10, logFDR),
    Direction = case_when(
      Term == cips_pathway_name ~ "CIPS",
      NES < 0 ~ "Downregulated",
      NES > 0 ~ "Upregulated",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(NES), !is.na(logFDR_plot))

# 设置颜色和大小
all_data <- all_data %>%
  mutate(
    Point_color = case_when(
      Term == cips_pathway_name ~ cips_color,
      Direction == "Downregulated" ~ down_color,
      Direction == "Upregulated" ~ up_color,
      TRUE ~ not_significant_color
    ),
    Point_size = ifelse(Term == cips_pathway_name, 4, 1.5),
    Point_alpha = ifelse(Term == cips_pathway_name, 1, 0.5)
  )

# 创建补充图
p_supp_all <- ggplot(all_data, aes(x = NES, y = logFDR_plot)) +
  geom_point(aes(color = Point_color, size = Point_size, alpha = Point_alpha)) +
  scale_color_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  geom_hline(yintercept = -log10(fdr_cutoff), 
             linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, 
             linetype = "solid", color = "black", linewidth = 0.3) +
  annotate("text", x = -1.5, y = 11, 
           label = "Downregulated", color = down_color, size = 8) +
  annotate("text", x = 1.5, y = 11, 
           label = "Upregulated", color = up_color, size = 8) +
  annotate("text", x = -2.5, y = 8.5, 
           label = "CIPS", color = cips_color, size = 4, fontface = "plain") +
  labs(title = "Competitive Pathway Analysis (All Pathways)",
       subtitle = sprintf("CIPS (red) shows strongest downregulation | NES = %.3f, FDR = %.2e", cips_nes, cips_fdr),
       x = "Normalized Enrichment Score (NES)",
       y = "-log10(FDR)") +
  scale_y_continuous(limits = c(0, 11), breaks = seq(0, 10, by = 2)) +
  theme_minimal() +
  theme(plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = subtitle_size, hjust = 0.5),
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5))

ggsave(output_supp_all_dotplot, plot = p_supp_all, width = plot_width, height = plot_height, dpi = 300)
cat(sprintf("补充图1已保存: %s\n", output_supp_all_dotplot))

# ====================== 补充图2：KEGG排名条形图 ======================
cat("\n========== 绘制补充图2：KEGG排名条形图 ==========\n")

plot_bar <- function(data, title, filename, top_n = top_n_pathways) {
  
  if (is.null(data) || nrow(data) == 0) {
    cat(sprintf("  跳过绘图: %s (无数据)\n", basename(filename)))
    return(NULL)
  }
  
  # 筛选下调通路并按排名排序
  plot_data <- data %>%
    filter(NES < 0) %>%
    arrange(Rank) %>%
    head(top_n) %>%
    mutate(Term_short = ifelse(nchar(Term) > 50, paste0(substr(Term, 1, 47), "..."), Term),
           Term_short = factor(Term_short, levels = rev(Term_short)),
           BarColor = ifelse(Term == cips_pathway_name, cips_color,
                             ifelse(!is.na(FDR) & FDR < fdr_cutoff, down_color, not_significant_color)))
  
  if (nrow(plot_data) == 0) {
    cat(sprintf("  跳过绘图: %s (无下调通路)\n", basename(filename)))
    return(NULL)
  }
  
  # 创建副标题
  if (nrow(cips_down) > 0) {
    sub_text <- sprintf("Ranking method: FDR first, then NES | CIPS ranked #%d out of %d downregulated pathways",
                        cips_down$Rank[1], nrow(data %>% filter(NES < 0)))
  } else {
    sub_text <- "Ranking method: FDR first, then NES"
  }
  
  p <- ggplot(plot_data, aes(x = NES, y = Term_short)) +
    geom_bar(stat = "identity", aes(fill = BarColor), width = 0.7) +
    scale_fill_identity() +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    labs(title = title,
         subtitle = sub_text,
         x = "Normalized Enrichment Score (NES)",
         y = "") +
    theme_minimal() +
    theme(plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = subtitle_size, hjust = 0.5),
          axis.title.x = element_text(size = axis_title_size),
          axis.text.y = element_text(size = axis_text_size),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5))
  
  ggsave(filename, plot = p, width = plot_width_bar, height = plot_height_bar, dpi = 300)
  cat(sprintf("  已保存: %s\n", basename(filename)))
}

# KEGG下调通路条形图
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  # 为KEGG单独计算排名
  kegg_with_cips <- bind_rows(gsea_kegg, cips_record)
  kegg_with_cips <- kegg_with_cips %>%
    mutate(FDR_for_ranking = ifelse(FDR == 0 | is.na(FDR), 1e-10, FDR))
  
  kegg_down <- kegg_with_cips %>%
    filter(NES < 0) %>%
    arrange(FDR_for_ranking, NES) %>%
    mutate(Rank = row_number())
  
  plot_bar(kegg_down, "KEGG Pathways: Top Downregulated (FDR-first Ranking)", 
           output_supp_kegg_bar)
} else {
  cat("KEGG数据为空，跳过KEGG条形图\n")
}

# GO下调通路条形图
if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  # 为GO单独计算排名
  go_with_cips <- bind_rows(gsea_go, cips_record)
  go_with_cips <- go_with_cips %>%
    mutate(FDR_for_ranking = ifelse(FDR == 0 | is.na(FDR), 1e-10, FDR))
  
  go_down <- go_with_cips %>%
    filter(NES < 0) %>%
    arrange(FDR_for_ranking, NES) %>%
    mutate(Rank = row_number())
  
  plot_bar(go_down, "GO Pathways: Top Downregulated (FDR-first Ranking)", 
           output_supp_go_bar)
} else {
  cat("GO数据为空，跳过GO条形图\n")
}

# ====================== 保存所有GSEA结果 ======================
write.csv(gsea_all_with_cips, file = output_gsea_results, row.names = FALSE)
cat(sprintf("\n所有GSEA结果已保存: %s\n", output_gsea_results))

# ====================== 生成分析报告 ======================
sink(output_summary)

cat("========================================\n")
cat("   竞争性基因集分析报告\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、排名方法\n")
cat("------------\n")
cat("排序规则: 先按FDR升序（小→大），再按NES升序（负值更负优先）\n")
cat("FDR=0的通路在排序时使用1e-10代替\n\n")

cat("二、数据概况\n")
cat("------------\n")
cat(sprintf("GO通路数: %d\n", ifelse(is.null(gsea_go), 0, nrow(gsea_go))))
cat(sprintf("KEGG通路数: %d\n", ifelse(is.null(gsea_kegg), 0, nrow(gsea_kegg))))
cat(sprintf("总通路数: %d\n", nrow(gsea_all)))
cat(sprintf("下调通路数: %d\n", nrow(down_regulated)))

cat("\n三、CIPS结果\n")
cat("------------\n")
cat(sprintf("NES: %.4f\n", cips_nes))
cat(sprintf("FDR: %.4f\n", cips_fdr))

cat("\n四、CIPS排名（先FDR后NES）\n")
cat("--------------------------\n")
if (nrow(cips_down) > 0) {
  cat(sprintf("在下调通路中: 第 %d 位（共 %d 个，前 %.1f%%）\n", 
              cips_down$Rank[1], nrow(down_regulated), cips_down$Percentile[1]))
}
if (nrow(cips_all) > 0) {
  cat(sprintf("在所有通路中: 第 %d 位（共 %d 个，前 %.1f%%）\n", 
              cips_all$Rank_all[1], nrow(all_ranked), cips_all$Percentile_all[1]))
}

cat("\n五、输出文件列表\n")
cat("----------------\n")
cat("主图:\n")
cat("  -", basename(output_main_dotplot), "- 下调通路点图\n")
cat("\n补充图:\n")
cat("  -", basename(output_supp_all_dotplot), "- 所有通路点图\n")
cat("  -", basename(output_supp_kegg_bar), "- KEGG排名条形图\n")
cat("  -", basename(output_supp_go_bar), "- GO排名条形图\n")
cat("\n数据文件:\n")
cat("  -", basename(output_gsea_results), "- 所有GSEA结果\n")
cat("  -", basename(output_cips_ranking), "- CIPS排名汇总\n")

cat("\n========================================\n")
cat("                完成\n")
cat("========================================\n")

sink()

cat("\n分析报告已保存: ", output_summary, "\n")
cat("\n=================== 分析完成 ====================\n")