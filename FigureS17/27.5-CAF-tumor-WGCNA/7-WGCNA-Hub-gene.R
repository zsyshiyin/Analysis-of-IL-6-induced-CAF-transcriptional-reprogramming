# ============================================================================
# WGCNA分析 - 模块5：Hub基因识别（完整版）
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
hub_dir <- file.path(output_dir, "hub_gene_analysis")
if (!dir.exists(hub_dir)) {
  dir.create(hub_dir, recursive = TRUE)
}

# ========================== 参数设置 =========================================
# Hub基因选择参数
top_n_hub_genes <- 30              # 每个模块选择的Top Hub基因数
top_n_network_genes <- 100         # 网络图中显示的Top基因数
correlation_threshold <- 0.5       # 相关性阈值
pvalue_threshold <- 0.05           # 显著性阈值

# 可视化参数
network_width <- 14
network_height <- 12
heatmap_width <- 12
heatmap_height <- 10

# 输出文件
output_rdata <- file.path(hub_dir, "hub_gene_analysis.RData")
output_hub_genes <- file.path(hub_dir, "hub_genes.csv")
output_module_network <- file.path(hub_dir, "module_network.graphml")
output_connectivity_plot <- file.path(hub_dir, "gene_connectivity.pdf")
output_hub_heatmap <- file.path(hub_dir, "hub_genes_heatmap.pdf")
output_hub_correlation <- file.path(hub_dir, "hub_genes_correlation.pdf")
output_summary <- file.path(hub_dir, "hub_gene_summary.txt")

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

library(WGCNA)
library(igraph)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(reshape2)

# 允许WGCNA使用多线程
enableWGCNAThreads()

cat("✓ 所有包已加载\n")

# ========================== 加载数据 ====================================
cat("\n========== 加载数据 ==========\n")

# 加载WGCNA结果
load(file.path(output_dir, "02_network_modules.RData"))
load(file.path(output_dir, "01_data_prepared.RData"))
load(file.path(output_dir, "03_module_enrichment.RData"))

cat(sprintf("Top模块: %s\n", paste(top_modules, collapse = ", ")))
cat(sprintf("表达矩阵维度: %d 基因 × %d 样品\n", nrow(expr_matrix_norm), ncol(expr_matrix_norm)))

# ========================== 计算基因连接度 ====================================
cat("\n========== 计算基因连接度 ==========\n")

# 计算每个基因的模块内连接度
intramodular_connectivity <- intramodularConnectivity(
  adj = adjacency,
  colors = moduleColors,
  scaleByMax = TRUE
)

# 添加基因名
intramodular_connectivity$Gene <- rownames(expr_matrix_norm)

cat("连接度计算完成\n")

# ========================== 计算模块隶属度（修复） ====================================
cat("\n========== 计算模块隶属度 ==========\n")

# 注意：expr_matrix_norm是基因×样品，MEs是样品×模块
# 需要转置expr_matrix_norm使其成为样品×基因
expr_matrix_t <- t(expr_matrix_norm)

# 计算每个基因与每个模块特征基因的相关性
# 结果矩阵：基因 × 模块
module_membership <- cor(expr_matrix_t, MEs, use = "pairwise.complete.obs")
colnames(module_membership) <- gsub("ME", "", colnames(module_membership))

cat(sprintf("模块隶属度矩阵维度: %d 基因 × %d 模块\n", 
            nrow(module_membership), ncol(module_membership)))

# ========================== 计算基因显著性 ====================================
cat("\n========== 计算基因显著性 ==========\n")

# 基因显著性 = 与模块特征基因相关性的绝对值
gene_significance <- abs(module_membership)

cat("基因显著性计算完成\n")

# ========================== 提取Top模块的Hub基因 ====================================
cat("\n========== 提取Top模块的Hub基因 ==========\n")

hub_genes_list <- list()
all_hub_genes <- data.frame()

for (mod in top_modules) {
  cat(sprintf("\n处理模块: %s\n", mod))
  
  # 获取该模块的基因
  mod_genes <- rownames(expr_matrix_norm)[moduleColors == mod]
  cat(sprintf("  模块基因数: %d\n", length(mod_genes)))
  
  # 获取该模块的连接度
  kWithin_col <- paste0("kWithin_", mod)
  
  if (kWithin_col %in% colnames(intramodular_connectivity)) {
    # 构建连接度数据框
    mod_connectivity <- data.frame(
      Gene = mod_genes,
      stringsAsFactors = FALSE
    )
    
    # 添加各项指标
    for (i in 1:nrow(mod_connectivity)) {
      gene <- mod_connectivity$Gene[i]
      
      # 模块内连接度
      mod_connectivity$kWithin[i] <- intramodular_connectivity[
        intramodular_connectivity$Gene == gene, kWithin_col]
      
      # 总连接度
      mod_connectivity$kTotal[i] <- intramodular_connectivity[
        intramodular_connectivity$Gene == gene, "kTotal"]
      
      # 模块外连接度
      kOut_col <- paste0("kOut_", mod)
      if (kOut_col %in% colnames(intramodular_connectivity)) {
        mod_connectivity$kOut[i] <- intramodular_connectivity[
          intramodular_connectivity$Gene == gene, kOut_col]
      }
      
      # 模块隶属度
      if (gene %in% rownames(module_membership)) {
        mod_connectivity$ModuleMembership[i] <- module_membership[gene, mod]
        mod_connectivity$GeneSignificance[i] <- abs(module_membership[gene, mod])
      }
    }
    
    # 按连接度排序
    mod_connectivity <- mod_connectivity %>%
      arrange(desc(kWithin)) %>%
      mutate(Rank = 1:n())
    
    # 选择Top Hub基因
    top_hubs <- head(mod_connectivity, min(top_n_hub_genes, nrow(mod_connectivity)))
    top_hubs$Module <- mod
    
    hub_genes_list[[mod]] <- top_hubs
    all_hub_genes <- rbind(all_hub_genes, top_hubs)
    
    cat(sprintf("  选择Top %d Hub基因\n", nrow(top_hubs)))
    cat("  Top 5 Hub基因:\n")
    for (i in 1:min(5, nrow(top_hubs))) {
      cat(sprintf("    %d. %s (kWithin=%.3f, MM=%.3f)\n", 
                  i, top_hubs$Gene[i], top_hubs$kWithin[i], 
                  top_hubs$ModuleMembership[i]))
    }
  }
}

# 保存Hub基因列表
if (nrow(all_hub_genes) > 0) {
  write.csv(all_hub_genes, output_hub_genes, row.names = FALSE)
  cat(sprintf("\n✓ Hub基因列表已保存: %s (%d 个基因)\n", 
              output_hub_genes, nrow(all_hub_genes)))
} else {
  cat("警告：未找到Hub基因\n")
}

# ========================== 计算Hub基因的表达相关性 ====================================
cat("\n========== 计算Hub基因的表达相关性 ==========\n")

if (nrow(all_hub_genes) > 0) {
  # 提取所有Hub基因
  hub_genes <- unique(all_hub_genes$Gene)
  
  # 提取Hub基因的表达矩阵
  hub_expr_matrix <- expr_matrix_norm[hub_genes, , drop = FALSE]
  
  # 计算Hub基因之间的相关性
  if (nrow(hub_expr_matrix) > 1) {
    hub_correlation <- cor(t(hub_expr_matrix), method = "spearman")
  } else {
    hub_correlation <- matrix(1, nrow = 1, ncol = 1)
    rownames(hub_correlation) <- hub_genes
    colnames(hub_correlation) <- hub_genes
  }
  
  cat(sprintf("Hub基因相关性计算完成 (%d × %d)\n", 
              nrow(hub_correlation), ncol(hub_correlation)))
}

# ========================== 图1：基因连接度分布图 ====================================
cat("\n========== 绘制基因连接度分布图 ==========\n")

# 准备连接度数据
connectivity_data <- data.frame(
  Gene = rownames(expr_matrix_norm),
  kTotal = intramodular_connectivity$kTotal,
  Module = moduleColors,
  stringsAsFactors = FALSE
)

# 为Top模块标记颜色
connectivity_data$IsTopModule <- connectivity_data$Module %in% top_modules

p1 <- ggplot(connectivity_data, aes(x = kTotal, fill = IsTopModule)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("TRUE" = "#C60036", "FALSE" = "#2C7BB6"),
                    labels = c("TRUE" = "Top Modules", "FALSE" = "Other Modules")) +
  labs(title = "Gene Connectivity Distribution",
       x = "Total Connectivity (kTotal)",
       y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top")

ggsave(output_connectivity_plot, p1, width = 8, height = 6, dpi = 300)
cat("✓ 基因连接度分布图已保存\n")

# ========================== 图2：Hub基因表达热图 ====================================
cat("\n========== 绘制Hub基因表达热图 ==========\n")

if (nrow(all_hub_genes) > 0 && nrow(hub_expr_matrix) > 0) {
  # Z-score标准化
  hub_expr_scaled <- t(scale(t(hub_expr_matrix)))
  hub_expr_scaled[is.na(hub_expr_scaled)] <- 0
  hub_expr_scaled[hub_expr_scaled > 3] <- 3
  hub_expr_scaled[hub_expr_scaled < -3] <- -3
  
  # 创建注释
  annotation_col <- data.frame(
    Group = sample_groups$Group,
    row.names = sample_groups$SampleID
  )
  
  # 为Hub基因添加模块注释
  annotation_row <- data.frame(
    Module = all_hub_genes$Module[match(rownames(hub_expr_scaled), all_hub_genes$Gene)],
    row.names = rownames(hub_expr_scaled)
  )
  
  # 颜色设置
  group_colors <- c("Control" = "#2E86AB", "DMEM" = "#1F8A4C", 
                    "I6" = "#C60036", "R2" = "#E4945A")
  
  # 为每个模块分配颜色
  unique_modules <- unique(all_hub_genes$Module)
  module_color_palette <- c("blue" = "#377EB8", "brown" = "#8B4513", 
                            "turquoise" = "#40E0D0", "yellow" = "#FFD700", 
                            "green" = "#4DAF4A", "red" = "#D7191C",
                            "black" = "#000000", "pink" = "#FF69B4",
                            "magenta" = "#FF00FF", "purple" = "#800080")
  module_colors_used <- module_color_palette[names(module_color_palette) %in% unique_modules]
  
  annotation_colors <- list(
    Group = group_colors,
    Module = module_colors_used
  )
  
  pdf(output_hub_heatmap, width = heatmap_width, height = heatmap_height)
  pheatmap(hub_expr_scaled,
           color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           annotation_row = annotation_row,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 6,
           fontsize_col = 8,
           main = paste("Hub Genes Expression Heatmap (", nrow(hub_expr_scaled), " genes)", sep = ""),
           border_color = NA)
  dev.off()
  cat("✓ Hub基因表达热图已保存\n")
}

# ========================== 图3：Hub基因相关性热图 ====================================
cat("\n========== 绘制Hub基因相关性热图 ==========\n")

if (exists("hub_correlation") && nrow(hub_correlation) > 1 && ncol(hub_correlation) > 1) {
  pdf(output_hub_correlation, width = heatmap_width, height = heatmap_height)
  pheatmap(hub_correlation,
           color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
           breaks = seq(-1, 1, length.out = 101),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = "Hub Genes Correlation Heatmap",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = FALSE)
  dev.off()
  cat("✓ Hub基因相关性热图已保存\n")
}

# ========================== 图4：模块网络图 ====================================
cat("\n========== 绘制模块网络图 ==========\n")

if (nrow(all_hub_genes) > 0) {
  for (mod in top_modules) {
    cat(sprintf("\n构建模块 %s 的网络图...\n", mod))
    
    # 获取该模块的Hub基因
    mod_hubs <- all_hub_genes[all_hub_genes$Module == mod, ]
    top_network_genes <- head(mod_hubs$Gene, min(top_n_network_genes, nrow(mod_hubs)))
    
    if (length(top_network_genes) >= 3) {
      # 提取这些基因的表达数据
      network_expr <- expr_matrix_norm[top_network_genes, , drop = FALSE]
      
      # 计算相关性矩阵
      if (nrow(network_expr) > 1) {
        cor_matrix <- cor(t(network_expr), method = "spearman")
        
        # 筛选显著相关的边
        edge_list <- data.frame()
        for (i in 1:(nrow(cor_matrix)-1)) {
          for (j in (i+1):nrow(cor_matrix)) {
            if (abs(cor_matrix[i, j]) >= correlation_threshold) {
              edge_list <- rbind(edge_list, data.frame(
                from = rownames(cor_matrix)[i],
                to = colnames(cor_matrix)[j],
                weight = cor_matrix[i, j],
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        if (nrow(edge_list) > 0) {
          # 创建网络图
          g <- graph_from_data_frame(edge_list, directed = FALSE)
          
          # 设置节点大小（基于连接度）
          degree_centrality <- degree(g)
          V(g)$size <- 5 + (degree_centrality / max(degree_centrality)) * 10
          
          # 设置节点颜色（基于模块内连接度）
          node_connectivity <- mod_hubs$kWithin[match(V(g)$name, mod_hubs$Gene)]
          if (length(unique(node_connectivity)) > 1) {
            V(g)$color <- colorRampPalette(c("#2C7BB6", "#D7191C"))(100)[
              as.numeric(cut(node_connectivity, breaks = 100))]
          } else {
            V(g)$color <- "#D7191C"
          }
          
          # 设置边的宽度和颜色
          E(g)$width <- abs(E(g)$weight) * 3
          E(g)$color <- ifelse(E(g)$weight > 0, "#D7191C", "#2C7BB6")
          
          # 保存为.graphml格式
          network_file <- file.path(hub_dir, paste0("module_", mod, "_network.graphml"))
          write.graph(g, file = network_file, format = "graphml")
          
          # 绘制网络图
          network_plot_file <- file.path(hub_dir, paste0("module_", mod, "_network.pdf"))
          pdf(network_plot_file, width = network_width, height = network_height)
          plot(g,
               layout = layout_with_fr(g),
               vertex.color = V(g)$color,
               vertex.size = V(g)$size,
               vertex.label = V(g)$name,
               vertex.label.cex = 0.7,
               vertex.label.dist = 0.5,
               edge.width = E(g)$width,
               edge.color = E(g)$color,
               main = paste("Module", mod, "Hub Gene Network"))
          legend("topright",
                 legend = c("Positive Correlation", "Negative Correlation"),
                 col = c("#D7191C", "#2C7BB6"),
                 lty = 1,
                 cex = 0.7)
          dev.off()
          cat(sprintf("  ✓ 模块 %s 网络图已保存 (%d 个节点, %d 条边)\n", 
                      mod, vcount(g), ecount(g)))
        }
      }
    }
  }
}

# ========================== 生成汇总报告 ====================================
cat("\n========== 生成汇总报告 ==========\n")

sink(output_summary)

cat("================================================================================\n")
cat("                      Hub基因识别分析报告\n")
cat("================================================================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、分析参数\n")
cat("--------------------------------------------------------------------------------\n")
cat(sprintf("  Top模块: %s\n", paste(top_modules, collapse = ", ")))
cat(sprintf("  每个模块Hub基因数: %d\n", top_n_hub_genes))
cat(sprintf("  相关性阈值: %.2f\n", correlation_threshold))

cat("\n二、Hub基因统计\n")
cat("--------------------------------------------------------------------------------\n")
for (mod in top_modules) {
  mod_hubs <- all_hub_genes[all_hub_genes$Module == mod, ]
  if (nrow(mod_hubs) > 0) {
    cat(sprintf("\n模块 %s:\n", mod))
    cat(sprintf("  Hub基因数: %d\n", nrow(mod_hubs)))
    cat("  Top 10 Hub基因:\n")
    for (i in 1:min(10, nrow(mod_hubs))) {
      cat(sprintf("    %d. %s (kWithin=%.3f, MM=%.3f)\n", 
                  i, mod_hubs$Gene[i], mod_hubs$kWithin[i], mod_hubs$ModuleMembership[i]))
    }
  }
}

cat("\n三、输出文件\n")
cat("--------------------------------------------------------------------------------\n")
cat("1.", basename(output_hub_genes), "- Hub基因列表\n")
cat("2.", basename(output_connectivity_plot), "- 基因连接度分布图\n")
cat("3.", basename(output_hub_heatmap), "- Hub基因表达热图\n")
if (exists("hub_correlation") && nrow(hub_correlation) > 1) {
  cat("4.", basename(output_hub_correlation), "- Hub基因相关性热图\n")
}
cat("5. module_*_network.pdf - 各模块网络图\n")
cat("6. module_*_network.graphml - 各模块网络文件（Cytoscape格式）\n")

sink()

cat("\n========================================\n")
cat("Hub基因识别完成！\n")
cat("输出目录:", hub_dir, "\n")
cat("========================================\n")