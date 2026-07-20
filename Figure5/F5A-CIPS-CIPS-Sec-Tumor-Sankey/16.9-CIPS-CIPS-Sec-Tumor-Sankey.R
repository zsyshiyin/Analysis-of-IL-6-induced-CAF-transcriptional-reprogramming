# ============================================================================
# 整合桑基图：上游TF → 分泌蛋白配体 → 受体 → 下游肿瘤TF
# 严格使用输入文件数据（无演示数据）
# ============================================================================

rm(list = ls())

# ========================== 工作路径设定 ====================================
work_dir <- "D:/zsy/SX/Fomal-final/16.9-CIPS-CIPS-Sec-Tumor-Sankey"
setwd(work_dir)
cat("当前工作路径:", getwd(), "\n")

# ========================== 创建输出目录 ====================================
output_dir <- "Integrated_Sankey_Visualization"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# ============================================================================
# 输入文件（请根据实际情况修改路径）【必须存在】
# ============================================================================
SECRETED_FILE <- "secreted_proteins.csv"              # 分泌蛋白列表
LIGAND_RECEPTOR_TF_FILE <- "ligand_receptor_tf.csv"   # 配体→受体→TF信号链
UPSTREAM_TF_FILE <- "upstream_tf_targets.csv"       # 上游TF→分泌蛋白调控关系

# ============================================================================
# 可视化参数（可自定义）
# ============================================================================
sankey_font_size <- 20
sankey_node_width <- 18
sankey_node_padding <- 5
sankey_height <- 1900
sankey_width <- 900
sankey_link_thickness <- 0.01       # 连线粗细缩放系数（越小越细，推荐0.3-0.7）

# 节点颜色
color_upstream_tf <- "#C60036"       # 上游调控TF（橙色）
color_ligand <- "#F39B7F"            # 配体（深红）
color_receptor <- "#2C7BB6"          # 受体（蓝色）
color_downstream_tf <- "#2AA7DE"     # 下游肿瘤活化TF（亮红）
color_placeholder <- "#E0E0E0"       # 占位节点（浅灰色）

# 筛选参数
TOP_UPSTREAM_TFS <- 15
TOP_LIGANDS <- 20
TOP_DOWNSTREAM_TFS <- 15

# ========================== 加载必要的包 ====================================
cat("\n========== 加载必要的包 ==========\n")

required_packages <- c("dplyr", "tidyr", "readr", "networkD3", "htmlwidgets")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
  cat(sprintf("✓ %s\n", pkg))
}

# ========================== 检查输入文件 ====================================
cat("\n========== 检查输入文件 ==========\n")

files_to_check <- c(SECRETED_FILE, LIGAND_RECEPTOR_TF_FILE, UPSTREAM_TF_FILE)
missing_files <- !file.exists(files_to_check)

if (any(missing_files)) {
  cat("\n错误：以下文件不存在:\n")
  for (f in files_to_check[missing_files]) cat("  -", f, "\n")
  stop("请修改输入文件路径")
}

# ========================== 读取数据 ====================================
cat("\n========== 读取数据 ==========\n")

secreted_data <- read.csv(SECRETED_FILE, stringsAsFactors = FALSE)
secreted_proteins <- unique(secreted_data[[1]])
cat(sprintf("✓ 分泌蛋白: %d 个\n", length(secreted_proteins)))

lr_tf_data <- read.csv(LIGAND_RECEPTOR_TF_FILE, stringsAsFactors = FALSE)
cat(sprintf("✓ 配体→受体→TF: %d 行\n", nrow(lr_tf_data)))

upstream_data <- read.csv(UPSTREAM_TF_FILE, stringsAsFactors = FALSE)
cat(sprintf("✓ 上游TF→靶基因: %d 行\n", nrow(upstream_data)))

# ========================== 名称清洗 ====================================
clean_gene_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x <- gsub("\\s+", "", x)
  return(x)
}

# ========================== 数据整合 ====================================
cat("\n========== 整合数据 ==========\n")

detect_columns <- function(data, possible_names) {
  col_names <- colnames(data)
  matched <- intersect(col_names, possible_names)
  if (length(matched) == 0) stop("无法识别列名: ", paste(col_names, collapse = ", "))
  return(matched[1])
}

# 上游数据
tf_col <- detect_columns(upstream_data, c("tf", "TF", "source", "from"))
target_col <- detect_columns(upstream_data, c("targets", "Targets", "target", "Target", "Target_Genes"))

upstream_expanded <- upstream_data %>%
  dplyr::select(upstream_tf = !!sym(tf_col), targets = !!sym(target_col)) %>%
  filter(!is.na(targets) & targets != "") %>%
  mutate(target_list = strsplit(as.character(targets), ";\\s*")) %>%
  tidyr::unnest(target_list) %>%
  mutate(
    ligand_raw = trimws(target_list),
    ligand = clean_gene_name(ligand_raw)
  ) %>%
  filter(ligand != "") %>%
  dplyr::select(upstream_tf, ligand, ligand_raw) %>%
  distinct()

# 下游数据
ligand_col <- detect_columns(lr_tf_data, c("ligand", "Ligand", "source", "from"))
receptor_col <- detect_columns(lr_tf_data, c("receptor", "Receptor", "target"))
tf_lr_col <- detect_columns(lr_tf_data, c("tf", "TF", "transcription_factor", "downstream_tf"))

lr_expanded <- lr_tf_data %>%
  dplyr::rename(
    ligand_raw = !!sym(ligand_col),
    receptor_raw = !!sym(receptor_col),
    downstream_tf_raw = !!sym(tf_lr_col)
  ) %>%
  filter(!is.na(ligand_raw) & !is.na(receptor_raw) & !is.na(downstream_tf_raw)) %>%
  mutate(
    ligand = clean_gene_name(ligand_raw),
    receptor = clean_gene_name(receptor_raw),
    downstream_tf = clean_gene_name(downstream_tf_raw)
  ) %>%
  distinct()

# 配体显示名称映射
ligand_display_map <- bind_rows(
  upstream_expanded %>% dplyr::select(ligand, ligand_raw),
  lr_expanded %>% dplyr::select(ligand, ligand_raw)
) %>% distinct(ligand, .keep_all = TRUE)

# ========================== 筛选节点 ====================================
cat("\n========== 筛选节点 ==========\n")

# 上游TF
upstream_degree <- upstream_expanded %>%
  group_by(upstream_tf) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))
top_upstream_tfs <- head(upstream_degree, TOP_UPSTREAM_TFS)$upstream_tf

# 下游TF
downstream_degree <- lr_expanded %>%
  group_by(downstream_tf) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))
top_downstream_tfs <- head(downstream_degree, TOP_DOWNSTREAM_TFS)$downstream_tf

# 配体（取并集）
upstream_ligands <- upstream_expanded %>%
  filter(upstream_tf %in% top_upstream_tfs) %>%
  pull(ligand) %>% unique()

downstream_ligands <- lr_expanded %>%
  filter(downstream_tf %in% top_downstream_tfs) %>%
  pull(ligand) %>% unique()

all_ligands <- union(upstream_ligands, downstream_ligands)

# 按连接数排序
ligand_degree <- bind_rows(
  upstream_expanded %>% filter(ligand %in% all_ligands) %>% group_by(ligand) %>% summarise(n = n(), .groups = 'drop'),
  lr_expanded %>% filter(ligand %in% all_ligands) %>% group_by(ligand) %>% summarise(n = n(), .groups = 'drop')
) %>%
  group_by(ligand) %>%
  summarise(total = sum(n), .groups = 'drop') %>%
  arrange(desc(total))

top_ligands <- head(ligand_degree, TOP_LIGANDS)$ligand
top_ligands_display <- ligand_display_map$ligand_raw[match(top_ligands, ligand_display_map$ligand)]
top_ligands_display <- ifelse(is.na(top_ligands_display), top_ligands, top_ligands_display)

cat(sprintf("上游TF: %d, 配体: %d, 下游TF: %d\n", 
            length(top_upstream_tfs), length(top_ligands), length(top_downstream_tfs)))

# ========================== 构建节点（关键：用前缀区分角色） ====================================
cat("\n========== 构建节点 ==========\n")

# 占位节点名称
PLACEHOLDER_UPSTREAM <- "Unknown TF"
PLACEHOLDER_RECEPTOR <- "Unknown Receptor"
PLACEHOLDER_DOWNSTREAM <- "Unknown TF"

nodes <- data.frame(
  name = character(),
  display_name = character(),
  type = character(),
  stringsAsFactors = FALSE
)

# 第1列：上游TF（前缀 UP_）
upstream_nodes <- data.frame(
  name = paste0("UP_", top_upstream_tfs),
  display_name = top_upstream_tfs,
  type = "Upstream_TF",
  stringsAsFactors = FALSE
)
upstream_nodes <- rbind(upstream_nodes, data.frame(
  name = "UP_PLACEHOLDER",
  display_name = PLACEHOLDER_UPSTREAM,
  type = "Placeholder",
  stringsAsFactors = FALSE
))
nodes <- bind_rows(nodes, upstream_nodes)

# 第2列：配体（前缀 LIG_）
ligand_nodes <- data.frame(
  name = paste0("LIG_", top_ligands),
  display_name = top_ligands_display,
  type = "Ligand",
  stringsAsFactors = FALSE
)
nodes <- bind_rows(nodes, ligand_nodes)

# 第3列：受体（前缀 REC_）
receptors <- lr_expanded %>%
  filter(ligand %in% top_ligands) %>%
  pull(receptor) %>% unique()
receptor_nodes <- data.frame(
  name = paste0("REC_", receptors),
  display_name = receptors,
  type = "Receptor",
  stringsAsFactors = FALSE
)
receptor_nodes <- rbind(receptor_nodes, data.frame(
  name = "REC_PLACEHOLDER",
  display_name = PLACEHOLDER_RECEPTOR,
  type = "Placeholder",
  stringsAsFactors = FALSE
))
nodes <- bind_rows(nodes, receptor_nodes)

# 第4列：下游TF（前缀 DOWN_）
downstream_tf_list <- lr_expanded %>%
  filter(ligand %in% top_ligands, downstream_tf %in% top_downstream_tfs) %>%
  pull(downstream_tf) %>% unique()
downstream_nodes <- data.frame(
  name = paste0("DOWN_", downstream_tf_list),
  display_name = downstream_tf_list,
  type = "Downstream_TF",
  stringsAsFactors = FALSE
)
downstream_nodes <- rbind(downstream_nodes, data.frame(
  name = "DOWN_PLACEHOLDER",
  display_name = PLACEHOLDER_DOWNSTREAM,
  type = "Placeholder",
  stringsAsFactors = FALSE
))
nodes <- bind_rows(nodes, downstream_nodes)

# 添加索引和颜色
nodes <- nodes %>%
  mutate(
    id = 0:(n() - 1),
    color_group = case_when(
      type == "Upstream_TF" ~ "Upstream_TF",
      type == "Ligand" ~ "Ligand",
      type == "Receptor" ~ "Receptor",
      type == "Downstream_TF" ~ "Downstream_TF",
      type == "Placeholder" ~ "Placeholder",
      TRUE ~ "Other"
    )
  )

cat(sprintf("总节点数: %d\n", nrow(nodes)))
print(table(nodes$type))

# 创建映射表
name_to_id <- setNames(nodes$id, nodes$name)

# ========================== 构建链接（使用带前缀的节点名） ====================================
cat("\n========== 构建链接 ==========\n")

links <- data.frame()

# 为每个展示的配体构建连接
for (lig in top_ligands) {
  lig_node <- paste0("LIG_", lig)
  
  # 1. 上游TF → 配体
  ut_list <- upstream_expanded %>%
    filter(ligand == lig, upstream_tf %in% top_upstream_tfs) %>%
    pull(upstream_tf)
  
  if (length(ut_list) > 0) {
    for (ut in ut_list) {
      links <- rbind(links, data.frame(
        source = name_to_id[paste0("UP_", ut)],
        target = name_to_id[lig_node],
        value = 1,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    links <- rbind(links, data.frame(
      source = name_to_id["UP_PLACEHOLDER"],
      target = name_to_id[lig_node],
      value = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # 2. 配体 → 受体
  rec_list <- lr_expanded %>%
    filter(ligand == lig) %>%
    pull(receptor) %>% unique()
  
  if (length(rec_list) > 0) {
    for (rec in rec_list) {
      rec_node <- paste0("REC_", rec)
      links <- rbind(links, data.frame(
        source = name_to_id[lig_node],
        target = name_to_id[rec_node],
        value = 1,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    links <- rbind(links, data.frame(
      source = name_to_id[lig_node],
      target = name_to_id["REC_PLACEHOLDER"],
      value = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  # 3. 受体 → 下游TF
  if (length(rec_list) > 0) {
    for (rec in rec_list) {
      rec_node <- paste0("REC_", rec)
      
      dt_list <- lr_expanded %>%
        filter(ligand == lig, receptor == rec, downstream_tf %in% top_downstream_tfs) %>%
        pull(downstream_tf) %>% unique()
      
      if (length(dt_list) > 0) {
        for (dt in dt_list) {
          links <- rbind(links, data.frame(
            source = name_to_id[rec_node],
            target = name_to_id[paste0("DOWN_", dt)],
            value = 1,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        links <- rbind(links, data.frame(
          source = name_to_id[rec_node],
          target = name_to_id["DOWN_PLACEHOLDER"],
          value = 1,
          stringsAsFactors = FALSE
        ))
      }
    }
  } else {
    links <- rbind(links, data.frame(
      source = name_to_id["REC_PLACEHOLDER"],
      target = name_to_id["DOWN_PLACEHOLDER"],
      value = 1,
      stringsAsFactors = FALSE
    ))
  }
}

# 聚合相同链接
links <- links %>%
  group_by(source, target) %>%
  summarise(value = sum(value), .groups = 'drop')

# 调整连线粗细
links$value <- links$value * sankey_link_thickness

cat(sprintf("总链接数: %d\n", nrow(links)))

# ========================== 创建桑基图 ====================================
cat("\n========== 生成桑基图 ==========\n")

color_scale <- paste0(
  'd3.scaleOrdinal()',
  '.domain(["Upstream_TF", "Ligand", "Receptor", "Downstream_TF", "Placeholder"])',
  '.range(["', color_upstream_tf, '", "', color_ligand, '", "', 
  color_receptor, '", "', color_downstream_tf, '", "', color_placeholder, '"])'
)

sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "display_name",
  NodeGroup = "color_group",
  colourScale = color_scale,
  fontSize = sankey_font_size,
  nodeWidth = sankey_node_width,
  nodePadding = sankey_node_padding,
  height = sankey_height,
  width = sankey_width,
  iterations = 0,
  sinksRight = TRUE
)

# 设置 Arial 字体
sankey_plot <- htmlwidgets::onRender(
  sankey_plot,
  '
  function(el, x) {
    d3.select(el).selectAll(".node text")
      .style("font-family", "Arial");
  }
  '
)

html_file <- file.path(output_dir, "Integrated_Sankey_With_Prefix.html")
saveWidget(sankey_plot, html_file, selfcontained = TRUE)
cat("✓ 桑基图已保存:", html_file, "\n")

# 保存数据
write.csv(nodes, file.path(output_dir, "Sankey_Nodes.csv"), row.names = FALSE)
write.csv(links, file.path(output_dir, "Sankey_Links.csv"), row.names = FALSE)

# ========================== 输出 PDF 和 PNG ====================================
cat("\n========== 导出 PDF 和 PNG ==========\n")

if (!requireNamespace("webshot", quietly = TRUE)) {
  install.packages("webshot")
  library(webshot)
  webshot::install_phantomjs()
}
library(webshot)

webshot(html_file, 
        file.path(output_dir, "Integrated_Sankey_Diagram.pdf"),
        vwidth = sankey_width, 
        vheight = sankey_height)

webshot(html_file, 
        file.path(output_dir, "Integrated_Sankey_Diagram.png"),
        vwidth = sankey_width, 
        vheight = sankey_height,
        zoom = 2)

cat("✓ PDF 和 PNG 已保存\n")

cat("\n完成！输出目录:", output_dir, "\n")


# ========================== 保存数据 ====================================
write.csv(nodes, file.path(output_dir, "Sankey_Nodes.csv"), row.names = FALSE)
write.csv(links, file.path(output_dir, "Sankey_Links.csv"), row.names = FALSE)

cat("\n完成！输出目录:", output_dir, "\n")