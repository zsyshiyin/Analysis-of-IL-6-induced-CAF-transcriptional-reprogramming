setwd("D:\\zsy\\SX\\Fomal SX\\diff-CAF-NF\\Diff-FAP")         # 设置工作目录
inputFile <- "singleGene.txt"                          # 输入文件名
yMin <- 0                                             # y轴最小值
yMax <- 15000                                         # y轴最大值
ySeg <- yMax * 0.94

library(limma)
library(beeswarm)

rt <- read.table(inputFile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
geneName <- colnames(rt)[1]

# 修改：交换组别顺序（改为"NF", "CAF"）
labels <- c("NF", "CAF")  # 交换组别顺序，现在NF在前，CAF在后
colnames(rt) <- c("expression", "Type")

# 为了确保正确的因子顺序，重新设置因子水平
rt$Type <- factor(rt$Type, levels = c("NF", "CAF"))

# 离散值处理方法：将所有超出阈值的值替换为下四分位数（Q1）
replace_with_q1 <- function(x) {
  # 使用标准阈值识别离散值
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  H <- 0.5 * IQR(x, na.rm = TRUE)  # 标准阈值：1.5倍IQR
  
  # 识别离散值
  lower_bound <- qnt[1] - H
  upper_bound <- qnt[2] + H
  outliers <- x < lower_bound | x > upper_bound
  
  # 将离散值替换为下四分位数（Q1）
  if(any(outliers)) {
    q1 <- quantile(x, probs = 0.25, na.rm = TRUE)
    x[outliers] <- q1
  }
  return(x)
}

# 应用离散值处理
rt_clean <- rt
rt_clean$expression <- ave(rt$expression, rt$Type, FUN = replace_with_q1)

# 计算统计检验
wilcoxTest <- wilcox.test(expression ~ Type, data = rt_clean)
wilcoxP <- wilcoxTest$p.value
pvalue <- signif(wilcoxP, 4)
pval <- ifelse(pvalue < 0.001, 
               format(signif(pvalue, 4), scientific = TRUE), 
               round(pvalue, 3))

# 计算各组中位数用于显示
medians <- aggregate(expression ~ Type, rt_clean, median)
# 确保中位数按照labels的顺序排列
medians <- medians[match(labels, medians$Type), ]

outFile <- paste0(geneName, ".pdf")
pdf(file = outFile, width = 10, height = 10)

# 调整边距以适应更大的字体
par(mar = c(5, 8, 4, 3), mgp = c(4, 1, 0))

# 修改：交换颜色顺序（现在NF用蓝色，CAF用红色）
col_fill <- c("#25377F4D", "#CA0E124D")  # 30%透明度，NF蓝色，CAF红色
col_border <- c("#25377F", "#CA0E12")    # NF深蓝色，CAF深红色

# 绘制箱线图 - 缩短箱型
boxplot(expression ~ Type, data = rt_clean, names = labels,
        ylab = paste0(geneName, " expression"),
        cex.main = 3,        # 加大标题字体
        cex.lab = 3,         # 加大坐标轴标签字体
        cex.axis = 2.5,        # 加大坐标轴刻度字体
        ylim = c(yMin, yMax),
        outline = FALSE, 
        col = col_fill, 
        border = col_border,
        boxwex = 0.38,          # 缩短箱型宽度
        staplewex = 0.3,       # 缩短须线宽度
        whisklty = 1, 
        whiskcol = col_border,
        lwd = 2)               # 加粗箱线图边框

# 添加蜂群点图（使用交换后的颜色）
beeswarm(expression ~ Type, data = rt_clean, 
         col = col_border, lwd = 0.1, pch = 16, 
         add = TRUE, corral = "wrap", cex = 1.2)

# 添加p值标注（位置不变）
segments(1, ySeg, 2, ySeg, lwd = 2)
segments(1, ySeg, 1, ySeg * 0.96, lwd = 2)
segments(2, ySeg, 2, ySeg * 0.96, lwd = 2)
text(1.5, ySeg * 1.05, 
     labels = paste0("p = ", pval), 
     cex = 2.2,                # 加大p值字体
     font = 1)                 # 加粗字体

# 中位数标注 - 优化后的位置
text(x = 1:2 + 0.18,          # 向右移动一些
     y = medians$expression,  
     labels = sprintf("%.2f", medians$expression), 
     pos = 4,                 # 文字位于指定位置的右侧
     cex = 2.2,               # 加大中位数标记字体
     font = 1,                # 加粗字体
     col = "black",
     xpd = TRUE)

dev.off()

# 输出统计信息
cat("\n========== 离散值替换为下四分位数(Q1)的结果 ==========\n")
cat("基因:", geneName, "\n")
cat("处理前样本数:", nrow(rt), "\n")
cat("处理后样本数:", nrow(rt_clean), "(所有样本保留)\n")
cat("处理方式: 使用1.5倍IQR识别离散值，替换为该组的下四分位数(Q1)\n\n")

# 统计每组的离散值数量和处理效果
for(type in levels(rt$Type)) {
  type_data_orig <- rt$expression[rt$Type == type]
  type_data_clean <- rt_clean$expression[rt_clean$Type == type]
  
  qnt <- quantile(type_data_orig, probs = c(0.25, 0.75))
  H <- 0.5 * IQR(type_data_orig)
  lower_bound <- qnt[1] - H
  upper_bound <- qnt[2] + H
  
  n_lower_outliers <- sum(type_data_orig < lower_bound)
  n_upper_outliers <- sum(type_data_orig > upper_bound)
  n_outliers <- n_lower_outliers + n_upper_outliers
  
  q1 <- quantile(type_data_orig, 0.25)
  q3 <- quantile(type_data_orig, 0.75)
  
  cat(sprintf("\n%s组:\n", type))
  cat(sprintf("  总样本数: %d\n", length(type_data_orig)))
  cat(sprintf("  离散值数量: %d (%.1f%%)\n", n_outliers, 100 * n_outliers/length(type_data_orig)))
  cat(sprintf("  原始中位数: %.2f\n", median(type_data_orig)))
  cat(sprintf("  处理后中位数: %.2f\n", median(type_data_clean)))
  cat(sprintf("  下四分位数(Q1): %.2f\n", q1))
  cat(sprintf("  上四分位数(Q3): %.2f\n", q3))
  cat(sprintf("  原始范围: %.2f - %.2f\n", min(type_data_orig), max(type_data_orig)))
  cat(sprintf("  处理后范围: %.2f - %.2f\n", min(type_data_clean), max(type_data_clean)))
}

cat("\n处理后的p值:", pval, "\n")
cat("==================================================\n")

# 保存处理后的数据
write.csv(rt_clean, paste0(geneName, "_outliers_to_q1.csv"), row.names = TRUE)