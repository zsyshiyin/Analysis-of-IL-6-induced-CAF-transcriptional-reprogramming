setwd("D:\\zsy\\SX\\DJY SX\\Gene-co-expression\\IL6-VEGFB")  #dir
rt=read.table("data.txt", sep="\t", header=T, check.names=F, row.names=1)
i=1  #gene 1
j=2  #gene 2 
x=as.numeric(rt[i,])
y=as.numeric(rt[j,])
xName=row.names(rt[i,])
yName=row.names(rt[j,])

# 计算Spearman相关系数和p值
cor_test_result = cor.test(x, y, method = "spearman")
R = round(cor_test_result$estimate, 3)
p_value = round(cor_test_result$p.value, 4)

# 设置更宽的PDF（默认width=7, height=7）
pdf(file="cor.pdf", width=7, height=7)  # 宽度10英寸，高度7英寸
z=lm(y~x)
plot(x, y, 
     type="p", pch=16, cex=1.5, col="#25377F",
     main=paste("Spearman's correlation = ", R, "\np-value = ", p_value, sep=""),
     xlab=paste(xName, "mRNA expression"), 
     ylab=paste(yName, "mRNA expression"),
     cex.main=1.8,     # 加大标题（Spearman和p值）字体
     cex.lab=1.8,      # 加大X/Y轴标签字体
     cex.axis=1.5,     # 加大坐标轴刻度字体
     mgp=c(2.5, 0.8, 0)) # 调整坐标轴标题位置（加大距离）
lines(x, fitted(z), col="#CA0E12", lwd=2)  # 加粗回归线
dev.off()