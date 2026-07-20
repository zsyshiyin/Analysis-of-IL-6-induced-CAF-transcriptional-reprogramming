######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("devtools")
#library(devtools)
#devtools::install_github("GfellerLab/EPIC")
#install_github('ebecht/MCPcounter',ref='master', subdir='Source')
#devtools::install_github('dviraran/xCell')

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#引用包
library(limma)
library(EPIC)
library(xCell)
library(MCPcounter)
library(estimate)

inputFile="mRNA.txt"     #表达数据文件
setwd("C:\\biowolf\\CAF\\12.tcgaCAF")     #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
rt=data[rowMeans(data)>0,]

#去除正常样品
group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
rt=rt[,group==0]
rt=t(rt)
row.names(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(rt))
rt=avereps(rt)
data=t(rt)

#输出整理后的表达数据文件
out=rbind(ID=colnames(data), data)
write.table(out,file="TCGA.symbol.txt",sep="\t",quote=F,col.names=F)

#输出TIDE输入文件
tide=log2(data+1)
tideMean=rowMeans(tide)
tide=tide-tideMean
tideOut=rbind(ID=colnames(tide), tide)
write.table(tideOut, file="TCGA.tide.txt", sep="\t", quote=F, col.names=F)

#运行EPIC
data=log2(data+1)
epic=EPIC(bulk=data)
CAF_EPIC=epic$cellFractions[,"CAFs"]
outTab=rbind(ID=colnames(epic$cellFractions), epic$cellFractions)
write.table(outTab, file="score.EPIC.txt", sep="\t", quote=F, col.names=F)

#运行MCPcounter
MCPcounter.estimate <- MCPcounter.estimate(data,
	featuresType="HUGO_symbols",
	probesets=read.table("MCPcounter.probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character"),
	genes=read.table("MCPcounter.genes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
)
CAF_MCPcounter=MCPcounter.estimate["Fibroblasts",]
outTab=rbind(ID=colnames(MCPcounter.estimate), MCPcounter.estimate)
write.table(outTab, file="score.MCPcounter.txt", sep="\t", quote=F, col.names=F)

#运行xCell
xCell=xCellAnalysis(data, rnaseq=TRUE, file.name="score.xCell.txt", parallel.sz=1)
CAF_xCell=xCell["Fibroblasts",]

#运行estimate包,得到肿瘤微环境的打分
filterCommonGenes(input.f="TCGA.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")
#对肿瘤微环境的打分进行整理, 输出每个样品的打分
scores=read.table("estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
StromalScore=scores[,"StromalScore"]
out=rbind(ID=colnames(scores), scores)
write.table(out, file="score.estimate.txt", sep="\t", quote=F, col.names=F)

#输出CAF打分结果
caf=cbind(CAF_EPIC, CAF_MCPcounter, CAF_xCell, StromalScore)
outTab=rbind(ID=colnames(caf), caf)
write.table(outTab, file="TCGA.CAF.txt", sep="\t", quote=F, col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

