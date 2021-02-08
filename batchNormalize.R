#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")
#biocLite("limma")

library(sva)
library(limma)
setwd("C:\\Users\\lexb4\\Desktop\\geoBatch\\06.batchNormalize")
rt=read.table("merge.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType=c(rep(1,50),rep(2,10))
modType=c(rep("normal",25),rep("tumor",25),rep("normal",5),rep("tumor",5))
mod = model.matrix(~as.factor(modType))
outTab=ComBat(data, batchType, mod, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)



###Video source: http://study.163.com/u/biowolf
######Video source: https://shop119322454.taobao.com
######速科生物: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio
######QQ群:  259208034
