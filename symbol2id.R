#install.packages("RSQLite")
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")

setwd("C:\\Users\\lexb4\\Desktop\\geoBatch\\09.symbol2id")

library("org.Hs.eg.db")
rt=read.table("symbol.txt",sep="\t",check.names=F,header=T)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)



###Video source: http://study.163.com/u/biowolf
######Video source: https://shop119322454.taobao.com
######�ٿ�����: http://www.biowolf.cn/
######�������䣺2740881706@qq.com
######����΢��: seqBio
######QQȺ:  259208034