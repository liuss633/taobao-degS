#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")

logFoldChange=1
adjustP=0.05

library(limma)
setwd("C:\\Users\\lexb4\\Desktop\\geoBatch\\07.diff")
rt=read.table("normalize.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
#rt=log2(rt+1)

#differential
modType=c(rep("normal",25),rep("tumor",25),rep("normal",5),rep("tumor",5))
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)

#write expression level of diff gene
hmExp=rt[rownames(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)

#volcano
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$adj.P.Val))
yMax=max(abs(allDiff$logFC))
plot(-log10(allDiff$adj.P.Val), allDiff$logFC, xlab="-log10(adj.P.Val)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC>logFoldChange)
points(-log10(diffSub$adj.P.Val), diffSub$logFC, pch=20, col="red",cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC<(-logFoldChange))
points(-log10(diffSub$adj.P.Val), diffSub$logFC, pch=20, col="green",cex=0.8)
abline(h=0,lty=2,lwd=3)
dev.off()



###Video source: http://study.163.com/u/biowolf
######Video source: https://shop119322454.taobao.com
######速科生物: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio
######QQ群:  259208034
