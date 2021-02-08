#install.packages("pheatmap")

setwd("C:\\Users\\lexb4\\Desktop\\geoBatch\\08.pheatmap")      #设置工作目录
rt=read.table("diffExp.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#rt=log2(rt+1)
#rt[rt>15]=15

library(pheatmap)
Geo=c(rep("GSE33335",50),rep("GSE56807",10))
Type=c(rep("N",25),rep("T",25),rep("N",5),rep("T",5))    #修改正常和癌症样品数目
names(Geo)=colnames(rt)
ann=cbind(Geo,Type)
ann=as.data.frame(ann)

tiff(file="heatmap.tiff",
       width = 20,            #图片的宽度
       height =25,            #图片的高度
       units ="cm",
       compression="lzw",
       bg="white",
       res=500)
pheatmap(rt, annotation=ann, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
#         scale="row",
         fontsize_row=3,
         fontsize_col=5)
dev.off()



###Video source: http://study.163.com/u/biowolf
######Video source: https://shop119322454.taobao.com
######速科生物: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio
######QQ群:  259208034
