#install.packages("pheatmap")

setwd("C:\\Users\\lexb4\\Desktop\\geoBatch\\08.pheatmap")      #���ù���Ŀ¼
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
Type=c(rep("N",25),rep("T",25),rep("N",5),rep("T",5))    #�޸������Ͱ�֢��Ʒ��Ŀ
names(Geo)=colnames(rt)
ann=cbind(Geo,Type)
ann=as.data.frame(ann)

tiff(file="heatmap.tiff",
       width = 20,            #ͼƬ�Ŀ���
       height =25,            #ͼƬ�ĸ߶�
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
######�ٿ�����: http://www.biowolf.cn/
######�������䣺2740881706@qq.com
######����΢��: seqBio
######QQȺ:  259208034