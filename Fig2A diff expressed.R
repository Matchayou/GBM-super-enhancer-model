rm(list=ls())

library(stringr)
library(limma)
library(edgeR)
library(clusterProfiler)

##edgeR
setwd("D:/SuperEnhancer/Fig2 TCGA TF overlap")
mRNA <- read.csv('TCGA-GBM.htseq_counts.tsv.gz', 
                 sep = "\t", header = TRUE, row.names = 1)
rownames(mRNA)=unlist(str_split(rownames(mRNA),"[.]",simplify=T))[,1]
tumor <- colnames(mRNA)[as.integer(substr(colnames(mRNA),14,15)) < 10]#TCGA编号0~9号是肿瘤样本
normal<- colnames(mRNA)[as.integer(substr(colnames(mRNA),14,15)) > 10]#TCGA编号10~19是正常对照
length(tumor) #168
length(normal) #5
normal_sample=mRNA[,normal]
tumor_sample=mRNA[,tumor]
data=cbind(normal_sample,tumor_sample)
exp=2^(data)-1
dim(exp)  #60488 genes
exp<-exp[rowSums((exp) > 1) >= 130,] #~75%
dim(exp) #25428 genes


group_list<- factor(c(rep('normal',length(normal)),rep('tumor',length(tumor))))
design<-model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exp)

DGElist <- DGEList( counts = exp, group = group_list)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]

d<- calcNormFactors(DGElist)
d<- estimateGLMCommonDisp(d, design)
d<- estimateGLMTrendedDisp(d, design)
d<- estimateGLMTagwiseDisp(d, design)

fit <- glmFit(d, design)
lrt<- glmLRT(fit, contrast = c(-1,1))
edgeR_DEG<- topTags(lrt, n = nrow(exp))
edgeR_DEG<- as.data.frame(edgeR_DEG)
edgeR_DEG$group='not-significant'
edgeR_DEG$group[edgeR_DEG$FDR<0.01 & edgeR_DEG$logFC>=1]='up-regulated'
edgeR_DEG$group[edgeR_DEG$FDR<0.01 & edgeR_DEG$logFC<=(-1)]='down-regulated'
table(edgeR_DEG$group)
#down-regulated not-significant    up-regulated 
#     3360           13503            3724 


#gene name transformation
gene<-bitr(rownames(edgeR_DEG),fromType = 'ENSEMBL',
           toType = c('SYMBOL','ENTREZID'),
           OrgDb = 'org.Hs.eg.db')
edgeR_DEG$ENSEMBL=rownames(edgeR_DEG)
edgeR_DEG=edgeR_DEG[,c(6,1,2,3,4,5)]
rownames(edgeR_DEG)=NULL
all_genes=merge(gene,edgeR_DEG,by='ENSEMBL')
dim(all_genes)
table(all_genes$group) 

#down-regulated not-significant    up-regulated 
#3036           11826            3109 
save(all_genes,file='all_genes.Rdata')


load("all_genes.Rdata")
all_genes$group<-factor(all_genes$group,
                        levels=c("up-regulated","not-significant","down-regulated"),
                        labels=c("Up-regulated","Not-significant","Down-regulated"))
library(ggplot2)
pdf("VolcanoNEW.pdf",width=5.5,height=4)
ggplot(all_genes,aes(x=logFC,y= -log10(FDR),color = group))+    
  geom_point(size=0.3)+                     
  labs(x="log2(FoldChange)",y="-log10(FDR)",color="Group")+  
  scale_color_manual(values =c("#BC3C28","grey","#0072B5"))+   
  geom_hline(yintercept=-log10(0.01),linetype=4)+            
  geom_vline(xintercept=c(-1,1),linetype=4)+ 
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))+
  theme_bw()
dev.off()  



