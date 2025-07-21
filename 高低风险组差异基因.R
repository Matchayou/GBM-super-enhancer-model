rm(list=ls())

library(stringr)
library(limma)
library(edgeR)
library(clusterProfiler)

##edgeR
setwd("D:/SuperEnhancer/Fig3 prognosis1")
load("mRNA_exp.Rdata")
mRNA <- as.data.frame(t(mRNA))
mRNA$sample <- rownames(mRNA)
riskscore <- read.table("riskcore.txt",,header=TRUE)
risk <-riskscore[,c(1,99)]
estdata <- merge(mRNA,risk,by = "sample")
high_sample <- estdata[estdata$risk == "high",]
high <- high_sample$sample
low_sample <- estdata[estdata$risk == "low",]
low <- low_sample$sample
setwd("D:/SuperEnhancer/Fig2 TCGA TF overlap")
mRNA <- read.csv('TCGA-GBM.htseq_counts.tsv.gz', 
                 sep = "\t", header = TRUE, row.names = 1)
high_sample=mRNA[,high]
low_sample=mRNA[,low]
data=cbind(high_sample,low_sample)
exp=2^(data)-1
dim(exp)  #60488 genes
exp<-exp[rowSums((exp) > 1) >= 130,] #~75%
dim(exp) #25060 genes

group_list<- factor(c(rep('high',83),rep('low',84)))
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
#187           200172              83 
rownames(edgeR_DEG)=unlist(str_split(rownames(edgeR_DEG),"[.]",simplify=T))[,1]


gene<-bitr(rownames(edgeR_DEG),fromType = 'ENSEMBL',
           toType = c('SYMBOL','ENTREZID'),
           OrgDb = 'org.Hs.eg.db')
edgeR_DEG$ENSEMBL=rownames(edgeR_DEG)
edgeR_DEG=edgeR_DEG[,c(7,1,2,3,4,5,6)]
rownames(edgeR_DEG)=NULL
all_genes=merge(gene,edgeR_DEG,by='ENSEMBL')
dim(all_genes)
table(all_genes$group) 
#down-regulated not-significant    up-regulated 
#172           17643                 72 

setwd("D:/SuperEnhancer/高低风险组差异基因")
save(all_genes,file='risk_genes.Rdata')
load("risk_genes.Rdata")
all_genes$group<-factor(all_genes$group,
                        levels=c("up-regulated","not-significant","down-regulated"),
                        labels=c("Up-regulated","Not-significant","Down-regulated"))
library(ggplot2)
pdf("Volcano_risk.pdf",width=5.5,height=4)
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

risk_exp_dif<-subset(all_genes,group=="Up-regulated"|group=="Down-regulated")
dim(risk_exp_dif) #244

table(risk_exp_dif$group)
#Up-regulated Not-significant  Down-regulated 
#72               0             172 

risk_expup <- subset(risk_exp_dif,group=="Up-regulated")
setwd("D:/SuperEnhancer/高低风险组差异基因")
write.csv(risk_expup,file="72uprisk.csv")
risk_expdown <- subset(risk_exp_dif,group=="Down-regulated")
write.csv(risk_expdown,file="172downrisk.csv")

GO_SEDGs<-enrichGO(SEDGs$ENTREZID,
                   OrgDb = 'org.Hs.eg.db',
                   keyType = "ENTREZID",
                   ont = "BP",pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",qvalueCutoff = 0.01,
                   minGSSize = 10,maxGSSize = 500,readable = TRUE)
dim(as.data.frame(GO_SEDGs))
write.csv(as.data.frame(GO_SEDGs),file="GO terms.csv")
pdf("GO BP top10new.pdf",height=7.68,width=7.14)
dotplot(GO_SEDGs, showCategory = 10,title="top 10 GO BP")+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=25))
dev.off()







