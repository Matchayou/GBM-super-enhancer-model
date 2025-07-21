
setwd("D:/SuperEnhancer/Fig1 SE overview/GENE_TO_ENHANCER")
files<-list.files(pattern="*GENE_TO_ENHANCER.txt")

#cancer 1-34
gene<-c()
for(i in 1:34)
{
  data<-read.table(files[i],header=TRUE)
  gene<-c(gene,data$GENE_NAME)
}
uni_gene_cancer<-unique(gene) 
length(uni_gene_cancer) #9036

#normal 35-37
gene<-c()
for(i in 35:37)
{
  data<-read.table(files[i],header=TRUE)
  gene<-c(gene,data$GENE_NAME)
}
uni_gene_normal<-unique(gene) 
length(uni_gene_normal) #4196

unique_all <- union(uni_gene_cancer,uni_gene_normal)
length(unique_all)

cancer_specific<-setdiff(uni_gene_cancer,uni_gene_normal)
length(cancer_specific)  #5240



load("D:/SuperEnhancer/Fig2 TCGA TF overlap/all_genes.Rdata")
head(all_genes)

SE_exp<-all_genes[all_genes$SYMBOL%in%unique_all,]
dim(SE_exp) #3549

SE_exp_dif<-subset(SE_exp,group=="up-regulated"|group=="down-regulated")
dim(SE_exp_dif) #1286

table(SE_exp_dif$group)
#down-regulated   up-regulated 
#416               870
SE_expup <- subset(SE_exp,group=="up-regulated")

setwd("D:/SuperEnhancer/Fig2 TCGA TF overlap")
write.csv(SE_expup,file="1154 SEs-related.csv")

#Venn
up_TCGA<-subset(all_genes,group=="up-regulated")$SYMBOL
down_TCGA<-subset(all_genes,group=="down-regulated")$SYMBOL
SE_related<-unique_all

length(up_TCGA)
length(down_TCGA)

library(venn)
library(VennDiagram)
library(RColorBrewer)

venn_list<-list(up_TCGA,SE_related)
names(venn_list)<-c("Up-regulated genes",
                    "SE-related genes")
mycolor=brewer.pal(3,"Set3")
venn(venn_list,
     zcolor=mycolor,  ##style??Ĭ????ɫ??bw??????ɫ
     opacity=0.5,#??????ɫ͸????
     box=F,      #?Ƿ????ӱ߿?
     ilcs=1.6,   #???ִ?С
     sncs=1.5)     #???????ִ?С








