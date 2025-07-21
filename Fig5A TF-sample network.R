
library(stringr)
setwd("D:/SuperEnhancer/Fig1 SE overview/GENE_TO_ENHANCER")
files<-list.files(pattern="*GENE_TO_ENHANCER.txt")
sample_name<-substring(files,1,5)

genes<-c("BEST3","CD180","DUSP6","G6PC3","HEY1","EN2","HOXB2","DLEU1","IRX5","LBH","ZEB1-AS1","LINC01265","AGAP2-AS1","POM121L9P")
network<-c()
net<-c()
for(j in 1:length(genes)){
  select_gene<-genes[j]
  for(i in 1:length(files))
    {
    data<-read.table(files[i],header=TRUE)
    gene<-data$GENE_NAME
    if(sum(gene%in%select_gene)!=0)
      net<-c(sample_name[i],select_gene)
    else net=NULL
    network<-rbind(network,net)
  }
}

head(network)
write.csv(network,
          file="D:/SuperEnhancer/Fig5A TF-sample network/network.csv")


setwd("D:/SuperEnhancer/Fig2 TCGA TF overlap")
#attribute
load("D:/SuperEnhancer/Fig2 TCGA TF overlap/all_genes.Rdata")
up_TCGA<-subset(all_genes,group=="up-regulated")
down_TCGA<-subset(all_genes,group=="down-regulated")
all_diff <- rbind(up_TCGA,down_TCGA)
gene_dif<-all_diff[all_diff$SYMBOL%in%genes,]
dim(gene_dif)
write.csv(gene_dif[,c(2,4,8)],
          file="D:/SuperEnhancer/Fig5A TF-sample network/TF attribute.csv")


