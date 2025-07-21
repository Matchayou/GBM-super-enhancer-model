

setwd("D:/SuperEnhancer/Fig1 SE overview/AnimalTFDB download TF and cofactor")
TF<-read.csv("Homo_sapiens_TF.csv")
TF_co<-read.csv("Homo_sapiens_TF_cofactors.csv")

TF<-TF[,c(2,3,6)]
colnames(TF)<-c("SYMBOL","ENSEMBL","ENTREZID")
TF_co<-TF_co[,c(2,3,5)]
colnames(TF_co)<-c("SYMBOL","ENSEMBL","ENTREZID")



setwd("D:/SuperEnhancer/Fig1 SE overview/GENE_TO_ENHANCER")
files<-list.files(pattern="*GENE_TO_ENHANCER.txt")

TF_num_all<-c()
TF_co_num_all<-c()
for(i in 1:length(files))
{
  data<-read.table(files[i],header=TRUE)
  gene<-data$GENE_NAME
  TF_num<-length(intersect(gene,TF$SYMBOL))
  TF_co_num<-length(intersect(gene,TF_co$SYMBOL))
  
  TF_num_all<-c(TF_num_all,TF_num)
  TF_co_num_all<-c(TF_co_num_all,TF_co_num)
}

MM<-data.frame(sample=substring(files,1,5),
               TF_num=TF_num_all,
               TF_co_num=TF_co_num_all)
library(tidyverse)
MM<-MM%>%pivot_longer(TF_num:TF_co_num,values_to="number",
                      names_to="var")
MM$var<-factor(MM$var,levels=c("TF_num","TF_co_num"),labels=c("TF","Cofactor"))

MM$sample[MM$sample=="N1_pe"]<-"Normal1"
MM$sample[MM$sample=="N2_pe"]<-"Normal2"
MM$sample[MM$sample=="N3_pe"]<-"Normal3"

library(ggplot2)
setwd("D:/SuperEnhancer/Fig1 SE overview")
pdf("Fig 1D Number of TF.pdf",height=7,width=5.5)
ggplot(MM,aes(x=sample,y=number,fill=var))+
  geom_bar(stat="identity",width=0.6)+
  scale_fill_manual(values=c("khaki1","lightpink"))+
  geom_text(aes(x=sample,y=number-15,label=number),position="stack")+
  labs(x="",
       y="Number")+
  coord_flip()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=12),
        axis.title = element_text(color = "black",size=15),
        legend.position = "bottom")
dev.off()