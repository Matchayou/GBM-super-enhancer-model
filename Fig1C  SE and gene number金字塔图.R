
setwd("D:/SuperEnhancer/Fig1 SE overview/GENE_TO_ENHANCER")
files<-list.files(pattern="*GENE_TO_ENHANCER.txt")
files
files <- files[c(1:34)]
gene_num<-c()
for(i in 1:length(files))
{
  data<-read.table(files[i],header=TRUE)
  gene_num[i]<-length(data$GENE_NAME)
}
gene_table<-data.frame(sample=substring(files,1,5),gene_num=gene_num)


setwd("D:/SuperEnhancer/Fig1 SE overview/SuperEnhancer")
SE_files<-list.files(pattern="*SuperEnhancers.table.txt")
SE_files
SE_files <- SE_files[c(1:34)]
SE_num<-c()
for(i in 1:length(SE_files))
{
  data<-read.table(SE_files[i],header=TRUE)
  SE_num[i]<-nrow(data)
}
SE_table<-data.frame(sample=substring(SE_files,1,5),SE_num=SE_num)
number_table<-merge(SE_table,gene_table,by="sample")

colnames(number_table)[c(2,3)]<-c("SE_number","SE_annotated_genes")



number_table$SE_number<- -number_table$SE_number
library(tidyverse)
table2<-number_table%>% 
  pivot_longer(cols=SE_number:SE_annotated_genes,
               names_to="var",values_to="number")

table2$sample[table2$sample=="N1_pe"]<-"Normal1"
table2$sample[table2$sample=="N2_pe"]<-"Normal2"
table2$sample[table2$sample=="N3_pe"]<-"Normal3"

p=ggplot(table2,aes(x=sample,y=number,fill=var))+
  geom_col()+
  coord_flip()+
  scale_fill_manual(values=c("pink1","lightskyblue1"))+
  scale_y_continuous(breaks = seq(-1500, 4000, 500), 
                     labels = as.character(abs(seq(-1500, 4000, 500))), 
                     limits = c(-1500, 4000))+  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=12),
        axis.title = element_text(color = "black",size=15)) +
  geom_hline(yintercept = 0, size = 0.4) 
  #annotate('text',label = 'SEs annotated genes', 32, 3300,size=4)+
  #annotate('text',label = 'SEs', 32, -1400,size=4)  


setwd("D:/SuperEnhancer/Fig1 SE overview")
pdf("Fig1C Number of SEs1.pdf",height=7,width=6)
p+geom_text(data=subset(table2,var=="SE_annotated_genes"),
            aes(x=sample,y=number-300,label=number))+
  geom_text(data=subset(table2,var=="SE_number"),
            aes(x=sample,y=number+200,label=abs(number)))+
  labs(x="Sample",y="Number")+
  theme(legend.position = "bottom")
dev.off()
