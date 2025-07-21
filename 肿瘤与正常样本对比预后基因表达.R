setwd("D:/SuperEnhancer/Fig3 prognosis")
load("SE_diff_exp.Rdata") 
tumor <- colnames(SE_exp)[as.integer(substr(colnames(SE_exp),14,15)) < 10]#TCGA编号0~9号是肿瘤样本
normal<- colnames(SE_exp)[as.integer(substr(colnames(SE_exp),14,15)) > 10]
normal_sample=SE_exp[,normal]
tumor_sample=SE_exp[,tumor]
data=cbind(normal_sample,tumor_sample)
genes <- c("BEST3","CD180","DUSP6","G6PC3","HEY1","EN2","HOXB2","DLEU1","IRX5","LBH","ZEB1-AS1","LINC01265","AGAP2-AS1","POM121L9P")
selected_rows <- data[rownames(data) %in% genes,]
selected_rows <- t(selected_rows)
library(reshape2)
df_long <- melt(selected_rows)

df_long$TYPE[df_long$Var1 %in% tumor] <- "tumor"
df_long$TYPE[df_long$Var1 %in% normal] <- "normal"

library(RColorBrewer)
library(viridis)
colors <- brewer.pal(9, "RdPu")
library("ggplot2")
library(ggpubr)
ggboxplot(df_long,
          x = "Var2",
          y = "value",
          color = "black",
          fill = "TYPE",
          xlab = "",
          ylab = "expression",
          main = "Genes expression by group") +
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",hide.ns = F) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))
