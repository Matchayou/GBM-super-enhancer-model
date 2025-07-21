
setwd("D:/SuperEnhancer/Fig3 prognosis1")
load("mRNA_exp.Rdata") #SE_exp
riskscore <- read.table("riskcore.txt",,header=TRUE)
group <- riskscore[,c(1,99)]
common_cols <- intersect(group$sample, colnames(mRNA))
expr_data <- mRNA[, common_cols]
rownames(expr_data) <- rownames(mRNA)
care_genes <- c("NCF2","MTHFS","DUSP6","G6PC3","EN2","HOXB2","DLEU1","LBH","ZEB1-AS1","LINC01265","AGAP2-AS1")
exp_data <- data.matrix(expr_data)


library(oncoPredict)
setwd("D:/SuperEnhancer/药物敏感性")
dir='./DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = exp_data,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' 
              )
setwd("D:/SuperEnhancer/药物敏感性")
drug_predict <- read.csv("DrugPredictions.csv",,header=TRUE)
drug_risk <- cbind(drug_predict,group)
library(ggpubr)
library(ggplot2)
# 定义画图函数

plot_drug_sensitivity <- function(df, drug_name, output_dir, test_method) {
  # 筛选指定药品的敏感性数据
  drug_data <- df[, c(drug_name, "risk")]
  # 计算高低风险组之间的p值
  risk_groups <- unique(drug_data[["risk"]])
  p_value <- NULL
  if (test_method == "t.test") {
    p_value <- t.test(subset(drug_data, drug_data$`risk` == risk_groups[1])[[drug_name]],
                      subset(drug_data, drug_data$`risk` == risk_groups[2])[[drug_name]])
    p_value <- p_value[["p.value"]]
  } else if (test_method == "wilcox.test") {
    p_value <- wilcox.test(subset(drug_data, drug_data$`risk` == risk_groups[1])[[drug_name]],
                           subset(drug_data, drug_data$`risk` == risk_groups[2])[[drug_name]])
    p_value <-p_value[["p.value"]]
  }
  
  # 转换为长格式
  drug_data_long <- reshape2::melt(drug_data, id.vars = "risk",
                                   variable.name = "Drug",
                                   value.name = "Sensitivity")
  # 画箱式图
  p <- ggplot(drug_data_long, aes(x = `risk`, y = Sensitivity, fill = `risk`)) +
    geom_boxplot() +
    labs(x = "Risk Group", y = paste("Drug Sensitivity of", drug_name))+
    stat_compare_means(label = "p.format",method = "wilcox.test",hide.ns = T)+
    theme_bw() +
    theme(panel.grid = element_blank())+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),  ##去掉背景网格
          text = element_text(size = 15),
          axis.text.x = element_text(face="bold",angle = 0,hjust = 1,color = 'black'),
          axis.title.x = element_blank(),
          legend.position =  "top",#绘图图例的位置
          legend.direction = "horizontal",
          legend.title =element_blank())
  # 如果p值小于0.05，则保存图像

  if (p_value < 0.05) {
    pv <- c()
    ggsave(filename = paste(output_dir, "/", drug_name, ".pdf", sep = ""),
           plot = p, dpi = 300)
  }
}

# 设置输出目录
output_dir <- "./drug_sensitivity_plot4"

# 创建输出目录（如果不存在）
dir.create(output_dir, showWarnings = FALSE)

# 循环画出每个药品的箱式图并保存图像
drug_names <- colnames(drug_risk)[2:199]
test_method <- "wilcox.test"  # 差异检验方法
for (drug_name in drug_names) {
  plot_drug_sensitivity(drug_risk, drug_name, output_dir, test_method)
}



plot_drug_t <- function(df, drug_name, output_dir, test_method) {
  pv <- c()
  # 筛选指定药品的敏感性数据
  drug_data <- df[, c(drug_name, "risk")]
  # 计算高低风险组之间的p值
  risk_groups <- unique(drug_data[["risk"]])
  p_value <- NULL
  if (test_method == "t.test") {
    p_value <- t.test(subset(drug_data, drug_data$`risk` == risk_groups[1])[[drug_name]],
                      subset(drug_data, drug_data$`risk` == risk_groups[2])[[drug_name]])
    p_value <- p_value[["p.value"]]
  } else if (test_method == "wilcox.test") {
    p_value <- wilcox.test(subset(drug_data, drug_data$`risk` == risk_groups[1])[[drug_name]],
                           subset(drug_data, drug_data$`risk` == risk_groups[2])[[drug_name]])
    p_value <-p_value[["p.value"]]
  }
  if (p_value < 0.05) {
    pv<-c(drug_name,p_value)
  }
  return(pv)
}
data_p<- c()
for (drug_name in drug_names) {
  value <- plot_drug_t(drug_risk, drug_name, output_dir, test_method)
  data_p <- rbind(data_p,value)
}
write.csv(data_p,file ="p.csv")

drug <- data_p[,1]
drugplot <- drug_predict[,drug]
drugplot <- cbind(drugplot,group)
plotdata <-pivot_longer(data = drugplot,
                        cols = 1:31,
                        names_to = "Drug",
                        values_to = "Sensitivity")

library("ggplot2")
library(ggpubr)
setwd("D:/SuperEnhancer/药物敏感性")
pdf("高低风险组药物敏感性差异1.pdf",width=30,height=8)
ggboxplot(plotdata,
          x = "Drug",
          y = "Sensitivity",
          color = "black",
          fill = "risk",
          xlab = "Drug type",
          ylab = "Drug Sensitivity",
          main = "Drug Sensitivity group by Risk",
          width = 0.7) +
  stat_compare_means(label = "p.signif",method = "wilcox.test",aes(group=risk),hide.ns = F) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()

expr_data <- t(expr_data)
expr_data <- as.data.frame(expr_data)
rownames(drug_predict) <- drug_predict$X
drug_predict <- drug_predict[,c(2:199)]

library(dplyr)
library(tidyr)
library(qvalue)
library(ggplot2)

# 计算药物敏感性和基因表达量的相关性
columns= c("gene name","drug name","cor","p_value")
drug_gene_cor = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(drug_gene_cor) = columns
# 存储所有药物和基因的相关性结果
for (drug_name in colnames(drug_predict)) {
  for(gene_name in care_genes){
  # 根据药物和基因的共同样本，将数据合并为一个数据框
    merged_data <- merge(drug_predict[,drug_name, drop = FALSE],
                         expr_data, by = 'row.names')
    colnames(merged_data)[1] <- 'Sample'
  # 计算相关系数和p值
    cor_result <- cor.test(merged_data[[2]], merged_data[,gene_name], method = 'spearman',exact=FALSE)
  # 存储结果
    cor_re <- c(gene_name,drug_name,cor = cor_result$estimate,p_value = cor_result$p.value)
    drug_gene_cor <- rbind(drug_gene_cor,cor_re)
  }
  }
colnames(drug_gene_cor) = columns
# 将结果转换为数据框，并添加FDR值
write.csv(drug_gene_cor,file ="drug gene cor.csv")


drug_gene_cor$p_value <- as.numeric(drug_gene_cor$p_value)
fdr_result <- qvalue(drug_gene_cor$p_value)$qvalues
drug_gene_cor$FDR <- fdr_result

# 选择FDR小于0.05的结果
significant_cor_df <- filter(drug_gene_cor, FDR < 0.05)
columns= c("gene","drug","cor","p_value","FDR")
colnames(significant_cor_df) = columns
significant_cor_df$cor <- as.numeric(significant_cor_df$cor)
significant_cor_df$FDR <- as.numeric(significant_cor_df$FDR)
cor_df <- filter(significant_cor_df, cor>0.01|cor<0.01)
cor_df <- filter(cor_df, FDR < 0.01)
cor_df$logFDR <- -log10(cor_df$FDR)
cor_df$score <- cor_df$logFDR*cor_df$cor
cor_df$score <- abs(cor_df$score)
rank <- cor_df[order(cor_df$score,decreasing = TRUE),]


# 按照药物名称进行分组，并计算每个药物的评分总数
drug_scores <- aggregate(rank$score, by=list(drug=rank$drug), FUN=sum)

# 按照评分总数进行排序，并选择评分总数前10的药物
top_drugs <- head(arrange(drug_scores, desc(x = x)), n = 30)
topdrug30 <- rank[rank$drug %in% top_drugs$drug, ]
write.csv(topdrug30,file="topdrug30.csv")

# 根据药物和基因的相关性系数和FDR值，画气泡图
ggplot(topdrug30, aes(x = drug, y = gene, size =FDR, color = cor)) +
  geom_point(alpha=1) +
  scale_size(range = c(10, 1)) +
  scale_color_gradient2(low = "purple", mid = "white", high =  "red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 15)) +
  labs(title = 'Correlation between Drug Sensitivity and Gene Expression',
       x = 'Drug', y = 'Gene', size = 'FDR', color = 'cor')+
  theme_bw() +
  theme(panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'))

ggplot(topdrug30, aes(x = gene, y = drug, size =FDR, color = cor)) +
  geom_point(alpha=1) +
  scale_size(range = c(10, 1)) +
  scale_color_gradient2(low = "purple", mid = "white", high =  "red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 15)) +
  labs(title = 'Correlation between Drug Sensitivity and Gene Expression',
       x = 'Gene', y = 'Drug', size = 'FDR', color = 'cor')+
  theme_bw() +
  theme(panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'))

