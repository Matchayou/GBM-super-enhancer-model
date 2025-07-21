
#270 diff SE TCGA exp FPKM
library(stringr)
setwd("D:/SuperEnhancer/Fig3 prognosis")
mRNA <- read.csv('TCGA-GBM.htseq_fpkm.tsv.gz', 
                 sep = "\t", header = TRUE, row.names = 1)
rownames(mRNA)=unlist(str_split(rownames(mRNA),"[.]",simplify=T))[,1]
dim(mRNA)
gene_trans<-bitr(rownames(mRNA),fromType = 'ENSEMBL',
                 toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
gene_trans <-unique(gene_trans)
table(duplicated(gene_trans$ENSEMBL))
which(duplicated(gene_trans$ENSEMBL))
unique_gene_trans<- gene_trans[!duplicated(gene_trans$ENSEMBL),]
unique_gene_trans<- unique_gene_trans[!duplicated(unique_gene_trans$SYMBOL),]
mRNA$ENSEMBL=rownames(mRNA)
mRNA=merge(unique_gene_trans,mRNA,by='ENSEMBL')
rownames(mRNA)<-mRNA$SYMBOL
mRNA <- mRNA[,-c(1,2)]
write.csv(mRNA,"mRNA_exp.csv")
save(mRNA,file="mRNA_exp.Rdata") 


SE_exp_dif<-read.csv("D:/SuperEnhancer/Fig2 TCGA TF overlap/1154 SEs-related.csv")
SE_exp<-mRNA[rownames(mRNA)%in% SE_exp_dif$ENSEMBL,]
dim(SE_exp) #1152 173

library(clusterProfiler)
gene_trans<-bitr(rownames(SE_exp),fromType = 'ENSEMBL',
                 toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
gene_trans <-unique(gene_trans)
table(duplicated(gene_trans$SYMBOL))
which(duplicated(gene_trans$SYMBOL))

#table(duplicated(gene_trans$ENSEMBL))
#which(duplicated(gene_trans$ENSEMBL))
#unique_gene_trans<- gene_trans[!duplicated(gene_trans$ENSEMBL),]
#sort(table(gene_trans$ENSEMBL),decreasing=TRUE)

SE_exp$ENSEMBL=rownames(SE_exp)
SE_exp=merge(gene_trans,SE_exp,by='ENSEMBL')
rownames(SE_exp)<-SE_exp$SYMBOL
SE_exp <- SE_exp[,-c(1,2)]
save(SE_exp,file="SE_diff_exp.Rdata") 



library(survival)
library(survminer)
library(survivalROC)
library(tibble)
library(stringr)

setwd("D:/SuperEnhancer/Fig3 prognosis")
data<-read.table("TCGA-GBM.survival.tsv",header=T)
data$sample<-gsub("-",".",data$sample)
head(data)
dim(data)  #649 samples


load("SE_diff_exp.Rdata") #SE_exp
SE_exp<-t(SE_exp)
dim(SE_exp)  #173 samples 

SE_exp<-rownames_to_column(as.data.frame(SE_exp),var="sample")
data<-merge(data[,c(1,2,4)],SE_exp,by="sample")
dim(data) #167 


#Univariate Cox
outTab=data.frame()
for(i in colnames(data[,4:ncol(data)])){
  cox <- coxph(Surv(data$OS.time, data$OS) ~ data[,i], data = data)
  a = summary(cox)
  b = as.data.frame(a$coefficients)[1,]
  c = as.data.frame(a$conf.int)[1,]
  outTab=rbind(outTab,cbind(b,c))
}
rownames(outTab)=colnames(data)[4:ncol(data)]
outTab_0.05=subset(outTab,outTab$`Pr(>|z|)`<0.05)
dim(outTab_0.05)  #236


##KM analysis
my.surv <- Surv(data$OS.time, data$OS)
log_rank_p <- apply(data[,4:ncol(data)], 2, function(values1){
  group=ifelse(values1>median(values1),'high','low')
  kmfit2 <- survfit(my.surv~group,data=data)
  data.survdiff=survdiff(my.surv~group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
})
log_rank_p <- log_rank_p[log_rank_p<0.05]
cc=names(log_rank_p)
length(cc)  #134


#overlap of Univariate Cox and KM survival analysis
gene2<-intersect(cc,rownames(outTab_0.05))
length(gene2) #90

#data<-cbind(data[,c(1:3)],data[,gene2])
data<-cbind(data[,c(1:3)],data[,rownames(outTab_0.05)])

###############################
#lasso
library(glmnet)
x=as.matrix(data[,-c(1:3)])
y=data.matrix(Surv(time = data$OS.time,event = data$OS))
fit=glmnet(x,y,family = 'cox',alpha = 1)
plot(fit,xvar='lambda',label=TRUE)
set.seed(1000)
lasso_fit=cv.glmnet(x,y,family='cox',type.measure = 'deviance')
plot(lasso_fit)
coefficient=coef(lasso_fit,s=lasso_fit$lambda.min)  #lasso_fit$lambda.min
active_index<-which(as.numeric(coefficient)!=0)
active_coefficient=as.numeric(coefficient)[active_index]
gene=rownames(coefficient)[active_index]
gene
length(gene) #26
##################################
#"SH3BP2"       "BEST3"        "POM121L9P"    "LOXL1"        "CD180"        "DUSP6"       
#"G6PC3"        "FCN3"         "OSMR"         "TIMP4"        "HEY1"         "EN2"         
#"RPL13"        "SCARA3"       "HOXB2"        "DLEU1"        "DLEU7-AS1"    "IRX5"        
#"UPP1"         "LBH"          "TSPAN4"       "ZEB1-AS1"     "CRNDE"        "LOC101927480"
#"LINC01265"    "AGAP2-AS1"   

data_ori =data
colnames(data)[184] <- 'DLEU7_AS1'
colnames(data)[228] <- "ZEB1_AS1"
colnames(data)[237] <- "AGAP2_AS1"

#colnames(data)[74] <- 'DLEU7_AS1'
#colnames(data)[92] <- "ZEB1_AS1"
#colnames(data)[96] <- "AGAP2_AS1"

#cox_multi1<-coxph(Surv(OS.time,OS)~AEBP1+NCF2+LOXL1+MTHFS+DUSP6+G6PC3+FCN3+OSMR+EN2+SCARA3+ISG20+HOXB2+DLEU1+DLEU7_AS1+UPP1+LBH+TSPAN4+ZEB1_AS1+CRNDE+LOC101927480+LINC01265+AGAP2_AS1,data = data)
cox_multi1<-coxph(Surv(OS.time,OS)~SH3BP2+BEST3+POM121L9P+LOXL1+CD180+DUSP6+G6PC3+FCN3+OSMR+TIMP4+HEY1+EN2+RPL13+SCARA3+HOXB2+DLEU1+DLEU7_AS1+IRX5+UPP1+LBH+TSPAN4+ZEB1_AS1+CRNDE+LOC101927480+LINC01265+AGAP2_AS1,data = data)
summary(cox_multi1)
cox1=step(cox_multi1,direction = 'both')


#c1("BEST3","CD180","DUSP6","G6PC3","HEY1","EN2","HOXB2","DLEU1","IRX5","LBH","ZEB1_AS1","LINC01265","AGAP2_AS1","POM121L9P")

riskScore=predict(cox1,type = 'risk',data=data)
risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))
data$riskScore=riskScore
data$risk=risk
data_risk =data
write.table(data_risk,file = "riskcore.txt",row.names = F,quote = F,sep = "\t")
save(data,file="dataprognosis.Rdata")

load("dataprognosis.Rdata")
#KM curve
sfit <- survfit(Surv(OS.time, OS)~risk, data=data)
ggsurv <- ggsurvplot(sfit, pval.method=T,
                     surv.median.line ="hv",cumevents=T,xlab = "Days",
                     legend = "top",title  = "combined-genes",
                     risk.table.y.text = F,conf.int=F,ncensor.plot = F,
                     risk.table = F, pval=TRUE,data=data)
pdf("KM curve trainnew1.pdf",width=5,height=5)
ggsurv
dev.off()

my.surv <- Surv(data$OS.time, data$OS)
log_rank_p <- apply(data[,4:ncol(data)], 2, function(values1){
  group=ifelse(values1>median(values1),'high','low')
  kmfit2 <- survfit(my.surv~group,data=data)
  data.survdiff=survdiff(my.surv~group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
})
log_rank_p <- log_rank_p[log_rank_p<0.05]
cc=names(log_rank_p)
length(cc)  #134


#DLEU1
DLEU1 <- data$DLEU1
risk_DLEU1 <- as.vector(ifelse(DLEU1>median(DLEU1),"high","low"))
data_DLEU1 =data
data_DLEU1$risk_DLEU1=risk_DLEU1
sfit <- survfit(Surv(OS.time, OS)~risk_DLEU1, data=data_DLEU1)
ggsurv <- ggsurvplot(sfit, pval.method=T,
                     surv.median.line ="hv",cumevents=T,xlab = "Days",
                     legend = "top",title  = "combined-genes",
                     risk.table.y.text = F,conf.int=F,ncensor.plot = F,
                     risk.table = F, pval=TRUE,data=data)
pdf("KM curve trainnewDLEU1.pdf",width=5,height=5)
ggsurv
dev.off()

#HEY1
HEY1 <- data$HEY1
risk_HEY1 <- as.vector(ifelse(HEY1>median(HEY1),"high","low"))
data_HEY1 =data
data_HEY1$risk_DLEU1=risk_HEY1
sfit <- survfit(Surv(OS.time, OS)~risk_HEY1, data=data_HEY1)
ggsurv <- ggsurvplot(sfit, pval.method=T,
                     surv.median.line ="hv",cumevents=T,xlab = "Days",
                     legend = "top",title  = "combined-genes",
                     risk.table.y.text = F,conf.int=F,ncensor.plot = F,
                     risk.table = F, pval=TRUE,data=data)
pdf("KM curve trainnewHEY1.pdf",width=5,height=5)
ggsurv
dev.off()

#ZEB1_AS1的KM曲线
ZEB1_AS1 <- data$ZEB1_AS1
risk_ZEB1_AS1 <- as.vector(ifelse(ZEB1_AS1>median(ZEB1_AS1),"high","low"))
data_ZEB1_AS1 =data
data_ZEB1_AS1$risk_ZEB1_AS1=risk_ZEB1_AS1
sfit <- survfit(Surv(OS.time, OS)~risk_ZEB1_AS1, data=data_ZEB1_AS1)
ggsurv <- ggsurvplot(sfit, pval.method=T,
                     surv.median.line ="hv",cumevents=T,xlab = "Days",
                     legend = "top",title  = "combined-genes",
                     risk.table.y.text = F,conf.int=F,ncensor.plot = F,
                     risk.table = F, pval=TRUE,data=data)
pdf("KM curve trainnewZEB1_AS1.pdf",width=5,height=5)
ggsurv
dev.off()

#BEST3
BEST3 <- data$BEST3
risk_BEST3 <- as.vector(ifelse(BEST3>median(BEST3),"high","low"))
data_BEST3 =data
data_BEST3$risk_BEST3=risk_BEST3
sfit <- survfit(Surv(OS.time, OS)~risk_BEST3, data=data_BEST3)
ggsurv <- ggsurvplot(sfit, pval.method=T,
                     surv.median.line ="hv",cumevents=T,xlab = "Days",
                     legend = "top",title  = "combined-genes",
                     risk.table.y.text = F,conf.int=F,ncensor.plot = F,
                     risk.table = F, pval=TRUE,data=data)
pdf("KM curve trainnewBEST3.pdf",width=5,height=5)
ggsurv
dev.off()





#ROC curve
pdf("train ROC1.pdf",width=5,height=5)
roc=survivalROC(Stime = data$OS.time, status=data$OS, 
                marker = data$riskScore,predict.time =365, method="KM")
roc1=survivalROC(Stime = data$OS.time, status=data$OS, 
                 marker = data$riskScore,predict.time =365*2, method="KM") #15/12
roc2=survivalROC(Stime = data$OS.time, status=data$OS, 
                 marker = data$riskScore,predict.time =365*3, method="KM") #18/12

plot(roc$FP, roc$TP, type="l",main="ROC Curve", col='#8B0000',xlim=c(0,1),bg="white", ylim=c(0,1),xlab="False positive rate", ylab="True positive rate",lwd = 2,cex.lab=1.2, cex.axis=1.2, font=2,font.lab=2)
par(new=TRUE)
plot(roc1$FP, roc1$TP, type="l", axes = FALSE,col='#00BFC4',xlab="", ylab="",lwd = 2)
par(new=TRUE)
plot(roc2$FP, roc2$TP, type="l", axes = FALSE,col='#F8766D',xlab="", ylab="",lwd = 2)
abline(0,1)
legend(x=0.4,y=0.4,
       c(paste("1 years AUC:",round(roc$AUC,3)),
         paste("2 years AUC:",round(roc1$AUC,3)),
         paste("3 years AUC:",round(roc2$AUC,3))),
       col=c('#8B0000','#00BFC4','#F8766D'),
       x.intersp=1, y.intersp=0.8,
       bty = 'n',text.font=2,lty=1,lwd=2,
       seg.len=1,cex=1.0)
dev.off()


# Forest plot多因素
#c1 <- c("NCF2","MTHFS","DUSP6","G6PC3","EN2","HOXB2"," DLEU1","LBH","ZEB1_AS1","LINC01265","AGAP2_AS1")
c1 <- c("BEST3","CD180","DUSP6","G6PC3","HEY1","EN2","HOXB2","DLEU1","IRX5","LBH","ZEB1-AS1","LINC01265","AGAP2-AS1","POM121L9P")
onecoxp <- outTab_0.05[rownames(outTab_0.05) %in% c1,]

multiCoxSum <- summary(cox1)
out_multi <- data.frame()
out_multi <- cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

out_multi <- as.data.frame(cbind(id=row.names(out_multi),out_multi)) 

out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
hz <- paste(round(out_multi$HR,3),
            "(",round(out_multi$HR.95L,3),
            "-",round(out_multi$HR.95H,3),")",sep = "")

tabletext <- cbind(c(NA,"Gene",out_multi$id),
                   c(NA,"Coefficient",round(out_multi$coef,3)),
                   c(NA,"P value",ifelse(out_multi$pvalue<0.001,"P < 0.001",round(out_multi$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))


library(forestplot)
pdf("Forestnew1.pdf",width=7,height=6)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,out_multi$HR),
           lower=c(NA,NA,out_multi$HR.95L), 
           upper=c(NA,NA,out_multi$HR.95H), 
           boxsize=0.3,lwd.ci=2,  
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=1,     
           colgap=unit(5,"mm"),   
           xticks = c(0.5, 1,1.5), 
           lwd.xaxis=1,            
           lineheight = unit(0.8,"cm"), 
           graphwidth = unit(.3,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), 
                           "17" = gpar(lwd=2, col="black")),  #"12"-nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),
          
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.5),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")         
dev.off()         

#forestplot 单因素
outTab=data.frame()
for(i in colnames(data[,4:ncol(data)])){
  cox <- coxph(Surv(data$OS.time, data$OS) ~ data[,i], data = data)
  a = summary(cox)
  b = as.data.frame(a$coefficients)[1,]
  c = as.data.frame(a$conf.int)[1,]
  outTab=rbind(outTab,cbind(b,c))
}
rownames(outTab)=colnames(data)[4:ncol(data)]
outTab_0.05=subset(outTab,outTab$`Pr(>|z|)`<0.05)
dim(outTab_0.05)  #236
c1 <- c("BEST3","CD180","DUSP6","G6PC3","HEY1","EN2","HOXB2","DLEU1","IRX5","LBH","ZEB1_AS1","LINC01265","AGAP2_AS1","POM121L9P")
onecoxp <- outTab_0.05[rownames(outTab_0.05) %in% c1,]
out_uni <- data.frame()
out_uni <- cbind(
  coef=onecoxp[,c(1)],
  HR=onecoxp[,c(2)],
  HR.95L=onecoxp[,c(8)],
  HR.95H=onecoxp[,c(9)],
  pvalue=onecoxp[,c(5)])
out_uni <- as.data.frame(cbind(id=row.names(onecoxp),out_uni)) 

out_uni[,2:ncol(out_uni)] <- as.numeric(unlist(out_uni[,2:ncol(out_uni)]))
hz <- paste(round(out_uni$HR,3),
            "(",round(out_uni$HR.95L,3),
            "-",round(out_uni$HR.95H,3),")",sep = "")

tabletext <- cbind(c(NA,"Gene",out_uni$id),
                   c(NA,"Coefficient",round(out_uni$coef,3)),
                   c(NA,"P value",ifelse(out_uni$pvalue<0.001,"P < 0.001",round(out_uni$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))

library(forestplot)
pdf("Forestuni.pdf",width=7,height=6)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,out_uni$HR),
           lower=c(NA,NA,out_uni$HR.95L), 
           upper=c(NA,NA,out_uni$HR.95H), 
           boxsize=0.3,lwd.ci=2,  
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=1,     
           colgap=unit(5,"mm"),   
           xticks = c(0.5, 1,1.5), 
           lwd.xaxis=1,            
           lineheight = unit(0.8,"cm"), 
           graphwidth = unit(.3,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), 
                           "17" = gpar(lwd=2, col="black")),  #"12"-nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),
           
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.5),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")         
dev.off()         




data=data[order(data$riskScore),]
p=ggplot(data = data)+
  geom_point(aes(x=seq(0:166),y=riskScore,color=risk))+
  scale_x_continuous(breaks = seq(0,166,50))+
  #scale_y_continuous(breaks = seq(0,7,1))+
  geom_vline(aes(xintercept=84),colour='#BB0000',linetype='dashed')
p1=p+theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank())+
  labs(x='Patients(increasing riskScore)',y='Risk Score')+
  theme(axis.line = element_line(size = 1,colour = 'black'))


data$status=ifelse(data$OS==1,'Dead','Alive')
cols <- c("Alive" = "#00BFC4", "Dead" = "#F8766D")
p2=ggplot()+geom_point(data = data,aes(x=seq(0:166),y=OS.time,color=status))+
  scale_x_continuous(breaks = seq(0,166,50))+
  scale_y_continuous(breaks = seq(0,3000,1000))+
  geom_vline(aes(xintercept=84),colour='#BB0000',linetype='dashed')+
  scale_color_manual(values = cols)
p3=p2+theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
  labs(x='Patients(increasing riskScore)',y='Survival Days')+
  theme(axis.line = element_line(size = 1,colour = 'black'))

#c("DUSP6","BEST3","ZEB1_AS1","HEY1","DLEU1","POM121L9P","LINC01265","IRX5","HOXB2","AGAP2_AS1","EN2","CD180","LBH","G6PC3")
gene=c("BEST3","ZEB1_AS1","HEY1","DLEU1","POM121L9P","LINC01265","IRX5","HOXB2","AGAP2_AS1","EN2","CD180","LBH","G6PC3","DUSP6")
library(pheatmap)
pheatmap_data=t(data[,gene])
p4=pheatmap(pheatmap_data,scale = 'row',
            fontsize_row = 11,
            color = colorRampPalette(c('blue','white','red'))(100),
            show_colnames = F,show_rownames = T,
            cluster_rows = F,cluster_cols = F)

library(cowplot)
p4=ggplotify::as.ggplot(p4)

pdf("risk score trainnew1.pdf")
plot_grid(p3,p1,p4,ncol = 1,axis = 'l',align = 'v')
dev.off()

DUSP6+PTPRN2+DLEU1+EN2+UPP1+LBH+TSPAN4+CRNDE

#nomogram
library(rms)
f_cph <- cph(Surv(OS.time, OS) ~ BEST3 + CD180 + DUSP6 + G6PC3 + HEY1 + EN2 + HOXB2 + DLEU1 + IRX5 + LBH + ZEB1_AS1 + LINC01265 + AGAP2_AS1 + POM121L9P,
             x=T, y=T, surv=T,
             data=data)
print(f_cph)

ddist <- datadist(data)
options(datadist='ddist')
med  <- Quantile(f_cph)
surv <- Survival(f_cph) 

pdf("nomogramnew.pdf",width=8, height=6)
plot(nomogram(f_cph, fun=list(function(x) surv(365, x),
                              function(x) surv(365*2, x),
                              function(x) surv(365*3, x)),
              funlabel=c("1-year Survival Probability", 
                         "2-year Survival Probability",
                         "3-year Survival Probability")))
dev.off()





#Univariate Cox model: KM curve
#gene=c("HDAC4","ZNF22","LBH","HOXB2","EN2","AEBP1")
cox<-coxph(Surv(OS.time,OS)~LBH,data = data)  
summary(cox)

riskScore=predict(cox,type = 'risk',data=data)
risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))
data$riskScore=riskScore
data$risk=risk


sfit <- survfit(Surv(OS.time, OS)~risk, data=data)
ggsurv <- ggsurvplot(sfit, pval.method=T,
                     surv.median.line ="hv",cumevents=T,xlab = "Days",
                     legend = "top",title  = "LBH",
                     risk.table.y.text = F,conf.int=F,ncensor.plot = F,
                     risk.table = F, pval=TRUE,data=data)
#pdf("KM curve LBH.pdf",width=5,height=5)
ggsurv
#dev.off()

#ROC curve
#pdf("ROC LBH.pdf",width=5,height=5)
roc=survivalROC(Stime = data$OS.time, status=data$OS, 
                marker = data$riskScore,predict.time =365, method="KM")
roc1=survivalROC(Stime = data$OS.time, status=data$OS, 
                 marker = data$riskScore,predict.time =365*2, method="KM")
roc2=survivalROC(Stime = data$OS.time, status=data$OS, 
                 marker = data$riskScore,predict.time =365*3, method="KM")

plot(roc$FP, roc$TP, type="l",main="LBH ROC Curve", col='#8B0000',xlim=c(0,1),bg="white", ylim=c(0,1),xlab="False positive rate", ylab="True positive rate",lwd = 2,cex.lab=1.2, cex.axis=1.2, font=2,font.lab=2)
par(new=TRUE)
plot(roc1$FP, roc1$TP, type="l", axes = FALSE,col='#00BFC4',xlab="", ylab="",lwd = 2)
par(new=TRUE)
plot(roc2$FP, roc2$TP, type="l", axes = FALSE,col='#F8766D',xlab="", ylab="",lwd = 2)
abline(0,1)
legend(x=0.4,y=0.4,
       c(paste("1 years AUC:",round(roc$AUC,3)),
         paste("2 years AUC:",round(roc1$AUC,3)),
         paste("3 years AUC:",round(roc2$AUC,3))),
       col=c('#8B0000','#00BFC4','#F8766D'),
       x.intersp=1, y.intersp=0.8,
       bty = 'n',text.font=2,lty=1,lwd=2,
       seg.len=1,cex=1.0)
#dev.off()





#phenotype
phenotype<-read.csv("TCGA-GBM.GDC_phenotype.tsv.gz", 
                    sep = "\t", header = TRUE, row.names = 1)
dim(phenotype) #671 94
phen<-phenotype[,c(1,50,19,30)]
phen$sample<-gsub("-",".",rownames(phen))

mer_data<-merge(phen,data[,c("sample","risk","riskScore","OS","OS.time")],by="sample")
dim(mer_data)
mydata <- mer_data
mydata[mydata == ""] <- NA
mydata <- na.omit(mydata)
#not significant
t.test(mer_data[mer_data$gender.demographic=="male","riskScore"],
       mer_data[mer_data$gender.demographic=="female","riskScore"])
chisq.test(table(mer_data$gender.demographic,mer_data$risk))

ggboxplot(mer_data,
          x = "gender.demographic",
          y = "riskScore",
          color = "black",
          fill = "gender.demographic",
          xlab = "",
          ylab = "riskScore",
          main = "riskScore by gender") +
  stat_compare_means(comparisons = list(c("female", "male")),methods="wilcox.test") +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))



cor.test(mer_data$riskScore,mer_data$age_at_initial_pathologic_diagnosis)

mer_data$age_group<-ifelse(mer_data$age_at_initial_pathologic_diagnosis>40,">40","<40")
dev.off()

ggboxplot(mer_data,
          x = "age_group",
          y = "riskScore",
          color = "black",
          fill = "age_group",
          xlab = "",
          ylab = "riskScore",
          main = "riskScore by age") +
  stat_compare_means(comparisons = list(c(">40", "<40")),methods="wilcox.test") +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))
boxplot(riskScore~age_group,data=mer_data)
chisq.test(table(mer_data$age_group,mer_data$risk))

#not significant
a <- table(mer_data$risk,mer_data$radiation_therapy)
chisq.test(table(mer_data$risk,mer_data$radiation_therapy)[,-1])

t.test(mer_data[mer_data$radiation_therapy=="YES","riskScore"],
       mer_data[mer_data$radiation_therapy=="NO","riskScore"])
new_data <- mer_data[!is.na(mer_data$radiation_therapy),]

ggboxplot(mydata,
          x = "radiation_therapy",
          y = "riskScore",
          color = "black",
          fill = "radiation_therapy",
          xlab = "",
          ylab = "riskScore",
          main = "riskScore by radiation_therapy") +
  stat_compare_means(comparisons = list(c("YES", "NO")),methods="wilcox.test") +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))



chisq.test(table(mer_data$risk,mer_data$karnofsky_performance_score))

cor.test(mer_data$riskScore,mer_data$karnofsky_performance_score)
mydata$KPS<-ifelse(mydata$karnofsky_score>55,">55","<55")
ggboxplot(mydata,
          x = "KPS",
          y = "riskScore",
          color = "black",
          fill = "KPS",
          xlab = "",
          ylab = "riskScore",
          main = "riskScore by KPS") +
  stat_compare_means(comparisons = list(c(">55", "<55")),methods="wilcox.test") +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))







colnames(mydata)[2] <- "age"
colnames(mydata)[3] <- "gender"
colnames(mydata)[4] <- "karnofsky_score"


mydata$radiation_therapy <- factor(ifelse(mydata$radiation_therapy == "YES", 1, 0))
#mydata$gender <- factor(ifelse(mydata$gender == "female", 1, 0))

outcli=data.frame()
for(i in colnames(mydata[,c(2,3,4,5,7)])){
  cox <- coxph(Surv(mydata$OS.time, mydata$OS) ~ mydata[,i], data = mydata)
  a = summary(cox)
  b = as.data.frame(a$coefficients)[1,]
  c = as.data.frame(a$conf.int)[1,]
  outcli=rbind(outcli,cbind(b,c))
}
rownames(outcli)=colnames(mydata)[c(2,3,4,5,7)]
outcli_0.05=subset(outcli,outcli$`Pr(>|z|)`<0.05)
dim(outcli)
cli_uniti <- data.frame()
cli_uniti <- cbind(
  coef=outcli[,"coef"],
  HR=outcli[,"exp(coef)"],
  HR.95L=outcli[,c(8)],
  HR.95H=outcli[,c(9)],
  pvalue=outcli[,"Pr(>|z|)"])
cli_uniti <- as.data.frame(cbind(id=row.names(outcli),cli_uniti)) 
rownames(cli_uniti) <- cli_uniti$id
cli_uniti[,2:ncol(cli_uniti)] <- as.numeric(unlist(cli_uniti[,2:ncol(cli_uniti)]))
hz <- paste(round(cli_uniti$HR,3),
            "(",round(cli_uniti$HR.95L,3),
            "-",round(cli_uniti$HR.95H,3),")",sep = "")
tabletext <- cbind(c(NA,"Gene",cli_uniti$id),
                   c(NA,"Coefficient",round(cli_uniti$coef,3)),
                   c(NA,"P value",ifelse(cli_uniti$pvalue<0.001,"P < 0.001",round(cli_uniti$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))

library(forestplot)
pdf("Forestcliuni.pdf",width=7,height=6)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,cli_uniti$HR),
           lower=c(NA,NA,cli_uniti$HR.95L), 
           upper=c(NA,NA,cli_uniti$HR.95H), 
           boxsize=0.2,lwd.ci=2,  
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=1,     
           colgap=unit(5,"mm"),   
           xticks = c(0, 0.5,1.0,1.5,2.0), 
           lwd.xaxis=1,            
           lineheight = unit(0.8,"cm"), 
           graphwidth = unit(.3,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "8" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black")),  #"12"-nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),
           
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.5),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio",
           )            
dev.off()         




#多因素Cox回归结果
res.cox <- coxph(Surv(OS.time,OS) ~age + gender+karnofsky_score+radiation_therapy+riskScore, data = mydata)
summary(res.cox)
#Surv(OS.time, OS) ~ riskScore + radiation_therapy

multiCoxcli <- summary(res.cox)
cli_multi <- data.frame()
cli_multi <- cbind(
  coef=multiCoxcli$coefficients[,"coef"],
  HR=multiCoxcli$conf.int[,"exp(coef)"],
  HR.95L=multiCoxcli$conf.int[,"lower .95"],
  HR.95H=multiCoxcli$conf.int[,"upper .95"],
  pvalue=multiCoxcli$coefficients[,"Pr(>|z|)"])

cli_multi <- as.data.frame(cbind(id=row.names(cli_multi),cli_multi)) 

cli_multi[,2:ncol(cli_multi)] <- as.numeric(unlist(cli_multi[,2:ncol(cli_multi)]))
hz <- paste(round(cli_multi$HR,3),
            "(",round(cli_multi$HR.95L,3),
            "-",round(cli_multi$HR.95H,3),")",sep = "")

tabletext <- cbind(c(NA,"Gene",cli_multi$id),
                   c(NA,"Coefficient",round(cli_multi$coef,3)),
                   c(NA,"P value",ifelse(cli_multi$pvalue<0.001,"P < 0.001",round(cli_multi$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))


library(forestplot)
pdf("Forestclimul.pdf",width=7,height=6)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,cli_multi$HR),
           lower=c(NA,NA,cli_multi$HR.95L), 
           upper=c(NA,NA,cli_multi$HR.95H), 
           boxsize=0.2,lwd.ci=2,  
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=1,     
           colgap=unit(5,"mm"),   
           xticks = c(0.5, 1,1.5), 
           lwd.xaxis=1,            
           lineheight = unit(0.8,"cm"), 
           graphwidth = unit(.3,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "8" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black")),  #"12"-nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),
           
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.5),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")         
dev.off()         






#结果不可靠

re.cox <- coxph(Surv(OS.time,OS) ~age+radiation_therapy+riskScore, data = mydata)
mydata$predicted <- exp(-predict(re.cox, type = "lp"))
library(rms)
sum_cph <- cph(Surv(OS.time,OS) ~ age+radiation_therapy+riskScore,
             x=T, y=T, surv=T,
             time.inc = 365*1,
             data=mydata)
summary(sum_cph)

ddist <- datadist(mydata)
options(datadist='ddist')
med  <- Quantile(sum_cph)
surv <- Survival(sum_cph) 

pdf("nomogramcli.pdf",width=8, height=6)
plot(nomogram(sum_cph, fun=list(function(x) surv(365, x),
                              function(x) surv(365*2, x),
                              function(x) surv(365*3, x)),
              funlabel=c("1-year Survival Probability", 
                         "2-year Survival Probability",
                         "3-year Survival Probability")))
dev.off()

mydata$predicted <- exp(-predict(sum_cph, type = "lp"))
sum_cph_1 <- cph(Surv(OS.time,OS) ~ age+radiation_therapy+riskScore,
               x=T, y=T, surv=T,
               time.inc = 365*1,
               data=mydata)
sum_cph_2 <- cph(Surv(OS.time,OS) ~ age+radiation_therapy+riskScore,
               x=T, y=T, surv=T,
               time.inc = 365*2,
               data=mydata)
sum_cph_3 <- cph(Surv(OS.time,OS) ~ age+radiation_therapy+riskScore,
                 x=T, y=T, surv=T,
                 time.inc = 365*3,
                 data=mydata)
cal <- calibrate(sum_cph_1, cmethod = "KM", method = "boot", u = 365*1, m = 40, B = 400)
cal_2 <- calibrate(sum_cph_2, cmethod = "KM", method = "boot", u = 365*2, m = 40, B = 400)
cal_3 <- calibrate(sum_cph_3, cmethod = "KM", method = "boot", u = 365*3, m = 40, B = 200)
plot(cal, lwd = 2, lty = 1, errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = "Nomogram-Predicted Probability of 1-Year OS", ylab = "Actual 1-Year OS (proportion)",
     col = c(rgb(192, 98, 83, maxColorValue = 255)), 
     subtitles = FALSE, 
     xlim = c(0,1), 
     ylim = c(0, 1),
     main = "Calibrate plot")
lines(cal[, c("mean.predicted", "KM")], 
      type = "l", 
      lwd = 2,
      col = c(rgb(192, 98,83, maxColorValue = 255)), 
      pch = 16)
abline(0, 1, lty = 3, lwd = 2, col = c(rgb(0, 118, 192, maxColorValue = 255)))

plot(cal_2, lwd = 2, lty = 1, errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = "Nomogram-Predicted Probability of 2-Year OS", ylab = "Actual 1-Year OS (proportion)",
     col = c(rgb(192, 98, 83, maxColorValue = 255)), 
     subtitles = FALSE, 
     xlim = c(0,1), 
     ylim = c(0, 1),
     main = "Calibrate plot")
lines(cal[, c("mean.predicted", "KM")], 
      type = "l", 
      lwd = 2,
      col = c(rgb(192, 98,83, maxColorValue = 255)), 
      pch = 16)
abline(0, 1, lty = 3, lwd = 2, col = c(rgb(0, 118, 192, maxColorValue = 255)))

plot(cal_3, lwd = 2, lty = 1, errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = "Nomogram-Predicted Probability of 3-Year OS", ylab = "Actual 1-Year OS (proportion)",
     col = c(rgb(192, 98, 83, maxColorValue = 255)), 
     subtitles = FALSE, 
     xlim = c(0,1), 
     ylim = c(0, 1),
     main = "Calibrate plot")
lines(cal[, c("mean.predicted", "KM")], 
      type = "l", 
      lwd = 2,
      col = c(rgb(192, 98,83, maxColorValue = 255)), 
      pch = 16)
abline(0, 1, lty = 3, lwd = 2, col = c(rgb(0, 118, 192, maxColorValue = 255)))


newdata <- mydata[60:120, ]  #用来校准的数据，这里从源数据中调取了部分
pred.lg <- predict(sum_cph, newdata)  #每位患者的风险评分
newdata$prob <- 1/(1 + exp(-pred.lg))  #将pred.lg做数据转化，数值更直观
prob <- newdata$prob
if (!require("PredictABEL")) {
  install.packages("PredictABEL")
}
library(PredictABEL)
plotCalibration(data = newdata, cOutcome = 8, predRisk = prob, groups = 10, plottitle = "Calibration plot")  #3为结局指标所在列数
