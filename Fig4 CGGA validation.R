rm(list=ls())
library(survival)
library(survminer)
library(survivalROC)
library(tibble)
library(stringr)

setwd('D:/SuperEnhancer/Fig4 external validation')

#CGGA325
CGGA<-read.table("CGGA.mRNAseq_325.RSEM-genes.20200506.txt",header=T) #325
dim(CGGA)
which(CGGA$Gene_Name == "ZEB1-AS1")
which(CGGA$Gene_Name == "AGAP2-AS1")
CGGA$Gene_Name[23636] <- "ZEB1_AS1"
CGGA$Gene_Name[1232] <- "AGAP2_AS1"

gene=c("DUSP6","BEST3","ZEB1_AS1","HEY1","DLEU1","POM121L9P","LINC01265","IRX5","HOXB2","AGAP2_AS1","EN2","CD180","LBH","G6PC3")

CGGA_gene<-CGGA[CGGA$Gene_Name%in%gene,]
dim(CGGA_gene)
rownames(CGGA_gene)<-CGGA_gene$Gene_Name
CGGA_gene<-CGGA_gene[,-1]
CGGA_gene<-t(CGGA_gene)
CGGA_gene<-rownames_to_column(as.data.frame(CGGA_gene),var="CGGA_ID")

survival_data<-read.csv("CGGA.mRNAseq_325_clinical.20200506.csv") #693
head(survival_data)
dim(survival_data)


survival_data_GBM<-subset(survival_data,Grade=="WHO IV"&PRS_type=="Primary")  
dim(survival_data_GBM)  #primary 85 samples
#survival_data_GBM<-subset(survival_data_GBM,OS>=90)

validate_data<-merge(survival_data_GBM[,c(1,7,8)],CGGA_gene,by="CGGA_ID")
colnames(validate_data)[c(2,3)]<-c("OS.time","OS")
validate_data<-na.omit(validate_data)

#cox1<-coxph(Surv(OS.time,OS)~HDAC4+ZNF22+LBH+HOXB2+EN2+AEBP1,data = validate_data)
cox1<-coxph(Surv(OS.time,OS)~BEST3 + CD180 + DUSP6 + G6PC3 + HEY1 + EN2 + HOXB2 + DLEU1 + IRX5 + LBH + ZEB1_AS1 + AGAP2_AS1 + POM121L9P,data = validate_data)

riskScore=predict(cox1,type = 'risk',data=validate_data)
risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))
validate_data$riskScore=riskScore
validate_data$risk=risk


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
pdf("Forestplot_CGGA325.pdf",width=7,height=6)
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
                           "16" = gpar(lwd=2, col="black")),  #"12"-nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),
           
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.5),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")         
dev.off()         











sfit <- survfit(Surv(OS.time, OS)~risk, data=validate_data)
ggsurv <- ggsurvplot(sfit, pval.method=T,surv.median.line ="hv",
                     cumevents=T,xlab = "Days",legend = "top",
                     title  = "combined-genes",risk.table.y.text = F,
                     conf.int=F,ncensor.plot = F,risk.table = F, 
                     pval=TRUE,data=validate_data)

pdf("validate CGGA325 survival curvenew.pdf",width=5,height=5)
ggsurv
dev.off()


pdf("test CGGA325 ROCnew1.pdf",width=5,height=5)
roc=survivalROC(Stime = validate_data$OS.time, status=validate_data$OS, 
                marker = validate_data$riskScore,predict.time =365, method="KM")
roc1=survivalROC(Stime = validate_data$OS.time, status=validate_data$OS, 
                 marker = validate_data$riskScore,predict.time =365*2, method="KM")
roc2=survivalROC(Stime = validate_data$OS.time, status=validate_data$OS, 
                 marker = validate_data$riskScore,predict.time =365*3, method="KM")

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


#phenotype
prog<-merge(validate_data[,c("CGGA_ID","risk")],survival_data,by="CGGA_ID")

table(prog$risk,prog$PRS_type)
chisq.test(table(prog$risk,prog$PRS_type))

table(prog$risk,prog$Gender)
chisq.test(table(prog$risk,prog$Gender))

table(prog$risk,prog$Radio_status..treated.1.un.treated.0.)
chisq.test(table(prog$risk,prog$Radio_status..treated.1.un.treated.0.))

table(prog$risk,prog$MGMTp_methylation_status)
chisq.test(table(prog$risk,prog$MGMTp_methylation_status))


table(prog$risk,prog$Histology)
chisq.test(table(prog$risk,prog$Histology))

table(prog$risk,prog$Grade)
chisq.test(table(prog$risk,prog$Grade))

table(prog$risk,prog$Chemo_status..TMZ.treated.1.un.treated.0.)
chisq.test(table(prog$risk,prog$Chemo_status..TMZ.treated.1.un.treated.0.))

chisq.test(table(prog$risk,prog$IDH_mutation_status))

table(prog$risk,prog$X1p19q_codeletion_status)
chisq.test(table(prog$risk,prog$X1p19q_codeletion_status))



table(prog$risk,prog$IDH_mutation_status)

#IDH_mutation_status
library(plyr)
a <- data.frame(table(prog$risk,prog$IDH_mutation_status))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a
pvalue <- chisq.test(table(prog$risk,prog$IDH_mutation_status))$p.value 

pdf("IDH_mutation_status325.pdf",width=4,height=5)
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack(),width=0.5)+
  scale_fill_manual(values = c("#DB423E","#008ECA"))+ #"gold"
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ 
  labs(x="Risk",y="Percent (%)",fill="IDH_mutation")+
  geom_text(aes(label=label,y=percent-2),
            position = position_stack(),size=6,color="black")+ #vjust=3
  annotate(geom = "text",
           cex=6,
           x=1.5, y=105, 
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", 
                                     paste0("= ",round(pvalue,3)))),
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
dev.off()
