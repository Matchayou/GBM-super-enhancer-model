setwd('D:/SuperEnhancer/Fig4 external validation/GSE7696')

library(GEOquery)
library(dplyr)

gset = getGEO(GEO='GSE7696', destdir=".",getGPL = F)
probe=read.csv('GPL570-55999.txt',header = T,sep = '\t',comment.char = "#")

e=gset[[1]]
exp=exprs(e)
pheno=pData(e)
expr=data.frame(exp)
expr$ID=rownames(expr)

##????????????
ids=probe[,c(1,11)]
ids = ids[-grep('///',ids$Gene.Symbol),]## һ??̽????Ӧ??????????ȥ??
ids=ids[!ids$Gene.Symbol=='',]##̽??û?ж?Ӧ??????ȥ??
exprSet = inner_join(ids,expr,by = 'ID')
library(limma)
exprSet= avereps(exprSet[,-c(1,2)],           
                 ID = exprSet$Gene.Symbol)## ????̽????Ӧһ????????ȡ??ֵ
which(rownames(exprSet) == "ZEB1-AS1")
which(rownames(exprSet) == "AGAP2-AS1")
rownames(exprSet)[18486] <- "ZEB1_AS1"
rownames(exprSet)[1988] <- "AGAP2_AS1"


#boxplot(exprSet)##????ͼ??һ?????ݷֲ?
save(exprSet,file = 'expr.Rdata')
save(ids,file = 'GPL570.Rdata')


##?????ٴ???Ϣ

type=pheno$characteristics_ch1.1
type=strsplit(type,split = ": ",fixed = T)
type=sapply(type,function(x){x[2]})

age=pheno$characteristics_ch1.3
age=strsplit(age,split = ": ",fixed = T)
age=sapply(age,function(x){x[2]})

gender=pheno$characteristics_ch1.4
gender=strsplit(gender,split = ": ",fixed = T)
gender=sapply(gender,function(x){x[2]})

treat=pheno$characteristics_ch1.5
treat=strsplit(treat,split = ": ",fixed = T)
treat=sapply(treat,function(x){x[2]})

OS=pheno$characteristics_ch1.6
OS=strsplit(OS,split = ": ",fixed = T)
OS=sapply(OS,function(x){x[2]})

OS.time=pheno$characteristics_ch1.7
OS.time=strsplit(OS.time,split = ": ",fixed = T)
OS.time=sapply(OS.time,function(x){x[2]})

pdata=data.frame(cbind(OS,OS.time,age,gender,type,treat))
rownames(pdata)=rownames(pheno)

#pdata<-subset(pdata,type=="GBM")
save(pdata,file = 'pdata.Rdata')




load('expr.Rdata')
load("pdata.Rdata")
###????ģ????֤

#TF and co-factor
gene=c("BEST3","CD180","DUSP6","G6PC3","HEY1","EN2","HOXB2","DLEU1","IRX5","LBH","ZEB1_AS1","LINC01265","AGAP2_AS1","POM121L9P")

exp=exprSet[rownames(exprSet)%in%gene,]
exp=data.frame(exp)
exp=t(exp)

sur=subset(pdata,!pdata$OS=='NA')
exp=data.frame(exp)
exp=exp[rownames(sur),]
data=cbind(sur,exp)

data$OS.time=as.numeric(data$OS.time)
data$OS=as.numeric(data$OS)
data$OS.time=30*data$OS.time

data=subset(data,type=="GBM")

#TF model
cox<-coxph(Surv(OS.time,OS)~BEST3+CD180+DUSP6+G6PC3+HEY1+EN2+HOXB2+DLEU1+IRX5+LBH+ZEB1_AS1+AGAP2_AS1+POM121L9P,data = data)

riskScore=predict(cox,type = 'risk',data=data)
risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))  #median
data$riskScore=riskScore
data$risk=risk

validate_data=data
sfit <- survfit(Surv(OS.time, OS)~risk, data=validate_data)
ggsurv <- ggsurvplot(sfit, pval.method=T,surv.median.line ="hv",
                     cumevents=T,xlab = "Days",legend = "top",
                     title  = "combined-genes",risk.table.y.text = F,
                     conf.int=F,ncensor.plot = F,risk.table = F, 
                     pval=TRUE,data=validate_data)
pdf("validate GSE7696 survival curvenew.pdf",width=5,height=5)
ggsurv
dev.off()


#ROC????
pdf("test GSE7696 ROCnew1.pdf",width=5,height=5)
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

validate_data$GEA_ID <- rownames(validate_data)
sur$GEA_ID <- rownames(sur)
prog<-merge(validate_data[,c("GEA_ID","risk")],sur,by="GEA_ID")
#IDH_mutation_status
library(plyr)
a <- data.frame(table(prog$risk,prog$IDH_mutation_status))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a
pvalue <- chisq.test(table(prog$risk,prog$IDH_mutation_status))$p.value 

#pdf("IDH_mutation_status325.pdf",width=4,height=5)
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
#dev.off()












