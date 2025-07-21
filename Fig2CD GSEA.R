library(clusterProfiler)
library(enrichplot)
library(GOplot)
library(DOSE)

setwd("D:/SuperEnhancer/Fig2 TCGA TF overlap")
SEDGs<-read.csv("1154 SEs-related.csv")
dim(SEDGs)  

##GO??KEGG analysis
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

GO_SEDGs_CC<-enrichGO(SEDGs$ENTREZID,
                   OrgDb = 'org.Hs.eg.db',
                   keyType = "ENTREZID",
                   ont = "CC",pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",qvalueCutoff = 0.01,
                   minGSSize = 10,maxGSSize = 500,readable = TRUE)
dim(as.data.frame(GO_SEDGs_CC))
write.csv(as.data.frame(GO_SEDGs_CC),file="GO CC terms.csv")

GO_SEDGs_CC <- read.csv("GO CC terms.csv")
pdf("GO CC top10new.pdf",height=5,width=4.5)
dotplot(GO_SEDGs_CC, showCategory = 10,title="top 10 GO CC")+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=25))
dev.off()


pdf("GO BP top10new.pdf",height=5,width=4.5)
dotplot(GO_SEDGs, showCategory = 10,title="top 10 GO BP")+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=25))
dev.off()

GO_SEDGs_MF<-enrichGO(SEDGs$ENTREZID,
                      OrgDb = 'org.Hs.eg.db',
                      keyType = "ENTREZID",
                      ont = "MF",pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",qvalueCutoff = 0.01,
                      minGSSize = 10,maxGSSize = 500,readable = TRUE)
dim(as.data.frame(GO_SEDGs_MF))
write.csv(as.data.frame(GO_SEDGs_MF),file="GO MF terms.csv")
pdf("GO MF top10new.pdf",height=5,width=4.5)
dotplot(GO_SEDGs_MF, showCategory = 10,title="top 10 GO MF")+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=25))
dev.off()


##or using codes
ego <- enrichGO(gene         = SEDGs$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

#cutoff=0.6
ego_sim<-simplify(GO_SEDGs, cutoff=0.4,by="p.adjust",select_fun=min) 
dim(as.data.frame(ego_sim))

pdf("GO terms simplify0.4new.pdf",height=5.5,width=8)
dotplot(ego_sim,title = 'GO BP',showCategory = 20)+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=70))
dev.off()
write.csv(as.data.frame(GO_SEDGs),file="GO slim cutoff termsnew.csv")




#KEGG pathway
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
KEGG_SEDGs<-enrichKEGG(SEDGs$ENTREZID,organism = "hsa",keyType = "kegg",
                      pvalueCutoff = 0.05,pAdjustMethod = "BH", minGSSize = 5,
                      maxGSSize = 500,qvalueCutoff = 0.01,use_internal_data = FALSE)


pdf("KEGG pathwaynew.pdf",height=8,width=7)
dotplot(KEGG_SEDGs,title = 'KEGG',showCategory = 39)+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=70))
dev.off()

write.csv(as.data.frame(KEGG_SEDGs),file="KEGG pathwaysnew.csv")

browseKEGG(KEGG_SEDGs, 'hsa05214')





x <- enrichDO(gene          = SEDGs$ENTREZID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)
dim(x)
write.csv(x,file="disease ontologynew.csv")
