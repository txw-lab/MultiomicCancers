
setwd('data/pancancer_atac/tumor/combined')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)
library(RColorBrewer)

###################
###diff peak
##################
#load data
fibro<-readRDS('cells/fibro/fibro_atac.rds')
meta<-fibro@meta.data
peak_fibro<-readRDS('cells/fibro/peaks_fibro.rds')
#plot
library(pheatmap)
peakdata<-fibro@assays$peaks@data
peakdata<-peakdata[unique(peak_fibro$gene),]
plotdata<-data.frame(matrix(0,nrow = nrow(peakdata),ncol = length(unique(peak_fibro$cluster))))
rownames(plotdata)<-rownames(peakdata)
colnames(plotdata)<-unique(peak_fibro$cluster)
for(i in colnames(plotdata)){ 
  tmp<-subset(meta,type==i)
  plotdata[,i]<-rowMeans(peakdata[,rownames(tmp)])
}
p<-pheatmap(t(scale(t(plotdata))),show_rownames = F,cluster_rows = F,cluster_cols = T)
ggsave('f5/peak_heatmap.pdf',p,width = 4,height = 4)

###################
###great
###################
library(ggplot2)
great_cc<-read.delim('cells/fibro/diffpeak/greatExportAll_cc.tsv',skip = 3)
great_ec<-read.delim('cells/fibro/diffpeak/greatExportAll_ec.tsv',skip = 3)
great_oc<-read.delim('cells/fibro/diffpeak/greatExportAll_oc.tsv',skip = 3)
p<-ggplot(data = great_ec[1:10,], 
       aes(x = Desc, y = ObsRegions,fill = BinomBonfP))+
  geom_bar(stat = "identity",width = 0.9)+ 
  coord_flip()+theme_bw()+ 
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ 
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ 
  theme(axis.title = element_text(size = 13), 
        axis.text = element_text(size = 11), 
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 11), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) 
ggsave('f5/go_barplot.pdf',p,width = 4,height = 4)

###################
###coverage plot
###################
fibro<-readRDS('cells/fibro/fibro_atac.rds')
Idents(fibro)<-'type'
p1<-CoveragePlot(fibro,'TNFSF11',extend.upstream = 5000,extend.downstream = 5000)
p2<-CoveragePlot(fibro,'TNFSF18',extend.upstream = 5000,extend.downstream = 5000)
p3<-CoveragePlot(fibro,'JAG2',extend.upstream = 5000,extend.downstream = 5000)
p4<-CoveragePlot(fibro,'POSTN',extend.upstream = 5000,extend.downstream = 5000)
ggsave('coverageplot.pdf',p1+p2+p3+p4,height = 8,width = 8)

###################
###footprint
###################
fibro<-readRDS('cells/fibro/fibro_atac.rds')
Idents(fibro)<-'type'
p<-PlotFootprint(fibro,features=c('MA1638.1','MA0878.1','MA0152.1'))
ggsave('footprint.pdf',p,height = 8,width = 8)
