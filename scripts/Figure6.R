
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
###umap rna
###################
immune.rna<-readRDS('cells/immune/immune_rna.rds')
data.rna<-data.frame(immune.rna@reductions$umap@cell.embeddings)
data.rna$celltype<-immune.rna$subtype
allcolour=c(colorRampPalette(brewer.pal(7,'Accent'))(10)[-2],'grey')
p<-ggplot(data.rna, aes(x = UMAP_1, 
                 y = UMAP_2, 
                 fill = celltype,
                 color = celltype)) +
  geom_point(size = 0.8) +
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = allcolour) +
  scale_color_manual(values = allcolour)
ggsave('f6/umap_rna.pdf',p,width = 8,height = 6)
#dotplot
immune.rna<-subset(immune.rna,subtype !='UD')
rna_marker=c('MS4A1','CD4','CD8A','APOE','KIT','FCN1','FCGR3A','JCHAIN','FOXP3')
rna_marker<-factor(rna_marker,levels = rna_marker)
p<-DotPlot(immune.rna,group.by = 'subtype',features = rna_marker,cols =c('lightgrey','red') )+ coord_flip() + RotatedAxis()
ggsave('f6/rna_dotplot.pdf',p,width = 10,height = 6)

###################
###umap atac
###################
immune.atac<-readRDS('cells/immune/immune_atac.rds')
data.atac<-data.frame(immune.atac@reductions$umap@cell.embeddings)
data.atac$celltype<-immune.atac$subtype
allcolour=c(colorRampPalette(brewer.pal(7,'Accent'))(10)[-2],'grey')
p<-ggplot(data.atac, aes(x = UMAP_1, 
                 y = UMAP_2, 
                 fill = celltype,
                 color = celltype)) +
  geom_point(size = 0.8) +
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = allcolour) +
  scale_color_manual(values = allcolour)
ggsave('f6/umap_atac.pdf',p,width = 8,height = 6)
#coverageplot
Idents(immune.atac)<-'subtype'
atac_marker=c('MS4A1','CD4','CD8A','APOE','KIT','FCN1','FCGR3A','JCHAIN','FOXP3')
for(i in atac_marker){
  p<-CoveragePlot(immune.atac,i,extend.upstream = 5000,extend.downstream = 5000)
  p<-p & scale_fill_manual(values =colorRampPalette(brewer.pal(7,'Accent'))(10)[-2] )
  ggsave(paste0('f6/',as.character(i),'.pdf'),p,width = 8,height = 6)
}

###################
###links barplot
###################
#load links data
ll<-readRDS('peakgenelinks/links_end.rds')
ll<-subset(ll,peak_promoter==1)
plotdf<-data.frame(matrix(0,nrow =length(unique(ll$gene_name)),ncol = 2 ))
colnames(plotdf)<-c('Gene','Linked_CRE')
plotdf$Gene<-unique(ll$gene_name)
for(i in 1:nrow(plotdf)){
  tmp<-subset(ll,gene_name==plotdf$Gene[i])
  plotdf[i,2]<-nrow(tmp)
}
#plot
library(ggplot2)
plotdata<-data.frame(table(plotdf$Linked_CRE))
plotdata$Var1<-as.numeric(plotdata$Var1)
plotdata$type<-ifelse(plotdata$Var1>=30,'30',plotdata$Var1)
plotdata$Freq[30:68]<-sum(plotdata$Freq[30:68])
plotdata<-plotdata[1:30,]
plotdata$type<-factor(plotdata$type,levels=c(1:30))
p<-ggplot(plotdata,aes(x=type,y=Freq))+geom_bar(stat='identity',fill='#000000')+
  theme_classic()+
  geom_vline(xintercept = 6,linetype='dashed',size=0.63)
ggsave('f6/link_barplot.pdf',height = 6,width = 8)

####################
###random background
####################
#load hic data
df_res<-readRDS('download/hic/res_hic.rds')
hic<-readRDS('download/hic/hic.rds')
colnames(df_res)<-colnames(hic)[12:28]
df_res<-df_res[,c(1:5,11,15,17)]
#plot
plotdf1<-data.frame(matrix(0,nrow = 100*8,ncol = 2))
colnames(plotdf1)<-c('celltype','percent')
plotdf1$celltype<-rep(colnames(df_res),each=100)
plotdf1$percent<-c(df_res[-1,]$Mon,df_res[-1,]$Mac0,df_res[-1,]$Mac1,df_res[-1,]$Mac2,
                   df_res[-1,]$Neu,df_res[-1,]$tCD4,df_res[-1,]$tCD8,df_res[-1,]$tB)
plotdf1$percent<-plotdf1$percent/53330
plotdf2<-data.frame(matrix(0,nrow = 8,ncol = 2))
colnames(plotdf2)<-c('celltype','percent')
plotdf2$celltype<-colnames(df_res)
plotdf2$percent<-as.numeric(df_res[1,])/53330
p<-ggplot(data=plotdf1,aes(x=celltype,y=percent))+
  geom_violin()+
  geom_point(data=plotdf2,aes(x=celltype,y=percent))+
  theme_bw()
ggsave('f6/links_background.pdf',p,height = 6,width = 8)

###################
###link plot
###################
#load links data
ll<-readRDS('peakgenelinks/links_end.rds')
ll<-subset(ll,peak_promoter==1)
ll_mrc1<-subset(ll,gene_name=='MRC1')
ll_trem2<-subset(ll,gene_name=='TREM2')
conns<-readRDS('cicero/immune/conns_600.rds')
conns<-subset(conns,abs(coaccess)>0.2)
peakanno<-readRDS('chipseeker/peakanno.rds')
peakanno<-subset(peakanno,type=='Promoter')
conns<-subset(conns,Peak2 %in% rownames(peakanno))
test.df <- data.frame(peakSite = conns$Peak2, 
                      peak2gene = peakanno[as.character(conns$Peak2),'SYMBOL'])
new_conns <- cbind(conns, test.df)
new_conns<-subset(new_conns,peak2gene %in% c('MRC1','TREM2'))
new_conns<-subset(new_conns,Peak1 %in% c(ll_trem2$peakName,ll_mrc1$peakName))
new_conns<-subset(new_conns,Peak2 %in% c('chr6-41161877-41162160','chr10-17808668-17809793'))
links <- ConnectionsToLinks(conns = new_conns)
Links(immune.atac) <- links
#plot
p1<-CoveragePlot(immune.atac,'TREM2',extend.upstream = 180000,extend.downstream = 150000)
p2<-CoveragePlot(immune.atac,'MRC1',extend.upstream = 70000,extend.downstream = 10000)
p1<-p1 & scale_fill_manual(values =colorRampPalette(brewer.pal(7,'Accent'))(10)[-2] )
p2<-p2 & scale_fill_manual(values =colorRampPalette(brewer.pal(7,'Accent'))(10)[-2] )
ggsave('f6/TREM2_coverageplot.pdf',p1,height = 6,width = 8)
ggsave('f6/MRC1_coverageplot.pdf',p2,height = 6,width = 8)