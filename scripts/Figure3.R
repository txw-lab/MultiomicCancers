
setwd('data/pancancer_atac/tumor/combined/')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggunchained)

#############
###cc
#############
#load data
cc<-readRDS('cc/chromvar/cc.rds')
cc_motifs<-c('ASCL1(var.2)','CDX1','FOXA1','GATA2','GRHL1','HNF1A','HNF4A','HOXA10','LEF1','TCF7','TEAD4',
             'ATF4','TWIST1','ZNF24','ZBTB26','ZBTB18',
             'SOX2','SOX8',
             'CEBPA','EHF','ELF3','IKZF1','SPI1','SPIB','SPIC','ZKSCAN5',
             'ELK4','RUNX3',
             'MSC','POU1F1','POU2F1',
             'PAX9'
)
cc_data<-cc@assays$chromvar@data
rownames(cc_data)<-cc@assays$peaks@motifs@motif.names[rownames(cc_data)]
cc_data<-cc_data[unique(cc_motifs),]
cc_df<-data.frame(matrix(0,nrow = nrow(cc_data),ncol = 7))
colnames(cc_df)<-c('Tumor cell','Fibroblast','Endothelial','Myeloid cell','T cell','Plasma','B cell')
rownames(cc_df)<-rownames(cc_data)
meta<-cc@meta.data
for(i in 1:7){
  tmp<-subset(meta,celltype==colnames(cc_df)[i])
  cc_df[,i]<-rowMeans(cc_data[,rownames(tmp)])
}
#plot
library(pheatmap)
p<-pheatmap(t(scale(t(cc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
         border_color = NA)
ggsave('f3/cc_heatmap.pdf',p,height = 7,width = 4)

##############
###ec
##############
#load data
ec<-readRDS('ec/chromvar/ec.rds')
ec_motifs<-c('ASCL1(var.2)','FOXA1','GRHL1','HNF1A','HNF4A','HOXA1','LEF1','NFIC','SOX2','TCF7','TEAD4','TFAP2A',
             'ARGFX','CDX1','HAND2','HOXA10','TWIST1','ZBTB18',
             'EBF1','EBF3','PBX2',
             'STAT1','STAT3',
             'CEBPA','EHF','ELF3','IKZF1','SPI1','SPIB','SPIC','ZKSCAN5',
             'ELK4','EOMES','RUNX2','RUNX3',
             'IRF4','POU1F1',
             'GATA1','GATA1::TAL1'
)
ec_data<-ec@assays$chromvar@data
rownames(ec_data)<-ec@assays$peaks@motifs@motif.names[rownames(ec_data)]
ec_data<-ec_data[unique(ec_motifs),]
ec_df<-data.frame(matrix(0,nrow = nrow(ec_data),ncol = 8))
colnames(ec_df)<-c('Tumor cell','Fibroblast','Myofibroblast','Endothelial','Myeloid  cell','T cell','B cell','Mast')
rownames(ec_df)<-rownames(ec_data)
meta<-ec@meta.data
for(i in 1:8){
  tmp<-subset(meta,celltype==colnames(ec_df)[i])
  ec_df[,i]<-rowMeans(ec_data[,rownames(tmp)])
}
#plot
library(pheatmap)
p<-pheatmap(t(scale(t(ec_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
         border_color = NA)
ggsave('f3/ec_heatmap.pdf',p,height = 7,width = 4)

###############
###oc
###############
#load data
oc<-readRDS('oc/chromvar/oc.rds')
oc_motifs<-c('ASCL1(var.2)','GRHL1','NFIC','PAX1','SIX1','SOX2','TEAD4','TFAP2A','ZEB1',
             'ARGFX','HAND2','TWIST1','ZBTB18',
             'EBF1','EBF3','NR3C1','NR4A1',
             'STAT1','STAT3',
             'CEBPA','EHF','ELF3','SPI1','SPIB','SPIC',
             'ELK4','EOMES','ETV1','RUNX2','RUNX3',
             'IRF4','POU1F1'
)
oc_data<-oc@assays$chromvar@data
rownames(oc_data)<-oc@assays$peaks@motifs@motif.names[rownames(oc_data)]
oc_data<-oc_data[unique(oc_motifs),]
oc_df<-data.frame(matrix(0,nrow = nrow(oc_data),ncol = 7))
colnames(oc_df)<-c('Tumor cell','Fibroblast','Myofibroblast','Endothelial','Myeloid cell','T cell','Plasma')
rownames(oc_df)<-rownames(oc_data)
meta<-oc@meta.data
for(i in 1:7){
  tmp<-subset(meta,celltype==colnames(oc_df)[i])
  oc_df[,i]<-rowMeans(oc_data[,rownames(tmp)])
}
#plot
library(pheatmap)
p<-pheatmap(t(scale(t(oc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
         border_color = NA)
ggsave('f3/oc_heatmap.pdf',p,height = 7,width = 4)

###################
###rcc
###################
#load data
rcc<-readRDS('rcc/chromvar/rcc.rds')
rcc_motifs<-c('HNF1A','HNF4A','NFIC','PAX3','POU4F1','TEAD4','VAX1','ZNF24','FOXA1',
              'EBF1','EBF3','MEF2A','MEF2B',
              'STAT1','STAT3','SOX2','SOX4',
              'ATF4','CEBPA','EHF','ELF3','SPI1','SPIB','SPIC',
              'EOMES','RUNX2','RUNX3','TBX1',
              'ASCL1','FIGLA','IRF4',
              'SNAI2','POU1F1','POU2F1',
              'GATA1','GATA1::TAL1'
)
rcc_data<-rcc@assays$chromvar@data
rownames(rcc_data)<-rcc@assays$peaks@motifs@motif.names[rownames(rcc_data)]
rcc_data<-rcc_data[unique(rcc_motifs),]
rcc_df<-data.frame(matrix(0,nrow = nrow(rcc_data),ncol = 8))
colnames(rcc_df)<-c('Tumor cell','Myofibroblast','Endothelial','Myeloid cell','T cell','Plasma','B cell','Mast')
rownames(rcc_df)<-rownames(rcc_data)
meta<-rcc@meta.data
for(i in 1:8){
  tmp<-subset(meta,celltype==colnames(rcc_df)[i])
  rcc_df[,i]<-rowMeans(rcc_data[,rownames(tmp)])
}
#plot
library(pheatmap)
p<-pheatmap(t(scale(t(rcc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
         border_color = NA)
ggsave('f3/rcc_heatmap.pdf',p,height = 7,width = 4)

#############
###atac umap
#############
#load data
colon<-readRDS('../../../colondata/schtan/atac/colon.rds')
Idents(colon)<-'celltype'
allcolour=c("#FEFB98","#C4AFCB",'#C91889',"#96BF9E","#FDCD8A",'#7FC97F',"#DBB6AF",'#D63048')
umap_cc = colon@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = colon$celltype) 
#plot
p<-ggplot(umap_cc,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 0.8 , alpha =1 )  +  
  scale_color_manual(values = allcolour)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=20), 
    legend.key.size=unit(1,'cm') ) +  
  guides(color = guide_legend(override.aes = list(size=5)))+ 
  geom_segment(aes(x = min(umap_cc$UMAP_1) , y = min(umap_cc$UMAP_2) ,
                   xend = min(umap_cc$UMAP_1) +3, yend = min(umap_cc$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap_cc$UMAP_1)  , y = min(umap_cc$UMAP_2)  ,
                   xend = min(umap_cc$UMAP_1) , yend = min(umap_cc$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap_cc$UMAP_1) +1.5, y = min(umap_cc$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap_cc$UMAP_1) -1, y = min(umap_cc$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
ggsave('f3/colon_atac_umap.pdf',p,height = 6,width = 8)

##########################
##coverageplot
##########################
#load data
colon<-readRDS('../../../colondata/schtan/atac/colon.rds')
colon$celltype<-factor(colon$celltype,levels = sort(unique(colon$celltype)))
Idents(colon)<-'celltype'
#plot
genes<-c('MS4A1','PECAM1',"EPCAM",'PDGFRA','PLP1','ITGAX','ACTA2','NRG1','CD247','LGR5')
allcolour=c("#FEFB98","#C4AFCB",'#C91889',"#96BF9E","#FDCD8A",'#7FC97F',"#DBB6AF",'#D63048')
for(i in genes){
  p<-CoveragePlot(colon,i,extend.upstream = 5000,extend.downstream = 5000)
  p<-p & scale_fill_manual(values = allcolour)
  ggsave(paste0('f3/',i,'.pdf'),p,width = 6,height = 4)
}

################
##footprint
################
epi<-readRDS('../../../colondata/tf/tumor/epi.rds')
p1<-PlotFootprint(epi, features = c('MA0867.2'))
ggsave('f3/sox4_footprint.pdf',p1,width = 6,height = 5)

p2<-PlotFootprint(epi, features = c('MA0768.1'))
ggsave('f3/lef1_footprint.pdf',p2,width =6,height = 5)

p3<-PlotFootprint(epi, features = c('MA0769.2'))
ggsave('f3/tcf7_footprint.pdf',p3,width = 6,height = 5)

p4<-PlotFootprint(epi, features = c('MA0809.2'))
ggsave('f3/tead4_footprint.pdf',p4,width = 6,height = 5)

##############
###chromvar
##############
#load data
epi<-readRDS('../../../colondata/chromvar/all/colon.rds')
epi<-subset(epi,celltype %in% c('Tumor cell','Epithelial'))
#plot
data<-epi@assays$chromvar@data[c('MA0867.2','MA0768.1',"MA0769.2",'MA0809.2'),]
rownames(data)<-c('SOX4','LEF1','TCF7','TEAD4')
epi$SOX4<-data['SOX4',]
epi$LEF1<-data['LEF1',]
epi$TCF7<-data['TCF7',]
epi$TEAD4<-data['TEAD4',]
meta<-epi@meta.data
df_sox4<-meta[,13:14]
df_sox4$TF<-'SOX4'
colnames(df_sox4)[2]<-'Chromvar_score'
df_lef1<-meta[,c(13,15)]
df_lef1$TF<-'LEF1'
colnames(df_lef1)[2]<-'Chromvar_score'
df_tcf7<-meta[,c(13,16)]
df_tcf7$TF<-'TCF7'
colnames(df_tcf7)[2]<-'Chromvar_score'
df_tead4<-meta[,c(13,17)]
df_tead4$TF<-'TEAD4'
colnames(df_tead4)[2]<-'Chromvar_score'
df<-rbind(df_lef1,df_tcf7,df_tead4,df_sox4)
df$TF<-factor(df$TF,levels = sort(unique(df$TF)))
p<-ggplot(df,aes(x = TF,y =Chromvar_score,fill=celltype)) +
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
        legend.position = 'top') +
  scale_fill_brewer(palette = 'Set1') +
  scale_fill_jco(name = '') +
  theme(panel.grid=element_blank())+
  stat_compare_means(aes(group=celltype),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif")
ggsave('f3/chromvar_tf.pdf',p,width = 8,height = 4)

################
###rna
################
#load data
colon<-readRDS('../../../colondata/cells/all/colon.rds')
epi<-subset(colon,celltype %in% c('Tumor cell','Epithelial'))
#plot
data<-epi@assays$RNA@counts[c('SOX4','LEF1','TCF7','TEAD4'),]
data<-t(scale(t(data)))
epi$SOX4<-data['SOX4',]
epi$LEF1<-data['LEF1',]
epi$TCF7<-data['TCF7',]
epi$TEAD4<-data['TEAD4',]
meta<-epi@meta.data
df_sox4<-meta[,34:35]
df_sox4$TF<-'SOX4'
colnames(df_sox4)[2]<-'Expression'
df_lef1<-meta[,c(34,36)]
df_lef1$TF<-'LEF1'
colnames(df_lef1)[2]<-'Expression'
df_tcf7<-meta[,c(34,37)]
df_tcf7$TF<-'TCF7'
colnames(df_tcf7)[2]<-'Expression'
df_tead4<-meta[,c(34,38)]
df_tead4$TF<-'TEAD4'
colnames(df_tead4)[2]<-'Expression'
df<-rbind(df_lef1,df_tcf7,df_tead4,df_sox4)
df$TF<-factor(df$TF,levels = sort(unique(df$TF)))
p<-ggplot(df,aes(x = TF,y =Expression,fill=celltype)) +
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
        legend.position = 'top') +
  scale_fill_brewer(palette = 'Set1') +
  scale_fill_jco(name = '') +
  theme(panel.grid=element_blank())+
  stat_compare_means(aes(group=celltype),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif")
ggsave('f3/rna_tf.pdf',p,width = 8,height = 4)