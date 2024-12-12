
setwd('tumor')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)
library(RColorBrewer)

###################
###tf heatmap
###################
#cc
cc<-readRDS('cc/chromvar/cc.rds')
cc_motifs<-c('ASCL1(var.2)','CDX1','FOXA1','GATA2','GRHL1','HNF1A','HNF4A','HOXA10','LEF1','TCF7','ZEB1',
             'TWIST1','TEAD4','NFATC2','ZBTB26','ZBTB18',
             'SOX2','SOX8',
             'CEBPA','EHF','ELF3','IKZF1','SPI1','SPIB','SPIC',
             'ELK4','RUNX3','ETV1','NFYA',
             'IRF2','POU1F1','POU2F1')
cc_data<-cc@assays$chromvar@data
rownames(cc_data)<-cc@assays$peaks@motifs@motif.names[rownames(cc_data)]
cc_data<-cc_data[unique(cc_motifs),]
cc_df<-data.frame(matrix(0,nrow = nrow(cc_data),ncol = 6))
colnames(cc_df)<-c('Tumor cell','Fibroblast','Endothelial','Myeloid cell','T cell','B cell')
rownames(cc_df)<-rownames(cc_data)
meta<-cc@meta.data
for(i in 1:6){
  tmp<-subset(meta,celltype==colnames(cc_df)[i])
  cc_df[,i]<-rowMeans(cc_data[,rownames(tmp)])
}
library(pheatmap)
p<-pheatmap(t(scale(t(cc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/cc_tf.pdf',p,height = 6,width = 2.5)

#ec
ec<-readRDS('ec/chromvar/ec.rds')
ec_motifs<-c('ASCL1(var.2)','FOXA1','GRHL1','HNF1A','HNF4A','HOXA1','LEF1','NFIC','SOX2','TCF7','TEAD4','TFAP2A',
             'ARGFX','CDX1','HAND2','HOXA10','TWIST1','ZBTB18',
             'EBF1','EBF3','PBX2',
             'STAT1','STAT3',
             'CEBPA','EHF','ELF3','IKZF1','SPI1','SPIB','SPIC','ZKSCAN5',
             'ELK4','EOMES','RUNX2','RUNX3',
             'IRF4','POU1F1',
             'GATA1','GATA1::TAL1')
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
library(pheatmap)
p<-pheatmap(t(scale(t(ec_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/ec_tf.pdf',p,height = 6,width = 2.5)

#oc
oc<-readRDS('oc/chromvar/oc.rds')
oc_motifs<-c('ASCL1(var.2)','GRHL1','NFIC','PAX1','SIX1','SOX2','TEAD4','TFAP2A','ZEB1',
             'ARGFX','HAND2','TWIST1','ZBTB18',
             'EBF1','EBF3','NR3C1','NR4A1',
             'STAT1','STAT3',
             'CEBPA','EHF','ELF3','SPI1','SPIB','SPIC',
             'ELK4','EOMES','ETV1','RUNX2','RUNX3',
             'IRF4','POU1F1')
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
library(pheatmap)
p<-pheatmap(t(scale(t(oc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/oc_tf.pdf',p,height = 6,width = 2.5)

#rcc
rcc<-readRDS('rcc/chromvar/rcc.rds')
rcc_motifs<-c('HNF1A','HNF4A','NFIA','PAX3','POU4F1','TEAD4','VAX1','ZNF24','ARGFX',
              'EBF1','EBF3',
              'STAT1','STAT3','SOX2','SOX4',
              'CEBPA','EHF','ELF3','SPI1','SPIB','SPIC',
              'EOMES','RUNX2','RUNX3','TBX1',
              'ASCL1(var.2)','FIGLA','IRF4','MYOD1','SNAI2','POU1F1','POU2F1')
rcc_data<-rcc@assays$chromvar@data
rownames(rcc_data)<-rcc@assays$peaks@motifs@motif.names[rownames(rcc_data)]
rcc_data<-rcc_data[unique(rcc_motifs),]
rcc_df<-data.frame(matrix(0,nrow = nrow(rcc_data),ncol = 6))
colnames(rcc_df)<-c('Tumor cell','Myofibroblast','Endothelial','Myeloid cell','T cell','Plasma')
rownames(rcc_df)<-rownames(rcc_data)
meta<-rcc@meta.data
for(i in 1:6){
  tmp<-subset(meta,celltype==colnames(rcc_df)[i])
  rcc_df[,i]<-rowMeans(rcc_data[,rownames(tmp)])
}
library(pheatmap)
p<-pheatmap(t(scale(t(rcc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/rcc_tf.pdf',p,height = 6,width = 2.5)

#plc
plc<-readRDS('plc/chromvar/plc.rds')
plc_motifs<-c('CUX1','FOXA1','GATA2','HNF1A','HNF4A','NEUROD1','PAX3','TCF7L1','TEAD4',
              'EBF1','EBF3',
              'SOX2','SOX8','SOX9',
              'CEBPA','EHF','ELF3','ETV4','IKZF1','SPI1','SPIB','SPIC',
              'EOMES','RUNX2','RUNX3',
              'ASCL1(var.2)','IRF2','MYF5','SNAI1','TCF3',
              'POU1F1','POU2F1','POU5F1','POU5F1B')
plc_data<-plc@assays$chromvar@data
rownames(plc_data)<-plc@assays$peaks@motifs@motif.names[rownames(plc_data)]
plc_data<-plc_data[unique(plc_motifs),]
plc_df<-data.frame(matrix(0,nrow = nrow(plc_data),ncol = 7))
colnames(plc_df)<-c('Tumor cell','Myofibroblast','Endothelial','Myeloid cell','T cell','Plasma','B cell')
rownames(plc_df)<-rownames(plc_data)
meta<-plc@meta.data
for(i in 1:7){
  tmp<-subset(meta,celltype==colnames(plc_df)[i])
  plc_df[,i]<-rowMeans(plc_data[,rownames(tmp)])
}
library(pheatmap)
p<-pheatmap(t(scale(t(plc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/plc_tf.pdf',p,height=6,width=2.5)

#bcc
bcc<-readRDS('bcc/chromvar/bcc.rds')
bcc_motifs<-c('ALX3','EMX1','GRHL1','HOXD3','LBX1','MEOX1','NFIA','ZNF24','PHOX2A','TEAD4','TFAP2A','TP53','VAX1',
              'CEBPA','EBF3','NFATC2','POU4F1','TWIST1',
              'EHF','ELF1','IKZF1','SPI1','SPIB','SPIC',
              'EOMES','MGA','TBR1','TBX1',
              'POU1F1','POU2F1','POU3F1','POU5F1',
              'IRF1')
bcc_data<-bcc@assays$chromvar@data
rownames(bcc_data)<-bcc@assays$peaks@motifs@motif.names[rownames(bcc_data)]
bcc_data<-bcc_data[unique(bcc_motifs),]
bcc_df<-data.frame(matrix(0,nrow = nrow(bcc_data),ncol = 9))
colnames(bcc_df)<-c('Tumor cell','Myofibroblast','Endothelial','Myeloid cell','NK','T cell','Plasma','B cell','Melanocyte')
rownames(bcc_df)<-rownames(bcc_data)
meta<-bcc@meta.data
for(i in 1:9){
  tmp<-subset(meta,celltype==colnames(bcc_df)[i])
  bcc_df[,i]<-rowMeans(bcc_data[,rownames(tmp)])
}
library(pheatmap)
p<-pheatmap(t(scale(t(bcc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/bcc_tf.pdf',p,height=6,width=2.5)

#lc
lc<-readRDS('lc/chromvar/lc.rds')
lc_motifs<-c('ASCL1(var.2)','FIGLA','GRHL1','FOXA1','HNF1A','MYOD1','NFIA','ZNF24','NKX2-2','TP53','ZEB1',
             'TEAD1','TEAD2','TEAD3','TEAD4','RUNX2',
             'SOX2','SOX8','SOX9','SOX13','SOX14',
             'CEBPA','ETV1','EHF','ELF1','IKZF1','SPI1','SPIB','SPIC',
             'IRF1','POU1F1','POU2F1','POU3F1','POU5F1')
lc_data<-lc@assays$chromvar@data
rownames(lc_data)<-lc@assays$peaks@motifs@motif.names[rownames(lc_data)]
lc_data<-lc_data[unique(lc_motifs),]
lc_df<-data.frame(matrix(0,nrow = nrow(lc_data),ncol = 8))
colnames(lc_df)<-c('Tumor cell','Fibroblast','Myofibroblast','Endothelial','Myeloid cell','Plasma','B cell','T cell')
rownames(lc_df)<-rownames(lc_data)
meta<-lc@meta.data
for(i in 1:8){
  tmp<-subset(meta,celltype==colnames(lc_df)[i])
  lc_df[,i]<-rowMeans(lc_data[,rownames(tmp)])
}
library(pheatmap)
p<-pheatmap(t(scale(t(lc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/lc_tf.pdf',p,height = 6,width = 2.5)

#bc
bc<-readRDS('bc/chromvar/bc.rds')
bc_motifs<-c('ASCL1(var.2)','FOXA1','GRHL1','MYOD1','SCRT1','TCF3','TFCP2','ZEB1',
             'NFATC2','TEAD4','TWIST1',
             'EBF1','SOX8','TP53','NFIA',
             'SOX8','SOX9','SOX15',
             'CEBPA','EHF','ELF1','ETV1','IKZF1','SPI1','USF1','ZKSCAN5',
             'IRF2','IRF4','IRF4','IRF8')
bc_data<-bc@assays$chromvar@data
rownames(bc_data)<-bc@assays$peaks@motifs@motif.names[rownames(bc_data)]
bc_data<-bc_data[unique(bc_motifs),]
bc_df<-data.frame(matrix(0,nrow = nrow(bc_data),ncol = 6))
colnames(bc_df)<-c('Tumor cell','Fibroblast','Myofibroblast','Endothelial','Myeloid cell','Lymphocyte')
rownames(bc_df)<-rownames(bc_data)
meta<-bc@meta.data
for(i in 1:6){
  tmp<-subset(meta,celltype==colnames(bc_df)[i])
  bc_df[,i]<-rowMeans(bc_data[,rownames(tmp)])
}
library(pheatmap)
p<-pheatmap(t(scale(t(bc_df))),show_colnames = T,cluster_rows = F,cluster_cols = F,
            border_color = NA)
ggsave('f5/bc_tf.pdf',p,height=6,width=2.5)

##################
###motif umap
##################
cc<-readRDS('cc/chromvar/cc.rds')
ec<-readRDS('ec/chromvar/ec.rds')
oc<-readRDS('oc/chromvar/oc.rds')
rcc<-readRDS('rcc/chromvar/rcc.rds')
bcc<-readRDS('bcc/chromvar/bcc.rds')
plc<-readRDS('plc/chromvar/plc.rds')
Idents(cc)<-'celltype'
Idents(ec)<-'celltype'
Idents(oc)<-'celltype'
Idents(rcc)<-'celltype'
Idents(bcc)<-'celltype'
Idents(plc)<-'celltype'
Idents(lc)<-'celltype'
Idents(bc)<-'celltype'
umap_cc = cc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = cc$celltype) 
umap_ec = ec@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = ec$celltype) 
umap_oc = oc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = oc$celltype) 
umap_rcc = rcc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = rcc$celltype) 
umap_bcc = bcc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = bcc$celltype) 
umap_plc = plc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = plc$celltype)
umap_bc = bc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = bc$celltype) 
umap_lc = lc@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = lc$celltype) 
allcolour_cc=c("#FEFB98","#C4AFCB","#96BF9E","#FDCD8A","#FEE491","#DBB6AF",'#D63048')
allcolour_ec=c("#FEFB98","#C4AFCB","#96BF9E",'#A75D2B',"#FDCD8A","#7FC97F","#DBB6AF",'#D63048')
allcolour_oc=c("#C4AFCB","#96BF9E","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_rcc=c("#FEFB98","#C4AFCB","#A75D2B","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_bcc=c("#FEFB98","#C4AFCB","#96BF9E","#A75D2B","#769AA8","#FDCD8A","#7FC97F","#C91889",'#FEE491',"#DBB6AF",'#D63048')
allcolour_plc=c("#FEFB98","#C4AFCB","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_lc=c("#FEFB98","#C4AFCB",'#96BF9E',"#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_bc=c("#C4AFCB",'#96BF9E','#F3BD92',"#FDCD8A","#7FC97F",'#D63048')

#plot
for(umap in c(umap_cc,umap_ec,umap_oc,umap_rcc,umap_plc,umap_bcc,umap_lc,umap_bc)){
  p<-ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
    geom_point(size = 1 , alpha =1 )  +  
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
    geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                     xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
    geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                     xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
             color="black",size = 3, fontface="bold" ) + 
    annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
             color="black",size = 3, fontface="bold" ,angle=90) +
    theme(legend.position = "none")
  ggsave('f5/Umap.atac.pdf',p,height = 6,width = 8)
}

###################
###tead rna
###################
plotdata<-data.frame(TF=c('TEAD1','TEAD2','TEAD3','TEAD4'),T_Value=c(99.321,30.383,47.559,34.272))
p1<-ggplot(data = plotdata, aes(x = TF, y = T_Value,fill=TF)) +
  geom_bar(stat = 'identity', width = .5,position = position_dodge(0.85))+
  theme_test()+
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme(legend.position = 'none')+
  labs(title = 'Gene Expression')
ggsave('f5/tead_rna.pdf',p1,width = 5,height = 2.5)

###################
###tead chromvar
###################
plotdata<-data.frame(TF=c('TEAD1','TEAD2','TEAD3','TEAD4'),T_Value=c(227.87,207.77,213.08,219.45))
p2<-ggplot(data = plotdata, aes(x = TF, y = T_Value,fill=TF)) +
  geom_bar(stat = 'identity', width = .5,position = position_dodge(0.85))+
  theme_test()+
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme(legend.position = 'none')+
  labs(title = 'Chromvar Score')
ggsave('f5/tead_chromvar.pdf',p2,width = 5,height = 2.5)

###################
###tead pathway
###################
bcc_res<-readRDS('combined/tf/all/bcc_res.rds')
cc_res<-readRDS('combined/tf/all/cc_res.rds')
ec_res<-readRDS('combined/tf/all/ec_res.rds')
oc_res<-readRDS('combined/tf/all/oc_res.rds')
plc_res<-readRDS('combined/tf/all/plc_res.rds')
rcc_res<-readRDS('combinedtf/all/rcc_res.rds')
lc_res<-readRDS('combinedtf/all/lc_res.rds')
bc_res<-readRDS('combined/tf/all/bc_res.rds')
plotdata<-data.frame(matrix(0,ncol = 5,nrow = 160))
colnames(plotdata)<-c('Pathway','TF','Cancer','N_gene','Pvalue')
plotdata$Pathway<-rep(c('hippo signaling','canonical Wnt signaling pathway',
                        'transforming growth factor beta receptor signaling pathway',
                        'epithelial cell development','cell-cell junction organization'),32)
plotdata$TF<-rep(rep(c(1,2,3,4),each=5),8)
plotdata$Cancer<-rep(c(1,2,3,4,5,6,7,8),each=20)
for(i in 1:160){
  tmp_cancer=list(bc_res,bcc_res,cc_res,ec_res,lc_res,oc_res,plc_res,rcc_res)[[plotdata[i,3]]]
  tmp_tf=tmp_cancer[[plotdata[i,2]]]
  tmp_res<-subset(tmp_tf,Description==plotdata[i,1])
  if(nrow(tmp_res)==1){
    plotdata[i,4]<-tmp_res[1,9]
    plotdata[i,5]<-tmp_res[1,5]
  }
}
plotdata$TF<-rep(rep(c('TEAD1','TEAD2','TEAD3','TEAD4'),each=5),8)
plotdata$Cancer<-rep(c('BC','BCC','CC','EC','LC','OC','PLC','RCC'),each=20)
plotdata$Pathway<-rep(c('Hippo signaling','Wnt signaling',
                        'TGF-β signaling',
                        'epithelial cell development','cell-cell junction organization'),32)
plotdata$N_gene<-ifelse(plotdata$N_gene==0,NA,plotdata$N_gene)
p<-ggplot(plotdata, aes(x = TF, y = Pathway, color = Pvalue)) + 
  geom_point(aes(size = N_gene)) + 
  scale_color_distiller(palette = "RdBu") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line=element_line(size = .3, colour="black")) + 
  facet_grid(~ Cancer, scales = "free_x",space = 'free')
ggsave('f5/tead_pathway.pdf',p,width = 10,height = 3)

###################
###tead footprinting
###################
combined<-readRDS('combined/tf/combined.rds')
meta<-combined@meta.data
meta<-subset(meta,!sample %in% paste0('brca',1:16))
combined<-combined[,rownames(meta)]
DefaultAssay(combined)<-'peaks'
Idents(combined)<-'celltype'
for(i in c('TEAD1','TEAD2','TEAD3','TEAD4')){
  p<- PlotFootprint(combined, features = i)
  ggsave(paste0('f5/',gsub('::','',i),'.pdf'),p,width = 8,height = 6)
}
