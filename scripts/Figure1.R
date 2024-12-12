
setwd('tumor')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)

###################
###umap atac
###################
#load data
cc.atac<-readRDS('cc/data/scatac/cc.rds')
Idents(cc.atac)<-'celltype'
ec.atac<-readRDS('ec/data/scatac/ec.rds')
Idents(ec.atac)<-'celltype'
oc.atac<-readRDS('oc/data/scatac/oc.rds')
Idents(oc.atac)<-'celltype'
rcc.atac<-readRDS('rcc/data/scatac/rcc.rds')
Idents(rcc.atac)<-'celltype'
bcc<-readRDS('bcc/data/scrna/bcc.rds')
Idents(bcc)<-'celltype'
plc<-readRDS('plc/data/scrna/plc.rds')
Idents(plc)<-'celltype'
lc.atac<-readRDS('lc/data/scatac/lc.rds')
Idents(lc.atac)<-'celltype'
bc.atac<-readRDS('bc/data/scatac/bc.rds')
Idents(bc.atac)<-'celltype'
#get embedding
umap_cc.atac = cc.atac@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = cc.atac$celltype) 
umap_ec.atac = ec.atac@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = ec.atac$celltype) 
umap_oc.atac = oc.atac@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = oc.atac$celltype) 
umap_rcc.atac = rcc.atac@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = rcc.atac$celltype) 
umap_bcc = bcc.atac@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = bcc.atac$celltype) 
umap_plc.atac = plc.atac@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = plc.atac$celltype) 
umap_lc.atac = lc.atac@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = lc.atac$celltype) 
umap_bc.atac = bc.atac@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = bc.atac$celltype) 
#color
allcolour_cc=c("#FEFB98","#C4AFCB","#96BF9E","#FDCD8A","#FEE491","#DBB6AF",'#D63048')
allcolour_ec=c("#FEFB98","#C4AFCB","#96BF9E",'#A75D2B',"#FDCD8A","#7FC97F","#DBB6AF",'#D63048')
allcolour_oc=c("#C4AFCB","#96BF9E","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_rcc=c("#FEFB98","#C4AFCB","#A75D2B","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_bcc=c("#FEFB98","#C4AFCB","#96BF9E","#A75D2B","#769AA8","#FDCD8A","#7FC97F","#C91889",'#FEE491',"#DBB6AF",'#D63048')
allcolour_plc=c("#FEFB98","#C4AFCB","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_lc=c("#FEFB98","#C4AFCB",'#96BF9E',"#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_bc=c("#C4AFCB",'#96BF9E','#F3BD92',"#FDCD8A","#7FC97F",'#D63048')
#plot
for(umap in c(umap_cc.atac,umap_ec.atac,umap_oc.atac,umap_rcc.atac,umap_plc.atac,umap_bcc.atac,umap_lc.atac,umap_bc.atac)){
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
  ggsave('f1/Umap.atac.pdf',p,height = 6,width = 8)
}

###################
###umap rna
###################
#load data
cc.rna<-readRDS('cc/data/scrna/cc.rds')
Idents(cc.rna)<-'celltype'
ec.rna<-readRDS('ec/data/scrna/ec.rds')
Idents(ec.rna)<-'celltype'
oc.rna<-readRDS('oc/data/scrna/oc.rds')
Idents(oc.rna)<-'celltype'
rcc.rna<-readRDS('rcc/data/scrna/rcc.rds')
Idents(rcc.rna)<-'celltype'
lc.rna<-readRDS('luad/data/scrna/lc.rds')
Idents(lc.rna)<-'celltype'
bc.rna<-readRDS('brca/data/scrna/bc.rds')
Idents(bc.rna)<-'celltype'
#get embedding
umap_cc.rna = cc.rna@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = cc.rna$celltype) 
umap_ec.rna = ec.rna@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = ec.rna$celltype) 
umap_oc.rna = oc.rna@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = oc.rna$celltype) 
umap_rcc.rna = rcc.rna@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = rcc.rna$celltype) 
umap_lc.rna = lc.rna@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = lc.rna$celltype) 
umap_bc.rna = bc.rna@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = bc.rna$celltype) 
#color
allcolour_cc=c("#FEFB98","#C4AFCB","#96BF9E","#FDCD8A","#DBB6AF",'#D63048')
allcolour_ec=c("#FEFB98","#C4AFCB","#96BF9E",'#A75D2B',"#FDCD8A","#7FC97F","#DBB6AF",'#D63048')
allcolour_oc=c("#C4AFCB","#96BF9E","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_rcc=c("#FEFB98","#C4AFCB","#A75D2B","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_lc=c("#FEFB98","#C4AFCB",'#96BF9E','#A75D2B',"#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_bc=c("#FEFB98","#C4AFCB",'#96BF9E',"#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
#plot
for(umap in c(umap_cc.atac,umap_ec.atac,umap_oc.atac,umap_rcc.atac,umap_lc.atac,umap_bc.atac)){
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
  ggsave('f1/Umap.rna.pdf',p,height = 6,width = 8)
}

###################
###cell percentage
###################
color<-c("#FEFB98","#C4AFCB",'#96BF9E','#F3BD92',"#A75D2B",'#769AA8',"#FDCD8A","#7FC97F","#C91889",'#FEE491',"#DBB6AF",'#D63048')
combined_atac<-readRDS('combined/scatac/combined.rds')
meta_atac<-combined_atac@meta.data
combined_rna<-readRDS('combined/scrna/combined_rna.rds')
meta_rna<-combined_rna@meta.data
plotdf<-data.frame(matrix(0,nrow = 12*8*2,ncol = 4))
colnames(plotdf)<-c('Seq','Tissue','Celltype','Percent')
plotdf$Seq<-c(rep('ATAC',96),rep('RNA',96))
plotdf$Tissue<-c(rep('BC',12),rep('BCC',12),rep('CC',12),rep('EC',12),rep('LC',12),rep('OC',12),rep('PLC',12),rep('RCC',12))
plotdf$Celltype<-rep(unique(meta_atac$celltype),8)
for(i in 1:192){
  if(plotdf$Seq[i]=='ATAC'){
    tmp<-subset(meta_atac,tissue==plotdf$Tissue[i])
    if(is.na(table(tmp$celltype)[plotdf$Celltype[i]])){
      next
    }else{
      plotdf$Percent[i]<-table(tmp$celltype)[plotdf$Celltype[i]]/nrow(tmp)
    }
  }else{
    tmp<-subset(meta_rna,tissue==plotdf$Tissue[i])
    if(is.na(table(tmp$celltype)[plotdf$Celltype[i]])){
      next
    }else{
      plotdf$Percent[i]<-table(tmp$celltype)[plotdf$Celltype[i]]/nrow(tmp)
    }
  }
}
plotdf$x<-ifelse(plotdf$Tissue=='BCC',2,ifelse(plotdf$Tissue=='CC',3,
                 ifelse(plotdf$Tissue=='EC',4,ifelse(plotdf$Tissue=='LC',5,
                 ifelse(plotdf$Tissue=='OC',6,ifelse(plotdf$Tissue=='PLC',7,
                 ifelse(plotdf$Tissue=='RCC',8,1)))))))
df1<-plotdf[1:96,]
df2<-plotdf[97:192,]
#plot
p<-ggplot()+
  geom_bar(data=df1,aes(x=x,y=Percent,fill=Celltype),stat = "identity",position = "stack",width=0.3)+
  geom_bar(data=df2,aes(x=x+0.3+0.05,y=Percent,fill=Celltype),stat="identity",position = "stack",width=0.3)+
  scale_x_continuous(breaks = c(1.15,2.15,3.15,4.15,5.15,6.15,7.15,8.15),
                     labels = c("BC","BCC","CC","EC",'LC',"OC",'PLC','RCC'))+
  scale_fill_manual(values = color)+
  theme_set(theme_bw())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black", size = 1, ))+
  xlab('')+
  guides(fill = guide_legend(title = NULL))+
  theme(axis.title.y = element_text(face = "bold", colour = 'black', size = 17))+
  theme(axis.text.x = element_text(size = 16, color = "black", face = "bold"))+
  theme(legend.text=element_text( size = 14))
ggsave('f1/celltype_percentage.pdf',p,height = 6,width = 8)

###################
###quality control
###################
combined_atac<-readRDS('combined/scatac/combined.rds')
Idents(combined_atac)<-'orig.ident'
p1<-FragmentHistogram(combined)
ggsave('f1/frag_plot.pdf',p1,height = 5,width = 5)
combined_atac<-TSSEnrichment(combined_atac,fast = FALSE)
p2<-TSSPlot(combined_atac)
ggsave('f1/tssplot.pdf',p2,height = 5,width = 5)

#####################
###umap_combined atac
#####################
combined_atac<-readRDS('combined/scatac/combined.rds')
umap_atac = combined_atac@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = combined_atac$celltype) 
color_atac<-c("#FEFB98","#C4AFCB",'#96BF9E','#F3BD92',"#A75D2B",'#769AA8',"#FDCD8A","#7FC97F","#C91889",'#FEE491',"#DBB6AF",'#D63048')

p_atac<-ggplot(umap_atac,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  
  scale_color_manual(values = color_atac)+
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
  geom_segment(aes(x = min(umap_atac$UMAP_1) , y = min(umap_atac$UMAP_2) ,
                   xend = min(umap_atac$UMAP_1) +3, yend = min(umap_atac$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap_atac$UMAP_1)  , y = min(umap_atac$UMAP_2)  ,
                   xend = min(umap_atac$UMAP_1) , yend = min(umap_atac$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap_atac$UMAP_1) +1.5, y = min(umap_atac$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap_atac$UMAP_1) -1, y = min(umap_atac$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) +
  theme(legend.position = "none")
ggsave('f1/umap_atac.pdf',p_atac,height =6 ,width = 6)

#####################
###umap_combined rna
#####################
combined_rna<-readRDS('combined/scrna/combined.rds')
color_rna<-c("#FEFB98","#C4AFCB",'#96BF9E',"#A75D2B",'#769AA8',"#FDCD8A","#7FC97F","#C91889",'#FEE491',"#DBB6AF",'#D63048')
umap_rna = combined_rna@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = combined_rna$celltype) 
p_rna<-ggplot(umap_rna,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 0.3 , alpha =1 )  +  
  scale_color_manual(values = color)+
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
  geom_segment(aes(x = min(umap_rna$UMAP_1) , y = min(umap_rna$UMAP_2) ,
                   xend = min(umap_rna$UMAP_1) +3, yend = min(umap_rna$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap_rna$UMAP_1)  , y = min(umap_rna$UMAP_2)  ,
                   xend = min(umap_rna$UMAP_1) , yend = min(umap_rna$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap_rna$UMAP_1) +1.5, y = min(umap_rna$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap_rna$UMAP_1) -1, y = min(umap_rna$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) +
  theme(legend.position = "none")
ggsave('f1/umap_rna.pdf',p_rna,height = 6,width = 6)
 
