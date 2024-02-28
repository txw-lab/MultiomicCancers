
setwd('data/pancancer_atac/tumor')
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
#color
allcolour_cc=c("#FEFB98","#C4AFCB","#96BF9E","#FDCD8A","#FEE491","#DBB6AF",'#D63048')
allcolour_ec=c("#FEFB98","#C4AFCB","#96BF9E",'#A75D2B',"#FDCD8A","#7FC97F","#DBB6AF",'#D63048')
allcolour_oc=c("#C4AFCB","#96BF9E","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_rcc=c("#FEFB98","#C4AFCB","#A75D2B","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
#plot
for(umap in c(umap_cc.atac,umap_ec.atac,umap_oc.atac,umap_rcc.atac)){
  p<-ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
    geom_point(size = 1 , alpha =1 )  +  
    scale_color_manual(values = allcolour)+
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标题
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background=element_rect(fill="white"))+
    theme(
      legend.title = element_blank(), #去掉legend.title 
      legend.key=element_rect(fill='white'), #
      legend.text = element_text(size=20), #设置legend标签的大小
      legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小
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
#color
allcolour_cc=c("#FEFB98","#C4AFCB","#96BF9E","#FDCD8A","#DBB6AF",'#D63048')
allcolour_ec=c("#FEFB98","#C4AFCB","#96BF9E",'#A75D2B',"#FDCD8A","#7FC97F","#DBB6AF",'#D63048')
allcolour_oc=c("#C4AFCB","#96BF9E","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
allcolour_rcc=c("#FEFB98","#C4AFCB","#A75D2B","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
#plot
for(umap in c(umap_cc.atac,umap_ec.atac,umap_oc.atac,umap_rcc.atac)){
  p<-ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
    geom_point(size = 1 , alpha =1 )  +  
    scale_color_manual(values = allcolour)+
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标题
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background=element_rect(fill="white"))+
    theme(
      legend.title = element_blank(), #去掉legend.title 
      legend.key=element_rect(fill='white'), #
      legend.text = element_text(size=20), #设置legend标签的大小
      legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小
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
color<-c("#FEFB98","#C4AFCB",'#96BF9E',"#A75D2B","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048')
combined_atac<-readRDS('combined/scatac/combined2.rds')
meta_atac<-combined_atac@meta.data
meta_atac$tissue<-'EC'
meta_atac$tissue<-ifelse(meta_atac$sample %in% c('CC1','CC2','CC3','CC4'),'CC',meta_atac$tissue)
meta_atac$tissue<-ifelse(meta_atac$sample %in% c('OC1','OC2'),'OC',meta_atac$tissue)
meta_atac$tissue<-ifelse(meta_atac$sample %in% c('RCC1','RCC2','RCC3'),'RCC',meta_atac$tissue)
combined<-readRDS('combined_rna.rds')
meta_rna<-combined@meta.data
plotdf<-data.frame(matrix(0,nrow = 9*4*2,ncol = 4))
colnames(plotdf)<-c('Seq','Tissue','Celltype','Percent')
plotdf$Seq<-c(rep('ATAC',36),rep('RNA',36))
plotdf$Tissue<-c(rep('CC',9),rep('EC',9),rep('OC',9),rep('RCC',9))
plotdf$Celltype<-rep(unique(meta_atac$celltype),4)
for(i in 1:72){
  if(plotdf$Seq[i]=='ATAC'){
    tmp<-subset(meta_atac,tissue==plotdf$Tissue[i])
    if(is.na(table(tmp$celltype)[plotdf$Celltype[i]])){
      next
    }else{
      plotdf$Percent[i]<-table(tmp$celltype)[plotdf$Celltype[i]]/nrow(tmp)
    }
  }else{
    tmp<-subset(meta_rna,orig.ident==plotdf$Tissue[i])
    if(is.na(table(tmp$celltype)[plotdf$Celltype[i]])){
      next
    }else{
      plotdf$Percent[i]<-table(tmp$celltype)[plotdf$Celltype[i]]/nrow(tmp)
    }
  }
  
}
plotdf$x<-ifelse(plotdf$Tissue=='CC',1,ifelse(plotdf$Tissue=='EC',2,ifelse(plotdf$Tissue=='OC',3,4)))
df1<-plotdf[1:36,]
df2<-plotdf[37:72,]
#plot
p<-ggplot()+
  geom_bar(data=df1,aes(x=x,y=Percent,fill=Celltype),stat = "identity",position = "stack",width=0.3)+
  geom_bar(data=df2,aes(x=x+0.3+0.05,y=Percent,fill=Celltype),stat="identity",position = "stack",width=0.3)+
  scale_x_continuous(breaks = c(1.15,2.15,3.15,4.15),
                     labels = c("CC","EC","OC",'RCC'))+
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
###peak annotation
###################
#load data
crc_peak<-readRDS('cc/chipseeker/peakanno.rds')
crc_peak<-data.frame(type=unique(crc_peak$type),counts=as.vector(table(crc_peak$type)[unique(crc_peak$type)]),tumor=rep('CC',7))
oc_peak<-readRDS('oc/chipseeker/peakanno.rds')
oc_peak<-data.frame(type=unique(oc_peak$type),counts=as.vector(table(oc_peak$type)[unique(oc_peak$type)]),tumor=rep('OC',7))
ec_peak<-readRDS('ec/chipseeker/peakanno.rds')
ec_peak<-data.frame(type=unique(ec_peak$type),counts=as.vector(table(ec_peak$type)[unique(ec_peak$type)]),tumor=rep('EC',7))
rcc_peak<-readRDS('rcc/chipseeker/peakanno.rds')
rcc_peak<-data.frame(type=unique(rcc_peak$type),counts=as.vector(table(rcc_peak$type)[unique(rcc_peak$type)]),tumor=rep('RCC',7))
combined_peak<-readRDS('chipseeker/peakanno.rds')
combined_peak<-data.frame(type=unique(combined_peak$type),counts=as.vector(table(combined_peak$type)[unique(combined_peak$type)]),tumor=rep('All_sample',7))
#plot
plotdf<-rbind(crc_peak,oc_peak,ec_peak,rcc_peak,combined_peak)
p<-ggplot(plotdf, aes
          (tumor, weight = counts, fill = type)) +
  #geom_hline(yintercept = seq(25, 100, 25), color = 'gray') +
  geom_bar(color = "white", width = .7, position = 'stack') +
  labs( y = 'Counts') +
  scale_fill_brewer(palette = "Set3")+
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
ggsave('f1/Chipseeker_anno.pdf',p,height = 6,width = 8)

###########################################
###cancer-associated celltype-specific DARs
###########################################
#load data
peaks_cc<-readRDS('../cc/data/scatac/peaks_crc.rds')
peaks_ec<-readRDS('../ec/data/scatac/peaks_ec.rds')
peaks_oc<-readRDS('../oc/data/scatac/peaks_oc.rds')
peaks_rcc<-readRDS('../rcc/data/scatac/peak_rcc.rds')
#celltype DARs
cc_list<-list()
for(i in unique(peaks_cc$cluster)){
  tmp<-subset(peaks_cc,cluster==i)
  chr<-c()
  start<-c()
  end<-c()
  for(j in tmp$gene){
    t<-strsplit(j,'-')[[1]]
    chr<-c(chr,t[1])
    start<-c(start,as.numeric(t[2]))
    end<-c(end,as.numeric(t[3]))
  }
  cc_list[[i]]<-GRanges(seqnames = Rle(chr, rep(1,length(chr))),
                        ranges = IRanges(start, end = end))
}
ec_list<-list()
for(i in unique(peaks_ec$cluster)){
  tmp<-subset(peaks_ec,cluster==i)
  chr<-c()
  start<-c()
  end<-c()
  for(j in tmp$gene){
    t<-strsplit(j,'-')[[1]]
    chr<-c(chr,t[1])
    start<-c(start,as.numeric(t[2]))
    end<-c(end,as.numeric(t[3]))
  }
  ec_list[[i]]<-GRanges(seqnames = Rle(chr, rep(1,length(chr))),
                        ranges = IRanges(start, end = end))
}
oc_list<-list()
for(i in unique(peaks_oc$cluster)){
  tmp<-subset(peaks_oc,cluster==i)
  chr<-c()
  start<-c()
  end<-c()
  for(j in tmp$gene){
    t<-strsplit(j,'-')[[1]]
    chr<-c(chr,t[1])
    start<-c(start,as.numeric(t[2]))
    end<-c(end,as.numeric(t[3]))
  }
  oc_list[[i]]<-GRanges(seqnames = Rle(chr, rep(1,length(chr))),
                        ranges = IRanges(start, end = end))
}
rcc_list<-list()
for(i in unique(peaks_rcc$cluster)){
  tmp<-subset(peaks_rcc,cluster==i)
  chr<-c()
  start<-c()
  end<-c()
  for(j in tmp$gene){
    t<-strsplit(j,'-')[[1]]
    chr<-c(chr,t[1])
    start<-c(start,as.numeric(t[2]))
    end<-c(end,as.numeric(t[3]))
  }
  rcc_list[[i]]<-GRanges(seqnames = Rle(chr, rep(1,length(chr))),
                         ranges = IRanges(start, end = end))
}
#cancer-associated
df_cc<-data.frame(celltype=names(cc_list),Count=rep(0,7),specific=rep(0,7),row.names = names(cc_list),percent=rep(0,7))
for(i in df_cc$celltype){
  tmp<-GRanges()
  if(i %in% names(ec_list)){
    tmp<-c(tmp,ec_list[[i]])
  }
  if(i %in% names(oc_list)){
    tmp<-c(tmp,oc_list[[i]])
  }
  if(i %in% names(rcc_list)){
    tmp<-c(tmp,rcc_list[[i]])
  }
  over<-findOverlaps(cc_list[[i]],tmp)
  df_cc[i,2]<-nrow(subset(peaks_cc,cluster==i))
  df_cc[i,3]<-df_cc[i,2]-length(unique(over@from))
  df_cc[i,4]<-df_cc[i,3]/df_cc[i,2]
}
df_ec<-data.frame(celltype=names(ec_list),Count=rep(0,8),specific=rep(0,8),row.names = names(ec_list),percent=rep(0,8))
for(i in df_ec$celltype){
  tmp<-GRanges()
  if(i %in% names(cc_list)){
    tmp<-c(tmp,cc_list[[i]])
  }
  if(i %in% names(oc_list)){
    tmp<-c(tmp,oc_list[[i]])
  }
  if(i %in% names(rcc_list)){
    tmp<-c(tmp,rcc_list[[i]])
  }
  over<-findOverlaps(ec_list[[i]],tmp)
  df_ec[i,2]<-nrow(subset(peaks_ec,cluster==i))
  df_ec[i,3]<-df_ec[i,2]-length(unique(over@from))
  df_ec[i,4]<-df_ec[i,3]/df_ec[i,2]
}
df_oc<-data.frame(celltype=names(oc_list),Count=rep(0,7),spocific=rep(0,7),row.names = names(oc_list),percent=rep(0,7))
for(i in df_oc$celltype){
  tmp<-GRanges()
  if(i %in% names(cc_list)){
    tmp<-c(tmp,cc_list[[i]])
  }
  if(i %in% names(ec_list)){
    tmp<-c(tmp,ec_list[[i]])
  }
  if(i %in% names(rcc_list)){
    tmp<-c(tmp,rcc_list[[i]])
  }
  over<-findOverlaps(oc_list[[i]],tmp)
  df_oc[i,2]<-nrow(subset(peaks_oc,cluster==i))
  df_oc[i,3]<-df_oc[i,2]-length(unique(over@from))
  df_oc[i,4]<-df_oc[i,3]/df_oc[i,2]
}
df_rcc<-data.frame(celltype=names(rcc_list),Count=rep(0,8),sprccific=rep(0,8),row.names = names(rcc_list),percent=rep(0,8))
for(i in df_rcc$celltype){
  tmp<-GRanges()
  if(i %in% names(cc_list)){
    tmp<-c(tmp,cc_list[[i]])
  }
  if(i %in% names(ec_list)){
    tmp<-c(tmp,ec_list[[i]])
  }
  if(i %in% names(oc_list)){
    tmp<-c(tmp,oc_list[[i]])
  }
  over<-findOverlaps(rcc_list[[i]],tmp)
  df_rcc[i,2]<-nrow(subset(peaks_rcc,cluster==i))
  df_rcc[i,3]<-df_rcc[i,2]-length(unique(over@from))
  df_rcc[i,4]<-df_rcc[i,3]/df_rcc[i,2]
}
#percentage of DARs
tumor<-mean(c(df_cc['Tumor cell',4],df_ec['Tumor cell',4],df_oc['Tumor cell',4],df_rcc['Tumor cell',4]))
t<-mean(c(df_cc['T cell',4],df_ec['T cell',4],df_oc['T cell',4],df_rcc['T cell',4]))
fibro<-mean(c(df_cc['Fibroblast',4],df_ec['Fibroblast',4],df_oc['Fibroblast',4]))
myo<-mean(c(df_rcc['Myofibroblast',4],df_ec['Myofibroblast',4],df_oc['Myofibroblast',4]))
b<-mean(c(df_cc['B cell',4],df_ec['B cell',4],df_rcc['B cell',4]))
endo<-mean(c(df_cc['Endothelial',4],df_ec['Endothelial',4],df_oc['Endothelial',4],df_rcc['Endothelial',4]))
pla<-mean(c(df_cc['Plasma',4],df_oc['Plasma',4],df_rcc['Plasma',4]))
mye<-mean(c(df_cc['Myeloid cell',4],df_ec['Myeloid cell',4],df_oc['Myeloid cell',4],df_rcc['Myeloid cell',4]))
mast<-mean(c(df_ec['Mast',4],df_rcc['Mast',4]))
#count of DARs
tumor_count<-mean(c(df_cc['Tumor cell',3],df_ec['Tumor cell',3],df_oc['Tumor cell',3],df_rcc['Tumor cell',3]))
t_count<-mean(c(df_cc['T cell',3],df_ec['T cell',3],df_oc['T cell',3],df_rcc['T cell',3]))
fibro_count<-mean(c(df_cc['Fibroblast',3],df_ec['Fibroblast',3],df_oc['Fibroblast',3]))
myo_count<-mean(c(df_rcc['Myofibroblast',3],df_ec['Myofibroblast',3],df_oc['Myofibroblast',3]))
b_count<-mean(c(df_cc['B cell',3],df_ec['B cell',3],df_rcc['B cell',3]))
endo_count<-mean(c(df_cc['Endothelial',3],df_ec['Endothelial',3],df_oc['Endothelial',3],df_rcc['Endothelial',3]))
pla_count<-mean(c(df_cc['Plasma',3],df_oc['Plasma',3],df_rcc['Plasma',3]))
mye_count<-mean(c(df_cc['Myeloid cell',3],df_ec['Myeloid cell',3],df_oc['Myeloid cell',3],df_rcc['Myeloid cell',3]))
mast_count<-mean(c(df_ec['Mast',3],df_rcc['Mast',3]))
#plot
library(ggplot2)
plotdata<-data.frame(celltype=c('Mast','Fibroblast','Tumor cell','Myofibroblast','Plasma','B cell','Endothelial','Myeloid cell','T cell'),
                     Percentage=c(mast,fibro,tumor,myo,pla,b,endo,mye,t),
                     Number=c(mast_count,fibro_count,tumor_count,myo_count,pla_count,b_count,endo_count,mye_count,t_count))
plotdata$celltype<-factor(plotdata$celltype,levels = plotdata$celltype)
p<-ggplot(plotdata)+
  geom_bar(aes(x=celltype,y=Number),stat = 'identity',fill='#F8766D')+
  geom_point(aes(x=celltype,y=Percentage*2500),size=3,shape=17)+
  geom_line(aes(x=celltype,y=Percentage*2500,group=1))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,size = 11))+
  scale_y_continuous(sec.axis=sec_axis(~.*0.0004, name="Percentage", labels=scales::percent))
ggsave('f1/CancerCelltypeDARs.pdf',p,height = 6,width = 8)