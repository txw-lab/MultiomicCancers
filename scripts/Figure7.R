
setwd('tumor')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)
library(RColorBrewer)
library(ChIPseeker)
library(IRanges)
library(GenomicRanges)

###################
###patient plots
###################
#patient1
combined1<-readRDS('cc/multiome/rds/combined.rds')
combined1$celltype<-ifelse(combined1$celltype=='Macrophage','Myeloid cell',combined1$celltype)
umap_p1 = combined1@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = combined1$celltype) 
allcolour1=c("#FEFB98","#C4AFCB",'#C91889',"#96BF9E",'#ADB5BD',"#FDCD8A",'#769AA8','#866148',"#DBB6AF",'#D63048')
p1<-ggplot(umap_p1,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  
  theme_test()+
  scale_color_manual(values = allcolour1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 2, fill = NA), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+theme(legend.position = "none")+
  labs(title = 'Patient1-scRNA')

umap_p2 = combined1@reductions$umap_atac@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = combined1$celltype)
p2<-ggplot(umap_p2,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  
  theme_test()+
  scale_color_manual(values = allcolour1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 2, fill = NA), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+theme(legend.position = "none")+
  labs(title = 'Patient1-scATAC')
grid.arrange(p1, p2,nrow=1)

#patient2
combined2<-readRDS('cc/multiome1/rds/combined.rds')
umap_p3 = combined2@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = combined2$subtype) 
allcolour2=c("#FEFB98","#C4AFCB",'#C91889',"#96BF9E",'#ADB5BD',"#FDCD8A",'#7FC97F','#866148',"#DBB6AF",'#D63048')
p3<-ggplot(umap_p3,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  
  theme_test()+
  scale_color_manual(values = allcolour2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 2, fill = NA), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+theme(legend.position = "none")+
  labs(title = 'Patient2-scRNA')

umap_p4 = combined2@reductions$umap_atac@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = combined2$subtype)
p4<-ggplot(umap_p4,aes(x= umap_atac_1 , y = umap_atac_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  
  theme_test()+
  scale_color_manual(values = allcolour1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 2, fill = NA), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+theme(legend.position = "none")+
  labs(title = 'Patient2-scATAC')

#patient3
atac3<-readRDS('cc/p3/data/scatac/p3.rds')
umap_p5 = atac3@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = atac3$celltype) 
allcolour3=c("#FEFB98",'#C91889',"#96BF9E","#FDCD8A",'#866148',"#DBB6AF",'#D63048')
p5<-ggplot(umap_p5,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  
  theme_test()+
  scale_color_manual(values = allcolour3)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 2, fill = NA), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+theme(legend.position = "none")+
  labs(title = 'Patient3-scATAC')

rna3<-readRDS('../cc/data/normal/p3/data/scrna/p3.rds')
umap_p6 = rna3@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = rna3$celltype)
allcolour4=c("#FEFB98",'#C91889',"#96BF9E","#FDCD8A","#DBB6AF",'#D63048')
p6<-ggplot(umap_p6,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  
  theme_test()+
  scale_color_manual(values = allcolour4)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 2, fill = NA), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+theme(legend.position = "none")+
  labs(title = 'Patient3-scRNA')

p<-grid.arrange(p1,p2,p3, p4,p6,p5,nrow=1)
ggsave('f7/umap_plots.pdf',p,height = 3,width = 16)

###################
###dot t value
###################
tumor1_rna<-readRDS('cc/multiome/rds/combined.rds')
tumor1_rna<-subset(tumor1_rna,celltype %in% c('Tumor cell','Epithelial'))
rna1<-data.frame(t(tumor1_rna@assays$RNA@data[c('CEBPG','LEF1','SOX4','TCF7','TEAD4'),]))
rna1[,'type']<-tumor1_rna$celltype

tumor2_rna<-readRDS('cc/multiome1/rds/combined.rds')
tumor2_rna<-subset(tumor2_rna,subtype %in% c('Tumor cell','Epithelial'))
rna2<-data.frame(t(tumor2_rna@assays$RNA@data[c('CEBPG','LEF1','SOX4','TCF7','TEAD4'),]))
rna2[,'type']<-tumor2_rna$subtype

tumor3_rna<-readRDS('cc/data/p3/data/scrna/p3.rds')
tumor3_rna<-subset(tumor3_rna,celltype %in% c('Tumor cell','Epithelial'))
rna3<-data.frame(t(tumor3_rna@assays$RNA@data[c('CEBPG','LEF1','SOX4','TCF7','TEAD4'),]))
rna3[,'type']<-tumor3_rna$celltype

tumor1_chrom<-readRDS('cc/multiome/chromvar/cc.rds')
tumor1_chrom<-subset(tumor1_chrom,celltype %in% c('Tumor cell','Epithelial'))
chrom1<-data.frame(t(tumor1_chrom@assays$chromvar@data[c('MA0838.1','MA0768.1','MA0867.2',"MA0769.2",'MA0809.2'),]))
colnames(chrom1)<-c('CEBPG','LEF1','SOX4','TCF7','TEAD4')
chrom1[,'type']<-tumor1_chrom$celltype

tumor2_chrom<-readRDS('cc/multiome1/chromvar/cc.rds')
tumor2_chrom<-subset(tumor2_chrom,subtype %in% c('Tumor cell','Epithelial'))
chrom2<-data.frame(t(tumor2_chrom@assays$chromvar@data[c('MA0838.1','MA0768.1','MA0867.2',"MA0769.2",'MA0809.2'),]))
colnames(chrom2)<-c('CEBPG','LEF1','SOX4','TCF7','TEAD4')
chrom2[,'type']<-tumor2_chrom$subtype

tumor3_chrom<-readRDS('cc/data/normal/p3/chromvar/cc.rds')
tumor3_chrom<-subset(tumor3_chrom,celltype %in% c('Tumor cell','Epithelial'))
chrom3<-data.frame(t(tumor3_chrom@assays$chromvar@data[c('MA0838.1','MA0768.1','MA0867.2',"MA0769.2",'MA0809.2'),]))
colnames(chrom3)<-c('CEBPG','LEF1','SOX4','TCF7','TEAD4')
chrom3[,'type']<-tumor3_chrom$celltype

df<-data.frame(matrix(0,nrow = 30,ncol = 7))
colnames(df)<-c('Patient','TF','Type','Mean_tumor','Mean_epi','T_value','P_value')
df$Patient<-rep(c('P1','P2','p3'),each=10)
df$TF<-rep(c('CEBPG','LEF1','SOX4','TCF7','TEAD4'),6)
df$Type<-rep(rep(c('Expression','Chromvar_score'),each=5),3)
for(i in 1:5){
  t<-t.test(subset(rna1,type=='Tumor cell')[,df[i,2]],subset(rna1,type=='Epithelial')[,df[i,2]])
  df[i,4]<-t$estimate[1]
  df[i,5]<-t$estimate[2]
  df[i,6]<-t$statistic
  df[i,7]<-t$p.value
}
for(i in 6:10){
  t<-t.test(subset(chrom1,type=='Tumor cell')[,df[i,2]],subset(chrom1,type=='Epithelial')[,df[i,2]])
  df[i,4]<-t$estimate[1]
  df[i,5]<-t$estimate[2]
  df[i,6]<-t$statistic
  df[i,7]<-t$p.value
}
for(i in 11:15){
  t<-t.test(subset(rna2,type=='Tumor cell')[,df[i,2]],subset(rna2,type=='Epithelial')[,df[i,2]])
  df[i,4]<-t$estimate[1]
  df[i,5]<-t$estimate[2]
  df[i,6]<-t$statistic
  df[i,7]<-t$p.value
}
for(i in 16:20){
  t<-t.test(subset(chrom2,type=='Tumor cell')[,df[i,2]],subset(chrom2,type=='Epithelial')[,df[i,2]])
  df[i,4]<-t$estimate[1]
  df[i,5]<-t$estimate[2]
  df[i,6]<-t$statistic
  df[i,7]<-t$p.value
}
for(i in 21:25){
  t<-t.test(subset(rna3,type=='Tumor cell')[,df[i,2]],subset(rna3,type=='Epithelial')[,df[i,2]])
  df[i,4]<-t$estimate[1]
  df[i,5]<-t$estimate[2]
  df[i,6]<-t$statistic
  df[i,7]<-t$p.value
}
for(i in 26:30){
  t<-t.test(subset(chrom3,type=='Tumor cell')[,df[i,2]],subset(chrom3,type=='Epithelial')[,df[i,2]])
  df[i,4]<-t$estimate[1]
  df[i,5]<-t$estimate[2]
  df[i,6]<-t$statistic
  df[i,7]<-t$p.value
}
df$logFC<-log2(df$Mean_tumor/df$Mean_epi)

p<-ggplot(df, aes(x = Type, y = TF, color = P_value)) + 
  geom_point(aes(size = T_value)) + 
  scale_color_distiller(palette = "RdBu") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line=element_line(size = .3, colour="black"))+
  facet_grid(~ Patient, scales = "free_x",space = 'free')
ggsave('f7/dot_t.pdf',p,height = 5,width = 3.5)

###################
###tfs footprinting
###################
#p1
tumor1<-readRDS('cc/multiome/tf/epi.rds')
p1<-PlotFootprint(tumor1, features = c("LEF1"))
p2<-PlotFootprint(tumor1, features = c("SOX4"))
p3<-PlotFootprint(tumor1, features = c("TCF7"))
p4<-PlotFootprint(tumor1, features = c("TEAD4"))
p5<-PlotFootprint(tumor1,features = c('CEBPG'))
ggsave('f7/lef1.pdf',p1,height = 4,width = 5)
ggsave('f7/sox4.pdf',p2,height = 4,width = 5)
ggsave('f7/tcf7.pdf',p3,height = 4,width = 5)
ggsave('f7/tead4.pdf',p4,height = 4,width = 5)
ggsave('f7/cebpg.pdf',p5,height = 4,width = 5)

#p2
tumor2<-readRDS('cc/multiome1/tf/epi.rds')
p1<-PlotFootprint(tumor2, features = c("LEF1"))
p2<-PlotFootprint(tumor2, features = c("SOX4"))
p3<-PlotFootprint(tumor2, features = c("TCF7"))
p4<-PlotFootprint(tumor2, features = c("TEAD4"))
p5<-PlotFootprint(tumor2,features = c('CEBPG'))
ggsave('f7/lef1_2.pdf',p1,height = 4,width = 5)
ggsave('f7/sox4_2.pdf',p2,height = 4,width = 5)
ggsave('f7/tcf7_2.pdf',p3,height = 4,width = 5)
ggsave('f7/tead4_2.pdf',p4,height = 4,width = 5)
ggsave('f7/cebpg_2.pdf',p5,height = 4,width = 5)

#p3
tumor3<-readRDS('cc/data/p3/tf/epi.rds')
p1<-PlotFootprint(tumor3, features = c("LEF1"))
p2<-PlotFootprint(tumor3, features = c("SOX4"))
p3<-PlotFootprint(tumor3, features = c("TCF7"))
p4<-PlotFootprint(tumor3, features = c("TEAD4"))
p5<-PlotFootprint(tumor3,features = c('CEBPG'))
ggsave('f7/lef1_3.pdf',p1,height = 4,width = 5)
ggsave('f7/sox4_3.pdf',p2,height = 4,width = 5)
ggsave('f7/tcf7_3.pdf',p3,height = 4,width = 5)
ggsave('f7/tead4_3.pdf',p4,height = 4,width = 5)
ggsave('f7/cebpg_3.pdf',p5,height = 4,width = 5)

#normdata
tf_p1<-GetFootprintData(tumor1,c('CEBPG','LEF1','SOX4','TCF7','TEAD4'))
tf_p11<-subset(tf_p1,class=='Observed')
tf_p12<-subset(tf_p1,class=='Observed' & abs(position)<=100)
tf_p2<-GetFootprintData(tumor2,c('CEBPG','LEF1','SOX4','TCF7','TEAD4'))
tf_p21<-subset(tf_p2,class=='Observed')
tf_p22<-subset(tf_p2,class=='Observed' & abs(position)<=100)
tf_p3<-GetFootprintData(tumor3,c('CEBPG','LEF1','SOX4','TCF7','TEAD4'))
tf_p31<-subset(tf_p3,class=='Observed')
tf_p32<-subset(tf_p3,class=='Observed' & abs(position)<=100)
df<-data.frame(matrix(0,nrow = 15,ncol=4))
colnames(df)<-c('Patient','TF','P_value(250bp)','P_value(100bp)')
df$Patient<-rep(c('p1','p2','p3'),each=5)
df$TF<-rep(c('CEBPG','LEF1','SOX4','TCF7','TEAD4'),3)
for(i in 1:5){
  t<-t.test(subset(tf_p11,feature==df[i,2] & group=='Tumor cell')[,4],
            subset(tf_p11,feature==df[i,2] & group=='Epithelial')[,4])
  df[i,3]<-t$p.value
  t<-t.test(subset(tf_p12,feature==df[i,2] & group=='Tumor cell')[,4],
            subset(tf_p12,feature==df[i,2] & group=='Epithelial')[,4])
  df[i,4]<-t$p.value
}
for(i in 6:10){
  t<-t.test(subset(tf_p21,feature==df[i,2] & group=='Tumor cell')[,4],
            subset(tf_p21,feature==df[i,2] & group=='Epithelial')[,4])
  df[i,3]<-t$p.value
  t<-t.test(subset(tf_p22,feature==df[i,2] & group=='Tumor cell')[,4],
            subset(tf_p22,feature==df[i,2] & group=='Epithelial')[,4])
  df[i,4]<-t$p.value
}
for(i in 11:15){
  t<-t.test(subset(tf_p31,feature==df[i,2] & group=='Tumor cell')[,4],
            subset(tf_p31,feature==df[i,2] & group=='Epithelial')[,4])
  df[i,3]<-t$p.value
  t<-t.test(subset(tf_p32,feature==df[i,2] & group=='Tumor cell')[,4],
            subset(tf_p32,feature==df[i,2] & group=='Epithelial')[,4])
  df[i,4]<-t$p.value
}
write.csv(df,'tables/normlized_tfs.csv',quote = F,row.names = F)
