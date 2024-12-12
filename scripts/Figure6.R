
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
###umap rna
###################
colon<-readRDS('cc/data/scrna/cc.rds')
colon$celltype<-factor(colon$celltype,levels = sort(unique(colon$celltype)))
Idents(colon)<-'celltype'
allcolour=c("#FEFB98","#C4AFCB",'#C91889',"#96BF9E",'#A75D2B',"#FDCD8A",'#7FC97F','#866148','#FEE491',"#DBB6AF",'#D63048','#666666')
umap_cc = colon@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = colon$celltype) 
p1<-ggplot(umap_cc,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
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
    legend.text = element_text(size=10), 
    legend.key.size=unit(0.6,'cm') ) +  
  guides(color = guide_legend(override.aes = list(size=4.5)))+ 
  geom_segment(aes(x = min(umap_cc$UMAP_1)-2 , y = min(umap_cc$UMAP_2)-2 ,
                   xend = min(umap_cc$UMAP_1) +1, yend = min(umap_cc$UMAP_2)-2 ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+
  geom_segment(aes(x = min(umap_cc$UMAP_1)-2  , y = min(umap_cc$UMAP_2)-2  ,
                   xend = min(umap_cc$UMAP_1)-2 , yend = min(umap_cc$UMAP_2) + 1),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap_cc$UMAP_1) -0.5, y = min(umap_cc$UMAP_2) -3, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) +
  annotate("text", x = min(umap_cc$UMAP_1) -3, y = min(umap_cc$UMAP_2)-0.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)
ggsave('f6/rna_umap.pdf',p1,width = 8,height = 6)

###################
###umap atac
###################
colon<-readRDS('cc/data/scatac/cc.rds')
Idents(colon)<-'celltype'
allcolour=c("#FEFB98","#C4AFCB",'#C91889',"#96BF9E","#FDCD8A",'#7FC97F','#FEE491',"#DBB6AF",'#D63048')
umap_cc = colon@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = colon$celltype) 
p2<-ggplot(umap_cc,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
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
    legend.text = element_text(size=10), 
    legend.key.size=unit(1,'cm') ) +
  guides(color = guide_legend(override.aes = list(size=5)))+ 
  geom_segment(aes(x = min(umap_cc$UMAP_1)-2 , y = min(umap_cc$UMAP_2)-2 ,
                   xend = min(umap_cc$UMAP_1) +1, yend = min(umap_cc$UMAP_2)-2 ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+
  geom_segment(aes(x = min(umap_cc$UMAP_1)-2  , y = min(umap_cc$UMAP_2)-2  ,
                   xend = min(umap_cc$UMAP_1)-2 , yend = min(umap_cc$UMAP_2) + 1),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap_cc$UMAP_1) -0.5, y = min(umap_cc$UMAP_2) -3, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) +
  annotate("text", x = min(umap_cc$UMAP_1) -3, y = min(umap_cc$UMAP_2)-0.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)
ggsave('f6/atac_umap.pdf',p2,width = 8,height = 6)

###################
###dot tfs
###################
diffmarker<-readRDS('cc/data/scrna/tumormarker.rds')
tfs<-c('CEBPG','LEF1','SOX4','TCF7','TEAD4')
plotdata<-data.frame(matrix(0,nrow = 10,ncol = 4))
colnames(plotdata)<-c('TF','Avg_log2FC','P.adj','Type')
plotdata$Type<-rep(c('RNA','Chromvar'),each=5)
plotdata$TF<-rep(tfs,2)
for(i in 1:5){
  plotdata[i,2]<-rnamarker[plotdata[i,1],2]
  plotdata[i,3]<-rnamarker[plotdata[i,1],5]
}
for(i in 6:10){
  plotdata[i,2]<-rnamarker[plotdata[i,1],2]
  plotdata[i,3]<-rnamarker[plotdata[i,1],5]
}
p<-ggplot(plotdata, aes(x = Type, y = TF, color = P.adj)) + 
  geom_point(aes(size = Avg_log2FC)) + 
  scale_color_distiller(palette = "RdBu") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line=element_line(size = .3, colour="black"))
ggsave('f6/dot_rna_chrom.pdf',p,width = 3,height = 4)

###################
###metaprogram tfs
###################
motif_data<-readRDS('cc/chromvar/cc_motif_data.rds')
meta1_data<-motif_data[unique(meta1$peakName),c('CEBPG','LEF1','SOX4','TCF7','TEAD4')]
meta2_data<-motif_data[unique(meta2$peakName),c('CEBPG','LEF1','SOX4','TCF7','TEAD4')]
meta3_data<-motif_data[unique(meta3$peakName),c('CEBPG','LEF1','SOX4','TCF7','TEAD4')]
meta4_data<-motif_data[unique(meta4$peakName),c('CEBPG','LEF1','SOX4','TCF7','TEAD4')]
colSums(meta1_data)/nrow(meta1_data)
colSums(meta2_data)/nrow(meta2_data)
colSums(meta3_data)/nrow(meta3_data)
colSums(meta4_data)/nrow(meta4_data)
plotdata<-data.frame(matrix(0,ncol = 4,nrow = 20))
colnames(plotdata)<-c('Program','TF','Percent')
plotdata$Program<-rep(c('metaProgram1','metaProgram2','metaProgram3','metaProgram4'),each=5)
plotdata$TF<-rep(c('CEBPG','LEF1','SOX4','TCF7','TEAD4'),4)
plotdata$TF<-factor(plotdata$TF,levels = c('CEBPG','LEF1','SOX4','TCF7','TEAD4'))
plotdata$Percent<-c(0.03651685 ,0.30617978 ,0.19101124 ,0.29119850 ,0.11985019 ,
                    0.03707224 ,0.19961977 ,0.15304183 ,0.16920152 ,0.13498099,
                    0.04454343 ,0.15144766 ,0.13140312 ,0.14699332 ,0.10244989,
                    0.03755869 ,0.19014085 ,0.16666667 ,0.18544601 ,0.13380282)
plotdata$value<-c('**','***','***','***','*',
                  '**','***','***','***','***',
                  '**','**','ns','***','ns',
                  '*','***','***','***','**')
p<-ggplot(data = plotdata, aes(x = Program, y = Percent, fill =  TF,label=value)) +
  geom_bar(stat = 'identity', width = .7,position = position_dodge(0.85)) +
  theme_set(theme_bw())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black", size = 1 ))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(11))+
  geom_text(aes(label = value), position = position_dodge(0.9), vjust = -0.8)+
  ylim(c(0,0.35))+
  theme(axis.text.x = element_text(angle=30,vjust = 0.5,hjust = 0.5,size = 8 ))
ggsave('f6/meta_tf.pdf',p,height = 3.5,width = 4.5)

###################
###footprinting
###################
#plot
cc<-readRDS('cc/data/tf/cc.rds')
DefaultAssay(cc)<-'peaks'
Idents(cc)<-'celltype'
for(i in c('CEBPG','LEF1','SOX4','TCF7','TEAD4')){
  p<- PlotFootprint(cc, features = i)
  ggsave(paste0('f6/',gsub('::','',i),'.pdf'),p,width = 8,height = 6)
}
#normalized value
cebpg<-readRDS('cc/data/tf/tf_cebpg.rds')
lef1<-readRDS('cc/data/tf/tf_lef1.rds')
sox4<-readRDS('cc/data/tf/tf_sox4.rds')
tcf7<-readRDS('cc/data/tf/tf_tcf7.rds')
tead4<-readRDS('cc/data/tf/tf_tead4.rds')
tf_p<-rbind(cebpg,lef1,sox4,tcf7,tead4)
tf_p<-subset(tf_p,class=='Observed' )
df1<-data.frame(matrix(0,nrow = 5,ncol=2))
colnames(df1)<-c('TF','P_value')
df1$TF<-c('CEBPG','LEF1','SOX4','TCF7','TEAD4')
for(i in 1:5){
  t<-wilcox.test(subset(tf_p,feature==df1[i,1] & group=='Tumor cell')[,4],
                 subset(tf_p,feature==df1[i,1] & group=='Epithelial')[,4])
  df1[i,2]<-t$p.value
}
write.csv(df1,'tables/tf_250.csv',quote = F,row.names = F)

tf_p<-rbind(cebpg,lef1,sox4,tcf7,tead4)
tf_p<-subset(tf_p,class=='Observed' & abs(position)<=100)
df2<-data.frame(matrix(0,nrow = 5,ncol=2))
colnames(df2)<-c('TF','P_value')
df2$TF<-c('CEBPG','LEF1','SOX4','TCF7','TEAD4')
for(i in 1:5){
  t<-wilcox.test(subset(tf_p,feature==df2[i,1] & group=='Tumor cell')[,4],
                 subset(tf_p,feature==df2[i,1] & group=='Epithelial')[,4])
  df2[i,2]<-t$p.value
}
write.csv(df2,'tables/tf_250.csv',quote = F,row.names = F)

###################
###cell similarity
###################
CancerCell<-readRDS('cc/nmf/CancerCell.rds')
CancerCell.four.matrix <- GetAssayData(CancerCell, slot = "scale.data" ,assay.type = "SCT")
cluster.order <- sapply(c("A001-C-007_", "CRC1_", "CRC3_"), function(x){
  index <- grep(x, colnames(CancerCell.four.matrix))
  patient.matrix <- CancerCell.four.matrix[,index]
  corrMatrix <- (1- cor(patient.matrix, method="pearson"))
  hc <- hclust(as.dist(corrMatrix), method="ward.D")
  row_dend <- dendsort(hc)
  return(colnames(patient.matrix)[row_dend$order])
})
cluster.order <- unlist(cluster.order)
index <- match(cluster.order, colnames(CancerCell.four.matrix))
CancerCell.four.matrix <- CancerCell.four.matrix[,index]
cell.similarity <- cor(CancerCell.four.matrix, method="pearson")
saveRDS(cell.similarity, file = "cell.similarity.rds")
diag(cell.similarity) <- 0
diag(cell.similarity) <- max(cell.similarity)
pdf("f6/cell.similarity.pdf")
row_split <- gsub("_\\d+", "", names(cluster.order))
column_split <- gsub("_\\d+", "", names(cluster.order))
Heatmap(cell.similarity, cluster_rows = F, border = TRUE, row_split = row_split, column_split = column_split, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), width = unit(12, "cm"), height = unit(12, "cm"), cluster_columns = F, show_column_names = F, show_row_names = F, use_raster = T, heatmap_legend_param = list(title = "Correlation coeficient"))
dev.off()

################################
###pathway plot
################################
res1<-readRDS('cc/nmf/res1_enrich.rds')
res2<-readRDS('cc/nmf/res2_enrich.rds')
res3<-readRDS('cc/nmf/res3_enrich.rds')
res4<-readRDS('cc/nmf/res4_enrich.rds')
path1<-c('GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY','GOBP_NEGATIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY','GOBP_CELL_CELL_SIGNALING_BY_WNT',
         'GOBP_CANONICAL_WNT_SIGNALING_PATHWAY','KEGG_PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM')
path2<-c('KEGG_ADHERENS_JUNCTION','GOBP_CELL_GROWTH','GOBP_CELL_JUNCTION_ASSEMBLY','GOBP_ACTIN_FILAMENT_BASED_MOVEMENT',
         'GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY')
path3<-c('REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION','KEGG_RIBOSOME','REACTOME_SELENOAMINO_ACID_METABOLISM',
         'GOBP_CYTOPLASMIC_TRANSLATION','REACTOME_TRANSLATION')
path4<-c('HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','GOBP_CHROMOSOME_SEGREGATION',
         'GOBP_DNA_REPLICATION','HALLMARK_MITOTIC_SPINDLE')
res1_path1<-res1[path1,]
res2_path2<-res2[path2,]
res3_path3<-res3[path3,]
res4_path4<-res4[path4,]
plotdata<-rbind(res1_path1,res2_path2,res3_path3,res4_path4)
plotdata$Count<-ifelse(plotdata$Count>=15,15,plotdata$Count)
plotdata$pro<-c(rep('metaProgram1',5),rep('metaProgram2',5),rep('metaProgram3',5),rep('metaProgram4',5))
plotdata$Description<-factor(plotdata$Description,levels = plotdata$Description)
p = ggplot(data = plotdata,aes(x = pro, y = Description))+
  geom_point(aes(size = Count,color = p.adjust))+
  theme_bw()+
  scale_colour_gradient(low = "green",high = "red")+
  scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
  labs(x = "metaProgram",y = "",title = "Dotplot",
       color = expression(p.adjust),size = "Count")+
  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 11),legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 60,hjust = 1))
ggsave('f6/metaprogram.pdf',p,width = 6.5,height = 6)
