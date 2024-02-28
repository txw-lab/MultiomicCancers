
setwd('data/pancancer_atac/tumor/combined')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)
library(RColorBrewer)

#######################
###cell similarity
#######################
#load data
CancerCell.four.matrix<-readRDS('NMF/CancerCell.four.matrix.rds')
cluster.order <- sapply(c("T1_", "T2_", "T3_", "T4_"), function(x){
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
diag(cell.similarity) <- 0
diag(cell.similarity) <- max(cell.similarity)
#plot
pdf("f4/cell.similarity.pdf")
row_split <- gsub("_\\d+", "", names(cluster.order))
column_split <- gsub("_\\d+", "", names(cluster.order))
Heatmap(cell.similarity, cluster_rows = F, border = TRUE, row_split = row_split, column_split = column_split, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), width = unit(12, "cm"), height = unit(12, "cm"), cluster_columns = F, show_column_names = F, show_row_names = F, use_raster = T, heatmap_legend_param = list(title = "Correlation coeficient"))
dev.off()

####################
###program cluster
####################
program.score<-readRDS('NMF/program.score.rds')
pdf("f4/programScore.hierarchical.clustering.pdf")
singleScore.cluster <- function(program.score, cluster.order, row_split){
  a <- t(program.score)
  index <- match(cluster.order, colnames(a))
  a <- a[,na.omit(index)]
  p1 <- Heatmap(a, cluster_rows = T, border = TRUE, row_split = row_split, clustering_distance_rows = "pearson", clustering_method_rows = "ward.D", cluster_column = F, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), cluster_columns = F, show_column_names = F, show_row_names = T, use_raster = T, width = unit(8, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Expression score"))
  print(p1)
}
res <- singleScore.cluster(program.score = T1.program.score, cluster.order = cluster.order, row_split = 5)
res <- singleScore.cluster(program.score = T2.program.score, cluster.order = cluster.order, row_split = 5)
res <- singleScore.cluster(program.score = T3.program.score, cluster.order = cluster.order, row_split = 5)
res <- singleScore.cluster(program.score = T4.program.score, cluster.order = cluster.order, row_split = 5)
Heatmap(mean.corrM, cluster_rows = T, cluster_columns = T, border = TRUE, row_split = 5, column_split = 5, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_column_names = T, show_row_names = T, use_raster = T, width = unit(6, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Mean correlation"))
Heatmap(mean.corrM, cluster_rows = T, cluster_columns = T, clustering_distance_rows = "pearson", clustering_method_rows = "ward.D", clustering_distance_columns = "pearson", clustering_method_columns = "ward.D", border = TRUE, row_split = 5, column_split = 5, row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_column_names = T, show_row_names = T, use_raster = T, width = unit(6, "cm"), height = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = "Mean correlation"))
dev.off()

################################
###metaprogram
################################
#load data
res1<-readRDS('NMF/res1_enrich.rds')
res2<-readRDS('NMF/res2_enrich.rds')
res3<-readRDS('NMF/res3_enrich.rds')
path1<-c('GOBP_ESTABLISHMENT_OF_CELL_POLARITY','GOBP_FOCAL_ADHESION_ASSEMBLY','GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY_PLANAR_CELL_POLARITY_PATHWAY',
         'BIOCARTA_CCR3_PATHWAY','GOBP_CELL_SUBSTRATE_JUNCTION_ORGANIZATION','GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY')
path2<-c('HALLMARK_G2M_CHECKPOINT','REACTOME_DNA_REPAIR','GOBP_DNA_REPLICATION','GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION',
         'GOBP_DNA_DEPENDENT_DNA_REPLICATION','KEGG_MISMATCH_REPAIR')
path3<-c('GOBP_CYTOPLASMIC_TRANSLATION','REACTOME_CELLULAR_RESPONSE_TO_STARVATION','GOBP_RIBOSOME_BIOGENESIS',
         'REACTOME_SIGNALING_BY_EGFR_IN_CANCER','REACTOME_REGULATION_OF_BACH1_ACTIVITY','REACTOME_SIGNALING_BY_VEGF')
res1_path1<-res1[path1,]
res2_path2<-res2[path2,]
res3_path3<-res3[path3,]
#plot
plotdata<-rbind(res1_path1,res2_path2,res3_path3)
plotdata$Count<-ifelse(plotdata$Count>=15,15,plotdata$Count)
plotdata$pro<-c(rep('metaProgram1',6),rep('metaProgram2',6),rep('metaProgram3',6))
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
        axis.text.x = element_text(angle = 90))
ggsave('f4/metaprogram.pdf',p,width = 8.5,height = 6)

#######################
###tf gene
#######################
#load data
meta.Signature<-readRDS('NMF/meta.Signature.rds')
df<-readRDS('peak_gene_link/archr/tumor_epi/df.rds')
tfs<-readRDS('diff/diff_tfs.rds')
tfs<-rownames(tfs)
#meta1
genes1<-intersect(names(df),meta.Signature[[1]])
subdf<-df[genes1]
tfs<-c('LEF1','SOX4','TCF7','TEAD4')
gg<-data.frame(matrix(0,nrow = length(subdf),ncol = length(tfs)))
colnames(gg)<-tfs
rownames(gg)<-names(subdf)
for(i in rownames(gg)){
  for(j in colnames(gg)){
    for(k in 1:length(unlist(subdf[[i]][2]))){
      if(grepl(j,unlist(subdf[[i]][2])[k])){
        gg[i,j]=1
        break
      }
    }
  }
}
# colSums(gg)
# LEF1  SOX4  TCF7 TEAD4 
# 69    61    73    59 
# LEF1      SOX4      TCF7     TEAD4 
# 0.8214286 0.7261905 0.8690476 0.7023810 

#meta2
genes2<-intersect(names(df),meta.Signature[[2]])
subdf<-df[genes2]
tfs<-c('LEF1','SOX4','TCF7','TEAD4')
gg<-data.frame(matrix(0,nrow = length(subdf),ncol = length(tfs)))
colnames(gg)<-tfs
rownames(gg)<-names(subdf)
for(i in rownames(gg)){
  for(j in colnames(gg)){
    for(k in 1:length(unlist(subdf[[i]][2]))){
      if(grepl(j,unlist(subdf[[i]][2])[k])){
        gg[i,j]=1
        break
      }
    }
  }
}
# colSums(gg)
# LEF1  SOX4  TCF7 TEAD4 
# 26    22    28    22 
# LEF1      SOX4      TCF7     TEAD4 
# 0.7222222 0.6111111 0.7777778 0.6111111 

#meta3
genes3<-intersect(names(df),meta.Signature[[3]])
subdf<-df[genes3]
tfs<-c('LEF1','SOX4','TCF7','TEAD4')
gg<-data.frame(matrix(0,nrow = length(subdf),ncol = length(tfs)))
colnames(gg)<-tfs
rownames(gg)<-names(subdf)
for(i in rownames(gg)){
  for(j in colnames(gg)){
    for(k in 1:length(unlist(subdf[[i]][2]))){
      if(grepl(j,unlist(subdf[[i]][2])[k])){
        gg[i,j]=1
        break
      }
    }
  }
}
# colSums(gg)
# LEF1  SOX4  TCF7 TEAD4 
# 32    31    39    25 
# LEF1      SOX4      TCF7     TEAD4 
# 0.6666667 0.6458333 0.8125000 0.5208333 

#plot
plotdata<-data.frame(matrix(0,ncol = 3,nrow = 12))
colnames(plotdata)<-c('Program','TF','Percent')
plotdata$Program<-rep(c('metaProgram1','metaProgram2','metaProgram3'),each=4)
plotdata$TF<-rep(c('LEF1','SOX4','TCF7','TEAD4'),3)
plotdata$TF<-factor(plotdata$TF,levels = c('LEF1','SOX4','TCF7','TEAD4'))
plotdata$Percent<-c(0.8214286,0.7261905,0.8690476,0.7023810 ,
                    0.7222222,0.6111111,0.7777778,0.6111111 ,
                    0.6666667,0.6458333,0.8125000,0.5208333)
p<-ggplot(data = plotdata, aes(x = Program, y = Percent, fill =  TF)) +
  geom_bar(stat = 'identity', width = .7,position = position_dodge(0.85)) +
  theme_set(theme_bw())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black", size = 1, ))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(11)) 
ggsave('f4/tf_percent.pdf',p,height = 3,width = 6)
