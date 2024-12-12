
setwd('tumor')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(ragg)
library(ggpp)

###################
###gpl plots
###################
bcc_ll<-readRDS('bcc/peakgenelinks/links_end.rds')
bcc_ll<-subset(bcc_ll,peak_promoter==1)
cc_ll<-readRDS('cc/peakgenelinks/links_end.rds')
cc_ll<-subset(cc_ll,peak_promoter==1)
ec_ll<-readRDS('ec/peakgenelinks/links_end.rds')
ec_ll<-subset(ec_ll,peak_promoter==1)
oc_ll<-readRDS('oc/peakgenelinks/links_end.rds')
oc_ll<-subset(oc_ll,peak_promoter==1)
plc_ll<-readRDS('plc/peakgenelinks/links_end.rds')
plc_ll<-subset(plc_ll,peak_promoter==1)
rcc_ll<-readRDS('rcc/peakgenelinks/links_end.rds')
rcc_ll<-subset(rcc_ll,peak_promoter==1)
lc_ll<-readRDS('lc/peakgenelinks/links_end.rds')
lc_ll<-subset(lc_ll,peak_promoter==1)
bc_ll<-readRDS('bc/peakgenelinks/links_end.rds')
bc_ll<-subset(bc_ll,peak_promoter==1)

peaks_bcc<-readRDS('bcc/data/scatac/peaks_bcc.rds')
peaks_bcc<-subset(peaks_bcc,gene %in% bcc_ll$peakName)
peaks_bc<-readRDS('bc/data/scatac/peaks_bc.rds')
peaks_bc<-subset(peaks_bc,gene %in% bc_ll$peakName)
peaks_cc<-readRDS('cc/data/scatac1/peaks_cc.rds')
peaks_cc<-subset(peaks_cc,gene %in% cc_ll$peakName)
peaks_plc<-readRDS('plc/data/scatac/peaks_plc.rds')
peaks_plc<-subset(peaks_plc,gene %in% plc_ll$peakName)
peaks_lc<-readRDS('lc/data/scatac/peaks_lc.rds')
peaks_lc<-subset(peaks_lc,gene %in% lc_ll$peakName)
peaks_rcc<-readRDS('rcc/data/scatac2/peaks_rcc.rds')
peaks_rcc<-subset(peaks_rcc,gene %in% rcc_ll$peakName)
peaks_ec<-readRDS('ec/data/scatac/peaks_ec.rds')
peaks_ec<-subset(peaks_ec,gene %in% ec_ll$peakName)
peaks_oc<-readRDS('oc/data/scatac/peaks_oc.rds')
peaks_oc<-subset(peaks_oc,gene %in% oc_ll$peakName)

df_list<-list()
for(i in 1:8){
  tumor<-c('BC','BCC','CC','EC','LC','OC','PLC','RCC')[i]
  peaks<-list(peaks_bc,peaks_bcc,peaks_cc,peaks_ec,peaks_lc,peaks_oc,peaks_plc,peaks_rcc)[[i]]
  tmp<-list()
  for(j in unique(peaks$cluster)){
    tmp[[j]]<-subset(peaks,cluster==j)[,'gene']
  }
  df_list[[tumor]]<-tmp
}

#correlation
meta<-readRDS('combined/conserved/meta.rds')
plotdata<-data.frame(matrix(0,nrow = 8,ncol = 3))
colnames(plotdata)<-c('Type','Cell','GPL')
plotdata$Type<-c('BCC','CC','EC','OC','PLC','RCC','LC','BC')
plotdata$Peak<-c(12214,27409,29250,13031,14231,26925,23868,30000)
plotdata$GPL<-c(11955,38292,39640,15202,10140,13368,18976,20109)
p<-ggscatter(plotdata, x = "Cell", y = "GPL",
             add = "reg.line", conf.int = TRUE,    
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson", 
           label.x = 15000, label.y = 50000)
ggsave('f4/correlation_gpl.pdf',p,width = 4,height = 3)

#barplot
bardata<-c(20109,11955,38292,39640,18976,15202,10140,13368)
names(bardata)<-c('BC','BCC','CC','EC','LC','OC','PLC','RCC')
p<-barplot(bardata,col = rgb(0.2,0.4,0.6,0.6),horiz = T,las=1)
ggsave('f4/barplot_gpl.pdf',p,width = 5,height = 4)

#cell type gpl
cells<-c('Tumor cell','T cell','Myeloid cell','Endothelial','Myofibroblast','Fibroblast',
         'B cell','Plasma')
plotdata<-data.frame(row.names = c('Cancer','Celltype','N'))
for(i in 1:8){
  x<-c('BCC','CC','EC','OC','PLC','RCC','LC','BC')[i]
  tmp<-df_list[[i]]
  for(j in names(tmp)){
    if(j %in% cells){
      tmp_plot<-data.frame(c(x,j,length(unique(tmp[[j]]))))
      plotdata<-cbind(plotdata,tmp_plot)
    }
  }
}
plotdata<-data.frame(t(plotdata))
plotdata$N<-as.numeric(plotdata$N)
p<-ggplot(plotdata,aes(x=Celltype,y=N))+
  geom_boxplot(aes(fill=Celltype),fill=NA)+
  geom_jitter(aes(color=Cancer),size=3)+
  scale_color_manual(values=c("#C91889","#C4AFCB",'#96BF9E',"#A75D2B",'#769AA8',"#FDCD8A","#D63048","#666666"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,size = 11))
ggsave('f4/celltype_gpl.pdf',p,width = 5,height = 3)

###################
###heatmap
###################
data<-readRDS('combined/conserved/data.rds')
p1<-data[[1]]
p2<-data[[2]]
p3<-data[[3]]
plotdata1<-data[[4]]
plotdata2<-data[[5]]
plotdata3<-data[[6]]
anno_col_all<-data[[7]]
plotdata_all<-rbind(plotdata1[p1$tree_row$order,rownames(anno_col_all)],
                    plotdata2[p2$tree_row$order,rownames(anno_col_all)],
                    plotdata3[p3$tree_row$order,rownames(anno_col_all)])
pall<-pheatmap(plotdata_all,cluster_cols = F,show_colnames = F,
                 show_rownames = F,cluster_rows = F,annotation_names_col = T,
                 annotation_col = anno_col_all,
                 gaps_row=c(),color = colorRampPalette(c('#48025b','#1f938a','#e7dd1b'))(100),
                 gaps_col = c(1884,3227,3291,3350,3771,3812,4094,4581,4969))
ggsave('f4/pheatmap_all.png',pall,width = 8,height = 10)

###################
###tumor great
###################
great_tumor<-read.delim('combined/conserved/diffpeak/tumor_greatExportAll.tsv',skip = 3)
great_tumor$Desc<-factor(great_tumor$Desc,levels = unique(great_tumor$Desc))
p<-ggplot(great_tumor[1:10,],aes(SetCov,Desc)) + 
  geom_point(aes(size=ObsRegions,color=BinomBonfP)) +
  scale_color_gradient(low="red",high = "green") + 
  labs(color="BinomBonfP",size="ObsRegions",
       x="SetCov",y="Top 10 GO terms",title="GREAT enrichment") + 
  theme_bw()
ggsave('f4/tumor_great.pdf',p,width = 5,height = 3.5)

###################
###overlapping paths
###################
tumor_res<-readRDS('combined/conserved/tumor_res_list.rds')
bcc_tumor<-subset(tumor_res[[1]],pvalue<0.05)
cc_tumor<-subset(tumor_res[[2]],pvalue<0.05)
ec_tumor<-subset(tumor_res[[3]],pvalue<0.05)
oc_tumor<-subset(tumor_res[[4]],pvalue<0.05)
plc_tumor<-subset(tumor_res[[5]],pvalue<0.05)
rcc_tumor<-subset(tumor_res[[6]],pvalue<0.05)
lc_tumor<-subset(tumor_res[[7]],pvalue<0.05)
bc_tumor<-subset(tumor_res[[8]],pvalue<0.05)
paths<-Reduce(intersect,list(bcc_tumor$ID,cc_tumor$ID,ec_tumor$ID,oc_tumor$ID,
                             plc_tumor$ID,rcc_tumor$ID,lc_tumor$ID,bc_tumor$ID))
circle_data<-data.frame(GO=rep(paths,8),
                        label=rep( 1:8,each=12),
                        Pvalue=0,
                        Count=0,
                        log10_Pvalue=0,
                        group=rep(c('A','B','C','D','E','F','G','H','I','J','K','L'),8))
for(i in 1:nrow(circle_data)){
  tmp<-list(bcc_tumor,cc_tumor,ec_tumor,oc_tumor,plc_tumor,rcc_tumor,lc_tumor,bc_tumor)[[circle_data[i,2]]]
  circle_data[i,3]<-tmp[circle_data[i,1],5]
  circle_data[i,4]<-tmp[circle_data[i,1],9]
}
circle_data$label<-rep( c('BCC','CC','EC','OC','PLC','RCC','LC','BC'),each=12)
datagroup <- circle_data$group %>% unique()
allplotdata <- tibble('group' = datagroup,
                      'label' = paste0('empty_individual_', seq_along(datagroup)),
                      'Pvalue' = 0,
                      'Count'=0) %>% 
  bind_rows(circle_data) %>% arrange(group) %>% mutate(xid = 1:n()) %>% 
  mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  mutate(angle = ifelse(angle < -90, angle+180, angle)) 
firstxid <- which(str_detect(allplotdata$label, pattern = "empty_individual")) 
segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$Count), to = max(allplotdata$Count), 10),
                 'coordytext' = as.character(round(coordylocation, 2)),
                 'x' = 1)
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)
p1<-ggplot() +   coord_polar() +  theme_void() +
  geom_bar(data = allplotdata, aes(x = xid, y = Count, fill = Pvalue), stat = 'identity') + 
  scale_fill_gradient2(low = "yellow",mid='#DBB6AF',high = "#4166AD",midpoint = 4.202505e-04)+
  geom_text(data = allplotdata %>% filter(!str_detect(label, pattern = "empty_individual")), 
            aes(x = xid, label = label, y = Count+10, angle = angle, hjust = hjust),
            color="black", fontface="bold",alpha=0.6, size=1.5) + 
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -5, yend=-5) + 
  geom_text(data = segment_data, aes(x = labelx, label = label), y = -10,size=3) + 
  geom_segment(data = griddata, 
               aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy),
               colour = "grey", alpha=0.8, size=0.6) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-30,100)) +
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color="black", size=2 , angle=0, fontface="bold") 
ggsave('f4/circleplot.pdf',p1,width = 9,height = 7)

table_plot<-data.frame(matrix(0,nrow = 12,ncol = 3))
table_plot$X1<-c('A','B','C','D','E','F','G','H','I','J','K','L')
table_plot$X2<-paths
table_plot$X3<-c('cell-cell junction organization','cell-cell junction assembly',
                 'apical junction assembly','tight junction organization',
                 'tight junction assembly','bicellular tight junction assembly',
                 'neural tube development','tube formation',
                 'epithelial tube formation','epithelial cell development',
                 'regulation of cell morphogenesis','cell differentiation involved in embryonic placenta development')
p2<-ggtexttable(table_plot, rows = NULL,cols = NULL,
                theme = ttheme(
                  colnames.style = colnames_style(color = "white", fill = "#8cc257")
                )
)
ggsave('f4/paths.pdf',p2,width = 6,height = 9)

###################
##motif enrich
###################
combined_atac<-readRDS('combined/chromvar/combined.rds')
peak_all<-StringToGRanges(rownames(combined_atac))
tumor_peak<-ChIPseeker::readPeakFile('combined/conserved/diffpeak/tumor_peak.bed')
tumor_peak<-peak_all[unique(findOverlaps(peak_all,tumor_peak)@from)]
t_peak<-ChIPseeker::readPeakFile('combined/conserved/diffpeak/t_peak.bed')
t_peak<-peak_all[unique(findOverlaps(peak_all,t_peak)@from)]
mye_peak<-ChIPseeker::readPeakFile('combined/conserved/diffpeak/mye_peak.bed')
mye_peak<-peak_all[unique(findOverlaps(peak_all,mye_peak)@from)]
fibro_peak<-ChIPseeker::readPeakFile('combined/conserved/diffpeak/fibro_peak.bed')
fibro_peak<-peak_all[unique(findOverlaps(peak_all,fibro_peak)@from)]
endo_peak<-ChIPseeker::readPeakFile('combined/conserved/diffpeak/endo_peak.bed')
endo_peak<-peak_all[unique(findOverlaps(peak_all,endo_peak)@from)]
myof_peak<-ChIPseeker::readPeakFile('combined/conserved/diffpeak/myof_peak.bed')
myof_peak<-peak_all[unique(findOverlaps(peak_all,myof_peak)@from)]
tumor_motif<-FindMotifs(combined_atac,GRangesToString(tumor_peak))
t_motif<-FindMotifs(combined_atac,GRangesToString(t_peak))
mye_motif<-FindMotifs(combined_atac,GRangesToString(mye_peak))
fibro_motif<-FindMotifs(combined_atac,GRangesToString(fibro_peak))
endo_motif<-FindMotifs(combined_atac,GRangesToString(endo_peak))
myof_motif<-FindMotifs(combined_atac,GRangesToString(myof_peak))
#plot
plot_tumor<-tumor_motif[1:15,]
plot_t<-t_motif[1:15,]
plot_mye<-mye_motif[1:15,]
plot_fibro<-fibro_motif[1:15,]
plot_endo<-endo_motif[1:15,]
plot_myof<-myof_motif[1:15,]
plotdata<-data.frame(matrix(0,nrow = 90,ncol = 4))
colnames(plotdata)<-c('Motif','CellType','Rank','P_value')
plotdata$Motif<-c(plot_tumor$motif.name,plot_t$motif.name,plot_mye$motif.name,
                  plot_fibro$motif.name,plot_endo$motif.name,plot_myof$motif.name)
plotdata$CellType<-rep(c('Tumor cell','T cell','Myeloid cell','Fibroblast','Endothelial','Myofibroblast'),each=15)
plotdata$Rank<-rep(1:15,6)
plotdata$P_value<-c(plot_tumor$pvalue,plot_t$pvalue,plot_mye$pvalue,
                    plot_fibro$pvalue,plot_endo$pvalue,plot_myof$pvalue)
plotdata$tf<-plotdata$Motif
plotdata$Motif<-paste0(plotdata$CellType,'_',plotdata$Motif)
labels<-plotdata$tf
names(labels)<-plotdata$Motif
plotdata$Motif<-factor(plotdata$Motif,levels = unique(plotdata$Motif))
p<-ggplot(plotdata, aes(x = Motif, y = Rank, color = P_value)) + 
  geom_point() + 
  scale_color_distiller(palette = "RdBu") + 
  #scale_size(limits = c(1, 6), range = c(0.5, 3.5)) + 
  theme_classic() + 
  scale_x_discrete("Motif")+
  scale_x_discrete("Motif", labels = labels)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.line=element_line(size = .3, colour="black")) + 
  facet_grid(~ CellType, scales = "free_x",space = 'free')
ggsave('f4/dot_motif.pdf',p,width = 13,height = 4)

###################
###klf4/fos/spib/elk4
###################
combined_atac<-readRDS('combined/scatac/combined.rds')
p1<-MotifPlot(combined_atac,c('MA0039.4'))
p2<-MotifPlot(combined_atac,c('MA0081.2'))
p3<-MotifPlot(combined_atac,c('MA0476.1'))
p4<-MotifPlot(combined_atac,c('MA0076.2'))
ggsave('f4/klf4_motif.pdf',p1,width = 4,height = 2)
ggsave('f4/spib_motif.pdf',p2,width = 4,height = 2)
ggsave('f4/fos_motif.pdf',p3,width = 4,height = 2)
ggsave('f4/elk4_motif.pdf',p4,width = 4,height = 2)

###################
###klf4-related pathways
###################
gene_go<-readRDS('combined/conserved/klf4_res.rds')
p<-ggplot(gene_go[c('GO:0010464','GO:2001053','GO:0097152',"GO:0002053",'GO:0010463','GO:0014031',
                    'GO:0048762','GO:0001837','GO:0010718','GO:0060485'),],aes(Count,Description)) + 
  geom_point(aes(size=Count,color=pvalue)) +
  scale_color_gradient(low="red",high = "green") + 
  labs(color="P_value",size="Count",
       x="Count",y="Top 10 GO terms",title="GO enrichment") + 
  theme_bw()+
  scale_y_discrete(labels=function(x) str_wrap(x, width=20))
ggsave('f4/klf4_emt.pdf',p,width = 3.5,height = 5)
