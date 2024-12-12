
setwd('tumor')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)

###################
###peak annotation
###################
#load data
cc_peak<-readRDS('cc/chipseeker/peakanno.rds')
cc_peak<-data.frame(type=unique(crc_peak$type),counts=as.vector(table(crc_peak$type)[unique(crc_peak$type)]),tumor=rep('CC',7))
oc_peak<-readRDS('oc/chipseeker/peakanno.rds')
oc_peak<-data.frame(type=unique(oc_peak$type),counts=as.vector(table(oc_peak$type)[unique(oc_peak$type)]),tumor=rep('OC',7))
ec_peak<-readRDS('ec/chipseeker/peakanno.rds')
ec_peak<-data.frame(type=unique(ec_peak$type),counts=as.vector(table(ec_peak$type)[unique(ec_peak$type)]),tumor=rep('EC',7))
rcc_peak<-readRDS('rcc/chipseeker/peakanno.rds')
rcc_peak<-data.frame(type=unique(rcc_peak$type),counts=as.vector(table(rcc_peak$type)[unique(rcc_peak$type)]),tumor=rep('RCC',7))
lc_peak<-readRDS('lc/chipseeker/peakanno.rds')
lc_peak<-data.frame(type=unique(lc_peak$type),counts=as.vector(table(lc_peak$type)[unique(lc_peak$type)]),tumor=rep('LC',7))
bc_peak<-readRDS('bc/chipseeker/peakanno.rds')
bc_peak<-data.frame(type=unique(bc_peak$type),counts=as.vector(table(bc_peak$type)[unique(bc_peak$type)]),tumor=rep('BC',7))
combined_peak<-readRDS('chipseeker/peakanno.rds')
combined_peak<-data.frame(type=unique(combined_peak$type),counts=as.vector(table(combined_peak$type)[unique(combined_peak$type)]),tumor=rep('All_sample',7))
#plot
plotdf<-rbind(bc_peak,bcc_peak,crc_peak,oc_peak,ec_peak,lc_peak,plc_peak,rcc_peak,combined_peak)
p<-ggplot(plotdf, aes
          (tumor, weight = counts, fill = type)) +
  geom_bar(color = "white", width = .7, position = 'stack') +
  labs( y = 'Counts') +
  scale_fill_brewer(palette = "Set3")+
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
ggsave('f2/Chipseeker_anno.pdf',p,height = 5.5,width = 3.5)

################
###peak barplot
################
peaks_bcc<-readRDS('bcc/data/scatac/peaks_bcc.rds')
peaks_cc<-readRDS('cc/data/scatac/peaks_cc.rds')
peaks_ec<-readRDS('ec/data/scatac/peaks_ec.rds')
peaks_oc<-readRDS('oc/data/scatac/peaks_oc.rds')
peaks_rcc<-readRDS('rcc/data/scatac/peak_rcc.rds')
peaks_plc<-readRDS('plc/data/scatac/peaks_plc.rds')
peaks_lc<-readRDS('lc/data/scatac/peak_lc.rds')
peaks_bc<-readRDS('bc/data/scatac/peak_bc.rds')
bcc_df<-data.frame(table(peaks_bcc$cluster))
bcc_df$Tumor<-'BCC'
cc_df<-data.frame(table(peaks_cc$cluster))
cc_df$Tumor<-'CC'
ec_df<-data.frame(table(peaks_ec$cluster))
ec_df$Tumor<-'EC'
oc_df<-data.frame(table(peaks_oc$cluster))
oc_df$Tumor<-'OC'
plc_df<-data.frame(table(peaks_plc$cluster))
plc_df$Tumor<-'PLC'
rcc_df<-data.frame(table(peaks_rcc$cluster))
rcc_df$Tumor<-'RCC'
lc_df<-data.frame(table(peaks_lc$cluster))
lc_df$Tumor<-'LC'
bc_df<-data.frame(table(peaks_bc$cluster))
bc_df$Tumor<-'BC'
plotdata<-rbind(bc_df,bcc_df,cc_df,ec_df,lc_df,oc_df,plc_df,rcc_df)
colnames(plotdata)<-c('Celltype','Number','Tumor')
p<-ggplot(plotdata,aes(x=Celltype,y=Number,fill=Tumor))+
  geom_bar(stat = 'identity')+
  theme_test()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1,size = 11)) + 
  facet_grid(~ Tumor, scales = "free_x",space = 'free')
ggsave('f2/Diff_peak.pdf',p,width = 12,height = 4)

###################
###coverageplot
###################
snp_cc<-readRDS('cc/gwas/download/snp.rds')
tumor_atac<-readRDS('cells/tumor/tumor.rds')
fibro_atac<-readRDS('cells/fibro/fibro.rds')
cc<-readRDS('cc/data/scatac/cc.rds')
Idents(tumor_atac)<-'tumorcell'
Idents(fibro_atac)<-'tumorcell'
Idents(cc)<-'celltype'
p1<-CoveragePlot(cc,'chr10-112518672-112519744',extend.downstream = 10000,extend.upstream = 10000)
p2<-CoveragePlot(cc,'chr3-64266065-64266806',extend.upstream = 10000,extend.downstream = 10000)
p3<-CoveragePlot(tumor_atac,'chr10-112518672-112519744',extend.downstream = 10000,extend.upstream = 10000)
p4<-CoveragePlot(fibro_atac,'chr3-64266065-64266806',extend.upstream = 10000,extend.downstream = 10000)
ggsave('f2/tcf7l2_cc.pdf',p1,width = 6,height = 4.5)
ggsave('f2/tcf7l2_tumor.pdf',p3,width = 6,height = 4.5)
ggsave('f2/prickle2_cc.pdf',p2,width = 6,height = 4.5)
ggsave('f2/prickle2_fibro.pdf',p4,width = 6,height = 4.5)

###################
###ldsc
###################
###cc
files<-list.files('/cc/gwas/ldsc')
df_enrich<-data.frame()
for(i in files){
  tmp_files<-list.files(paste0('cc/gwas/ldsc/',i))
  tmp_files<-tmp_files[grep('results',tmp_files)]
  if(grepl('results',i)){
    tmp<-read.delim(paste0('cc/gwas/ldsc/',i))
    tmp_df<-data.frame(cell=tmp[,1],disease=rep(strsplit(i,'\\.')[[1]][1],nrow(tmp)),pvalue=tmp[,4])
    df_enrich<-rbind(df_enrich,tmp_df)}
}
plotdata<-df_enrich %>% mutate(text = case_when(
  pvalue < 0.0001 ~ "****",pvalue<0.001 ~ "***",pvalue<0.01 ~"**",pvalue<0.05 ~"*"))
plotdata<-plotdata[plotdata$cell!='all',]
plotdata$log10P<- -log10(plotdata$pvalue)
p1<-ggplot(plotdata,aes(y=disease,x=log10P,color=cell))+
  geom_point()+
  geom_vline(xintercept = 1.30103,linetype='dashed')+
  scale_color_manual(values=c("#FEFB98","#C4AFCB","#96BF9E","#FDCD8A","#FEE491","#DBB6AF",'#D63048'))+
  theme_test()
ggsave('f2/cc_ldsc.pdf',p1,width =4.3 ,height = 3)

###oc
files<-list.files('oc/gwas/ldsc')
df_enrich<-data.frame()
for(i in files){
  tmp_files<-list.files(paste0('oc/gwas/ldsc/',i))
  tmp_files<-tmp_files[grep('results',tmp_files)]
  if(grepl('results',i)){
    tmp<-read.delim(paste0('oc/gwas/ldsc/',i))
    tmp_df<-data.frame(cell=tmp[,1],disease=rep(strsplit(i,'\\.')[[1]][1],nrow(tmp)),pvalue=tmp[,4])
    df_enrich<-rbind(df_enrich,tmp_df)}
}
plotdata<-df_enrich %>% mutate(text = case_when(
  pvalue < 0.0001 ~ "****",pvalue<0.001 ~ "***",pvalue<0.01 ~"**",pvalue<0.05 ~"*"))
plotdata<-plotdata[plotdata$cell!='all',]
plotdata$log10P<- -log10(plotdata$pvalue)
p2<-ggplot(plotdata,aes(y=disease,x=log10P,color=cell))+
  geom_point()+
  geom_vline(xintercept = 1.30103,linetype='dashed')+
  scale_color_manual(values=c("#C4AFCB","#96BF9E","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048'))+
  theme_test()
ggsave('f2/oc_ldsc.pdf',p2,width = 4.3,height = 3)

###bc
files<-list.files('bc/gwas/ldsc')
df_enrich<-data.frame()
for(i in files){
  tmp_files<-list.files(paste0('bc/gwas/ldsc/',i))
  tmp_files<-tmp_files[grep('results',tmp_files)]
  if(grepl('results',i)){
    tmp<-read.delim(paste0('bc/gwas/ldsc/',i))
    tmp_df<-data.frame(cell=tmp[,1],disease=rep(strsplit(i,'\\.')[[1]][1],nrow(tmp)),pvalue=tmp[,4])
    df_enrich<-rbind(df_enrich,tmp_df)}
}
plotdata<-df_enrich %>% mutate(text = case_when(
  pvalue < 0.0001 ~ "****",pvalue<0.001 ~ "***",pvalue<0.01 ~"**",pvalue<0.05 ~"*"))
plotdata<-plotdata[plotdata$cell!='all',]
plotdata$log10P<- -log10(plotdata$pvalue)
p3<-ggplot(plotdata,aes(y=disease,x=log10P,color=cell))+
  geom_point()+
  geom_vline(xintercept = 1.30103,linetype='dashed')+
  scale_color_manual(values=c("#C4AFCB","#96BF9E","#FDCD8A","#7FC97F",'#FEE491',"#DBB6AF",'#D63048'))+
  theme_test()
ggsave('f2/bc_ldsc.pdf',p3,width = 4.3,height = 3)
