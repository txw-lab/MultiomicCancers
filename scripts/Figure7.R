
setwd('data/pancancer_atac/tumor')
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

###############################
###gwas snp intersect with peak
###############################
library(ggforce)
big.circle<-data.frame(x=0,y=0,radius=1)
small.circle<-data.frame(x=0,y=-0.2,radius=0.5)
p1<-ggplot()+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=0.5,data=big.circle,fill='#F8766D')+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=1,data=small.circle,fill='#F8766D')+
  geom_text(aes(x=0,y=-0.2,label='219'))+
  geom_text(aes(x=0,y=0.6,label='1197'))+
  theme_void()+
  coord_fixed()
p2<-ggplot()+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=0.5,data=big.circle,fill='#7CAE00')+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=1,data=small.circle,fill='#7CAE00')+
  geom_text(aes(x=0,y=-0.2,label='1'))+
  geom_text(aes(x=0,y=0.6,label='117'))+
  theme_void()+
  coord_fixed()
p3<-ggplot()+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=0.5,data=big.circle,fill='#00BFC4')+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=1,data=small.circle,fill='#00BFC4')+
  geom_text(aes(x=0,y=-0.2,label='60'))+
  geom_text(aes(x=0,y=0.6,label='506'))+
  theme_void()+
  coord_fixed()
p4<-ggplot()+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=0.5,data=big.circle,fill='#C77CFF')+
  geom_circle(aes(x0=x,y0=y,r=radius),alpha=1,data=small.circle,fill='#C77CFF')+
  geom_text(aes(x=0,y=-0.2,label='4'))+
  geom_text(aes(x=0,y=0.6,label='45'))+
  theme_void()+
  coord_fixed()
ggsave('f7/Venn_plots.pdf',p1+p2+p3+p4,height = 10,width = 10)

###################
###cc ldsc
###################
#load ldsc result
files<-list.files('cc/gwas/enrichresult')
df_enrich<-data.frame()
for(i in files){
  tmp_files<-list.files(paste0('cc/gwas/enrichresult/',i))
  tmp_files<-tmp_files[grep('results',tmp_files)]
  for(j in tmp_files){
    tmp<-read.delim(paste0('cc/gwas/enrichresult/',i,'/',j))
    info<-strsplit(j,'\\.')[[1]]
    tmp_df<-data.frame(cell=info[1],disease=info[2],enrichment=tmp[1,5],Pvalue=tmp[1,7])
    df_enrich<-rbind(df_enrich,tmp_df)}
}
#plot
library(dplyr)
library(ggplot2)
plotdata<-df_enrich %>% mutate(text = case_when(
  Pvalue<0.0001 ~ "****",Pvalue<0.001 ~ "***",Pvalue<0.01 ~"**",Pvalue<0.05 ~"*"))
plotdata$enrichment<-ifelse(plotdata$enrichment<(-50),-50,plotdata$enrichment)
plotdata$enrichment<-ifelse(plotdata$enrichment>50,50,plotdata$enrichment)
p1<-ggplot(plotdata, aes(cell, disease)) + 
  geom_tile(aes(fill = Pvalue), colour = "black", size = 0.2)+
  scale_fill_gradient2(low = '#1A5592',mid = "white",high = "#B83D3D") + 
  geom_text(aes(label=text),col ="black",size = 5) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(size = 10,color = 'black',angle = 90, hjust = 1), 
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.position = 'top') + 
  scale_y_discrete(position = "right")
ggsave('f7/cc_ldsc.pdf',p1,height = 6,width = 8)

###################
###oc ldsc
###################
#load ldsc result
files<-list.files('oc/gwas/enrichresult')
df_enrich<-data.frame()
for(i in files){
  tmp_files<-list.files(paste0('oc/gwas/enrichresult/',i))
  tmp_files<-tmp_files[grep('results',tmp_files)]
  for(j in tmp_files){
    tmp<-read.delim(paste0('oc/gwas/enrichresult/',i,'/',j))
    padj<-p.adjust(tmp$Enrichment_p,method = 'fdr')
    info<-strsplit(j,'\\.')[[1]]
    tmp_df<-data.frame(cell=info[1],disease=info[2],enrichment=tmp[1,5],Pvalue=tmp[1,7])
    df_enrich<-rbind(df_enrich,tmp_df)}
}
#plot
library(dplyr)
library(ggplot2)
plotdata<-df_enrich %>% mutate(text = case_when(
  Pvalue<0.0001 ~ "****",Pvalue<0.001 ~ "***",Pvalue<0.01 ~"**",Pvalue<0.05 ~"*"))
plotdata$enrichment<-ifelse(plotdata$enrichment<(-50),-50,plotdata$enrichment)
plotdata$enrichment<-ifelse(plotdata$enrichment>50,50,plotdata$enrichment)
p2<-ggplot(plotdata, aes(cell, disease)) + 
  geom_tile(aes(fill = Pvalue), colour = "black", size = 0.2)+
  scale_fill_gradient2(low = '#1A5592',mid = "white",high = "#B83D3D") + 
  geom_text(aes(label=text),col ="black",size = 5) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(size = 10,color = 'black',angle = 90, hjust = 1), 
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.position = 'top') + 
  scale_y_discrete(position = "right")
ggsave('f7/oc_ldsc.pdf',p2,height = 6,width = 8)

###################
###links cc
###################
#load data
cc<-readRDS('../cc/data/scatac/cc.rds')
cc$celltype<-factor(cc$celltype,levels = sort(unique(cc$celltype)))
Idents(cc)<-'celltype'
x<-c('chr16-68778175-68779190','chr13-110422090-110422624')
#links
conns<-readRDS('../cc/cicero/conns_600.rds')
conns<-subset(conns,abs(coaccess)>0.2)
peakanno<-readRDS('../cc/chipseeker/peakanno.rds')
peakanno<-subset(peakanno,type=='Promoter')
conns<-subset(conns,Peak2 %in% rownames(peakanno))
test.df <- data.frame(peakSite = conns$Peak2, 
                      peak2gene = peakanno[as.character(conns$Peak2),'SYMBOL'])
new_conns <- cbind(conns, test.df)
new_conns<-subset(new_conns,Peak1 %in% x)
links <- ConnectionsToLinks(conns = new_conns)
Links(cc) <- links
#plot
gr_cdh1<-StringToGRanges('chr16-68778175-68779190')
gr_cdh1$color<-'yellow'
p1<-CoveragePlot(cc,'chr16-68778175-68779190',extend.upstream = 100000,extend.downstream = 100000,region.highlight =gr_cdh1 )
ggsave('f7/CDH1_links.pdf',p1,height = 6,width = 8)
gr_col4a2<-StringToGRanges('chr13-110422090-110422624')
gr_col4a2$color<-'yellow'
p2<-CoveragePlot(cc,'chr13-110422090-110422624',extend.upstream = 150000,extend.downstream = 50000,region.highlight =gr_col4a2 )
ggsave('f7/COL4A2_links.pdf',p2,height = 6,width = 8)

#manhattan plot
colon_snp<-read.delim('cc/gwas/gwasdata/2020nc_cc_GCST90011811_buildGRCh37.tsv')
ref_38<-read.delim('.cc/gwas/gwasdata/liftover/2020nc_hg19Tohg38.bed')
colon_snp<-colon_snp[ref_38[,4],]
colon_snp$base_pair_location<-ref_38[,2]
#chr16-68778175-68779190  
#rs9939049
plotdata<-subset(colon_snp,chromosome==16 & base_pair_location>=68678175 & base_pair_location<=68879190)
plotdf<-data.frame(matrix(NA,nrow = 201016,ncol = 2))
plotdf$X2<-1
plotdf$X1<-68678175:68879190
for(i in 1:nrow(plotdf)){
  if(plotdf[i,1] %in% plotdata$base_pair_location){
    tmp<-subset(plotdata,base_pair_location==plotdf[i,1])
    plotdf[i,2]<-tmp[1,7]
    next
  }
}
plotdf[,'logPvalue']<-(-log(plotdf$X2))
plotdf$label = ifelse(plotdf$X1==68778398, 'rs9939049',"")
plotdf$shape<-'a'
plotdf$shape<-ifelse(plotdf$X1==68778398,'b',plotdf$shape)
p<-ggplot(plotdf,aes(X1,logPvalue,color=logPvalue))+geom_point(shape=20)+
  scale_color_gradient(low="green",high = "#990000")+
  theme_test()
ggsave('f7/cdh1_man.pdf',p,width = 18,height = 3)

#chr13-110422090-110422624  
#rs7993934
plotdata<-subset(colon_snp,chromosome==13 & base_pair_location>=110272090 & base_pair_location<=110472624)
plotdf<-data.frame(matrix(NA,nrow = 200535,ncol = 2))
plotdf$X2<-1
plotdf$X1<-110272090:110472624
for(i in 1:nrow(plotdf)){
  if(plotdf[i,1] %in% plotdata$base_pair_location){
    tmp<-subset(plotdata,base_pair_location==plotdf[i,1])
    plotdf[i,2]<-tmp[1,7]
    next
  }
}
plotdf[,'logPvalue']<-(-log(plotdf$X2))
plotdf$label = ifelse(plotdf$X1==110422568, 'rs7993934',"")
plotdf$shape<-'a'
plotdf$shape<-ifelse(plotdf$X1==110422568,'b',plotdf$shape)
p<-ggplot(plotdf,aes(X1,logPvalue,color=logPvalue))+geom_point(shape=20)+
  scale_color_gradient(low="green",high = "#990000")+
  theme_test()
p+geom_text_repel(data = plotdf, aes(x = X1,y = logPvalue, label = label),
                  size = 4,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",
                  show.legend = FALSE,max.overlaps = 1000)
ggsave('f7/col4a2_man.pdf',p,width = 18,height = 3)

###################
###links oc
###################
#load data
oc<-readRDS('oc/data/scatac/oc.rds')
DefaultAssay(oc)<-'peaks'
oc$celltype<-factor(oc$celltype,levels = sort(unique(oc$celltype)))
Idents(oc)<-'celltype'
x<-c('chr2-176177677-176179143','chr11-10778352-10778942')
#links
conns<-readRDS('../oc/cicero/conns_600.rds')
conns<-subset(conns,abs(coaccess)>0.2)
peakanno<-readRDS('../oc/chipseeker/peakanno.rds')
peakanno<-subset(peakanno,type=='Promoter')
conns<-subset(conns,Peak2 %in% rownames(peakanno))
test.df <- data.frame(peakSite = conns$Peak2, 
                      peak2gene = peakanno[as.character(conns$Peak2),'SYMBOL'])
new_conns <- cbind(conns, test.df)
new_conns<-subset(new_conns,Peak1 %in% x)
links <- ConnectionsToLinks(conns = new_conns)
Links(oc) <- links
#plot
gr_haglr<-StringToGRanges('chr2-176177677-176179143')
gr_haglr$color<-'yellow'
CoveragePlot(oc,'chr2-176177677-176179143',extend.upstream = 50000,extend.downstream = 150000,region.highlight =gr_haglr )
ggsave('f7/HAGLR_links.pdf',p1,height = 6,width = 8)
gr_ctr9<-StringToGRanges('chr11-10778352-10778942')
gr_ctr9$color<-'yellow'
CoveragePlot(oc,'chr11-10778352-10778942',extend.upstream = 100000,extend.downstream = 100000,region.highlight =gr_ctr9 )
ggsave('f7/CTR9_links.pdf',p1,height = 6,width = 8)

#manhattan plot
oc_snp<-read.delim('../oc/gwas/gwasdata/GCST90041902_20001_1039.v1.0.fastGWA.tsv')
ref_38<-read.delim('../oc/gwas/gwasdata/liftover/2020nc_hg19Tohg38.bed')
oc_snp<-oc_snp[ref_38[,4],]
oc_snp$base_pair_location<-ref_38[,3]
#chr11-10778352-10778942 
#rs141131642
plotdata<-subset(oc_snp,chromosome==11 & base_pair_location>=10678352 & base_pair_location<=10878942)
plotdf<-data.frame(matrix(NA,nrow = 200591,ncol = 2))
plotdf$X2<-1
plotdf$X1<-10678352:10878942
for(i in 1:nrow(plotdf)){
  if(plotdf[i,1] %in% plotdata$base_pair_location){
    tmp<-subset(plotdata,base_pair_location==plotdf[i,1])
    plotdf[i,2]<-tmp[1,13]
    next
  }
}
plotdf[,'logPvalue']<-(-log(plotdf$X2))
plotdf$label = ifelse(plotdf$X1==10778733, 'rs141131642',"")
plotdf$shape<-'a'
plotdf$shape<-ifelse(plotdf$X1==10778733,'b',plotdf$shape)
p<-ggplot(plotdf,aes(X1,logPvalue,color=logPvalue))+geom_point(shape=20)+
  scale_color_gradient(low="green",high = "#990000")+
  theme_test()
p<-p+geom_text_repel(data = plotdf, aes(x = X1,y = logPvalue, label = label),
                     size = 4,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"),
                     segment.color = "black",
                     show.legend = FALSE,max.overlaps = 1000)
ggsave('f7/ctr9_man.pdf',p,width = 18,height = 3)

#chr2-176177677-176179143
#rs6755777
plotdata<-subset(oc_snp,chromosome==2 & base_pair_location>=176077677 & base_pair_location<=176279143)
plotdf<-data.frame(matrix(NA,nrow = 201467,ncol = 2))
plotdf$X2<-1
plotdf$X1<-176077677:176279143
for(i in 1:nrow(plotdf)){
  if(plotdf[i,1] %in% plotdata$base_pair_location){
    tmp<-subset(plotdata,base_pair_location==plotdf[i,1])
    plotdf[i,2]<-tmp[1,13]
    next
  }
}
plotdf[,'logPvalue']<-(-log(plotdf$X2))
plotdf$label = ifelse(plotdf$X1==176178499, 'rs6755777',"")
plotdf$shape<-'a'
plotdf$shape<-ifelse(plotdf$X1==176178499,'b',plotdf$shape)
p<-ggplot(plotdf,aes(X1,logPvalue,color=logPvalue))+geom_point(shape=20)+
  scale_color_gradient(low="green",high = "#990000")+
  theme_test()
p<-p+geom_text_repel(data = plotdf, aes(x = X1,y = logPvalue, label = label),
                  size = 4,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",
                  show.legend = FALSE,max.overlaps = 1000)
ggsave('f7/haglr_man.pdf',p,width = 18,height = 3)