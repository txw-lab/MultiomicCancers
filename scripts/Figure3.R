
setwd('tumor')
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

###################
#coverageplot MYC
###################
tumor<-readRDS('combined/cells/tumor/tumor.rds')
tumor$tumorcell<-factor(tumor$tumorcell,levels = sort(unique(tumor$tumorcell)))
Idents(tumor)<-'tumorcell'
peaks<-readRDS('combined/cells/tumor/peaks_tumor.rds')
peakanno<-readRDS('chipseeker/peakanno.rds')
peaks$closegene<-peakanno[peaks$gene,'SYMBOL']
peaks$type<-peakanno[peaks$gene,'type']
peaks<-subset(peaks,p_val_adj<0.05)
gr<-StringToGRanges(c('chr8-128531703-128542048','chr8-127870000-127883526',
                      'chr8-127387527-127460527','chr8-128154301-128178044'))
gr$color<-'yellow'
p<-CoveragePlot(tumor,'MYC',extend.upstream = 350000,extend.downstream = 800000,
                region.highlight = gr)
ggsave('f3/MYC.pdf',p,width = 6,height = 5)

gr1<-StringToGRanges(c('chr8-128531703-128542048'))
gr1$color<-'yellow'
p1<-CoveragePlot(tumor,'chr8-128531703-128542048',extend.upstream = 20000,extend.downstream = 20000,region.highlight = gr1)
gr2<-StringToGRanges(c('chr8-127870000-127883526'))
gr2$color<-'yellow'
p2<-CoveragePlot(tumor,'chr8-127870000-127883526',extend.upstream = 20000,extend.downstream = 20000,region.highlight = gr2)
gr3<-StringToGRanges(c('chr8-127387527-127460527'))
gr3$color<-'yellow'
p3<-CoveragePlot(tumor,'chr8-127387527-127460527',extend.upstream = 30000,extend.downstream = 30000,region.highlight = gr3)
gr4<-StringToGRanges(c('chr8-128154301-128178044'))
gr4$color<-'yellow'
p4<-CoveragePlot(tumor,'chr8-128154301-128178044',extend.upstream = 30000,extend.downstream = 30000,region.highlight = gr4)
ggsave('f3/MYC1.pdf',p1,width = 2.5,height = 5)
ggsave('f3/MYC2.pdf',p2,width = 2.5,height = 5)
ggsave('f3/MYC3.pdf',p3,width = 2.5,height = 5)
ggsave('f3/MYC4.pdf',p4,width = 2.5,height = 5)

###################
#coverageplot NDRG1
###################
tumor<-readRDS('combined/cells/tumor/tumor.rds')
tumor$tumorcell<-factor(tumor$tumorcell,levels = sort(unique(tumor$tumorcell)))
Idents(tumor)<-'tumorcell'
gr<-StringToGRanges(c('chr8-133255001-133263233','chr8-133244401-133251916','chr8-133302001-133304600','chr8-133214802-133219600','chr8-133320000-133321000','chr8-133287484-133298561'))
gr$color<-'yellow'
p<-CoveragePlot(tumor,'NDRG1',extend.upstream = 40000,extend.downstream = 40000,region.highlight = gr)
ggsave('f3/NDRG1.pdf',p,width = 6,height = 5)

###################
###diff peak
###################
tumor<-readRDS('combined/cells/tumor/tumor.rds')
Idents(tumor)<-'tumorcell'
peaks_tumor<-FindMarkers(tumor,ident.1 = 'CC', test.use = 'LR',
                         latent.vars = 'nCount_peaks',logfc.threshold = 0)
peaks_tumor<-subset(peaks_tumor,avg_log2FC<0 & p_val_adj<0.05)
peaks_tumor<-peaks_tumor[order(peaks_tumor$avg_log2FC,decreasing = F),]
peaks_tumor$rank<-c(1:nrow(peaks_tumor))
plot(x=peaks_tumor$rank,y=-peaks_tumor$avg_log2FC)
peaks_tumor$peak<-NA
peaks_tumor$peak[41]<-rownames(peaks_tumor)[41]
p<-ggplot(peaks_tumor,aes(x=rank,y=-avg_log2FC,))+
  geom_point(shape=1)+
  theme_classic()+
  annotate('text',x=41,y=0.3424377,label='chr8-133320109-133321082')
ggsave('f3/CC_downpeak.pdf',p,height = 4,width = 5.5)

################
##NDRG1 plot
################
gtf_file = 'E:/gene_id_file/Homo_sapiens.GRCh38.104.chr.gtf'
gene_anno <- rtracklayer::readGFF(gtf_file)
# rename some columns to match requirements
gene_anno$chromosome <- paste0('chr',gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
gene_anno<-subset(gene_anno,type=='exon')
gene_anno<-gene_anno[!is.na(gene_anno$symbol),]
gene_anno<-gene_anno[!is.na(gene_anno$transcript),]
gene_anno<-gene_anno[!is.na(gene_anno$gene),]

chia_conns_bc<-data.frame(Peak1 = c("chr8-133320333-133320960"), 
                          Peak2 = c("chr8-133296788-133298035"),
                          coaccess = c(0.03))
chia_conns_bcc<-data.frame(Peak1 = c("chr8-133320333-133320960"), 
                          Peak2 = c("chr8-133296788-133298035"),
                          coaccess = c(0.24))
chia_conns_ec<-data.frame(Peak1 = c("chr8-133320333-133320960"), 
                          Peak2 = c("chr8-133296788-133298035"),
                          coaccess = c(0.32))
chia_conns_lc<-data.frame(Peak1 = c("chr8-133320333-133320960"), 
                          Peak2 = c("chr8-133296788-133298035"),
                          coaccess = c(0.58))
chia_conns_oc<-data.frame(Peak1 = c("chr8-133320333-133320960"), 
                          Peak2 = c("chr8-133296788-133298035"),
                          coaccess = c(0.62))
chia_conns_plc<-data.frame(Peak1 = c("chr8-133320333-133320960"), 
                          Peak2 = c("chr8-133296788-133298035"),
                          coaccess = c(0.26))
chia_conns_rcc<-data.frame(Peak1 = c("chr8-133320333-133320960"), 
                          Peak2 = c("chr8-133296788-133298035"),
                          coaccess = c(0.88))
p_bc<-plot_connections(chia_conns_bc,'chr8',133275000,133325000,gene_model = gene_anno,
                   collapseTranscripts = "longest",include_axis_track = T,
                   gene_model_color = 'black',
                   coaccess_cutoff = 0)
p_bcc<-plot_connections(chia_conns_bc,'chr8',133275000,133325000,gene_model = gene_anno,
                       collapseTranscripts = "longest",include_axis_track = T,
                       gene_model_color = 'black',
                       coaccess_cutoff = 0)
p_ec<-plot_connections(chia_conns_bcc,'chr8',133275000,133325000,gene_model = gene_anno,
                       collapseTranscripts = "longest",include_axis_track = T,
                       gene_model_color = 'black',
                       coaccess_cutoff = 0)
p_lc<-plot_connections(chia_conns_lc,'chr8',133275000,133325000,gene_model = gene_anno,
                       collapseTranscripts = "longest",include_axis_track = T,
                       gene_model_color = 'black',
                       coaccess_cutoff = 0)
p_oc<-plot_connections(chia_conns_oc,'chr8',133275000,133325000,gene_model = gene_anno,
                       collapseTranscripts = "longest",include_axis_track = T,
                       gene_model_color = 'black',
                       coaccess_cutoff = 0)
p_plc<-plot_connections(chia_conns_plc,'chr8',133275000,133325000,gene_model = gene_anno,
                       collapseTranscripts = "longest",include_axis_track = T,
                       gene_model_color = 'black',
                       coaccess_cutoff = 0)
p_rcc<-plot_connections(chia_conns_rcc,'chr8',133275000,133325000,gene_model = gene_anno,
                       collapseTranscripts = "longest",include_axis_track = T,
                       gene_model_color = 'black',
                       coaccess_cutoff = 0)
ggsave('f3/ndrg1_bc.pdf',p_bc,height = 5,width = 5)
ggsave('f3/ndrg1_bcc.pdf',p_bcc,height = 5,width = 5)
ggsave('f3/ndrg1_ec.pdf',p_ec,height = 5,width = 5)
ggsave('f3/ndrg1_lc.pdf',p_lc,height = 5,width = 5)
ggsave('f3/ndrg1_oc.pdf',p_oc,height = 5,width = 5)
ggsave('f3/ndrg1_plc.pdf',p_plc,height = 5,width = 5)
ggsave('f3/ndrg1_rcc.pdf',p_rcc,height = 5,width = 5)
