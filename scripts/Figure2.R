
setwd('data/pancancer_atac/tumor/combined/')
set.seed(1001)

#load package
library(Seurat)
library(Signac)
library(ggplot2)
library(tidydr)
library(dplyr)

#######################
###umap tumorcell
########################
#load data
tumor<-readRDS('cells/tumor/tumor.rds')
p<-DimPlot(tumor,group.by = 'tumorcell')
ggsave(p,'f2/TumorUmap.pdf',height = 6,width = 8)

################
###peak heatmap
################
#load data
tumor<-readRDS('cells/tumor/tumor.rds')
peaks_tumor<-readRDS('cells/tumor/peaks_tumor.rds')
peaks_tumor<-subset(peaks_tumor,p_val_adj<0.05)
meta<-tumor@meta.data
#plot
library(pheatmap)
peakdata<-tumor@assays$peaks@data
peakdata<-peakdata[unique(peaks_tumor$gene),]
plotdata<-data.frame(matrix(0,nrow = nrow(peakdata),ncol = length(unique(peaks_tumor$cluster))))
rownames(plotdata)<-rownames(peakdata)
colnames(plotdata)<-unique(peaks_tumor$cluster)
for(i in colnames(plotdata)){ 
  tmp<-subset(meta,tumorcell==i)
  plotdata[,i]<-rowMeans(peakdata[,rownames(tmp)])
}
p<-pheatmap(t(scale(t(plotdata))),show_rownames = F,cluster_rows = F,cluster_cols = T)
ggsave('f2/peak_heatmap.pdf',p,width = 4,height = 4)

###################
###coverageplot
###################
#load data
peaks<-readRDS('cells/tumor/peaks_tumor.rds')
peakanno<-readRDS('chipseeker/peakanno.rds')
peaks$closegene<-peakanno[peaks$gene,'SYMBOL']
peaks$type<-peakanno[peaks$gene,'type']
peaks<-subset(peaks,p_val_adj<0.05)
#plot
gr<-StringToGRanges(c('chr8-128531703-128542048','chr8-127875000-127878526','chr8-127389527-127458527'))
gr$color<-'yellow'
p<-CoveragePlot(tumor,'MYC',extend.upstream = 500000,extend.downstream = 900000,
                region.highlight = gr)
ggsave('f2/myc.pdf',p,width = 6,height = 4)
gr<-StringToGRanges(c('chr8-133255001-133263233','chr8-133244401-133251916','chr8-133301001-133304600','chr8-133214802-133219600','chr8-133320000-133321000','chr8-133287484-133300561'))
gr$color<-'yellow'
p<-CoveragePlot(tumor,'NDRG1',extend.upstream = 40000,extend.downstream = 40000,region.highlight = gr)
ggsave('f2/ndrg1.pdf',p,width = 6,height = 4)

#####################
###ndrg1 link plot
#####################
#load data
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
#cicero
conns<-readRDS('cicero/conns_600.rds')
conns<-subset(conns,abs(coaccess)>0.2)
peakanno<-readRDS('chipseeker/peakanno.rds')
peakanno<-subset(peakanno,type=='Promoter')
conns<-subset(conns,Peak2 %in% rownames(peakanno))
test.df <- data.frame(peakSite = conns$Peak2, 
                      peak2gene = peakanno[as.character(conns$Peak2),'SYMBOL'])
new_conns <- cbind(conns, test.df)
chia_conns <-  data.frame(Peak1 = c("chr8-133320367-133321081"), 
                          Peak2 = c("chr8-133298299-133298499"),
                          coaccess = c(0.3203461))
#plot
p<-plot_connections(chia_conns,'chr8',133275000,133325000,gene_model = gene_anno,
                 collapseTranscripts = "longest",include_axis_track = T,
                 gene_model_color = 'black',
                 coaccess_cutoff = 0)
ggsave('f2/cancer.pdf',p,width = 6,height = 6)