setwd("C:/Users/yuchaoj/Dropbox/ATAC_multi/data_script/")
setwd("~/Dropbox/ATAC_multi/data_script/")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(cowplot)
library(aricode)
library(PMA)
library(cluster)
library(lisi)
library(FNN)
library(clustree)
library(MAESTRO)
library(slingshot)
library(irlba)
library(cisTopic)
library(umap)
library(clevr)
library(mogsa)
source('0_utility.R')

data='bmmc'
data='pbmc'
data='brain'
data='human'

for(data in c('bmmc','pbmc','brain','human')){
  cat(data,'\n\n\n\n\n')
  
  load(paste0(data,'_processed.rda'))
  eval(parse(text=paste('obj=',data,';','rm(',data,')')))
  
  #########################################################
  ##  Peak Accessibility: LSI
  #########################################################
  
  set.seed(1234)
  DefaultAssay(obj)='ATAC'
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:31, verbose = FALSE, graph.name=c('lsi_nn','lsi_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'lsi_snn')
  
  obj$signac_clusters=obj$seurat_clusters
  p1=DimPlot(object = obj, group.by='signac_clusters', label = TRUE, reduction = 'lsi.umap') + NoLegend()+ggtitle('Peak accessibility LSI')
  p2=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'lsi.umap') + NoLegend()+ggtitle('Peak accessibility LSI')
  p1+p2
  
  temp=get.accuracy(pred.cluster = obj$signac_clusters, real.cluster = obj$predicted.id,
                    title='Signac', reduction.mat=obj@reductions$lsi.umap@cell.embeddings)
  metric=temp$metric
  output=t(as.matrix(metric))
  rownames(output)='Peak_LSI'
  
  #########################################################
  ##  Peak Accessibility: LDA
  #########################################################
  
  DefaultAssay(obj)='ATAC'
  obj <- FindNeighbors(object = obj, reduction = 'lda', dims = 1:30, verbose = FALSE, graph.name=c('lda_nn','lda_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'lda_snn')
  
  obj$cistopic_clusters=obj$seurat_clusters
  p1=DimPlot(object = obj, group.by='cistopic_clusters', label = TRUE, reduction = 'lda.umap') + NoLegend()+ggtitle('Peak accessibility LDA')
  p2=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'lda.umap') + NoLegend()+ggtitle('Peak accessibility LDA')
  p1+p2
  
  temp=get.accuracy(pred.cluster = obj$cistopic_clusters, real.cluster = obj$predicted.id,
                    title='cisTopic', reduction.mat=obj@reductions$lda.umap@cell.embeddings)
  metric=temp$metric
  output=rbind(output, Peak_LDA=metric)
  
  # #########################################################
  # ##  Motif
  # #########################################################
  # 
  # DefaultAssay(obj)='MOTIF'
  # obj <- FindNeighbors(object = obj, reduction = 'motif.pca', dims = 1:30, verbose = FALSE, graph.name=c('motif_nn','motif_snn'))
  # obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'motif_snn')
  # 
  # obj$motif_clusters = obj$seurat_clusters
  # p1=DimPlot(object = obj, group.by='motif_clusters', label = TRUE, reduction = 'motif.umap') + NoLegend()+ggtitle('Motif score')
  # p2=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'motif.umap') + NoLegend()+ggtitle('Motif score')
  # p1+p2
  # 
  # temp=get.accuracy(pred.cluster = obj$motif_clusters, real.cluster = obj$predicted.id,
  #                   title='Motif deviation', reduction.mat=obj@reductions$motif.umap@cell.embeddings)
  # metric=temp$metric
  # output=rbind(output, Motif=metric)
  
  #########################################################
  ##  Gene Activity
  #########################################################
  
  DefaultAssay(obj)='ACTIVITY'
  obj <- FindNeighbors(object = obj, reduction = 'activity.pca', dims = 1:30, verbose = FALSE, graph.name=c('activity_nn','activity_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'activity_snn')
  
  obj$activity_clusters = obj$seurat_clusters
  p1=DimPlot(object = obj, group.by='activity_clusters', label = TRUE, reduction = 'activity.umap') + NoLegend()+ggtitle('Gene activity')
  p2=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'activity.umap') + NoLegend()+ggtitle('Gene activity')
  p1+p2
  
  temp=get.accuracy(pred.cluster = obj$activity_clusters, real.cluster = obj$predicted.id,
                    title='Gene activity', reduction.mat=obj@reductions$activity.umap@cell.embeddings)
  metric=temp$metric
  output=rbind(output, GeneActivity=metric)
  
  #########################################################
  ##  Consensus PCA
  #########################################################
  
  obj <- GetConsensusPCA(obj, reduction.list=list('lsi','lda','activity.pca'),
                         dims.list = list(2:31, 1:30, 1:30),
                         reduction.name = 'consensus.pca',
                         reduction.key = "consensuspca_")
  
  obj <- RunUMAP(object = obj, reduction = 'consensus.pca', dims = 1:30, reduction.name = 'consensus.umap', 
                 reduction.key='consensusumap_', verbose = FALSE)
  
  DefaultAssay(obj)='ATAC'
  obj <- FindNeighbors(object = obj, reduction = 'consensus.pca', dims = 1:30, verbose = FALSE, graph.name=c('consensuspca_nn','consensuspca_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name='consensuspca_snn')
  
  obj$consensus_clusters=obj$seurat_clusters
  p1=DimPlot(object = obj, group.by='consensus_clusters', label = TRUE, reduction = 'consensus.umap') + 
    NoLegend()+ggtitle('Consensus PCA')
  p2=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'consensus.umap') +
    NoLegend()+ggtitle('Consensus PCA')
  p1+p2
  
  temp=get.accuracy(pred.cluster = obj$consensus_clusters, real.cluster = obj$predicted.id,
                    title='Consensus PCA', reduction.mat=obj@reductions$consensus.umap@cell.embeddings)
  metric=temp$metric
  output=rbind(output, ConsensusPCA=metric)
  
  #########################################################
  ##  Multiple CCA
  #########################################################
  
  obj <- GetMultiCCA(obj, reduction.list=list('lsi', 'lda','activity.pca'),
                     dims.list = list(2:31, 1:30, 1:30),
                     reduction.name = 'multicca.pca',
                     reduction.key = "multiccapca_")
  
  obj <- RunUMAP(obj, dims = 1:30, reduction='multicca.pca', 
                 reduction.name='multicca.umap', reduction.key = 'multiccaumap_', verbose = FALSE)
  
  obj <- FindNeighbors(object = obj, reduction = 'multicca.pca', dims = 1:30, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name='multicca_snn')
  
  obj$cc_clusters = obj$seurat_clusters
  p1=DimPlot(object = obj, label = TRUE, group.by='cc_clusters',reduction = 'multicca.umap') + NoLegend()+ggtitle('MultiCCA')
  p2=DimPlot(object = obj, group.by='predicted.id',label = TRUE, reduction = 'multicca.umap') + NoLegend() +ggtitle('MultiCCA')
  p1+p2
  
  temp=get.accuracy(pred.cluster = obj$cc_clusters, real.cluster = obj$predicted.id,
                    title='MultiCCA', reduction.mat=obj@reductions$multicca.umap@cell.embeddings)
  metric=temp$metric
  output=rbind(output, MultiCCA=metric)
  
  #########################################################
  ##  Weighted Nearest Neighbor
  #########################################################
  
  DefaultAssay(obj)='ATAC'
  obj <- FindMultiModalNeighbors(obj,
                                 reduction.list=list('lsi','lda','activity.pca'),
                                 dims.list = list(2:31, 1:30, 1:30), verbose=FALSE)
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                 reduction.key = "wnnumap_", verbose = FALSE) # Run UMAP based on the knn distance
  obj <- FindClusters(object = obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  obj$wnn_clusters=obj$seurat_clusters
  
  p1=DimPlot(object = obj, label = TRUE, group.by='wnn_clusters',reduction = 'wnn.umap') + NoLegend()+ggtitle('WNN')
  p2=DimPlot(object = obj, group.by='predicted.id',label = TRUE, reduction = 'wnn.umap') + NoLegend() +ggtitle('WNN')
  p1+p2
  
  temp=get.accuracy(pred.cluster = obj$wnn_clusters, real.cluster = obj$predicted.id,
                    title='WNN', reduction.mat=obj@reductions$wnn.umap@cell.embeddings)
  metric=temp$metric
  output=rbind(output, WNN=metric)
  
  print(output)
  cat('\n\n\n')
  write.csv(output, file=paste0('metrics/metrics_output_',data,'_no_motif.csv'), row.names = T)
  
  #########################################################
  ##  clustree to determine the number of clusters / resolution
  #########################################################
  
  resolution.range <- seq(from = 0.2, to = 1.2, by = 0.2)
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, resolution = resolution.range, graph.name='consensuspca_snn')
  obj.cpca.clustree=clustree(obj, prefix = 'consensuspca_snn_res.')
  p1=obj.cpca.clustree
  # The stability index from the SC3 package (Kiselev et al. 2017) measures the stability of clusters across resolutions.
  # The stability index is automatically calculated when a clustering tree is built.
  obj.cpca.clustree=clustree(obj, prefix = 'consensuspca_snn_res.', node_colour = "sc3_stability")
  p2=obj.cpca.clustree
  print(p1+p2)
  
  resolution.range <- seq(from = 0.2, to = 1.2, by = 0.2)
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, resolution = resolution.range, graph.name='multicca_snn')
  obj.mcca.clustree=clustree(obj, prefix = 'multicca_snn_res.')
  p1=obj.mcca.clustree
  obj.mcca.clustree=clustree(obj, prefix = 'multicca_snn_res.', node_colour = "sc3_stability")
  p2=obj.mcca.clustree
  p1+p2
  
  resolution.range <- seq(from = 0.2, to = 1.2, by = 0.2)
  obj <- FindClusters(object = obj, algorithm = 1, resolution = resolution.range, verbose = FALSE, graph.name = "wsnn")
  obj.wnn.clustree=clustree(obj, prefix = 'wsnn_res.')
  p1=obj.wnn.clustree
  obj.wnn.clustree=clustree(obj, prefix = 'wsnn_res.', node_colour = "sc3_stability")
  p2=obj.wnn.clustree
  p1+p2
  
  #########################################################
  ##  Slingshot
  #########################################################
  
  SlingshotLineages <- GetSlingshot(obj, reduction.to.construct = 'consensus.pca',
                                    reduction.to.plot='consensus.umap',
                                    cluster='consensus_clusters', 
                                    predicted.id='predicted.id',
                                    generate.plot = TRUE)
  
}


# Compare results with and without motif
setwd("~/Dropbox/ATAC_multi/data_script/")
data='pbmc'
results=read.csv(paste0('metrics/metrics_output_',data,'.csv'), row.names = 1)
results=results[-(1:4),]
results.nomotif=read.csv(paste0('metrics/metrics_output_',data,'_no_motif.csv'), row.names = 1)
results.nomotif=results.nomotif[-(1:3),]
rownames(results.nomotif)=paste0(rownames(results.nomotif),'_no_motif')
results.all.temp=rbind(results, results.nomotif)
results.all=results.all.temp

data='brain'
results=read.csv(paste0('metrics/metrics_output_',data,'.csv'), row.names = 1)
results=results[-(1:4),]
results.nomotif=read.csv(paste0('metrics/metrics_output_',data,'_no_motif.csv'), row.names = 1)
results.nomotif=results.nomotif[-(1:3),]
rownames(results.nomotif)=paste0(rownames(results.nomotif),'_no_motif')
results.all.temp=rbind(results, results.nomotif)
results.all=cbind(results.all, results.all.temp)

data='bmmc'
results=read.csv(paste0('metrics/metrics_output_',data,'.csv'), row.names = 1)
results=results[-(1:4),]
results.nomotif=read.csv(paste0('metrics/metrics_output_',data,'_no_motif.csv'), row.names = 1)
results.nomotif=results.nomotif[-(1:3),]
rownames(results.nomotif)=paste0(rownames(results.nomotif),'_no_motif')
results.all.temp=rbind(results, results.nomotif)
results.all=cbind(results.all, results.all.temp)

data='human'
results=read.csv(paste0('metrics/metrics_output_',data,'.csv'), row.names = 1)
results=results[-(1:4),]
results.nomotif=read.csv(paste0('metrics/metrics_output_',data,'_no_motif.csv'), row.names = 1)
results.nomotif=results.nomotif[-(1:3),]
rownames(results.nomotif)=paste0(rownames(results.nomotif),'_no_motif')
results.all.temp=rbind(results, results.nomotif)
results.all=cbind(results.all, results.all.temp)

# Reorder the rows
results.all=results.all[c(1,4,2,5,3,6),]

pdf(file='pdf/metrics_motif_no_motif.pdf', width=12, height=4)
par(mar=c(11.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,4))
for(i in 1:4){
  if(colnames(results.all)[i] == 'cLISI'){
    ylabc=paste(colnames(results.all)[i], '(smaller = better)')
  } else{
    ylabc=paste(colnames(results.all)[i], '(larger = better)')
  }
  boxplot(t(results.all[,seq(i,ncol(results.all),4)]),
          las=2,
          ylab=ylabc,
          col=rep(c("grey47","grey77"),3),
          main=colnames(results.all)[i])
}
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,1))


