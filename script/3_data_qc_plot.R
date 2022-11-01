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

output=matrix(nrow=4, ncol=5)
rownames(output)=c('bmmc','pbmc','brain','human')
colnames(output)=c('# Cells',	'# peaks',	'Median # reads per cell',
                   'Median frac cells per peak',	'# Cell types')
                   

for(data in c('bmmc','pbmc','brain','human')){
  cat(data,'\n')
  
  load(paste0(data,'_processed_analyzed.rda'))
  
  temp=obj@assays$ATAC@counts
  temp=temp>0 # Binarize the count matrix
  
  # Number of cells
  n.cell=ncol(temp)
  
  # Number of peaks
  n.peak=nrow(temp)
  
  # ATAC reads per cell
  median.reads.per.cell=round(median(colSums(temp)))
  
  # Cells per peak
  median.cells.per.peak=round(median(rowSums(temp)/ncol(temp)), 3)
  
  # Number of cell types
  n.celltype=length(unique(obj$predicted.id))
  
  output[data,]=c(n.cell, n.peak, median.reads.per.cell, median.cells.per.peak, n.celltype)
  
  gc()
}

output
write.csv(output, file='output.csv')



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

library(patchwork)

# Generate dimension reduction
for(data in c('bmmc','pbmc','brain','human')){
  cat(data,'\n')
  
  load(paste0(data,'_processed_analyzed.rda'))
  
  p1=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'lsi.umap') + NoLegend()+ggtitle('Peak accessibility LSI')
  p2=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'lda.umap') + NoLegend()+ggtitle('Peak accessibility LDA')
  p3=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'motif.umap') + NoLegend()+ggtitle('Motif score')
  p4=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'activity.umap') + NoLegend()+ggtitle('Gene activity')
  p5=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'consensus.umap') + NoLegend()+ggtitle('Consensus PCA')
  p6=DimPlot(object = obj, group.by='predicted.id',label = TRUE, reduction = 'multicca.umap') + NoLegend() +ggtitle('MultiCCA')
  p7=DimPlot(object = obj, group.by='predicted.id',label = TRUE, reduction = 'wnn.umap') + NoLegend() +ggtitle('WNN')
  p8=DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = 'lsi.umap') +ggtitle('Place Holder for Legend')
  p=p1+p2+p3+p4+p5+p6+p7+p8+ plot_layout(ncol = 4)
  ggsave(filename=paste0('pdf/umap_',data,'.pdf'), plot = p, width=20, height=10)
}

# Generate heatmap of confusion matrix
for(data in c('bmmc','pbmc','brain','human')){
  cat(data,'\n')
  
  load(paste0(data,'_processed_analyzed.rda'))
  
  p1=get.accuracy(pred.cluster = obj$signac_clusters, real.cluster = obj$predicted.id,
                  title='Peak accessibility LSI', reduction.mat=obj@reductions$lsi.umap@cell.embeddings)$p
  p2=get.accuracy(pred.cluster = obj$cistopic_clusters, real.cluster = obj$predicted.id,
                  title='Peak accessibility LDA', reduction.mat=obj@reductions$lda.umap@cell.embeddings)$p
  p3=get.accuracy(pred.cluster = obj$motif_clusters, real.cluster = obj$predicted.id,
                  title='Motif score', reduction.mat=obj@reductions$motif.umap@cell.embeddings)$p
  p4=get.accuracy(pred.cluster = obj$activity_clusters, real.cluster = obj$predicted.id,
                  title='Gene activity', reduction.mat=obj@reductions$activity.umap@cell.embeddings)$p
  p5=get.accuracy(pred.cluster = obj$consensus_clusters, real.cluster = obj$predicted.id,
                  title='Consensus PCA', reduction.mat=obj@reductions$consensus.umap@cell.embeddings)$p
  p6=get.accuracy(pred.cluster = obj$cc_clusters, real.cluster = obj$predicted.id,
                  title='MultiCCA', reduction.mat=obj@reductions$multicca.umap@cell.embeddings)$p
  p7=get.accuracy(pred.cluster = obj$wnn_clusters, real.cluster = obj$predicted.id,
                  title='WNN', reduction.mat=obj@reductions$wnn.umap@cell.embeddings)$p

  p=p1+p2+p3+p4+p5+p6+p7+ plot_layout(ncol = 4)
  if(data=='pbmc'){
    ggsave(filename=paste0('pdf/heatmap_confusion_',data,'.pdf'), plot = p, width=28, height=10)
  } else{
    ggsave(filename=paste0('pdf/heatmap_confusion_',data,'.pdf'), plot = p, width=25, height=10)
  }

}

# Generate clustree output
for(data in c('bmmc','pbmc','brain','human')){
  cat(data,'\n')
  load(paste0(data,'_processed_analyzed.rda'))
  
  resolution.range <- seq(from = 0.2, to = 1.2, by = 0.2)
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, resolution = resolution.range, graph.name='consensuspca_snn')
  obj.cpca.clustree=clustree(obj, prefix = 'consensuspca_snn_res.')
  p1=obj.cpca.clustree
  # The stability index from the SC3 package (Kiselev et al. 2017) measures the stability of clusters across resolutions.
  # The stability index is automatically calculated when a clustering tree is built.
  obj.cpca.clustree=clustree(obj, prefix = 'consensuspca_snn_res.', node_colour = "sc3_stability")
  p2=obj.cpca.clustree
  print(p1+p2)
  p=p1+p2
  ggsave(filename=paste0('pdf/clustree_cpca_',data,'.pdf'), plot = p, width=12, height=6)
}


# Generate slingshot output
for(data in c('bmmc','pbmc','brain','human')){
  cat(data,'\n')
  load(paste0(data,'_processed_analyzed.rda'))
  
  pdf(file=paste0('pdf/slingshot_cpca_',data,'.pdf'), width=10, height=3.5)
  SlingshotLineages <- GetSlingshot(obj, reduction.to.construct = 'consensus.pca',
                                    reduction.to.plot='consensus.umap',
                                    cluster='consensus_clusters', 
                                    predicted.id='predicted.id',
                                    generate.plot = TRUE)
  dev.off()
}
