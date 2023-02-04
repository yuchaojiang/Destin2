setwd("C:/Users/yuchaoj/Dropbox/ATAC_multi/data_script/")
# setwd("~/Dropbox/ATAC_multi/data_script/")

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

for (seed in 1:20){
  cat('Seed',seed,':\t')
  
  for(data in c('bmmc','pbmc','brain','human')){
    cat(data,'...\n\n')
    
    load(paste0(data,'_processed.rda'))
    eval(parse(text=paste('obj=',data,';','rm(',data,')')))
    
    #########################################################
    ##  Peak Accessibility: LSI
    #########################################################
    
    set.seed(1234)
    DefaultAssay(obj)='ATAC'
    obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:31, verbose = FALSE, graph.name=c('lsi_nn','lsi_snn'))
    obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'lsi_snn', random.seed = seed)
    
    obj$signac_clusters=obj$seurat_clusters

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
    obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'lda_snn', random.seed = seed)
    
    obj$cistopic_clusters=obj$seurat_clusters
    
    temp=get.accuracy(pred.cluster = obj$cistopic_clusters, real.cluster = obj$predicted.id,
                      title='cisTopic', reduction.mat=obj@reductions$lda.umap@cell.embeddings)
    metric=temp$metric
    output=rbind(output, Peak_LDA=metric)
    
    #########################################################
    ##  Motif
    #########################################################
    
    DefaultAssay(obj)='MOTIF'
    obj <- FindNeighbors(object = obj, reduction = 'motif.pca', dims = 1:30, verbose = FALSE, graph.name=c('motif_nn','motif_snn'))
    obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'motif_snn', random.seed = seed)
    
    obj$motif_clusters = obj$seurat_clusters

    temp=get.accuracy(pred.cluster = obj$motif_clusters, real.cluster = obj$predicted.id,
                      title='Motif deviation', reduction.mat=obj@reductions$motif.umap@cell.embeddings)
    metric=temp$metric
    output=rbind(output, Motif=metric)
    
    #########################################################
    ##  Gene Activity
    #########################################################
    
    DefaultAssay(obj)='ACTIVITY'
    obj <- FindNeighbors(object = obj, reduction = 'activity.pca', dims = 1:30, verbose = FALSE, graph.name=c('activity_nn','activity_snn'))
    obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'activity_snn', random.seed = seed)
    
    obj$activity_clusters = obj$seurat_clusters

    temp=get.accuracy(pred.cluster = obj$activity_clusters, real.cluster = obj$predicted.id,
                      title='Gene activity', reduction.mat=obj@reductions$activity.umap@cell.embeddings)
    metric=temp$metric
    output=rbind(output, GeneActivity=metric)
    
    #########################################################
    ##  Consensus PCA
    #########################################################
    
    obj <- GetConsensusPCA(obj, reduction.list=list('lsi','lda','activity.pca','motif.pca'),
                           dims.list = list(2:31, 1:30, 1:30, 1:30),
                           reduction.name = 'consensus.pca',
                           reduction.key = "consensuspca_")
    
    obj <- RunUMAP(object = obj, reduction = 'consensus.pca', dims = 1:30, reduction.name = 'consensus.umap', 
                   reduction.key='consensusumap_', verbose = FALSE)
    
    DefaultAssay(obj)='ATAC'
    obj <- FindNeighbors(object = obj, reduction = 'consensus.pca', dims = 1:30, verbose = FALSE, graph.name=c('consensuspca_nn','consensuspca_snn'))
    obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name='consensuspca_snn', random.seed = seed)
    
    obj$consensus_clusters=obj$seurat_clusters

    temp=get.accuracy(pred.cluster = obj$consensus_clusters, real.cluster = obj$predicted.id,
                      title='Consensus PCA', reduction.mat=obj@reductions$consensus.umap@cell.embeddings)
    metric=temp$metric
    output=rbind(output, ConsensusPCA=metric)
    
    #########################################################
    ##  Multiple CCA
    #########################################################
    
    obj <- GetMultiCCA(obj, reduction.list=list('lsi', 'lda','activity.pca','motif.pca'),
                       dims.list = list(2:31, 1:30, 1:30, 1:30),
                       reduction.name = 'multicca.pca',
                       reduction.key = "multiccapca_")
    
    obj <- RunUMAP(obj, dims = 1:30, reduction='multicca.pca', 
                   reduction.name='multicca.umap', reduction.key = 'multiccaumap_', verbose = FALSE)
    
    obj <- FindNeighbors(object = obj, reduction = 'multicca.pca', dims = 1:30, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
    obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name='multicca_snn', random.seed = seed)
    
    obj$cc_clusters = obj$seurat_clusters

    temp=get.accuracy(pred.cluster = obj$cc_clusters, real.cluster = obj$predicted.id,
                      title='MultiCCA', reduction.mat=obj@reductions$multicca.umap@cell.embeddings)
    metric=temp$metric
    output=rbind(output, MultiCCA=metric)
    
    #########################################################
    ##  Weighted Nearest Neighbor
    #########################################################
    
    DefaultAssay(obj)='ATAC'
    obj <- FindMultiModalNeighbors(obj,
                                   reduction.list=list('lsi','lda','activity.pca','motif.pca'),
                                   dims.list = list(2:31, 1:30, 1:30, 1:30), verbose=FALSE)
    obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                   reduction.key = "wnnumap_", verbose = FALSE) # Run UMAP based on the knn distance
    obj <- FindClusters(object = obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE, random.seed = seed)
    obj$wnn_clusters=obj$seurat_clusters
    
    temp=get.accuracy(pred.cluster = obj$wnn_clusters, real.cluster = obj$predicted.id,
                      title='WNN', reduction.mat=obj@reductions$wnn.umap@cell.embeddings)
    metric=temp$metric
    output=rbind(output, WNN=metric)
    
    print(output)
    cat('\n\n\n')
    save(output, file=paste0('random_seeds/metrics_output_',data,'_',seed,'.rda'))
  }
}



setwd("~/Dropbox/ATAC_multi/data_script/")

results.all=vector(mode = "list", length = 20)

for(seed in 1:20){
  data='pbmc'
  load(paste0('random_seeds/metrics_output_',data,'_',seed,'.rda'))
  results=output
  
  data='brain'
  load(paste0('random_seeds/metrics_output_',data,'_',seed,'.rda'))
  results=cbind(results, output)
  
  data='bmmc'
  load(paste0('random_seeds/metrics_output_',data,'_',seed,'.rda'))
  results=cbind(results, output)
  
  data='human'
  load(paste0('random_seeds/metrics_output_',data,'_',seed,'.rda'))
  results=cbind(results, output)

  results.all[[seed]]=results
}

median=results.all[[1]]
median[1:nrow(median),1:ncol(median)]=NA
lower=upper=output=median

for(i in 1:nrow(median)){
  for(j in 1:ncol(median)){
    data=c()
    for(seed in 1:20){
      data=c(data, results.all[[seed]][i,j])
    }
    median[i,j]=format(round(median(data),3), nsmall=3)
    lower[i,j]=format(round(quantile(data,0.05),3), nsmall=3)
    upper[i,j]=format(round(quantile(data,0.95),3), nsmall=3)
    output[i,j]=paste0(median[i,j], '; [',lower[i,j],', ', upper[i,j],']')
  }
}

write.csv(output, file='output_CI.csv')
