GetConsensusPCA=function(obj, reduction.list,
                         dims.list,
                         reduction.name=NULL,
                         reduction.key=NULL,
                         n.consensus.pc=NULL){
  if(is.null(reduction.name))     reduction.name = 'consensus.pca'
  if(is.null(reduction.key)) reduction.key = 'consensuspca_'
  if(is.null(n.consensus.pc)) n.consensus.pc=30
    
  K=length(reduction.list) # Number of modalities
  mod1.pc=obj@reductions[[reduction.list[[1]]]]@cell.embeddings[,dims.list[[1]]]
  pc.all=mod1.pc
  for(k in 2:K){
    modk.pc=obj@reductions[[reduction.list[[k]]]]@cell.embeddings[,dims.list[[k]]]
    pc.all=cbind(pc.all, modk.pc)
  }
  
  pc.all=apply(pc.all, 2, scale)
  rownames(pc.all)=colnames(obj)
  consensus.pca=pc.all%*%(svd(pc.all)$v[,1:n.consensus.pc])
  colnames(consensus.pca)=paste0('consensuspca_',1:n.consensus.pc)
  
  obj[["consensus.pca"]] <- CreateDimReducObject(embeddings = consensus.pca, 
                                                 key = "consensuspca_", assay = "ATAC")
  
  return(obj)
}


GetMultiCCA=function(obj, reduction.list,
                         dims.list,
                         reduction.name=NULL,
                         reduction.key=NULL,
                         n.cca=NULL){
  if(is.null(reduction.name))     reduction.name = 'multicca.pca'
  if(is.null(reduction.key)) reduction.key = 'multiccapca_'
  if(is.null(n.cca)) n.cca=30
  
  K=length(reduction.list) # Number of modalities
  mod1.pc=obj@reductions[[reduction.list[[1]]]]@cell.embeddings[,dims.list[[1]]]
  pc.all.list=list(t(mod1.pc))
  for(k in 2:K){
    modk.pc=obj@reductions[[reduction.list[[k]]]]@cell.embeddings[,dims.list[[k]]]
    pc.all.list[[k]]=t(modk.pc)
  }
  cc_pca <- mbpca(x = pc.all.list, 
                method = "blockScore", ncomp = n.cca, verbose = FALSE, 
                moa = FALSE, scale = TRUE)$t
  rownames(cc_pca)=colnames(obj)
  colnames(cc_pca)=paste0('ccpca_',1:ncol(cc_pca))
  
  obj[["multicca.pca"]] <- CreateDimReducObject(embeddings = cc_pca, 
                                                key = "multiccapca_", assay = "ATAC")
  return(obj)
}


GetSlingshot=function(obj, reduction.to.construct, reduction.to.plot, cluster,
                      predicted.id=NULL,
                      generate.plot=NULL){
  if(is.null(plot)) plot=TRUE
  set.seed(1)
  dimPlot <- obj@reductions[[reduction.to.plot]]@cell.embeddings
  dimConstruct <- obj@reductions[[reduction.to.construct]]@cell.embeddings
  clustering <-  factor(obj@meta.data[,cluster])
  # We use the PC to infer the lineages
  # But use the UMAP to show the cluster centroids and to
  # show the lineages in the UMAP plot
  lineagesConstruct <- SlingshotDataSet(getLineages(data = dimConstruct, 
                                             clusterLabels = clustering))
  lineagesPlot <- SlingshotDataSet(getLineages(data = dimPlot, 
                                               clusterLabels = clustering))
  lineagesConstruct@reducedDim=lineagesPlot@reducedDim
  SlingshotLineages=lineagesConstruct
  if(generate.plot){
    
    if(!is.null(predicted.id)){
      par(mfrow=c(1,3))
      celltypes <- factor(obj$predicted.id)
      p <-   DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = reduction.to.plot)
      pbuild <- ggplot_build(p)$data[[1]]
      
      plot(dimPlot[, 1:2], col = pbuild$colour, cex = 0.5, pch = 16,
           xlab=paste(reduction.to.plot, '1'), ylab=paste(reduction.to.plot,'2'))
      for (i in levels(celltypes)) {
        text(mean(dimPlot[celltypes == i, 1]), mean(dimPlot[celltypes == i, 2]), labels = i, font = 2)
      }
      title('Predicted cell types')
    } else{
      par(mfrow=c(1,2))
    }
    
    p <-   DimPlot(object = obj, group.by=cluster, label = TRUE, reduction = reduction.to.plot) 
    pbuild <- ggplot_build(p)$data[[1]]
    plot(dimPlot[, 1:2], col = pbuild$colour, cex = 0.5, pch = 16,
         xlab=paste(reduction.to.plot, '1'), ylab=paste(reduction.to.plot,'2'))
    for (i in levels(clustering)) {
      text(mean(dimPlot[clustering == i, 1]), mean(dimPlot[clustering == i, 2]), labels = i, font = 2)
    }
    title('Cell clusters')
    
    plot(dimPlot[, 1:2], col = pbuild$colour, cex = 0.5, pch = 16,
         xlab=paste(reduction.to.plot, '1'), ylab=paste(reduction.to.plot,'2'))
    lines(SlingshotLineages, lwd = 1, col = "black", cex=1)
    title('Slingshot trajectory')
    par(mfrow=c(1,1))
  }
  return(SlingshotLineages)
}


get.accuracy=function(pred.cluster, real.cluster, x.lab=NULL, y.lab=NULL, title=NULL, reduction.mat){
  if(is.null(x.lab)) x.lab='ATAC clusters'
  if(is.null(y.lab)) y.lab='Cell types transferred from RNA'
  if(is.null(title)) title=''
  AccuracyTable=table(pred.cluster, real.cluster)
  pred.label= colnames(AccuracyTable)[apply(AccuracyTable,1, which.max)]
  ari = signif(ARI(pred.cluster, real.cluster),4)
  ami = signif(AMI(pred.cluster, real.cluster),4)
  hscore = signif(homogeneity(real.cluster, pred.cluster),4)
  
  if(!any(is.na(reduction.mat))){
    meta.data=data.frame(pred.cluster, real.cluster)
    clisi = as.matrix(compute_lisi(reduction.mat, meta.data, 'real.cluster'))
    clisi = signif(mean(clisi),4)
    
    # Order the table to be plotted in the heatmap
    AccuracyTable=AccuracyTable[order(pred.label),]
    AccuracyTable=AccuracyTable[,order(colnames(AccuracyTable))]
    predictions <- AccuracyTable/rowSums(AccuracyTable) # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)
    p <- ggplot(predictions, aes(pred.cluster, real.cluster, fill = Freq)) + 
      geom_tile() + scale_fill_gradient(name = "Fraction of cells", low = "#ffffc8", high = "#7d0025") + 
      xlab(x.lab) + ylab(y.lab) +
      theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      ggtitle(paste0(title, ':\nARI = ', ari,'; AMI = ', ami, '; H-score = ', hscore, '; cLISI = ',clisi))+
      theme(plot.title = element_text(size = 10, face = "bold"))
    print(p)
    metric=c(ari, ami, hscore, clisi)
    names(metric)=c('ARI', 'AMI', 'Homogeneity', 'cLISI')
  } else{
    metric=c(ari, ami, hscore, NA, NA)
    names(metric)=c('ARI', 'AMI', 'Homogeneity', 'cLISI')
  }
  return(list(AccuracyTable=AccuracyTable, metric=metric, p=p))
}

