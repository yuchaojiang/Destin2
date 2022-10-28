#' Consensus PCA of Seurat
#'
#' This function takes a Seurat object and conducts consensus principal component analysis (CPCA) as joint dimension reduction
#'
#' @param obj Seurat object
#' @param reduction.list A list of of dimension reduction objects within obj
#' @param dims.list A list of which dimensions to use for each combined dimension reduction object
#' @param assay A character string specifying the assay used to calculate this dimensional reduction
#' @param reduction.name A character string providing the name of the CPCA dimension reduction
#' @param reduction.key A character string that acts as a prefix for the CPCA dimension reduction
#' @param n.consensus.pc Number of principal components (PC's) to store
#'
#' @return A Seurat object containing CPCA joint dimension reduction in slot "consensus.pca"
#'
#' @importFrom SeuratObject CreateDimReducObject
#' @export
GetConsensusPCA=function(obj, reduction.list,
                         dims.list,
                         assay = NULL,
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

  obj[["consensus.pca"]] <- SeuratObject::CreateDimReducObject(embeddings = consensus.pca,
                                                 key = "consensuspca_", assay = assay)

  return(obj)
}

#' MultiCCA of Seurat
#'
#' This function takes a Seurat object and conducts multiple canonical correlation analysis (MCCA) as joint dimension reduction
#'
#' @param obj Seurat object
#' @param reduction.list A list of of dimension reduction objects within obj
#' @param dims.list A list of which dimensions to use for each combined dimension reduction object
#' @param assay A character string specifying the assay used to calculate this dimensional reduction
#' @param reduction.name A character string providing the name of the MCCA dimension reduction
#' @param reduction.key A character string that acts as a prefix for the MCCA dimension reduction
#' @param n.consensus.pc Number of principal components (PC's) to store
#'
#' @return A Seurat object containing MCCA joint dimension reduction in slot "multicca.pca"
#'
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom mogsa mbpca
#' @export
GetMultiCCA=function(obj, reduction.list,
                     dims.list,
                     assay = NULL,
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
  cc_pca <- mogsa::mbpca(x = pc.all.list,
                  method = "blockScore", ncomp = n.cca, verbose = FALSE,
                  moa = FALSE, scale = TRUE)$t
  rownames(cc_pca)=colnames(obj)
  colnames(cc_pca)=paste0('ccpca_',1:ncol(cc_pca))

  obj[["multicca.pca"]] <- SeuratObject::CreateDimReducObject(embeddings = cc_pca,
                                                key = "multiccapca_", assay = assay)
  return(obj)
}

#' Get Slingshot Lineages of Seurat object
#'
#' This function takes a Seurat object and displays trajectory inference using Slingshot
#'
#' @param obj A Seurat object
#' @param reduction.to.construct A character string specifying dimension reduction object used to infer lineages
#' @param reduction.to.plot A character string specifying dimension reduction object to cluster centroids and plot lineages
#' @param cluster A vector specifying cluster identity for each sample
#' @param predicted.id A character string specifying (transferred) cell labels
#' @param generate.plot If true returns Slingshot plots based on cluster and predicted.id
#'
#' @return Slingshot Lineages
#'
#' @importFrom SeuratObject CreateDimReducObject
#' @import Seurat
#' @import slingshot
#' @import ggplot2
#' @export
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
  lineagesConstruct <- slingshot::SlingshotDataSet(slingshot::getLineages(data = unname(dimConstruct),
                                                    clusterLabels = clustering))
  lineagesPlot <- slingshot::SlingshotDataSet(slingshot::getLineages(data = unname(dimPlot),
                                               clusterLabels = clustering))
  lineagesConstruct@reducedDim=lineagesPlot@reducedDim
  SlingshotLineages=lineagesConstruct
  if(generate.plot){

    if(!is.null(predicted.id)){
      par(mfrow=c(1,3))
      celltypes <- factor(obj$predicted.id)
      p <-   Seurat::DimPlot(object = obj, group.by='predicted.id', label = TRUE, reduction = reduction.to.plot)
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

    p <-   Seurat::DimPlot(object = obj, group.by=cluster, label = TRUE, reduction = reduction.to.plot)
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

