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
library(slingshot)
library(irlba)
library(cisTopic)
library(umap)
library(clevr)
library(mogsa)
library(motifmatchr)
library(chromVAR)

source('0_utility.R')

data='bmmc'
data='pbmc'
data='brain'
data='human'

output=matrix(nrow=4, ncol=7)
rownames(output)=c('bmmc','pbmc','brain','human')
colnames(output)=c('lsi','lda','activity','motif','cpca','mcca','wnn')

for(data in rownames(output)){
  cat(data,'\n\n')
  
  load(paste0(data,'_processed.rda'))
  eval(parse(text=paste('obj=',data,';','rm(',data,')')))
  
  if(data=='human'){
    ref_genome='hg19'
    path = 'Human/fragment.txt.gz'
  } else if(data=='bmmc'){
    ref_genome='hg19'
    path = "BMMC/GSM4138889_scATAC_BMMC_D6T1.fragments.tsv.gz"
  } else if(data=='pbmc'){
    ref_genome='hg19'
    path = 'PBMC/atac_v1_pbmc_10k_fragments.tsv.gz'
  } else if (data=='brain'){
    ref_genome='mm10'
    path = 'Adult_mouse_brain/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz'
  }
  
  # Update the path to the fragment files
  frags <- Fragments(obj)  # get list of fragment objects
  Fragments(obj) <- NULL  # remove fragment information from assay
  frags[[1]] <-  UpdatePath(frags[[1]], new.path = path)
  Fragments(obj) <- frags # assign updated list back to the object
  
  if(ref_genome=='hg38'){
    library(EnsDb.Hsapiens.v86) # hg38
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome=BSgenome.Hsapiens.UCSC.hg38
    ensdb=EnsDb.Hsapiens.v86
  } else if (ref_genome=='hg19'){
    library(EnsDb.Hsapiens.v75) # hg19
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome=BSgenome.Hsapiens.UCSC.hg19
    ensdb=EnsDb.Hsapiens.v75
  } else if (ref_genome=='mm10'){
    library(EnsDb.Mmusculus.v79) # mm10
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome=BSgenome.Mmusculus.UCSC.mm10
    ensdb=EnsDb.Mmusculus.v79
  }
  
  
  
  #############################################
  #############################################
  ##
  ##              Peak matrix processing
  ##
  #############################################
  #############################################
  start_time <- Sys.time()
  DefaultAssay(obj)='ATAC'
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = 'q0')
  obj <- RunSVD(obj) # This by default returns an lsi reduction
  
  obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                 reduction.key='lsiumap_', verbose = FALSE)
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
  obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
  end_time <- Sys.time()  
  output[data,'lsi']=round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "mins")),2)
  
  #############################################
  #############################################
  ##
  ##              Gene activity processing
  ##
  #############################################
  #############################################
  
  # To create a gene activity matrix, we extract gene coordinates 
  # and extend them to include the 2 kb upstream region (as promoter 
  # accessibility is often correlated with gene expression). We then 
  # count the number of fragments for each cell that map to each of 
  # these regions, using the using the FeatureMatrix() function. 
  start_time <- Sys.time()
  gene.activities <- GeneActivity(obj)
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  obj[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
  DefaultAssay(obj) <- "ACTIVITY"
  obj <- NormalizeData(
    object = obj,
    assay = 'ACTIVITY',
    normalization.method = 'LogNormalize',
    scale.factor = median(obj$nCount_ACTIVITY)
  )
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = rownames(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj), 
                reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, reduction='activity.pca', 
                 reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)
  rm(gene.activities)
  end_time <- Sys.time()  
  output[data,'activity']=round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "mins")),2)

  #############################################
  #############################################
  ##
  ##              Motif score
  ##
  #############################################
  #############################################
  start_time <- Sys.time()
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
  # add motif information
  DefaultAssay(obj)='ATAC'
  obj <- AddMotifs(
    object = obj,
    genome = genome,
    pfm = pfm
  )
  rm(pfm)
  
  # Computing motif activities
  obj <- RunChromVAR(
    object = obj,
    genome = genome,
    new.assay.name = 'MOTIF'
  )
  
  DefaultAssay(obj) <- 'MOTIF'
  #Using PCA dimensional reduction on motif modality and put it into Seurat Object
  pca_motif<-prcomp(t(obj@assays[["MOTIF"]]@data))$x[,1:50]
  obj[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
  rm(pca_motif)
  obj <- RunUMAP(obj, dims = 1:20, reduction='motif.pca', 
                 reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)
  end_time <- Sys.time()  
  output[data,'motif']=round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "mins")),2)
  
  #############################################
  #############################################
  ##
  ##              cisTopic LDA
  ##
  #############################################
  #############################################
  start_time <- Sys.time()
  DefaultAssay(obj) <- 'ATAC'
  temp = obj@assays$ATAC@counts
  rownames(temp)=GRangesToString(obj@assays$ATAC@ranges, sep=c(':','-'))
  cisTopicobj = createcisTopicObject(count.matrix = temp)
  rm(temp)
  
  cisTopicObject <- runWarpLDAModels(cisTopicobj, topic=c(10,20,40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
  cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
  cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')
  
  obj[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
  obj[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")
  end_time <- Sys.time()  
  output[data,'lda']=round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "mins")),2)
  
  
  #########################################################
  ##  Peak Accessibility: LSI
  #########################################################
  
  set.seed(1234)
  DefaultAssay(obj)='ATAC'
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:31, verbose = FALSE, graph.name=c('lsi_nn','lsi_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'lsi_snn')
  
  obj$signac_clusters=obj$seurat_clusters
  
  
  #########################################################
  ##  Peak Accessibility: LDA
  #########################################################
  
  DefaultAssay(obj)='ATAC'
  obj <- FindNeighbors(object = obj, reduction = 'lda', dims = 1:30, verbose = FALSE, graph.name=c('lda_nn','lda_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'lda_snn')
  
  obj$cistopic_clusters=obj$seurat_clusters
  
  
  #########################################################
  ##  Motif
  #########################################################
  
  DefaultAssay(obj)='MOTIF'
  obj <- FindNeighbors(object = obj, reduction = 'motif.pca', dims = 1:30, verbose = FALSE, graph.name=c('motif_nn','motif_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'motif_snn')
  
  obj$motif_clusters = obj$seurat_clusters
  

  #########################################################
  ##  Gene Activity
  #########################################################
  
  DefaultAssay(obj)='ACTIVITY'
  obj <- FindNeighbors(object = obj, reduction = 'activity.pca', dims = 1:30, verbose = FALSE, graph.name=c('activity_nn','activity_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name = 'activity_snn')
  
  obj$activity_clusters = obj$seurat_clusters
  
  
  #########################################################
  ##  Consensus PCA
  #########################################################
  start_time <- Sys.time()
  obj <- GetConsensusPCA(obj, reduction.list=list('lsi','lda','activity.pca','motif.pca'),
                  dims.list = list(2:31, 1:30, 1:30, 1:30),
                  reduction.name = 'consensus.pca',
                  reduction.key = "consensuspca_")
  
  obj <- RunUMAP(object = obj, reduction = 'consensus.pca', dims = 1:30, reduction.name = 'consensus.umap', 
                 reduction.key='consensusumap_', verbose = FALSE)
  
  DefaultAssay(obj)='ATAC'
  obj <- FindNeighbors(object = obj, reduction = 'consensus.pca', dims = 1:30, verbose = FALSE, graph.name=c('consensuspca_nn','consensuspca_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name='consensuspca_snn')
  
  obj$consensus_clusters=obj$seurat_clusters

  end_time <- Sys.time()  
  output[data,'cpca']=round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "mins")),2)
  
  #########################################################
  ##  Multiple CCA
  #########################################################
  start_time <- Sys.time()
  obj <- GetMultiCCA(obj, reduction.list=list('lsi', 'lda','activity.pca','motif.pca'),
                         dims.list = list(2:31, 1:30, 1:30, 1:30),
                         reduction.name = 'multicca.pca',
                         reduction.key = "multiccapca_")
  
  obj <- RunUMAP(obj, dims = 1:30, reduction='multicca.pca', 
                 reduction.name='multicca.umap', reduction.key = 'multiccaumap_', verbose = FALSE)
  
  obj <- FindNeighbors(object = obj, reduction = 'multicca.pca', dims = 1:30, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
  obj <- FindClusters(object = obj, algorithm = 1, verbose = FALSE, graph.name='multicca_snn')
  
  obj$cc_clusters = obj$seurat_clusters
  end_time <- Sys.time()  
  output[data,'mcca']=round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "mins")),2)
  
  
  #########################################################
  ##  Weighted Nearest Neighbor
  #########################################################
  start_time <- Sys.time()
  DefaultAssay(obj)='ATAC'
  obj <- FindMultiModalNeighbors(obj,
                                 reduction.list=list('lsi','lda','activity.pca','motif.pca'),
                                 dims.list = list(2:31, 1:30, 1:30, 1:30), verbose=FALSE)
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                 reduction.key = "wnnumap_", verbose = FALSE) # Run UMAP based on the knn distance
  obj <- FindClusters(object = obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  obj$wnn_clusters=obj$seurat_clusters
  end_time <- Sys.time()  
  output[data,'wnn']=round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "mins")),2)
  
  print(output)
  cat('\n\n')
}
save(output, file='running_time_output.rda')

write.csv(output, file='running_time_output.csv')
