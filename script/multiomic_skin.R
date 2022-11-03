#set working directory as needed

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
library(MAESTRO)
library(cisTopic)
library(mogsa)
library(BSgenome.Mmusculus.UCSC.mm10)
library(irlba)

skin=readRDS('skin.chromvar.rds')

DefaultAssay(skin) <- 'ACTIVITY'

skin <- RunPCA(skin, features = VariableFeatures(skin), 
                     reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)

DefaultAssay(skin) <- 'chromvar'


#Using PCA dimensional reduction on motif modality and put it into Seurat Object
pca_motif<-prcomp(t(skin@assays[["chromvar"]]@data))$x[,1:50]
skin[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "chromvar")
rm(pca_motif)
skin <- RunUMAP(skin, dims = 1:20, reduction='motif.pca', 
                      reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)


#############################################
#############################################
##
##              cisTopic LDA
##
#############################################
#############################################

DefaultAssay(skin) <- 'ATAC'
temp = skin@assays$ATAC@counts
rownames(temp)=GRangesToString(skin@assays$ATAC@ranges, sep=c(':','-'))
cisTopicskin = createcisTopicObject(count.matrix = temp)
rm(temp)

cisTopicObject <- runWarpLDAModels(cisTopicskin, topic=c(5, 10, 20, 25, 30, 40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

dim(skin)
dim(cisTopicObject@selected.model$topics)
dim(cisTopicObject@selected.model$document_expects)

skin[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
skin[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")


#########################################
#########################################
##
##          Joint Dimension Reduction
##
#########################################
#########################################

skin <- FindMultiModalNeighbors(
  skin, reduction.list = list("lsi", "activity.pca", "motif.pca", "lda"), 
  dims.list = list(2:31, 1:30, 1:30, 1:30))


skin <- RunUMAP(skin, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
skin <- FindClusters(skin, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


skin_mcca <- mbpca(x = list(t(skin[["lsi"]]@cell.embeddings[, 2:31]), 
                                  t(skin[["activity.pca"]]@cell.embeddings[, 1:30]),
                                  t(skin[["motif.pca"]]@cell.embeddings[, 1:30]), 
                                  t(skin[["lda"]]@cell.embeddings[, 1:30])),
                         method = "blockScore", verbose = T, ncomp = 30, moa = F, scale = T)

skin_cpca <- prcomp(cbind(skin[["lsi"]]@cell.embeddings[, 2:31], 
                                skin[["activity.pca"]]@cell.embeddings[, 1:30],
                                skin[["motif.pca"]]@cell.embeddings[, 1:30], 
                                skin[["lda"]]@cell.embeddings[, 1:30]), scale = T)$x[, 1:30]


rownames(skin_mcca$t) <- colnames(skin)


skin[["CPCA"]] <- CreateDimReducObject(embeddings = skin_cpca, key = "CPCA_", assay = "ATAC")
skin[["MCCA"]] <- CreateDimReducObject(embeddings = skin_mcca$t, key = "MCCA_", assay = "ATAC")

save(skin, file = 'skin_processed.rda')

