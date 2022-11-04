# set working directory as needed

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
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(irlba)

set.seed(1234)

ref_genome='hg38'

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


# load the RNA and ATAC data
counts <- Read10X_h5("./pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "./pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotationsg
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
multi_pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
multi_pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)



DefaultAssay(multi_pbmc) <- "ATAC"

multi_pbmc <- NucleosomeSignal(multi_pbmc)
multi_pbmc <- TSSEnrichment(multi_pbmc)

# filter out low quality cells
multi_pbmc <- subset(
  x = multi_pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)



#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################

DefaultAssay(multi_pbmc) <- "ATAC"
multi_pbmc <- RunTFIDF(multi_pbmc)
multi_pbmc <- FindTopFeatures(multi_pbmc, min.cutoff = 5)
multi_pbmc <- RunSVD(multi_pbmc)

DepthCor(multi_pbmc) ## 1st pc very correlated with sequencing depth, as expected

multi_pbmc <- RunUMAP(object = multi_pbmc, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                reduction.key='lsiumap_', verbose = FALSE)

#############################################
#############################################
##
##              Gene activity processing
##
#############################################
#############################################


gene.activities <- GeneActivity(multi_pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
multi_pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(multi_pbmc) <- "ACTIVITY"
multi_pbmc <- NormalizeData(
  object = multi_pbmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(multi_pbmc$nCount_ACTIVITY)
)
multi_pbmc <- FindVariableFeatures(multi_pbmc, selection.method = "vst", nfeatures = 2000)
multi_pbmc <- ScaleData(multi_pbmc, features = rownames(multi_pbmc))
multi_pbmc <- RunPCA(multi_pbmc, features = VariableFeatures(multi_pbmc), 
               reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
multi_pbmc <- RunUMAP(multi_pbmc, dims = 1:30, reduction='activity.pca', 
                reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)


#############################################
#############################################
##
##              Motif score
##
#############################################
#############################################

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

main.chroms <- standardChromosomes(genome)
keep.peaks <- which(as.character(seqnames(granges(multi_pbmc))) %in% main.chroms)
multi_pbmc[["ATAC"]] <- subset(multi_pbmc[["ATAC"]], features = rownames(multi_pbmc[["ATAC"]])[keep.peaks])

rm(main.chroms)
rm(keep.peaks)


# add motif information
DefaultAssay(multi_pbmc)='ATAC'
multi_pbmc <- AddMotifs(
  object = multi_pbmc,
  genome = genome,
  pfm = pfm
)
rm(pfm)

# Computing motif activities
multi_pbmc <- RunChromVAR(
  object = multi_pbmc,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(multi_pbmc) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
pca_motif<-prcomp(t(multi_pbmc@assays[["MOTIF"]]@data))$x[,1:50]
multi_pbmc[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
rm(pca_motif)
multi_pbmc <- RunUMAP(multi_pbmc, dims = 1:20, reduction='motif.pca', 
                reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)


#############################################
#############################################
##
##              cisTopic LDA
##
#############################################
#############################################

DefaultAssay(multi_pbmc) <- 'ATAC'
temp = multi_pbmc@assays$ATAC@counts
rownames(temp)=GRangesToString(multi_pbmc@assays$ATAC@ranges, sep=c(':','-'))
cisTopicmulti_pbmc = createcisTopicObject(count.matrix = temp)
rm(temp)

cisTopicObject <- runWarpLDAModels(cisTopicmulti_pbmc, topic=c(5, 10, 20, 25, 30, 40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

dim(multi_pbmc)
dim(cisTopicObject@selected.model$topics)
dim(cisTopicObject@selected.model$document_expects)

multi_pbmc[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
multi_pbmc[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")

save(multi_pbmc, file='multi_pbmc_processed.rda')


#########################################
#########################################
##
##          Joint Dimension Reduction
##
#########################################
#########################################

multi_pbmc <- FindMultiModalNeighbors(
  multi_pbmc, reduction.list = list("lsi", "activity.pca", "motif.pca", "lda"), 
  dims.list = list(2:31, 1:30, 1:30, 1:30))


multi_pbmc <- RunUMAP(multi_pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
multi_pbmc <- FindClusters(multi_pbmc, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


multi_pbmc_mcca <- mbpca(x = list(t(multi_pbmc[["lsi"]]@cell.embeddings[, 2:31]), 
                                   t(multi_pbmc[["activity.pca"]]@cell.embeddings[, 1:30]),
                                   t(multi_pbmc[["motif.pca"]]@cell.embeddings[, 1:30]), 
                                   t(multi_pbmc[["lda"]]@cell.embeddings[, 1:30])),
                          method = "blockScore", verbose = T, ncomp = 30, moa = F, scale = T)

multi_pbmc_cpca <- prcomp(cbind(multi_pbmc[["lsi"]]@cell.embeddings[, 2:31], 
                                 multi_pbmc[["activity.pca"]]@cell.embeddings[, 1:30],
                                 multi_pbmc[["motif.pca"]]@cell.embeddings[, 1:30], 
                                 multi_pbmc[["lda"]]@cell.embeddings[, 1:30]), scale = T)$x[, 1:30]


rownames(multi_pbmc_mcca$t) <- colnames(multi_pbmc)


multi_pbmc[["CPCA"]] <- CreateDimReducObject(embeddings = multi_pbmc_cpca, key = "CPCA_", assay = "ATAC")
multi_pbmc[["MCCA"]] <- CreateDimReducObject(embeddings = multi_pbmc_mcca$t, key = "MCCA_", assay = "ATAC")



#########################################
#########################################
#### RNA Preprocessing & Reduction
#########################################
#########################################

DefaultAssay(multi_pbmc) <- 'RNA'

multi_pbmc <- NormalizeData(multi_pbmc)
multi_pbmc <- FindVariableFeatures(multi_pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(multi_pbmc)
multi_pbmc <- ScaleData(multi_pbmc, features = all.genes)
multi_pbmc <- RunPCA(multi_pbmc, features = VariableFeatures(object = multi_pbmc))

save(multi_pbmc, file='multi_pbmc_processed.rda')





