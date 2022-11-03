## set working directory as needed

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
library(fields)
library(kableExtra)

set.seed(1234)

ref_genome='mm10'

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
counts <- Read10X_h5("./e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
fragpath <- "./e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz"

# get gene annotationsg
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
multi_brain <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
multi_brain[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

rm(counts)

DefaultAssay(multi_brain) <- "ATAC"

multi_brain <- NucleosomeSignal(multi_brain)
multi_brain <- TSSEnrichment(multi_brain)

VlnPlot(
  object = multi_brain,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
multi_brain <- subset(
  x = multi_brain,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 50000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
)


#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################

DefaultAssay(multi_brain) <- "ATAC"
multi_brain <- RunTFIDF(multi_brain)
multi_brain <- FindTopFeatures(multi_brain, min.cutoff = 5)
multi_brain <- RunSVD(multi_brain)

DepthCor(multi_brain) ## 1st pc very correlated with sequencing depth, as expected

multi_brain <- RunUMAP(object = multi_brain, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                   reduction.key='lsiumap_', verbose = FALSE)

#############################################
#############################################
##
##              Gene activity processing
##
#############################################
#############################################


gene.activities <- GeneActivity(multi_brain)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
multi_brain[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(multi_brain) <- "ACTIVITY"
multi_brain <- NormalizeData(
  object = multi_brain,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(multi_brain$nCount_ACTIVITY)
)
multi_brain <- FindVariableFeatures(multi_brain, selection.method = "vst", nfeatures = 2000)
multi_brain <- ScaleData(multi_brain, features = rownames(multi_brain))
multi_brain <- RunPCA(multi_brain, features = VariableFeatures(multi_brain), 
                  reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
multi_brain <- RunUMAP(multi_brain, dims = 1:30, reduction='activity.pca', 
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

std.chroms <- standardChromosomes(genome)
keep.seq <- which(as.character(seqnames(granges(multi_brain))) %in% std.chroms)
multi_brain[["ATAC"]] <- subset(multi_brain[["ATAC"]], features = rownames(multi_brain[["ATAC"]])[keep.seq])

rm(main.chroms)
rm(keep.peaks)


# add motif information
DefaultAssay(multi_brain)='ATAC'
multi_brain <- AddMotifs(
  object = multi_brain,
  genome = genome,
  pfm = pfm
)
rm(pfm)

# Computing motif activities
multi_brain <- RunChromVAR(
  object = multi_brain,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(multi_brain) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
pca_motif<-prcomp(t(multi_brain@assays[["MOTIF"]]@data))$x[,1:50]
multi_brain[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
rm(pca_motif)
multi_brain <- RunUMAP(multi_brain, dims = 1:20, reduction='motif.pca', 
                   reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)


#############################################
#############################################
##
##              cisTopic LDA
##
#############################################
#############################################

DefaultAssay(multi_brain) <- 'ATAC'
temp = multi_brain@assays$ATAC@counts
rownames(temp)=GRangesToString(multi_brain@assays$ATAC@ranges, sep=c(':','-'))
cisTopicmulti_brain = createcisTopicObject(count.matrix = temp)
rm(temp)

cisTopicObject <- runWarpLDAModels(cisTopicmulti_brain, topic=c(5, 10, 20, 25, 30, 40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

dim(multi_brain)
dim(cisTopicObject@selected.model$topics)
dim(cisTopicObject@selected.model$document_expects)

multi_brain[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
multi_brain[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")


###################


#########################################
#########################################
##
##          Joint Dimension Reduction
##
#########################################
#########################################

multi_brain <- FindMultiModalNeighbors(
  multi_brain, reduction.list = list("lsi", "activity.pca", "motif.pca", "lda"), 
  dims.list = list(2:31, 1:30, 1:30, 1:30))


multi_brain <- RunUMAP(multi_brain, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
multi_brain <- FindClusters(multi_brain, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


multi_brain_mcca <- mbpca(x = list(t(multi_brain[["lsi"]]@cell.embeddings[, 2:31]), 
                             t(multi_brain[["activity.pca"]]@cell.embeddings[, 1:30]),
                             t(multi_brain[["motif.pca"]]@cell.embeddings[, 1:30]), 
                             t(multi_brain[["lda"]]@cell.embeddings[, 1:30])),
                    method = "blockScore", verbose = T, ncomp = 30, moa = F, scale = T)

multi_brain_cpca <- prcomp(cbind(multi_brain[["lsi"]]@cell.embeddings[, 2:31], 
      multi_brain[["activity.pca"]]@cell.embeddings[, 1:30],
      multi_brain[["motif.pca"]]@cell.embeddings[, 1:30], 
      multi_brain[["lda"]]@cell.embeddings[, 1:30]), scale = T)$x[, 1:30]


rownames(multi_brain_mcca$t) <- colnames(multi_brain)


multi_brain[["CPCA"]] <- CreateDimReducObject(embeddings = multi_brain_cpca, key = "CPCA_", assay = "ATAC")
multi_brain[["MCCA"]] <- CreateDimReducObject(embeddings = multi_brain_mcca$t, key = "MCCA_", assay = "ATAC")



#############################################################
### Process RNA for FOSCTTNN
#############################################################
DefaultAssay(multi_brain) <- 'RNA'

multi_brain <- NormalizeData(multi_brain)
multi_brain <- FindVariableFeatures(multi_brain, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(multi_brain)
multi_brain <- ScaleData(multi_brain, features = all.genes)
multi_brain <- RunPCA(multi_brain, features = VariableFeatures(object = multi_brain))

save(multi_brain, file='multi_brain_processed.rda')








