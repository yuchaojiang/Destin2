setwd("C:/Users/yuchaoj/Dropbox/ATAC_multi/data_script/")
setwd("~/Dropbox/ATAC_multi/data_script/")

ref_genome='hg19'

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
library(umap)


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

set.seed(1234)

#############################################
#############################################
##
##              Pre-processing
##
#############################################
#############################################

#Read in the raw counts data of scATAC-seq data
counts <- Read10X_h5(filename = "PBMC/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")

#Read in the metadata for cells
metadata <- read.csv(
  file = "PBMC/atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

#Create Seurat Object for analysis
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = ref_genome,
  fragments = 'PBMC/atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = metadata
)

pbmc[['ATAC']]
rm(counts); rm(chrom_assay); rm(metadata)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(pbmc) <- annotations
rm(annotations)

#############################################
#############################################
##
##              Quality controls
##
#############################################
#############################################

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc

save(pbmc, file='pbmc_raw.rda')

load('pbmc_raw.rda')

#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################
DefaultAssay(pbmc)='ATAC'
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc) # This by default returns an lsi reduction

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                reduction.key='lsiumap_', verbose = FALSE)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(pbmc@assays$ATAC@data) # peak data
dim(pbmc@reductions$lsi@cell.embeddings) # peak lsi reduction

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

gene.activities <- GeneActivity(pbmc)
# This can be substituted with MAESTRO's regulatory potential model


# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(pbmc) <- "ACTIVITY"
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_ACTIVITY)
)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), 
               reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction='activity.pca', 
                reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)
DimPlot(object = pbmc, label = TRUE, reduction = 'activity.umap') + NoLegend()

dim(pbmc@assays$ACTIVITY@scale.data) # Scaled gene activity
dim(pbmc@reductions$activity.pca@cell.embeddings) # Gene activity pc


#############################################
#############################################
##
##              Integrating with scRNA-seq
##
#############################################
#############################################

# Load the pre-processed scRNA-seq data for PBMCs
# https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds 
pbmc_rna <- readRDS("PBMC/pbmc_10k_v3.rds")
DefaultAssay(pbmc)='ACTIVITY'
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

table(pbmc$seurat_clusters, pbmc$predicted.id)
# Remove cluster 14 with lower QC
pbmc <- subset(pbmc, idents = 14, invert = TRUE)
pbmc$seurat_clusters=droplevels(pbmc$seurat_clusters)

# Remove cells with low label-transfer scores
hist(pbmc$prediction.score.max)
pbmc=subset(pbmc, cells=which(pbmc$prediction.score.max>0.8))
hist(pbmc$prediction.score.max)

table(pbmc$predicted.id)
pbmc=subset(pbmc, cells=which(pbmc$predicted.id!='Dendritic cell'))
table(pbmc$predicted.id)


#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################
DefaultAssay(pbmc)='ATAC'
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc) # This by default returns an lsi reduction

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                reduction.key='lsiumap_', verbose = FALSE)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE, reduction = 'lsi.umap') + NoLegend()
DimPlot(object = pbmc, group.by='predicted.id', label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(pbmc@assays$ATAC@data) # peak data
dim(pbmc@reductions$lsi@cell.embeddings) # peak lsi reduction

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
gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(pbmc) <- "ACTIVITY"
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_ACTIVITY)
)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), 
               reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction='activity.pca', 
                reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)
DimPlot(object = pbmc, label = TRUE, reduction = 'activity.umap') + NoLegend()
DimPlot(object = pbmc, group.by='predicted.id', label = TRUE, reduction = 'activity.umap') + NoLegend()

dim(pbmc@assays$ACTIVITY@scale.data) # Scaled gene activity
dim(pbmc@reductions$activity.pca@cell.embeddings) # Gene activity pc

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

# add motif information
DefaultAssay(pbmc)='ATAC'
pbmc <- AddMotifs(
  object = pbmc,
  genome = genome,
  pfm = pfm
)
rm(pfm)

# Computing motif activities
pbmc <- RunChromVAR(
  object = pbmc,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(pbmc) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
pca_motif<-prcomp(t(pbmc@assays[["MOTIF"]]@data))$x[,1:50]
pbmc[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
rm(pca_motif)
pbmc <- RunUMAP(pbmc, dims = 1:20, reduction='motif.pca', 
                reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)
DimPlot(object = pbmc, label = TRUE, reduction = 'motif.umap') + NoLegend()
DimPlot(object = pbmc, group.by='predicted.id',label = TRUE, reduction = 'motif.umap') + NoLegend()

dim(pbmc@assays$MOTIF@data) # Scaled motif scores
dim(pbmc@reductions$motif.pca@cell.embeddings) # Motif pc


#############################################
#############################################
##
##              cisTopic LDA
##
#############################################
#############################################

DefaultAssay(pbmc) <- 'ATAC'
temp = pbmc@assays$ATAC@counts
rownames(temp)=GRangesToString(pbmc@assays$ATAC@ranges, sep=c(':','-'))
cisTopicpbmc = createcisTopicObject(count.matrix = temp)
rm(temp)

cisTopicObject <- runWarpLDAModels(cisTopicpbmc, topic=c(10,20,40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

dim(pbmc)
dim(cisTopicObject@selected.model$topics)
dim(cisTopicObject@selected.model$document_expects)

pbmc[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
pbmc[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")

save(pbmc, file='pbmc_processed.rda')

