setwd("C:/Users/yuchaoj/Dropbox/ATAC_multi/data_script/")
setwd("~/Dropbox/ATAC_multi/data_script/")

ref_genome='mm10'

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
counts <- Read10X_h5(filename = "Adult_mouse_brain/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5")

#Read in the metadata for cells
metadata <- read.csv(
  file = "Adult_mouse_brain/atac_v1_adult_brain_fresh_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

#Create Seurat Object for analysis
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = ref_genome,
  fragments = 'Adult_mouse_brain/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

brain <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = metadata
)

brain[['ATAC']]
rm(counts); rm(chrom_assay); rm(metadata)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(brain) <- annotations
rm(annotations)

#############################################
#############################################
##
##              Quality controls
##
#############################################
#############################################

# compute nucleosome signal score per cell
brain <- NucleosomeSignal(object = brain)

# compute TSS enrichment score per cell
brain <- TSSEnrichment(object = brain, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
brain$pct_reads_in_peaks <- brain$peak_region_fragments / brain$passed_filters * 100
brain$blacklist_ratio <- brain$blacklist_region_fragments / brain$peak_region_fragments

VlnPlot(
  object = brain,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

brain <- subset(
  x = brain,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
brain

save(brain, file='brain_raw.rda')

load('brain_raw.rda')

#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################
DefaultAssay(brain)='ATAC'
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(brain) # This by default returns an lsi reduction

brain <- RunUMAP(object = brain, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                reduction.key='lsiumap_', verbose = FALSE)
brain <- FindNeighbors(object = brain, reduction = 'lsi', dims = 2:30)
brain <- FindClusters(object = brain, verbose = FALSE, algorithm = 3)
DimPlot(object = brain, label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(brain@assays$ATAC@data) # peak data
dim(brain@reductions$lsi@cell.embeddings) # peak lsi reduction

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

gene.activities <- GeneActivity(brain)
# # This can be substituted with MAESTRO's regulatory potential model
# # python3.8 -m pip install tables
# # python3.8 -m pip install h5py
# # python3.8 -m pip install scipy
# # python3.8 -m pip install pandas
# reticulate::use_python('/Users/yuchaojiang/pythonProject/bin/python', required = TRUE)
# gene.activities <- ATACCalculateGenescore(brain@assays$ATAC@counts, organism = "GRCm38")

# add the gene activity matrix to the Seurat object as a new assay and normalize it
brain[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(brain) <- "ACTIVITY"
brain <- NormalizeData(
  object = brain,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_ACTIVITY)
)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <- ScaleData(brain, features = rownames(brain))
brain <- RunPCA(brain, features = VariableFeatures(brain), 
                reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
brain <- RunUMAP(brain, dims = 1:30, reduction='activity.pca', 
                 reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)
DimPlot(object = brain, label = TRUE, reduction = 'activity.umap') + NoLegend()

dim(brain@assays$ACTIVITY@scale.data) # Scaled gene activity
dim(brain@reductions$activity.pca@cell.embeddings) # Gene activity pc


#############################################
#############################################
##
##              Integrating with scRNA-seq
##
#############################################
#############################################

# Load the pre-processed scRNA-seq data for brains
# https://signac-objects.s3.amazonaws.com/allen_brain.rds
brain_rna <- readRDS("Adult_mouse_brain/allen_brain.rds")
DefaultAssay(brain)='ACTIVITY'
transfer.anchors <- FindTransferAnchors(
  reference = brain_rna,
  query = brain,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = brain_rna$subclass,
  weight.reduction = brain[['lsi']],
  dims = 2:30
)

brain <- AddMetaData(object = brain, metadata = predicted.labels)
plot1 <- DimPlot(
  object = brain_rna,
  group.by = 'subclass',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = brain,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

table(brain$seurat_clusters, brain$predicted.id)

# Remove cells with low label-transfer scores
hist(brain$prediction.score.max)
brain=subset(brain, cells=which(brain$prediction.score.max>0.8))
hist(brain$prediction.score.max)


#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################
DefaultAssay(brain)='ATAC'
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(brain) # This by default returns an lsi reduction

brain <- RunUMAP(object = brain, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                reduction.key='lsiumap_', verbose = FALSE)
brain <- FindNeighbors(object = brain, reduction = 'lsi', dims = 2:30)
brain <- FindClusters(object = brain, verbose = FALSE, algorithm = 3)
DimPlot(object = brain, label = TRUE, reduction = 'lsi.umap') + NoLegend()
DimPlot(object = brain, group.by='predicted.id', label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(brain@assays$ATAC@data) # peak data
dim(brain@reductions$lsi@cell.embeddings) # peak lsi reduction

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
gene.activities <- GeneActivity(brain)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
brain[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(brain) <- "ACTIVITY"
brain <- NormalizeData(
  object = brain,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_ACTIVITY)
)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <- ScaleData(brain, features = rownames(brain))
brain <- RunPCA(brain, features = VariableFeatures(brain), 
                reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
brain <- RunUMAP(brain, dims = 1:30, reduction='activity.pca', 
                 reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)
DimPlot(object = brain, label = TRUE, reduction = 'activity.umap') + NoLegend()
DimPlot(object = brain, group.by='predicted.id', label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(brain@assays$ACTIVITY@scale.data) # Scaled gene activity
dim(brain@reductions$activity.pca@cell.embeddings) # Gene activity pc

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
DefaultAssay(brain)='ATAC'
brain <- AddMotifs(
  object = brain,
  genome = genome,
  pfm = pfm
)
rm(pfm)

# Computing motif activities
brain <- RunChromVAR(
  object = brain,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(brain) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
pca_motif<-prcomp(t(brain@assays[["MOTIF"]]@data))$x[,1:50]
brain[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
rm(pca_motif)
brain <- RunUMAP(brain, dims = 1:20, reduction='motif.pca', 
                 reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)
DimPlot(object = brain, label = TRUE, reduction = 'motif.umap') + NoLegend()
DimPlot(object = brain, group.by='predicted.id', label = TRUE, reduction = 'motif.umap') + NoLegend()

dim(brain@assays$MOTIF@data) # Scaled motif scores
dim(brain@reductions$motif.pca@cell.embeddings) # Motif pc

#############################################
#############################################
##
##              cisTopic LDA
##
#############################################
#############################################

DefaultAssay(brain) <- 'ATAC'
temp = brain@assays$ATAC@counts
rownames(temp)=GRangesToString(brain@assays$ATAC@ranges, sep=c(':','-'))
cisTopicbrain = createcisTopicObject(count.matrix = temp)
rm(temp)

cisTopicObject <- runWarpLDAModels(cisTopicbrain, topic=c(10,20,40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

dim(brain)
dim(cisTopicObject@selected.model$topics) # 40 topics
dim(cisTopicObject@selected.model$document_expects)

brain[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
brain[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")

save(brain, file='brain_processed.rda')
