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

load("Human/gse.rda")

#Read in the raw counts data of scATAC-seq data
counts = gse@assays$peaks@counts

#Read in the metadata for cells
metadata <- gse@meta.data

#Create Seurat Object for analysis
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = ref_genome,
  fragments = 'Human/fragment.txt.gz',
  min.cells = 10,
  min.features = 200
)

human <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = metadata
)


human[['ATAC']]
rm(counts); rm(chrom_assay); rm(metadata)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(human) <- annotations
rm(annotations)

#############################################
#############################################
##
##              Quality controls
##
#############################################
#############################################

# compute nucleosome signal score per cell
human <- NucleosomeSignal(object = human)

# compute TSS enrichment score per cell
human <- TSSEnrichment(object = human, fast = TRUE)

# add blacklist ratio and fraction of reads in peaks
human$pct_reads_in_peaks <- human$frip * 100 #frip = fraction of reads in peaks
human$blacklist_ratio <- human$blacklist_fraction 
human$peak_region_fragments <- human$nCount_ATAC

VlnPlot(
  object = human,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

human <- subset(
  x = human,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
human

save(human, file='human_raw.rda')

load('human_raw.rda')

#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################
DefaultAssay(human)='ATAC'
human <- RunTFIDF(human)
human <- FindTopFeatures(human, min.cutoff = 'q0')
human <- RunSVD(human) # This by default returns an lsi reduction

human <- RunUMAP(object = human, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                 reduction.key='lsiumap_', verbose = FALSE)
human <- FindNeighbors(object = human, reduction = 'lsi', dims = 2:30)
human <- FindClusters(object = human, verbose = FALSE, algorithm = 3)
DimPlot(object = human, label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(human@assays$ATAC@data) # peak data
dim(human@reductions$lsi@cell.embeddings) # peak lsi reduction

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

gene.activities <- GeneActivity(human)
# # This can be substituted with MAESTRO's regulatory potential model
# # python3.8 -m pip install tables
# # python3.8 -m pip install h5py
# # python3.8 -m pip install scipy
# # python3.8 -m pip install pandas
# reticulate::use_python('/Users/yuchaojiang/pythonProject/bin/python', required = TRUE)
# gene.activities <- ATACCalculateGenescore(human@assays$ATAC@counts, organism = "GRCm38")

# add the gene activity matrix to the Seurat object as a new assay and normalize it
human[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(human) <- "ACTIVITY"
human <- NormalizeData(
  object = human,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(human$nCount_ACTIVITY)
)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)
human <- ScaleData(human, features = rownames(human))
human <- RunPCA(human, features = VariableFeatures(human), 
                reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
human <- RunUMAP(human, dims = 1:30, reduction='activity.pca', 
                 reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)
DimPlot(object = human, label = TRUE, reduction = 'activity.umap') + NoLegend()

dim(human@assays$ACTIVITY@scale.data) # Scaled gene activity
dim(human@reductions$activity.pca@cell.embeddings) # Gene activity pc


plot1 <- DimPlot(
  object = human,
  group.by = 'tissue',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('Tissue')

plot2 <- DimPlot(
  object = human,
  group.by = 'seurat_clusters',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

table(human$seurat_clusters, human$tissue)

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
DefaultAssay(human)='ATAC'
human <- AddMotifs(
  object = human,
  genome = genome,
  pfm = pfm
)
rm(pfm)

# Computing motif activities
human <- RunChromVAR(
  object = human,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(human) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
pca_motif<-prcomp(t(human@assays[["MOTIF"]]@data))$x[,1:50]
human[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
rm(pca_motif)
human <- RunUMAP(human, dims = 1:20, reduction='motif.pca', 
                 reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)
DimPlot(object = human, label = TRUE, reduction = 'motif.umap') + NoLegend()
#DimPlot(object = human, group.by='predicted.id', label = TRUE, reduction = 'motif.umap') + NoLegend()

dim(human@assays$MOTIF@data) # Scaled motif scores
dim(human@reductions$motif.pca@cell.embeddings) # Motif pc

#############################################
#############################################
##
##              cisTopic LDA
##
#############################################
#############################################

DefaultAssay(human) <- 'ATAC'
temp = human@assays$ATAC@counts
rownames(temp)=GRangesToString(human@assays$ATAC@ranges, sep=c(':','-'))
cisTopichuman = createcisTopicObject(count.matrix = temp)
rm(temp)

cisTopicObject <- runWarpLDAModels(cisTopichuman, topic=c(10,20,40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = human@meta.data)

cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

dim(human)
dim(cisTopicObject@selected.model$topics) # 40 topics
dim(cisTopicObject@selected.model$document_expects)

human[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
human[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")

human$predicted.id=factor(human$tissue)
save(human, file='human_processed.rda')

