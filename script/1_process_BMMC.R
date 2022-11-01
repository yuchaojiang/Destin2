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
library(SummarizedExperiment)
granja.atac.se=readRDS('BMMC/scATAC-Healthy-Hematopoiesis-191120_cellxpeak.rds')

colData(granja.atac.se)
dim(colData(granja.atac.se))
colData(granja.atac.se)$BioClassification

toplot=colData(granja.atac.se)
toplot=as.data.frame(toplot)
ggplot(toplot, aes(x=UMAP1, y=UMAP2)) + geom_point(aes(color=BioClassification), size=0.2)


metadata = toplot; rm(toplot)
counts=granja.atac.se@assays$data@listData$counts
rownames(counts)=names(granja.atac.se@rowRanges)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("_", "_"),
  genome = ref_genome,
  min.cells = 10,
  min.features = 200
)

granja.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data=metadata
)
rm(counts); rm(chrom_assay); rm(granja.atac.se); rm(metadata)
granja.atac
granja.atac[['ATAC']]

# Focus on BMMC only
table(granja.atac$orig.ident)
table(granja.atac$Group)
granja.atac=granja.atac[,granja.atac$Group=='BMMC_D6T1']
# Remove unknown cell types
granja.atac=granja.atac[,!grepl('Unk', granja.atac$BioClassification)]
table(granja.atac$BioClassification)
granja.atac=granja.atac[rowSums(granja.atac)>400,]

fragments <- CreateFragmentObject(
  path = "BMMC/GSM4138889_scATAC_BMMC_D6T1.fragments.tsv.gz",
  cells = granja.atac$Barcode,
  validate.fragments = TRUE,
  verbose = TRUE
)
granja.atac=SetAssayData(granja.atac, slot = "fragments", new.data = fragments)
bmmc=granja.atac
rm(granja.atac)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(bmmc) <- annotations
rm(annotations)

#############################################
#############################################
##
##              Quality controls
##
#############################################
#############################################

# compute nucleosome signal score per cell
bmmc <- NucleosomeSignal(object = bmmc)

# compute TSS enrichment score per cell
bmmc <- TSSEnrichment(object = bmmc, fast = FALSE)

VlnPlot(
  object = bmmc,
  features = c('nCount_ATAC', 'nFeature_ATAC',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

bmmc <- subset(
  x = bmmc,
  subset = nCount_ATAC < 25000 &
    nFeature_ATAC < 15000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
)
bmmc

save(bmmc, file='bmmc_raw.rda')

load('bmmc_raw.rda')

#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################
DefaultAssay(bmmc)='ATAC'
bmmc <- RunTFIDF(bmmc)
bmmc <- FindTopFeatures(bmmc, min.cutoff = 'q0')
bmmc <- RunSVD(bmmc) # This by default returns an lsi reduction

bmmc <- RunUMAP(object = bmmc, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                reduction.key='lsiumap_', verbose = FALSE)
bmmc <- FindNeighbors(object = bmmc, reduction = 'lsi', dims = 2:30)
bmmc <- FindClusters(object = bmmc, verbose = FALSE, algorithm = 3)
DimPlot(object = bmmc, label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(bmmc@assays$ATAC@data) # peak data
dim(bmmc@reductions$lsi@cell.embeddings) # peak lsi reduction

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

gene.activities <- GeneActivity(bmmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
bmmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(bmmc) <- "ACTIVITY"
bmmc <- NormalizeData(
  object = bmmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(bmmc$nCount_ACTIVITY)
)
bmmc <- FindVariableFeatures(bmmc, selection.method = "vst", nfeatures = 2000)
bmmc <- ScaleData(bmmc, features = rownames(bmmc))
bmmc <- RunPCA(bmmc, features = VariableFeatures(bmmc), 
               reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
bmmc <- RunUMAP(bmmc, dims = 1:30, reduction='activity.pca', 
                reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)
DimPlot(object = bmmc, label = TRUE, reduction = 'activity.umap') + NoLegend()

dim(bmmc@assays$ACTIVITY@scale.data) # Scaled gene activity
dim(bmmc@reductions$activity.pca@cell.embeddings) # Gene activity pc


#############################################
#############################################
##
##              Integrating with scRNA-seq
##
#############################################
#############################################

bmmc$predicted.id=bmmc$BioClassification
table(bmmc$seurat_clusters, bmmc$predicted.id)

table(bmmc$predicted.id)
bmmc=subset(bmmc, cells=which(bmmc$predicted.id!='10_cDC' &
                                bmmc$predicted.id!='18_Plasma' &
                                bmmc$predicted.id!='21_CD4.N2' &
                                bmmc$predicted.id!='23_CD8.EM' &
                                bmmc$predicted.id!='03_Late.Eryth'))

temp=sapply(strsplit(bmmc$predicted.id,'_'), "[[", 2)
temp[temp=='CLP.1']='CLP'
temp[temp=='CLP.2']='CLP'
temp[temp=='CD14.Mono.1']='CD14.Mono'
temp[temp=='CD14.Mono.2']='CD14.Mono'
temp[temp=='CD4.N1']='CD4.N'
bmmc$predicted.id=temp

#############################################
#############################################
##
##              Peak matrix processing
##
#############################################
#############################################
DefaultAssay(bmmc)='ATAC'
bmmc <- RunTFIDF(bmmc)
bmmc <- FindTopFeatures(bmmc, min.cutoff = 'q0')
bmmc <- RunSVD(bmmc) # This by default returns an lsi reduction

bmmc <- RunUMAP(object = bmmc, reduction = 'lsi', dims = 2:30, reduction.name = 'lsi.umap', 
                reduction.key='lsiumap_', verbose = FALSE)
bmmc <- FindNeighbors(object = bmmc, reduction = 'lsi', dims = 2:30)
bmmc <- FindClusters(object = bmmc, verbose = FALSE, algorithm = 3)
DimPlot(object = bmmc, label = TRUE, reduction = 'lsi.umap') + NoLegend()
DimPlot(object = bmmc, group.by='predicted.id', label = TRUE, reduction = 'lsi.umap') + NoLegend()

dim(bmmc@assays$ATAC@data) # peak data
dim(bmmc@reductions$lsi@cell.embeddings) # peak lsi reduction

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
gene.activities <- GeneActivity(bmmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
bmmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(bmmc) <- "ACTIVITY"
bmmc <- NormalizeData(
  object = bmmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(bmmc$nCount_ACTIVITY)
)
bmmc <- FindVariableFeatures(bmmc, selection.method = "vst", nfeatures = 2000)
bmmc <- ScaleData(bmmc, features = rownames(bmmc))
bmmc <- RunPCA(bmmc, features = VariableFeatures(bmmc), 
               reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
bmmc <- RunUMAP(bmmc, dims = 1:30, reduction='activity.pca', 
                reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

rm(gene.activities)
DimPlot(object = bmmc, label = TRUE, reduction = 'activity.umap') + NoLegend()
DimPlot(object = bmmc, group.by='predicted.id', label = TRUE, reduction = 'activity.umap') + NoLegend()

dim(bmmc@assays$ACTIVITY@scale.data) # Scaled gene activity
dim(bmmc@reductions$activity.pca@cell.embeddings) # Gene activity pc

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
DefaultAssay(bmmc)='ATAC'
bmmc <- AddMotifs(
  object = bmmc,
  genome = genome,
  pfm = pfm
)
rm(pfm)

# Computing motif activities
bmmc <- RunChromVAR(
  object = bmmc,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(bmmc) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
pca_motif<-prcomp(t(bmmc@assays[["MOTIF"]]@data))$x[,1:50]
bmmc[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
rm(pca_motif)
bmmc <- RunUMAP(bmmc, dims = 1:20, reduction='motif.pca', 
                reduction.name='motif.umap', reduction.key = 'motifumap_', verbose = FALSE)
DimPlot(object = bmmc, label = TRUE, reduction = 'motif.umap') + NoLegend()
DimPlot(object = bmmc, group.by='predicted.id',label = TRUE, reduction = 'motif.umap') + NoLegend()

dim(bmmc@assays$MOTIF@data) # Scaled motif scores
dim(bmmc@reductions$motif.pca@cell.embeddings) # Motif pc


#############################################
#############################################
##
##              cisTopic LDA
##
#############################################
#############################################

DefaultAssay(bmmc) <- 'ATAC'
temp = bmmc@assays$ATAC@counts
rownames(temp)=GRangesToString(bmmc@assays$ATAC@ranges, sep=c(':','-'))
cisTopicbmmc = createcisTopicObject(count.matrix = temp)
rm(temp)

cisTopicObject <- runWarpLDAModels(cisTopicbmmc, topic=c(10,20,40), seed=987, nCores=6, iterations = 500, addModels=FALSE)
cisTopicObject <- selectModel(cisTopicObject, type = 'derivative')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

dim(bmmc)
dim(cisTopicObject@selected.model$topics)
dim(cisTopicObject@selected.model$document_expects)

bmmc[['lda']] <- CreateDimReducObject(embeddings = t(cisTopicObject@selected.model$document_expects), key = "lda_", assay = "ATAC")
bmmc[["lda.umap"]] <- CreateDimReducObject(embeddings = cisTopicObject@dr$cell$Umap, key = "ldaumap_", assay = "ATAC")

save(bmmc, file='bmmc_processed.rda')
