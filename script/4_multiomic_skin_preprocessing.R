
ct=read.table('celltype_skin.txt', header = T, sep='\t')
# first column: atac barcode
# second column: rna barcode
# third column: annotated cell type from the paper

rna.count=read.table('GSM4156608_skin.late.anagen.rna.counts.txt.gz', header = T, row.names = 1)
rna.count[1:5, 1:5]

# match the rna.count with the cell type labels
any(is.na(match(ct$rna.bc, colnames(rna.count))))
rna.count=rna.count[,match(ct$rna.bc, colnames(rna.count))]

# The barcode file is the same as the annotation ct; omitted
peak.barcodes <- scan('GSM4156597_skin.late.anagen.barcodes.txt.gz', what="")
all(ct$rna.bc==peak.barcodes)
rm(peak.barcodes)

# Read in the peak matrix: stored as Market matrix
# https://atlas.gs.washington.edu/mouse-atac/docs/
library(Matrix)
library(GenomicRanges)
peak.count= readMM("GSM4156597_skin.late.anagen.counts.txt.gz")
dim(peak.count)
peak.bed=read.delim('GSM4156597_skin.late.anagen.peaks.bed.gz', header=FALSE)
peak.granges=GRanges(seqnames=peak.bed$V1, ranges=IRanges(st=peak.bed$V2, end=peak.bed$V3))
rm(peak.bed)
rownames(peak.count)=paste(peak.granges)
colnames(peak.count)=ct$rna.bc # we will use the rna barcode from here and on

dim(ct) # barcode and cell type
dim(peak.count) # peak count
length(peak.granges) # GRanges for peaks
dim(rna.count) # rna count


library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse skin
library(dplyr)
library(ggplot2)

# Create Seurat object
skin <- CreateSeuratObject(counts = rna.count)
skin[["percent.mt"]] <- PercentageFeatureSet(skin, pattern = "^MT-")
rm(rna.count)

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.use <- seqnames(peak.granges) %in% standardChromosomes(peak.granges)
peak.count <- peak.count[as.vector(grange.use), ]
peak.granges <- peak.granges[as.vector(grange.use)]
rm(grange.use)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

rm(peak.granges)
peak.count=peak.count[rowSums(peak.count)>10,]
# Leave the fragment file out for now.
# frag.file <- "GSM4156599_skin.atac.fragments.bed.gz"
chrom_assay <- CreateChromatinAssay(
  counts = peak.count,
  sep = c(":", "-"),
  genome = 'mm10',
  # fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
skin[["ATAC"]] <- chrom_assay
skin@assays$RNA
skin@assays$ATAC
rm(peak.count)
save(skin, file='skin_cluster.rda')


library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse skin
library(dplyr)
library(ggplot2)
load('skin_cluster.rda')

## perform basic qc based on the number of detected molecules for each modality
## as well as mitochondrial percentage
VlnPlot(skin, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
skin <- subset(
  x = skin,
  subset = nCount_ATAC < 25000 &
    nCount_ATAC > 500 &
    nCount_RNA < 10000 &
    nCount_RNA > 500 &
    percent.mt < 20
)
VlnPlot(skin, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

## perform pre-processing and dimensional reduction on both assays independently
## using standard approaches for rna and atac-seq data
DefaultAssay(skin) <- "RNA"
skin <- SCTransform(skin, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(skin) <- "ATAC"
skin <- RunTFIDF(skin)
skin <- FindTopFeatures(skin, min.cutoff = 'q0')
skin <- RunSVD(skin)
skin <- RunUMAP(skin, reduction = 'lsi', dims = 2:50,
                reduction.name = "umap.atac",
                reduction.key = "atacUMAP_")

## calculate a wnn graph, representing a weighted combination of rna and atac-seq modalities
## use this graph for umap visualization and clustering
skin <- FindMultiModalNeighbors(skin, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
skin <- RunUMAP(skin, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
skin <- FindClusters(skin, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
table(skin$seurat_clusters)

p1 <- DimPlot(skin, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(skin, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(skin, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# We will now get the gene activity matrix
DefaultAssay(skin) <- "ATAC"
# gene.activities <- GeneActivity(skin) # This would not run because it requires fragment file to be loaded.

# Below is a short-cut but is not very precise: we only aggregate reads in the peak regions
# But there are reads off peak regions but within gene bodies that are defined
# Will load the fragment file in the second script.
annotation <- Annotation(object = skin)
transcripts <- Signac:::CollapseToLongestTranscript(ranges = annotation)
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"]
transcripts <- Extend(x = transcripts, upstream = 2000, 
                      downstream = 0)
gene.activities=matrix(nrow=length(transcripts), ncol=ncol(skin))
rownames(gene.activities)=transcripts$gene_name
colnames(gene.activities)=colnames(skin)
peak.count=skin@assays$ATAC@counts
peak.granges=skin@assays$ATAC@ranges
for(i in 1:nrow(gene.activities)){
  if(i %%200 ==0) cat(i,' ')
  peak.index=which(countOverlaps(peak.granges, transcripts[i])>0)
  gene.activities[i,]=apply(peak.count[peak.index,, drop=FALSE],2, sum)
}
gene.activities=gene.activities[apply(gene.activities,1,sum)>10,]
gene.activities=gene.activities[!duplicated(rownames(gene.activities)),] # remove duplicate genes
gene.activities=gene.activities[-which(rownames(gene.activities)==''),]

skin[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)

# Perform normalization/scaling of the gene activity matrix
DefaultAssay(skin) <- "ACTIVITY"
skin <- NormalizeData(skin)
skin <- FindVariableFeatures(skin)
all.genes <- rownames(skin)
skin <- ScaleData(skin, features = all.genes)
skin@assays$ACTIVITY@scale.data[1:5,1:5]
save(skin, file='skin.rda')


wnn.resolution=5
skin <- FindClusters(skin, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
p5=DimPlot(skin, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
        label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(skin@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()
p5

wnn.resolution=15
skin <- FindClusters(skin, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
p15=DimPlot(skin, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
            label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(skin@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()
p15

wnn.resolution=25
skin <- FindClusters(skin, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
p25=DimPlot(skin, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
            label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(skin@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()
p25


p5|p15|p25



dir.out <- "./"
dir.in <- "./"

library(Seurat)
library(Signac)

## add a fragment object
# read in the skin seurat object
file <- "skin.rda"
path <- paste(dir.in, file, sep="")
load(path)

# read in the barcode mapping
file <- 'celltype_skin.txt'
path <- paste(dir.in, file, sep="")
ct <- read.table(path, header = T, sep='\t')
dim(ct)
# [1] 34774     3

# rename cells
head(x = colnames(x = skin))
# [1] "R1.01.R2.01.R3.06.P1.55" "R1.01.R2.03.R3.68.P1.55"
# [3] "R1.01.R2.05.R3.15.P1.53" "R1.01.R2.05.R3.40.P1.55"
# [5] "R1.01.R2.05.R3.49.P1.55" "R1.01.R2.06.R3.14.P1.55"
length(colnames(skin))
# [1] 29308
table(!is.na(match(colnames(skin), ct$rna.bc)))
# 
# TRUE
# 29308
atac.bc <- ct$atac.bc[match(colnames(skin), ct$rna.bc)]
# atac.bc <- ct$atac.bc[!is.na(match(ct$rna.bc, colnames(skin)))]
skin <- RenameCells(object = skin, new.names = atac.bc)
head(x = colnames(x = skin))
# [1] "R1.01.R2.01.R3.06.P1.07" "R1.01.R2.03.R3.68.P1.07"
# [3] "R1.01.R2.05.R3.15.P1.05" "R1.01.R2.05.R3.40.P1.07"
# [5] "R1.01.R2.05.R3.49.P1.07" "R1.01.R2.06.R3.14.P1.07"

# create a fragment object
file <- "skin.atac.fragments.tsv.gz"
path <- paste(dir.in, file, sep="")
frags <- CreateFragmentObject(path)
DefaultAssay(skin) <- "ATAC"
# data <- GetAssayData(skin, slot = "data")
# chrom_assay <- skin[["ATAC"]] 
# Fragments(chrom_assay) 
Fragments(skin) <- frags
class(Fragments(skin)[[1]])
# [1] "Fragment"
# attr(,"package")
# [1] "Signac"

file <- "skin.frag.rds"
path <- paste(dir.in, file, sep="")
if (!file.exists(path)) saveRDS(skin.frag, file=path)

# # sanity check (see Figure 2H in Ma et al. (2020))
# # Gapdh
# fig <- paste("covplot_Gapdh.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Gapdh",
#   annotation = T,
#   peaks = T,
#   extend.downstream = 10000
# )
# dev.off()
# 
# # Krt5
# fig <- paste("covplot_Krt5.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Krt5",
#   annotation = T,
#   peaks = T,
#   extend.downstream = 10000
# )
# dev.off()
# 
# # Krt15
# fig <- paste("covplot_Krt15.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Krt15",
#   annotation = T,
#   peaks = T,
#   extend.upstream = 10000
# )
# dev.off()
# 
# # Lgr5
# fig <- paste("covplot_Lgr5.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Lgr5",
#   annotation = T,
#   peaks = T,
#   extend.upstream = 10000
# )
# dev.off()



library(Seurat)
library(Signac)

skin=readRDS('skin.frag.rds')

## run chromVAR
DefaultAssay(skin) <- "ATAC"
levels(skin)
#  [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
# [16] "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27"
celltype.names <- levels(skin)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(skin) <- "ATAC"
# https://satijalab.org/signac/articles/motif_vignette.html
# get a list of motif position frequency matrices from the jaspar database
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=10090
# 10090 is the species ID for mouse
# 9606 is the species ID for human
# https://satijalab.org/signac/articles/motif_vignette.html uses 9606 for mouse
# and this link gives the justification: https://github.com/timoast/signac/issues/58 

pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(skin), 
                                  pwm = pwm_set, 
                                  genome = 'mm10', 
                                  use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
skin <- SetAssayData(skin, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
skin <- RunChromVAR( # note that this step can take 30-60 minutes
  object = skin,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(skin, file='skin.chromvar.rds')


