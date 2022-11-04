palettecb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


library(gridExtra)
library(gtable)
library(grid)
library(rdist)
library(FNN)
library(tidyverse)
library(Hmisc)

### Note that throughout titles of tables/graphs should be edited depending on dataset used


# For this data, MCCA and CPCA do not have the corresponding UMAP slot.
obj <- RunUMAP(obj, dims = 1:30, reduction='MCCA', 
               reduction.name='MCCA.umap', reduction.key = 'MCCAumap_', verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, reduction='CPCA', 
               reduction.name='CPCA.umap', reduction.key = 'CPCAumap_', verbose = FALSE)

# The parts below should apply to all datasets, so long as they are processed in the same format

#### X - Matrix in ATAC space (ex: cell by motif matrix) or matrix of distances
#### neighbor_indices are n x 2 matrix of with first column indicating cell and second indicating index of that cells nearest neighbor 
#### in RNA space
####
#### Note that for this to work cells must be in the same order for ATAC and RNA data
#### dist is parameter specifying whether X is either matrix in ATAC space or matrix of distances (defaults to false)
calcFOSCTTNN <- function(X, neighbor_indices, dist = FALSE) {
  X_distances <- X
  if(!dist) {
    X_distances <- as.matrix(rdist(X))
  }
  ## nearest neighbors from RNA 
  matching_dist <- X_distances[neighbor_indices]
  FOSCTTNN <- (rowSums(X_distances <= matching_dist) - 1) / (length(matching_dist) - 1)
  return(FOSCTTNN)
}


###############################
###############################
#### FOSCTTNN
###############################
###############################

DefaultAssay(obj) <- 'RNA' ## Change to SCT for skin dataset
obj <- FindNeighbors(obj, features = VariableFeatures(object = obj), return.neighbor = TRUE)
obj_rna_indices <- obj[["SCT.nn"]]@nn.idx[, 1:2]

obj_df <- data.frame(Method = c(rep("Peak_LSI", ncol(obj)), 
                                rep("Peak_LDA", ncol(obj)),
                                rep("Motif", ncol(obj)), 
                                rep("GeneActivity", ncol(obj)),
                                rep("ConsensusPCA", ncol(obj)), 
                                rep("MultiCCA", ncol(obj)),
                                rep("WNN", ncol(obj))), 
                     FOSCTTNN = c(calcFOSCTTNN(obj@reductions$lsi.umap@cell.embeddings, obj_rna_indices),
                                 calcFOSCTTNN(obj@reductions$lda.umap@cell.embeddings, obj_rna_indices),
                                 calcFOSCTTNN(obj@reductions$motif.umap@cell.embeddings, obj_rna_indices),
                                 calcFOSCTTNN(obj@reductions$activity.umap@cell.embeddings, obj_rna_indices),
                                 calcFOSCTTNN(obj@reductions$CPCA.umap@cell.embeddings, obj_rna_indices),
                                 calcFOSCTTNN(obj@reductions$MCCA.umap@cell.embeddings, obj_rna_indices),
                                 calcFOSCTTNN(obj@reductions$wnn.umap@cell.embeddings, obj_rna_indices)))

obj_df$Method=factor(obj_df$Method, levels=unique(obj_df$Method)) # Change the order of the methods

ggplot(obj_df, aes(x = Method, y = FOSCTTNN, fill = Method)) + geom_boxplot() + ggtitle("Dataset FOSCTNN") + 
  scale_fill_manual(values = palettecb) + theme_gray(base_size = 15) + theme(axis.text.x = element_text(size = 9.5))

ggsave("data_FOSCTTNN_boxplot.pdf", device = "pdf")

write.csv(obj_df, "data_FOSCTTNN.csv")


################## 
### Table Creation
##################
my_obj_df <- obj_df %>% group_by(Method) %>% dplyr::summarize(
  q1 = quantile(FOSCTTNN, 0.25),
  median = median(FOSCTTNN),
  GMD = GiniMd(FOSCTTNN),
  mean = mean(FOSCTTNN),
  sd = sd(FOSCTTNN),
  q3 = quantile(FOSCTTNN, 0.75)) 

my_obj_df[, -1] <- my_obj_df[, -1] %>% round(3)


my_obj_df <- tableGrob(my_obj_df, rows = NULL)
title <- textGrob("ata FOSCTTNN", gp = gpar(fontsize = 18))
padding <- unit(.75,"line")
my_obj_df <- gtable_add_rows(
  my_obj_df, heights = grobHeight(title) + padding, pos = 0
)
my_obj_df <- gtable_add_grob(
  my_obj_df, list(title),
  t = 1, l = 1, r = ncol(my_obj_df)
)

#Export to pdf
pdf('data_FOSCTTNN_tbl.pdf',width = 10)
grid.newpage()
grid.draw(my_obj_df)
dev.off()

###############################
###############################
#### Agreement
###############################
###############################
# X is embedding from various methods from ATAC
# obj_rna_indices: the nearest neighbor indices for RNA (not including the cell itself)
# k.param: how many neighbors?
calcAgreement <- function(X, obj_rna_indices, k.param) {
  obj_atac_indices <- get.knn(X, k= k.param)$nn.index
  agreement=rep(NA, nrow(obj_rna_indices))
  for(i in 1:length(agreement)){
    agreement[i]=length(intersect(obj_rna_indices[i,], obj_atac_indices[i,]))/ncol(obj_rna_indices)
  }
  return(mean(agreement))
}

obj <- FindNeighbors(obj, features = VariableFeatures(object = obj), return.neighbor = T)

obj_rna_indices <- obj[["RNA.nn"]]@nn.idx[,-1] ## Change to SCT.nn for skin dataset

obj_df <- data.frame( k = c(50, 100, 150, 200),
                     
                     Peak_LSI = c(calcAgreement(obj@reductions$lsi.umap@cell.embeddings, obj_rna_indices, k.param=50),
                             calcAgreement(obj@reductions$lsi.umap@cell.embeddings, obj_rna_indices, k.param=100),
                             calcAgreement(obj@reductions$lsi.umap@cell.embeddings, obj_rna_indices, k.param=150),
                             calcAgreement(obj@reductions$lsi.umap@cell.embeddings, obj_rna_indices, k.param=200)),
                     
                     Peak_LDA = c(calcAgreement(obj@reductions$lda.umap@cell.embeddings, obj_rna_indices, k.param=50),
                              calcAgreement(obj@reductions$lda.umap@cell.embeddings, obj_rna_indices, k.param=100),
                              calcAgreement(obj@reductions$lda.umap@cell.embeddings, obj_rna_indices, k.param=150),
                              calcAgreement(obj@reductions$lda.umap@cell.embeddings, obj_rna_indices, k.param=200)),
                     
                     Motif = c(calcAgreement(obj@reductions$motif.umap@cell.embeddings, obj_rna_indices, k.param=50),
                               calcAgreement(obj@reductions$motif.umap@cell.embeddings, obj_rna_indices, k.param=100),
                               calcAgreement(obj@reductions$motif.umap@cell.embeddings, obj_rna_indices, k.param=150),
                               calcAgreement(obj@reductions$motif.umap@cell.embeddings, obj_rna_indices, k.param=200)),
                     
                     GeneActivity = c(calcAgreement(obj@reductions$activity.umap@cell.embeddings, obj_rna_indices, k.param=50),
                                       calcAgreement(obj@reductions$activity.umap@cell.embeddings, obj_rna_indices, k.param=100),
                                       calcAgreement(obj@reductions$activity.umap@cell.embeddings, obj_rna_indices, k.param=150),
                                       calcAgreement(obj@reductions$activity.umap@cell.embeddings, obj_rna_indices, k.param=200)),
                     

                     
                     ConsensusPCA = c(calcAgreement(obj@reductions$CPCA.umap@cell.embeddings, obj_rna_indices, k.param=50),
                              calcAgreement(obj@reductions$CPCA.umap@cell.embeddings, obj_rna_indices, k.param=100),
                              calcAgreement(obj@reductions$CPCA.umap@cell.embeddings, obj_rna_indices, k.param=150),
                              calcAgreement(obj@reductions$CPCA.umap@cell.embeddings, obj_rna_indices, k.param=200)), 
                     
                     MultiCCA = c(calcAgreement(obj@reductions$MCCA.umap@cell.embeddings, obj_rna_indices, k.param=50),
                              calcAgreement(obj@reductions$MCCA.umap@cell.embeddings, obj_rna_indices, k.param=100),
                              calcAgreement(obj@reductions$MCCA.umap@cell.embeddings, obj_rna_indices, k.param=150),
                              calcAgreement(obj@reductions$MCCA.umap@cell.embeddings, obj_rna_indices, k.param=200)), 
                     
                     WNN = c(calcAgreement(obj@reductions$wnn.umap@cell.embeddings, obj_rna_indices, k.param=50),
                             calcAgreement(obj@reductions$wnn.umap@cell.embeddings, obj_rna_indices, k.param=100),
                             calcAgreement(obj@reductions$wnn.umap@cell.embeddings, obj_rna_indices, k.param=150),
                             calcAgreement(obj@reductions$wnn.umap@cell.embeddings, obj_rna_indices, k.param=200)))

write.csv(obj_df, "data_agreement.csv")

library(reshape2)
obj_long <- melt(obj_df, id = "k")
names(obj_long) <- c("k", "Method", "Agreement")
ggplot(data = obj_long, aes(x = k, y = Agreement, color = Method)) + geom_line() + ggtitle("Data Agreement") + 
  scale_color_manual(values = palette) + theme_gray(base_size = 15)
  
ggsave("data_agreement.pdf", device = "pdf")


my_obj_df <- tableGrob(obj_df %>% round(3), rows = NULL)
title <- textGrob("Data Agreement", gp = gpar(fontsize = 18))
padding <- unit(.75,"line")
my_obj_df <- gtable_add_rows(
  my_obj_df, heights = grobHeight(title) + padding, pos = 0
)
my_obj_df <- gtable_add_grob(
  my_obj_df, list(title),
  t = 1, l = 1, r = ncol(my_obj_df)
)

#Export to pdf
pdf('data_agreement_tbl.pdf',width = 10)
grid.newpage()
grid.draw(my_obj_df)
dev.off()


