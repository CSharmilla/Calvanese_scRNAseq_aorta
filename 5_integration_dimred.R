library(here)
library(scran)
library(scater)
library(batchelor)

if (!file.exists(here("output", "mnn.out.RDS"))){
  
  
  all.sce.norm <- readRDS(here("output", "all.sce.norm.RDS"))
  all.var <- readRDS(here("output", "all.var.RDS"))
  all.hvgs <- readRDS(here("output", "all.hvgs.RDS"))
  
  
  # 
  # ### dimensionality reduction
  # ### PCA
  # ### double checking
  # all(names(all.sce.norm) == names(all.hvgs))
  # 
  # ### using technical variation to find threshold for useful PCs
  # set.seed(111001001)
  # all.sce.norm <- mapply(FUN=denoisePCA, x=all.sce.norm, technical=all.var, subset.row = all.hvgs, name = "pca_denoise")
  # plot_list <- list()
  # for (n in names(all.sce.norm)){
  #   plot_list[[n]] <- plotReducedDim(all.sce.norm[[n]], colour_by="sum", dimred = "pca_denoise") + ggtitle(paste0("PCA denoise ", n))
  #   str(reducedDim(all.sce.norm[[n]], "pca_denoise")) # compactly display structure
  # }
  # patchwork::wrap_plots(plot_list)
  # 
  # ### plotting TSNE
  # set.seed(00101001101)
  # ### When dimred is specified, no additional feature selection or standardization is performed. This means that any settings of ntop, subset_row and scale are ignored.
  # all.sce.norm <- lapply(all.sce.norm, runTSNE, dimred="pca_denoise", perplexity=90)
  # plots <- list()
  # for (i in names(all.sce.norm)){
  #   plots[[i]] <- plotReducedDim(all.sce.norm[[i]], dimred ="TSNE", colour_by = "sum") + ggtitle(i)
  # }
  # patchwork::wrap_plots(plots)
  # 
  # saveRDS(all.sce.norm, file = here("output", "sce.norm.dimred.RDS"))
  
  
  # integration
  # common universe of features (genes)
  universe <- Reduce(intersect, lapply(all.sce.norm, rownames))
  length(universe)
  
  
  # subsetting the sce objects to only include common genes
  all.sce.com <- lapply(all.sce.norm, "[", i=universe,)
  all.var.com <- lapply(all.var, "[", i=universe,)
  
  
  
  # tidying up a bit
  for (n in names(all.sce.com)){
    rowData(all.sce.com[[n]]) <- rowData(all.sce.com[[n]])[,c("ID", "Symbol", "new_names")]
    colData(all.sce.com[[n]]) <- colData(all.sce.com[[n]])[,c("label", "sum", "subsets_Mito_percent",
                                                              "subsets_Ribo_percent")]
    reducedDims(all.sce.com[[n]]) <- NULL
  }
  
  
  # rescale each batch to adjust for differences in sequencing depth between batches
  renorm.sce <- do.call(multiBatchNorm, all.sce.com)
  
  
  # perform feature selection by averaging the variance components across all batches with the combineVar() function
  combined.var <- do.call(combineVar, all.var.com)
  combined.hvgs <- getTopHVGs(combined.var, prop = 0.15) #combined.var$bio > 0
  
  # unlikely but check anyway
  for (i in names(renorm.sce)){
    if(any(duplicated(colnames(renorm.sce[[i]])))){
      print("duplicated barcodes")
    } else{
      print("no duplicated barcodes")
    }
  }
  
  ### looking at data before correction
  # Synchronizing the metadata for cbind()ing.
  uncorrected <- renorm.sce
  
  # lapply() always returns a list, whereas sapply() tries to simplify the 
  # result into a vector or matrix
  # check if all rowdatas are equal
  #all(lapply(lapply(uncorrected, rowData), identical, rowData(uncorrected[[1]])))
  
  if (all(sapply(lapply(uncorrected, rowData), identical, rowData(uncorrected[[1]])))) {
    uncorrected <- do.call(cbind, uncorrected)
  } else{print("rowdatas are not in sync")}
  
  # fix duplicated barcodes in combined data
  if(any(duplicated(colnames(uncorrected)))){
    colnames(uncorrected) <- paste0(colnames(uncorrected), sep="_", uncorrected$label)
  } else{
    print("no duplicated barcodes")
  }

  coldata <- data.frame(colData(uncorrected))
  rowdata <- data.frame(rowData(uncorrected))
  
  uncorrected$label <- factor(uncorrected$label)
  
  # Using RandomParam() as it is more efficient for file-backed matrices.
  set.seed(0010101010)
  uncorrected <- runPCA(uncorrected, subset_row=combined.hvgs,
                        BSPARAM=BiocSingular::RandomParam())
  
  uncorrected <- denoisePCA(uncorrected, technical=combined.var, subset.row = combined.hvgs, name = "pca_denoise")
  
  # plotting tsne to see how the cells cluster
  set.seed(1111001)
  uncorrected <- runTSNE(uncorrected, dimred="pca_denoise")
  uncorrected <- runUMAP(uncorrected, dimred="pca_denoise", min_dist = 0.08)
  
  #tsne
  plotTSNE(uncorrected, colour_by="label", point_alpha = 0.2) + ggtitle("before correction") # alot of batch effect
  #umap
  plotUMAP(uncorrected, colour_by="label", point_alpha = 0.2) + ggtitle("before correction") 
  
  # correcting data
  set.seed(1000101001)
  
  # At the highest level, we have the correctExperiments() function that wraps the batchCorrect() function. 
  # This is intended for the specific case where correction is applied on data in existing SingleCellExperiment objects and we wish to 
  # preserve the pre-existing data and metadata from those objects in the corrected output.
  
  # if a merged dataset is submitted then do specify batch 
  mnn.out <- correctExperiments(uncorrected, batch = uncorrected$label, subset.row=combined.hvgs, 
                                correct.all=TRUE, PARAM = FastMnnParam(d=50, prop.k=0.05, BSPARAM=BiocSingular::RandomParam(),
                                                                       merge.order=list(
                                                                         "CS14_wk4", 
                                                                         list("CS15_wk5", "CS15_wk5.5"))))
                                                                         #"week6")))
  mnn.out
  
  # plotting 
  set.seed(0010101010)
  mnn.out <- runTSNE(mnn.out, dimred="corrected")
  mnn.out <- runUMAP(mnn.out, dimred="corrected", min_dist = 0.08) # minimum distance between points in low-dimensional space
  # controls how tightly UMAP clumps points together, with low values leading to more tightly packed embeddings. Larger values 
  # of min_dist will make UMAP pack points together more loosely, focusing instead on the preservation of the broad 
  # topological structure.

  # plotting
  plotUMAP(mnn.out, colour_by="label", point_alpha = 0.3) + ggtitle("UMAP after correction")
  
  plotTSNE(mnn.out, colour_by="label", point_alpha = 0.3) + ggtitle("TSNE after correction")
  plotUMAP(mnn.out, colour_by="sum", point_alpha = 0.3) + ggtitle("UMAP after correction")
  plotTSNE(mnn.out, colour_by="subsets_Ribo_percent", point_alpha = 0.3) + ggtitle("TSNE after correction")
  
  

  
  # saving
  saveRDS(mnn.out, file = here("output", "mnn.out.RDS"))
  saveRDS(uncorrected, file = here("output", "uncorrected.RDS"))
  
} else {
  mnn.out <- readRDS(here("output", "mnn.out.RDS"))
}



