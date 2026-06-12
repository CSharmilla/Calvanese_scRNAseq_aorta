library(here)
library(MatrixGenerics)
library(scran)
library(scater)

sce <- readRDS(here("output", "sce.RDS"))

matLog2 <- function(spmat, scale = FALSE, scaleFactor = 1e6) {
  
  
  if (scale == TRUE) {
    spmat <- t( t(spmat) / colSums(spmat)) * scaleFactor
  }
  
  if (is(spmat, "sparseMatrix")) {
    matsum <- summary(spmat)
    
    logx <- log2(matsum$x + 1)
    
    logmat <- sparseMatrix(i = matsum$i, j = matsum$j,
                           x = logx, dims = dim(spmat),
                           dimnames = dimnames(spmat))
  } else {
    logmat <- log2(spmat + 1)
  }
  
  
  return(logmat)
  
}


set.seed(100)
all.clust.sce <- lapply(sce, quickCluster)
all.sce.norm <- list()
plot <- hist <- list()
for(n in names(sce)){
  all.sce.norm[[n]] <- computeSumFactors(sce[[n]], cluster=all.clust.sce[[n]], min.mean=0.1)
  print(summary(sizeFactors(all.sce.norm[[n]])))
  
  all.sce.norm[[n]] <- logNormCounts(all.sce.norm[[n]])
  
  # plotting sf normalised data
  hist(log10(all.sce.norm[[n]]$sizeFactor), xlab="Log10[Size factor]", col='grey80')
  
  ## plotting deconvolution size factors and size factors 
  sf <- librarySizeFactors(all.sce.norm[[n]])
  deconv <- calculateSumFactors(all.sce.norm[[n]], cluster=all.clust.sce[[n]])
  plot(sf, deconv, xlab="Library size factor",
       ylab="Deconvolution size factor", log='xy', pch=16, col=all.clust.sce[[n]])
  abline(a=0, b=1, col="red")
  # log transform raw counts also
  assay(all.sce.norm[[n]], "logcounts_raw") <- matLog2(counts(all.sce.norm[[n]]))            
  
}
all.sce.norm



### did normalisation work?

### looking at variance explained

### before normalisation
### Compute, for each gene, the percentage of variance that is explained by one or more variables of interest.
var <- lapply(all.sce.norm, scater::getVarianceExplained, exprs_values = "logcounts_raw", variables = c("sum", "detected"))

### after normalisation
var_after <- lapply(all.sce.norm, getVarianceExplained, variables = c("sum", "detected"))

# Plot explanatory variables ordered by percentage of variance explained - before and after normalisation
p <- list()
for (x in names(var_after)){
  p[[1]] <- plotExplanatoryVariables(var[[x]]) + ggtitle(paste0("before normalisation for ",x))
  p[[2]] <- plotExplanatoryVariables(var_after[[x]]) + ggtitle(paste0("after normalisation for ",x))
  print(patchwork::wrap_plots(p))
}

### now looking at pca
### before
raw_pca <- lapply(all.sce.norm, runPCA, exprs_values = "logcounts_raw")
### after
all.sce.norm <- lapply(all.sce.norm, runPCA, name = "pca")
plots <- list()
for (n in names(var_after)){
  plots[[1]] <- plotPCA(raw_pca[[n]], colour_by="sum") + ggtitle(paste0("Before normalisation ", n))
  plots[[2]] <- plotReducedDim(all.sce.norm[[n]], colour_by="sum", dimred = "pca") + ggtitle(paste0("after normalisation ", n))
  print(patchwork::wrap_plots(plots))
}


saveRDS(all.sce.norm, file = here("output", "all.sce.norm.RDS"))
