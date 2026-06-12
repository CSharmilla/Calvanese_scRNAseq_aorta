library(here)
library(scater)
library(AnnotationHub)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(scuttle)
library(patchwork)



sample.list <- readRDS(here("output", "sample.list"))

# converting to singlecellexperiment objects
sce.list <- lapply(sample.list, function(x){SingleCellExperiment(assays = list(counts = x))})


# adding rowdata, converting ensembl ids to gene names
for (i in names(sce.list)){
  print(i)
  rowData(sce.list[[i]])["ID"] <- genes <- rownames(sce.list[[i]])
  # useful https://jorainer.github.io/ensembldb/reference/EnsDb-AnnotationDbi.html
  gene_annot <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", 
                       column = c("SYMBOL"), multiVals = "first")

  if (all(rownames(gene_annot) %in% genes)){
    rowData(sce.list[[i]]) <- merge(rowData(sce.list[[i]]), gene_annot, by = 0, sort = FALSE, all.x = TRUE)
    rowData(sce.list[[i]]) <- rowData(sce.list[[i]])[-1]
    colnames(rowData(sce.list[[i]])) <- c("ID", "Symbol")
  }else{
    print("all genes were not annotated")
  }
  
  # sanity check
  if (all(rownames(sce.list[[i]]) == rowData(sce.list[[i]])[,"ID"])){
    new_names <- uniquifyFeatureNames(rownames(sce.list[[i]]), rowData(sce.list[[i]])[, "Symbol"])
    rownames(sce.list[[i]]) <- new_names
    rowData(sce.list[[i]])$new_names <- new_names
  }else{
    print("rowdata and rownames are out of sync")
  }
}
#a <- data.frame(rowData(sce.list[[3]]))
# QC will be done for each batch separately

### adding tags for each sample to keep them separate
for (i in names(sce.list)){
  colData(sce.list[[i]])$label <- i
}

### find undetected genes
#detected_genes <- lapply(sce.list, function(x) rowSums(counts(x))>0)
#lapply(detected_genes, table)

# making sure we keep genes that are detected in at least one dataset
# initialise
gene_names <- rownames(sce.list[[1]])

# Initialising a logical vector of FALSE
genes_detected <- rep(FALSE, length(gene_names))
names(genes_detected) <- gene_names

# Loop over each sce object
for (sce in sce.list) {
  # Check which genes are detected 
  detected_in_sce <- rowSums(counts(sce))>0 
  
  # Update the overall detection vector
  genes_detected <- genes_detected | detected_in_sce
}


### filter out undetected genes
det_sce <- list()
for (i in names(sce.list)){
  det_sce[[i]] <- sce.list[[i]][genes_detected, ]
}

### adding cell - level QC metrics
det_sce_qc <-  lapply(det_sce, function(a){addPerCellQC(a, subsets=list(Mito=grep("MT-", rowData(a)$Symbol), Ribo=grep("^RP[SL]", rowData(a)$Symbol)))})
### summary of QC results
lapply(det_sce_qc, function(x){summary(colData(x)$sum)})
lapply(det_sce_qc, function(x){summary(colData(x)$detected)})
lapply(det_sce_qc, function(x){summary(colData(x)$subsets_Mito_percent)})
lapply(det_sce_qc, function(x){summary(colData(x)$subsets_Ribo_percent)})

### adaptive thresholds
### library size
low_lib_size <- lapply(det_sce_qc, function(a){isOutlier(a$sum, log = TRUE, type = "lower")})
lapply(low_lib_size, table) 
lapply(low_lib_size, function(a){attr(a, "thresholds")})
for (i in names(det_sce_qc)){
  colData(det_sce_qc[[i]])$low_lib_size <- low_lib_size[[i]]
  print(plotColData(det_sce_qc[[i]],
                    x="label",
                    y="sum",
                    colour_by = "low_lib_size") +
          scale_y_log10() +
          labs(y = "Total count", title = paste0("Total count for ", i)) +
          guides(colour=guide_legend(title="Discarded")))
}

### number of genes
low_n_features <- lapply(det_sce_qc, function(a){isOutlier(a$detected, log = TRUE, type = "lower")})
lapply(low_n_features, table) 
lapply(low_n_features, function(a){attr(a, "thresholds")})
for (i in names(det_sce_qc)){
  colData(det_sce_qc[[i]])$low_n_features <- low_n_features[[i]]
  print(plotColData(det_sce_qc[[i]],
                    x="label",
                    y="detected",
                    colour_by = "low_n_features") +
          scale_y_log10() +
          labs(y = "Genes detected", title = paste0("Genes detected for ", i)) +
          guides(colour=guide_legend(title="Discarded")))
}

### mito genes
high_Mito_percent <- lapply(det_sce_qc, function(a){isOutlier(a$subsets_Mito_percent, type = "higher")})
lapply(high_Mito_percent, table) 
lapply(high_Mito_percent, function(a){attr(a, "thresholds")})
for (i in names(det_sce_qc)){
  colData(det_sce_qc[[i]])$high_Mito_percent <- high_Mito_percent[[i]]
  print(plotColData(det_sce_qc[[i]],
                    x="label",
                    y="subsets_Mito_percent",
                    colour_by = "high_Mito_percent") +
          labs(y = "Percentage mitochondrial UMIs", title = paste0("Mitochondrial UMIs ", i)) +
          guides(colour=guide_legend(title="Discarded")))
}



### summary of discarded cells
for (i in names(det_sce_qc)){
  print(data.frame(`Library Size` = sum(low_lib_size[[i]]),
                   `Genes detected` = sum(low_n_features[[i]]),
                   `Mitochondrial UMIs` = sum(high_Mito_percent[[i]]),
                   Total = sum(low_lib_size[[i]] | low_n_features[[i]] | high_Mito_percent[[i]])))
}


### adding a discard column to each dataset individually
for (i in names(det_sce_qc)){
  colData(det_sce_qc[[i]])$discard <- with(colData(det_sce_qc[[i]]), ifelse(low_lib_size == TRUE | low_n_features == TRUE | (low_lib_size == TRUE & high_Mito_percent == TRUE), TRUE, FALSE))
  
  
  plots <- list()
  # plots[[1]] <- plotColData(det_sce_qc[[i]],
  #                           x="sum",
  #                           y="subsets_Mito_percent",
  #                           colour_by="discard")
  # plots[[2]] <- plotColData(det_sce_qc[[i]],
  #                           x="detected",
  #                           y="subsets_Mito_percent",
  #                           colour_by="discard")
  plots[[3]] <- plotColData(det_sce_qc[[i]],
                            x="sum",
                            y="detected",
                            colour_by="discard")
  #print(patchwork::wrap_plots(plots) + patchwork::plot_annotation(title = paste0("plots for ", i)))
}

### filter out poor quality cells
sce_fil <- lapply(det_sce_qc, function(a){a[,!a$discard]})


### recalculating QC metrics
### look at coldata
lapply(sce_fil, colData)
sce <- sce_fil
for (i in names(sce)){
  colData(sce[[i]]) <- colData(sce[[i]])[,1:3]
}
lapply(sce, colData)
lapply(sce, rowData)
sce <-  lapply(sce, function(a){addPerCellQC(a, subsets=list(Mito=grep("MT-", rowData(a)$Symbol), Ribo=grep("^RP[SL]", rowData(a)$Symbol)))})
lapply(sce, colData)

saveRDS(sce, file = here("output", "sce.RDS"))


