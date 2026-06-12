library(here)
library(scater)
library(slingshot)


subcluster_list <- readRDS(here("output", "subcluster_list.RDS"))

k <- "cluster_k10" 
plotTSNE(sce.ann1.aec.sub, colour_by = k, point_size=1, text_by = k, text_size =5)

# with slingshot

# fit a single principal curve to yields a pseudotime ordering of cells based on their relative positions when projected onto the curve
# as issue as https://github.com/kstreet13/slingshot/issues/87
# throws an error when using all 50 dimensions. 1:49 works fine which is very strange. 

if (!file.exists(here("output", "sce.sling.RDS"))){
  celltype <- "hema_2"
  cluster <- "subcluster"
  start.clus <- "hema_2.5"
  end.clus <- NULL
  sce.sling <- slingshot(subcluster_list[[celltype]], 
                         reducedDim=reducedDims(subcluster_list[[celltype]])$corrected[,1:49], 
                         cluster = cluster,
                         start.clus=start.clus,
                         end.clus=end.clus)
                         
  
  # investigating
  # getting an error when using mnn corrected reddim
  # table(duplicated(reducedDims(sce.ann1.aec.sub)$corrected))
  # ncol(reducedDims(sce.ann1.aec.sub)$corrected)
  # # check whether each cluster's covariance matrix is invertible 
  # for(cl in unique(sce.ann1.aec.sub$comb_ann1)){
  #   solve(cov(reducedDims(sce.ann1.aec.sub[,sce.ann1.aec.sub$comb_ann1 == cl])$corrected[,1:49]))
  # }
  # 
  # cl <- "endo_1_NA"
  # x <- reducedDims(sce.ann1.aec.sub[,sce.ann1.aec.sub$comb_ann1 == cl])$corrected
  # heatmap(cor(x), symm = TRUE, Rowv=NA, Colv = NA)
  
  
  saveRDS(sce.sling, here("output", "sce.sling.RDS"))
  
}else{
  sce.sling <- readRDS(here("output", "sce.sling.RDS"))
  
}


# looking at lineages
slingLineages(sce.sling) 


# creating matrix of pseudotimes
pseudo.paths <- slingPseudotime(sce.sling)
head(pseudo.paths) 

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have similar pseudo-time values in 
# all paths anyway, so taking the rowMeans is not particularly controversial.
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Need to loop over the paths and add each one separately.
gg <- plotTSNE(sce.sling, colour_by=I(shared.pseudo), text_by = cluster)
embedded <- embedCurves(sce.sling, "TSNE")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=TSNE1, y=TSNE2), size=1.2)
  print(gg)
}

gg

# using monocle

