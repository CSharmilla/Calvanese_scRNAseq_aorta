library(here)
library(scater)
library(scran)
library(patchwork)
library(dplyr)


mnn.out.clus <- readRDS(here("output", "mnn.out.clus.RDS"))
k <- "cluster_k10" 

plotTSNE(mnn.out.clus, colour_by = k, point_size=1, text_by = k, text_size =5,text_colour = "red")


# creating markers list
markers_list <- list(endothelial = c("CDH5", "CD34", "KDR","TIE1", "APLNR", "NRP2", "GJA5", "CXCR4", "NR2F2"),
                     hsc = c("RUNX1", "SPN", "PTPRC", "MLLT3", "MYB", "GFI1",
                                      "GATA2","HOXA9", "MECOM", "HLF", "SPINK2", 
                                      "CD34", "THY1", "ACE"),
                     hematopoietic = c("RUNX1", "SPN", "PTPRC"),
                     stroma = c("COL1A1", "PDGFRA", "CXCL12", "POSTN", "PTN"),
                     erythroid = c("HBE1", "HBZ", "GYPA"), 
                     proliferative = c("MKI67", "TOP2A", "AURKB"), 
                     epithelial = c("EPCAM", "AFP"), 
                     fibroblasts = c("PAX1", "SOX9", "HAND1", "NKX2-3", "DCN", "ALDH1A2", 
                                     "COL14A1", "CRABP1", "LUM", "PAX3", "TBX5", "HAND2"),
                     EHT_genes = c("JUN",
                                   "FOSB", 
                                   "NR4A1",
                                   "KLF2",
                                   "REL", 
                                   "HES1",
                                   "EGR3",
                                   "IRF1",
                                   "ATF3",
                                   "BACH2",
                                   "RARG",
                                   "SOX17",
                                   "FOXC2",
                                   "EBF1",
                                   "RORA",
                                   "SNAI1",
                                   "HOXA3",
                                   "ZBTB16",
                                   "YAP1",
                                   "WWTR1", "CSMD2", "RAB27B"
                     ),
                     Venous_EC = c("EPHB4","LEFTY1", "LEFTY2", "FLT4", "NRP2", "NR2F2", "EMCN", "APLNR"))

# plot markers TSNE
plot_list <- lapply(markers_list[["endothelial"]],
                    function(x)plotTSNE(mnn.out.clus, colour_by = x, text_by = k, text_colour = "red") + ggtitle(x))

patchwork::wrap_plots(plot_list) +  patchwork::plot_annotation(title = "endothelial")



plot_list <- lapply(markers_list[["hsc"]],
                    function(x)plotTSNE(mnn.out.clus, colour_by = x, text_by = k, text_colour = "red") + ggtitle(x))

patchwork::wrap_plots(plot_list) +  patchwork::plot_annotation(title = "hsc")

plotTSNE(mnn.out.clus, colour_by = "label", point_size=1, text_by = k, text_size =5)

# markers dot plots
lapply(names(markers_list), function(x){plotDots(mnn.out.clus, features = markers_list[[x]], group = k) + ggtitle(x)})




# add annotation to sce - reading in ann file
annotation.1 <- read.csv(here("data", "annotation1.csv"))
ann.vec <- annotation.1$celltype
names(ann.vec) <- annotation.1$cluster_k10

## adding to coldata in a new sce
# creating new
sce.ann.1 <- mnn.out.clus
coldata <- colData(sce.ann.1) %>% as.data.frame %>% mutate(annotation_1=recode(cluster_k10, !!!ann.vec)) %>% as("DataFrame")
# sanity check
all(colnames(sce.ann.1) == rownames(coldata))
colData(sce.ann.1) <- coldata

# plotting
plotTSNE(sce.ann.1, colour_by = "annotation_1", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")

plotTSNE(sce.ann.1, colour_by = "SPN", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "PTPRC", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "NR2F2", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "APLNR", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "GJA4", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "GJA5", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "ETS1", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "IL33", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "PLXNA4", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")
plotTSNE(sce.ann.1, colour_by = "RUNX1", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")

saveRDS(sce.ann.1, here("output", "sce.ann.1.RDS"))


# creating some dot plots of genes of interest
plotDots(sce.ann.1, features = markers_list[["EHT_genes"]], group = "annotation_1") + ggtitle("EHT_genes") + theme(axis.text.x = element_text(angle = 45))

plotTSNE(sce.ann.1, colour_by = "CSMD2", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")

plotTSNE(sce.ann.1, colour_by = "RUNX1", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")

# dot plots for sub clustered sce
sce.ann1.sub <- readRDS(here("output", "sce.ann1.sub.RDS"))
lapply(names(markers_list), function(x){plotDots(sce.ann1.sub, features = markers_list[[x]], group = "annotation_1") + ggtitle(x) + theme(axis.text.x = element_text(angle = 45))}) 


