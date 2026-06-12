library(here)
library(EnhancedVolcano)
library(scran)



sce.ann.1 <- readRDS(here("output", "sce.ann.1.RDS"))



summed <- aggregateAcrossCells(sce.ann.1, id=colData(sce.ann.1)[,c("annotation_1", "label")])

comparisons <- list(hsc_vs_aec = list(comp = c("HSC", "aEC"),
                                      coef = "annotation_1HSC"),
                    hsc_vs_hema1 = list(comp = c("HSC", "hema_1"),
                                        coef = "annotation_1HSC"), 
                    hema1_vs_hema2 = list(comp = c("hema_1", "hema_2"),
                                          coef = "annotation_1hema_1"))

for (i in names(comparisons)){
  #i <- "hsc_vs_hema1"
  summed.sub <- summed[,summed$annotation_1 %in% comparisons[[i]][["comp"]]]
  summed.sub$annotation_1 <- factor(summed.sub$annotation_1)
  summed.sub$batch <- factor(summed.sub$batch)
  
  # taking the second of comp
  summed.sub$annotation_1 <- relevel(summed.sub$annotation_1, ref=comparisons[[i]][["comp"]][2])
  print(model.matrix(~batch + annotation_1, data=colData(summed.sub)))
  
  between.res <- pseudoBulkDGE(summed.sub,
                               label=rep("dummy", ncol(summed.sub)),
                               design=~batch + annotation_1, 
                               coef=comparisons[[i]][["coef"]])[[1]]
  
  print(table(Sig=between.res$FDR <= 0.05, Sign=sign(between.res$logFC)))
  
  res <- data.frame(between.res)
  assign(paste0(i,"res"), res)
  
  p <- EnhancedVolcano(between.res,
                  lab = rownames(between.res),
                  x = 'logFC',
                  y = 'PValue',
                  title = i,
                  subtitle = NULL,
                  pCutoff = 0.05,
                  FCcutoff = 0.5, 
                  labSize = 3.5)
  
  print(p)
}



plotTSNE(sce.ann.1, colour_by = "IL33", point_size=1, text_by ="annotation_1", text_size =5, text_colour = "red")

