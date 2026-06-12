library(here)
library(EnhancedVolcano)
library(scran)



sce.ann.sub.1 <- subcluster_list$aEC #readRDS(here("output", "sce.ann.sub.1.RDS"))
plotTSNE(sce.ann.sub.1, colour_by = "aec.sub", point_size=1)

a <- data.frame(colData(summed2))

summed2 <- aggregateAcrossCells(sce.ann.sub.1, id=colData(sce.ann.sub.1)[,c("aec.sub", "label")])

comparisons <- list(aec_2_vs_aec_3 = list(comp = c("aec_2", "aec_3"),
                                          coef = "aec.subaec_2"),
                    pre_HE_vs_aec_2 = list(comp = c("pre_HE", "aec_2"),
                                           coef = "aec.subpre_HE"), 
                    pre_HE_vs_aec_3 = list(comp = c("pre_HE", "aec_3"),
                                           coef = "aec.subpre_HE"),
                    HE_vs_pre_HE = list(comp = c("HE", "pre_HE"),
                                        coef = "aec.subHE"))

for (i in names(comparisons)){
  #i <- "hsc_vs_hema1"
  summed.sub <- summed2[,summed2$aec.sub %in% comparisons[[i]][["comp"]]]
  summed.sub$aec.sub <- factor(summed.sub$aec.sub)
  summed.sub$batch <- factor(summed.sub$batch)
  
  # taking the second of comp
  summed.sub$aec.sub <- relevel(summed.sub$aec.sub, ref=comparisons[[i]][["comp"]][2])
  print(model.matrix(~batch + aec.sub, data=colData(summed.sub)))
  
  between.res <- pseudoBulkDGE(summed.sub,
                               label=rep("dummy", ncol(summed.sub)),
                               design=~batch + aec.sub, 
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



plotTSNE(sce.ann1.sub, colour_by = "FREM1", point_size=1, text_by ="aec.sub", text_size =5, text_colour = "red")
plotTSNE(subcluster_list$aEC, colour_by = "PCDH7", point_size=1, text_by ="aec_ann", text_size =5, text_colour = "red")
plotTSNE(sce.ann.sub.1, colour_by = "PCDH7", point_size=1, text_by ="aec.sub", text_size =5, text_colour = "red")

