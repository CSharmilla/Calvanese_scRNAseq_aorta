library(here)
library(EnhancedVolcano)
library(scran)



sce.ann1.sub <- readRDS(here("output", "sce.ann1.sub.RDS"))
plotTSNE(sce.ann1.sub, colour_by = "aec.sub", point_size=1)


summed <- aggregateAcrossCells(sce.ann1.sub, id=colData(sce.ann1.sub)[,c("HSC.sub", "label")])

# comparisons <- list(aec_2_vs_aec_3 = list(comp = c("aec_2", "aec_3"),
#                                       coef = "aec.subaec_2"),
#                     pre_HE_vs_aec_2 = list(comp = c("pre_HE", "aec_2"),
#                                         coef = "aec.subpre_HE"), 
#                     pre_HE_vs_aec_3 = list(comp = c("pre_HE", "aec_3"),
#                                           coef = "aec.subpre_HE"),
#                     HE_vs_pre_HE = list(comp = c("HE", "pre_HE"),
#                                            coef = "aec.subHE"))

comparisons <- list(HSC.2_vs_HSC.1 = list(comp = c("HSC.2", "HSC.1"),
                                          coef = "HSC.subHSC.2"))

for (i in names(comparisons)){
  #i <- "hsc_vs_hema1"
  summed.sub <- summed[,summed$HSC.sub %in% comparisons[[i]][["comp"]]]
  summed.sub$HSC.sub <- factor(summed.sub$HSC.sub)
  summed.sub$batch <- factor(summed.sub$batch)
  
  # taking the second of comp
  summed.sub$HSC.sub <- relevel(summed.sub$HSC.sub, ref=comparisons[[i]][["comp"]][2])
  print(model.matrix(~batch + HSC.sub, data=colData(summed.sub)))
  
  between.res <- pseudoBulkDGE(summed.sub,
                               label=rep("dummy", ncol(summed.sub)),
                               design=~batch + HSC.sub, 
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



plotTSNE(sce.ann1.sub, colour_by = "TOP2A", point_size=1, text_by ="HSC.sub", text_size =5, text_colour = "red")
plotTSNE(sce.ann1.sub, colour_by = "DLL4", point_size=1, text_by ="aec.sub", text_size =5, text_colour = "red")


