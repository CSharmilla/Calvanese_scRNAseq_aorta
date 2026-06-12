library(here)
library(bluster)
library(scran)
library(dplyr)


sce.ann.1 <- readRDS(here("output", "sce.ann.1.RDS"))

if (!file.exists(here("output", "subcluster_list.RDS"))){
  subcluster_list <- quickSubCluster(sce.ann.1, sce.ann.1$annotation_1, 
                                     prepFUN=function(x) {
                                       set.seed(0010101010)
                                       dec <- modelGeneVarByPoisson(x)
                                       top <- getTopHVGs(dec, n=5000)
                                       denoisePCA(x, technical=dec,
                                                  subset.row=top, name = "pca_denoise")
                                     },
                                     clusterFUN=function(x) {
                                       set.seed(001010)
                                       clusterRows(reducedDim(x, "pca_denoise"), 
                                                   NNGraphParam(k = 30, cluster.fun="louvain"))
                                     })
  saveRDS(subcluster_list, file = here("output", "subcluster_list.RDS"))
} else{
  subcluster_list <- readRDS(here("output", "subcluster_list.RDS"))
  
}


# renaming aEC
plotTSNE(subcluster_list$aEC, colour_by = "subcluster", point_size=1)

plotTSNE(subcluster_list$aEC, colour_by = "IL33", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$aEC, colour_by = "ALDH1A1", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$aEC, colour_by = "KCNK17", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$aEC, colour_by = "RUNX1", point_size=1, text_by = "subcluster", text_colour = "red")

aec_cluster <- c(paste0("aEC.", 1:4))
aec.ann <- c("HE", "aec_2", "aec_3", "pre_HE")
names(aec.ann) <- aec_cluster
subcluster_list$aEC$subcluster <- factor(subcluster_list$aEC$subcluster)
aec.coldata <- colData(subcluster_list$aEC) %>% as.data.frame %>% mutate(aec.sub=recode_factor(subcluster, !!!aec.ann)) %>% as("DataFrame")

all(colnames(subcluster_list$aEC) == rownames(aec.coldata))
colData(subcluster_list$aEC) <- aec.coldata

sub_cells <- colnames(subcluster_list$aEC)
sub_labels <- subcluster_list$aEC$aec.sub

full <- rep(NA, ncol(sce.ann.1))
names(full) <- colnames(sce.ann.1)

full[sub_cells] <- as.character(sub_labels)

sce.ann1.sub <- sce.ann.1
sce.ann1.sub$aec.sub <- full
plotTSNE(sce.ann1.sub, colour_by = "aec.sub", point_size=1)
plotTSNE(subcluster_list$aEC, colour_by = "DLL4", point_size=1, text_by = "aec.sub", text_colour = "red")

# combining subclusters with the full annotation column
#sce.ann1.aec.sub$comb_ann1 <- paste0(sce.ann1.aec.sub$annotation_1, "_", sce.ann1.aec.sub$aec.sub)
#plotTSNE(sce.ann1.aec.sub, colour_by = "comb_ann1", point_size=1)

saveRDS(sce.ann1.sub, file = here("output", "sce.ann1.sub.RDS"))



# renaming hema_2
plotTSNE(subcluster_list$hema_2, colour_by = "subcluster", point_size=1)

plotTSNE(subcluster_list$hema_2, colour_by = "SPN", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$hema_2, colour_by = "PTPRC", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$hema_2, colour_by = "KCNK17", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$hema_2, colour_by = "RUNX1", point_size=1, text_by = "subcluster", text_colour = "red")

# hema_2_cluster <- c(paste0("hema_2.", 1:4))
# hema_2.ann <- c("HE", "hema_2_2", "hema_2_3", "pre_HE")
# names(hema_2.ann) <- hema_2_cluster
# subcluster_list$hema_2$subcluster <- factor(subcluster_list$hema_2$subcluster)
# hema_2.coldata <- colData(subcluster_list$hema_2) %>% as.data.frame %>% mutate(hema_2.sub=recode_factor(subcluster, !!!hema_2.ann)) %>% as("DataFrame")
# 
# all(colnames(subcluster_list$hema_2) == rownames(hema_2.coldata))
# colData(subcluster_list$hema_2) <- hema_2.coldata

sub_cells <- colnames(subcluster_list$hema_2)
sub_labels <- subcluster_list$hema_2$subcluster

full <- rep(NA, ncol(sce.ann.1))
names(full) <- colnames(sce.ann.1)

full[sub_cells] <- as.character(sub_labels)

#sce.ann1.sub <- sce.ann.1
sce.ann1.sub$hema_2.sub <- full
plotTSNE(sce.ann1.sub, colour_by = "hema_2.sub", point_size=1)
plotTSNE(subcluster_list$hema_2, colour_by = "DLL4", point_size=1, text_by = "hema_2.sub", text_colour = "red")

# combining subclusters with the full annotation column
#sce.ann1.aec.sub$comb_ann1 <- paste0(sce.ann1.aec.sub$annotation_1, "_", sce.ann1.aec.sub$aec.sub)
#plotTSNE(sce.ann1.aec.sub, colour_by = "comb_ann1", point_size=1)

saveRDS(sce.ann1.sub, file = here("output", "sce.ann1.sub.RDS"))


# renaming HSC
plotTSNE(subcluster_list$HSC, colour_by = "subcluster", point_size=1)

plotTSNE(subcluster_list$HSC, colour_by = "SPN", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$HSC, colour_by = "PTPRC", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$HSC, colour_by = "KCNK17", point_size=1, text_by = "subcluster", text_colour = "red")
plotTSNE(subcluster_list$HSC, colour_by = "RUNX1", point_size=1, text_by = "subcluster", text_colour = "red")

# HSC_cluster <- c(paste0("HSC.", 1:4))
# HSC.ann <- c("HE", "HSC_2", "HSC_3", "pre_HE")
# names(HSC.ann) <- HSC_cluster
# subcluster_list$HSC$subcluster <- factor(subcluster_list$HSC$subcluster)
# HSC.coldata <- colData(subcluster_list$HSC) %>% as.data.frame %>% mutate(HSC.sub=recode_factor(subcluster, !!!HSC.ann)) %>% as("DataFrame")
# 
# all(colnames(subcluster_list$HSC) == rownames(HSC.coldata))
# colData(subcluster_list$HSC) <- HSC.coldata

sub_cells <- colnames(subcluster_list$HSC)
sub_labels <- subcluster_list$HSC$subcluster

full <- rep(NA, ncol(sce.ann.1))
names(full) <- colnames(sce.ann.1)

full[sub_cells] <- as.character(sub_labels)

#sce.ann1.sub <- sce.ann.1
sce.ann1.sub$HSC.sub <- full
plotTSNE(sce.ann1.sub, colour_by = "HSC.sub", point_size=1)
plotTSNE(subcluster_list$HSC, colour_by = "DLL4", point_size=1, text_by = "HSC.sub", text_colour = "red")

# combining subclusters with the full annotation column
#sce.ann1.aec.sub$comb_ann1 <- paste0(sce.ann1.aec.sub$annotation_1, "_", sce.ann1.aec.sub$aec.sub)
#plotTSNE(sce.ann1.aec.sub, colour_by = "comb_ann1", point_size=1)

saveRDS(sce.ann1.sub, file = here("output", "sce.ann1.sub.RDS"))





