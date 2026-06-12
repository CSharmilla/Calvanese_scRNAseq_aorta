library(here)
library(scran)



all.sce.norm <- readRDS(here("output", "all.sce.norm.RDS"))

### variance modelling
### with density weights
### variance modelling - by quantifying technical noise
all.var <- lapply(all.sce.norm, modelGeneVarByPoisson)


### plotting
par(mfrow=c(3,4)) ### this works with graphics package stuff (like plot function)
for (n in names(all.var)){
  var <- all.var[[1]]
  var <- var[order(var$bio, decreasing=TRUE),]
  fit <- metadata(var)
  graphics::plot(var$mean, var$total, xlab="mean of log-expression", ylab="variance of log-expression", main = n) 
  curve(fit$trend(x), col="blue", add = TRUE, lwd = 3)
}

### get HVGs
all.hvgs <- lapply(all.var, getTopHVGs, prop=0.15)
### look at how many are there
lapply(all.hvgs, length)
### subset HVGs
for (n in names(all.sce.norm)){
  rowSubset(all.sce.norm[[n]]) <- all.hvgs[[n]]
}

saveRDS(all.sce.norm, file = here("output", "all.sce.norm.RDS"))
saveRDS(all.var, file = here("output", "all.var.RDS"))
saveRDS(all.hvgs, file = here("output", "all.hvgs.RDS"))

