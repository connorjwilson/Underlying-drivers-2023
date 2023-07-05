########################
# fragrate-div-analysis
########################

library(openCR)
library(layeranalyzer)

# Load time bins
stages <- read.csv("data/stages.csv")

###############
# Import data
###############

# Import fragmentation index
frag_index <- read.csv("data/fragmentation_data.csv")
fragts_in <- layer.data.series(time.points=as.numeric(frag_index$time),
                               value.points=frag_index$fragmentation.index,
                               name="fragts_in")
fragin_struct <- layer.series.structure(fragts_in, numlayers=1)

# Import invert data
load("data/results/pradel_results.RData")
load("data/results/pradel_stdev.RData")

origts <- layer.data.series(time.points=(-stages$bottom[6:95]),
                            value.points=pradel_results[[1]][3:92],
                            std.dev=pradel_stdev[[1]][3:92],
                            name="origts")

extts <- layer.data.series(time.points=(-stages$top[4:93]),
                           value.points=pradel_results[[2]][1:90],
                           std.dev=pradel_stdev[[2]][1:90],
                           name="extts")

sampts <- layer.data.series(time.points=(-stages$mid[5:94]),
                            value.points=pradel_results[[3]][2:91],
                            std.dev=pradel_stdev[[3]][2:91],
                            name="sampts")

origts_struct <- layer.series.structure(origts, numlayers=2)
extts_struct <- layer.series.structure(extts, numlayers=2)
sampts_struct <- layer.series.structure(sampts, numlayers=2)

#########################
# layeranalyzer analysis
#########################

results_iorigN2 <- list()
for (i in 1:1) {
  orig_null <- layer.analyzer(fragin_struct, origts_struct, do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_caus <- layer.analyzer(fragin_struct, origts_struct, causal=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_corr <- layer.analyzer(fragin_struct, origts_struct, corr=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_caus2 <- layer.analyzer(fragin_struct, origts_struct, causal=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_corr2 <- layer.analyzer(fragin_struct, origts_struct, corr=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_models <- list(orig_null, orig_caus, orig_corr, orig_caus2, orig_corr2)
  
  results_iorigN2[[i]] <- orig_models
  
  rm(orig_null, orig_caus, orig_corr, orig_models, orig_caus2, orig_corr2)
  #save(results_iorigN2, file="data/results/results_iorigN2.Rdata") if testing stability, can save results after every run
}



results_iextN2 <- list()
for (i in 1:1) {
  ext_null <- layer.analyzer(fragin_struct, extts_struct, do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_caus <- layer.analyzer(fragin_struct, extts_struct, causal=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_corr <- layer.analyzer(fragin_struct, extts_struct, corr=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_caus2 <- layer.analyzer(fragin_struct, extts_struct, causal=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_corr2 <- layer.analyzer(fragin_struct, extts_struct, corr=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_models <- list(ext_null, ext_caus, ext_corr, ext_caus2, ext_corr2)
  
  results_iextN2[[i]] <- ext_models
  
  rm(ext_null, ext_caus, ext_corr, ext_caus2, ext_corr2, ext_models)
  #save(results_iextN2, file="data/results/results_iextN2.Rdata")
}

results_isampN2 <- list()
for (i in 1:1) {
  samp_null <- layer.analyzer(fragin_struct, sampts_struct, do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_caus <- layer.analyzer(fragin_struct, sampts_struct, causal=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_corr <- layer.analyzer(fragin_struct, sampts_struct, corr=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_caus2 <- layer.analyzer(fragin_struct, sampts_struct, corr=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_corr2 <- layer.analyzer(fragin_struct, sampts_struct, corr=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_models <- list(samp_null, samp_caus, samp_corr, samp_caus2, samp_corr2)
  
  results_isampN2[[i]] <- samp_models
  
  rm(samp_null, samp_caus, samp_corr, samp_saus2, camp_corr2, samp_models)
  #save(results_isampN2, file="data/results/results_isampN2.Rdata")
}

# View model weights
for (i in 1:10) { print(compare.layered(results_iext[[i]], ML.IC="AICc"))}


index_orig <- results_iorigN2[[1]]
save(index_orig, file="data/results/index_orig.RData")
index_ext <- results_iextN2[[1]]
save(index_ext, file="data/results/index_ext.RData")
index_samp <- results_isampN2[[1]]
save(index_samp, file="data/results/index_samp.RData")
