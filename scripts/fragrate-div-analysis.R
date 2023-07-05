########################
# fragrate-div-analysis
########################

library(openCR)
library(divDyn)
library(layeranalyzer)

# Load time bins

stages <- read.csv("data/stages.csv")

###############
# Import data
###############

# Import smoothed fragmentation rate
frag_data_s <- read.csv2("data/smoothed_frag_rate.csv")
fragts_sm <- layer.data.series(time.points=as.numeric(frag_data_s$time),
                               value.points=frag_data_s$value,
                               name="fragts_sm")
fragsm_struct <- layer.series.structure(fragts_sm, numlayers=1)

frag_data_ss <- frag_data_s[seq(from=17, by=20, length.out=533),]
fragts_ss <- layer.data.series(time.points=as.numeric(frag_data_ss$time),
                               value.points=frag_data_ss$value,
                               name="fragts_ss")

fragss_struct <- layer.series.structure(fragts_ss, numlayers=1)

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

results_origN2 <- list()
for (i in 1:1) {
  orig_null <- layer.analyzer(fragss_struct, origts_struct, do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_caus <- layer.analyzer(fragss_struct, origts_struct, causal=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_corr <- layer.analyzer(fragss_struct, origts_struct, corr=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_caus2 <- layer.analyzer(fragss_struct, origts_struct, causal=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_corr2 <- layer.analyzer(fragss_struct, origts_struct, corr=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  orig_models <- list(orig_null, orig_caus, orig_corr, orig_caus2, orig_corr2)
  
  results_origN2[[i]] <- orig_models
  
  rm(orig_null, orig_caus, orig_corr, orig_caus2, orig_corr2, orig_models)
  #save(results_origN2, file="data/results/results_origN2.Rdata") if testing stability, can save results after every run
}



results_extN2 <- list()
for (i in 1:1) {
  ext_null <- layer.analyzer(fragss_struct, extts_struct, do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_caus <- layer.analyzer(fragss_struct, extts_struct, causal=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_corr <- layer.analyzer(fragss_struct, extts_struct, corr=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_caus2 <- layer.analyzer(fragss_struct, extts_struct, causal=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_corr2 <- layer.analyzer(fragss_struct, extts_struct, corr=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  ext_models <- list(ext_null, ext_caus, ext_corr, ext_caus2, ext_corr2)
  
  results_extN2[[i]] <- ext_models
  
  rm(ext_null, ext_caus, ext_corr, ext_caus2, ext_corr2, ext_models)
  #save(results_extN2, file="data/results/results_extN2.Rdata")
}

results_sampN2 <- list()
for (i in 1:1) {
  samp_null <- layer.analyzer(fragss_struct, sampts_struct, do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_caus <- layer.analyzer(fragss_struct, sampts_struct, causal=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_corr <- layer.analyzer(fragss_struct, sampts_struct, corr=matrix(c(1,1,2,1)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_caus2 <- layer.analyzer(fragss_struct, sampts_struct, causal=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_corr2 <- layer.analyzer(fragss_struct, sampts_struct, corr=matrix(c(1,1,2,2)), do.maximum.likelihood=TRUE, num.MCMC=100000, maximum.likelihood.numstart=10000)
  samp_models <- list(samp_null, samp_caus, samp_corr, samp_caus2, samp_corr2)
  
  results_sampN2[[i]] <- samp_models
  
  rm(samp_null, samp_caus, samp_corr, samp_caus2, samp_corr2, samp_models)
  #save(results_sampN2, file="data/results/results_sampN2.Rdata")
}

rate_orig <- results_origN2[[1]]
save(rate_orig, file="data/results/rate_orig.RData")
rate_ext <- results_extN2[[1]]
save(rate_ext, file="data/results/rate_ext.RData")
rate_samp <- results_sampN2[[1]]
save(rate_samp, file="data/results/rate_samp.RData")

