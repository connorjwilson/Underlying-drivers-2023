library(layeranalyzer)

# Load time bins

stages <- read.csv("data/stages.csv")

###############
# Import data
###############

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


#########################
# layeranalyzer analysis
#########################


#########################
# Origination analysis:
#########################


orig <- traverse.standalone.layered(origts,
  max.layers=2, talkative=T,allow.one.feedback.loop=T, 
  just.stationary=F,no.rw=F, time.integrals.possible=F,                                    
  allow.deterministic.layers=T,do.maximum.likelihood=T, 
  maximum.likelihood.numstart=1000, num.MCMC=1000, spacing=10, 
  burnin=2000, num.temp = 4)
# Remove linear trend models:
orig.nolin=list(orig[[1]],orig[[2]],orig[[4]],orig[[5]],
  orig[[6]],orig[[7]],orig[[8]],orig[[9]])

compare.layered(orig.nolin)
#          weight=-0.5*AIC Post. Prob.(%)
#Model   1       -112.0541        0.02269
#Model   2       -137.2967        0.00000
#Model   3       -103.7218       94.28720
#Model   4       -106.6933        4.82992
#Model   5       -111.6208        0.03499
#Model   6       -112.5777        0.01344
#Model   7       -108.4767        0.81176
#Model   8       -139.2663        0.00000

summary.layered(orig.nolin[[3]])
#                  ML estimate Bayesian Lower 95% Bayesian #Upper 95%
#mu_origts      -3.705739   -3.642893   -3.272083
#dt_origts_1    1.375067    0.027746    2.491004
#sigma_origts_1 0.766006    0.030247    2.749660
#dt_origts_2    83.561837   0.070320    9.903549
#sigma_origts_2 0.005593    0.008374    4.453098
#init_or_l_s0   -1.508269   -1.614393   -1.403538
#init_or_l2_s0  -2.112977   -70.592484  102.961823

save(orig.nolin, file="single_orig.RData")


# Reconstruct both layers:
orig.struct=layer.series.structure(origts, numlayers=2,
   init.0=TRUE)
orig.mod=layer.analyzer(orig.struct,    
     do.maximum.likelihood=TRUE,
     maximum.likelihood.numstart=1000,
     num.MCMC=1000, spacing=10, burnin=2000, num.temp = 4,
     smoothing.specs=
     list(do.smoothing=TRUE,smoothing.time.diff=1,
          smoothing.start=NULL,smoothing.end=NULL,
          num.smooth.per.mcmc=10,
     do.return.smoothing.samples=FALSE)  )

save(orig.mod, file="data/results/layer_reconstructions/orig.mod.Rdata")





#########################
# extinction analysis:
#########################


ext <- traverse.standalone.layered(extts, max.layers=2,
   talkative=T,allow.one.feedback.loop=T, just.stationary=F,
   no.rw=F, time.integrals.possible=F,
   allow.deterministic.layers=T, do.maximum.likelihood=T, 
   maximum.likelihood.numstart=1000, 
   num.MCMC=1000, spacing=10, burnin=2000, 
   num.temp = 4)
# Remove linear trend models:
ext.nolin=list(ext[[1]],ext[[2]],ext[[4]],ext[[5]],
  ext[[6]],ext[[7]],ext[[8]],ext[[9]])

compare.layered(ext.nolin)
#          weight=-0.5*AIC Post. Prob.(%)
#Model   1       -132.1773       33.41984
#Model   2       -168.3211        0.00000
#Model   3       -135.1171        1.76701
#Model   4       -135.8273        0.86858
#Model   5       -134.1773        4.52276
#Model   6       -135.1773        1.66388
#Model   7       -131.6302       57.75793
#Model   8       -170.0030        0.00000

summary.layered(ext.nolin[[7]])
# ML estimate Bayesian Lower 95% Bayesian Upper 95%
#mu_extts      -8.356207  -9.664402    10.212970
#dt_extts_1    0.017442   0.006677     0.773023
#sigma_extts_1 9.146166   1.278735     14.468617
#sigma_extts_2 0.029938   0.014742     0.096901
#init_ex_l_s0  -4.336193  -4.725734    -3.957657
#init_ex_l_s0  -3.000669  -3.582532    -2.113587

# Trond: 2 layers with RW bottom layer. 

save(ext.nolin, file="single_ext.RData")

# Reconstruct both layers:
ext.struct=layer.series.structure(extts, numlayers=2,
   no.pull = TRUE, init.0=TRUE)
ext.mod=layer.analyzer(ext.struct,    
     do.maximum.likelihood=TRUE,
     maximum.likelihood.numstart=1000,
     num.MCMC=1000, spacing=10, burnin=2000, num.temp = 4,
     silent.mode=FALSE,
     smoothing.specs=
     list(do.smoothing=TRUE,smoothing.time.diff=0.1,
          smoothing.start=-600,smoothing.end=0,
          num.smooth.per.mcmc=10,
     do.return.smoothing.samples=FALSE)  )

save(ext.mod, file="data/results/layer_reconstructions/ext.mod.Rdata")

# Trond - interpretation: Measurement noise probably
# underestimated. First layer is there to add the extra noise
# necessary to explain the variation in the measurements that
# is lacking in the measurement standard deviation. 
# The second layer is probably the extinction process itself.




#########################
# Sampling analysis:
#########################


samp <- traverse.standalone.layered(sampts, max.layers=2,
   talkative=T, allow.one.feedback.loop=T, just.stationary=F,
   no.rw=F, time.integrals.possible=F,
   allow.deterministic.layers=T, do.maximum.likelihood=T, 
   maximum.likelihood.numstart=1000, 
   num.MCMC=1000, spacing=10, burnin=2000, 
   num.temp = 4)
# Remove linear trend models:
samp.nolin=list(samp[[1]],samp[[2]],samp[[4]],samp[[5]],
  samp[[6]],samp[[7]],samp[[8]],samp[[9]])
compare.layered(samp.nolin)
#          weight=-0.5*AIC Post. Prob.(%)
#Model   1       -84.34685        2.64828
#Model   2       -96.79618        0.00001
#Model   3       -81.31273       55.03815
#Model   4       -82.35619       19.38636
#Model   5       -86.40285        0.33889
#Model   6       -87.45627        0.11818
#Model   7       -82.20857       22.47014
#Model   8       -98.22418        0.00000

summary.layered(samp.nolin[[3]])
#   ML estimate Bayesian Lower95% Bayesian Upper 95%
#mu_sampts     -2.143221 -2.513120 -1.826869
#dt_sampts_1   1.350923  0.226192  17.540403
#sigma_s_1     0.523857  0.029248  0.950380
#dt_sampts_2   42.036325 0.889329  74.239554
#sigma_s_2     0.115950  0.012924   1.301279
#init_s_l2_s0  -0.618456 -1.225231  -0.054444
#init_s_l2_s0  -1.584356 -19.095603 49.310120

save(samp.nolin, file="single_samp.RData")


# Reconstruct both layers:
samp.struct=layer.series.structure(sampts, numlayers=2,
   no.pull = TRUE, init.0=TRUE)
samp.mod=layer.analyzer(samp.struct,    
     do.maximum.likelihood=TRUE,
     maximum.likelihood.numstart=1000,
     num.MCMC=1000, spacing=10, burnin=2000, num.temp = 4,
     silent.mode=FALSE,
     smoothing.specs=
     list(do.smoothing=TRUE,smoothing.time.diff=0.1,
          smoothing.start=-600,smoothing.end=0,
          num.smooth.per.mcmc=10,
     do.return.smoothing.samples=FALSE)  )

save(samp.mod, file="data/results/layer_reconstructions/samp.mod.Rdata")




