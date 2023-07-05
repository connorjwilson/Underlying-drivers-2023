# Read normalized fragmentation index

# logit_frag <- log(frag_data$fragmentation.index/(1-frag_data$fragmentation.index))

frag=read.table("data/lfrag_norm.txt",sep=",",header=T)

# Put into layeranalyzer context:
library(layeranalyzer)

lfrag=layer.data.series(frag$time, frag$trans.frag.index,
                        "frag")

# Restrict characteristic time and expected value
pr1=layer.prior(mu=c(-0.0000002,0.0000002),init=c(-0.0000002,0.0000002),
                dt=c(2,5),sigma=c(1e-6,1000),
                obs=c(0.01,0.1),lin=c(-0.01,0.01), beta=c(-1000,1000))

frag.struct=layer.series.structure(lfrag, numlayers=2,
                                   time.integral=c(1), prior=pr1) 

# Test:
test1=layer.analyzer(frag.struct, 
                     silent.mode=F, talkative.burnin=T, 
                     smoothing.specs=
                       list(do.smoothing=TRUE,smoothing.time.diff=0.05,
                            smoothing.start=NULL,smoothing.end=NULL,
                            num.smooth.per.mcmc=10,   
                            do.return.smoothing.samples=TRUE) )

summary.layered(test1)

smoothed_frag_rate <- cbind(test1$process.time.points, test1$process.mean)
write.csv(smoothed_frag_rate, file="smoothed_frag_rate.csv")