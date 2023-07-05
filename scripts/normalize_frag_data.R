# Normalizes fragmentation index


# Read data:
#lfrag.orig=read.table("data/logit_frag_data.csv",sep=",",header=T)

frag_data <- read.csv("data/fragmentation_data.csv")
frag_data$logit.fragmentation.index <- log(frag_data$fragmentation.index/(1-frag_data$fragmentation.index))
lfrag.orig <- frag_data

# Order chronologically:
lfrag.orig=lfrag.orig[order(lfrag.orig$time),]





# Check normality
qqnorm(lfrag.orig$logit.fragmentation.index)
qqline(lfrag.orig$logit.fragmentation.index)
# Samples seem to have lighter tails than the 
# normal distribution! (I.e. seems more restricted to be 
# close to the mean than the normal distribution is. 
# Usually data is more heavy-tailed.)

shapiro.test(lfrag.orig$logit.fragmentation.index)
# p-value = 1.134e-11




# Transform the data:

# Store it in an easyly accessible variable. 
# We need to reference this quite a few times.
y=lfrag.orig$logit.fragmentation.index

# Make cumulative density approximation:
# First, determine accuracy from the range of the data.
scale=10^ceiling(log10(max(y)-min(y)))
# Density approximation from R's kernel density estimation method:
dd=density(y, n=100001, from=min(y)-scale, to=max(y)+scale, adjust=0.2) 

hist(y,freq=F)
lines(dd)
# Seems to also have some signs of multimodality
# Could be too detailed, but it seems to need this detail
# level to pass the Shapiro test.

# Find cumulative distribution:
dd$F=cumsum(dd$y)*(dd$x[2]-dd$x[1])
dd$F=dd$F/max(dd$F)

# Transform values to something standard normalized-like:
trans=function(x)
  qnorm(dd$F[max(which(dd$x<x))])

# Transform values back again:
invtrans=function(x)
  dd$x[min(which(dd$F>pnorm(x)))]

# Transform data:
x=0*y
for(i in 1:length(x))
  x[i]=trans(y[i])

shapiro.test(x)
# p-value = 0.06882

# Put into a new table:
lfrag.new=lfrag.orig[,3:5]
lfrag.new$trans.frag.index=x

# Write transformed data to file:
write.table(lfrag.new,file="data/lfrag_norm.txt",
  col.names=T,row.names=F,sep=",")

