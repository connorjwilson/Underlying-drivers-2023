#########################
###
### CMR-analysis.R
###
#########################

library(openCR)

stages <- read.csv("data/stages.csv")

#################################
# Import and prep occurrence data
#################################

occ_data <- read.csv("data/occ_data.csv", header=TRUE, skip=19)
occ_data$mid_ma <- rowMeans(occ_data[,15:16])



# Prepare occurrence data for conversion to capture history
breaks <- c(-stages$bottom, 0)
stages$bins <- cut(-(stages$mid), breaks=breaks)
stages$mid <- -(stages$mid)
breaks <- c(-stages$bottom, 0)
occ_data$bins <- cut(-(occ_data$mid_ma), breaks=breaks)
stages_1 <- subset(stages, select=c(stage, stg, bins))
occ_data <- merge(occ_data, stages_1, by="bins")

# Create detection record
detection <- as.data.frame(rep(1, nrow(occ_data)))
colnames(detection) <- "Session"
detection$ID <- occ_data$genus
detection$Occasion <- occ_data$stg
detection$Detector <- 1
detection$Occasion <- detection$Occasion - 3
detection <- detection[-which(detection$ID==""),]

# Create trap record
dummytraps <- data.frame(row.names=list(1))
dummytraps$Detector <- 1
dummytraps$x <- 1
dummytraps$y <- 1
dummytraps <- read.traps(data=dummytraps, detector="multi")


#########################
# Run Pradel CMR analysis
#########################

# Create capture history
capt_matrix <- suppressWarnings(make.capthist(detection, dummytraps, fmt="trapID"))

# Run CMR analysis
pradel_1 <- (openCR.fit(capt_matrix, type="Pradelg", model=list(p~t, phi~t, gamma~t)))
save(pradel_1, file="data/pradel_1.RData")

pradel1_est <- predict(pradel_1)

#######################
# Plot CMR results
#######################

# Extinction
plot(-stages$top[4:nrow(stages)], 1-pradel1_est$phi$estimate, main="Extinction probability", type="b")

plot(-stages$top[4:nrow(stages)], -log(pradel1_est$phi$estimate)/stages$dur[4:nrow(stages)], type="b", main="Extinction rate", log="y")
arrows(-stages$top[4:nrow(stages)], -log(pradel1_est$phi$lcl)/stages$dur[4:nrow(stages)], -stages$top[4:nrow(stages)], -log(pradel1_est$phi$ucl)/stages$dur[4:nrow(stages)], code=3, angle=90, length=0.1, col="black", lty=2)


# Origination
plot(-stages$top[4:nrow(stages)], 1-pradel1_est$gamma$estimate, main="Origination probability", type="b")

plot(-stages$top[4:nrow(stages)], -log(pradel1_est$gamma$estimate)/stages$dur[4:nrow(stages)], type="b", main="Origination rate", log="y")
arrows(-stages$top[4:nrow(stages)], -log(pradel1_est$gamma$lcl)/stages$dur[4:nrow(stages)], -stages$top[4:nrow(stages)], -log(pradel1_est$gamma$ucl)/stages$dur[4:nrow(stages)], code=3, angle=90, length=0.1, col="black", lty=2)


# Sampling
plot(-stages$top[4:nrow(stages)], pradel1_est$p$estimate, main="Sampling probability", type="b")

plot(-stages$top[4:nrow(stages)], -log(1-pradel1_est$p$estimate)/stages$dur[4:nrow(stages)], type="b", main="Sampling rate", log="y")
arrows(-stages$top[4:nrow(stages)], -log(1-pradel1_est$p$lcl)/stages$dur[4:nrow(stages)], -stages$top[4:nrow(stages)], -log(1-pradel1_est$p$ucl)/stages$dur[4:nrow(stages)], code=3, angle=90, length=0.1, col="black", lty=2)

###############
# Save results
###############

pradel1_est$gamma$estimate[c(1,2)] <- NA
pradel1_est$phi$estimate[c(91,92)] <- NA
pradel1_est$p$estimate[c(1,length(pradel1_est$p$estimate))] <- NA

pradel1_est$gamma$lcl[c(1,2)] <- NA
pradel1_est$phi$lcl[c(91,92)] <- NA
pradel1_est$p$lcl[c(1,length(pradel1_est$p$estimate))] <- NA
pradel1_est$gamma$ucl[c(1,2)] <- NA
pradel1_est$phi$ucl[c(91,92)] <- NA
pradel1_est$p$ucl[c(1,length(pradel1_est$p$estimate))] <- NA

# Estimates
pradel_results <- list(log(-log(pradel1_est$gamma$estimate)/stages$tdur[3:(nrow(stages)-1)]), log(-log(pradel1_est$phi$estimate)/stages$tdur[4:nrow(stages)]), log(-log(1-pradel1_est$p$estimate)/stages$dur[4:nrow(stages)]))

plot(-stages$bottom[4:nrow(stages)], pradel_results[[1]], type="b", main="Origination", xlab="", ylab="Log origination rate")
arrows(-stages$bottom[4:nrow(stages)], log(-log(pradel1_est$gamma$lcl)/stages$tdur[3:(nrow(stages)-1)]), -stages$bottom[4:nrow(stages)], log(-log(pradel1_est$gamma$ucl)/stages$tdur[3:(nrow(stages)-1)]), code=3, angle=90, length=0.1, col="black", lty=2)


plot(-stages$top[4:nrow(stages)], pradel_results[[2]], type="b", main="Extinction", xlab="", ylab="Log extinction rate")
arrows(-stages$top[4:nrow(stages)], log(-log(pradel1_est$phi$lcl)/stages$tdur[4:(nrow(stages))]), -stages$top[4:nrow(stages)], log(-log(pradel1_est$phi$ucl)/stages$tdur[4:(nrow(stages))]), code=3, angle=90, length=0.1, col="black", lty=2)


plot(-stages$mid[4:nrow(stages)], pradel_results[[3]], type="b", main="Sampling", xlab="Mya", ylab="Log sampling rate")
arrows(-stages$mid[4:nrow(stages)], log(-log(1-pradel1_est$p$lcl)/stages$dur[4:(nrow(stages))]), -stages$mid[4:nrow(stages)], log(-log(1-pradel1_est$p$ucl)/stages$dur[4:(nrow(stages))]), code=3, angle=90, length=0.1, col="black", lty=2)


save(pradel_results, file="data/results/pradel_results.RData")

# Errors
pradel_stdev <- list( -(log(-log(pradel1_est$gamma$ucl)/stages$tdur[3:(nrow(stages)-1)])-log(-log(pradel1_est$gamma$lcl)/stages$tdur[3:(nrow(stages)-1)]))/3.92, 
                      -(log(-log(pradel1_est$phi$ucl)/stages$tdur[4:nrow(stages)])-log(-log(pradel1_est$phi$lcl)/stages$tdur[4:nrow(stages)]))/3.92, 
                      (log(-log(1-pradel1_est$p$ucl)/stages$dur[4:nrow(stages)])-log(-log(1-pradel1_est$p$lcl)/stages$dur[4:nrow(stages)]))/3.92)

save(pradel_stdev, file="data/results/pradel_stdev.RData")
