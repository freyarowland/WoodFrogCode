# Code for Bayesian variable selection
# written by F.E. Rowland
# Last edit January 2022

# load library
library(BayesVarSel)

# load data
alldata2 <- read.csv("woodfrogdata.csv", header = TRUE)

# compute weighted neighborhood effect from two years prior
alldata2$wEggs_t2 <- NA
alldata2$Pond.ID <- as.character(alldata2$Pond.ID)
for (i in 3:dim(alldata2)[1]) {
  if (alldata2$Pond.ID[i] == alldata2$Pond.ID[i - 2]) {
    alldata2$wEggs_t2[i] <- alldata2$wEggs[i - 2]
  } else{
    alldata2$wEggs_t2[i] <- NA
  }
}

# reduce data to only complete observations
red_data <- na.omit(alldata2)

# population growth BVS ----
growthBVS <- Bvs(formula = GrowthRate ~
                   RASYdens_t2 +
                   Depth +
                   Canopy +
                   spfaPDSI +
                   spfaPDSIlag +
                   days_thawed +
                   wEggs_t2,
                 data = red_data)
summary(growthBVS)

# plot of posterior dimension probabilities
plot(growthBVS)

# egg mass count BVS ----
countBVS <- Bvs(formula = Avg.RASY.Count ~
                  RASYdens_t2 +
                  Depth +
                  Canopy +
                  spfaPDSI +
                  spfaPDSIlag +
                  days_thawed +
                  wEggs_t2,
                data = red_data)
summary(countBVS)

# plot of posterior dimension probabilties
plot(countBVS)


# Calculate the most probable models ---
# computes it based on all variables instead of single
# inclusion probabilities
growth.Bvs <- Bvs(formula=
                   GrowthRate ~
                   RASYdens_t2 +
                   Depth +
                   Canopy +
                   spfaPDSI +
                   spfaPDSIlag +
                   days_thawed +
                   wEggs_t2,
                 data=red_data,
                 prior.betas="robust",
                 prior.models="Constant")


count.Bvs <- Bvs(formula=
                    Avg.RASY.Count ~
                    RASYdens_t2 +
                    Depth +
                    Canopy +
                    spfaPDSI +
                    spfaPDSIlag +
                    days_thawed +
                    wEggs_t2,
                  data=red_data,
                  prior.betas="robust",
                  prior.models="Constant")