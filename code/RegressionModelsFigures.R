# Code for partial pooling Bayesian models
# written by F.E. Rowland
# Last edit January 2022

# load data
alldata2 <- read.csv("woodfrogdata.csv", header = TRUE)
summary(alldata2)

# libraries
library(rstanarm)
library(tidyverse)
library(tidyr)
library(dplyr)
library(viridis)
library(ggplot2)
library(scales)
library(readr)
library(gridExtra)
library(grid)
library(colorspace)
library(lattice)
library(data.table)
library(tidyselect)
library(MASS)
library(bayesplot)
library(tibble)
library(tidybayes)
library(lattice)
library(extrafont)
library(modelr)
library(parallel)

# Pond.ID as factor
alldata2$Pond.ID <- as.factor(alldata2$Pond.ID)

# how many cores for running models?
detectCores()

# count vs. weighted eggs ----
# recalc weighted eggs t-2
alldata2$wEggs_t2 <- NA
for(i in 3:dim(alldata2)[1]){
  if(alldata2$Pond.ID[i]==alldata2$Pond.ID[i-2]){
    alldata2$wEggs_t2[i] <- alldata2$wEggs[i-2]
  } else{
    alldata2$wEggs_t2[i] <- NA
  }
}



# fit model
wEggs_Cmod <- stan_glmer(Avg.RASY.Count ~ wEggs_t2 + (wEggs_t2|Pond.ID),
                         iter = 10000,
                         na.action = "na.omit",
                         data = alldata2,
                         family = gaussian(),
                         cores = 3,
                         seed = 12345)
summary(wEggs_Cmod, digits = 5)

# leave-one-out cross validation
loo_mod <- loo(wEggs_Cmod, k_threshold = 0.7) # All are okay.

# check model fit on Shiny app
launch_shinystan(wEggs_Cmod) # model fit okay!

# model saved
saveRDS(wEggs_Cmod, "wEggs_logCmod.rds")

# load Rdata instead of re-running model
wEggs_Cmod <- readRDS("wEggs_Cmod.rds")

summary(wEggs_Cmod, digits = 4)

# get coefficients
coefs <- coef(wEggs_Cmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# join coefficients with the dataframe for graphing
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# credible interval summaries
ci95 <- posterior_interval(wEggs_Cmod, prob = 0.95, pars = "wEggs_t2")
ci50 <- posterior_interval(wEggs_Cmod, prob = 0.50, pars = "wEggs_t2")
round(ci95, 4)
round(ci50, 6)

# overall plot
wEggs_Counts <- alldata3 %>%
  ggplot(aes(x = wEggs_t2, y = Avg.RASY.Count), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 77.8856, slope = -0.0255, col = "black", lwd = 1) +
  xlab(expression(Neighborhood~competition[t - 2])) +
  ylab(expression(Egg~mass~count[t])) +
  theme_bw(base_size = 12)


# growth rate vs. weighted eggs ----

# partial pooling model
wEggs_Gmod <- stan_glmer(growthrate2 ~ wEggs_t2 + (wEggs_t2|Pond.ID),
                         iter = 10000,
                         na.action = "na.omit",
                         data = alldata2,
                         family = gaussian(),
                         cores = 3,
                         seed = 12345)
summary(wEggs_Gmod, digits = 5)
loo_mod <- loo(wEggs_Gmod) # checks model fit. All are okay.
launch_shinystan(wEggs_Gmod) # model fit okay!

# model saved
saveRDS(wEggs_Gmod, "wEggs_Gmod.rds")

# load model data so no need to rerun
wEggs_Gmod <- readRDS("wEggs_Gmod.rds")

# get coefficients
coefs <- coef(wEggs_Gmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# join coefficients to dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# credible interval summaries
ci95 <- posterior_interval(wEggs_Gmod, prob = 0.95, pars = "wEggs_t2")
ci50 <- posterior_interval(wEggs_Gmod, prob = 0.50, pars = "wEggs_t2")
round(ci95, 4)
round(ci50, 6)

# overall plot
wEggs_Growth <- alldata3 %>%
  ggplot(aes(x = wEggs_t2, y = growthrate2), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5) +
  geom_abline(intercept = 0.02458, slope = -0.00019, col = "black", lwd = 1) +
  xlab(expression(Neighborhood~competition[t - 2])) +
  ylab(expression(Population~growth~rate[t])) +
  scale_x_continuous(limits = c(0, 1500)) +
  theme_bw(base_size = 12)

# count vs. Avg PDSI ----

# fit model
PDSI_Cmod <- stan_glmer(Avg.RASY.Count ~ spfaPDSI + (spfaPDSI|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        seed = 12345)
summary(PDSI_Cmod, digits = 5)

# leave-one-out cross validation
loo_mod <- loo(PDSI_Cmod, k_threshold = 0.7) # checks model fit.

# Shiny app for model fit checks
launch_shinystan(PDSI_Cmod) # model fit good

# model saved
saveRDS(PDSI_Cmod, "PDSI_logCmod.rds")

# load model data if do not want to rerun
PDSI_Cmod <- readRDS("PDSI_Cmod.rds")

# get coefficients
coefs <- coef(PDSI_Cmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# join coefficients with dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# credible intervals
ci95 <- posterior_interval(PDSI_Cmod, prob = 0.95, pars = "spfaPDSI")
round(ci95, 4)

# overall plot
PDSI_Counts <- alldata3 %>%
  ggplot(aes(x = spfaPDSI, y = Avg.RASY.Count), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 71.9667, slope = 4.1812, col = "black", lwd = 1) +
  xlab(expression(Drought~severity~(Mar~-~Nov)[t])) +
  ylab(expression(Egg~mass~count[t])) +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_bw(base_size = 12)


# growth rate vs. avg PDSI ----

# fit model
PDSI_Gmod <- stan_glmer(growthrate2 ~ spfaPDSI + (spfaPDSI|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        seed = 12345)
summary(PDSI_Gmod, digits = 5)
# loo_mod <- loo(PDSI_Gmod, k_threshold = 0.7) # checks model fit. All are okay.
# launch_shinystan(PDSI_Gmod) # model fit okay!

# model saved
saveRDS(PDSI_Gmod, "PDSI_Gmod.rds")

# load model data
PDSI_Gmod <- readRDS("PDSI_Gmod.rds")

# get coefficients
coefs <- coef(PDSI_Gmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# append coefficients to dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# get the 95% credible intervals
ci95 <- posterior_interval(PDSI_Gmod, prob = 0.95, pars = "spfaPDSI")
round(ci95, 4)

# overall plot
PDSI_Growth <- alldata3 %>%
  ggplot(aes(x = spfaPDSI, y = growthrate2), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5) +
  geom_abline(intercept = -0.00241, slope = 0.01595, col = "black", lwd = 1) +
  xlab(expression(Drought~severity~(Mar~-~Nov)[t])) +
  ylab(expression(Population~growth~rate[t])) +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_bw(base_size = 12)


# count vs. winter thaw ----

# fit model
thaw_Cmod <- stan_glmer(Avg.RASY.Count ~ days_thawed.y + (days_thawed.y|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        prior = normal(0, 2, autoscale = TRUE),
                        prior_intercept = normal(0, 5, autoscale = TRUE),
                        cores = 3,
                        seed = 12345)
summary(thaw_Cmod, digits = 4)
loo_mod <- loo(thaw_Cmod) # checks model fit. All are okay.
launch_shinystan(thaw_Cmod) # model fit okay!

# model saved
saveRDS(thaw_Cmod, "thaw_Cmod.rds")
thaw_Cmod <- readRDS("thaw_Cmod.rds")

# get coefficients
coefs <- coef(thaw_Cmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# add coefficients to dataframe
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# 95% credible interval
ci95 <- posterior_interval(thaw_Cmod, prob = 0.95, pars = "days_thawed.y")
round(ci95, 4)

# overall plot
thaw_Counts <- alldata3 %>%
  ggplot(aes(x = days_thawed.y, y = Avg.RASY.Count), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 180.5958, slope = -0.6998, col = "black", lwd = 1) +
  xlab(expression(Winter~thaw~days[t-1])) +
  ylab(expression(Egg~mass~count[t])) +
  # scale_x_continuous(limits = c(-4, 4)) +
  theme_bw(base_size = 12)


# growth rate vs. winter thaw ----

# fit model
thaw_Gmod <- stan_glmer(growthrate2 ~ days_thawed.y + (days_thawed.y|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        seed = 12345)
summary(thaw_Gmod, digits = 5)

# get 95% credible interval
ci95 <- posterior_interval(thaw_Gmod, prob = 0.95, pars = "days_thawed.y")
round(ci95, 4)

# check model fit
# loo_mod <- loo(PDSIlag_Gmod, k_threshold = 0.7) # checks model fit. All are okay.
# launch_shinystan(PDSIlag_Gmod) # model fit okay!

# model saved
saveRDS(thaw_Gmod, "thaw_Gmod.rds")

# model load
thaw_Gmod <- readRDS("thaw_Gmod.rds")

# get coefficients
coefs <- coef(thaw_Gmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# add coefficients to dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# verall plot
thaw_Growth <- alldata3 %>%
  ggplot(aes(x = days_thawed.y, y = growthrate2), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 0.26996, slope = -0.00177, col = "black", lwd = 1) +
  xlab(expression(Winter~thaw~days[t-1])) +
  ylab(expression(Population~growth~rate[t])) +
  theme_bw(base_size = 12)



# count vs. PDSI lag ----

# fit model
PDSIlag_Cmod <- stan_glmer(Avg.RASY.Count ~ spfaPDSIlag + (spfaPDSIlag|Pond.ID),
                           iter = 10000,
                           na.action = "na.omit",
                           data = alldata2,
                           family = gaussian(),
                           prior = normal(0, 2, autoscale = TRUE),
                           prior_intercept = normal(0, 5, autoscale = TRUE),
                           cores = 3,
                           seed = 12345)
summary(PDSIlag_Cmod, digits = 5)

# check model fit
loo_mod <- loo(PDSIlag_Cmod, k_threshold = 0.7) # All are okay. 1 problematic, but was okay.
launch_shinystan(PDSIlag_Cmod) # model fit okay!

# model saved
saveRDS(PDSIlag_Cmod, "PDSIlag_Cmod.rds")

# model data load
PDSIlag_Cmod <- readRDS("PDSIlag_Cmod.rds")

# get coefficients
coefs <- coef(PDSIlag_Cmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# add coefficients to dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# get 95% credible inverval
ci95 <- posterior_interval(PDSIlag_Cmod, prob = 0.95, pars = "spfaPDSIlag")
round(ci95, 5)

# overall plot
PDSIlag_Counts <- alldata3 %>%
  ggplot(aes(x = spfaPDSIlag, y = Avg.RASY.Count), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 72.32583, slope = 3.02360, col = "black", lwd = 1) +
  xlab(expression(Drought~severity~(Mar~-~Nov)[t-2])) +
  ylab(expression(Egg~mass~count[t])) +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_bw(base_size = 12)

# growth rate vs. PDSI lag ----

# fit model
PDSIlag_Gmod <- stan_glmer(growthrate2 ~ spfaPDSIlag + (spfaPDSIlag|Pond.ID),
                           iter = 10000,
                           na.action = "na.omit",
                           data = alldata2,
                           family = gaussian(),
                           cores = 3,
                           seed = 12345)
summary(PDSIlag_Gmod, digits = 5)

# get 95% credible interval
ci95 <- posterior_interval(PDSIlag_Gmod, prob = 0.95, pars = "spfaPDSIlag")
round(ci95, 4)

# check model fit
loo_mod <- loo(PDSIlag_Gmod, k_threshold = 0.7) # checks model fit. All are okay.
launch_shinystan(PDSIlag_Gmod) # model fit okay!

# model saved
saveRDS(PDSIlag_Gmod, "PDSIlag_Gmod.rds")

# load model data
PDSIlag_Gmod <- readRDS("PDSIlag_Gmod.rds")

# get coefficients
coefs <- coef(PDSIlag_Gmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# append coefficients to dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# overall plot
PDSIlag_Growth <- alldata3 %>%
  ggplot(aes(x = spfaPDSIlag, y = growthrate2), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = -0.00500, slope = -0.00526, col = "black", lwd = 1) +
  xlab(expression(Drought~severity~(Mar~-~Nov[t-2]))) +
  ylab(expression(Population~growth~rate[t])) +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_bw(base_size = 12)

# count vs. pond competition ----

# fit model
Dens_Cmod <- stan_glmer(Avg.RASY.Count ~ RASYdens_t2 + (RASYdens_t2|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        seed = 12345)
summary(Dens_Cmod, digits = 5)

# check model fit
loo_mod <- loo(Dens_Cmod) # All are okay.
launch_shinystan(Dens_Cmod) # model fit okay!

# model data saved
saveRDS(Dens_Cmod, "Dens_Cmod.rds")

# model data load
Dens_Cmod <- readRDS("Dens_Cmod.rds")

# get coefficients
coefs <- coef(Dens_Cmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# get 95% credible interval
posterior_interval(Dens_Cmod, regex_pars = "Slope")

# append data to df for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# get 95% credible interval
ci95 <- posterior_interval(Dens_Cmod, prob = 0.95, pars = "RASYdens_t2")
round(ci95, 4)

# overall plot
Dens_Counts <- alldata3 %>%
  ggplot(aes(x = RASYdens_t2, y = Avg.RASY.Count), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 53.4540, slope = 186.4976, col = "black", lwd = 1) +
  xlab(expression(Pond~competition[t-2])) +
  ylab(expression(Egg~mass~count[t])) +
  #scale_x_continuous(limits = c(-4, 4)) +
  theme_bw(base_size = 12)


# growth rate vs. competition ----

# fit model
Dens_Gmod <- stan_glmer(growthrate2 ~ RASYdens_t2 + (RASYdens_t2|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        seed = 12345)
summary(Dens_Gmod, digits = 5)

# get 95% credible interval
ci95 <- posterior_interval(Dens_Gmod, prob = 0.95, pars = "RASYdens_t2")
round(ci95, 4)

# check model fit
loo_mod <- loo(PDSIlag_Gmod) # checks model fit. All are okay.
launch_shinystan(Dens_Gmod) # model fit okay!

# model data saved
saveRDS(Dens_Gmod, "Dens_Gmod.rds")

# model data load
Dens_Gmod <- readRDS("Dens_Gmod.rds")

# get coefficients
coefs <- coef(Dens_Gmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# append coefficients to df for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# overall plot
Dens_Growth <- alldata3 %>%
  ggplot(aes(x = RASYdens_t2, y = growthrate2), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 0.25682, slope = -5.71335, col = "black", lwd = 1) +
  xlab(expression(Pond~competition[t-2])) +
  ylab(expression(Population~growth~rate[t])) +
  theme_bw(base_size = 12)


# count vs. depth ----

# fit model
Depth_Cmod <- stan_glmer(Avg.RASY.Count ~ Depth + (Depth|Pond.ID),
                         iter = 10000,
                         na.action = "na.omit",
                         data = alldata2,
                         family = gaussian(),
                         prior = normal(0, 2, autoscale = TRUE),
                         prior_intercept = normal(0, 5, autoscale = TRUE),
                         cores = 3,
                         seed = 12345)
summary(Depth_Cmod, digits = 4)

# get 95% credible intervals
ci95 <- posterior_interval(Depth_Cmod, prob = 0.95, pars = "Depth")
round(ci95, 4)

# check model fit
loo_mod <- loo(Depth_Cmod, k_threshold = 0.7) # 3 suspicious. All are okay.
launch_shinystan(Depth_Cmod) # model fit okay!

# model data saved
saveRDS(Depth_Cmod, "Depth_logCmod.rds")

# model data load
Depth_Cmod <- readRDS("Depth_Cmod.rds") # also Depth_Cmod.rds for non transformed

# get coefficients
coefs <- coef(Depth_Cmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# append coefficients to df for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# overall plot
Depth_Counts <- alldata3 %>%
  ggplot(aes(x = Depth, y = Avg.RASY.Count), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 44.7179, slope = 0.4461, col = "black", lwd = 1) +
  xlab(expression(Pond~depth[t])) +
  ylab(expression(Egg~mass~count[t])) +
  theme_bw(base_size = 12)


# growth rate vs. depth ----

# fit model
Depth_Gmod <- stan_glmer(growthrate2 ~ Depth + (Depth|Pond.ID),
                         iter = 10000,
                         na.action = "na.omit",
                         data = alldata2,
                         family = gaussian(),
                         cores = 3,
                         seed = 12345)
summary(Depth_Gmod, digits = 5)

# get 95% credible interval
ci95 <- posterior_interval(Depth_Gmod, prob = 0.95, pars = "Depth")
round(ci95, 4)

# check model fit
loo_mod <- loo(Depth_Gmod, k_threshold = 0.7) # checks model fit. All are okay.
launch_shinystan(Depth_Gmod) # model fit okay!

# model data saved
saveRDS(Depth_Gmod, "Depth_Gmod.rds")

# model data load
Depth_Gmod <- readRDS("Depth_Gmod.rds")

# get coefficients
coefs <- coef(Depth_Gmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# append coefficients to dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# overall plot
Depth_Growth <- alldata3 %>%
  ggplot(aes(x = Depth, y = growthrate2), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = -0.03124, slope = 0.00057, col = "black", lwd = 1) +
  xlab(expression(Pond~depth[t])) +
  ylab(expression(Population~growth~rate[t])) +
  theme_bw(base_size = 12)


# count vs. Canopy ----

# fit model
Canopy_Cmod <- stan_glmer(Avg.RASY.Count ~ canopy + (canopy|Pond.ID),
                          iter = 10000,
                          na.action = "na.omit",
                          data = alldata2,
                          prior = normal(0, 2, autoscale = TRUE),
                          prior_intercept = normal(0, 5, autoscale = TRUE),
                          family = gaussian(),
                          cores = 3,
                          seed = 12345)
summary(Canopy_Cmod, digits = 4)

# get 95% credible intervals
ci95 <- posterior_interval(Canopy_Cmod, prob = 0.95, pars = "canopy")
round(ci95, 4)

# check model fit
loo_mod <- loo(Canopy_Cmod, k_threshold = 0.7) # checks model fit. 1 suspicious points. All are okay.
launch_shinystan(Canopy_Cmod) # model fit okay!

# model data saved
saveRDS(Canopy_Cmod, "Canopy_logCmod.rds")

# model data loaded
Canopy_Cmod <- readRDS("Canopy_Cmod.rds")

# get coefficients
coefs <- coef(Canopy_Cmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# append coefficients to df for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# overall plot
Canopy_Counts <- alldata3 %>%
  ggplot(aes(x = canopy, y = Avg.RASY.Count), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = 40.5616, slope = 0.6650, col = "black", lwd = 1) +
  xlab(expression(Pond~canopy~(GSF))) +
  ylab(expression(Egg~mass~count[t])) +
  #scale_x_continuous(limits = c(-4, 4)) +
  theme_bw(base_size = 12)


# growth rate vs. depth ----

# fit model
Canopy_Gmod <- stan_glmer(growthrate2 ~ canopy + (canopy|Pond.ID),
                          iter = 10000,
                          na.action = "na.omit",
                          data = alldata2,
                          family = gaussian(),
                          cores = 3,
                          seed = 12345)
summary(Canopy_Gmod, digits = 4)

# get 95% credible intervals
ci95 <- posterior_interval(Canopy_Gmod, prob = 0.95, pars = "canopy")
round(ci95, 4)

# check model fit
loo_mod <- loo(Depth_Gmod, k_threshold = 0.7) # checks model fit. All are okay.
launch_shinystan(Depth_Gmod) # model fit okay!

# model data saved
saveRDS(Canopy_Gmod, "Canopy_Gmod.rds")

# model data loaded
Canopy_Gmod <- readRDS("Canopy_Gmod.rds")

# get coefficients
coefs <- coef(Canopy_Gmod)$Pond.ID %>%
  rownames_to_column("Pond.ID")
names(coefs) <- c("Pond.ID", "Intercept", "Slope")

# append coefficients to dataframe for plotting
alldata3 <- left_join(alldata2, coefs, by="Pond.ID")

# overall plot
Canopy_Growth <- alldata3 %>%
  ggplot(aes(x = canopy, y = growthrate2), alldata3) +
  geom_point(pch = 21, fill = "grey80", alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, slope = Slope), data = alldata3, color = "grey80", lwd = 0.5, alpha = 0.5) +
  geom_abline(intercept = -0.03124, slope = 0.00057, col = "black", lwd = 1) +
  xlab(expression(Pond~canopy~(GSF))) +
  ylab(expression(Population~growth~rate[t])) +
  theme_bw(base_size = 12)



## multipanel egg count figure ----
library("ggpubr")

allcounts <- ggarrange(
  wEggs_Counts,
  PDSI_Counts +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  PDSIlag_Counts,
  thaw_Counts +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  Dens_Counts,
  Depth_Counts +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  Canopy_Counts ,
  nrow = 4, ncol = 2,
  font.label = list(size = 12),
  align = "v",
  label.x = 0.05,
  labels = "auto")


ggsave(allcounts,
       dpi = 600,
       filename = "counts_figure.pdf",
       height = 9,
       width = 6)


counts_figure <- ggarrange(Dens_Counts,
                           wEggs_Counts,
                           PDSI_Counts,
                           PDSIlag_Counts,
                           thaw_Counts,
                           Depth_Counts,
                           Canopy_Counts,
                           font.label = list(size = 12),
                           labels = c("a", "b", "c", "d", "e", "f", "g"),
                           ncol = 2, nrow = 4)


# multipanel population growth rate figure ----

allgrowth <- ggarrange(
  wEggs_Growth + ylab(expression(Pop.~growth~rate[t])),
  PDSI_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  PDSIlag_Growth + ylab(expression(Pop.~growth~rate[t])),
  thaw_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  Dens_Growth + ylab(expression(Pop.~growth~rate[t])),
  Depth_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  Canopy_Growth+ ylab(expression(Pop.~growth~rate[t])) ,
  nrow = 4, ncol = 2,
  font.label = list(size = 12),
  align = "v",
  label.x = 0.05,
  labels = "auto")


ggsave(allgrowth,
       dpi = 600,
       filename = "growth_figure_v5.pdf",
       height = 9,
       width = 6)