# load data
alldata2 <- read.csv("woodfrogdata.csv", header = TRUE)
summary(alldata2)

# libraries
library(rstanarm)
library(tidybayes)
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(ggplot2)
library(cowplot)
library(rstan)
library(RColorBrewer)



# compute neighborhood effect (wEggs) two years prior ----
# add in wEggs(t-2)
alldata2$wEggst2 <- NA
alldata2$Pond.ID <- as.character(alldata2$Pond.ID)
for (i in 3:dim(alldata2)[1]) {
  if (alldata2$Pond.ID[i] == alldata2$Pond.ID[i - 2]) {
    alldata2$wEggst2[i] <- alldata2$wEggs[i - 2]
  } else{
    alldata2$wEggst2[i] <- NA
  }
}

# egg mass count models ----

# set priors
my_priors <- normal(location = (0, 0), scale = (2.5, 2.5))

# normal model
w_eggs_Cmod <- stan_glmer(Avg.RASY.Count ~ wEggst2 + (wEggst2|Pond.ID),
                       iter = 10000,
                       na.action = "na.omit",
                       data = alldata2,
                       family = gaussian(),
                       cores = 3,
                       prior = my_priors,
                       seed = 12345)
summary(w_eggs_Cmod, digits = 6)
ci90 <- posterior_interval(w_eggs_Cmod, prob = 0.90, pars = "wEggst2")
round(ci90, 6)

# log-linear model
alldata2$logAvg.RASY.Count <- log(alldata2$Avg.RASY.Count + 1)

w_eggs_Cmod <- stan_glmer(logAvg.RASY.Count ~ wEggst2 + (wEggst2|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        prior = my_priors,
                        seed = 12345)
summary(w_eggs_Cmod, digits = 6)

draws <- as.data.frame(w_eggs_Cmod)
colnames(draws)[1:2] <- c("a", "b")

wEggs_Counts +
  geom_abline(data = draws, aes(intercept = a, slope = b),
              color = "skyblue", size = 0.2, alpha = 0.25) +
  geom_abline(intercept = coef(w_eggs_Cmod)[1], slope = coef(w_eggs_Cmod)[2],
              color = "skyblue4", size = 1)

# plot it
wEggs_Counts <- alldata2 %>%
  data_grid(wEggst2 = seq_range(wEggst2, n = 101)) %>%
  add_fitted_draws(w_eggs_Cmod, n = 100) %>%
  ggplot(aes(x = wEggst2, y = logAvg.RASY.Count)) +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = .value, group = paste(.draw)), alpha = .1, color = "lightblue") +
  theme_bw(base_size = 20) +
  xlab("Neighborhood competition(t-2)") +
  ylab("Egg mass count") +
  geom_abline(intercept = coef(w_eggs_Cmod)[1], slope = coef(w_eggs_Cmod)[2],
              color = "skyblue4", size = 0.5) +
  ylim(0, 1000) +
  xlim(0, 1500)

wEggs_Counts <- alldata2 %>%
  data_grid(wEggst2 = seq_range(wEggst2, n = 51)) %>%
  add_fitted_draws(w_eggs_Cmod) %>%
  ggplot(aes(x = wEggst2, y = logAvg.RASY.Count)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Neighborhood~competition[t-2])) +
  ylab("Egg mass count") +
  theme(legend.position = "none")

# density(t-2) model
dens_Cmod <- stan_glmer(logAvg.RASY.Count ~ RASYdens_t2 + (1|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        prior = my_priors,
                        seed = 12345)
loo.mod <- loo(dens_Cmod)
summary(dens_Cmod, digits = 6)
ci95 <- posterior_interval(dens_Cmod, prob = 0.95, pars = "RASYdens_t2")
round(ci95, 6)
saveRDS(dens_Cmod, "LogCountdensmod.rds")


Dens_Counts <- alldata2 %>%
  data_grid(RASYdens_t2 = seq_range(RASYdens_t2, n = 51)) %>%
  add_fitted_draws(dens_Cmod) %>%
  ggplot(aes(x = RASYdens_t2, y = logAvg.RASY.Count)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Pond~competition[t-2])) +
  ylab("Egg mass count") +
  theme(legend.position = "none")

# depth model
depth_Cmod <- stan_glmer(Avg.RASY.Count ~ Depth + (Depth|Pond.ID),
                      iter = 10000,
                      na.action = "na.omit",
                      data = alldata2,
                      family = gaussian(),
                      cores = 3,
                      prior = my_priors,
                      seed = 12345)
loo.mod <- loo(depth_Cmod)
summary(depth_Cmod, digits = 6)
ci90 <- posterior_interval(depth_Cmod, prob = 0.90, pars = "Depth")
round(ci90, 6)


Depth_Counts <- alldata2 %>%
  data_grid(Depth = seq_range(Depth, n = 51)) %>%
  add_fitted_draws(depth_Cmod) %>%
  ggplot(aes(x = Depth, y = Avg.RASY.Count)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Pond~depth~(cm))) +
  ylab("Egg mass count") +
  theme(legend.position = "none")


# canopy model
canopy_Cmod <- stan_glmer(Avg.RASY.Count ~ canopy + (canopy|Pond.ID),
                       iter = 10000,
                       na.action = "na.omit",
                       data = alldata2,
                       family = gaussian(),
                       cores = 3,
                       prior = my_priors,
                       seed = 12345)
loo.mod <- loo(canopy_Cmod, k_threshold = 0.7)
summary(canopy_Cmod, digits = 4)
ci90 <- posterior_interval(canopy_Cmod, prob = 0.90, pars = "canopy")
round(ci90, 6)


Canopy_Counts <- alldata2 %>%
  data_grid(canopy = seq_range(canopy, n = 51)) %>%
  add_fitted_draws(canopy_Cmod) %>%
  ggplot(aes(x = canopy, y = Avg.RASY.Count)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Pond~canopy~(GSF))) +
  ylab("Egg mass count") +
  theme(legend.position = "none")

# winter thaw model
winter_Cmod <- stan_glmer(Avg.RASY.Count ~ days_thawed.y + (days_thawed.y|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        prior = my_priors,
                        seed = 12345)
loo.mod <- loo(winter_Cmod)
summary(winter_Cmod, digits = 4)
ci90 <- posterior_interval(winter_Cmod, prob = 0.90, pars = "days_thawed.y")
round(ci90, 4)


Thaw_Counts <- alldata2 %>%
  data_grid(days_thawed.y = seq_range(days_thawed.y, n = 51)) %>%
  add_fitted_draws(winter_Cmod) %>%
  ggplot(aes(x = days_thawed.y, y = Avg.RASY.Count)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value), alpha = 0.5) +
  theme_bw(base_size = 12) +
  xlab(expression(Winter~thaw~days[t-1])) +
  ylab("Egg mass count") +
  theme(legend.position = "none")

# PDSI model
PDSI_Cmod <- stan_glmer(Avg.RASY.Count ~ spfaPDSI + (spfaPDSI|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        prior = my_priors,
                        seed = 12345)
loo.mod <- loo(PDSI_Cmod)
summary(PDSI_Cmod, digits = 4)
ci90 <- posterior_interval(PDSI_Cmod, prob = 0.90, pars = "spfaPDSI")
round(ci90, 4)


PDSI_Counts <- alldata2 %>%
  data_grid(spfaPDSI = seq_range(spfaPDSI, n = 51)) %>%
  add_fitted_draws(PDSI_Cmod) %>%
  ggplot(aes(x = spfaPDSI, y = Avg.RASY.Count)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Avg.~PDSI~(March-Nov))) +
  ylab("Egg mass count") +
  theme(legend.position = "none")

# PDSI lag model
PDSIlag_Cmod <- stan_glmer(Avg.RASY.Count ~ spfaPDSIlag + (spfaPDSIlag|Pond.ID),
                      iter = 10000,
                      na.action = "na.omit",
                      data = alldata2,
                      family = gaussian(),
                      cores = 3,
                      prior = my_priors,
                      seed = 12345)
loo.mod <- loo(PDSIlag_Cmod)
summary(PDSIlag_Cmod, digits = 4)
ci90 <- posterior_interval(PDSIlag_Cmod, prob = 0.90, pars = "spfaPDSIlag")
round(ci90, 4)


PDSIlag_Counts <- alldata2 %>%
  data_grid(spfaPDSIlag = seq_range(spfaPDSIlag, n = 51)) %>%
  add_fitted_draws(PDSIlag_Cmod) %>%
  ggplot(aes(x = spfaPDSIlag, y = Avg.RASY.Count)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Avg.~PDSI~(March-Nov[t-2]))) +
  ylab("Egg mass count") +
  theme(legend.position = "none")

# multipanel egg count figure ----
library("ggpubr")

allcounts <- ggarrange(
          wEggs_Counts,
          PDSI_Counts +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ),
          PDSIlag_Counts +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ),
          Thaw_Counts +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ),
          Dens_Counts,
          Depth_Counts +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ),
          Canopy_Counts +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ),
          nrow = 2, ncol = 4,
          font.label = list(size = 12),
          align = "v",
          common.legend = TRUE,
          legend = "bottom",
          label.x = 0.05,
          labels = "auto")

ggexport(width = 1000, height = 500, allcounts, filename = "counts_figure.jpg")

counts_figure <- ggarrange(Dens_Counts,
                             wEggs_Counts,
                             PDSI_Counts,
                             PDSIlag_Counts,
                             Thaw_Counts,
                             Depth_Counts,
                             Canopy_Counts,
                             font.label = list(size = 12),
                             labels = c("a", "b", "c", "d", "e", "f", "g"),
                             ncol = 2, nrow = 4) %>%
ggexport(allcounts, filename = "counts_figure.jpg")



############################
####### growth rates ----

# growth rate panels ----

# weighted eggs ----
# add in wEggs(t-2)
alldata2$wEggst2 <- NA
alldata2$Pond.ID <- as.character(alldata2$Pond.ID)
for (i in 3:dim(alldata2)[1]) {
  if (alldata2$Pond.ID[i] == alldata2$Pond.ID[i - 2]) {
    alldata2$wEggst2[i] <- alldata2$wEggs[i - 2]
  } else{
    alldata2$wEggst2[i] <- NA
  }
}

# model
w_eggs_Gmod <- stan_glmer(growthrate2 ~ wEggst2 + (wEggst2|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        prior = my_priors,
                        seed = 12345)
summary(w_eggs_Gmod, digits = 6)
ci90 <- posterior_interval(w_eggs_Gmod, prob = 0.90, pars = "wEggst2")
round(ci90, 6)

draws <- as.data.frame(w_eggs_Gmod)
colnames(draws)[1:2] <- c("a", "b")



wEggs_Growth <- alldata2 %>%
  data_grid(wEggst2 = seq_range(wEggst2, n = 51)) %>%
  add_fitted_draws(w_eggs_Gmod) %>%
  ggplot(aes(x = wEggst2, y = growthrate2)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Neighborhood~competition[t-2])) +
  ylab("Population growth rate") +
  theme(legend.position = "none")

# density(t-2) model
dens_Gmod <- stan_glmer(growthrate2 ~ RASYdens_t2 + (RASYdens_t2|Pond.ID),
                      iter = 10000,
                      na.action = "na.omit",
                      data = alldata2,
                      family = gaussian(),
                      cores = 3,
                      prior = my_priors,
                      seed = 12345)
loo.mod <- loo(dens_Gmod)
summary(dens_Gmod, digits = 6)
ci90 <- posterior_interval(dens_Gmod, prob = 0.90, pars = "RASYdens_t2")
round(ci90, 6)


Dens_Growth <- alldata2 %>%
  data_grid(RASYdens_t2 = seq_range(RASYdens_t2, n = 51)) %>%
  add_fitted_draws(dens_Gmod) %>%
  ggplot(aes(x = RASYdens_t2, y = growthrate2)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Pond~competition[t-2])) +
  ylab("Population growth rate") +
  theme(legend.position = "none")

# depth model
depth_Gmod <- stan_glmer(growthrate2 ~ Depth + (Depth|Pond.ID),
                       iter = 10000,
                       na.action = "na.omit",
                       data = alldata2,
                       family = gaussian(),
                       cores = 3,
                       prior = my_priors,
                       seed = 12345)
loo.mod <- loo(depth_Gmod)
summary(depth_Gmod, digits = 6)
ci90 <- posterior_interval(depth_Gmod, prob = 0.90, pars = "Depth")
round(ci90, 6)


Depth_Growth <- alldata2 %>%
  data_grid(Depth = seq_range(Depth, n = 51)) %>%
  add_fitted_draws(depth_Gmod) %>%
  ggplot(aes(x = Depth, y = growthrate2)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Pond~depth~(cm))) +
  ylab("Population growth rate") +
  theme(legend.position = "none")


# canopy model
canopy_Gmod <- stan_glmer(growthrate2 ~ canopy + (canopy|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        prior = my_priors,
                        seed = 12345)
loo.mod <- loo(canopy_Gmod)
summary(canopy_Cmod, digits = 6)
ci90 <- posterior_interval(canopy_Gmod, prob = 0.90, pars = "canopy")
round(ci90, 6)


Canopy_Growth <- alldata2 %>%
  data_grid(canopy = seq_range(canopy, n = 51)) %>%
  add_fitted_draws(canopy_Gmod) %>%
  ggplot(aes(x = canopy, y = growthrate2)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Pond~canopy~(GSF))) +
  ylab("Population growth rate") +
  theme(legend.position = "none")

# winter thaw model
winter_Gmod <- stan_glmer(growthrate2 ~ days_thawed.y + (days_thawed.y|Pond.ID),
                        iter = 10000,
                        na.action = "na.omit",
                        data = alldata2,
                        family = gaussian(),
                        cores = 3,
                        prior = my_priors,
                        seed = 12345)
loo.mod <- loo(winter_Gmod)
summary(winter_Gmod, digits = 5)
ci90 <- posterior_interval(winter_Gmod, prob = 0.90, pars = "days_thawed.y")
round(ci90, 6)


Thaw_Growth <- alldata2 %>%
  data_grid(days_thawed.y = seq_range(days_thawed.y, n = 51)) %>%
  add_fitted_draws(winter_Gmod) %>%
  ggplot(aes(x = days_thawed.y, y = growthrate2)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value), alpha = 0.5) +
  theme_bw(base_size = 12) +
  xlab(expression(Winter~thaw~days[t-1])) +
  ylab("Population growth rate") +
  theme(legend.position = "none")

# PDSI model
PDSI_Gmod <- stan_glmer(growthrate2 ~ spfaPDSI + (spfaPDSI|Pond.ID),
                      iter = 10000,
                      na.action = "na.omit",
                      data = alldata2,
                      family = gaussian(),
                      cores = 3,
                      prior = my_priors,
                      seed = 12345)
loo.mod <- loo(PDSI_Gmod)
summary(PDSI_Gmod, digits = 5)
ci90 <- posterior_interval(PDSI_Gmod, prob = 0.90, pars = "spfaPDSI")
round(ci90, 5)


PDSI_Growth <- alldata2 %>%
  data_grid(spfaPDSI = seq_range(spfaPDSI, n = 51)) %>%
  add_fitted_draws(PDSI_Gmod) %>%
  ggplot(aes(x = spfaPDSI, y = growthrate2)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Avg.~PDSI~(March-Nov))) +
  ylab("Population growth rate") +
  theme(legend.position = "none")

# PDSI lag model
PDSIlag_Gmod <- stan_glmer(growthrate2 ~ spfaPDSIlag + (spfaPDSIlag|Pond.ID),
                         iter = 10000,
                         na.action = "na.omit",
                         data = alldata2,
                         family = gaussian(),
                         cores = 3,
                         prior = my_priors,
                         seed = 12345)
loo.mod <- loo(PDSIlag_Gmod)
summary(PDSIlag_Gmod, digits = 5)
ci90 <- posterior_interval(PDSIlag_Gmod, prob = 0.90, pars = "spfaPDSIlag")
round(ci90, 5)


PDSIlag_Growth <- alldata2 %>%
  data_grid(spfaPDSIlag = seq_range(spfaPDSIlag, n = 51)) %>%
  add_fitted_draws(PDSIlag_Gmod) %>%
  ggplot(aes(x = spfaPDSIlag, y = growthrate2)) +
  scale_fill_brewer(palette = "Blues") +
  geom_point(data = alldata2, pch = 21, fill = "grey80", alpha = 0.5) +
  stat_lineribbon(aes(y = .value)) +
  theme_bw(base_size = 12) +
  xlab(expression(Avg.~PDSI~(March-Nov[t-2]))) +
  ylab("Population growth rate") +
  theme(legend.position = "none")

# multipanel egg count figure ----
library("ggpubr")

allgrowth <- ggarrange(
  wEggs_Growth,
  PDSI_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  PDSIlag_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  Thaw_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  Dens_Growth,
  Depth_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  Canopy_Growth +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank() ),
  nrow = 2, ncol = 4,
  font.label = list(size = 12),
  align = "v",
  common.legend = TRUE,
  legend = "bottom",
  label.x = 0.05,
  labels = "auto")

ggexport(width = 1000, height = 500, allcounts, filename = "counts_figure.jpg")

counts_figure <- ggarrange(Dens_Counts,
                           wEggs_Counts,
                           PDSI_Counts,
                           PDSIlag_Counts,
                           Thaw_Counts,
                           Depth_Counts,
                           Canopy_Counts,
                           font.label = list(size = 12),
                           labels = c("a", "b", "c", "d", "e", "f", "g"),
                           ncol = 2, nrow = 4) %>%
  ggexport(allcounts, filename = "counts_figure.jpg")
