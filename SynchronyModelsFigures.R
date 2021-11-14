# Code to compute spatial and temporal synchrony
# written by F.E. Rowland
# January 2021

# libraries
library(ncf)
library(ape)
library(cowplot)
library(reshape2)
library(viridis)
library(ggrepel)
library(grid)
library(gridExtra)
library(ggpubr)
library(lmap)

# get distance matrix between ponds ----
ponds2$name <- ponds2$Pond.ID
#coords <- ponds[c(1:16, 18:23, 25:32, 34:50, 52, 56:58, 60:62, 64:66, 69:71), c(5:7)] # take out the ponds with no data
coords <- ponds2[ , c(1, 5:6)] # subset Pond.ID, lat, lon

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.

  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }

  n.geopoints <- nrow(df.geopoints)

  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints

  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})

  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name

  return(mat.distances)
}

distance.mat.m <- GeoDistanceInMetresMatrix(coords)

# inverse the distance matrix (so closer = smaller)
inv.dist <- 1/distance.mat.m

# have diag spots be zero, since that is the pond effect on itself
diag(inv.dist) <- 0

# use ncf package and spline.correlog to detect synchrony ----

# data
countmat <- read.csv("wide_eggs.csv", header = TRUE)
countmat_med <- read.csv("wide_eggs_med.csv", header = TRUE)
coords <- read.csv("pond_coordinates3.csv", header = TRUE)
coords2 <- subset(coords, eliminate == 0)
x <- coords2$lat
y <- coords2$lon
z <- countmat[,-1]
z.med <- countmat_med[,-1]


# Spatial correlogram
fit1 <- correlog.nc(x, y, z, w = NULL, increment = 0.5, resamp = 5000,
                    na.rm = TRUE, latlon = TRUE, quiet = FALSE)
summary(fit1)
plot(fit1)

# with median distance only
fit1 <- correlog.nc(x, y, z.med, w = NULL, increment = 0.5, resamp = 5000,
                    na.rm = TRUE, latlon = TRUE, quiet = FALSE)
summary(fit1)
plot(fit1)

# noncentered (Mantel) correlogram shows 2016 as outlier
fit1 <- correlog.nc(x = x, y = y, z = z, increment = 0.5, na.rm = TRUE, resamp = 1000)  ## Not run: plot(fit1)

# plot spline correlelogram USE ----
fit2 <- spline.correlog(x = x, y = y, z = z, latlon = TRUE,
                        resamp = 5000, na.rm = TRUE)
plot(fit2, text = TRUE)
summary(fit2)

# better plot of spline
x <- as.vector(fit2$real$predicted$x)
y <- as.vector(fit2$real$predicted$y)
newdata <- data.frame(x, y)

lCI <- fit2$boot$boot.summary$predicted$y["0.025",]
uCI <- fit2$boot$boot.summary$predicted$y["0.975",]
xx <- c(x,rev(x))
yy <- c(lCI,rev(uCI))

shapedata <- data.frame(lCI, uCI)

# Create spatial correlation plot in ggplot (use in pub) ----
correlog.plot <- ggplot(aes(x = x, y = y), data = newdata) +
  xlim(0, 6) +
  ylim(-1, 1) +
  xlab("Distance (km)") +
  ylab("Moran's I") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw(base_size = 20) +
  geom_ribbon(data = shapedata, aes(ymin = lCI, ymax = uCI), alpha = 0.5, fill = "grey50") +
  geom_line(data = newdata, aes(x = x, y = y))

# spatial plot
spatial.plot(x, y, z, ctr = TRUE)


# Calculate Moran's I ----

# make sure names match up for calcu
# eggs3$Pond.ID[!(eggs3$Pond.ID %in% inv.dist$name)]

# calculate distances other ways
pond.dist <- as.matrix(dist(cbind(ponds2$lon, ponds2$lat)))
pond.dist.inv <- 1/pond.dist
diag(pond.dist.inv) <- 0
pond.dist.inv[is.infinite(pond.dist.inv)] <- 0

eggs_pond <- eggs3 %>%
  dplyr::group_by(Pond.ID) %>%
  dplyr::summarize(avg = mean(Avg.RASY.Count, na.rm = TRUE),
                   sd = sd(Avg.RASY.Count, na.rm = TRUE),
                   CV = sd(Avg.RASY.Count, na.rm = TRUE)/mean(Avg.RASY.Count, na.rm = TRUE)*100)

# calculate Moran I
Moran.I(eggs_pond$avg, pond.dist.inv, scaled = TRUE, na.rm = TRUE) # p = 0.0148, Moran's I = 0.0707 (weak)
Moran.I(eggs_pond$sd, pond.dist.inv, scaled = TRUE, na.rm = TRUE) # p = 0.0024, Moran's I = 0.0927 (weak)
Moran.I(eggs_pond$CV, pond.dist.inv, na.rm = TRUE) # p = 0.7409, Moran's I = -0.0274 (none)
Moran.I(eggs_pond$avg_gr, pond.dist.inv, na.rm = TRUE) # p = 0.938, Moran's I = -0.0187 (none)

# calculate Moran's I by Avg Egg masses per year ----
eggs_red <- eggs3[,c(1,3,5,8)]

# fill in missing datas
eggs_red2 <- complete(eggs_red, Year=2000:2020,
                      fill = list(Avg.RASY.Count = NA))

# change data format from long to wide
egg_wide2 <- dcast(eggs_red2, Pond.ID ~ Year, value.var="Avg.RASY.Count", mean, na.rm = TRUE)

egg_wide2[egg_wide2 =="NaN"] <- "NA"
egg_wide3 <- egg_wide2 %>% mutate_if(is.character,as.numeric)

Moran2000 <- Moran.I(egg_wide3$"2000", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2001 <- Moran.I(egg_wide3$"2001", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2002 <- Moran.I(egg_wide3$"2002", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2003 <- Moran.I(egg_wide3$"2003", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2004 <- Moran.I(egg_wide3$"2004", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2005 <- Moran.I(egg_wide3$"2005", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2006 <- Moran.I(egg_wide3$"2006", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2007 <- Moran.I(egg_wide3$"2007", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2008 <- Moran.I(egg_wide3$"2008", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2009 <- Moran.I(egg_wide3$"2009", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2010 <- Moran.I(egg_wide3$"2010", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2011 <- Moran.I(egg_wide3$"2011", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2012 <- Moran.I(egg_wide3$"2012", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2013 <- Moran.I(egg_wide3$"2013", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2014 <- Moran.I(egg_wide3$"2014", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2015 <- Moran.I(egg_wide3$"2015", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2016 <- Moran.I(egg_wide3$"2016", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2017 <- Moran.I(egg_wide3$"2017", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2018 <- Moran.I(egg_wide3$"2018", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2019 <- Moran.I(egg_wide3$"2019", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2020 <- Moran.I(egg_wide3$"2020", pond.dist.inv, scaled = TRUE, na.rm = TRUE)

# manually entered data from above because... Luddite
Moran_sum <- read.csv("Moran_I.csv", header = TRUE)

# Bonferroni adjustement
Moran_sum$Bonf.p <- p.adjust(Moran_sum$p.value, method="bonferroni")

# Holm-Bonferroni method
Moran_sum$HBonf.p <- p.adjust(Moran_sum$p.value, method = "holm")

round_function <- function(x, digits) {
  # round all numeric variables
  # x: data frame
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

Moran_round <- round_function(Moran_sum, 3)

# create plot of Moran's I over time ----

Moran_plot <- ggplot(aes(
  x = Year,
  y = Moran.obs,
  group = p.value,
  label = p.value
), data = Moran_sum) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Moran's I") +
  ylim(-1, 1) +
  geom_point(
    aes(
      x = Year,
      y = Moran.obs,
      fill = ifelse(p.value < 0.05, '≤ 0.05', '> 0.05')
    ),
    pch = 21,
    size = 4
  ) +
  scale_fill_manual(values = c("grey50", "red")) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")
# labs(fill = "P-value")

## make multipanel plot of Moran's I and ncf plot ----

autocorr_figure <- ggarrange(Moran_plot, correlog.plot,
                             font.label = list(size = 20),
                             labels = c("a", "b"),
                             ncol = 1, nrow = 2,
                             align = "hv",
                             width = 400, height = 600) %>%
  ggexport(filename = "autocorr_figure.jpeg")
autocorr_figure

# combine Moran plot and correlog plot
grid.arrange(Moran_plot, correlog.plot, ncol = 1)

# calculate Moran's I by *egg mass density* per year ----

eggs3$dens <- (eggs3$Avg.RASY.Count/eggs3$Area)
eggs_red <- eggs3[,c(1,3,5,8,13)]

# make long format for density
egg_wide2 <- dcast(eggs_red, Pond.ID ~ Year, value.var=c("dens"), mean, na.rm = TRUE)

egg_wide2[egg_wide2 =="NaN"] <- "NA"
egg_wide3 <- egg_wide2 %>% mutate_if(is.character,as.numeric)

Moran2000 <- Moran.I(egg_wide3$"2000", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2001 <- Moran.I(egg_wide3$"2001", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2002 <- Moran.I(egg_wide3$"2002", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2003 <- Moran.I(egg_wide3$"2003", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2004 <- Moran.I(egg_wide3$"2004", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2005 <- Moran.I(egg_wide3$"2005", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2006 <- Moran.I(egg_wide3$"2006", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2007 <- Moran.I(egg_wide3$"2007", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2008 <- Moran.I(egg_wide3$"2008", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2009 <- Moran.I(egg_wide3$"2009", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2010 <- Moran.I(egg_wide3$"2010", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2011 <- Moran.I(egg_wide3$"2011", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2012 <- Moran.I(egg_wide3$"2012", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2013 <- Moran.I(egg_wide3$"2013", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2014 <- Moran.I(egg_wide3$"2014", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2015 <- Moran.I(egg_wide3$"2015", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2016 <- Moran.I(egg_wide3$"2016", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2017 <- Moran.I(egg_wide3$"2017", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2018 <- Moran.I(egg_wide3$"2018", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2019 <- Moran.I(egg_wide3$"2019", pond.dist.inv, scaled = TRUE, na.rm = TRUE)
Moran2020 <- Moran.I(egg_wide3$"2020", pond.dist.inv, scaled = TRUE, na.rm = TRUE)

# manually entered data from above because... Luddite
Moran_sum <- read.csv("Moran_I.csv", header = TRUE)

round_function <- function(x, digits) {
  # round all numeric variables
  # x: data frame
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

Moran_round <- round_function(Moran_sum, 3)

# create plot of Moran's I density over time ----
dev.off()
MoranDens_plot <- ggplot(aes(x = Year, y = dens_Moran.obs, group = dens_p.value, label = dens_p.value), data = Moran_sum) +
  geom_point(aes(x = Year, y = dens_Moran.obs, fill = dens_p.value), pch = 21, size = 4) +
  ylab("Egg Density Moran's I") +
  scale_fill_viridis_c(breaks = seq(0, 1, 0.2),
                       option = "magma",
                       direction = -1,
                       limits = c(0, 1)) +
  theme_bw(base_size = 16) +
  ylim(-0.11, 0.11) +
  geom_label_repel(
    data = subset(Moran_sum, dens_p.value <= 0.10),
    segment.color = "grey50",
    segment.size = 0.2,
    box.padding = 1,
    direction = "y") +
  labs(fill = "P-value")

Moran_plot <- ggplot(aes(x = Year, y = Moran.obs, group = p.value, label = p.value), data = Moran_sum) +
  geom_point(aes(x = Year, y = Moran.obs, fill = p.value), pch = 21, size = 4) +
  ylab("Moran's I") +
  ylim(-0.11, 0.11) +
  scale_fill_viridis_c(breaks = seq(0, 1, 0.2),
                       option = "magma",
                       direction = -1,
                       limits = c(0, 1)) +
  theme_bw(base_size = 16) +
  xlab(NULL) +
  geom_label_repel(
    data = subset(Moran_round, p.value <= 0.10),
    segment.color = "grey50",
    segment.size = 0.2,
    box.padding = 1,
    direction = "y") +
  labs(fill = "P-value")


columngrid <- plot_grid(
  Moran_plot + theme(legend.position = "none"),
  MoranDens_plot + theme(legend.position = "none"),
  labels = c('A', 'B'),
  ncol = 2,
  align="hv"
)

legend <- get_legend(
  # create some space to the left of the legend
  Moran_plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
plot_grid(columngrid, legend, rel_widths = c(2, .4))



