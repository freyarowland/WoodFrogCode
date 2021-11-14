# Code to compute spatial and temporal synchrony
# written by F.E. Rowland
# January 2021

# libraries
library(ncf)
library(ape)

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
