# This script produces maps displaying Yale Myers Forest wood frog ponds.
# Written by A. Z. Andis Arietta on 04 March 2021

library(tidyverse)
library(ggmap)
library(rgdal)
library(raster)
library(GGally)
library(network)

# Import pond coordinates
pond_coords <- read.csv("./YM_Pond_Coords.csv")

# Get a basemap image from stamen engine
basemap <- get_stamenmap(bbox = c(left = -72.19, bottom = 41.9, right = -72.08, top = 42.0), zoom = 12, mayptype = "terrain")

# Import shapefile for Yale Myers Forest boundary (dowloaded from YMF website)
YMF_shape <- readOGR(dsn = "./YMF_shapefiles", layer = "YMF_Boundary_Full") # Read in the shapefile
YMF_shape_ll <- spTransform(YMF_shape, "+proj=longlat +datum=WGS84") # Correct the coordinate system from UTM to latlong
YMF_shape_df <- fortify(YMF_shape_ll) # Coerce into a dataframe to play nicely with ggplot

# Plot the basemap with YMF boundary
YMFmap <- ggmap(basemap) +
  geom_polygon(data = YMF_shape_df, aes(x = long, y = lat, group = group), fill = NA, col = "grey40", size = 1)

# Create the network edges. To do this, we create a distance matrix, exclude edges more than 500m, and then convert this to a network for the graph
Points <- as.matrix(pond_coords[,3:2]) # get pond coordinates into matrix
distmat <- pointDistance(Points, Points, lonlat = TRUE, allpairs = TRUE) # create a full pair-wise matrix
rownames(distmat) <- pond_coords$Name
colnames(distmat) <- pond_coords$Name
xy <- t(combn(colnames(distmat), 2)) # Get the top triangle of the matrix
DistEdges <- data.frame(xy, dist = distmat[xy]) %>% # Make the pairwise values with distances into a dataframe
  filter(dist <= 500) %>% # Exclude pairs with more than 500m distance
  dplyr::select(-dist)

distnet <- network(DistEdges) # Convert pair-wise distance dataframe into a network object

row.names(pond_coords) <- pond_coords$Name # index the pond coordinate data by pond names

distnet %v% "lat" <- pond_coords[network.vertex.names(distnet), "Lat"] # Assign the coodinates to the network edges
distnet %v% "lon" <- pond_coords[network.vertex.names(distnet), "Long"] # Assign the coodinates to the network edges

### Create the maps
## Map nodes and edges
ggnetworkmap(YMFmap, distnet, size = 0, great.circles = FALSE, segment.color = "darkred") +
    geom_point(data = pond_coords, aes(y = Lat, x = Long), pch = 21, col = "grey20", size = 2, fill = "grey50")

## Map nodes only
YMFmap +
  geom_point(data = pond_coords, aes(y = Lat, x = Long), pch = 21, col = "grey20", size = 2, fill = "grey50")

## Overview map
library(sf)
library(rnaturalearth)
library(ggspatial)

world <- ne_countries(scale = "medium", returnclass = "sf") # Get country boundaries
us_states <- map_data("state") # Get state boundaries

# Plot the overview map with bounding box for the forest inset.
ggplot(data = world) +
  # geom_sf(color = "grey20", fill = "grey80") +
  geom_polygon(data = us_states, aes(x = long, y = lat, group = group), color = "grey20", fill = "grey80") +
  annotation_scale(location = "tl", width_hint = 0.25) +
  annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.00, "in"), pad_y = unit(0.35, "in"), style = north_arrow_minimal) +
  coord_sf(xlim = c(-74.6,-68.9), ylim = c(40.1,43.8), expand = F) +
  geom_point(data = pond_coords, aes(y = Lat, x = Long)) +
  geom_rect(xmin = -72.19, xmax = -72.08, ymin = 41.9, ymax = 42.0, col = "firebrick", alpha = 1, size = 1) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
(-74.6+-68.9)/2
(40.1+43.8)/2


# Helpful sites:
# https://cran.microsoft.com/snapshot/2016-01-19/web/packages/GGally/vignettes/ggnetworkmap.html
# https://rstudio-pubs-static.s3.amazonaws.com/298685_7ca83d01093a4ed79479197945d3783d.html
# https://www.r-bloggers.com/2018/05/three-ways-of-visualizing-a-graph-on-a-map/
# https://www.littlemissdata.com/blog/maps
