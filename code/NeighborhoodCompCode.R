# Code for calculated weighted neighborhood competition based
# on egg mass counts of ponds within 500 meters

# Code by F.E. Rowland
# Last edit January 2022

# three steps
# 1) make all ponds > 500 m away NA
# 2) calculate exp(-θd(ij))*areaj
# 3) sum together

# libraries
library(reshape2)
library(stringr)
library(dplyr)


# load distance matrix (see spatial synchrony code)
distance.mat.m <- read.csv("distance.mat.csv", header = TRUE)

# eliminate ponds > 500 m away
distance.mat.500 <- ifelse(test = distance.mat.m <= 500, distance.mat.m, NA)
diag(distance.mat.m) <- NA

# create for loop for connectivity ----

conn.mat <- matrix(data = NA, nrow = length(distance.mat.500[,2]), ncol = length(distance.mat.500[,2]))
for (r in 1:nrow(distance.mat.500)) {
  for (c in 1:ncol(distance.mat.500)){
    conn.mat[r, c] <- exp(-(1/500)*distance.mat.500[r, c])
  }
}

# Count data - turn all ponds to 1 ----
conn.mat <- ifelse(test = distance.mat.m <= 500, 1, NA)
diag(conn.mat) <- NA

# Nearest neighbor (nn) measurment "nn_mat"
# return only smallest value in each column ----
# load data
nn.mat <- read.csv("nn_mat.csv", header = TRUE)
rownames(nn.mat) <- rownames(distance.mat.m)
colnames(nn.mat) <- colnames(distance.mat.m)

# add in row and column names ----
rownames(conn.mat) <- rownames(distance.mat.500)
colnames(conn.mat) <- colnames(distance.mat.500)

# make diagonal NA
diag(conn.mat) <- NA

# multiply entry by pond average egg masses
conn.df <- data.frame(conn.mat)
nn.df <- data.frame(nn.mat)



# changing data so column is RASY count within that year ----

# fill in missing dates with NA
eggs_long <- complete(eggs2, Year=2000:2020,
                      fill = list(Avg.RASY.Count = NA))

# recast data from long to wide format
egg_wide2 <- dcast(eggs_long, Pond.ID ~ Year, value.var="Avg.RASY.Count", mean, na.rm = TRUE)

# some clean up
egg_wide2[egg_wide2 =="NaN"] <- "NA"
egg_wide3 <- egg_wide2 %>% mutate_if(is.character,as.numeric)

# example of how to use
# Moran2000 <- Moran.I(egg_wide3$"2000", pond.dist.inv, scaled = TRUE, na.rm = TRUE)

# make different years of egg mass counts
eggsavg <- eggs11 %>%
  dplyr::group_by(Pond.ID) %>%
  dplyr::summarize(avg = mean(Avg.RASY.Count, na.rm = TRUE))

# multiple connectivity matrix * RASY counts ----
conn_x_eggs <- eggsavg$avg * conn.df
conn_2000_eggs <- egg_wide3$"2000" * conn.df
conn_2001_eggs <- egg_wide3$"2001" * conn.df
conn_2002_eggs <- egg_wide3$"2002" * conn.df
conn_2003_eggs <- egg_wide3$"2003" * conn.df
conn_2004_eggs <- egg_wide3$"2004" * conn.df
conn_2005_eggs <- egg_wide3$"2005" * conn.df
conn_2006_eggs <- egg_wide3$"2006" * conn.df
conn_2007_eggs <- egg_wide3$"2007" * conn.df
conn_2008_eggs <- egg_wide3$"2008" * conn.df
conn_2009_eggs <- egg_wide3$"2009" * conn.df
conn_2010_eggs <- egg_wide3$"2010" * conn.df
conn_2011_eggs <- egg_wide3$"2011" * conn.df
conn_2012_eggs <- egg_wide3$"2012" * conn.df
conn_2013_eggs <- egg_wide3$"2013" * conn.df
conn_2014_eggs <- egg_wide3$"2014" * conn.df
conn_2015_eggs <- egg_wide3$"2015" * conn.df
conn_2016_eggs <- egg_wide3$"2016" * conn.df
conn_2017_eggs <- egg_wide3$"2017" * conn.df
conn_2018_eggs <- egg_wide3$"2018" * conn.df
conn_2019_eggs <- egg_wide3$"2019" * conn.df
conn_2020_eggs <- egg_wide3$"2020" * conn.df

conn_2000_eggs <- egg_wide3$"2000" * nn.df
conn_2001_eggs <- egg_wide3$"2001" * nn.df
conn_2002_eggs <- egg_wide3$"2002" * nn.df
conn_2003_eggs <- egg_wide3$"2003" * nn.df
conn_2004_eggs <- egg_wide3$"2004" * nn.df
conn_2005_eggs <- egg_wide3$"2005" * nn.df
conn_2006_eggs <- egg_wide3$"2006" * nn.df
conn_2007_eggs <- egg_wide3$"2007" * nn.df
conn_2008_eggs <- egg_wide3$"2008" * nn.df
conn_2009_eggs <- egg_wide3$"2009" * nn.df
conn_2010_eggs <- egg_wide3$"2010" * nn.df
conn_2011_eggs <- egg_wide3$"2011" * nn.df
conn_2012_eggs <- egg_wide3$"2012" * nn.df
conn_2013_eggs <- egg_wide3$"2013" * nn.df
conn_2014_eggs <- egg_wide3$"2014" * nn.df
conn_2015_eggs <- egg_wide3$"2015" * nn.df
conn_2016_eggs <- egg_wide3$"2016" * nn.df
conn_2017_eggs <- egg_wide3$"2017" * nn.df
conn_2018_eggs <- egg_wide3$"2018" * nn.df
conn_2019_eggs <- egg_wide3$"2019" * nn.df
conn_2020_eggs <- egg_wide3$"2020" * nn.df

connect_final <- data.frame(dplyr::as_tibble(rowSums(x = conn_x_eggs, na.rm = TRUE)))
connect_2000 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2000_eggs, na.rm = TRUE)))
connect_2001 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2001_eggs, na.rm = TRUE)))
connect_2002 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2002_eggs, na.rm = TRUE)))
connect_2003 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2003_eggs, na.rm = TRUE)))
connect_2004 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2004_eggs, na.rm = TRUE)))
connect_2005 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2005_eggs, na.rm = TRUE)))
connect_2006 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2006_eggs, na.rm = TRUE)))
connect_2007 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2007_eggs, na.rm = TRUE)))
connect_2008 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2008_eggs, na.rm = TRUE)))
connect_2009 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2009_eggs, na.rm = TRUE)))
connect_2010 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2010_eggs, na.rm = TRUE)))
connect_2011 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2011_eggs, na.rm = TRUE)))
connect_2012 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2012_eggs, na.rm = TRUE)))
connect_2013 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2013_eggs, na.rm = TRUE)))
connect_2014 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2014_eggs, na.rm = TRUE)))
connect_2015 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2015_eggs, na.rm = TRUE)))
connect_2016 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2016_eggs, na.rm = TRUE)))
connect_2017 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2017_eggs, na.rm = TRUE)))
connect_2018 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2018_eggs, na.rm = TRUE)))
connect_2019 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2019_eggs, na.rm = TRUE)))
connect_2020 <- data.frame(dplyr::as_tibble(rowSums(x = conn_2020_eggs, na.rm = TRUE)))

# connect pond.ID to data ----
connect_2000$Pond.ID <- coords$name
connect_2001$Pond.ID <- coords$name
connect_2002$Pond.ID <- coords$name
connect_2003$Pond.ID <- coords$name
connect_2004$Pond.ID <- coords$name
connect_2005$Pond.ID <- coords$name
connect_2006$Pond.ID <- coords$name
connect_2007$Pond.ID <- coords$name
connect_2008$Pond.ID <- coords$name
connect_2009$Pond.ID <- coords$name
connect_2010$Pond.ID <- coords$name
connect_2011$Pond.ID <- coords$name
connect_2012$Pond.ID <- coords$name
connect_2013$Pond.ID <- coords$name
connect_2014$Pond.ID <- coords$name
connect_2015$Pond.ID <- coords$name
connect_2016$Pond.ID <- coords$name
connect_2017$Pond.ID <- coords$name
connect_2018$Pond.ID <- coords$name
connect_2019$Pond.ID <- coords$name
connect_2020$Pond.ID <- coords$name

# make this into one big dataframe to merge with other data ----
dd1 <- left_join(connect_2000, connect_2001, by = "Pond.ID")
dd2 <- left_join(dd1, connect_2002, by = "Pond.ID")
dd3 <- left_join(dd2, connect_2003, by = "Pond.ID")
dd4 <- left_join(dd3, connect_2004, by = "Pond.ID")
dd5 <- left_join(dd4, connect_2005, by = "Pond.ID")
dd5n <- dd5 %>%
  dplyr::rename(
    w2000 = value.x,
    w2001 = value.y,
    w2002 = value.x.x,
    w2003 = value.y.y,
    w2004 = value.x.x.x,
    w2005 = value.y.y.y
  )
dd6 <- left_join(dd5n, connect_2006, by = "Pond.ID")
dd7 <- left_join(dd6, connect_2007, by = "Pond.ID")
dd8 <- left_join(dd7, connect_2008, by = "Pond.ID")
dd9 <- left_join(dd8, connect_2009, by = "Pond.ID")
dd10 <- left_join(dd9, connect_2010, by = "Pond.ID")
dd11 <- left_join(dd10, connect_2011, by = "Pond.ID")
dd11n <- dd11 %>%
  dplyr::rename(
    w2006 = value.x,
    w2007 = value.y,
    w2008 = value.x.x,
    w2009 = value.y.y,
    w2010 = value.x.x.x,
    w2011 = value.y.y.y
  )
dd12 <- left_join(dd11n, connect_2012, by = "Pond.ID")
dd13 <- left_join(dd12, connect_2013, by = "Pond.ID")
dd14 <- left_join(dd13, connect_2014, by = "Pond.ID")
dd15 <- left_join(dd14, connect_2015, by = "Pond.ID")
dd16 <- left_join(dd15, connect_2016, by = "Pond.ID")
dd17 <- left_join(dd16, connect_2017, by = "Pond.ID")
dd17n <- dd17 %>%
  dplyr::rename(
    w2012 = value.x,
    w2013 = value.y,
    w2014 = value.x.x,
    w2015 = value.y.y,
    w2016 = value.x.x.x,
    w2017 = value.y.y.y
  )
dd18 <- left_join(dd17n, connect_2018, by = "Pond.ID")
dd19 <- left_join(dd18, connect_2019, by = "Pond.ID")
dd20 <- left_join(dd19, connect_2020, by = "Pond.ID")
dd20n <- dd20 %>%
  dplyr::rename(
    w2018 = value.x,
    w2019 = value.y,
    w2020 = value
  )

# recast from long to wide ----
dd_final <- melt(dd20n, id.vars = "Pond.ID",
                 variable.name = "wYear",
                 value.name = "nnEggs")
head(dd_final)



# extract year
numextract <- function(string){
  str_extract(string, "\\-*\\d+\\.*\\d*")
}
dd_final$Year <- numextract(dd_final$wYear)
wDDfinal$Year <- numextract(wDDfinal$wYear)
nn_final$Year <- numextract(nn_final$wYear)

# save data

# weighted density dependence
write.csv(x = wDDfinal, file="wDDdata.csv")

# nearest neighbor
write.csv(nn_final, "nn_final.csv")
