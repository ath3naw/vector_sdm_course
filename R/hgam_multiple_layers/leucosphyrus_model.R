rm(list = ls())
library(tidyverse)
library(terra)
library(geodata)
library(gratia)
library(dplyr)
library(mgcv)
library(sp)
source("R/hgam_multiple_layers/functions_multiple.R")

leuco_data <- read_csv("data/tabular/hgam_multiple_layers/Leucosphyrus_Data.csv")

leuco_data <- leuco_data |> select(-any_of(c(
  "DNA/No DNA",
  "Coordinate system",
  "Coordinates Inferred",
  "Geographical\nregion",
  "GenBank\naccession\nnumbers",
  "Source",
  "Geographical\nregion"
  )))

leuco_data <- leuco_data |> rename(lon_ddm="DDM (Lon)",
                                   lat_ddm="DDM (Lat)", 
                                   lat2=Latitude,
                                   lon2=longitude,
                                   lon= "Longitude (GPS)",
                                   lat= "Latitude (GPS)")

chd <- "°"
chm <- "'"
chs <- "\""
leuco_data$lon_ddm <- fix_missing(leuco_data$lon_ddm)
leuco_data$lat_ddm <- fix_missing(leuco_data$lat_ddm)
leuco_data$lon_ddm <- sapply(leuco_data$lon_ddm, function(x) {
  if (is.na(x)) {
    NA_real_
  } else {
    as.numeric(char2dms(x, chd = chd, chm = chm, chs = chs))
  }
})
leuco_data$lat_ddm <- sapply(leuco_data$lat_ddm, function(x) {
  if (is.na(x)) {
    NA_real_
  } else {
    as.numeric(char2dms(x, chd = chd, chm = chm, chs = chs))
  }
})
leuco_data$lat2 <- sapply(leuco_data$lat2, average_coords)
leuco_data$lon2 <- sapply(leuco_data$lon2, average_coords)
leuco_data$lon <- sapply(leuco_data$lon, clean_coords)
leuco_data$lat <- sapply(leuco_data$lat, clean_coords)
leuco_data$lon <- as.numeric(leuco_data$lon)
leuco_data$lat <- as.numeric(leuco_data$lat)

leuco_data <- leuco_data %>%
  mutate(lonx = coalesce(lon_ddm, lon2, lon), laty = coalesce(lat_ddm, lat2, lat)) %>%
  select(-c(lon_ddm,lon2,lon,lat_ddm,lat2,lat))

leuco_data <- leuco_data[!is.na(leuco_data$lonx),]
leuco_data$Number <- case_when(
  leuco_data$Number == "Found" ~ 1,
  leuco_data$Number=="Nothing" ~ 0,
  !is.na(as.numeric(leuco_data$Number)) & as.numeric(leuco_data$Number)>0 ~ 1,
  TRUE ~ 0
  )
leuco_data$Number <- ifelse(as.numeric(leuco_data$Number)>0, 1, 0)
library(geodata)
wrld <- world(path=".")
plot(wrld, xlim=c(-10, 100), ylim=c(-40, 80), col="light yellow", border="light gray")
points(leuco_data$lonx, leuco_data$laty, col="red", pch=20)
acv <- vect(leuco_data, geom=c("lonx", "laty"), crs="+proj=longlat +datum=WGS84")
ovr <- extract(acv, wrld)
cntr <- ovr$NAME_0
i <- which(is.na(cntr))
i
j <- which(cntr != acv$Country[j])
plot(acv)
lines(wrld,col='blue',lwd=2)

#coords <- vect(df_clean, geom=c("lon", "lat"), crs="+proj=longlat +datum=WGS84")
 
 