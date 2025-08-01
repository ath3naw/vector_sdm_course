rm(list = ls())
library(tidyverse)
library(terra)
library(geodata)
library(gratia)
library(dplyr)
library(mgcv)
library(sp)
source("R/leucosphyrus_model/leuco_functions.R")

leuco_data <- read_csv("data/Leucosphyrus_Data.csv")

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
leuco_data$Species <- str_to_title(leuco_data$Species)
leuco_data$Species[leuco_data$Species=="Leucosphyrus Group"] <- "Leucosphyrus Subgroup"
library(geodata)
wrld <- world(path=".")
plot(wrld, xlim=c(50, 155), ylim=c(-20, 40), col="light yellow", border="light gray")
points(leuco_data$lonx, leuco_data$laty, col="red", pch=20)
leuco_vector <- vect(leuco_data, geom=c("lonx", "laty"), crs="+proj=longlat +datum=WGS84")
mis_coords <- extract(wrld, leuco_vector)
cntr <- mis_coords$NAME_0
i <- which(is.na(cntr))
i
plot(wrld, xlim=c(50, 155), ylim=c(-20, 40), col="light yellow", border="light gray")
points(leuco_data$lonx[i], leuco_data$laty[i], col="red", pch=20)
off <- cbind(mis_coords[i,], leuco_vector$Country[i])

# is just a little off, should be an island
j <- which(cntr != leuco_vector$Country)
m <- cbind(cntr[j], leuco_vector$Country[j])
colnames(m) <- c("polygons", "data")
m

leuco_data$Species <- str_to_title(leuco_data$Species)

library(tidyverse)
library(terra)
library(geodata)



# get extent map for our basic analyses

leuco_dist_map <- gadm(
  country = c("THA", "MM", "IN", "KHM", "IDN", "VN", "MYS","LA", "CN", "PH", "NPL", "BTN", "BD"),
  level = 0,
  path = "data/downloads"
) # slow to run though not large

plot(leuco_dist_map)

# let's download some data

# bioclimactic variables from worldclim
# https://worldclim.org/data/bioclim.html
# bio <- worldclim_global(var = "bio", res = 0.5, path = "data/downloads")
# se_asia_extent <- ext(85, 140, -15, 35)


library(terra)

# Path to the folder with the .tif files
bioclim_path <- "grids"

# Stack all 19 .tif files
bioclim_stack <- rast(file.path(bioclim_path, paste0("wc2.1_30s_bio_", 1:19, ".tif")))

seasia_ext <- ext(60, 150, -15, 55)  # adjust as needed

bioclim_leuco_map <- crop(bioclim_stack, seasia_ext)
plot(bioclim_leuco_map[[1]])


# generate mask of our area of interest so we can use it to process other data into this shape
leuco_map_mask <- bioclim_leuco_map[[1]] %>%
  mask(leuco_dist_map) * 0 + 1

plot(leuco_map_mask)


# process bioclim into shape
bc_leuco_map <- mask(
  x = bioclim_leuco_map,
  mask = leuco_map_mask
)

# # process travel data
# # this is a very large raster, we speed up the operation by cropping then masking
# # try looking at the result of only cropping to understand why both steps are needed
# travel_kenya <- crop(
#   x = travel,
#   y = kenya_mask
# ) %>%
#   mask(
#     mask = kenya_mask
#   )
# 
# # let's have a look at it



# 
# # now let's create a bias layer out of our travel data
# # it should scale from 0-1 where 0 is hard to get to and 1 is easy
# 
# plot(travel_kenya)
# 
# rescale_travel <- travel_kenya
# 
# rescale_travel[] <- rescale_travel[]/ max(rescale_travel[], na.rm = TRUE)
# 
# rescale_travel[] <- 1 - rescale_travel[]
# 
# plot(rescale_travel)
# 
# 

# now let's write our key rasters to disk
# NB the above download functions will only need to download once;
# once those objects are saved in that location, they will be re-read from disk

terra::writeRaster(
  x = leuco_map_mask,
  filename = "data/grids/leuco_map_mask.tif",
  overwrite=TRUE
)

terra::writeRaster(
  x = bc_leuco_map,
  filename = "data/grids/bc_leuco_map.tif",
  overwrite=TRUE
)

# terra::writeRaster(
#   x = rescale_travel,
#   filename = "data/grids/rescale_travel.tif"
# )

leuco_map_mask <- terra::rast("data/grids/leuco_map_mask.tif")
bc_leuco_map <- terra::rast("data/grids/bc_leuco_map.tif")

plot(bc_leuco_map)
plot(leuco_map_mask)

# Specify covariates used for the model
# BIO1 - Annual Mean Temperature
ttemp <- bc_leuco_map[[1]]
# BIO12 = Annual Precipitation
pprec <- bc_leuco_map[[12]]

covs <- c(ttemp, pprec)
names(covs) <- c("ttemp", "pprec")

writeRaster(covs,
            "data/grids/covs_hgam.tif",
            overwrite = TRUE)

pa_tabular <- leuco_data %>%
  group_by(lonx, laty) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  select(site_id, lonx, laty, Species, Number) |>
  mutate(not_complex = ifelse(Species %in% c("Leucosphyrus Subgroup", "Dirus Complex", "Leucosphyrus Group"), 0, 1)) |>
  mutate(dirus = ifelse(Species %in% c("Dirus Complex", "Anopheles Dirus", "Anopheles Cracens", "Anopheles Scanloni", "Anopheles Baimaii", "Anopheles Elegans", "Anopheles Takasagoensis", "Anopheles Nemophilous"), 1, 0)) |>
  mutate(leucosphyrus = ifelse(Species %in% c("Leucosphyrus Complex", "Anopheles Leucosphyrus", "Anopheles Latens", "Anopheles Introlatus", "Anopheles Balabacensis", "Anopheles Baisasi"), 1, 0))

max_id <- max(pa_tabular$site_id)
asia_vec_data <- read_csv("data/tabular/asia_vec_data.csv")
pa_tabular2 <- asia_vec_data %>%
  group_by(lonx, laty) %>%
  mutate(site_id = cur_group_id()+max_id, Number = 1) %>%
  ungroup() %>%
  select(site_id, lonx, laty, Species, Number) |>
  mutate(not_complex = ifelse(Species %in% c("Leucosphyrus Subgroup", "Dirus Complex", "Leucosphyrus Group"), 0, 1)) |>
  mutate(dirus = ifelse(Species %in% c("Dirus Complex", "Anopheles Dirus", "Anopheles Cracens", "Anopheles Scanloni", "Anopheles Baimaii", "Anopheles Elegans", "Anopheles Takasagoensis", "Anopheles Nemophilous"), 1, 0)) |>
  mutate(leucosphyrus = ifelse(Species %in% c("Leucosphyrus Complex", "Anopheles Leucosphyrus", "Anopheles Latens", "Anopheles Introlatus", "Anopheles Balabacensis", "Anopheles Baisasi"), 1, 0))

#could do this or not
pa_tabular <- rbind(pa_tabular,pa_tabular2)

pa_tabular <- pa_tabular %>% rename(pa = Number, sp = Species)
covs_vals <- extract(x=covs, y=pa_tabular |> dplyr::select(lonx, laty))
pa_model_data <- cbind(pa_tabular, covs_vals) |> 
  dplyr::select(-ID)
pa_model_data <- pa_model_data %>%
  distinct(site_id, sp, pa, .keep_all = TRUE)
pa_model_data$sp <- as.factor(pa_model_data$sp)

write.csv(pa_model_data,
          file = "data/tabular/hgam_leuco_data_nobias.csv",
          row.names = FALSE)

full_formula <- "pa ~ te(ttemp, pprec, bs=c('tp', 'tp')) + 
                      te(ttemp, pprec, bs=c('tp', 'tp'), by=leucosphyrus) +
                      te(ttemp, pprec, bs=c('tp', 'tp'), by=dirus) +
                      te(ttemp, pprec, sp, bs=c('tp', 'tp', 're'),
                      by=not_complex)"
full_formula <- as.formula(full_formula)
full_formula

mos_modGS <- gam(formula=full_formula, 
                 data = pa_model_data, 
                 optimizer=c("outer", "bfgs"),
                 family = "binomial", 
                 method = "REML")

# group model
covs$leucosphyrus <- 0
covs$dirus <- 0
covs$not_complex <- 0
covs$sp <- "Leucosphyrus Subgroup"

pred_pa_modGS_group <- sdm_predict(
  model = mos_modGS,
  covariates = covs
)

par(mfrow=c(1,1))
plot(pred_pa_modGS_group, main="predicted_dist - unbiased", range=c(0,1))
points(pa_model_data$lonx, pa_model_data$laty, bg=pa_model_data$pa, pch=21)

# complex model
n_cp <- 2
cp_names <- c("Leucosphyrus Complex", "Dirus Complex")
cp <- c("leucosphyrus", "dirus")
pred_pa_modGS_complex <- rast(rep(leuco_map_mask, n_cp))
covs$not_complex <- 0
par(mfrow=c(1,2))
# Plot all complexes
for(i in 1:n_cp){
  covs[["leucosphyrus"]] <- 0  
  covs[["dirus"]] <- 0  
  
  # Set current complex to 1
  covs[[cp[i]]] <- 1
  covs$sp <- cp_names[i]
  
  pred_pa_modGS_complex[[i]] <- sdm_predict(
    model = mos_modGS,
    covariates = covs
  )
  plot(pred_pa_modGS_complex[[i]], main = paste(cp_names[i], "- unbiased"), range=c(0,1))
}

# species model
covs$not_complex <- 1
species <- c("Anopheles Leucosphyrus", "Anopheles Latens", "Anopheles Introlatus", "Anopheles Balabacensis", "Anopheles Baisasi")
num_species <- length(species)
pred_pa_modGS_sp <- rast(rep(leuco_map_mask, num_species))
spcp <- data.frame(
  cp=c(1,1,1,1,1,2,2,2,2,2,2,2), sp=c("Anopheles Leucosphyrus", "Anopheles Latens",
                                      "Anopheles Introlatus", "Anopheles Balabacensis",
                                      "Anopheles Baisasi", "Anopheles Dirus", "Anopheles Cracens", "Anopheles Scanloni",
                                      "Anopheles Baimaii", "Anopheles Elegans", "Anopheles Takasagoensis", "Anopheles Nemophilous")
)
# Plot all species
par(mfrow=c(2,3))
for(i in 1:num_species){
  covs$sp <- species[i]
  comp <- spcp[i,1]
  covs[["leucosphyrus"]] <- 0  
  covs[["dirus"]] <- 0 
  # Set current complex to 1
  covs[[cp[comp]]] <- 1
  pred_pa_modGS_sp[[i]] <- sdm_predict(
    model = mos_modGS,
    covariates = covs
  )
  plot(pred_pa_modGS_sp[[i]], main = paste(species[i], "Complex", spcp[i,1], "- unbiased"), range=c(0,1))
}

par(mfrow=c(1,2))

partial_response_plot(
  model = mos_modGS,
  data = pa_model_data,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS,
  data = pa_model_data,
  var = "pprec",
  # scale = "link"
  scale = "response"
)

round(k.check(mos_modGS),2)