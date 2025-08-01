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

# crop to the southeast asia area
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

# # process travel data, right now in kenya, can change as needed
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

# write rasters to save
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

# write in the maps
leuco_map_mask <- terra::rast("data/grids/leuco_map_mask.tif")
bc_leuco_map <- terra::rast("data/grids/bc_leuco_map.tif")

# inspect them
plot(bc_leuco_map)
plot(leuco_map_mask)

# Specify covariates used for the model (only temp and prec)
# BIO1 - Annual Mean Temperature
ttemp <- bc_leuco_map[[1]]
# BIO12 = Annual Precipitation
pprec <- bc_leuco_map[[12]]

covs <- c(ttemp, pprec)
names(covs) <- c("ttemp", "pprec")

# save covariates
writeRaster(covs,
            "data/grids/covs_hgam.tif",
            overwrite = TRUE)

# prepare data and convert to same format as model data
pa_tabular <- leuco_data %>%
  group_by(lonx, laty) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  select(site_id, lonx, laty, Species, Number) |>
  mutate(not_complex = ifelse(Species %in% c("Leucosphyrus Subgroup", "Dirus Complex", "Leucosphyrus Group"), 0, 1)) |>
  mutate(dirus = ifelse(Species %in% c("Dirus Complex", "Anopheles Dirus", "Anopheles Cracens", "Anopheles Scanloni", "Anopheles Baimaii", "Anopheles Elegans", "Anopheles Takasagoensis", "Anopheles Nemophilous"), 1, 0)) |>
  mutate(leucosphyrus = ifelse(Species %in% c("Leucosphyrus Complex", "Anopheles Leucosphyrus", "Anopheles Latens", "Anopheles Introlatus", "Anopheles Balabacensis", "Anopheles Baisasi"), 1, 0))

# find max number of sites collected
max_id <- max(pa_tabular$site_id)
# read in leuco data (po points)
asia_vec_data <- read_csv("data/tabular/asia_vec_data.csv")
# convert to tabular format
pa_tabular2 <- asia_vec_data %>%
  group_by(lonx, laty) %>%
  mutate(site_id = cur_group_id()+max_id, Number = 1) %>%
  ungroup() %>%
  select(site_id, lonx, laty, Species, Number) |>
  mutate(not_complex = ifelse(Species %in% c("Leucosphyrus Subgroup", "Dirus Complex", "Leucosphyrus Group"), 0, 1)) |>
  mutate(dirus = ifelse(Species %in% c("Dirus Complex", "Anopheles Dirus", "Anopheles Cracens", "Anopheles Scanloni", "Anopheles Baimaii", "Anopheles Elegans", "Anopheles Takasagoensis", "Anopheles Nemophilous"), 1, 0)) |>
  mutate(leucosphyrus = ifelse(Species %in% c("Leucosphyrus Complex", "Anopheles Leucosphyrus", "Anopheles Latens", "Anopheles Introlatus", "Anopheles Balabacensis", "Anopheles Baisasi"), 1, 0))

#could do this or not, depending on if you want po data in
pa_tabular <- rbind(pa_tabular,pa_tabular2)

# convert the headings, etc. to be more readable
pa_tabular <- pa_tabular %>% rename(pa = Number, sp = Species)
covs_vals <- extract(x=covs, y=pa_tabular |> dplyr::select(lonx, laty))

# add to make model data for presence-absence
pa_model_data <- cbind(pa_tabular, covs_vals) |> 
  dplyr::select(-ID)
pa_model_data <- pa_model_data %>%
  distinct(site_id, sp, pa, .keep_all = TRUE)
pa_model_data$sp <- as.factor(pa_model_data$sp)

# save it
write.csv(pa_model_data,
          file = "data/tabular/hgam_leuco_data_nobias.csv",
          row.names = FALSE)
