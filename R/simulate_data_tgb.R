# preparing data part 2

# should run without part 1 if data/grids/ contains the below rasters.
# as an alternative to running script preparing_data_1.R, grids can be
# downloaded into data/grids from:
# https://doi.org/10.26188/21630146.v1 
# you can do this manually or work out how to in R!

# load packages
library(tidyverse)
library(terra)
source("R/functions.R")

# read our rasters in
mad_mask <- terra::rast("data/grids/mad_mask.tif")
bc_mad <- terra::rast("data/grids/bc_mad.tif")
rescale_travel <- terra::rast("data/grids/rescale_travel.tif")
bias <- rescale_travel ^ 2

# add fake target species data - make some that are widespread, some that are
# more restricted (add that information as a column). We can be confident that
# all have the same sampling bias as the target species (fake species made by
# the user)

# sample only biased ones from rescale_travel

# sample biased and environmentally specific ones from rescale_travel
# and some bioclim variables

par(mfrow=c(1,1))

# make two sets of fake species, with different numbers of points
n_species_each <- 10 # 10 different species
widespread_species_n_points <- 10 + rpois(n_species_each, rlnorm(n_species_each, 5, 1))
focal_species_n_points <- 10 + rpois(n_species_each, rlnorm(n_species_each, 5, 1))

# simulate all the widespread ones just from the bias layer
widespread_all_points <- terra::spatSample(bias,
                                           sum(widespread_species_n_points),
                                           method = "weights",
                                           replace = TRUE,
                                           na.rm = TRUE,
                                           as.points = TRUE)

plot(bias) # see overall bias
points(widespread_all_points, cex = 0.2)

# label them randomly according to different species
widespread_df <- widespread_all_points %>%
  crds() %>%
  as_tibble() %>%
  mutate(
    species_id = rep(seq_len(n_species_each),
                     widespread_species_n_points),
    type = "widespread"
  )


# now do focal ones, each time selecting a different bioclim layer and including it in the sampling

focal_df <- tibble(x = numeric(0),
                   y = numeric(0),
                   species_id = integer(0),
                   type = character(0))

for(i in seq_len(n_species_each)) {
  # selecting 1 of the bioclim layers (dim(bc_mad)[3])
  which_bioclim <- sample.int(dim(bc_mad)[3], 1)
  layer <- bc_mad[[which_bioclim]] * bias
  
  layer[] <- layer[] - min(layer[], na.rm = TRUE)
  layer[] <- layer[]/ max(layer[], na.rm = TRUE)
  
  # simulate all the focal ones from different bias layers
  one_focal_species_points <- terra::spatSample(layer,
                                                focal_species_n_points[i],
                                                method = "weights",
                                                replace = TRUE,
                                                na.rm = TRUE,
                                                as.points = TRUE) 
  
  one_focal_species_df <- one_focal_species_points %>%
    crds() %>%
    as_tibble() %>%
    mutate(
      species_id = i,
      type = "focal"
    )
  
  focal_df <- rbind(focal_df, one_focal_species_df)
  
}

species_df <- rbind(focal_df, widespread_df)

write_csv(
  x = species_df,
  file = "data/tabular/species_df.csv"
)

