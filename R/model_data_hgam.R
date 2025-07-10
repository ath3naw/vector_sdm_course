rm(list = ls())
library(tidyverse)
library(terra)
library(geodata)
library(gratia)
library(dplyr)
library(mgcv)
source("R/functions.R")

prob_pres <- terra::rast("data/grids/spec_prob_pres_hglm.tif")
covs <- terra::rast("data/grids/covariates.tif")
mad_mask <- terra::rast("data/grids/mad_mask.tif")
bias <- terra::rast("data/grids/bias.tif")
spbias <- terra::rast("data/grids/sp_bias.tif") # can comment out whichever one

pa_tab <- read_csv("data/tabular/hglm_pa_tab_data_med.csv")
pa_tab_biased <- read_csv("data/tabular/hglm_pa_tab_data_med_biased.csv")

# calculating group probability of presence from individual species
group_prob_abs <- app(1-prob_pres, fun=prod, na.rm=TRUE)
group_prob_pres <- 1-group_prob_abs
n_sp <- nlyr(prob_pres)

pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_23.csv")
pa_model_data$sp <- as.factor(pa_model_data$sp)
pa_model_data_biased <- read_csv("data/tabular/hglm_pa_data_med_bias_23.csv")
pa_model_data_biased$sp <- as.factor(pa_model_data_biased$sp)
pa_model_data_spbiased <- read_csv("data/tabular/hglm_pa_data_med_spbias_23.csv")
pa_model_data_spbiased$sp <- as.factor(pa_model_data_spbiased$sp)
po_model_data <- read_csv("data/tabular/presence_only_data_hgam.csv")
po_model_data$sp <- as.factor(po_model_data$sp)


# unbiased
mos_modGS <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                   t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                      by=not_complex), 
                 data = pa_model_data, 
                 family = "binomial", 
                 method = "REML")

covs$not_complex <- 0
covs$sp <- "x"
# specify the species
pred_pa_modGS_group <- sdm_predict(
  model = mos_modGS,
  covariates = covs
)
par(mfrow=c(1,1))
plot(pred_pa_modGS_group)
par(mfrow=c(1,2))
plot(group_prob_pres, main="Group Prob of Pres", range=c(0,1))
plot(pred_pa_modGS_group, main="predicted_dist - unbiased", range=c(0,1))

pred_pa_modGS <- rast(rep(mad_mask, n_sp))
covs$not_complex <- 1
par(mfrow=c(4,2))
# Plot all species
for(letter in letters[1:n_sp]){
  covs$sp <- letter
  i <- match(letter, letters)
  pred_pa_modGS[[i]] <- sdm_predict(
    model = mos_modGS,
    covariates = covs
  )
  plot(prob_pres[[i]], main = "True Prob of Pres", range=c(0,1))
  plot(pred_pa_modGS[[i]], main = paste("Species", letter, "- unbiased"), range=c(0,1))
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


# biased
mos_modGS_biased <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                          t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                             by=not_complex), 
                        data = pa_model_data_biased, 
                        family = "binomial", 
                        method = "REML")

mos_modGS_spbiased <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                          t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                             by=not_complex), 
                        data = pa_model_data_spbiased, 
                        family = "binomial", 
                        method = "REML")

covs$not_complex <- 0
covs$sp <- "x"
# specify the species
pred_pa_modGS_group_biased <- sdm_predict(
  model = mos_modGS_biased,
  covariates = covs
)
pred_pa_modGS_group_spbiased <- sdm_predict(
  model = mos_modGS_spbiased,
  covariates = covs
)

par(mfrow=c(2,2))
plot(group_prob_pres, main="Group Prob of Pres", range=c(0,1))
plot(pred_pa_modGS_group, main="predicted_dist - unbiased", range=c(0,1))
plot(pred_pa_modGS_group_biased, main="predicted_dist - travel biased", range=c(0,1))
plot(pred_pa_modGS_group_spbiased, main="predicted_dist - species biased", range=c(0,1))

pred_pa_modGS_biased <- rast(rep(mad_mask, n_sp))
pred_pa_modGS_spbiased <- rast(rep(mad_mask, n_sp))
covs$not_complex <- 1
# Plot all species
for(letter in letters[1:n_sp]){
  covs$sp <- letter
  i <- match(letter, letters)
  pred_pa_modGS_biased[[i]] <- sdm_predict(
    model = mos_modGS_biased,
    covariates = covs
  )
  pred_pa_modGS_spbiased[[i]] <- sdm_predict(
    model = mos_modGS_spbiased,
    covariates = covs
  )
  par(mfrow=c(2,2))
  plot(prob_pres[[i]], main = "True Prob of Pres", range=c(0,1))
  plot(pred_pa_modGS[[i]], main = paste("Species", letter, "- unbiased"), range=c(0,1))
  plot(pred_pa_modGS_biased[[i]], main = paste("Species", letter, "- travel biased"), range=c(0,1))
  #points(pa_model_data_biased$x, pa_model_data_biased$y, pch=21, bg=pa_tab_biased[[letter]])
  plot(pred_pa_modGS_spbiased[[i]], main = paste("Species", letter, "- species biased"), range=c(0,1))
}

par(mfrow=c(1,1))

plot(bias)
points(pa_model_data_biased$x, pa_model_data_biased$y, pch=20)

par(mfrow=c(2,2))
partial_response_plot(
  model = mos_modGS_biased,
  data = pa_model_data_biased,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_biased,
  data = pa_model_data_biased,
  var = "pprec",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_spbiased,
  data = pa_model_data_spbiased,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_spbiased,
  data = pa_model_data_spbiased,
  var = "pprec",
  # scale = "link"
  scale = "response"
)

round(k.check(mos_modGS_biased),2)

# presence-only model, unbiased
po_model_data_rbg <- read_csv("data/tabular/presence_only_data_rbg_hgam.csv")
po_model_data_rbg$sp <- as.factor(po_model_data_rbg$sp)
mos_modGS_rbg <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                   t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                      by=not_complex), 
                 data = po_model_data_rbg, 
                 family = "binomial", 
                 method = "REML")

covs$not_complex <- 0
covs$sp <- "x"
# specify the species
pred_po_modGS_rbg_group <- sdm_predict(
  model = mos_modGS_rbg,
  covariates = covs
)
par(mfrow=c(1,1))
plot(pred_po_modGS_rbg_group)
par(mfrow=c(1,2))
plot(group_prob_pres, main="Group Prob of Pres", range=c(0,1))
plot(pred_po_modGS_rbg_group, main="predicted_dist - unbiased")

pred_po_modGS_rbg <- rast(rep(mad_mask, n_sp))
covs$not_complex <- 1
# Plot all species
for(letter in letters[1:n_sp]){
  covs$sp <- letter
  i <- match(letter, letters)
  pred_po_modGS_rbg[[i]] <- sdm_predict(
    model = mos_modGS_rbg,
    covariates = covs
  )
  par(mfrow=c(1,2))
  plot(prob_pres[[i]], main = "True Prob of Pres")
  plot(pred_po_modGS_rbg[[i]], main = paste("Species", letter, "- unbiased PO"))
}

par(mfrow=c(1,2))

partial_response_plot(
  model = mos_modGS_rbg,
  data = po_model_data_rbg,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_rbg,
  data = po_model_data_rbg,
  var = "pprec",
  # scale = "link"
  scale = "response"
)

round(k.check(mos_modGS_rbg),2)

## presence-only data, bias in sampling
# biased data with random background
po_model_data_biased_rbg <- read_csv("data/tabular/presence_only_data_biased_rbg_hgam.csv")
po_model_data_biased_rbg$sp <- as.factor(po_model_data_biased_rbg$sp)
po_model_data_spbiased_rbg <- read_csv("data/tabular/presence_only_data_spbiased_rbg_hgam.csv")
po_model_data_spbiased_rbg$sp <- as.factor(po_model_data_spbiased_rbg$sp)

mos_modGS_biased_rbg <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                       t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                          by=not_complex), 
                     data = po_model_data_biased_rbg, 
                     family = "binomial", 
                     method = "REML")
mos_modGS_spbiased_rbg <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                              t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                 by=not_complex), 
                            data = po_model_data_spbiased_rbg, 
                            family = "binomial", 
                            method = "REML")

covs$not_complex <- 0
covs$sp <- "x"
# specify the species
pred_po_modGS_biased_rbg_group <- sdm_predict(
  model = mos_modGS_biased_rbg,
  covariates = covs
)
pred_po_modGS_spbiased_rbg_group <- sdm_predict(
  model = mos_modGS_spbiased_rbg,
  covariates = covs
)

par(mfrow=c(2,2))
plot(group_prob_pres, main="Group Prob of Pres")
plot(pred_po_modGS_rbg_group, main="predicted_dist - unbiased")
plot(pred_po_modGS_biased_rbg_group, main="predicted_dist - travel biased")
plot(pred_po_modGS_spbiased_rbg_group, main="predicted_dist - species biased")

pred_po_modGS_biased_rbg <- rast(rep(mad_mask, n_sp))
pred_po_modGS_spbiased_rbg <- rast(rep(mad_mask, n_sp))
covs$not_complex <- 1
# Plot all species
for(letter in letters[1:n_sp]){
  covs$sp <- letter
  i <- match(letter, letters)
  pred_po_modGS_biased_rbg[[i]] <- sdm_predict(
    model = mos_modGS_biased_rbg,
    covariates = covs
  )
  pred_po_modGS_spbiased_rbg[[i]] <- sdm_predict(
    model = mos_modGS_spbiased_rbg,
    covariates = covs
  )
  par(mfrow=c(2,2))
  plot(prob_pres[[i]], main = "True Prob of Pres")
  plot(pred_po_modGS_rbg[[i]], main = paste("Species", letter, "- unbiased PO"))
  plot(pred_po_modGS_biased_rbg[[i]], main = paste("Species", letter, "- travel biased PO"))
  plot(pred_po_modGS_spbiased_rbg[[i]], main = paste("Species", letter, "- species biased PO"))
}

par(mfrow=c(3,2))
partial_response_plot(
  model = mos_modGS_rbg,
  data = po_model_data_rbg,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_rbg,
  data = po_model_data_rbg,
  var = "pprec",
  # scale = "link"
  scale = "response"
)

partial_response_plot(
  model = mos_modGS_biased_rbg,
  data = po_model_data_biased_rbg,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_biased_rbg,
  data = po_model_data_biased_rbg,
  var = "pprec",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_spbiased_rbg,
  data = po_model_data_spbiased_rbg,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_spbiased_rbg,
  data = po_model_data_spbiased_rbg,
  var = "pprec",
  # scale = "link"
  scale = "response"
)
round(k.check(mos_modGS_rbg),2)
round(k.check(mos_modGS_biased_rbg),2)
round(k.check(mos_modGS_spbiased_rbg),2)

## biased background with biased data
po_model_data_biased_bbg <- read_csv("data/tabular/presence_only_data_biased_bbg_hgam.csv")
po_model_data_biased_bbg$sp <- as.factor(po_model_data_biased_bbg$sp)
po_model_data_spbiased_bbg <- read_csv("data/tabular/presence_only_data_spbiased_bbg_hgam.csv")
po_model_data_spbiased_bbg$sp <- as.factor(po_model_data_spbiased_bbg$sp)

mos_modGS_biased_bbg <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                              t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                 by=not_complex), 
                            data = po_model_data_biased_bbg, 
                            family = "binomial", 
                            method = "REML")
mos_modGS_spbiased_bbg <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                                t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                   by=not_complex), 
                              data = po_model_data_spbiased_bbg, 
                              family = "binomial", 
                              method = "REML")

covs$not_complex <- 0
covs$sp <- "x"
# specify the species
pred_po_modGS_biased_bbg_group <- sdm_predict(
  model = mos_modGS_biased_bbg,
  covariates = covs
)
pred_po_modGS_spbiased_bbg_group <- sdm_predict(
  model = mos_modGS_spbiased_bbg,
  covariates = covs
)

par(mfrow=c(2,2))
plot(group_prob_pres, main="Group Prob of Pres")
plot(pred_po_modGS_rbg_group, main = "predicted_dist - travel biased")
plot(pred_po_modGS_biased_bbg_group, main="predicted_dist - travel biased")
plot(pred_po_modGS_spbiased_bbg_group, main="predicted_dist - species biased")

pred_po_modGS_biased_bbg <- rast(rep(mad_mask, n_sp))
pred_po_modGS_spbiased_bbg <- rast(rep(mad_mask, n_sp))
covs$not_complex <- 1
# Plot all species
for(letter in letters[1:n_sp]){
  covs$sp <- letter
  i <- match(letter, letters)
  pred_po_modGS_biased_bbg[[i]] <- sdm_predict(
    model = mos_modGS_biased_bbg,
    covariates = covs
  )
  pred_po_modGS_spbiased_bbg[[i]] <- sdm_predict(
    model = mos_modGS_spbiased_bbg,
    covariates = covs
  )
  par(mfrow=c(2,2))
  plot(prob_pres[[i]], main = "True Prob of Pres")
  plot(pred_po_modGS_rbg[[i]], main = paste("Species", letter, "- unbiased PO"))
  plot(pred_po_modGS_biased_bbg[[i]], main = paste("Species", letter, "- travel biased PO"))
  plot(pred_po_modGS_spbiased_bbg[[i]], main = paste("Species", letter, "- species biased PO"))
}

par(mfrow=c(3,2))
partial_response_plot(
  model = mos_modGS_rbg,
  data = po_model_data_rbg,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_rbg,
  data = po_model_data_rbg,
  var = "pprec",
  # scale = "link"
  scale = "response"
)

partial_response_plot(
  model = mos_modGS_biased_bbg,
  data = po_model_data_biased_bbg,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_biased_bbg,
  data = po_model_data_biased_bbg,
  var = "pprec",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_spbiased_bbg,
  data = po_model_data_spbiased_bbg,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_spbiased_bbg,
  data = po_model_data_spbiased_bbg,
  var = "pprec",
  # scale = "link"
  scale = "response"
)
round(k.check(mos_modGS_rbg),2)
round(k.check(mos_modGS_biased_bbg),2)
round(k.check(mos_modGS_spbiased_bbg),2)

## bias in complex data
pa_model_data_cbias <- read_csv("data/tabular/hglm_pa_data_max_cbias_23.csv")
pa_model_data_cbias$sp <- as.factor(pa_model_data_cbias$sp)
mos_modGS_cbias <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                              t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                 by=not_complex), 
                            data = pa_model_data_cbias, 
                            family = "binomial", 
                            method = "REML")

covs$not_complex <- 0
covs$sp <- "x"
# specify the species
pred_po_modGS_cbias_group <- sdm_predict(
  model = mos_modGS_cbias,
  covariates = covs
)
par(mfrow=c(1,1))
plot(pred_po_modGS_cbias_group)
par(mfrow=c(1,2))
plot(group_prob_pres, main="Group Prob of Pres")
plot(pred_po_modGS_cbias_group, main="predicted_dist - biased complex")

pred_po_modGS_cbias <- rast(rep(mad_mask, n_sp))
covs$not_complex <- 1
# Plot all species
for(letter in letters[1:n_sp]){
  covs$sp <- letter
  i <- match(letter, letters)
  pred_po_modGS_cbias[[i]] <- sdm_predict(
    model = mos_modGS_cbias,
    covariates = covs
  )
  par(mfrow=c(1,2))
  plot(prob_pres[[i]], main = "True Prob of Pres")
  plot(pred_po_modGS_cbias[[i]], main = paste("Species", letter, "- biased complex"))
}

par(mfrow=c(1,2))

partial_response_plot(
  model = mos_modGS_cbias,
  data = pa_model_data_cbias,
  var = "ttemp",
  # scale = "link"
  scale = "response"
)
partial_response_plot(
  model = mos_modGS_cbias,
  data = pa_model_data_cbias,
  var = "pprec",
  # scale = "link"
  scale = "response"
)

round(k.check(mos_modGS_cbias),2)

# # GI model
# pa_practice <- pa_model_data |>
#   mutate(sp_indiv= ifelse(not_complex == 1, as.character(sp), NA)) |>
#   mutate(sp_indiv = factor(sp_indiv))
# mos_modGI <- gam(pa ~ s(ttemp, bs="tp") +
#                    s(ttemp, by=sp_indiv, bs="tp") +
#                    s(sp, by= not_complex, bs="re"),
#                  data=pa_practice,
#                  family = "binomial",
#                  method="REML")
# 
# 
# # Use one consistent level for each factor
# covs$sp_indiv <- factor(levels(pa_practice$sp_indiv)[1],
#                         levels = levels(pa_practice$sp_indiv))
# 
# covs$sp <- factor(levels(pa_practice$sp)[1],
#                   levels = levels(pa_practice$sp))
# 
# covs$not_complex <- 0  # or 0 if more appropriate
# 
# # Custom predict function to extract only the s(ttemp) term
# predict_s_ttemp <- function(model, newdata) {
#   predict(model, newdata = newdata, type = "terms", terms = "s(ttemp)")
# }
# 
# # Apply prediction function to raster
# pred_pa_modGI <- predict(covs, mos_modGI, fun = predict_s_ttemp)
# 
# plot(ttemp_effect_raster, main = "Effect of ttemp (s(ttemp))")
# par(mfrow=c(1,2))
# plot(group_prob_pres, main="True prob of pres")
# plot(pred_pa_modGI, main="predicted_distribution")
# 
# covs$not_complex <- 1
# for(letter in letters[1:n_sp]){
#   covs$sp <- letter
#   pred_pa_modGI <- sdm_predict(
#     model = mos_modGI,
#     covariates = covs
#   )
#   par(mfrow=c(1,2))
#   plot(prob_pres[[match(letter, letters)]], main = "True Distribution")
#   plot(pred_pa_modGI, main = paste("Species ", letter))
# }
# 
# 
# 
# # Custom function to only return global smooth
# predict_global_smooth <- function(model, data, ...) {
#   # Prepare data for prediction
#   data$sp <- factor("x", levels = levels(pa_practice$sp))
#   data$not_complex <- 0
#   data$sp_indiv <- NA
#   data$sp_indiv <- factor(data$sp_indiv, levels = levels(pa_practice$sp_indiv))
#   
#   # Predict with excluded smooth
#   predict(model, newdata = data, type = "response")
# }
# 
# pred_raster_global <- terra::predict(
#   object = covs,
#   model = mos_modGI,
#   fun = predict_global_smooth,
#   type = "response"
# )
# 
# pred_pa_modGI <- sdm_predict(
#   model = mos_modGI,
#   covariates = covs
# )
# plot(pred_pa_modGI)
# # mos_modGI <- gam(pa ~ sp*not_complex +
# #                    te(ttemp, pprec, bs=c("tp","tp")) +
# #                    te(ttemp, pprec, by=sp, bs=c("tp","tp")),
# #                  data=pa_model_data,
# #                  family="binomial",
# #                  method="REML")
# 
