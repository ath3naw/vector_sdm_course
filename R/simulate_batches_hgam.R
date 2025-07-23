rm(list=ls())
library(tidyverse)
library(terra)
library(geodata)
library(gratia)
library(dplyr)
library(mgcv)
source("R/functions.R")

par(mfrow = c(1,1))
covs <- terra::rast("data/grids/covariates.tif")
mad_mask <- terra::rast("data/grids/mad_mask.tif")
bc_mad <- terra::rast("data/grids/bc_mad.tif")
bias <- terra::rast("data/grids/bias.tif")

n_sp <- 10
max_catch_size <- 10  # was 10

# group level model
#ttemp was 0.02, ttemp2=-0.05, ttemp test was ttemp=0, ttemp2 = -0.2/0.3, tiso=0.001
beta_group <- c(ttemp = 0.02, ttemp2 = -0.1, 
                #tiso = 0.001, 
                #tseas = -0.001, 
                #twet = 0.08, twet2 = -0.02, 
                pprec = 0.001, pprec2 = -0.000004,
                #pseas = 0.04, pseas2 = -0.0007, pwarm = 0.0002, pwarm2 = -0.000001, 
                int=0.2) # int was 0.03

cons_group <- c(ttemp = 22, 
                #twet = 24, 
                pprec = 1400)#, pseas = 80, pwarm=1000)

group_abund <- exp(beta_group["int"]+beta_group["ttemp"]*covs$ttemp+beta_group["ttemp2"]*(covs$ttemp-cons_group["ttemp"])^2+
                     #beta_group["tiso"]*covs$tiso+
                     #beta_group["tseas"]*covs$tseas+
                     #beta_group["twet"]*covs$twet+beta_group["twet2"]*(covs$twet-cons_group["twet"])^2+
                     beta_group["pprec"]*covs$pprec+beta_group["pprec2"]*(covs$pprec-cons_group["pprec"])^2)#+
#beta_group["pseas"]*covs$pseas+beta_group["pseas2"]*(covs$pseas-cons_group["pseas"])^2+
#beta_group["pwarm"]*covs$pwarm+beta_group["pwarm2"]*(covs$pwarm-cons_group["pwarm"])^2)
plot(group_abund)
rel_group_abund <- rescale_abundance(group_abund)
names(rel_group_abund) <- "relative_abundance"
plot(rel_group_abund, main="Relative Abundance")
group_prob_pres <- probability_of_presence(rel_group_abund,
                                           max_catch_size)
plot(group_prob_pres, main="Prob of Presence")
writeRaster(group_prob_pres,
            "data/grids/group_prob_pres_hglm.tif",
            overwrite = TRUE)
# number of simulations
n <- 10
cor_unbiased <- vector("list", length=n)
cor_allbiased <- vector("list", length=n)
cor_po_unbiased <- vector("list", length=n)
cor_po_allbiased <- vector("list", length=n)
all_bias <- mad_mask

for(x in 1:n){
  # species-level deviations
  species_data <- data.frame(
    species = factor(paste("sp", 1:n_sp, sep="")),
    ttemp = rnorm(n_sp, 0, 0.015), # was 0.008
    ttemp2 = rnorm(n_sp, 0, 0.07), # was 0.015
    #tiso = rnorm(n_sp, 0, 0.003),# was 0.003
    #tseas = rnorm(n_sp, 0, 0.0025),
    #twet = rnorm(n_sp, 0, 0.03),
    #twet2 = rnorm(n_sp, 0, 0.02),
    pprec = rnorm(n_sp, 0, 0.0008),
    pprec2 = rnorm(n_sp, 0, 0.0000009),
    #pseas = rnorm(n_sp, 0, 0.02),
    #pseas2 = rnorm(n_sp, 0, 0.0005),
    #pwarm = rnorm(n_sp, 0, 0.0006),
    #pwarm2 = rnorm(n_sp, 0, 0.0000008),
    int = rnorm(n_sp, 0, 0.1) # was 0.02
  )
  species_abund <- exp(beta_group["int"]+species_data$int+(beta_group["ttemp"]+species_data$ttemp)*covs$ttemp+
                         (beta_group["ttemp2"]+species_data$ttemp2)*(covs$ttemp-cons_group["ttemp"])^2+
                         #(beta_group["tiso"]+species_data$tiso)*covs$tiso+
                         #(beta_group["tseas"]+species_data$tseas)*covs$tseas+
                         #(beta_group["twet"]+species_data$twet)*covs$twet+(beta_group["twet2"]+species_data$twet2)*(covs$twet-cons_group["twet"])^2+
                         (beta_group["pprec"]+species_data$pprec)*covs$pprec+(beta_group["pprec2"]+species_data$pprec2)*(covs$pprec-cons_group["pprec"])^2)#+
  #(beta_group["pseas"]+species_data$pseas)*covs$pseas+(beta_group["pseas2"]+species_data$pseas2)*(covs$pseas-cons_group["pseas"])^2+
  #(beta_group["pwarm"]+species_data$pwarm)*covs$pwarm+(beta_group["pwarm2"]+species_data$pwarm2)*(covs$pwarm-cons_group["pwarm"])^2)

  rel_species_abund <- rast(lapply(species_abund, rescale_abundance))
  
  prob_pres <- probability_of_presence(rel_species_abund,
                                       max_catch_size)
  plot(prob_pres, main="Prob of Presence")
  writeRaster(prob_pres,
              "data/grids/spec_prob_pres_hglm.tif",
              overwrite = TRUE)
  
  covs_bounds <- data.frame(
    variable = names(covs),
    min = covs_bounds <- data.frame(
      variable = names(covs),
      min = global(covs, "min", na.rm = TRUE)[, 1],
      max = global(covs, "max", na.rm = TRUE)[, 1]
    )
  )
  
  covs_means <- global(covs, "mean", na.rm=TRUE)[, 1]
  covs_means <- as.data.frame(t(covs_means))
  names(covs_means) <- names(covs)
  
  means <- global(rel_species_abund, "mean", na.rm = TRUE)[, 1]
  i <- which.min(means)
  sp_bias <- prob_pres[[i]]
  writeRaster(sp_bias,
              "data/grids/sp_bias.tif",
              overwrite = TRUE)
  all_bias <- bias*sp_bias
  
  n_samples <- 300
  pa_tab <- generate_data_tabular(n_samples, mad_mask, prob_pres=prob_pres)
  write.csv(pa_tab,
            file = "data/tabular/hglm_pa_tab_data_med.csv",
            row.names = FALSE)
  
  # ## 5/6 complex
  # n_cp <- round(n_samples*5/6)
  # pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  # pa_model_data
  # 
  # # save these items
  # write.csv(pa_model_data,
  #           file = "data/tabular/hglm_pa_data_med_nobias_56.csv",
  #           row.names = FALSE)
  
  ## 2/3 complex
  n_cp <- round(n_samples*2/3)
  pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  
  # save these items
  write.csv(pa_model_data,
            file = "data/tabular/hglm_pa_data_med_nobias_23.csv",
            row.names = FALSE)
  
  
  # ## 1/3 complex
  # n_cp <- round(n_samples*1/3)
  # pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  # 
  # # save these items
  # write.csv(pa_model_data,
  #           file = "data/tabular/hglm_pa_data_med_nobias_13.csv",
  #           row.names = FALSE)
  # 
  # ## 1/6 complex
  # n_cp <- round(n_samples*1/6)
  # pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  # 
  # # save these items
  # write.csv(pa_model_data,
  #           file = "data/tabular/hglm_pa_data_med_nobias_16.csv",
  #           row.names = FALSE)
  
  # with bias
  pa_tab <- generate_data_tabular(n_samples, all_bias, prob_pres=prob_pres, weighted=TRUE)
  write.csv(pa_tab,
            file = "data/tabular/hglm_pa_tab_data_med_allbiased.csv",
            row.names = FALSE)
  
  plot(all_bias, main="Biased Sample Locations")
  points(pa_tab$x, pa_tab$y, pch = 16)
  
  # number of complex data points
  n_cp <- round(n_samples*2/3)
  pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  
  write.csv(pa_model_data,
            file = "data/tabular/hglm_pa_data_med_allbias_23.csv",
            row.names = FALSE)
  
  # n_cp <- round(n_samples*5/6)
  # pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  # 
  # write.csv(pa_model_data,
  #           file = "data/tabular/hglm_pa_data_med_allbias_56.csv",
  #           row.names = FALSE)
  # 
  # n_cp <- round(n_samples*1/3)
  # pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  # 
  # write.csv(pa_model_data,
  #           file = "data/tabular/hglm_pa_data_med_allbias_13.csv",
  #           row.names = FALSE)
  # 
  # n_cp <- round(n_samples*1/6)
  # pa_model_data <- generate_model_data(n_samples, n_cp, pa_tab)
  # 
  # write.csv(pa_model_data,
  #           file = "data/tabular/hglm_pa_data_med_allbias_16.csv",
  #           row.names = FALSE)
  
  ### presence-only data
  ## unbiased data
  pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_23.csv")
  occurrence_coords <- pa_model_data[pa_model_data$pa==1,]
  plot(prob_pres[[1]])
  points(occurrence_coords$x, occurrence_coords$y, pch=16)
  # presence_only
  write.csv(
    x = occurrence_coords,
    file = "data/tabular/presence_only_data_hgam.csv",
    row.names = FALSE
  )
  
  n_bg_points <- 300
  random_bg <- random_locations(mad_mask,
                                n_bg_points)
  site_id <- seq(from=max(occurrence_coords$site_id)+1, length.out=nrow(random_bg))
  plot(mad_mask)
  points(random_bg)
  po_tab <- tibble()
  for(letter in c(letters[1:n_sp], "x")){
    po_tab <- rbind(po_tab, bind_cols(
      site_id=site_id,
      as_tibble(crds(random_bg)),
      tibble(
        pa=0,
        sp=rep(letter, length.out=nrow(random_bg)),
        not_complex = ifelse(sp == "x", 0, 1)
      )
    )
    )
  }
  
  covs_vals <- extract(x=covs, y=crds(random_bg)) |> select(-any_of(c("not_complex", "sp")))
  covs_vals <- bind_rows(replicate(n_sp+1, covs_vals, simplify=FALSE))
  po_model_data <- rbind(occurrence_coords, bind_cols(po_tab, covs_vals)) |>
    mutate(sp = as.factor(sp))
  
  write.csv(
    x = po_model_data,
    file = "data/tabular/presence_only_data_rbg_hgam.csv",
    row.names = FALSE
  )
  
  pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_allbias_23.csv")
  occurrence_coords <- pa_model_data[pa_model_data$pa==1,]
  
  # presence_only
  write.csv(
    x = occurrence_coords,
    file = "data/tabular/presence_only_data_allbiased_hgam.csv",
    row.names = FALSE
  )
  
  po_model_data <- rbind(occurrence_coords, bind_cols(po_tab, covs_vals)) |>
    mutate(sp = as.factor(sp))
  
  write.csv(
    x = po_model_data,
    file = "data/tabular/presence_only_data_allbiased_rbg_hgam.csv",
    row.names = FALSE
  )
  
  
  pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_23.csv")
  pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  write.csv(
    x = pa_model_data,
    file = "data/tabular/hglm_pa_data_med_nobias_23_nocp.csv",
    row.names = FALSE
  )
  
  # pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_56.csv")
  # pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  # write.csv(
  #   x = pa_model_data,
  #   file = "data/tabular/hglm_pa_data_med_nobias_56_nocp.csv",
  #   row.names = FALSE
  # )
  
  pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_allbias_23.csv")
  pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  write.csv(
    x = pa_model_data,
    file = "data/tabular/hglm_pa_data_med_allbias_23_nocp.csv",
    row.names = FALSE
  )
  # # varying quality (fraction of complex)
  # # 5/6 complex
  # pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_56.csv")
  # pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  # write.csv(
  #   x = pa_model_data,
  #   file = "data/tabular/hglm_pa_data_med_nobias_56_nocp.csv",
  #   row.names = FALSE
  # )
  # 
  # # 1/3 complex
  # pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_13.csv")
  # pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  # write.csv(
  #   x = pa_model_data,
  #   file = "data/tabular/hglm_pa_data_med_nobias_13_nocp.csv",
  #   row.names = FALSE
  # )
  # 
  # # 1/6 complex
  # pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_16.csv")
  # pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  # write.csv(
  #   x = pa_model_data,
  #   file = "data/tabular/hglm_pa_data_med_nobias_16_nocp.csv",
  #   row.names = FALSE
  # )
  # 
  # # varying quality, biased
  # # 5/6 complex
  # pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_allbias_56.csv")
  # pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  # write.csv(
  #   x = pa_model_data,
  #   file = "data/tabular/hglm_pa_data_med_allbias_56_nocp.csv",
  #   row.names = FALSE
  # )
  # 
  # # 1/3 complex
  # pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_allbias_13.csv")
  # pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  # write.csv(
  #   x = pa_model_data,
  #   file = "data/tabular/hglm_pa_data_med_allbias_13_nocp.csv",
  #   row.names = FALSE
  # )
  # 
  # # 1/6 complex
  # pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_allbias_16.csv")
  # pa_model_data <- pa_model_data[pa_model_data$sp != "x",]
  # write.csv(
  #   x = pa_model_data,
  #   file = "data/tabular/hglm_pa_data_med_allbias_16_nocp.csv",
  #   row.names = FALSE
  # )
  
  # presence-only data with species-only data
  # unbiased
  po_model_data <- read_csv("data/tabular/presence_only_data_rbg_hgam.csv")
  po_model_data <- po_model_data[po_model_data$sp != "x",]
  write.csv(
    x = po_model_data,
    file = "data/tabular/presence_only_data_rbg_hgam_nocp.csv",
    row.names = FALSE
  )
  # biased
  po_model_data <- read_csv("data/tabular/presence_only_data_allbiased_rbg_hgam.csv")
  po_model_data <- po_model_data[po_model_data$sp != "x",]
  write.csv(
    x = po_model_data,
    file = "data/tabular/presence_only_data_allbiased_rbg_hgam_nocp.csv",
    row.names = FALSE
  )
  
  group_prob_abs <- app(1-prob_pres, fun=prod, na.rm=TRUE)
  group_prob_pres <- 1-group_prob_abs
  
  ### presence-absence data!
  pa_model_data <- read_csv("data/tabular/hglm_pa_data_med_nobias_23.csv")
  pa_model_data$sp <- as.factor(pa_model_data$sp)
  pa_model_data_nocp <- read_csv("data/tabular/hglm_pa_data_med_nobias_23_nocp.csv")
  pa_model_data_nocp$sp <- as.factor(pa_model_data_nocp$sp)
  
  ## unbiased
  mos_modGS <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                     t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                        by=not_complex), 
                   data = pa_model_data, 
                   family = "binomial", 
                   method = "REML")
  mos_modGS_nocp <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                          t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                             by=not_complex), 
                        data = pa_model_data_nocp, 
                        family = "binomial", 
                        method = "REML")
  mos_modGS_nogp <- gam(pa ~ t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                               by=not_complex), 
                          data = pa_model_data_nocp, 
                          family = "binomial", 
                          method = "REML")
  
  pred_pa_modGS <- rast(rep(mad_mask, n_sp))
  pred_pa_modGS_nocp <- rast(rep(mad_mask, n_sp))
  pred_pa_modGS_nogp <- rast(rep(mad_mask, n_sp))
  covs$not_complex <- 1
  cor_unbiased[[x]] <- data.frame(CP = 0, NoCP=0, NoGP=0)
  for(letter in letters[1:n_sp]){
    covs$sp <- letter
    i <- match(letter, letters)
    pred_pa_modGS[[i]] <- sdm_predict(
      model = mos_modGS,
      covariates = covs
    )
    pred_pa_modGS_nocp[[i]] <- sdm_predict(
      model = mos_modGS_nocp,
      covariates = covs
    )
    pred_pa_modGS_nogp[[i]] <- sdm_predict(
      model = mos_modGS_nogp,
      covariates = covs
    )
    par(mfrow=c(2,2))
    cor_unbiased[[x]][i, 1] <- compute_cor(prob_pres[[i]], pred_pa_modGS[[i]])
    cor_unbiased[[x]][i, 2] <- compute_cor(prob_pres[[i]], pred_pa_modGS_nocp[[i]])
    cor_unbiased[[x]][i, 3] <- compute_cor(prob_pres[[i]], pred_pa_modGS_nogp[[i]])
  }
  
  # biased
  pa_model_data_allbiased <- read_csv("data/tabular/hglm_pa_data_med_allbias_23.csv")
  pa_model_data_allbiased$sp <- as.factor(pa_model_data_allbiased$sp)
  pa_model_data_allbiased_nocp <- read_csv("data/tabular/hglm_pa_data_med_allbias_23_nocp.csv")
  pa_model_data_allbiased_nocp$sp <- as.factor(pa_model_data_allbiased_nocp$sp)
  
  mos_modGS_allbiased <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                               t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                  by=not_complex), 
                             data = pa_model_data_allbiased, 
                             family = "binomial", 
                             method = "REML")
  mos_modGS_allbiased_nocp <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                                    t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                       by=not_complex), 
                                  data = pa_model_data_allbiased_nocp, 
                                  family = "binomial", 
                                  method = "REML")
  mos_modGS_allbiased_nogp <- gam(pa ~ t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                       by=not_complex), 
                                  data = pa_model_data_allbiased_nocp, 
                                  family = "binomial", 
                                  method = "REML")
  
  pred_pa_modGS_allbiased <- rast(rep(mad_mask, n_sp))
  pred_pa_modGS_allbiased_nocp <- rast(rep(mad_mask, n_sp))
  pred_pa_modGS_allbiased_nogp <- rast(rep(mad_mask, n_sp))
  covs$not_complex <- 1
  cor_allbiased[[x]] <- data.frame(CP = 0, NoCP=0, NoGP=0)
  for(letter in letters[1:n_sp]){
    covs$sp <- letter
    i <- match(letter, letters)
    pred_pa_modGS_allbiased[[i]] <- sdm_predict(
      model = mos_modGS_allbiased,
      covariates = covs
    )
    pred_pa_modGS_allbiased_nocp[[i]] <- sdm_predict(
      model = mos_modGS_allbiased_nocp,
      covariates = covs
    )
    pred_pa_modGS_allbiased_nogp[[i]] <- sdm_predict(
      model = mos_modGS_allbiased_nogp,
      covariates = covs
    )
    par(mfrow=c(2,2))
    cor_allbiased[[x]][i, 1] <- compute_cor(prob_pres[[i]], pred_pa_modGS_allbiased[[i]])
    cor_allbiased[[x]][i, 2] <- compute_cor(prob_pres[[i]], pred_pa_modGS_allbiased_nocp[[i]])
    cor_allbiased[[x]][i, 3] <- compute_cor(prob_pres[[i]], pred_pa_modGS_allbiased_nogp[[i]])
  }
  
  ### presence-only data!
  po_model_data <- read_csv("data/tabular/presence_only_data_rbg_hgam.csv")
  po_model_data$sp <- as.factor(po_model_data$sp)
  po_model_data_nocp <- read_csv("data/tabular/presence_only_data_rbg_hgam_nocp.csv")
  po_model_data_nocp$sp <- as.factor(po_model_data_nocp$sp)
  
  ## unbiased
  mos_modGS_rbg <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                         t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                            by=not_complex), 
                       data = po_model_data, 
                       family = "binomial", 
                       method = "REML")
  mos_modGS_rbg_nocp <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                              t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                 by=not_complex), 
                            data = po_model_data_nocp, 
                            family = "binomial", 
                            method = "REML")
  mos_modGS_rbg_nogp <- gam(pa ~ t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                 by=not_complex), 
                            data = po_model_data_nocp, 
                            family = "binomial", 
                            method = "REML")
  
  pred_po_modGS <- rast(rep(mad_mask, n_sp))
  pred_po_modGS_nocp <- rast(rep(mad_mask, n_sp))
  pred_po_modGS_nogp <- rast(rep(mad_mask, n_sp))
  covs$not_complex <- 1
  cor_po_unbiased[[x]] <- data.frame(CP = 0, NoCP=0, NoGP=0)
  for(letter in letters[1:n_sp]){
    covs$sp <- letter
    i <- match(letter, letters)
    pred_po_modGS[[i]] <- sdm_predict(
      model = mos_modGS_rbg,
      covariates = covs
    )
    pred_po_modGS_nocp[[i]] <- sdm_predict(
      model = mos_modGS_rbg_nocp,
      covariates = covs
    )
    pred_po_modGS_nogp[[i]] <- sdm_predict(
      model = mos_modGS_rbg_nogp,
      covariates = covs
    )
    par(mfrow=c(2,2))
    cor_po_unbiased[[x]][i, 1] <- compute_cor_po(prob_pres[[i]], pred_po_modGS[[i]])
    cor_po_unbiased[[x]][i, 2] <- compute_cor_po(prob_pres[[i]], pred_po_modGS_nocp[[i]])
    cor_po_unbiased[[x]][i, 3] <- compute_cor_po(prob_pres[[i]], pred_po_modGS_nogp[[i]])
  }
  # biased po data
  po_model_data_allbiased <- read_csv("data/tabular/presence_only_data_allbiased_rbg_hgam.csv")
  po_model_data_allbiased$sp <- as.factor(po_model_data_allbiased$sp)
  po_model_data_allbiased_nocp <- read_csv("data/tabular/presence_only_data_allbiased_rbg_hgam_nocp.csv")
  po_model_data_allbiased_nocp$sp <- as.factor(po_model_data_allbiased_nocp$sp)
  
  mos_modGS_allbiased_rbg <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                                   t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                      by=not_complex), 
                                 data = po_model_data_allbiased, 
                                 family = "binomial", 
                                 method = "REML")
  mos_modGS_allbiased_rbg_nocp <- gam(pa ~ te(ttemp, pprec, bs=c("tp", "tp")) + 
                                        t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                           by=not_complex), 
                                      data = po_model_data_allbiased_nocp, 
                                      family = "binomial", 
                                      method = "REML")
  mos_modGS_allbiased_rbg_nogp <- gam(pa ~ t2(ttemp, pprec, sp, bs=c("tp", "tp", "re"),
                                           by=not_complex), 
                                      data = po_model_data_allbiased_nocp, 
                                      family = "binomial", 
                                      method = "REML")
  
  pred_po_modGS_allbiased <- rast(rep(mad_mask, n_sp))
  pred_po_modGS_allbiased_nocp <- rast(rep(mad_mask, n_sp))
  pred_po_modGS_allbiased_nogp <- rast(rep(mad_mask, n_sp))
  covs$not_complex <- 1
  cor_po_allbiased[[x]] <- data.frame(CP = 0, NoCP=0, NoGP=0)
  for(letter in letters[1:n_sp]){
    covs$sp <- letter
    i <- match(letter, letters)
    pred_po_modGS_allbiased[[i]] <- sdm_predict(
      model = mos_modGS_allbiased_rbg,
      covariates = covs
    )
    pred_po_modGS_allbiased_nocp[[i]] <- sdm_predict(
      model = mos_modGS_allbiased_rbg_nocp,
      covariates = covs
    )
    pred_po_modGS_allbiased_nogp[[i]] <- sdm_predict(
      model = mos_modGS_allbiased_rbg_nogp,
      covariates = covs
    )
    par(mfrow=c(2,2))
    cor_po_allbiased[[x]][i, 1] <- compute_cor_po(prob_pres[[i]], pred_po_modGS_allbiased[[i]])
    cor_po_allbiased[[x]][i, 2] <- compute_cor_po(prob_pres[[i]], pred_po_modGS_allbiased_nocp[[i]])
    cor_po_allbiased[[x]][i, 3] <- compute_cor_po(prob_pres[[i]], pred_po_modGS_allbiased_nogp[[i]])
  }
}
cor_unbiased <- bind_rows(cor_unbiased, .id="column_label") %>% select(-column_label)
cor_allbiased <- bind_rows(cor_allbiased, .id="column_label") %>% select(-column_label)
cor_po_unbiased <- bind_rows(cor_po_unbiased, .id="column_label") %>% select(-column_label)
cor_po_allbiased <- bind_rows(cor_po_allbiased, .id="column_label") %>% select(-column_label)

write.csv(cor_unbiased,
          file = "data/tabular/cor_unbiased_23.csv",
          row.names = FALSE)
write.csv(cor_allbiased,
          file = "data/tabular/cor_allbiased_23.csv",
          row.names = FALSE)
write.csv(cor_po_unbiased,
          file = "data/tabular/cor_po_unbiased_23.csv",
          row.names = FALSE)
write.csv(cor_po_allbiased,
          file = "data/tabular/cor_po_allbiased_23.csv",
          row.names = FALSE)
cor_unbiased_avg <- cor_unbiased %>%
  colMeans
cor_unbiased_avg
cor_allbiased_avg <- cor_allbiased %>%
  colMeans
cor_allbiased_avg
cor_po_unbiased_avg <- cor_po_unbiased %>%
  colMeans
cor_po_unbiased_avg
cor_po_allbiased_avg <- cor_po_allbiased %>%
  colMeans
cor_po_allbiased_avg

# calculate overall how many of each model was the winner
max_col_unbiased <- max.col(cor_unbiased)
max_unbiased <- table(max_col_unbiased)
names(max_unbiased) <- c("CP", "NoCP", "NoGP")
max_unbiased

max_col_allbiased <- max.col(cor_allbiased)
max_allbiased <- table(max_col_allbiased)
names(max_allbiased) <- c("CP", "NoCP", "NoGP")
max_allbiased

max_col_po_unbiased <- max.col(cor_po_unbiased)
max_po_unbiased <- table(max_col_po_unbiased)
names(max_po_unbiased) <- c("CP", "NoCP", "NoGP")
max_po_unbiased

max_col_po_allbiased <- max.col(cor_po_allbiased)
max_po_allbiased <- table(max_col_po_allbiased)
names(max_po_allbiased) <- c("CP", "NoCP", "NoGP")
max_po_allbiased

par(mfrow=c(1,1))
plot(cor_unbiased$CP, pch=16)
points(cor_unbiased$NoCP, col="red", pch=16)
points(cor_unbiased$NoGP, col="blue", pch=16)
legend("bottomleft", legend=c("CP", "NoCP", "NoGP"),
       col=c("black", "red", "blue"), pch=16, cex=0.8)

plot(cor_allbiased$CP, pch=16)
points(cor_allbiased$NoCP, col="red", pch=16)
points(cor_allbiased$NoGP, col="blue", pch=16)
legend("bottomleft", legend=c("CP", "NoCP", "NoGP"),
       col=c("black", "red", "blue"), pch=16, cex=0.8)

plot(cor_po_unbiased$CP, pch=16)
points(cor_po_unbiased$NoCP, col="red", pch=16)
points(cor_po_unbiased$NoGP, col="blue", pch=16)
legend("bottomleft", legend=c("CP", "NoCP", "NoGP"),
       col=c("black", "red", "blue"), pch=16, cex=0.8)

plot(cor_po_allbiased$CP, pch=16)
points(cor_po_allbiased$NoCP, col="red", pch=16)
points(cor_po_allbiased$NoGP, col="blue", pch=16)
legend("bottomleft", legend=c("CP", "NoCP", "NoGP"),
       col=c("black", "red", "blue"), pch=16, cex=0.8)
