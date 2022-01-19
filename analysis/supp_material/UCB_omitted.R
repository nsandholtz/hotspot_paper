setwd("~/hotspot_paper")

library(tidyverse)
library(circular)
library(pracma)
library(data.table)
library(doParallel)

# Source utilities
source("./analysis/utils.R")
source("./analysis/constants.R")

# Read in argmax grids

acquisition_grid = readRDS(file = "./analysis/section_4/model_output/acquisition_grid.rds")

ARGS = commandArgs(trailingOnly = TRUE)

par_index = as.numeric(ARGS[1]) # parameter index
n_dat = as.numeric(ARGS[2]) # n_obs

# Simulation study params -------------------------------------------------

n_sim = 1000
# n_dat = 10
# par_index = 2

# Fix scale parameter at mean value from subject fits

mean_scale = readRDS(file = "./analysis/section_4/model_output/acquisition_fits.rds") %>%
  with(mean(unlist(lapply(., function(x) round(x$MAP$scale_val, digits = 2)))))
scale_val_sim = scale_vals[which.min(abs(scale_vals - mean_scale))]

min_r1 = range(r1_grid)[1]
max_r1 = range(r1_grid)[2]

get_coverage_UCB_omit = function(par_index_, scale_val_sim_, n_dat){
  
  acq_types = c("PI", "EI", "UCB")
  
  # generate data
  r1_vals = runif(n_dat, min_r1, max_r1)
  acquisition_surfaces = matrix(NA, n_dat, angle_resolution)
  for (i in 1:n_dat) {
    # Surrogate inference
    inferred_surface = infer_surface_1(
      x1 = 20,
      y1 = 0,
      r1 = r1_vals[i],
      r_sig = rew_sig,
      beta_Sig = diag(1, nrow = 2),
      init_score = 0,
      projection_grid = optim_grid
    )
    # acquisition
    acquisition_surfaces[i, ] = acquire_ucb(
      post_pred_mu = inferred_surface$post_pred_mu,
      post_pred_sig = inferred_surface$post_pred_sig,
      quant = par_vals[[3]][par_index_]
    )
  }
  which_mode = rbinom(n_dat, 1, .5) + 1
  sign_mode = c(-1, 1)[which_mode]
  acquisition_curve = cbind.data.frame(r1 = r1_vals,
                                       angle = theta_indexes[apply(acquisition_surfaces, 1, which.max)] * sign_mode)
  observed_angles = suppressWarnings(circular::rwrappedcauchy(
    n_dat,
    mu = acquisition_curve$angle,
    rho = exp(-scale_val_sim_)
  )) %>%
    ifelse(. > pi, . - 2 * pi, .)

  # Make inference ----------------------------------------------------------
  
  prob_old = prob_new = rep(.5, length(r1_vals))
  iter = 1
  while ((iter == 1 | sum(abs(prob_new - prob_old)) >= .001) & iter < 10) {
    prob_old = prob_new
    
    full_fits = list()
    for (i in 1:3) {
      if(i == 3) next
      # LOOP OVER PI, EI, UCB
      full_fits[[i]] = sym_wc_log_lik_over_grid(
        subject_rewards = r1_vals,
        subject_targets = observed_angles,
        r1_grid = r1_grid,
        acquisition_grid_ = acquisition_grid[[i]],
        par_vals_ = par_vals[[i]],
        scale_vals_ = scale_vals,
        weights_ = prob_new
      ) %>%
        mutate(acq_type = acq_types[i])
    }
    full_fits = rbindlist(full_fits) %>%
      mutate(
        likelihood = exp(log_lik),
        prior_prob = 1 / nrow(.),
        # Discrete uniform priors
        post_prob_un = likelihood * prior_prob,
        log_post_prob_un = log_lik + log(prior_prob),
      )
    
    log_norm_const = matrixStats::logSumExp(full_fits$log_post_prob_un)
    
    full_fits = full_fits %>%
      mutate(
        log_post_prob_norm = log_post_prob_un - log_norm_const,
        post_prob_norm = exp(log_post_prob_norm)
      ) %>%
      arrange(desc(log_post_prob_un)) %>%
      mutate(cum_post_prob = cumsum(post_prob_norm))
    
    MAP_iter = full_fits[which.max(full_fits$log_post_prob_un),]
    map_acq_type = ifelse(MAP_iter$acq_type == "PI",
                          1,
                          ifelse(MAP_iter$acq_type == "EI", 2, 3))
    
    map_acquisition_curve = pracma::interp2(
      r1_grid,
      par_vals[[map_acq_type]],
      Z = acquisition_grid[[map_acq_type]],
      xp = r1_vals,
      yp = rep(MAP_iter$par_val, length(r1_vals))
    )
    
    h_theta_1 = (1 / (2 * pi)) * sinh(MAP_iter$scale_val) / (cosh(MAP_iter$scale_val) - cos(observed_angles - map_acquisition_curve))
    h_theta_2 = (1 / (2 * pi)) * sinh(MAP_iter$scale_val) / (cosh(MAP_iter$scale_val) - cos(observed_angles + map_acquisition_curve))
    prob_new = h_theta_1 / (h_theta_1 + h_theta_2)
    iter = iter + 1
  }
  
  PI_in_CI = full_fits %>%
    slice(-c((min(which(cum_post_prob >= .95))+1):n())) %>%
    with(if_else("PI" %in% unique(acq_type), 1, 0))

  EI_in_CI = full_fits %>%
    slice(-c((min(which(cum_post_prob >= .95))+1):n())) %>%
    with(if_else("EI" %in% unique(acq_type), 1, 0))
  return(c(PI_in_CI, EI_in_CI))
}


# Do simulation -----------------------------------------------------------

#get_coverage_EI(15, scale_val_sim, n_dat = 100)

n.cores <- parallel::detectCores() 
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

UCB_omitted <- foreach(
  i = 1:n_sim,
  .combine = 'rbind',
  .packages = c("tidyverse", "data.table")
) %dopar% {
  get_coverage_UCB_omit(par_index, scale_val_sim, n_dat)
}

parallel::stopCluster(cl = my.cluster)

print(apply(UCB_omitted, 2, mean, na.rm = T)) 


saveRDS(UCB_omitted, file = paste0("./analysis/supp_material/simulation_output/UCB_omitted_",
                                  par_index,"_",
                                  n_dat,
                                  ".rds"))
