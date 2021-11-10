library(tidyverse)
library(circular)
library(pracma)
library(data.table)
library(doParallel)

n.cores <- parallel::detectCores() - 2

# Source utilities
source("./analysis/utils.R")
source("./analysis/constants.R")

# ACQUISITION GRID --------------------------------------------------------

acquisition_grid_threshold = readRDS(file = "./analysis/section_5/model_output/acquisition_grid_threshold.rds")
# STRUCTURE acquisition_grid_threshold[[thresh_val]][acq_family][par_value, r1_val]

# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

# Remove "burn-in" for each subject
first_n = 20
ok_sessions = session_dat %>%
  filter(round_index > first_n,
                target_dist > session_max_dist)
ok_event_dat = event_dat %>%
  filter(session_id %in% ok_sessions$session_id) 


# AUGMENTED FIT --------------------------------------------------------------

acquisition_fits_threshold = list()
acq_types = c("PI", "EI", "UCB")

for(my_sub in unique(event_dat$alt_id)){
  cat(my_sub, "\r")
  
  # Get the rewards from move 1
  rewards_1 = ok_event_dat %>%
    filter(alt_id == my_sub,
                  move == 1) %>%
    pull(delta_score) 
  
  # Calculate the targets in absolute value radians
  
  targets = ok_event_dat %>%
    mutate(lead_rad_relative = lead(rad_relative)) %>%
    filter(alt_id == my_sub,
           move == 1) %>%
    pull(lead_rad_relative)

  # Make training and test sets
  
  set.seed(my_sub)
  training = sample(1:length(rewards_1), 
                    size = round(length(rewards_1)*.8),
                    replace = F)

  full_fits_threshold = list()
  for(i in 1:3){ # LOOP OVER PI, EI, UCB
    
    my.cluster <- parallel::makeCluster(
      n.cores, 
      type = "PSOCK"
    )
    doParallel::registerDoParallel(cl = my.cluster)
    
    full_fits_threshold[[i]] <- foreach(
      j = 1:length(acquisition_grid_threshold),
      .combine = 'rbind',
      .packages = c("tidyverse", "data.table", "pracma", "circular")
    ) %dopar% {
      sym_wc_log_lik_over_grid(
        subject_rewards = rewards_1[training],
        subject_targets = targets[training],
        r1_grid = r1_grid,
        acquisition_grid_ = acquisition_grid_threshold[[j]][[i]],
        par_vals_ = par_vals[[i]],
        scale_vals_ = scale_vals) %>%
        mutate(threshold_val = threshold_vals[j],
               acq_type = acq_types[i])
    }
    parallel::stopCluster(cl = my.cluster)
  }
  full_fits_threshold = rbindlist(full_fits_threshold) %>%
      mutate(
        likelihood = exp(log_lik),
        prior_prob = 1 / nrow(.), # Discrete uniform priors
        post_prob_un = likelihood * prior_prob,
        log_post_prob_un = log_lik + log(prior_prob)
      )
  
  log_norm_const = matrixStats::logSumExp(full_fits_threshold$log_post_prob_un)
  
  full_fits_threshold = full_fits_threshold %>%
    mutate(log_post_prob_norm = log_post_prob_un - log_norm_const,
           post_prob_norm = exp(log_post_prob_norm)) 
  
  # "marginalize" over scale parameter and sort
  model_sorter = full_fits_threshold %>%
    group_by(acq_type, par_val, threshold_val) %>%
    summarise(marg_like = sum(likelihood),
              marg_post_norm = sum(post_prob_norm)) %>%
    ungroup() %>%
    arrange(desc(marg_post_norm)) %>%
    mutate(cum_post_prob = cumsum(marg_post_norm)) 

  acquisition_fits_threshold[[my_sub]] = list()
  acquisition_fits_threshold[[my_sub]]$posterior = model_sorter[model_sorter$marg_post_norm > 1e-10, ]
  acquisition_fits_threshold[[my_sub]]$MAP = full_fits_threshold[which.max(full_fits_threshold$log_post_prob_un),]
  
  # Out of sample log like
  map_acq_type = ifelse(acquisition_fits_threshold[[my_sub]]$MAP$acq_type == "PI",
                        1, ifelse(acquisition_fits_threshold[[my_sub]]$MAP$acq_type == "EI", 2, 3))
  map_threshold_index = which(threshold_vals == acquisition_fits_threshold[[my_sub]]$MAP$threshold_val)
  oos_acquisition_curve = pracma::interp2(
    r1_grid,
    par_vals[[map_acq_type]],
    Z = acquisition_grid_threshold[[map_threshold_index]][[map_acq_type]],
    xp = rewards_1[-training],
    yp = rep(acquisition_fits_threshold[[my_sub]]$MAP$par_val, length(rewards_1[-training]))
  )
  oos_log_lik = round(sym_wc_log_lik(data_vec = targets[-training],
                                     peak_vec = oos_acquisition_curve,
                                     scale = acquisition_fits_threshold[[my_sub]]$MAP$scale_val), digits = 2)
  acquisition_fits_threshold[[my_sub]]$out_log_lik = oos_log_lik 
}

saveRDS(acquisition_fits_threshold, file = "./analysis/section_5/model_output/acquisition_fits_threshold.rds")


