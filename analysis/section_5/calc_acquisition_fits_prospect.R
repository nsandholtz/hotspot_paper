library(tidyverse)
library(circular)
library(pracma)
library(data.table)

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

# PROSPECT FIT --------------------------------------------------------------

acquisition_fits_prospect = list()
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
  
  # condition on negative / positive delta_score on move 1
  neg_inds = which(rewards_1 < 0)
  
  neg_train_inds = neg_inds[neg_inds %in% training]
  neg_test_inds = neg_inds[!(neg_inds %in% training)]
  
  pos_inds = which(rewards_1 >= 0)
  
  pos_train_inds = pos_inds[pos_inds %in% training]
  pos_test_inds = pos_inds[!(pos_inds %in% training)]
  
  
  full_fits_threshold_neg = list()
  full_fits_threshold_pos = list()
  for(i in 1:3){ # LOOP OVER PI, EI, UCB
    for(j in 1:length(acquisition_grid_threshold)){
      aug_fit_neg = sym_wc_log_lik_over_grid(
        subject_rewards = rewards_1[neg_train_inds],
        subject_targets = targets[neg_train_inds],
        r1_grid = r1_grid,
        acquisition_grid_ = acquisition_grid_threshold[[j]][[i]],
        par_vals_ = par_vals[[i]],
        scale_vals_ = scale_vals) %>%
        mutate(threshold_val = threshold_vals[j],
               acq_type = acq_types[i])
      aug_fit_pos = sym_wc_log_lik_over_grid(
        subject_rewards = rewards_1[pos_train_inds],
        subject_targets = targets[pos_train_inds],
        r1_grid = r1_grid,
        acquisition_grid_ = acquisition_grid_threshold[[j]][[i]],
        par_vals_ = par_vals[[i]],
        scale_vals_ = scale_vals) %>%
        mutate(threshold_val = threshold_vals[j],
               acq_type = acq_types[i])
      if(j == 1){
        full_fits_threshold_neg[[i]] = aug_fit_neg
        full_fits_threshold_pos[[i]] = aug_fit_pos
      } else {
        full_fits_threshold_neg[[i]] = rbind(
          full_fits_threshold_neg[[i]],
          aug_fit_neg
        )
        full_fits_threshold_pos[[i]] = rbind(
          full_fits_threshold_pos[[i]],
          aug_fit_pos
        )
      }
    }
  }
  full_fits_threshold_neg = rbindlist(full_fits_threshold_neg) 
  full_fits_threshold_pos = rbindlist(full_fits_threshold_pos) 
  
  full_fits_prospect = full_fits_threshold_neg %>%
      left_join(full_fits_threshold_pos, by = c("par_val", 
                                                "scale_val", 
                                                "acq_type")) %>%
      mutate(log_lik = log_lik.x + log_lik.y,
             prior_prob = 1 / nrow(.), # Discrete uniform priors
             likelihood = exp(log_lik),
             post_prob_un = likelihood * prior_prob,
             log_post_prob_un = log_lik + log(prior_prob)) %>%
      select(acq_type, par_val, scale_val, threshold_val.x, 
             threshold_val.y, log_lik, prior_prob, likelihood, 
             post_prob_un, log_post_prob_un)
  
  log_norm_const = matrixStats::logSumExp(full_fits_prospect$log_post_prob_un)
  
  full_fits_prospect = full_fits_prospect %>%
    mutate(log_post_prob_norm = log_post_prob_un - log_norm_const,
           post_prob_norm = exp(log_post_prob_norm)) 
  
  # "marginalize" over scale parameter and sort
  model_sorter = full_fits_prospect %>%
    group_by(acq_type, par_val, threshold_val.x, threshold_val.y) %>%
    summarise(marg_like = sum(likelihood),
              marg_post_norm = sum(post_prob_norm)) %>%
    ungroup() %>%
    arrange(desc(marg_post_norm)) %>%
    mutate(cum_post_prob = cumsum(marg_post_norm)) 
  
  acquisition_fits_prospect[[my_sub]] = list()
  acquisition_fits_prospect[[my_sub]]$posterior = model_sorter[model_sorter$marg_post_norm > 1e-10, ]
  acquisition_fits_prospect[[my_sub]]$MAP = full_fits_prospect[which.max(full_fits_prospect$log_post_prob_un),]
  

  # Out of sample log like
  map_acq_type = ifelse(acquisition_fits_prospect[[my_sub]]$MAP$acq_type == "PI",
                        1, ifelse(acquisition_fits_prospect[[my_sub]]$MAP$acq_type == "EI", 2, 3))
  map_threshold_index_neg = which(threshold_vals == acquisition_fits_prospect[[my_sub]]$MAP$threshold_val.x)
  map_threshold_index_pos = which(threshold_vals == acquisition_fits_prospect[[my_sub]]$MAP$threshold_val.y)
  
  oos_acquisition_curve_neg = pracma::interp2(
    c(r1_grid[r1_grid < 0] ,0),
    par_vals[[map_acq_type]],
    Z = acquisition_grid_threshold[[map_threshold_index_neg]][[map_acq_type]][,r1_grid <= 0],
    xp = rewards_1[neg_test_inds],
    yp = rep(acquisition_fits_prospect[[my_sub]]$MAP$par_val, length(neg_test_inds))
  )
  oos_acquisition_curve_pos = pracma::interp2(
    r1_grid[r1_grid >= 0],
    par_vals[[map_acq_type]],
    Z = acquisition_grid_threshold[[map_threshold_index_pos]][[map_acq_type]][, r1_grid >= 0],
    xp = rewards_1[pos_test_inds],
    yp = rep(acquisition_fits_prospect[[my_sub]]$MAP$par_val, length(pos_test_inds))
  )
  
  neg_oos_log_lik = round(sym_wc_log_lik(data_vec = targets[neg_test_inds],
                                         peak_vec = oos_acquisition_curve_neg,
                                         scale = acquisition_fits_prospect[[my_sub]]$MAP$scale_val), digits = 2)  
  pos_oos_log_lik = round(sym_wc_log_lik(data_vec = targets[pos_test_inds],
                                         peak_vec = oos_acquisition_curve_pos,
                                         scale = acquisition_fits_prospect[[my_sub]]$MAP$scale_val), digits = 2)
  oos_log_lik = neg_oos_log_lik + pos_oos_log_lik
  acquisition_fits_prospect[[my_sub]]$out_log_lik = oos_log_lik 
}

saveRDS(acquisition_fits_prospect, file = "./analysis/section_5/model_output/acquisition_fits_prospect.rds")


