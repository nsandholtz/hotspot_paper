library(tidyverse)
library(circular)
library(pracma)
library(data.table)

# Source utilities
source("./analysis/utils.R")
source("./analysis/constants.R")

# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

# Read in acquisition fits 
acquisition_fits_prospect = readRDS(file = "./analysis/section_5/model_output/acquisition_fits_prospect.rds")

# Remove "burn-in" for each subject
first_n = 20
ok_sessions = session_dat %>%
  filter(round_index > first_n,
         target_dist > session_max_dist)
ok_event_dat = event_dat %>%
  filter(session_id %in% ok_sessions$session_id) 

# 
table_df = data.frame(matrix(NA,28,0))

# Subject rounds
table_df$n_rounds = ok_sessions %>%
  group_by(alt_id) %>%
  summarise(sub_rounds = n()) %>%
  pull(sub_rounds)

# MAP estimates
table_df$map_family = unlist(lapply(acquisition_fits_prospect, function(x) x$MAP$acq_type))
table_df$map_exp_param = unlist(lapply(acquisition_fits_prospect, function(x) x$MAP$par_val))
table_df$map_thresh_neg = unlist(lapply(acquisition_fits_prospect, function(x) pi/2 - x$MAP$threshold_val.x))
table_df$map_thresh_pos = unlist(lapply(acquisition_fits_prospect, function(x) pi/2 - x$MAP$threshold_val.y))
table_df$map_scale = unlist(lapply(acquisition_fits_prospect, function(x) round(x$MAP$scale_val, digits = 2)))
table_df$w_est = unlist(lapply(acquisition_fits_prospect, function(x) round(x$w_est, digits = 2)))

# Change PI at 0 to UCB at 0.5
change_inds = which(table_df$map_family == "PI" & table_df$map_exp_param == 0)
table_df$map_family[change_inds] = "UCB"
table_df$map_exp_param[change_inds] = .5


# Task Performance --------------------------------------------------------

#
table_df$mean_m2_score = ok_event_dat %>% 
  filter(move == 2) %>%
  group_by(alt_id) %>%
  summarise(mean_score = mean(delta_score)) %>% 
  pull(mean_score)

#
table_df$mean_m3_score = ok_event_dat %>% 
  filter(move == 3) %>%
  group_by(alt_id) %>%
  summarise(mean_score = mean(delta_score)) %>% 
  pull(mean_score)

table_df$mean_tot_score = ok_event_dat %>% 
  filter(move >= 2) %>%
  group_by(alt_id, session_id) %>%
  summarise(tot_score = sum(delta_score, na.rm = T)) %>% 
  group_by(alt_id) %>%
  summarise(mean_tot_score = mean(tot_score, na.rm = T)) %>%
  pull(mean_tot_score)

table_print = table_df %>%
  select(map_family, map_exp_param,
         map_thresh_neg, map_thresh_pos,
         map_scale,
         w_est,
         mean_m2_score,
         mean_m3_score,
         mean_tot_score
         )
xtable::xtable(table_print)

# Other performance summaries ---------------------------------------------

table_df$mean_exploit = ok_event_dat %>% 
  filter(move > 3) %>%
  group_by(alt_id) %>%
  summarise(mean_angle = mean(abs(rad_relative))) %>% 
  pull(mean_angle)

table_df$avg_exploit = table_df %>%
  mutate(avg_exploit = (map_thresh_neg + map_thresh_pos)/2) %>%
  pull(avg_exploit)

cor(table_df$avg_exploit, table_df$mean_m2_score)
cor(table_df$avg_exploit, table_df$mean_m3_score)
cor(table_df$mean_exploit, table_df$mean_tot_score)

