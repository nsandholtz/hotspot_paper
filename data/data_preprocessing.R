# Libraries ----
library(tidyverse)

# Source helper files ----
source("./analysis/utils.R")


# Read in raw data --------------------------------------------------------

# raw data parameters

x_origin = 540
y_origin = 540
ring_radius = 450

setwd("./data/raw_data")

session_col_names <- c('session_id',
                       'user_id',
                       'avg_score',
                       'final_score',
                       'max_score',
                       'time_start',
                       'total_time',
                       'num_moves')
event_col_names <- c('event_id',
                     'session_id',
                     'move',
                     'time',
                     'x_loc',
                     'y_loc',
                     'score',
                     'target_x',
                     'target_y',
                     'ring_center_x',
                     'ring_center_y',
                     'ring_radius')

event_files <- list.files(pattern = "eventrec")
session_files <- list.files(pattern = "sesrec")

session_dat <- lapply(session_files, read_csv, col_names = session_col_names) %>% 
  bind_rows() %>%
  as.data.frame() %>% 
  mutate(num_moves = num_moves - 1) %>% # Re-index moves such that the initial starting point of each round is "move 0"
  select(session_id, user_id, num_moves, max_score)

event_dat <- lapply(event_files, read_csv, col_names = event_col_names) %>% 
  bind_rows() %>%
  as.data.frame() %>%
  mutate(x_loc = x_loc - x_origin, # set origin to (0,0)
         y_loc = y_loc - y_origin,
         target_x = target_x - x_origin,
         target_y = target_y - y_origin,
         move = move - 1, # Re-index moves such that the initial starting point of each round is "move 0"
         ) %>%
  select(session_id, event_id, move, x_loc, y_loc, score, target_x, target_y) %>%
  left_join((session_dat %>%
               select(session_id,
                      user_id, 
                      max_score,
                      num_moves)), 
            by = "session_id")

session_dat <- session_dat %>%
  left_join(event_dat %>% 
              filter(move == 0) %>%
              select(session_id, target_x, target_y),
            by = "session_id")


# Creating new variables --------------------------------------------------


# * new vars for session_dat ------------------------------------------------

session_dat$target_dist <- to_polar(session_dat$target_x, 
                                            session_dat$target_y)$dist
session_dat$target_rad <- to_polar(session_dat$target_x, 
                                    session_dat$target_y)$rad

# Get reward gradients
session_dat$reward_gradient <- NA
for(sess in session_dat$session_id) {
  this_session_dat = event_dat %>%
    filter(session_id == sess)
  session_dat$reward_gradient[session_dat$session_id == sess] <- 
    (session_dat$max_score[session_dat$session_id == sess] - 
       this_session_dat$score[this_session_dat$move == 0])/
    session_dat$target_dist[session_dat$session_id == sess]
}

session_dat$session_max_dist = (session_dat$num_moves) * 20
session_dat$init_score = event_dat %>% 
  filter(move == 0) %>%
  pull(score)


# * new vars for event_dat ------------------------------------------------

# Getting each move's distance to the hotspot
event_dat$dist_hotspot <- to_polar((event_dat$x_loc - event_dat$target_x),
                         (event_dat$y_loc - event_dat$target_y))$dist

# Each move / hotspot locaction in polar coordinates

event_dat$dist <- to_polar(event_dat$x_loc,event_dat$y_loc)$dist
event_dat$rad <- to_polar(event_dat$x_loc,event_dat$y_loc)$rad

event_dat$target_dist <- to_polar(event_dat$target_x, event_dat$target_y)$dist
event_dat$target_rad <- to_polar(event_dat$target_x, event_dat$target_y)$rad

# rotate data such that move 2 falls on direction of radian = 0

event_dat$rot_rad <- NA
event_dat$target_rot_rad <- NA
for(id in unique(event_dat$user_id)){
  this_sub <- filter(event_dat, user_id == id)
  for(sess_id in unique(this_sub$session_id)){
    round_obs <- which(event_dat$user_id == id & event_dat$session_id == sess_id & event_dat$move >= 0)
    event_dat$rot_rad[round_obs[-1]] <- minuspi_to_pi(event_dat$rad[round_obs[-1]] - event_dat$rad[round_obs[-1]][1])
    event_dat$target_rot_rad[round_obs] <- minuspi_to_pi(event_dat$target_rad[round_obs] - event_dat$rad[round_obs[-1]][1])
  }
}

# calculate changes in dist / angle relative to  to previous move 

event_dat$dist_relative <- to_polar((event_dat$x_loc - lag(event_dat$x_loc)),
                                 (event_dat$y_loc - lag(event_dat$y_loc)))$dist
event_dat$dist_relative[event_dat$move == 0] <- NA

event_dat$rad_relative <- NA
for (id in unique(event_dat$user_id)) {
  this_sub <- filter(event_dat, user_id == id)
  for (sess_id in unique(this_sub$session_id)) {
    round_obs <-
      which(event_dat$user_id == id & event_dat$session_id == sess_id)
    event_dat$rad_relative[round_obs][2] <- 0
    for (i in 3:length(round_obs)) {
      event_dat$rad_relative[round_obs][i] <-
        minuspi_to_pi(to_polar(
          event_dat$x_loc[round_obs][i],
          event_dat$y_loc[round_obs][i],
          origin_x = event_dat$x_loc[round_obs][i - 1],
          origin_y = event_dat$y_loc[round_obs][i - 1]
        )[2] -
        to_polar(
          event_dat$x_loc[round_obs][i - 1],
          event_dat$y_loc[round_obs][i - 1],
          origin_x = event_dat$x_loc[round_obs][i - 2],
          origin_y = event_dat$y_loc[round_obs][i - 2]
        )[2]
        )
    }
  }
}

# change in score

event_dat$delta_score <- c(0,diff(event_dat$score))
event_dat$delta_score[c(1,which(diff(event_dat$move) < 0)+1)] <- NA


# Alternate subject labels ----

label_tool = cbind.data.frame(user_id = event_dat$user_id %>% unique(),
                              alt_id = 1:length(event_dat$user_id %>% unique()))

event_dat = event_dat %>% 
  left_join(label_tool, by = "user_id")

session_dat = session_dat %>% 
  left_join(label_tool, by = "user_id")


# Round indexes ----
session_dat = session_dat %>% 
  group_by(alt_id) %>% 
  mutate(round_index=1:n()) %>%
  select(alt_id, user_id, round_index, session_id, num_moves, 
         target_x, target_y, target_dist, target_rad, 
         reward_gradient, session_max_dist, max_score, init_score)

event_dat = event_dat %>% 
  left_join((session_dat %>%
              ungroup() %>%
               select(session_id, 
                      round_index)),
            by = "session_id") %>%
  select(alt_id, user_id, round_index, num_moves, event_id, session_id, move, 
         x_loc, y_loc, dist, rad, rot_rad, dist_relative, rad_relative,
         score, delta_score, dist_hotspot, max_score, target_x, target_y, 
         target_dist, target_rad, target_rot_rad)

plot(event_dat %>%
       filter(alt_id == 28, move == 1) %>%
       pull(delta_score),
     event_dat %>%
       filter(alt_id == 28, move == 2) %>%
       pull(rad_relative))


# Save data ---- 
setwd("~/Research/hotspot_paper")
saveRDS(session_dat,
        file = "./data/session_dat.rds")
saveRDS(event_dat,
        file = "./data/event_dat.rds")

