##### Libraries and Data ----
library(dplyr)
library(readr)

# Source helper files ----
source("~/Dropbox/Luke_Research/human_acquisition/manuscript/section_1/utils.R")

# Set script parameters
my_experiment <- "experiment_2"
setwd(paste0("~/Dropbox/Luke_Research/human_acquisition/data/",my_experiment,"/Rdata"))
event_col_names <- c('event_id',
                     'session_id',
                     'trial_number',
                     'time',
                     'x_loc',
                     'y_loc',
                     'score',
                     'target_x',
                     'target_y',
                     'ring_center_x',
                     'ring_center_y',
                     'ring_radius')
session_col_names <- c('session_id',
                       'user_id',
                       'avg_score',
                       'final_score',
                       'max_score',
                       'time_start',
                       'total_time',
                       'number_tries')

event_files <- list.files(pattern = "eventrec")
session_files <- list.files(pattern = "sesrec")

#### Preprocessing ----

full_events <- lapply(event_files, read_csv, col_names = F) %>% bind_rows()
full_session <- lapply(session_files, read_csv, col_names = F) %>% bind_rows()
event_dat <- data.frame(full_events)
names(event_dat) <- event_col_names
session_dat <- data.frame(full_session)
names(session_dat) <- session_col_names

event_dat = event_dat %>% 
  left_join((session_dat %>%
               select(session_id,
                      user_id, 
                      max_score,
                      number_tries)), 
            by = "session_id")

# Getting each move's distance to the hotspot
event_dat$dist_hotspot <- to_polar((event_dat$x_loc - event_dat$target_x),
                         (event_dat$y_loc - event_dat$target_y))$dist

# Polar Coordinates -------------------------------------------------------

event_dat$dist <- to_polar((event_dat$x_loc - event_dat$ring_center_x),
                         (event_dat$y_loc - event_dat$ring_center_y))$dist
event_dat$rad <- to_polar((event_dat$x_loc - event_dat$ring_center_x),
                          (event_dat$y_loc - event_dat$ring_center_y),
                          two_pi = T)$rad
event_dat$target_dist <- to_polar((event_dat$target_x - event_dat$ring_center_x),
                                  (event_dat$target_y - event_dat$ring_center_y))$dist
event_dat$target_rad <- to_polar((event_dat$target_x - event_dat$ring_center_x),
                          (event_dat$target_y - event_dat$ring_center_y),
                          two_pi = T)$rad

# Rotated Move 2 ----------------------------------------------------------
# rotating axis for move two and all moves after that

event_dat$rot_rad <- NA
event_dat$target_rot_rad <- NA
for(id in unique(event_dat$user_id)){
  this_sub <- filter(event_dat, user_id == id)
  for(sess_id in unique(this_sub$session_id)){
    round_obs <- which(event_dat$user_id == id & event_dat$session_id == sess_id & event_dat$trial_number >= 1)
    event_dat$rot_rad[round_obs[-1]] <- event_dat$rad[round_obs[-1]] - event_dat$rad[round_obs[-1]][1]
    event_dat$target_rot_rad[round_obs] <- event_dat$target_rad[round_obs] - event_dat$rad[round_obs[-1]][1]
  }
}

# corrections for greater than 2pi?
event_dat$rot_rad <- ifelse(event_dat$rot_rad > pi, 
                          event_dat$rot_rad - 2*pi,
                          event_dat$rot_rad)
event_dat$rot_rad <- ifelse(event_dat$rot_rad < -pi,
                            event_dat$rot_rad + 2*pi,
                            event_dat$rot_rad)
event_dat$target_rot_rad <- ifelse(event_dat$target_rot_rad < -pi,
                            event_dat$target_rot_rad + 2*pi,
                            event_dat$target_rot_rad)
event_dat$rot_deg <- rad2deg(event_dat$rot_rad)




# Relative to previous move ------------------------------------------------

event_dat$dist_local <- to_polar((event_dat$x_loc - lag(event_dat$x_loc)),
                                 (event_dat$y_loc - lag(event_dat$y_loc)))$dist
event_dat$dist_local[event_dat$trial_number == 1] <- NA

event_dat$rad_local <- NA
for (id in unique(event_dat$user_id)) {
  this_sub <- filter(event_dat, user_id == id)
  for (sess_id in unique(this_sub$session_id)) {
    round_obs <-
      which(event_dat$user_id == id & event_dat$session_id == sess_id)
    event_dat$rad_local[round_obs][2] <- 0
    for (i in 3:length(round_obs)) {
      event_dat$rad_local[round_obs][i] <-
        to_polar(
          event_dat$x_loc[round_obs][i],
          event_dat$y_loc[round_obs][i],
          origin_x = event_dat$x_loc[round_obs][i - 1],
          origin_y = event_dat$y_loc[round_obs][i - 1],
          two_pi = T
        )[2] -
        to_polar(
          event_dat$x_loc[round_obs][i - 1],
          event_dat$y_loc[round_obs][i - 1],
          origin_x = event_dat$x_loc[round_obs][i - 2],
          origin_y = event_dat$y_loc[round_obs][i - 2],
          two_pi = T
        )[2]
    }
  }
}

event_dat$rad_local <- ifelse(event_dat$rad_local > pi, 
                              event_dat$rad_local - 2*pi,
                              event_dat$rad_local)
event_dat$rad_local <- ifelse(event_dat$rad_local < -pi,
                              event_dat$rad_local + 2*pi,
                              event_dat$rad_local)

# Getting incremental scores
event_dat$score_inc <- c(0,diff(event_dat$score))
event_dat$score_inc[c(1,which(diff(event_dat$trial_number) < 0)+1)] <- NA

# Get reward gradients
session_dat$reward_gradient <- NA
for(sess in session_dat$session_id) {
  this_session_dat = event_dat %>%
    filter(session_id == sess)
  session_dat$reward_gradient[session_dat$session_id == sess] <- 
    (session_dat$max_score[session_dat$session_id == sess] - 
       this_session_dat$score[this_session_dat$trial_number == 1])/
    (this_session_dat$dist_hotspot[this_session_dat$trial_number == 1])
}

#hist(session_dat$reward_gradient)

session_dat$dist_hotspot = event_dat %>% 
  filter(trial_number == 1) %>%
  pull(dist_hotspot)
session_dat$round_max_dist = (session_dat$number_tries - 1) * 20
session_dat$init_score = event_dat %>% 
  filter(trial_number == 1) %>%
  pull(score)
session_dat$target_x = event_dat %>% 
  filter(trial_number == 1) %>%
  pull(target_x)
session_dat$target_y = event_dat %>% 
  filter(trial_number == 1) %>%
  pull(target_y)

# Alternate subject labels ----

label_tool = cbind.data.frame(user_id = event_dat$user_id %>% unique(),
                              alt_id = 1:length(event_dat$user_id %>% unique()))

event_dat = event_dat %>% 
  left_join(label_tool, by = "user_id")

session_dat = session_dat %>% 
  left_join(label_tool, by = "user_id")


# Round indexes
session_dat = session_dat %>% group_by(alt_id) %>% mutate(round_index=1:n())

event_dat = event_dat %>% 
  left_join((session_dat %>%
              ungroup() %>%
               select(session_id, 
                      round_index)),
            by = "session_id")

#### Save data ---- 
saveRDS(session_dat,
        file = paste0("~/Dropbox/Luke_Research/human_acquisition/data/",
                      my_experiment,
                      "/session_dat.rds"))
saveRDS(event_dat,
        file = paste0("~/Dropbox/Luke_Research/human_acquisition/data/",
                      my_experiment,
                      "/event_dat.rds"))

