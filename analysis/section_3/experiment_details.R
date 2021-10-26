# Experiment details

#### Libraries and Data ----
library(dplyr)

source("./analysis/utils.R")

# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

(session_dat %>%
  filter(target_dist <= session_max_dist) %>%
    nrow(.)) / nrow(session_dat)

mean(session_dat$session_max_dist <= session_dat$target_dist)


