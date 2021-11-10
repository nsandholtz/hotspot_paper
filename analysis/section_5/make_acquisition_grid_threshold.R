library(tidyverse)

# Source utilities and constants
source("./analysis/utils.R")
source("./analysis/constants.R")

# Read in acquisition surfaces
acquisition_surfaces = readRDS(file = "./analysis/section_4/model_output/acquisition_surfaces.rds")


# Threshold ---------------------------------------------------------------

acquisition_grid_threshold = list()

for(k in 1:length(threshold_vals)){
  cat(k, "\n")
  acquisition_grid_threshold[[k]] = list(
    PI = matrix(NA, par_resolution, length(r1_grid)),
    EI = matrix(NA, par_resolution, length(r1_grid)),
    UCB = matrix(NA, par_resolution, length(r1_grid))
  )
  for (j in 1:par_resolution) {
    # xi, quant
    for (i in 1:3) {
      # (pi, ei, ucb)
      for (h in 1:length(r1_grid)) {
        lower_lim_index = which.min(abs(theta_indexes - (pi/2 - threshold_vals[k])))
        upper_lim_index = which.min(abs(theta_indexes - (pi/2 + threshold_vals[k])))
        acquisition_threshold_values = acquisition_surfaces[[j]][[i]][h, ]
        acquisition_threshold_values[!(1:angle_resolution %in% lower_lim_index:upper_lim_index)] = -10000 # use this instead of min
        acquisition_grid_threshold[[k]][[i]][j, h] = theta_indexes[which.max(acquisition_threshold_values)]
      }
    }
  }
}

saveRDS(acquisition_grid_threshold, 
        file = "./analysis/section_5/model_output/acquisition_grid_threshold.rds")

