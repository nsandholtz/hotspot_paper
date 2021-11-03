library(tidyverse)

# Source utilities
source("./analysis/utils.R")
source("./analysis/constants.R")

acquisition_surfaces = list()

for(i in 1:length(r1_grid)) {
  cat(i, "\n")
  
  inferred_surface = infer_surface_1(
    x1 = 20,
    y1 = 0,
    r1 = r1_grid[i],
    r_sig = rew_sig,
    beta_Sig = diag(1, nrow = 2),
    init_score = 0,
    projection_grid = optim_grid
  )
  
  # CALCULATE ACQUISITION SURFACES ---
  for (j in 1:par_resolution) {
    if (i == 1) {
      acquisition_surfaces[[j]] = list(
        pi_surface = matrix(NA, length(r1_grid), angle_resolution),
        ei_surface = matrix(NA, length(r1_grid), angle_resolution),
        ucb_surface = matrix(NA, length(r1_grid), angle_resolution)
      )
    }
    for(k in 1:3){
      if(k %in% c(1,2)){
        acquisition_surfaces[[j]][[k]][i, ] = acquire_pi(
          post_pred_mu = inferred_surface$post_pred_mu,
          post_pred_sig = inferred_surface$post_pred_sig,
          r1 = r1_grid[i],
          xi_val = par_vals[[k]][j],
          init_score = 0
        )
      } else {
        acquisition_surfaces[[j]][[k]][i, ] = acquire_ucb(
          post_pred_mu = inferred_surface$post_pred_mu,
          post_pred_sig = inferred_surface$post_pred_sig,
          quant = par_vals[[k]][j]
        )
      }
    }
  }
}

# * # FIND THE MAXES ----

acquisition_grid = list(
  PI = matrix(NA, par_resolution, length(r1_grid)),
  EI = matrix(NA, par_resolution, length(r1_grid)),
  UCB = matrix(NA, par_resolution, length(r1_grid))
)

for(j in 1:par_resolution) {
  acquisition_maxes = list()
  for (i in 1:3) {
    acquisition_maxes[[i]] = cbind.data.frame(r1 = r1_grid,
                                              angle = theta_indexes[apply(acquisition_surfaces[[j]][[i]], 1, which.max)])
    acquisition_grid[[i]][j, ] = acquisition_maxes[[i]]$angle
  }
}

saveRDS(acquisition_surfaces, 
        file = "./analysis/section_4/model_output/acquisition_surfaces.rds")
saveRDS(acquisition_grid, 
        file = "./analysis/section_4/model_output/acquisition_grid.rds")


