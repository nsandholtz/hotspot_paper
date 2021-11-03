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

# Read in acquisition grids
acquisition_grid = readRDS(file = "./analysis/section_4/model_output/acquisition_grid.rds")

# Read in acquisition fits 
acquisition_fits = readRDS(file = "./analysis/section_4/model_output/acquisition_fits.rds")

# Remove "burn-in" for each subject
first_n = 20
ok_sessions = session_dat %>%
  dplyr::filter(round_index > first_n,
                target_dist > session_max_dist)
ok_event_dat = event_dat %>%
  dplyr::filter(session_id %in% ok_sessions$session_id) 


# Plotting ----------------------------------------------------------------

model_cols = c("green", "blue", "red")

pdf(file = "./figures/section_4/example_acquisition_fits.pdf",
    width = 8, height = 16/5)
par(mfrow = c(1,3))
par(mar = c(5,.3,3,.3),
    oma = c(0,6,0,.1))

iter = 1
for(my_sub in c(21, 24, 28)) {
  cat(my_sub,"\r")
  # Get the rewards from move 1
  rewards_1 = event_dat %>%
    dplyr::filter(session_id %in% ok_sessions$session_id,
                  alt_id == my_sub,
                  move == 1) %>%
    dplyr::pull(delta_score) 
  
  # Calculate the targets in absolute value radians
  
  targets = event_dat %>%
    dplyr::mutate(lead_rad_relative = lead(rad_relative)) %>%
    filter(session_id %in% ok_sessions$session_id,
           alt_id == my_sub,
           move == 1) %>%
    pull(lead_rad_relative)
  
  # MAP curve on grid
  map_curve_on_grid =  pracma::interp2(
    r1_grid,
    par_vals[[acquisition_fits[[my_sub]]$MAP$acq_type]],
    Z = acquisition_grid[[acquisition_fits[[my_sub]]$MAP$acq_type]],
    xp = r1_grid,
    yp = rep(acquisition_fits[[my_sub]]$MAP$par_val,
             length(r1_grid))
  )
  
# * Credible Intervals on MAP curve ----
  
posterior_sample_indices = sample(1:nrow(acquisition_fits[[my_sub]]$posterior),
                             size = 100000,
                             prob = acquisition_fits[[my_sub]]$posterior$marg_post_norm,
                             replace = T)
num_curves = ifelse(max(posterior_sample_indices) > 2, max(posterior_sample_indices), 3)

unique_posterior_acq_curves = list()
for(i in 1:num_curves){
  unique_posterior_acq_curves[[i]] = pracma::interp2(
    r1_grid,
    par_vals[[acquisition_fits[[my_sub]]$posterior$acq_type[i]]],
    Z = acquisition_grid[[acquisition_fits[[my_sub]]$posterior$acq_type[i]]],
    xp = r1_grid,
    yp = rep(acquisition_fits[[my_sub]]$posterior$par_val[i],
             length(r1_grid))
  )
}

posterior_weights = rep(1, num_curves)
my_inds = as.numeric(names(table(posterior_sample_indices)))
posterior_weights[my_inds] = posterior_weights[my_inds] + table(posterior_sample_indices)

cred_ints = mapply(function(x,y, ...) do.call("rbind", replicate(x, y, simplify = FALSE)), 
       posterior_weights,
       unique_posterior_acq_curves) %>%
  lapply(., as.data.frame) %>%
  rbindlist(.) %>%
  apply(., MARGIN = 2, FUN = quantile, probs = c(.025, .975))

# * Prediction intervals ----------------------------------------------

my_HDIs = get_curve_HDIs(
  r1_grid = r1_grid,
  map_curve_on_grid = map_curve_on_grid,
  map_scale = acquisition_fits[[my_sub]]$MAP$scale_val,
  num_samp = 1000000,
  from_to_limit = pi,
  kern_width = .15
)

# PLOT
acq_type_ind = ifelse(acquisition_fits[[my_sub]]$MAP$acq_type == "PI", 
                      1,
                      ifelse(acquisition_fits[[my_sub]]$MAP$acq_type == "EI",
                             2,
                             3))
title_help = c("PI", "EI", "UCB")
title_quote = c(bquote(.(title_help[acq_type_ind]) ~ "," ~ xi["PI"] ~ " = " ~ .(round(acquisition_fits[[my_sub]]$MAP$par_val, digits = 2))),
                bquote(.(title_help[acq_type_ind]) ~ "," ~ xi["EI"] ~ " = " ~ .(round(acquisition_fits[[my_sub]]$MAP$par_val, digits = 2))),
                bquote(.(title_help[acq_type_ind]) ~ "," ~ hat(p) ~ " = " ~ .(round(acquisition_fits[[my_sub]]$MAP$par_val, digits = 2))))
# * Plot Subject Fit ----
plot(rewards_1, targets,
     pch = 16,
     cex.lab = 1.5,
     ylab = expression(theta[2]),
     xlab = expression(paste(Delta,r[1])), 
     ylim = c(-pi,pi),
     xlim = range(r1_grid),
     axes = F, type = "n")
if(iter == 1) {
  axis(2, col = "gray", las = 1, at = c(-pi, -3 * pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3 * pi/4, pi),
     labels = c(expression(-pi), expression(-3 * pi/4), expression(-pi/2), 
                expression(-pi/4), 0, expression(pi/4), expression(pi/2), 
                expression(3 * pi/4), expression(pi)))
  mtext(expression(theta[2]),side = 2, outer = T, line = 3, cex = 1)
}
axis(1, col = "gray", at = seq(-30, 30, by = 10))
box(col = "gray", lwd = 2)
abline(h = c(-pi, -3 * pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3 * pi/4, pi),
       col = "gray",
       lty = "dotted")
abline(v = seq(-30, 30, by = 5),
       col = "gray",
       lty = "dotted")

# TEXT
title(paste("Subject", my_sub, "MAP Acquisition Curve"), line = .5, cex.main = 1)

# Plot prediction intervals
unimodal_ind = match(NA,my_HDIs[,3])
if(!is.na(unimodal_ind) & unimodal_ind != 1){
  x_poly = c(r1_grid, 
             rev(r1_grid[unimodal_ind:length(r1_grid)]),
             rev(r1_grid[1:(unimodal_ind - 1)]),
             r1_grid[1:(unimodal_ind - 1)],
             rev(r1_grid[1:(unimodal_ind - 1)]))   
  y_poly = c(my_HDIs[,1], 
             rev(my_HDIs[unimodal_ind:length(r1_grid), 2]),
             rev(my_HDIs[1:(unimodal_ind - 1), 4]),
             my_HDIs[1:(unimodal_ind - 1), 3],
             rev(my_HDIs[1:(unimodal_ind - 1), 2]))
  
  polygon(x = x_poly, y = y_poly, 
          col = gray(.75, alpha = .5),
          border = NA)
} else {
  x_poly_upper = c(r1_grid, 
                   rev(r1_grid))   
  y_poly_upper = c(my_HDIs[,3], 
                   rev(my_HDIs[, 4]))
  polygon(x = x_poly_upper, y = y_poly_upper, 
          col = gray(.75, alpha = .5),
          border = NA)
  
  x_poly_lower = c(r1_grid, 
                   rev(r1_grid))   
  y_poly_lower = c(my_HDIs[,1], 
                   rev(my_HDIs[, 2]))
  polygon(x = x_poly_lower, y = y_poly_lower, 
          col = gray(.75, alpha = .5),
          border = NA)
}

# Plot cred ints
x_poly_upper = c(r1_grid, 
                 rev(r1_grid))   
y_poly_upper = c(cred_ints[1, r1_grid < 0],
                 cred_ints[2, r1_grid >= 0],
                 rev(cred_ints[1, r1_grid >= 0]),
                 rev(cred_ints[2, r1_grid < 0]))
polygon(x = x_poly_upper, y = y_poly_upper, 
        col = gray(.25, alpha = .5),
        border = NA)

x_poly_lower = c(r1_grid, 
                 rev(r1_grid))   
y_poly_lower = c(-cred_ints[1, r1_grid < 0],
                 -cred_ints[2, r1_grid >= 0],
                 -rev(cred_ints[1, r1_grid >= 0]),
                 -rev(cred_ints[2, r1_grid < 0]))

polygon(x = x_poly_lower, y = y_poly_lower, 
        col = gray(.25, alpha = .5),
        border = NA)

# plot raw data
points(rewards_1, targets,
       pch = 16,
       cex = .4)
# Plot map curves

lines(r1_grid, 
      map_curve_on_grid, 
      lwd = 1,
      lty = 1,
      col = ifelse(acquisition_fits[[my_sub]]$MAP$acq_type == "PI", 
                   model_cols[1],
                   ifelse(acquisition_fits[[my_sub]]$MAP$acq_type == "EI",
                          model_cols[2],
                          model_cols[3])))
lines(r1_grid, 
      -map_curve_on_grid, 
      lwd = 1,
      lty = 1,
      col = ifelse(acquisition_fits[[my_sub]]$MAP$acq_type == "PI", 
                   model_cols[1],
                   ifelse(acquisition_fits[[my_sub]]$MAP$acq_type == "EI",
                          model_cols[2],
                          model_cols[3])))

text(2,-pi,
     bquote(.(title_quote[[acq_type_ind]])),
     pos = 4,cex = 1
)
iter = iter + 1
}
dev.off()
