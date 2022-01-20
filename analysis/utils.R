
# UTILITIES -------------------------------------------------------

`%>%` = dplyr::`%>%`

# Wrap angles into [0, 2*pi]
zero_to_2pi <- function(x)
{
  x %% (2*pi)
}

# makes a single angle in [-pi, pi]
minuspi_to_pi <- function(x) {
  y <- (x + pi) %% (2*pi)
  y_neg <- which(y < 0)
  y[y_neg] <- 2*pi + y[y_neg]
  y - pi
}  

# Convert cartesian (x,y) and an origin, into radians and radius
to_polar <- function(x,y,origin_x = 0, origin_y = 0, two_pi = T){
  dist <- sqrt((x-origin_x)^2 + (y-origin_y)^2)
  if(two_pi == T){
    rad <- zero_to_2pi(atan2(y-origin_y,x-origin_x))
  } else {
    rad <- minuspi_to_pi(atan2(y-origin_y,x-origin_x))
  }
  if(length(x) > 1) return(data.frame(dist,rad))
  if(length(x) == 1) return(c(dist,rad))
}

# Calculate reward of a location, given gradient, initial score, and hotspot
calc_reward <- function(loc_x, 
                   loc_y,
                   hotspot_x,
                   hotspot_y, 
                   grad,
                   R_0, 
                   origin_x = 0, 
                   origin_y = 0){
  my_reward <- R_0 + grad*(sqrt((origin_x-hotspot_x)^2 + (origin_y-hotspot_y)^2) - 
                                  sqrt((loc_x-hotspot_x)^2 + (loc_y-hotspot_y)^2))
  return(my_reward)
}

# convert (radius, radians) to (x,y)
to_cart <- function(r,theta){
  x <- r*cos(theta)
  y <- r*sin(theta)
  return(cbind.data.frame(x,y))
}

# For a given x,y pair (rotated) and hotspot (rotated), find the optimal 
# coordinates of the next move (i.e. directly toward the hotspot)
get_opt_xy = function(x_rot, y_rot, hotspot_x_rot, hotspot_y_rot, r = 20){
  opt_angle = to_polar(x = hotspot_x_rot,
                       y = hotspot_y_rot,
                       origin_x = x_rot,
                       origin_y = y_rot)[2]
  trans_x = to_cart(r = r, theta = opt_angle)$x + x_rot
  trans_y = to_cart(r = r, theta = opt_angle)$y + y_rot
  return(c(trans_x, trans_y))
}


# INFERENCE FUNCTIONS -----------------------------------------------------

infer_surface_1 <- function(x1 = 20, y1 = 0,
                            r1,
                            r_sig, 
                            beta_Sig,
                            init_score,
                            beta_mu = c(0,0),
                            projection_grid){
  X = t(c(x1, y1))
  post_Prec <- (1/r_sig^2) * t(X) %*% X + solve(beta_Sig)
  post_Sig <- solve(post_Prec)
  post_mu <- post_Sig %*% ((1/r_sig^2) * t(X) %*% (r1 - init_score) + solve(beta_Sig) %*% beta_mu)
  
  projection_grid$post_pred_mu <- apply(projection_grid[,c("x", "y")], 1,
                                        function(x) x %*% post_mu + init_score)
  projection_grid$post_pred_sig <- apply(projection_grid[,c("x", "y")], 1,
                                         function(x) sqrt(t(x) %*% post_Sig %*% x + r_sig^2))
  return(projection_grid)
}  

infer_surface_n <- function(X,
                            r,
                            r_sig, 
                            beta_Sig,
                            init_score,
                            beta_mu = c(0,0),
                            projection_grid){
  
  post_Prec <- (1/r_sig^2) * t(X) %*% X + solve(beta_Sig)
  post_Sig <- solve(post_Prec)
  post_mu <- post_Sig %*% ((1/r_sig^2) * t(X) %*% (r - init_score) + solve(beta_Sig) %*% beta_mu)
  
  projection_grid$post_pred_mu <- apply(projection_grid[,c("x", "y")], 1,
                                        function(x) x %*% post_mu + init_score)
  projection_grid$post_pred_sig <- apply(projection_grid[,c("x", "y")], 1,
                                         function(x) sqrt(t(x) %*% post_Sig %*% x + r_sig^2))
  return(projection_grid)
}  



# ACQUISITON FUNCTIONS ----------------------------------------------------

acquire_pi = function(post_pred_mu, 
                      post_pred_sig,
                      r1,
                      xi_val,
                      init_score){
  lambda <- (post_pred_mu - max(init_score, r1) - xi_val)/post_pred_sig
  pi_surface <- pnorm(lambda)
  return(pi_surface)
}

acquire_ei = function(post_pred_mu, 
                      post_pred_sig,
                      r1,
                      xi_val,
                      init_score){
  lambda <- (post_pred_mu - max(init_score, r1) - xi_val)/post_pred_sig
  ei_surface <- post_pred_sig * (pnorm(lambda)*lambda + dnorm(lambda))
  return(ei_surface)
}

acquire_ucb = function(post_pred_mu, 
                       post_pred_sig,
                       quant
){
  ucb_surface <- qnorm(p = quant, 
                       post_pred_mu, 
                       post_pred_sig)
  return(ucb_surface)
}


acquire_pi_curve = function(x1 = 20, y1 = 0,
                            r1_vals,
                            r_sig, 
                            beta_Sig,
                            init_score = 200,
                            beta_mu = c(0,0),
                            projection_grid,
                            xi_val){
  mult_max = NA
  pi_curve = NA
  opt_angle = NA
  for(i in 1:length(r1_vals)){
    post_pred = infer_surface_1(x1 = x1,
                                y1 = y1,
                                r1 = r1_vals[i], 
                                r_sig = r_sig, 
                                beta_Sig = beta_Sig, 
                                beta_mu = beta_mu,
                                init_score = init_score, 
                                projection_grid = projection_grid)
    pi_surface = acquire_pi(post_pred_mu = post_pred$post_pred_mu,
                            post_pred_sig = post_pred$post_pred_sig, 
                            r1 = r1_vals[i], 
                            xi_val = xi_val, 
                            init_score = init_score)
    mult_max[i] = ifelse(length(which(pi_surface == max(pi_surface))) > 1, 1, 0)
    pi_max_ind = which(pi_surface == max(pi_surface))[1]
    opt_angle[i] = to_polar(projection_grid$x[pi_max_ind],
                            projection_grid$y[pi_max_ind], 
                            origin_x = x1, 
                            origin_y = y1)[2] 
  }
  return(data.frame(r1_vals, 
                    opt_angle,
                    mult_max))
}

acquire_ei_curve = function(x1 = 20, y1 = 0,
                            r1_vals,
                            r_sig, 
                            beta_Sig,
                            init_score = 200,
                            beta_mu = c(0,0),
                            projection_grid,
                            xi_val){
  mult_max = NA
  ei_curve = NA
  opt_angle = NA
  for(i in 1:length(r1_vals)){
    post_pred = infer_surface_1(x1 = x1,
                                y1 = y1,
                                r1 = r1_vals[i], 
                                r_sig = r_sig, 
                                beta_Sig = beta_Sig, 
                                beta_mu = beta_mu,
                                init_score = init_score, 
                                projection_grid = projection_grid)
    ei_surface = acquire_ei(post_pred_mu = post_pred$post_pred_mu,
                            post_pred_sig = post_pred$post_pred_sig, 
                            r1 = r1_vals[i], 
                            xi_val = xi_val, 
                            init_score = init_score)
    mult_max[i] = ifelse(length(which(ei_surface == max(ei_surface))) > 1, 1, 0)
    ei_max_ind = which(ei_surface == max(ei_surface))[1]
    opt_angle[i] = to_polar(projection_grid$x[ei_max_ind],
                            projection_grid$y[ei_max_ind], 
                            origin_x = x1, 
                            origin_y = y1)[2] 
  }
  return(data.frame(r1_vals, 
                    opt_angle,
                    mult_max))
}

acquire_ucb_curve = function(x1 = 20, y1 = 0,
                             r1_vals,
                             r_sig, 
                             beta_Sig,
                             init_score = 200,
                             beta_mu = c(0,0),
                             projection_grid,
                             quant_val){
  mult_max = NA
  ucb_curve = NA
  opt_angle = NA
  for(i in 1:length(r1_vals)){
    post_pred = infer_surface_1(x1 = x1,
                                y1 = y1,
                                r1 = r1_vals[i], 
                                r_sig = r_sig, 
                                beta_Sig = beta_Sig, 
                                beta_mu = beta_mu,
                                init_score = init_score, 
                                projection_grid = projection_grid)
    ucb_surface = acquire_ucb(post_pred_mu = post_pred$post_pred_mu,
                              post_pred_sig = post_pred$post_pred_sig, 
                              quant = quant_val)
    mult_max[i] = ifelse(length(which(ucb_surface == max(ucb_surface))) > 1, 1, 0)
    ucb_max_ind = which(ucb_surface == max(ucb_surface))[1]
    opt_angle[i] = to_polar(projection_grid$x[ucb_max_ind],
                            projection_grid$y[ucb_max_ind], 
                            origin_x = x1, 
                            origin_y = y1)[2] 
  }
  return(data.frame(r1_vals, 
                    opt_angle,
                    mult_max))
}



# INVERSE UTILS -----------------------------------------------------------

sym_wc_log_lik4 = function(data_vec, 
                           peak_vec, 
                           scale,
                           weights){
  h_theta_1 = (1/(2*pi)) * sinh(scale) /(cosh(scale) - cos(data_vec - peak_vec))
  h_theta_2 = (1/(2*pi)) * sinh(scale) /(cosh(scale) - cos(data_vec + peak_vec))
  
  sum(log((weights * h_theta_1 + (1 - weights) * h_theta_2)))
}

sym_wc_log_lik_over_grid = function(subject_rewards,
                            subject_targets,
                            r1_grid,
                            acquisition_grid_,
                            par_vals_,
                            scale_vals_, 
                            weights_){
  par_l = length(par_vals_)
  scale_l = length(scale_vals_)
  num_rounds = length(subject_rewards)
  output = matrix(NA, par_l * scale_l, 3)
  for(i in 1:par_l){
    acquisition_curve = pracma::interp2(
      r1_grid,
      par_vals_,
      Z = acquisition_grid_,
      xp = subject_rewards,
      yp = rep(par_vals_[i], num_rounds)
    )
    for(j in 1:length(scale_vals_)){
      output[scale_l * (i-1) + j, 1] = par_vals_[i]
      output[scale_l * (i-1) + j, 2] = scale_vals_[j]
      output[scale_l * (i-1) + j, 3] = sym_wc_log_lik4(data_vec = subject_targets,
                                                       peak_vec = acquisition_curve,
                                                       scale = scale_vals_[j],
                                                       weights = weights_)
    }
  }
df_output = data.frame(par_val = output[,1],
                           scale_val = output[,2],
                           log_lik = output[,3])
  return(df_output)
}

# Getting HDI intervals for the acquisition curves
# requires some creativity
get_curve_HDIs = function(r1_grid,
                          map_curve_on_grid,
                          map_scale,
                          num_samp = 100000,
                          from_to_limit = pi,
                          kern_width = .15) {
  
  full_ints = matrix(NA, nrow = length(r1_grid), 4)
  # Estimate density
  for (i in 1:length(r1_grid)) {
    try_again = TRUE
    iter = 1
    while(try_again & iter < 10){
      density_samples = suppressWarnings(
        circular::rwrappedcauchy(num_samp,
                                 mu = map_curve_on_grid[i],
                                 rho = exp(-map_scale))
      ) %>%
        ifelse(. > pi, . - 2 * pi, .)
      density_samples_reflected = c(density_samples, density_samples * -1)
      kde_estimate = density(
        density_samples_reflected,
        from = -1 * from_to_limit,
        to = from_to_limit,
        width = kern_width
      )
      hdi_ints = HDInterval::hdi(kde_estimate,
                                 allowSplit = T)
      if(length(hdi_ints) > 2){
        crap_inds = which(hdi_ints[,2] - hdi_ints[,1] < .2)
      } else {
        crap_inds = integer(0)
      }
      if(length(crap_inds) > 0){
        hdi_ints = hdi_ints[-crap_inds,]
      } 
      if(length(hdi_ints) == 2){
        full_ints[i, 1:2] = hdi_ints
        try_again = FALSE
      } else if(length(hdi_ints) == 4){
        full_ints[i, 1:2] = hdi_ints[1,]
        full_ints[i, 3:4] = hdi_ints[2,]
        try_again = FALSE
      } 
      else {
        print("Potential error")
        iter = iter + 1
      }
    }
  }
  return(full_ints)
}

