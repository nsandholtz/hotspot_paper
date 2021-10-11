
# Utility Functions -------------------------------------------------------

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

# Convert (x,y) and an Origin, into radians and radius
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

# Convert radians to degrees
rad2deg <- function(rad) {(rad * 180) / (pi)}

# convert (radius, radians) to (x,y)
to_cart <- function(r,theta){
  x <- r*cos(theta)
  y <- r*sin(theta)
  return(cbind.data.frame(x,y))
}


get_opt_angle = function(x_rot, y_rot, hotspot_x_rot, hotspot_y_rot){
  opt_angle = to_polar(x = hotspot_x_rot,
                       y = hotspot_y_rot,
                       origin_x = r_rot,
                       origin_y = y_rot)[2]
  return(opt_angle)
}

get_opt_xy = function(x_rot, y_rot, hotspot_x_rot, hotspot_y_rot, r = 20){
  opt_angle = to_polar(x = hotspot_x_rot,
                       y = hotspot_y_rot,
                       origin_x = x_rot,
                       origin_y = y_rot)[2]
  trans_x = to_cart(r = r, theta = opt_angle)$x + x_rot
  trans_y = to_cart(r = r, theta = opt_angle)$y + y_rot
  return(c(trans_x, trans_y))
}

# Function to compute log-likelihood
aquisition_log_lik = function(data_vec, 
                              mean_vec, 
                              sig2,
                              neg_side = F){
  if(neg_side == FALSE){
    sum(VGAM::dfoldnorm(data_vec, 
                        mean = as.numeric(mean_vec), 
                        sd = sig2,
                        log = T))
  } else {
    sum(VGAM::dfoldnorm(pi - data_vec, 
                        mean = pi - as.numeric(mean_vec), 
                        sd = sig2,
                        log = T))
  }
}



# FASTER VERSIONS ---------------------------------------------------------

infer_surface_1 <- function(x1 = 20, y1 = 0,
                            r1,
                            r_sig, 
                            beta_Sig,
                            init_score = 200,
                            beta_mu = c(0,0),
                            projection_grid){
  
  #projection_grid$x = projection_grid$x + x1
  #projection_grid$y = projection_grid$y + y1
  
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
                            init_score = 200,
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

acquire_pi = function(post_pred_mu, 
                      post_pred_sig,
                      r1,
                      xi_val,
                      init_score = 200){
  lambda <- (post_pred_mu - max(init_score, r1) - xi_val)/post_pred_sig
  pi_surface <- pnorm(lambda)
  return(pi_surface)
}

acquire_ei = function(post_pred_mu, 
                      post_pred_sig,
                      r1,
                      xi_val,
                      init_score = 200){
  lambda <- (post_pred_mu - max(init_score, r1) - xi_val)/post_pred_sig
  ei_surface <- post_pred_sig * (pnorm(lambda)*lambda + dnorm(lambda))
  return(ei_surface)
}

acquire_ucb = function(post_pred_mu, 
                       post_pred_sig,
                       #r1,
                       #init_score = 200,
                       quant
){
  #ucb_surface <- post_pred_mu + kappa * post_pred_sig
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

# SLOWER VERSIONS ---------------------------------------------------------



# Posterior pred surface after move 1:
surface_move1 <- function(x1 = 20, y1 = 0,
                          r1,
                          r_sig, 
                          beta_Sig,
                          init_score = 200,
                          beta_mu = c(0,0),
                          click_radius = 20,
                          length_of_grid = 45,
                          projection_grid = expand.grid(x = seq(-click_radius, click_radius,
                                                                length.out = length_of_grid),
                                                        y = seq(-click_radius, click_radius,
                                                                length.out = length_of_grid))
) {
  
  X = t(c(x1, y1))
  post_Prec <- (1/r_sig^2) * t(X) %*% X + solve(beta_Sig)
  post_Sig <- solve(post_Prec)
  post_mu <- post_Sig %*% ((1/r_sig^2) * t(X) %*% (r1 - init_score) + solve(beta_Sig) %*% beta_mu)
  
  # Project the posterior predictive mean onto grid
  
  polar_projection_grid <- cbind.data.frame(projection_grid, 
                                            to_polar(projection_grid$x, 
                                                     projection_grid$y)) %>% 
    dplyr::filter(dist <= click_radius+.1) %>%
    dplyr::mutate(x = x + x1,
                  y = y + y1)
  polar_projection_grid$post_pred_mu <- apply(polar_projection_grid[,c("x","y")], 1,
                                              function(x) x %*% post_mu + init_score)
  polar_projection_grid$post_pred_sig <- apply(polar_projection_grid[,c("x","y")], 1,
                                               function(x) sqrt(t(x) %*% post_Sig %*% x + r_sig^2))
  return(polar_projection_grid)
}

# Calculate PI, EI, and UCB surfaces:

acq_surf <- function(x1 = 20, y1 = 0,
                     r1,
                     r_sig, 
                     beta_Sig,
                     init_score = 200,
                     xi_val = 0,
                     kappa = 0,
                     omega = 0,
                     beta_mu = c(0,0),
                     click_radius = 20,
                     length_of_grid = 45,
                     projection_grid = expand.grid(x = seq(-click_radius, click_radius,
                                                           length.out = length_of_grid),
                                                   y = seq(-click_radius, click_radius,
                                                           length.out = length_of_grid))) {
  surf_inf <- surface_move1(x1 = x1, y1 = y1,
                            r1 = r1, 
                            r_sig = r_sig[1],
                            beta_Sig = beta_Sig,
                            init_score = init_score,
                            beta_mu = beta_mu,
                            click_radius = click_radius,
                            length_of_grid = length_of_grid,
                            projection_grid = projection_grid)
  lambda <- (surf_inf$post_pred_mu - max(init_score, r1) - xi_val)/surf_inf$post_pred_sig
  penalty <- (40 - apply(surf_inf[,c("x","y")],
                         1,
                         FUN = function(x) dist(rbind(x,c(0,0)))))/40
  surf_inf$prob_imp <- pnorm(lambda) - omega * penalty
  surf_inf$exp_imp <- surf_inf$post_pred_sig * (pnorm(lambda)*lambda + dnorm(lambda)) - omega * penalty
  surf_inf$ucb_imp <- surf_inf$post_pred_mu + kappa * surf_inf$post_pred_sig - omega * penalty
  return(surf_inf)
}


# Functions to pass to optim to find optimal tuning params (xi, kappa)

pi_acq_optim = function(my_pars, rewards, r_sig, optim_grid, targets){
  pi_model_optimums <- NA
  for(i in 1:length(rewards)){
    a_surfaces = acq_surf(x1 = 20, y1 = 0,
                          r1 = rewards[i], 
                          r_sig = r_sig,
                          beta_Sig = diag(2, nrow = 2),
                          xi_val = my_pars[1],
                          omega = my_pars[2],
                          init_score = 200,
                          projection_grid = optim_grid)
    pi_model_optimums[i] = a_surfaces$rad[which.max(a_surfaces$prob_imp)]
  }
  pi_xi_error_sums = sum(abs(pi_model_optimums - targets))
  return(pi_xi_error_sums)
}

ei_acq_optim = function(my_pars, rewards, r_sig, optim_grid, targets){
  ei_model_optimums <- NA
  for(i in 1:length(rewards)){
    a_surfaces = acq_surf(x1 = 20, y1 = 0,
                          r1 = rewards[i], 
                          r_sig = r_sig,
                          beta_Sig = diag(2, nrow = 2),
                          xi_val = my_pars[1],
                          omega = my_pars[2],
                          init_score = 200,
                          projection_grid = optim_grid)
    ei_model_optimums[i] = a_surfaces$rad[which.max(a_surfaces$exp_imp)]
  }
  ei_xi_error_sums = sum(abs(ei_model_optimums - targets))
  return(ei_xi_error_sums)
}

ucb_acq_optim = function(my_pars,
                         rewards, 
                         r_sig, 
                         optim_grid, 
                         targets){
  ucb_model_optimums <- NA
  for(i in 1:length(rewards)){
    a_surfaces = acq_surf(x1 = 20, y1 = 0,
                          r1 = rewards[i], 
                          r_sig = r_sig,
                          beta_Sig = diag(2, nrow = 2),
                          kappa = my_pars[1],
                          omega = my_pars[2],
                          init_score = 200,
                          projection_grid = optim_grid)
    ucb_model_optimums[i] = a_surfaces$rad[which.max(a_surfaces$ucb_imp)]
  }
  ucb_kappa_error_sums = sum(abs(ucb_model_optimums - targets))
  return(ucb_kappa_error_sums)
}

subject_opt_curve_finder = function(rewards,
                                    targets,
                                    PI = T, EI = T, UCB = T){
  
  # Restrict grid to outer edge of circle 
  
  optim_grid = to_cart(r = 20, theta = seq(0, pi, length.out = 180))
  
  # Select grid search values for tuning parameters
  # to pass to optim as an initial value
  xi_vals <- seq(0,50, length.out = 5)
  kappa_vals <- seq(0,3, length.out = 5)
  omega_vals <- 0 # always fix this initial value to 0
  
  # Find good initial value for tuning parameter
  cat("Finding good initial values\n")
  PI_error_sums <- NA
  EI_error_sums <- NA
  UCB_error_sums <- NA
  
  for(j in 1:length(xi_vals)){
    PI_model_optimums <- NA
    EI_model_optimums <- NA
    UCB_model_optimums <- NA
    
    for(i in 1:length(rewards)){
      a_surfaces = acq_surf(x1 = 20, y1 = 0,
                            r1 = rewards[i], 
                            r_sig = r_sig,
                            beta_Sig = diag(2, nrow = 2),
                            xi_val = xi_vals[j],
                            kappa = kappa_vals[j],
                            omega = omega_vals,
                            init_score = 200,
                            projection_grid = optim_grid)
      PI_model_optimums[i] = a_surfaces$rad[which.max(a_surfaces$prob_imp)]
      EI_model_optimums[i] = a_surfaces$rad[which.max(a_surfaces$exp_imp)]
      UCB_model_optimums[i] = a_surfaces$rad[which.max(a_surfaces$ucb_imp)]
    }
    PI_error_sums[j] = sum(abs(PI_model_optimums - targets))
    EI_error_sums[j] = sum(abs(EI_model_optimums - targets))
    UCB_error_sums[j] = sum(abs(UCB_model_optimums - targets))
  }
  initial_PI = xi_vals[which.min(PI_error_sums)]
  initial_EI = xi_vals[which.min(EI_error_sums)]
  initial_UCB = kappa_vals[which.min(UCB_error_sums)]
  
  optim_PI = optim_EI = optim_UCB = NA
  if(PI){
    cat("Optimizing PI tuning paramaters\n")
    optim_PI =  suppressWarnings(optim(par = c(initial_PI, omega_vals), 
                                       fn = pi_acq_optim,
                                       rewards = rewards,
                                       r_sig = r_sig,
                                       optim_grid = optim_grid,
                                       targets = targets)$par)
  }
  if(EI){
    cat("Optimizing EI tuning paramaters\n")
    optim_EI =  suppressWarnings(optim(par = c(initial_EI, omega_vals), 
                                       fn = ei_acq_optim,
                                       rewards = rewards,
                                       r_sig = r_sig,
                                       optim_grid = optim_grid,
                                       targets = targets)$par)
  }
  if(UCB){
    cat("Optimizing UCB tuning paramaters\n")
    optim_UCB =  suppressWarnings(optim(par = c(initial_UCB, omega_vals), 
                                        fn = ucb_acq_optim,
                                        rewards = rewards,
                                        r_sig = r_sig,
                                        optim_grid = optim_grid,
                                        targets = targets)$par)
  }
  return(list(optim_PI,
              optim_EI,
              optim_UCB))
}

subject_opt_curve_plotter = function(rewards,
                                     targets,
                                     subject_label,
                                     r_sig,
                                     PI_optims, 
                                     EI_optims, 
                                     UCB_optims,
                                     reward_resolution = seq(-90,90, 
                                                             length.out = 300) + 200) {
  plot(0,0,type="n",
       xlim=c(-90, 90),ylim=c(0,180),
       xlab = expression(paste(Delta," reward")),
       ylab = "Angle of second move",
       main = paste("Subject", 1,"\nBest Fit Explore-Exploit Curves"),
       axes = F)
  abline(h = c(0,45,90,135,180),col = "gray",lty = "dotted")
  grid(ny = 0)
  axis(1)
  axis(2, at = c(0,45,90,135,180))
  points(rewards - 200,
         rad2deg(targets),
         pch=16)
  #abline(h = 0, lty = 2)
  
  optim_grid = to_cart(r = 20, theta = seq(0, pi, length.out = 180))
  
  pi_optim_curve <- NA
  ei_optim_curve <- NA
  ucb_optim_curve <- NA
  
  if(!is.null(PI_optims)){
    for(i in 1:length(reward_resolution)){
      a_surfaces = acq_surf(x1 = 20, y1 = 0,
                            r1 = reward_resolution[i], 
                            r_sig = r_sig,
                            beta_Sig = diag(2, nrow = 2),
                            xi_val = PI_optims[1],
                            omega = PI_optims[2],
                            init_score = 200,
                            projection_grid = optim_grid)
      pi_optim_curve[i] = a_surfaces$rad[which.max(a_surfaces$prob_imp)]
    }
    lines(reward_resolution - 200,
          rad2deg(pi_optim_curve),
          col = "green4",
          lwd = 2)
  }
  if(!is.null(EI_optims)){
    for(i in 1:length(reward_resolution)){
      a_surfaces = acq_surf(x1 = 20, y1 = 0,
                            r1 = reward_resolution[i], 
                            r_sig = r_sig,
                            beta_Sig = diag(2, nrow = 2),
                            xi_val = EI_optims[1],
                            omega = EI_optims[2],
                            init_score = 200,
                            projection_grid = optim_grid)
      ei_optim_curve[i] = a_surfaces$rad[which.max(a_surfaces$exp_imp)]
    }
    lines(reward_resolution - 200,
          rad2deg(ei_optim_curve),
          col = "blue",
          lwd = 2)
  }
  if(!is.null(UCB_optims)){
    for(i in 1:length(reward_resolution)){
      a_surfaces = acq_surf(x1 = 20, y1 = 0,
                            r1 = reward_resolution[i], 
                            r_sig = r_sig,
                            beta_Sig = diag(2, nrow = 2),
                            kappa = UCB_optims[1],
                            omega = UCB_optims[2],
                            init_score = 200,
                            projection_grid = optim_grid)
      ucb_optim_curve[i] = a_surfaces$rad[which.max(a_surfaces$ucb_imp)]
    }
    lines(reward_resolution - 200,
          rad2deg(ucb_optim_curve),
          col = "red",
          lwd = 2)
  }
  legend("topright",
         c(as.expression(bquote(.("PI:") ~ xi ~ .(paste0("= ",round(PI_optims[1], digits = 2),",")) ~ omega ~ .(paste("=",round(PI_optims[2], digits = 2))))),
           as.expression(bquote(.("EI:") ~ xi ~ .(paste0("= ",round(EI_optims[1], digits = 2),",")) ~ omega ~ .(paste("=",round(EI_optims[2], digits = 2))))),
           as.expression(bquote(.("UCB:") ~ kappa ~ .(paste0("= ",round(UCB_optims[1], digits = 2),",")) ~ omega ~ .(paste("=",round(UCB_optims[2], digits = 2)))))),
         col = c("green4",  "blue", "red"),
         pch = 15)
}



# Calculate the maximum possible score for a round, given the origin, 
# reward gradient, number of moves, and the hotspot location
max_score = function(hotspot_x, hotspot_y, grad,
                     origin_x = 960, origin_y = 540, 
                     max_moves = 11
){
  
  hotspot_dist <- sqrt((origin_x-hotspot_x)^2 + (origin_y-hotspot_y)^2)
  if(hotspot_dist <= 20*max_moves){
    max_x <- hotspot_x
    max_y <- hotspot_y
    max_reward <- reward(max_x,max_y,
                         hotspot_x = hotspot_x,
                         hotspot_y = hotspot_y,
                         origin_x = origin_x,
                         origin_y = origin_y,
                         grad = grad)
    return(c(max_x, max_y, max_reward))
  } else {
    theta <- atan2(hotspot_y-origin_y,hotspot_x-origin_x)
    if((hotspot_y > origin_y) & (hotspot_x > origin_x)){
      max_y <- sin(theta)*(20*max_moves) + origin_y
      max_x <- cos(theta)*(20*max_moves) + origin_x
    } else if((hotspot_y > origin_y) & (hotspot_x < origin_x)) {
      max_y <- sin(pi-theta)*(20*max_moves) + origin_y
      max_x <- -cos(pi-theta)*(20*max_moves) + origin_x
    } else if((hotspot_y < origin_y) & (hotspot_x < origin_x)) {
      max_y <- sin(-pi-theta)*(20*max_moves) + origin_y
      max_x <- -cos(-pi-theta)*(20*max_moves) + origin_x
    } else {
      max_y <- sin(theta)*(20*max_moves) + origin_y
      max_x <- cos(theta)*(20*max_moves) + origin_x
    }
    
    max_reward <- reward(max_x,max_y,
                         hotspot_x = hotspot_x,
                         hotspot_y = hotspot_y,
                         origin_x = origin_x,
                         origin_y = origin_y,
                         grad = grad)
    return(c(max_x, max_y, max_reward))
  }
}

# Calculates the reward of a mouse click, given
# the gradient, origin, and hotspot location.
reward <- function(mouse_x, mouse_y, hotspot_x, hotspot_y, 
                   R_0 = 200, grad = runif(1,15,75), 
                   origin_x = 960, origin_y = 540){
  my_reward <- R_0 + (grad/20)*(sqrt((origin_x-hotspot_x)^2 + (origin_y-hotspot_y)^2) - 
                                  sqrt((mouse_x-hotspot_x)^2 + (mouse_y-hotspot_y)^2))
  return(my_reward)
}


# REWARD SURFACE METROPOLIS ESTIMATION --------------------------------------------

sample_post = function(reward,
                       sig_r,
                       init_beta = c(0, 1),
                       sig_beta,
                       gam_a = 9,
                       gam_b = 4,
                       N_iter) {
  # * Initialize ----
  
  N_acc <- 0 # Number of accepted candidates
  
  out <- matrix(NA, N_iter, 3)
  init_lik <- dgamma(x = sqrt(sum(init_beta ^ 2)),
                     shape = gam_a,
                     rate = gam_b, log = T) + dnorm(x = reward,
                                           mean = 200 + init_beta[1] * 20,
                                           sd = sig_r, log = T)
  out[1, ] <- c(init_beta,
                init_lik)
  
  # ITERATE METROPOLIS
  
  # * Do metropolis ----
  
  for (i in 2:N_iter) {
    out[i, ] <- out[i - 1, ] # Initialize at previous value
    
    # beta candidate
    beta_cand = mvtnorm::rmvnorm(n = 1, mean = out[i, 1:2], sigma = sig_beta)
    if (beta_cand[2] >= 0) {
      accept <- (
        dgamma(
          x = sqrt(sum(beta_cand ^ 2)),
          shape = gam_a,
          rate = gam_b,
          log = T
        ) + dnorm(
          x = reward,
          mean = 200 + beta_cand[1] * 20,
          sd = sig_r,
          log = T
        )
      ) - out[i, 3]
      u <- runif(1, 0, 1)
      if (log(u) < accept) {
        out[i, 1:2] <- beta_cand
        out[i, 3] <- (
          dgamma(
            x = sqrt(sum(beta_cand ^ 2)),
            shape = gam_a,
            rate = gam_b,
            log = T
          ) + dnorm(
            x = reward,
            mean = 200 + beta_cand[1] * 20,
            sd = sig_r,
            log = T
          )
        )
        N_acc[1] <- N_acc[1] + 1
      }
    }
  }
  return(list(N_acc = N_acc,
              out = out))
}


# Plotting Functions ----------------------------------------------------

# Plot a round from the raw data:
plot_round <- function(subject,
                       round,
                       dat = all_dat,
                       info = dat_info){
  sub_info <- dplyr::filter(info,subj == subject & block_ct == round)
  sub_data <- dplyr::filter(dat,subj == subject & block_ct == round)
  plot(1,2,type="n",xlim=c(960-450,960+450),ylim=c(540-450,540+450),asp = 1,
       axes = F,main=paste("Subject: ",subject," of ",max(info$subj),"\n Round: ",round," of ",max(dat$block_ct[dat$subj == subject]),sep = ""),ylab = "",xlab = "")
  plotrix::draw.circle(960,540,450)
  points(sub_info$hotspot_cx, sub_info$hotspot_cy,col="red",pch=19,cex=.4)
  points(sub_info$hotspot_cx, sub_info$hotspot_cy,col="red",cex=1.5)
  points(sub_info$hotspot_cx, sub_info$hotspot_cy,col="red",cex=2.5)
  points(sub_data$cursor_cx,sub_data$cursor_cy,type="l",col="darkgray")
  points(sub_data$cursor_cx,sub_data$cursor_cy,pch=19,cex=.5,col="darkgray")
  points(sub_data$cursor_cx[1],sub_data$cursor_cy[1],pch=19,cex=.5,col="blue")
}


plot_fit = function(subject,
                    subject_dat,
                    fit_obj,
                    burn,
                    r1,
                    tuning_mins = c(0,0,.5),
                    tuning_maxs = c(25,25,.975),
                    tuning_vals,
                    tuning_grid,
                    block_ct_rid = 10,
                    trace_plot = T){
  
  # Get the rewards from move 1
  sub_rewards = subject_dat %>%
    dplyr::filter(subj == subject,
                  block_ct > block_ct_rid,
                  dots_shown <= 275,
                  dots_shown >= 125,
                  trial == 2) %>%
    dplyr::pull(dots_shown)
  
  # Calculate the targets in absolute value radians
  
  target_inds = subject_dat %>%
    with(which(subj == subject &
                 block_ct > block_ct_rid &
                 dots_shown <= 275 &
                 dots_shown >= 125 &
                 trial == 2))
  
  sub_targets = subject_dat[target_inds + 1, ] %>% 
    dplyr::mutate(abs_rad = abs(rot_rad)) %>%
    dplyr::pull(abs_rad)
  
  # par(mfrow = c(1,2))
  # plot(fit_obj$out[,7], type = "l",
  #      main = paste0(expression("Mixing Plot\n", 
  #                               sigma)),
  #      ylab = expression(sigma))
  # plot(fit_obj$out[,8], type = "l",
  #      main = paste0("Mixing Plot\n",
  #                    "Log-likelihood"),
  #      ylab = "Log-likelihood")
  # 
  title_help <- c("PI", "EI", "UCB")
  # par(mfrow = c(2,3))
  # if(trace_plot){
  #   # Trace chains
  #   for(i in 1:3){
  #     plot(fit_obj$out[, i],
  #          type = "l",
  #          ylab = ifelse(i == 3, expression(kappa), expression(xi)),
  #          main = paste("Mixing Plot\n", title_help[i], "Tuning Parameter"))
  #   }
  #   for(i in 1:3){
  #     plot(fit_obj$out[, i+3],
  #          type = "l",
  #          ylab = ifelse(i == 3, expression(rho[kappa]), expression(rho[xi])),
  #          main = paste("Mixing Plot\n", title_help[i], "Mixture Component"))
  #   }
  # }
  
  # for(i in 1:3){
  #   plot(density(fit_obj$out[burn:nrow(fit_obj$out), i], 
  #                from = tuning_mins[i], to = tuning_maxs[i]),
  #        main = paste("MCMC Density\n", title_help[i], "Tuning Parameter"),
  #        xlab = ifelse(i == 3, expression(kappa), expression(xi)))
  #   abline(v = mean(fit_obj$out[burn:nrow(fit_obj$out), i]), 
  #          col = "red", 
  #          lty = 2)
  # }
  # for(i in 1:3){
  #   plot(density(fit_obj$out[burn:nrow(fit_obj$out), i+3], from = 0, to = 1),
  #        main = paste("MCMC Density\n", title_help[i], "Mixture Component"),
  #        xlab = ifelse(i == 3, expression(rho[kappa]), expression(rho[xi])),
  #        xlim = c(0,1))
  #   abline(v = mean(fit_obj$out[burn:nrow(fit_obj$out), i+3]), 
  #          col = "red",
  #          lty = 2)
  # }
  #par(mfrow = c(2,1))
  
  plot(density(fit_obj$out[burn:nrow(fit_obj$out), 8]),
       main = paste("Subject", subject, "Mixture Model Log Likelihood"),
       xlab = "Log Likelihood")
  
  post_mean_vec = apply(fit_obj$out[burn:nrow(fit_obj$out), ], 2, mean)
  posterior_tuning_mean <- matrix(NA, 3, length(sub_rewards))
  
  for(k in 1:3){
    posterior_tuning_mean[k,] =  pracma::interp2(r1, 
                                                 tuning_vals[[k]], 
                                                 Z = tuning_grid[[k]],
                                                 xp = sub_rewards, yp = rep(post_mean_vec[k], 
                                                                            length(sub_rewards))) 
  }
  
  norm_post_mixture = post_mean_vec[4:6]/sum(post_mean_vec[4:6])
  posterior_mixture_mean = c(norm_post_mixture %*% posterior_tuning_mean)
  
  # Calculate CIs 
  posterior_mixture_variates <- matrix(NA, 
                                       nrow(fit_obj$out),
                                       length(sub_rewards))
  
  for(i in 1:nrow(fit_obj$out)){
    tuning_variates <- matrix(NA, 3, length(sub_rewards))
    for(k in 1:3){
      tuning_variates[k,] =  pracma::interp2(r1, 
                                             tuning_vals[[k]], 
                                             Z = tuning_grid[[k]],
                                             xp = sub_rewards,
                                             yp = rep(fit_obj$out[i, k], 
                                                      length(sub_rewards))) 
    }
    posterior_mixture_variates[i,] = c(fit_obj$out[i, 4:6] %*% tuning_variates)
  }
  lower_bound = apply(posterior_mixture_variates[-c(1:burn),], 2, quantile, probs = .025)
  upper_bound = apply(posterior_mixture_variates[-c(1:burn),], 2, quantile, probs = .975)
  
  plot(sub_rewards, sub_targets,
       pch = 16,
       cex = .7,
       main = paste("Subject", subject, "\n",
                    "Observed Data and Best Fit Mixture Aquisition Function"),
       ylab = "Move 2 Angle (in radians)",
       xlab = "Move 1 Reward", 
       ylim = c(0,pi),
       axes = F)
  axis(2, at = c(0, pi/4, pi/2, 3 * pi/4, pi),
       labels = c(0, expression(pi/4), expression(pi/2), 
                  expression(3 * pi/4), expression(pi)))
  axis(1, at = seq(120, 280, by = 5))
  box()
  lines(sub_rewards[order(sub_rewards)], 
        posterior_mixture_mean[order(sub_rewards)], 
        lwd = 3,
        lty = 1,
        col = "red")
  lines(sub_rewards[order(sub_rewards)], 
        lower_bound[order(sub_rewards)], 
        lwd = 2,
        lty = 2,
        col = "blue")
  lines(sub_rewards[order(sub_rewards)], 
        upper_bound[order(sub_rewards)], 
        lwd = 2,
        lty = 2,
        col = "blue")
  abline(v = seq(120, 280, by = 5), 
         lty = 2,
         col = "lightgray")
  abline(h =  c(0, pi/4, pi/2, 3 * pi/4, pi),
         lty = 2,
         col = "lightgray")
  abline(v = 200, lty = 2)
  abline(h = c(0, pi))
}

