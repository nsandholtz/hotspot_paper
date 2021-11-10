# Set global parameters
rew_sig = 0.01

# Set parameter bounds / resolutions --------------------------------------

pi_xi_max = 30
ei_xi_max = 30
quant_max = .99

par_resolution = 61
par_vals = list(PI = seq(0, pi_xi_max, length.out = par_resolution),
                EI = seq(0, ei_xi_max, length.out = par_resolution),
                UCB = seq(.5, quant_max, length.out = par_resolution))

# move 2 angle resolution
angle_resolution = 180
optim_grid = to_cart(r = 20, theta = seq(0, pi, length.out = angle_resolution))
optim_grid$x = optim_grid$x + 20
theta_indexes = seq(0, pi, length.out = angle_resolution)

# reward resolution
score_resolution =  .5
r1_grid <- seq(-34, 34, by = score_resolution)

# wrapped cauchy scale parameter resolution
scale_resolution = 61
scale_vals = seq(.01, pi/4, length.out = scale_resolution)

# Exploration threshold parameter tau
threshold_resolution = 30
threshold_vals = seq(0,pi/2, length.out = threshold_resolution)
