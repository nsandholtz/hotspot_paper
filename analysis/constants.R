# Set global parameters
rew_sig = 0.01

# Set parameter resolutions / bounds --------------------------------------

par_resolution = 61
score_resolution = .5
pi_xi_max = 30
ei_xi_max = 30
quant_max = .99

par_vals = list(PI = seq(0, pi_xi_max, length.out = par_resolution),
                EI = seq(0, ei_xi_max, length.out = par_resolution),
                UCB = seq(.5, quant_max, length.out = par_resolution))

# other parameters

angle_resolution = 180
optim_grid = to_cart(r = 20, theta = seq(0, pi, length.out = angle_resolution))
optim_grid$x = optim_grid$x + 20
theta_indexes = seq(0, pi, length.out = angle_resolution)

r1_grid <- seq(-34, 34, by = score_resolution)
scale_vals = seq(.01, pi/4, length.out = par_resolution)
