library(dplyr)
library(ggforce)
library(reshape2)

# Source utilities
source("./analysis/utils.R")


# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

# Example round -----------------------------------------------------------

subject_dat = session_dat %>% 
  group_by(alt_id) %>% summarise(num_sessions = n())
tens = which(session_dat$num_moves == 10)

example_session = session_dat[tens[4],]
example_events = event_dat %>%
  filter(session_id == example_session$session_id) %>%
  mutate(for_text = paste0("Move ",move,": ",round(score))) %>%
  left_join((example_session %>% 
               ungroup() %>%
               select(session_id)), by = "session_id")

# Parameters for plots

r_sig <- .01
xi_val = 1
quant_val = .95

m1_score_range = event_dat %>% 
  filter(move == 1) %>%
  with(round(range(delta_score)))

delta_r = seq(m1_score_range[1], m1_score_range[2], by = 1)
angle_resolution = 90

# normalize first move to (20,0)
optim_grid = to_cart(r = 20, theta = seq(0, pi, length.out = angle_resolution))
optim_grid$x = optim_grid$x + 20
theta_indexes = seq(0, pi, length.out = angle_resolution)

# Initialize matrices for acquisition surfaces
acquisition_grid = list()
acquisition_grid$PI_grid = matrix(NA, length(delta_r), nrow(optim_grid))
acquisition_grid$EI_grid = matrix(NA, length(delta_r), nrow(optim_grid))
acquisition_grid$UCB_grid = matrix(NA, length(delta_r), nrow(optim_grid))


for(i in 1:length(delta_r)){
  inferred_surface = infer_surface_1(x1 = 20, y1 = 0,
                                     r1 = delta_r[i] + example_events$score[1],
                                     r_sig = r_sig,
                                     beta_Sig = diag(1, nrow = 2),
                                     init_score = example_events$score[1],
                                     projection_grid = optim_grid)
  acquisition_grid$PI_grid[i,] = acquire_pi(post_pred_mu = inferred_surface$post_pred_mu,
                                        post_pred_sig = inferred_surface$post_pred_sig, 
                                        r1 = delta_r[i] + example_events$score[1], 
                                        xi_val = xi_val, 
                                        init_score = example_events$score[1])
  acquisition_grid$EI_grid[i,] = acquire_ei(post_pred_mu = inferred_surface$post_pred_mu,
                                            post_pred_sig = inferred_surface$post_pred_sig, 
                                            r1 = delta_r[i] + example_events$score[1], 
                                            xi_val = xi_val, 
                                            init_score = example_events$score[1])
  acquisition_grid$UCB_grid[i,] = acquire_ucb(post_pred_mu = inferred_surface$post_pred_mu,
                                            post_pred_sig = inferred_surface$post_pred_sig,
                                            quant = quant_val)
}


# GET maxes of acquisition surfaces for each value of delta_r 

acquisition_maxes = list()
for(i in 1:3) {
  acquisition_maxes[[i]] = cbind.data.frame(
    delta_r = delta_r,
    angle = theta_indexes[apply(acquisition_grid[[i]], 1, which.max)])
}



# Plotting ----------------------------------------------------------------


# Melt into data.frames for plotting  
acquisition_melt = list()
acquisition_melt_full = list()

for(i in 1:3) {  
  acquisition_melt[[i]] = reshape2::melt(acquisition_grid[[i]]) %>%
    dplyr::left_join(dplyr::tibble(Var1 = 1:length(delta_r),
                                   delta_r = delta_r),
                     by = "Var1") %>%
    dplyr::left_join(dplyr::tibble(Var2 = 1:length(theta_indexes),
                                   angle = theta_indexes),
                     by = "Var2") %>%
    dplyr::select(delta_r, angle, value)
  
  #  Symmetry across x-axis
  acquisition_melt_full[[i]] = acquisition_melt[[i]] %>%
    dplyr::filter(angle != 0) %>%
    dplyr::mutate(angle = -angle) %>%
    rbind(acquisition_melt[[i]])
}

acquisition_plots = list()

plot_titles = c("Move 2 Boundary PI",
                "Move 2 Boundary EI",
                "Move 2 Boundary UCB ")
sub_titles = c(bquote(xi["PI"] ~ " = " ~ .(xi_val)),
               bquote(xi["EI"] ~ " = " ~ .(xi_val)),
               bquote(p ~ " = " ~ .(quant_val)))
legend_titles = c("Prob. of Improvement", "Expected Improvement", "Upper Conf. Bound")

for(i in 1:3){
  acquisition_plots[[i]] = ggplot2::ggplot(data = acquisition_melt_full[[i]], 
                                           ggplot2::aes(x = delta_r, y = angle, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                  name = legend_titles[i]) + 
    ggplot2::geom_line(data = acquisition_maxes[[i]],
                       color = "green3",
                       size = .75,
                       ggplot2::aes(fill = NULL)) + 
    ggplot2::geom_line(data = acquisition_maxes[[i]] %>%
                         dplyr::mutate(angle = -angle),
                       color = "green3",
                       size = .75,
                       ggplot2::aes(fill = NULL)) + 
    ggplot2::geom_segment(x=5, xend = 5,
                 y = -pi, yend = pi,
                 color='black') + 
    ggplot2::geom_point(aes(x=delta_r, y=angle, fill = NULL), 
                        data = acquisition_maxes[[i]] %>%
                          filter(delta_r == 5), 
                        color = "green2",
                        shape = "*",
                        size = 10) +
    ggplot2::geom_point(aes(x=delta_r, y=angle, fill = NULL), 
                        data = acquisition_maxes[[i]] %>%
                          dplyr::mutate(angle = -angle) %>%
                          filter(delta_r == 5), 
                        color = "green2",
                        shape = "*",
                        size = 10) +
    ggplot2::scale_x_continuous(name=expression(paste(Delta, r[1])),
                                breaks=c(-30, -20, -10, 0, 10, 20, 30)) + 
    ggplot2::scale_y_continuous(name=expression(theta[2]),
                                breaks=c(-pi, -pi/2, 0, pi/2, pi),
                                labels = c(expression(-pi), 
                                           expression(-pi/2),
                                           0, 
                                           expression(pi/2), 
                                           expression(pi))) +
    ggplot2::theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(angle = 45),
          legend.key.width=unit(1,"cm")) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    ggplot2::labs(title = plot_titles[i],
         subtitle = bquote(.(sub_titles[[i]])))
}

pdf(file = "./figures/section_3/acquisition_surfaces_w_curves.pdf",
    width = 8, height = 4)
cowplot::plot_grid(acquisition_plots[[1]], 
                   acquisition_plots[[2]],
                   acquisition_plots[[3]],
                   nrow = 1)
dev.off()

