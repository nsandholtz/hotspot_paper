library(tidyverse)
library(ggforce)

# Source utilities
source("./analysis/utils.R")
source("./analysis/constants.R")

# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

# Global parameters
xi_val = 1
quant_val = .95

acquisition_grid = list()
acquisition_grid$PI = matrix(, length(r1_grid), nrow(optim_grid))
acquisition_grid$EI = matrix(, length(r1_grid), nrow(optim_grid))
acquisition_grid$UCB = matrix(, length(r1_grid), nrow(optim_grid))


for(i in 1:length(r1_grid)){
  inferred_surface = infer_surface_1(x1 = 20, y1 = 0,
                                     r1 = r1_grid[i] + 93,
                                     r_sig = rew_sig,
                                     beta_Sig = diag(1, nrow = 2),
                                     init_score = 93,
                                     projection_grid = optim_grid)
  acquisition_grid$PI[i,] = acquire_pi(post_pred_mu = inferred_surface$post_pred_mu,
                                            post_pred_sig = inferred_surface$post_pred_sig, 
                                            r1 = r1_grid[i] + 93, 
                                            xi_val = xi_val, 
                                            init_score = 93)
  acquisition_grid$EI[i,] = acquire_ei(post_pred_mu = inferred_surface$post_pred_mu,
                                            post_pred_sig = inferred_surface$post_pred_sig, 
                                            r1 = r1_grid[i] + 93, 
                                            xi_val = xi_val, 
                                            init_score = 93)
  acquisition_grid$UCB[i,] = acquire_ucb(post_pred_mu = inferred_surface$post_pred_mu,
                                              post_pred_sig = inferred_surface$post_pred_sig,
                                              quant = quant_val)
}


acquisition_grid_tilde = acquisition_grid

tau = pi/4
for(i in 1:3){
  neg_col_inds = which((theta_indexes > pi - tau))
  pos_col_inds = which((theta_indexes < tau))
  acquisition_grid_tilde[[i]][r1_grid < 0, neg_col_inds] = min(acquisition_grid[[i]])
  acquisition_grid_tilde[[i]][r1_grid >= 0, pos_col_inds] = min(acquisition_grid[[i]])
}
# GET MAXES ---------------------------------------------------------------

acquisition_maxes = list()
for(i in 1:3) {
  acquisition_maxes[[i]] = cbind.data.frame(
    r1 = r1_grid,
    angle = theta_indexes[apply(acquisition_grid_tilde[[i]], 1, which.max)])
}

acquisition_maxes_orig = list()
for(i in 1:3) {
  acquisition_maxes_orig[[i]] = cbind.data.frame(
    r1 = r1_grid,
    angle = theta_indexes[apply(acquisition_grid[[i]], 1, which.max)])
}

acquisition_melt = list()
acquisition_melt_full = list()

for(i in 1:3) {  
  # Melt into data.frames
  acquisition_melt[[i]] = reshape2::melt(acquisition_grid_tilde[[i]]) %>%
    dplyr::left_join(dplyr::tibble(Var1 = 1:length(r1_grid),
                                   r1 = r1_grid),
                     by = "Var1") %>%
    dplyr::left_join(dplyr::tibble(Var2 = 1:length(theta_indexes),
                                   angle = theta_indexes),
                     by = "Var2") %>%
    dplyr::select(r1, angle, value)
  
  # Mirror image across x-axis
  acquisition_melt_full[[i]] = acquisition_melt[[i]] %>%
    dplyr::filter(angle != 0) %>%
    dplyr::mutate(angle = -angle) %>%
    rbind(acquisition_melt[[i]])
}


acquisition_plots = list()

plot_titles = c(bquote("Move 2 Boundary " ~ widetilde("PI")), 
                bquote("Move 2 Boundary " ~ widetilde("EI")),
                bquote("Move 2 Boundary " ~ widetilde("UCB")))
sub_titles = c(bquote(xi["PI"] ~ " = " ~ .(xi_val) ~ ", " ~ tau ~ " = " ~ pi/4),
               bquote(xi["EI"] ~ " = " ~ .(xi_val) ~ ", " ~ tau ~ " = " ~ pi/4),
               bquote(p ~ " = " ~ .(quant_val) ~ ", " ~ tau ~ " = " ~ pi/4))
legend_titles = c(bquote("Augmented " ~ widetilde("PI")), bquote("Augmented " ~ widetilde("EI")), bquote("Augmented " ~ widetilde("UCB")))

for(i in 1:3){
  acquisition_plots[[i]] = ggplot2::ggplot(data = acquisition_melt_full[[i]], 
                                           ggplot2::aes(x = r1, y = angle, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                  name = bquote(.(legend_titles[[i]]))) + 
    ggplot2::geom_line(data = acquisition_maxes[[i]],
                       color = "green3",
                       size = 1,
                       ggplot2::aes(fill = NULL)) + 
    ggplot2::geom_line(data = acquisition_maxes[[i]] %>%
                         dplyr::mutate(angle = -angle),
                       color = "green3",
                       size = 1,
                       ggplot2::aes(fill = NULL)) + 
    ggplot2::geom_line(data = acquisition_maxes_orig[[i]],
                       color = "green3",
                       lty = 2,
                       size = .75,
                       ggplot2::aes(fill = NULL)) + 
    ggplot2::geom_line(data = acquisition_maxes_orig[[i]] %>%
                         dplyr::mutate(angle = -angle),
                       color = "green3",
                       lty = 2,
                       size = .75,
                       ggplot2::aes(fill = NULL)) + 
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
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 8),
                   plot.margin = ggplot2::margin(1, .1, .1, 1, "cm")
    ) + 
    ggplot2::labs(title = bquote(.(plot_titles[[i]])),
                  subtitle = bquote(.(sub_titles[[i]]))) +
    ggplot2::theme(
                   panel.border = element_rect(
                     colour = "gray",
                     fill = NA,
                     size = 1
                   ),
                   plot.margin = margin(5.5, -2.5, 5.5, 5.5, "pt")) + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(angle = 0),
          legend.key.width=unit(1,"cm")) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    ggplot2::annotate("text", x = 0, y = 7*pi/8+.1, label = '}', size = 6) +
    ggplot2::annotate("text", x = 3, y = 7*pi/8, label = deparse(bquote(tau)), size = 4, parse = T) +
    ggplot2::annotate("text", x = -1, y = 1*pi/8+.1, label = '{', size = 6) +
    ggplot2::annotate("text", x = -4, y = 1*pi/8, label = deparse(bquote(tau)), size = 4, parse = T) +
    {if(i != 1 )  theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())}+
    {if(i != 2 )  theme(axis.title.x=element_text(color = gray(1,0)))} + 
    {if(i == 3 )  theme( plot.margin = margin(5.5, 7, 5.5, 5.5, "pt"))}
}

# Save pdf
# pdf(file = "./figures/section_5/aug_acquisition_surfaces_w_curves.pdf",
#     width = 8, height = 4.35)
# cowplot::plot_grid(acquisition_plots[[1]],
#                    acquisition_plots[[2]],
#                    acquisition_plots[[3]], rel_widths = c(1.175,.95,1),
#                    nrow = 1)
# dev.off()

