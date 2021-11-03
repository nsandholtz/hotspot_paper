library(tidyverse)

# Source utilities
source("./analysis/utils.R")

# Set global parameters
beta_Sig = diag(1, nrow = 2)
rew_sig = .01
xi_vals = c(0,10,20,30)
quant_vals = c(.5, .65, .85, .99)

r1 = seq(-33,33, by = .5)
angle_resolution = 180*2
optim_grid = to_cart(r = 20, 
                     theta = seq(0, pi, length.out = angle_resolution)) %>%
  mutate(x = x + 20)



# Calculate Acquisition Curves

pi_curves = matrix(NA, length(r1), length(xi_vals))
ei_curves = matrix(NA, length(r1), length(xi_vals))
ucb_curves = matrix(NA, length(r1), length(quant_vals))
for(i in 1:length(xi_vals)){
  pi_curves[,i] = acquire_pi_curve(
    r1_vals = r1,
    r_sig = rew_sig,
    beta_Sig = beta_Sig,
    init_score = 0,
    xi_val = xi_vals[i],
    projection_grid = optim_grid
  )$opt_angle
  
  ei_curves[,i] = acquire_ei_curve(
    r1_vals = r1,
    r_sig = rew_sig,
    beta_Sig = beta_Sig,
    init_score = 0,
    xi_val = xi_vals[i],
    projection_grid = optim_grid
  )$opt_angle
  
  ucb_curves[,i] = acquire_ucb_curve(
    r1_vals = r1,
    r_sig = rew_sig,
    beta_Sig = beta_Sig,
    init_score = 0,
    quant_val = quant_vals[i],
    projection_grid = optim_grid
  )$opt_angle
}


all_curves = list(pi_curves,
                  ei_curves,
                  ucb_curves)

# PLOT ----

opt_curve_plots = list()
my_palette = list(colorRampPalette(c("black", "green")),
                  colorRampPalette(c("black", "blue")),
                  colorRampPalette(c("black", "red")))
plot_titles = c("Candidate PI Curves",
                "Candidate EI Curves",
                "Candidate UCB Curves")
legend_titles = c(bquote(xi[PI] ~ " : "),
                  bquote(xi[EI] ~ " : "),
                  bquote(p ~ " : "))
legend_names = list(xi_vals,
                    xi_vals,
                    quant_vals)


for(i in 1:3) {
  if(i != 3){
    par_vals = xi_vals
  } else {
    par_vals = quant_vals
  }
  plot_dat = cbind.data.frame(delta_r = rep(r1, length(par_vals)),
                              Exploration = rep(par_vals, each = length(r1)),
                              angle = c(all_curves[[i]]))
  opt_curve_plots[[i]] = ggplot2::ggplot(data = plot_dat, ggplot2::aes(x = delta_r,
                                                                       y = angle,
                                                                       group = Exploration)) +
    ggplot2::geom_line(ggplot2::aes(color = as.factor(Exploration))) +
    ggplot2::geom_line(data = plot_dat %>%
                         dplyr::mutate(angle = -angle),
                       ggplot2::aes(color = as.factor(Exploration))) +
    ggplot2::scale_color_manual(values = my_palette[[i]](length(par_vals)),
                                name = bquote(.(legend_titles[[i]])),
                                labels = legend_names[[i]]) +
    ggplot2::scale_x_continuous(name = expression(paste(Delta, r[1])),
                                breaks = c(-30,-20,-10, 0, 10, 20, 30)) +
    ggplot2::scale_y_continuous(
      name = expression(theta[2]),
      breaks = c(-pi,-pi / 2, 0, pi / 2, pi),
      labels = c(
        expression(-pi),
        expression(-pi / 2),
        0,
        expression(pi / 2),
        expression(pi)
      )
    ) +
    ggplot2::ggtitle(plot_titles[i]) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10)) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.size = unit(.1, 'cm'),
      plot.margin = ggplot2::margin(1, .1, .1, 1, "cm")
    ) + 
    ggplot2::theme(legend.position="bottom",
                   panel.border = element_rect(
                     colour = "gray",
                     fill = NA,
                     size = 1
                   ),
                   plot.margin = margin(5.5, -2.5, 5.5, 5.5, "pt")) + 
    {if(i != 1 )  theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())}+
    {if(i != 2 )  theme(axis.title.x=element_text(color = gray(1,0)))} + 
    {if(i == 3 )  theme( plot.margin = margin(5.5, 7, 5.5, 5.5, "pt"))}+
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 5)))
}

pdf(file = "./figures/section_4/sample_acquisition_curves.pdf",
    width = 8, height = 18/5)
cowplot::plot_grid(opt_curve_plots[[1]], 
                   opt_curve_plots[[2]],
                   opt_curve_plots[[3]],rel_widths = c(1.175,.95,1),
                   nrow = 1)  
dev.off()

# reshape
# sample_curves_df = rbind(pivot_longer(data.frame(pi_curves, 
#                                                  acq_type = "PI",
#                                                  delta_rew = r1), 
#                                       cols = c("X1", "X2", "X3", "X4"), 
#                                       names_to = "par_val", 
#                                       values_to = "acq_val") %>%
#                            arrange(par_val),
#                          pivot_longer(data.frame(ei_curves, 
#                                                  acq_type = "EI",
#                                                  delta_rew = r1), 
#                                       cols = c("X1", "X2", "X3", "X4"),
#                                       names_to = "par_val", 
#                                       values_to = "acq_val") %>%
#                            arrange(par_val),
#                          pivot_longer(data.frame(ucb_curves, 
#                                                  acq_type = "UCB",
#                                                  delta_rew = r1), 
#                                       cols = c("X1", "X2", "X3", "X4"),
#                                       names_to = "par_val", 
#                                       values_to = "acq_val") %>%
#                            arrange(par_val))







