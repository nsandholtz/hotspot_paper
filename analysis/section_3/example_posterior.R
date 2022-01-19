# Posterior predictive surfaces for experiment 2

library(dplyr)
library(ggplot2)
library(ggforce)

# Source utilities
source("./analysis/utils.R")
source("./analysis/constants.R")

# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

r1 <- seq(-30,30, length.out = 300)
# Example round -----------------------------------------------------------

subject_dat = session_dat %>% 
  group_by(alt_id) %>% summarise(num_sessions = n())
tens = which(session_dat$num_moves == 10)

example_session = session_dat[tens[4],]
example_events = event_dat %>%
  filter(session_id == example_session$session_id) %>%
  left_join((example_session %>% 
               ungroup() %>%
               select(session_id)), by = "session_id")

# Rotate data
example_events$x_rot = to_cart(example_events$dist, example_events$rot_rad)$x
example_events$y_rot = to_cart(example_events$dist, example_events$rot_rad)$y
example_events$x_rot[1] = 0
example_events$y_rot[1] = 0
example_events$target_x_rot = to_cart(example_events$target_dist, example_events$target_rot_rad)$x
example_events$target_y_rot = to_cart(example_events$target_dist, example_events$target_rot_rad)$y

my_mesh <- expand.grid(x = seq(0, 40,length.out = 200),
                       y = seq(-20, 20,length.out = 200))

my_mesh_rad <- cbind.data.frame(my_mesh,
                                to_polar(my_mesh$x,
                                         my_mesh$y,
                                         origin_x = 20, 
                                         origin_y = 0))

circ_mesh_rad = dplyr::filter(my_mesh_rad,dist <= 20)

# Infer the reward surface using the bayesian linear model
inferred_surface = infer_surface_1(x1 = 20, y1 = 0,
                                   r1 = example_events$score[2],
                                   r_sig = rew_sig,
                                   beta_Sig = diag(1, nrow = 2),
                                   init_score = example_events$score[1],
                                   projection_grid = circ_mesh_rad[,c("x", "y")])

circ_mesh_rad$pp_mu = inferred_surface$post_pred_mu
circ_mesh_rad$pp_sig = inferred_surface$post_pred_sig


# Plot predictive dist ----------------------------------------------------

circle3 <- data.frame(
  x = 20,
  y = 0,
  r = 20
)

local_surf = ggplot2::ggplot(circ_mesh_rad, aes(x, y)) + 
  ggplot2::geom_raster(aes(fill = pp_mu), interpolate=F) + 
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                name = "Reward"
                                #breaks = seq(370,430, by = 10)
                                ) +
  geom_circle(aes(x0=x, y0=y, r=r), data = circle3, color = "black", inherit.aes = F) +
  ggplot2::geom_point(aes(x=x_rot, y=y_rot), 
                      data = example_events[1:2,], 
                      color = "gray",
                      size = 2) +
  ggplot2::geom_label(aes(x=x_rot, y=y_rot),
                      data = example_events[1:2,],
                      hjust = 0,
                      label = c(deparse(bquote(r[.(example_events$move[1])] ~ " = " ~ .(round(example_events$score[1])))),
                                deparse(bquote(r[.(example_events$move[2])] ~ " = " ~ .(round(example_events$score[2]))))),
                      parse = T,
                      nudge_x = 2,
                      nudge_y = -5,
                      color = "black",
                      label.padding = unit(.15, "lines"),
                      size = 3) +
  coord_equal() + 
  theme_minimal() +  ggplot2::theme(panel.border = element_rect(
    colour = "gray",
    fill = NA,
    size = 1
  ),
  plot.margin = margin(5.5, -2.5, 5.5, 5.5, "pt")) + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45),
        legend.key.width=unit(1,"cm")) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Posterior Predictive Mean")


local_sd = ggplot2::ggplot(circ_mesh_rad, aes(x, y)) + 
  ggplot2::geom_raster(aes(fill = pp_sig), interpolate=F) + 
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                name = "Reward"
                                #breaks = seq(370,430, by = 10)
  ) +
  geom_circle(aes(x0=x, y0=y, r=r), data = circle3, color = "black", inherit.aes = F) +
  ggplot2::geom_point(aes(x=x_rot, y=y_rot), 
                      data = example_events[1:2,], 
                      color = "gray",
                      size = 2) +
  ggplot2::geom_label(aes(x=x_rot, y=y_rot),
                      data = example_events[1:2,],
                      hjust = 0,
                      label = c(deparse(bquote(r[.(example_events$move[1])] ~ " = " ~ .(round(example_events$score[1])))),
                                deparse(bquote(r[.(example_events$move[2])] ~ " = " ~ .(round(example_events$score[2]))))),
                      parse = T,
                      nudge_x = 2,
                      nudge_y = -5,
                      color = "black",
                      label.padding = unit(.15, "lines"),
                      size = 3) +
  coord_equal() + 
  theme_minimal() +  
  ggplot2::theme(panel.border = element_rect(
    colour = "gray",
    fill = NA,
    size = 1
  ),
  plot.margin = margin(5.5, 7, 14.5, 5.5, "pt")) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45),
        legend.key.width=unit(1,"cm")) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Posterior Predictive SD")

pdf(file = "./figures/section_3/example_posterior.pdf",
    width = 5.25, height = 4)
cowplot::plot_grid(local_surf, local_sd, rel_widths = c(1.175,1), nrow = 1)
dev.off()
