library(dplyr)
library(ggforce)

# Source utilities
source("./analysis/utils.R")
source("./analysis/constants.R")

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

# Calculate optimal angle

m2_opt_xy = get_opt_xy(x_rot = 20, 
                       y_rot = 0, 
                       hotspot_x_rot = to_cart(example_events$target_dist[1], example_events$target_rot_rad[1])$x,
                       hotspot_y_rot = to_cart(example_events$target_dist[1], example_events$target_rot_rad[1])$y)

# Calculate acquisition surfaces -----

# Set acquisition params 
xi_val = 1
quant_val = .95

circ_mesh_rad$pi_surf = acquire_pi(post_pred_mu = circ_mesh_rad$pp_mu,
                                   post_pred_sig = circ_mesh_rad$pp_sig, 
                                   r1 = example_events$score[2], 
                                   xi_val = xi_val, 
                                   init_score = example_events$score[1])
circ_mesh_rad$ei_surf = acquire_ei(post_pred_mu = circ_mesh_rad$pp_mu,
                                   post_pred_sig = circ_mesh_rad$pp_sig, 
                                   r1 = example_events$score[2], 
                                   xi_val = xi_val, 
                                   init_score = example_events$score[1])
circ_mesh_rad$ucb_surf = acquire_ucb(post_pred_mu = circ_mesh_rad$pp_mu,
                                     post_pred_sig = circ_mesh_rad$pp_sig, 
                                     quant = quant_val
)

# Get Maximums -----

pi_max = circ_mesh_rad[which(circ_mesh_rad$pi_surf == max(circ_mesh_rad$pi_surf, na.rm = TRUE)),] %>%
  filter(x == max(x)) %>%
  mutate(ab_y = abs(y)) %>%
  filter(ab_y == min(ab_y))
ei_max = circ_mesh_rad[which(circ_mesh_rad$ei_surf == max(circ_mesh_rad$ei_surf, na.rm = TRUE)),]
ucb_max = circ_mesh_rad[which(circ_mesh_rad$ucb_surf == max(circ_mesh_rad$ucb_surf, na.rm = TRUE)),]


# PI plot -----------------------------------------------------------------

circle0 <- data.frame(
  x = 0,
  y = 0,
  r = 20
)

circle1 <- data.frame(
  x = 20,
  y = 0,
  r = 20
)


local_pi = ggplot2::ggplot(circ_mesh_rad, aes(x, y)) + 
  geom_circle(aes(x0=x, y0=y, r=r), 
              data = circle0, 
              color = "black", 
              linetype = "dashed",
              inherit.aes = F) +
  ggplot2::geom_raster(aes(fill = pi_surf), 
                       interpolate=F,
                       alpha = 9/10) + 
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                name = "Prob. of Improvement") +
  geom_circle(aes(x0=x, y0=y, r=r), data = circle1, color = "black", inherit.aes = F) +
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
  ggplot2::geom_point(aes(x=x, y=y), 
                      data = pi_max, 
                      color = "green2",
                      shape = "*",
                      size = 10) +
  coord_equal() + 
  theme_minimal() +
  ggplot2::theme(panel.border = element_rect(
                   colour = "gray",
                   fill = NA,
                   size = 1
                 ),
                 plot.margin = margin(5.5, -2.5, 5.5, 5.5, "pt")) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(angle = 45),
        legend.key.width=unit(1,"cm")) +
  theme(axis.title.x=element_text(color = gray(1,0))) + 
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  labs(title="Move 2 PI Surface",
       subtitle=bquote(xi["PI"] ~ " = " ~ .(xi_val)))


# EI plot -----------------------------------------------------------------
  
local_ei = ggplot2::ggplot(circ_mesh_rad, aes(x, y)) + 
  geom_circle(aes(x0=x, y0=y, r=r), 
              data = circle0, 
              color = "black", 
              linetype = "dashed",
              inherit.aes = F) +
  ggplot2::geom_raster(aes(fill = ei_surf), 
                       interpolate=F,
                       alpha = 9/10) + 
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                name = "Expected Improvement"
  ) +
  geom_circle(aes(x0=x, y0=y, r=r), data = circle1, color = "black", inherit.aes = F) +
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
  ggplot2::geom_point(aes(x=x, y=y), 
                      data = ei_max, 
                      color = "green2",
                      shape = "*",
                      size = 10) +
  ggplot2::geom_point(aes(x=x, y=-y), 
                      data = ei_max, 
                      color = "green2",
                      shape = "*",
                      size = 10) +
  coord_equal() + 
  theme_minimal() +
  ggplot2::theme(panel.border = element_rect(
    colour = "gray",
    fill = NA,
    size = 1
  ),
  plot.margin = margin(5.5, -2.5, 12.5, 5.5, "pt")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(angle = 45),
        legend.key.width=unit(1,"cm"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  labs(title="Move 2 EI Surface",
       subtitle=bquote(xi["EI"] ~ " = " ~ .(xi_val)))


# UCB_plot ----------------------------------------------------------------


local_ucb = ggplot2::ggplot(circ_mesh_rad, aes(x, y)) + 
  geom_circle(aes(x0=x, y0=y, r=r), 
              data = circle0, 
              color = "black", 
              linetype = "dashed",
              inherit.aes = F) +
  ggplot2::geom_raster(aes(fill = ucb_surf), 
                       interpolate=F,
                       alpha = 9/10) + 
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                name = "Upper Conf. Bound"
  ) +
  geom_circle(aes(x0=x, y0=y, r=r), data = circle1, color = "black", inherit.aes = F) +
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
  ggplot2::geom_point(aes(x=x, y=y), 
                      data = ucb_max, 
                      color = "green2",
                      shape = "*",
                      size = 10) +
  coord_equal() + 
  theme_minimal() +
  ggplot2::theme(legend.position="bottom",
                 panel.border = element_rect(
                   colour = "gray",
                   fill = NA,
                   size = 1
                 )) +
 theme(axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())+
theme(axis.title.x=element_text(color = gray(1,0))) +
theme( plot.margin = margin(5.5, 7, 5.5, 5.5, "pt")) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(angle = 45),
        legend.key.width=unit(1,"cm")) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  labs(title="Move 2 UCB Surface",
       subtitle=bquote(p ~ " = " ~ .(quant_val)))

# Plot for manuscript ----

# pdf(file = "./figures/section_3/example_acquisition_m2.pdf",
#     width = 8, height = 18/5)
# cowplot::plot_grid(local_pi, local_ei, local_ucb, rel_widths = c(1.175,.92,1), nrow = 1)
# dev.off()
