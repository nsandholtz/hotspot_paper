library(dplyr)
library(ggforce)
source("~/Dropbox/Luke_Research/human_acquisition/manuscript/section_1/utils.R")

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
  mutate(for_text = NA) %>%
  left_join((example_session %>% 
               ungroup() %>%
               select(session_id)), by = "session_id")

# Move 2 ------------------------------------------------------------------


# Set acquisition params ---------------
r_sig <- .01
xi_val = 1
quant_val = .95


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
                                   r_sig = r_sig,
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


# Move 3  -----------------------------------------------------------------

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

m2_star = rbind(pi_max[1, c("x", "y")],
                ei_max[1, c("x", "y")],
                ucb_max[which(ucb_max$y <= 0), c("x", "y")])
m2_star[1,] = c(40,0)
# Need to manually change to c(40, 0)

# GET HYPOTHETICAL MOVE 2

sub_titles = c(bquote(xi["PI"] ~ " = " ~ .(xi_val)),
               bquote(xi["EI"] ~ " = " ~ .(xi_val)),
               bquote(p ~ " = " ~ .(quant_val)))

m3_acq_plot = list() 

for(i in 1:3){
  
  m3_opt_xy = get_opt_xy(x_rot = m2_star$x[i], 
                         y_rot = m2_star$y[i],
                         hotspot_x_rot = example_events$target_x_rot[1],
                         hotspot_y_rot = example_events$target_y_rot[1])
  
  move_2_rew = calc_reward(loc_x = m2_star$x[i], 
                         loc_y = m2_star$y[i],
                         hotspot_x = example_events$target_x_rot[1],
                         hotspot_y = example_events$target_y_rot[1],
                         grad = example_session$reward_gradient,
                         R_0 = example_events$score[1])
  m3_mesh <- expand.grid(x = seq(m2_star$x[i] - 20, 
                                    m2_star$x[i] + 20,
                                    length.out = 200),
                            y = seq(m2_star$y[i] - 20, 
                                    m2_star$y[i] + 20,
                                    length.out = 200))
  m3_mesh_rad <- cbind.data.frame(m3_mesh,
                                  to_polar(m3_mesh$x,
                                           m3_mesh$y,
                                           origin_x = m2_star$x[i], 
                                           origin_y = m2_star$y[i]))
  m3_circ_mesh_rad = dplyr::filter(m3_mesh_rad,dist <= 20)
  
  m3_inferred_surface = infer_surface_n(X = as.matrix(rbind(c(20,0),m2_star[i,])),
                                        r = c(example_events$score[2], move_2_rew),
                                        r_sig = r_sig,
                                        beta_Sig = diag(1, nrow = 2),
                                        init_score = example_events$score[1],
                                        projection_grid = m3_circ_mesh_rad[,c("x", "y")])
  
  m3_circ_mesh_rad$pp_mu = m3_inferred_surface$post_pred_mu
  m3_circ_mesh_rad$pp_sig = m3_inferred_surface$post_pred_sig
  
  # Get M3 Acquisition surfaces -----
  
  if(i == 1){
    m3_circ_mesh_rad$acq_surf = acquire_pi(post_pred_mu = m3_circ_mesh_rad$pp_mu,
                                          post_pred_sig = m3_circ_mesh_rad$pp_sig, 
                                          r1 = c(example_events$score[2], move_2_rew), 
                                          xi_val = xi_val, 
                                          init_score = example_events$score[1])
  } else if(i == 2){
    m3_circ_mesh_rad$acq_surf = acquire_ei(post_pred_mu = m3_circ_mesh_rad$pp_mu,
                                          post_pred_sig = m3_circ_mesh_rad$pp_sig, 
                                          r1 = c(example_events$score[2], move_2_rew), 
                                          xi_val = xi_val, 
                                          init_score = example_events$score[1])
  } else {
    m3_circ_mesh_rad$acq_surf = acquire_ucb(post_pred_mu = m3_circ_mesh_rad$pp_mu,
                                            post_pred_sig = m3_circ_mesh_rad$pp_sig, 
                                            quant = quant_val)
  }

  m3_max = m3_circ_mesh_rad[which(m3_circ_mesh_rad$acq_surf == max(m3_circ_mesh_rad$acq_surf, na.rm = TRUE)),] %>%
    filter(x == max(x)) %>%
    mutate(ab_y = abs(y)) %>%
    filter(ab_y == min(ab_y))
  
  # M3 plots -----------------------------------------------------------------
  
  circle_4 <- data.frame(
    x = m2_star$x[i],
    y = m2_star$y[i],
    r = 20
  )
  
  m3_dat = example_events[1:2,] %>% 
    select(x_rot, y_rot, score) %>%
    rbind.data.frame(data.frame(x_rot = m2_star$x[i],
                                y_rot = m2_star$y[i],
                                score = move_2_rew))
  
  
  m3_acq_plot[[i]] = ggplot2::ggplot(m3_circ_mesh_rad, aes(x, y)) + 
    geom_circle(aes(x0=x, y0=y, r=r), 
                data = circle0, 
                color = "black", 
                linetype = "dashed",
                inherit.aes = F) +
    geom_circle(aes(x0=x, y0=y, r=r), 
                data = circle1, 
                color = "black", 
                linetype = "dashed",
                inherit.aes = F) +
    ggplot2::geom_raster(aes(fill = acq_surf), 
                         interpolate=F,
                         alpha = 9/10) + 
    ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                  name = ifelse(i == 1, "Prob. of Improvement",
                                                ifelse(i == 2, "Expected Improvement", "Upper Conf. Bound"))
    ) +
    geom_circle(aes(x0=x, y0=y, r=r), 
                data = circle_4, 
                color = "black", 
                inherit.aes = F) +
    ggplot2::geom_point(aes(x=x_rot, y=y_rot), 
                        data = m3_dat, 
                        color = "gray",
                        size = 2) +
    ggplot2::geom_label(aes(x=x_rot, y=y_rot),
                        data = m3_dat,
                        hjust = 0,
                        label = c(deparse(bquote(r[.(example_events$move[1])] ~ " = " ~ .(round(m3_dat$score[1])))),
                                  deparse(bquote(r[.(example_events$move[2])] ~ " = " ~ .(round(m3_dat$score[2])))),
                                  deparse(bquote(r[.(example_events$move[3])] ~ " = " ~ .(round(m3_dat$score[3]))))),
                        parse = T,
                        nudge_x = 2,
                        nudge_y = 5,
                        color = "black",
                        label.padding = unit(.15, "lines"),
                        size = 3) +
    ggplot2::geom_point(aes(x=x, y=y), 
                        data = m3_max, 
                        color = "green2",
                        shape = "*",
                        size = 10) +
    ggplot2::geom_point(aes(x=x, y=y), 
                        data = data.frame(x = m3_opt_xy[1],
                                          y = m3_opt_xy[2]), 
                        color = "magenta",
                        size = 2) +
    xlim(-20,60) +
    ylim(-40,40) +
    coord_equal() + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(angle = 45),
          legend.key.width=unit(1,"cm")) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    labs(title= ifelse(i == 1,"Move 3 PI Surface",
                       ifelse(i == 2, "Move 3 EI Surface", 
                              "Move 3 UCB Surface")),
         subtitle = bquote(.(sub_titles[[i]])))
}

cowplot::plot_grid(m3_acq_plot[[1]],
                   m3_acq_plot[[2]],
                   m3_acq_plot[[3]], nrow = 1)


pdf(file = "./figures/section_3/example_acquisition_m3.pdf",
    width = 10, height = 5)
cowplot::plot_grid(
                   m3_acq_plot[[1]], m3_acq_plot[[2]], m3_acq_plot[[3]], 
                   nrow = 1)
dev.off()
