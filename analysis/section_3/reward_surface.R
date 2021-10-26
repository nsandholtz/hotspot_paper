# Visualizing the reward surface
library(ggforce)

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

# Global reward surface ---------------------------------------------------

# mesh over Reward Surface
my_mesh1 <- expand.grid(x = seq(-450, 450, length.out = 200),
                        y = seq(-450, 450, length.out = 200))
my_mesh1 <- cbind.data.frame(my_mesh1,
                             to_polar(my_mesh1$x,
                                      my_mesh1$y,
                                      origin_x = 0, 
                                      origin_y = 0))

# Calculate the reward over the grid
my_mesh1$rew <- calc_reward(loc_x = my_mesh1$x, 
                       loc_y = my_mesh1$y,
                       hotspot_x = example_session$target_x,
                       hotspot_y = example_session$target_y,
                       R_0 = example_session$init_score,
                       grad = example_session$reward_gradient)

# Get rid of points outside the task region:
my_resh1 = dplyr::filter(my_mesh1,dist <= 450)

# task region
circle1 <- data.frame(
  x = 0,
  y = 0,
  r = 450
)

# click region
circle2 <- data.frame(
  x = 0,
  y = 0,
  r = 20
)

# plot surface
global_surf = ggplot(my_resh1, aes(x, y)) + 
  geom_raster(aes(fill = rew),interpolate=F) + 
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                name = "Reward"
                                #breaks = c(-400, -200, 0, 200, 400, 600, 800)
                                ) +
  geom_circle(aes(x0=x, y0=y, r=r), data = circle1, inherit.aes = F) +
  geom_circle(aes(x0=x, y0=y, r=r), data = circle2, color = "green4", inherit.aes = F) +
  coord_equal() + 
  annotate("text", x = -175, y = 0, label = "Move 1\nboundary") + 
  #theme_minimal() +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(angle = 45),
        legend.key.width=unit(1,"cm")) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  ggtitle("Global Objective Function")



# Local reward surface ----------------------------------------------------

# Put a mesh over the coordinates of interst.
my_mesh2 <- expand.grid(x = seq(-20, 20,length.out = 200),
                        y = seq(-20, 20,length.out = 200))
my_mesh2 <- cbind.data.frame(my_mesh2,
                             to_polar(my_mesh2$x,
                                      my_mesh2$y,
                                      origin_x = 0, origin_y = 0))

# Calculate the  reward over the grid
my_mesh2$rew <- calc_reward(loc_x = my_mesh2$x, 
                            loc_y = my_mesh2$y,
                            hotspot_x = example_session$target_x,
                            hotspot_y = example_session$target_y,
                            R_0 = example_session$init_score,
                            grad = example_session$reward_gradient)

# Get rid of points outside the click region:
my_resh2 = dplyr::filter(my_mesh2,dist <= 20)

circle3 <- data.frame(
  x = 0,
  y = 0,
  r = 20
)

# Plot local surface
local_surf = ggplot(my_resh2, 
                    aes(x, y)) + 
  geom_raster(aes(fill = rew), 
              interpolate=F) + 
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                name = "Reward (rescaled) "
                                #breaks = seq(370,430, by = 10)
                                ) +
  geom_circle(aes(x0=x, 
                  y0=y,
                  r=r), 
              data = circle3, 
              color = "green4", 
              inherit.aes = F, 
              size = 2) +
  coord_equal() + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(angle = 45),
        legend.key.width=unit(1,"cm")) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  ggtitle("Move 1 Objective Function")


# Save plots
pdf(file = "./figures/section_3/reward_surface.pdf",
    width = 10*.8, height = 6*.8)
cowplot::plot_grid(global_surf, local_surf, nrow = 1)
dev.off()


