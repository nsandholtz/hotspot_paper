# Raw move 2 angles

# Source utilities
library(tidyverse)
source("./analysis/utils.R")

# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

# Remove "burn-in" for each subject
first_n = 20

ok_sessions = session_dat %>%
  dplyr::filter(round_index > first_n,
         dist_hotspot > round_max_dist) 


subjects_to_plot = c(21, 24, 28)

example_dat = event_dat %>%
  #dplyr::filter(session_id %in% ok_sessions$session_id) %>%
  dplyr::mutate(lag_delta_score = dplyr::lag(delta_score)) %>%
  dplyr::filter(alt_id %in% subjects_to_plot) %>%
  filter(move == 2) %>%
  dplyr::select(lag_delta_score, rad_relative, alt_id)

facet_names <- c(
  `21` = "Subject 21",
  `24` = "Subject 24",
  `28` = "Subject 28"
)

rotated_move_2 = ggplot2::ggplot(data = example_dat,
                                  ggplot2::aes(x = lag_delta_score,
                                               y = rad_relative)) +
  ggplot2::geom_point(size = .5) +
  ggplot2::scale_x_continuous(name = expression(paste(Delta, r[1])),
                              breaks = c(-30, -20, -10, 0, 10, 20, 30)) +
  ggplot2::scale_y_continuous(
    name = expression(theta[2]),
    breaks = c(-pi, -pi / 2, 0, pi / 2, pi),
    labels = c(
      expression(-pi),
      expression(-pi / 2),
      0,
      expression(pi / 2),
      expression(pi)
    )
  ) +
  facet_wrap( ~ alt_id, labeller = as_labeller(facet_names)) + 
  ggplot2::theme_minimal() +
  # FOR WNAR VERSION
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    panel.border = element_rect(
      colour = "gray",
      fill = NA,
      size = 1
    ),
    strip.text.x = element_text(size = 10)
    #panel.spacing = unit(2.5, "lines")
  ) 
rotated_move_2
  
pdf(file = "./figures/section_4/rotated_move_2_example.pdf",
      width = 8, height = 16/5)
  rotated_move_2
  dev.off()




