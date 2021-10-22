#### Libraries and Data ----

library(plotrix)
library(dplyr)
library(ggplot2)

source("./analysis/utils.R")

# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")


# Example round --------------------------------------------------------

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

circle_example <- data.frame(
  x = example_events$x_loc[1],
  y = example_events$y_loc[1],
  r = 450
)

example_round = ggplot2::ggplot(example_events, aes(x_loc, y_loc, label = for_text)) + 
  ggforce::geom_circle(aes(x0=x, y0=y, r=r), data = circle_example, color = "black", fill = gray(.75,.5), inherit.aes = F) +
  ggplot2::geom_point(size = .75) + 
  ggrepel::geom_label_repel(size = 8,
                            label.padding = .5,
                            fill = "white",
                            direction = "y", nudge_x = 220) +
  ggplot2::geom_point(aes(x=target_x, y=target_y), data = example_events[1,], color = "red", size = 1) +
  ggforce::geom_circle(aes(x0=target_x, y0=target_y, r = 15), data = example_events[1,], color = "red", inherit.aes = F) +
  ggforce::geom_circle(aes(x0=target_x, y0=target_y, r = 30), data = example_events[1,], color = "red", inherit.aes = F) +
  ggforce::geom_circle(aes(x0=target_x, y0=target_y, r = 45), data = example_events[1,], color = "red", inherit.aes = F) +
  ggrepel::geom_label_repel(aes(x=target_x, y=target_y, label = paste("Hotspot:",round(max_score))), data = example_events[1,],
                            color = "black", 
                            size = 8,
                            fill = "white",
                            nudge_y = 100,
                            nudge_x = -50,
                            label.padding = .5) +
  #ggplot2::geom_point(aes(x_loc, y_loc), data = example_events[1,], color = "blue", size = .75) + 
  coord_equal() + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
  #ggtitle(paste("Example Round\n Subject:", example_info$alt_id))
  #ggtitle(paste("Example Round"))
print(example_round)


pdf(file = "./figures/section_2/example_round.pdf",
    width = 8, height = 8)
example_round
dev.off()

