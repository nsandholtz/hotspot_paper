library(dplyr)
library(ggplot2)

# Source utilities
source("./analysis/utils.R")


# LOAD DATA ---------------------------------------------------------------

event_dat <- readRDS("./data/event_dat.rds")
session_dat <- readRDS("./data/session_dat.rds")

# Remove "burn-in" for each subject
first_n = 20

ok_sessions = session_dat %>%
  dplyr::filter(round_index > first_n,
                target_dist > session_max_dist)
ok_event_dat = event_dat %>%
  dplyr::filter(session_id %in% ok_sessions$session_id) 


# Do subjects on average get move 3 correct? -----------------------------

ok_event_dat = ok_event_dat %>%
  mutate(pos_rot_rad = ifelse(rot_rad < 0, 
                              rot_rad + 2*pi,
                              rot_rad),
         pos_target_rot_rad = ifelse(target_rot_rad < 0, 
                                     target_rot_rad + 2*pi,
                                     target_rot_rad),
         opt_dist = minuspi_to_pi(pos_rot_rad - pos_target_rot_rad))

my_cust_lab = expression(pi, 7*pi/8, 3*pi/4, 5 * pi/8, 
             pi / 2, 3 * pi / 8, pi / 4,  pi / 8, 0, - pi / 8, 
           -pi / 4, -3*pi / 8, -pi / 2, -5*pi / 8, -3 * pi / 4, -7 * 
             pi/8, -pi)

my_mids = (round(seq(-pi, pi, by = pi/16), digits = 2) + lead(round(seq(-pi, pi, by = pi/16), digits = 2)))/2 

for_gg = data.frame(move = integer(),
                    ind = integer(),
                    counts = integer(),
                    dens = double())

for(i in 1:3) {
  my_hist = ok_event_dat %>%
    filter(move == i) %>%
    with(hist(
      opt_dist,
      breaks = seq(-pi, pi, by = pi / 16),
      prob = F
    ))
  
  for_gg = rbind(
    for_gg,
    data.frame(
      move = i,
      ind = my_mids[1:length(my_hist$counts)],
      counts = my_hist$counts,
      dens = my_hist$density
    )
  )
}
  
facet_names <- c(
  `1` = "Move 1 angle",
  `2` = "Move 2 angle",
  `3` = "Move 3 angle"
)

polar_hist = ggplot(for_gg, aes(x=ind, y=dens)) +
  geom_bar(stat='identity') + theme_light() +
  theme(axis.title.y=element_text(angle=0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank()) + 
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1)) + 
  coord_polar(start = 3*pi/2, direction = 1) + 
  scale_x_continuous("", limits = c(-pi, pi),
                     breaks = seq(-pi, pi, by = pi/8),
                     labels = my_cust_lab) +
  scale_y_continuous(expression(paste("P(", theta, ")")),
                     limits=c(0,.4), breaks = seq(0,1,by = .1)) + 
  geom_hline(yintercept=seq(0,.4,by = .1), 
             color = "lightgray", size=.15) + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.border = element_blank()) + 
  facet_wrap(~ move, labeller = as_labeller(facet_names))  + 
  theme(
    strip.text.x = element_text(
      size = 12, color = "black"
    ),
    strip.background = element_rect(
      color="black", fill=gray(.9), linetype="solid"
    )
  )


# pdf(file = "./figures/section_3/angle_to_hotspot.pdf",
#     width = 10, height = 4)
# polar_hist
# dev.off()

