# Header ------------------------------------------------------------------

library(ggforce)
library(gridExtra)
library(grid)
library(dplyr)

# Source utility functions
source("./analysis/utils.R")

# Move 1 diagram ----

scores = data.frame(
  x = c(0,0),
  y = c(0, 20),
  for_text = c("Move 0: 93","Move 1: 98") 
)

m1_bound <- data.frame(
  x0 = 0,
  y0 = 0,
  r = 20
)

origin = data.frame(
  x = 0,
  y = 0
)

arrow_segs <- data.frame(x1 = c(0,0,0,0,0,0,0,0) , 
                         x2 = c(12,12,-12,-12,12*sqrt(2),0,-12*sqrt(2),0), 
                         y1 = c(0,0,0,0,0,0,0,0), 
                         y2 = c(12,-12,12,-12,0,12*sqrt(2),0,-12*sqrt(2)))

move_1 = ggplot() +
  geom_circle(aes(x0 = x0, y0 = y0, r = r),
              data = m1_bound,
              color = NA,
              fill = gray(.75,.5)) +
  ggtitle("Move 1") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(.3,"cm"), 
                             type = "closed"),
               size = 1,
               linetype = 2,
               color = "green2",
               data = arrow_segs) + 
  geom_point(aes(x = x, y = y), 
             data = origin,
             color = "black",
             size = 16) + 
  annotate(geom="text", 
           x=0, y=0, 
           label="0",
           color="white",
           size = 10
  ) + 
ggplot2::geom_label(aes(x=x, y=y, 
                        label = for_text), 
                    data = scores[1,],
                    hjust = .5, nudge_y = -5,
                    color = "black", 
                    size = 4) +
  coord_fixed()
move_1

# Move 2 diagram ----

m1 = data.frame(
  x = 0,
  y = 20
)

m2_bound <- data.frame(
  x0 = 0,
  y0 = 20,
  r = 20
)

two_prime = data.frame(
  x = c(-20, 0),
  y = c(20, 40)
)

two_prime_half = data.frame(
  x = c(20/sqrt(2)),
  y = c(20+20/sqrt(2))
)

m2_arrow_segs1 <- data.frame(x1 = c(0,0), 
                            x2 = c(-16,12), 
                            y1 = c(20,20), 
                            y2 = c(20,20))
m2_arrow_segs2 <- data.frame(x1 = c(0,20/sqrt(2)), 
                            x2 = c(0,20/sqrt(2)), 
                            y1 = c(20, 20), 
                            y2 = c(36, 30))

move_2 = ggplot() +
  geom_circle(aes(x0 = x0, y0 = y0, r = r),
              data = m1_bound,
              color = NA,
              fill = gray(.75,.5)) +
  geom_circle(aes(x0 = x0, y0 = y0, r = r),
              data = m2_bound,
              color = NA,
              fill = gray(.75,.5)) +
  ggtitle("Move 2") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(.3,"cm"), 
                             type = "closed"),
               size = 1,
               linetype = 2,
               color = "green2",
               data = m2_arrow_segs1) + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(.3,"cm"), 
                             type = "closed"),
               size = 1,
               linetype = 2,
               color = "blue",
               data = m2_arrow_segs2) + 
  geom_point(aes(x = x, y = y), 
             data = origin,
             color = "black",
             size = 8) + 
  geom_point(aes(x = x, y = y), 
             data = m1,
             color = "black",
             size = 8) +
  geom_point(aes(x = x, y = y), 
             data = two_prime,
             color = c("green2","blue"),
             size = 8) +
  geom_point(aes(x, y), 
             data = two_prime_half,
             shape="\u25D6", 
             colour="blue", 
             size=15) +
  geom_point(aes(x, y), 
             data = two_prime_half,
             shape="\u25D7", 
             colour="green2", 
             size=15) +
  annotate(geom="text", 
           x=0, y=0, 
           label="0",
           color="white",
           size = 5
  ) + 
  annotate(geom="text", 
           x=0, y=20, 
           label="1",
           color="white",
           size = 5
  ) +    
  annotate(geom="text", 
           x=-20, y=20, 
           label="2a",
           color="white",
           size = 4
  ) + 
  annotate(geom="text", 
           x=0, y=40, 
           label="2b",
           color="white",
           size = 4
  ) +    
  annotate(geom="text", 
           x=20/sqrt(2), y=20/sqrt(2)+20, 
           label="2c",
           color="white",
           size = 4
  ) +    
  ggplot2::geom_label(aes(x=x, y=y, 
                          label = for_text), 
                      data = scores,
                      hjust = .5, nudge_y = -5,
                      color = "black", 
                      size = 3) +
  coord_fixed()
move_2

# Requires X quartz
# cairo_pdf("./figures/section_2/Move_1_diagram.pdf",height = 4, family="Arial Unicode MS")
# grid.arrange(move_1, nrow = 1)
# dev.off()
# 
# cairo_pdf("./figures/section_2/Move_2_diagram.pdf",height = 4, family="Arial Unicode MS")
# grid.arrange(move_2, nrow = 1)
# dev.off()

