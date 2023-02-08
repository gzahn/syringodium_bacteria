library(tidyverse)
# library(imgpalr)
# library(colorblindr) # https://github.com/clauswilke/colorblindr

# pal.discrete <- 
#   imgpalr::image_pal(file="./media/palletepics.png",
#                    type = "qual",
#                    n = 12, 
#                    plot = TRUE,
#                    saturation = c(.5,1),
#                    brightness = c(0,1),
#                    seed = 123)

pal.discrete <- 
  c("#41132D","#944F1C","#13600E","#2A5C12","#10411D","#5A5F13",
    "#057BCC","#793025","#219CAE","#172543","#7B681E","#1B554C")

# palette_plot(pal.discrete) %>% 
#   cvd_grid()

ggplot2::theme_set(theme_minimal() +
                     theme(axis.text = element_text(face='bold'),
                           legend.title = element_text(face='bold')))





