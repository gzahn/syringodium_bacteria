library(tidyverse)
library(imgpalr)
library(colorblindr)

# Creating a custom color palette from an image
# pal.discrete <-
#   image_pal(file="./media/palletepics.png",
#           type = 'qual',
#           n = 12,
#           plot = TRUE,
#           saturation = c(.5, 1),
#           brightness = c(0, 1),
#           seed = 123)

pal.discrete <- 
  c("#13600E", "#2A5C12", "#10411D", "#219CAE", "#5A5F13", "#057BCC", "#944F1C", "#41132D", "#793025","#172543", "#7B681E", "#1B554C")

palette_plot(pal.discrete,label_size = 0) %>% 
   cvd_grid()

pal.discrete2 <- 
  c("#359B73", # 1
    "#F0E442", # 2
    "#FFB2FD", # 3
    "#FF7F52", # 4
    "#064273", # 5
    "#000000", # 6
    "#D55E00", # 7
    "#9F0162", # 8
    "#FFD3CD", # 9
    "#6A0213", # 10
    "#00AF8E", # 11
    "#004002"
    )

palette_plot(pal.discrete2, label_size = 0) %>% 
  cvd_grid()

theme_set(theme_minimal() +
            theme(axis.text = element_text(face = 'bold'),
                  legend.title = element_text(face = 'bold')))

