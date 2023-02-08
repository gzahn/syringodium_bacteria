library(tidyverse)
# library(imgpalr)
# library(colorblindr) # https://github.com/clauswilke/colorblindr
# library(colorBlindness)
library(ggplot2)
# library(reshape2)
#library(imgpalr)
#library(colorblindr) # https://github.com/clauswilke/colorblindr

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

#pal.discrete <- 
  #c("#172543", "#d05f0a", "#6FC9A6", "#0E494E", "#219CAE", "#057BCC", 
          #"#793025", "#65A00D", "#2D4828", "#f2a102", "#ebc765")

#palette_plot(pal.discrete) %>% 
#  cvd_grid()

ggplot2::theme_set(theme_minimal() +
                     theme(axis.text = element_text(face='bold'),
                           legend.title = element_text(face='bold')))




m_palette <-c("#40B0A0", "#EC6508", "#7715A3", "#0464DB", "#029AA0", "#0D6963",
      "#296B05", "#B70B46", "#DD303C", "#90B711", "#F89607", "#DBB678")



# creating a plot to test
# mat <- matrix(1:81, nrow = 9, ncol = 9)

#mat1 <- melt(t(mat[9:1, ]))
#len <- length(y)
#mat1$v2 <- cut(mat1$value,
               breaks = seq(0,ceiling(81/len)*len, 
                            length.out = len+1))
#ht <- ggplot(mat1) + 
#  geom_tile(aes(x=Var1, y=Var2, fill=v2)) + 
#  scale_fill_manual(values=y) + 
#  theme_bw()
# check the plot by CVD simulator
