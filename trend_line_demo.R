

# code to generate a trend line to add to a population trajectory plot

library(tidyverse)

# # my locally saved version of the annual indices associated with population goals
# all_inds <- readRDS("data/all_socb_goal_indices.rds")
#
# # selecting just American Tree Sparrow
# demo_inds <- all_inds %>%
#   filter(resultsCode == "CBC",
#          speciesCode == "amtspa")
# # saving just these indices for later sharing
# write_csv(demo_inds,"data/amtspa_goal_indices.csv")

#loading American Tree Sparrow CBC indices associated with goal
demo_inds <- read_csv("data/amtspa_goal_indices.csv")

#extracting the trend value associated with these indices
trend <- demo_inds %>%
  select(trnd) %>%
  distinct() %>% # removing repeated values
  mutate(trnd = log((trnd/100)+1)) %>%   # transforming %/year trend to log-scale change
  unlist()


## function to generate a trend line that represents the estimated trend
## calculated on the log-transformed indices and then
## transformed back to the original index scale
trend_line <- function(indices,#for this function indices is a vector of index values
                       years, #for this function year is a vector of year values
                       slope = trend){ # for this function slope is a single real value

  # identify the middle year of the available time-series
  # used to center the year values and ensure that the mid-point of the line
  # falls at the center of the time-series
  mid_year <- round(median(years))

  # log-transforming the annual indices ("index")
  l_ind <- log(indices)

  # vertical centering of the trend-line as the mean of the log-transformed indices
  y_center <- mean(l_ind) #intercept in equation below

  # generating the line using y = mx + b formula plus exponential re-transformation
  line = exp(slope*(years-mid_year) + y_center) #

  return(line)

}
# applying the above function to the index values.
demo_inds <- demo_inds %>%
  mutate(vis_line = trend_line(index,
                               year,
                               trend))


## setting an upper limit to the demo plot below, cutting off the extreme upper limits
y_upper_limit <- max(demo_inds$index)*1.1

plot_demo <- ggplot(data = demo_inds,
                    aes(x = year,y = index))+
  geom_point()+
  geom_errorbar(aes(ymin = indexLowerCI,
                    ymax = indexUpperCI),
                alpha = 0.3,
                width = 0)+
  geom_line(aes(y = vis_line))+
  coord_cartesian(ylim = c(0,y_upper_limit))+ # zooms in on the points,
  theme_bw()

pdf("demo_trend_line_plot.pdf",
    height = 5)
plot_demo
dev.off()




