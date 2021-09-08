library(tidyverse)
library(fs)

load("C:/Users/bayerdm/Downloads/dataForDamon.RData")

dir_create("data")
data.for.Damon %>%
  as_tibble() %>%
  mutate(y = as.integer(y),
         hispanic = hispanic == "Yes") %>%
  write_csv("data/kalish_data.csv")
