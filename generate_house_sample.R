library(tidyverse)

set.seed(200)
house_pop <- tibble(house_size = c(sample(7, size = 408663, replace = T, prob = c(36198, 44742, 19337, 16262, 7446, 2919, 1546)), 2, 1)) %>%
  mutate(group = row_number())
n_houses <- nrow(house_pop)

ICC <- .1
prevalence <- 0.01

set.seed(200)
many_house_scenario <- house_pop %>%
  mutate(house_cases = round(rbeta(n_houses, (1 - ICC) / ICC * prevalence, (1 - ICC) / ICC * (1 - prevalence)) * house_size))

ICC <- .91
prevalence <- 0.01
set.seed(200)
few_house_scenario <- house_pop %>%
  mutate(house_cases = round(rbeta(n_houses, (1 - ICC) / ICC * prevalence, (1 - ICC) / ICC * (1 - prevalence)) * house_size))

generate_house_sample <- function(scenario, total_houses_sampled) {
  scenario %>%
    sample_n(total_houses_sampled, replace = F) %>%
    select(group, size = house_size, cases = house_cases)
}

generate_house_sample(many_house_scenario, 1000)
generate_house_sample(few_house_scenario, 1000)
