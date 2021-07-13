library(tidyverse)

state_pop <- tibble(
  area = c(
    "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
    "Connecticut", "Delaware", "District of Columbia", "Florida",
    "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa",
    "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts",
    "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana",
    "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico",
    "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma",
    "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
    "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington",
    "West Virginia", "Wisconsin", "Wyoming"
  ),
  population = c(
    5024279, 733391, 7151502, 3011524, 39538223, 5773714, 3605944,
    989948, 689545, 21538187, 10711908, 1455271, 1839106, 12812508,
    6785528, 3190369, 2937880, 4505836, 4657757, 1362359, 6177224,
    7029917, 10077331, 5706494, 2961279, 6154913, 1084225, 1961504,
    3104614, 1377529, 9288994, 2117522, 20201249, 10439388, 779094,
    11799448, 3959353, 4237256, 13002700, 1097379, 5118425, 886667,
    6910840, 29145505, 3271616, 643077, 8631393, 7705281, 1793716,
    5893718, 576851
  )
) %>%
  mutate(state = row_number()) %>%
  select(state, population, area)

few_state_scenario <- state_pop %>%
  mutate(population_cases = case_when(
    area %in% c("California", "Texas", "Florida", "New York") ~ round(population * 0.015),
    # area %in% c("California", "Texas", "Florida", "New York") ~ round(population * 0.005),
    TRUE ~ 0
  ))


ICC <- 0.01
prevalence <- 0.005

set.seed(200)
many_state_scenario <- state_pop %>%
  mutate(., population_cases = round(rbeta(nrow(.), (1 - ICC) / ICC * prevalence, (1 - ICC) / ICC * (1 - prevalence)) * population))

scenario <- many_state_scenario



generate_state_sample <- function(scenario, m, n_tilde) {
  scenario %>%
    sample_n(size = m, replace = T, weight = population) %>%
    mutate(., sample_cases = rbinom(n = nrow(.), size = n_tilde, prob = population_cases / population)) %>%
    mutate(sample_size = n_tilde) %>%
    select(state, population, sample_cases, sample_size)
}




# group_by(group = state, population) %>%
#   summarize(cases = sum(sample_cases),
#             size = n_tilde * n(),
#             .groups = "drop") %>%
#   mutate(group_weight = population / sum(population)) %>%
#   select(-population)


# generate_state_sample(many_state_scenario, m = 12000, n_tilde = 1)
# generate_state_sample(many_state_scenario, 12000)
# generate_state_sample(few_state_scenario, 12000)
