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
  mutate(state = row_number(),
         pop_weight = population / sum(population)) %>%
  select(state, population, pop_weight, area)

few_state_scenario <- state_pop %>%
  mutate(population_cases = case_when(
    area %in% c("California", "Texas", "Florida", "New York") ~ round(population * 0.015),
    TRUE ~ 0
  ))

ICC <- 0.01
prevalence <- 0.005

set.seed(200)
many_state_scenario <- state_pop %>%
  mutate(., population_cases = round(rbeta(nrow(.), (1 - ICC) / ICC * prevalence, (1 - ICC) / ICC * (1 - prevalence)) * population))


generate_state_sample <- function(scenario, total_people_sampled) {
  n_states <- nrow(scenario)

  sample(n_states, size = total_people_sampled, replace = T, prob = scenario[["population"]]) %>%
    enframe(name = NULL, value = "state") %>%
    count(state, name = "n_sampled") %>%
    left_join(scenario, by = "state") %>%
    mutate(sample_cases = rbinom(n = n_states, size = n_sampled, prob = population_cases / population)) %>%
    select(
      group = state,
      size = n_sampled,
      cases = sample_cases,
      group_weight = pop_weight
    )
}

generate_state_sample(many_state_scenario, 12000)
generate_state_sample(few_state_scenario, 12000)
