library(tidyverse)
simulate_weights <- function(n_weights, coef_var) {
  # mu <- 1 / n_weights
  # sigma2 <- (coef_var / n_weights)^2
  # alpha <- ((1 - mu) / sigma2 - 1 / mu) * mu^2
  # beta <- alpha * (1 / mu - 1)

  alpha <- 1 / coef_var^2 - 1 / (n_weights * coef_var^2) - 1 / n_weights
  beta <- alpha * (n_weights - 1)

  x <- rbeta(n_weights, alpha, beta)
  if(sum(is.na(x))) return(NA)
  x / sum(x)
}

# We want to evaluate along 3 axis:

# Overall prevalence should always be low (say 0.5%)
# Weights should span a range of relevant coefficients of variation (say 10%-200%)
# Group-prevalences
# Prevalence can be allocated differently concenrated in highly weighted groups, low weighted groups,
population_prevalence <- 0.005
n_weights <- 8000
set.seed(200)
tmp <-
  crossing(prop_groups_with_cases = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 1),
           kind_of_groups_with_cases = c("mostly high weight", "any weight", "mostly low weight"),
           coefficient_of_variation = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10)) %>%
  # slice(rep(1:n(), 10)) %>%
  mutate(w = map(coefficient_of_variation, ~simulate_weights(n_weights = n_weights, coef_var = .))) %>%
  mutate(groups_with_prevalence = pmap(list(prop_groups_with_cases = prop_groups_with_cases,
                                            kind_of_groups_with_cases = kind_of_groups_with_cases,
                                            w = w),
                                       function(prop_groups_with_cases, kind_of_groups_with_cases, w) {
                                         case_when(
                                           kind_of_groups_with_cases == "mostly high weight" ~ sample(1:n_weights, size = round(prop_groups_with_cases * n_weights), replace = F, prob = w),
                                           kind_of_groups_with_cases == "any weight" ~ sample(1:n_weights, size = round(prop_groups_with_cases * n_weights), replace = F),
                                           kind_of_groups_with_cases == "mostly low weight" ~ sample(1:n_weights, size = round(prop_groups_with_cases * n_weights), replace = F, prob = 1 / (w + .Machine$double.xmin))
                                         )}
  )) %>%
  mutate(group_prevalence = map2_dbl(w, groups_with_prevalence, ~population_prevalence / sum(.x[.y]))) %>%
  filter(group_prevalence < 1) %>%
  mutate(theta = map2(groups_with_prevalence, group_prevalence, ~{
    theta <- numeric(n_weights)
    theta[.x] <- .y
    theta}))


tmp %>% mutate(theta_w_cor = map2_dbl(theta, w, ~cor(.x, .y))) %>%
  select(-groups_with_prevalence, -group_prevalence, -theta) %>%
  arrange(desc(theta_w_cor))
