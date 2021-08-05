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
  crossing(prop_groups_with_cases = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75),
           kind_of_groups_with_cases = c("mostly high weight", "any weight", "mostly low weight"),
           coefficient_of_variation = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10)) %>%
  bind_rows(crossing(prop_groups_with_cases = 1,
                     kind_of_groups_with_cases = "any weight",
                     coefficient_of_variation = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10))) %>%
  mutate(design = 1:n()) %>%
  slice(rep(1:n(), each = 10)) %>%
  group_by(design) %>%
  mutate(replication = row_number()) %>%
  ungroup() %>%
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
    theta})) %>%
  mutate(mean_weight = map_dbl(w, mean),
         var_weight = map_dbl(w, var)) %>%
  mutate(coef_var = sqrt(var_weight) / mean_weight,
         theta_w_cor = map2_dbl(theta, w, ~cor(.x, .y)))



crossing(prop_groups_with_cases = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75),
         kind_of_groups_with_cases = c("mostly high weight", "any weight", "mostly low weight"),
         coefficient_of_variation = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10)) %>%
  bind_rows(crossing(prop_groups_with_cases = 1,
                     kind_of_groups_with_cases = "any weight",
                     coefficient_of_variation = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10))) %>%
  mutate(design = 1:n()) %>%
  filter(design %in% which(!(1:171 %in% tmp$design)))




tmp %>%
  count(design) %>%
  arrange(n)


tmp %>%
  ggplot(aes(coef_var, theta_w_cor, color = kind_of_groups_with_cases)) +
  geom_point() +
  cowplot::theme_cowplot()

# You can't have a high coefficient of variation and maintain low theat_w correlation
# or even be able to generate datasets with infection propotion < 1.
# at least with 8000 weights.
tmp %>%
  count(kind_of_groups_with_cases, sign(theta_w_cor), coefficient_of_variation) %>%
  View()

tmp %>%
  filter(is.na(theta_w_cor)) %>%
  pull(prop_groups_with_cases)

tmp %>%
  mutate(mean_weight = map2_dbl(groups_with_prevalence, w, ~(mean(.y[.x]))),
         var_weight = map2_dbl(groups_with_prevalence, w, ~(var(.y[.x])))) %>%
  arrange(mean_weight) %>%
  mutate(rank = dense_rank(mean_weight)) %>%
  group_by(kind_of_groups_with_cases) %>%
  summarize(median(rank))
  tail()
  pull(kind_of_groups_with_cases) %>%
  cat(sep = "\n")

tmp

tmp %>% mutate(theta_w_cor = map2_dbl(theta, w, ~cor(.x, .y))) %>%
  select(-groups_with_prevalence, -group_prevalence, -theta) %>%
  arrange(desc(theta_w_cor))
