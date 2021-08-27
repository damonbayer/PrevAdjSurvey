library(tidyverse)
library(furrr)
source("code/WprevSeSp_SRS.R")
n_partitions <- 500
# Simulate Weights
simulate_weights <- function(n_groups, coef_var) {
  # mu <- 1 / n_groups
  # sigma2 <- (coef_var / n_groups)^2
  # alpha <- ((1 - mu) / sigma2 - 1 / mu) * mu^2
  # beta <- alpha * (1 / mu - 1)

  alpha <- 1 / coef_var^2 - 1 / (n_groups * coef_var^2) - 1 / n_groups
  beta <- alpha * (n_groups - 1)

  x <- rbeta(n_groups, alpha, beta)
  if(sum(is.na(x))) return(NA)
  x / sum(x)
}

set.seed(200)
weight_candidates <-
  crossing(n_groups = c(50, 8000),
           # target_coef_var = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 100, 1000)) %>%
           target_coef_var = seq(0.01, 5, length.out = 500)) %>%
  mutate(design = 1:n()) %>%
  # slice(rep(1:n(), each = 100)) %>%
  group_by(design) %>%
  mutate(replication = row_number()) %>%
  ungroup() %>%
  mutate(weights = future_map2(.options = furrr_options(seed = 200),
                               .x = n_groups,
                               .y = target_coef_var,
                               .f = ~simulate_weights(n_groups = .x, coef_var = .y)))
cat("weights simulated\n")


experimental_design <-
  weight_candidates %>%
  mutate(obs_coef_var = map_dbl(weights, ~(sd(.) / mean(.)))) %>%
  mutate(coef_var_rel_error = abs(obs_coef_var - target_coef_var) / target_coef_var) %>%
  group_by(design) %>%
  filter(coef_var_rel_error == min(coef_var_rel_error)) %>%
  select(design, n_groups, coef_var = obs_coef_var, weights) %>%
  mutate(weights = map(weights, ~sort(., decreasing = T))) %>%
  mutate(tests_per_group = case_when(
    n_groups == 50 ~ 200,
    n_groups == 8000 ~ 1
  )) %>%
  crossing(
    n_tested_for_sensitivity = 60,
    n_tested_for_specificity = 300,
    prevalence = c(0.005, 0.05),
    sensitivity = 0.95,
    specificity = seq(0.75, 1, by = 0.05),
    prop_groups_with_prev = c(0.02, 0.05, 0.10, 0.25, 0.5, 0.75),
    group_distribution = c("high", "uniform", "low")) %>%
  mutate(n_groups_with_prev = round(prop_groups_with_prev * n_groups)) %>%
  filter(n_groups_with_prev > 0) %>%
  mutate(groups_with_prev = pmap(
    list(n_groups = n_groups,
         n_groups_with_prev = n_groups_with_prev,
         group_distribution = group_distribution),
    function(n_groups, n_groups_with_prev, group_distribution) {
      case_when(
        n_groups == n_groups_with_prev ~ list(1:n_groups),
        n_groups == n_groups_with_prev + 1 ~ list((1:n_groups)[-round(n_groups / 2)]),
        group_distribution == "high" ~ list(head(1:n_groups, n_groups_with_prev)),
        group_distribution == "uniform" ~ list(as.integer(round(seq(1, n_groups, length.out = n_groups_with_prev + 2)[c(-1, -(n_groups_with_prev + 2))]))),
        group_distribution == "low" ~ list(tail(1:n_groups, n_groups_with_prev))) %>%
        unlist()
    })) %>%
  mutate(group_prev = pmap_dbl(list(prev = prev,
                                    weights = weights,
                                    groups_with_prev = groups_with_prev),
                               function(prev, weights, groups_with_prev)
                                 prev / sum(weights[groups_with_prev]))) %>%
  filter(group_prev <= 1) %>% # impossible to achieve desired population prevalence when only a small number of groups have cases
  mutate(design = 1:n()) %>%
  mutate(., partition = ceiling(row_number() / (nrow(.) / n_partitions)))


calculate_interval_results <- function(weights,
                                       n_groups,
                                       n_groups_with_prev,
                                       groups_with_prev,
                                       group_prev,
                                       tests_per_group,
                                       n_tested_for_sensitivity,
                                       n_tested_for_specificity,
                                       sensitivity,
                                       specificity) {
  # Working on this
  true_test_results <- integer(n_groups)
  positive_tests[groups_with_prev] <- rbinom(n = n_groups_with_prev, size = tests_per_group, prob = group_prev)
  sample_positiviity <- positive_tests / tests_per_group

  # Prepare for methods
  p_hat <- sum(sample_positiviity * weights)
  var_hat_p_hat <- sum((weights / tests_per_group)^2 * positive_tests)
  x <- positive_tests
  w <- weights / tests_per_group

  # Apply each method to simulated data
  result_wspoissonTest <- wspoissonTest(x = x, w = w, midp = F)
  result_wspoissonTest_midp <- wspoissonTest(x = x, w = w, midp = T)
  result_AC <- AC_method(p_hat = p_hat, var_hat_p_hat = var_hat_p_hat, n_psu = n_groups, n_ssu = tests_per_group, adjusted = F)
  result_CP <- CP_method(p_hat = p_hat, var_hat_p_hat = var_hat_p_hat, n_psu = n_groups, n_ssu = tests_per_group, adjusted = F)

  list(conf_int_wspoissonTest = as.vector(result_wspoissonTest$conf.int),
       conf_int_wspoissonTest_midp = as.vector(result_wspoissonTest_midp$conf.int),
       conf_int_AC = as.vector(result_AC$conf.int),
       conf_int_CP = as.vector(result_CP$conf.int))
}

# assume no relationship between spec and sens and weight
tmp <- experimental_design %>% slice(1)

weights <- tmp$weights[[1]]


weights <- tmp$weights[[1]]
n_groups <- tmp$n_groups[[1]]
n_groups_with_prev <- tmp$n_groups_with_prev[[1]]
groups_with_prev <- tmp$groups_with_prev[[1]]
group_prev <- tmp$group_prev[[1]]
tests_per_group <- tmp$tests_per_group[[1]]
n_tested_for_sensitivity <- tmp$n_tested_for_sensitivity[[1]]
n_tested_for_specificity <- tmp$n_tested_for_specificity[[1]]
sensitivity <- tmp$sensitivity[[1]]
specificity <- tmp$specificity[[1]]
