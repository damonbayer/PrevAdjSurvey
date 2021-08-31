args <- commandArgs(trailingOnly=TRUE)
sjob <- as.integer(args[1])
cat("sjob is ", sjob, "\n")

suppressMessages({
  library(tidyverse)
  library(asht)
  # library(parallel)
  # library(foreach)
  # library(doRNG)
  # library(doParallel)
  library(furrr)
  library(fs)
  library(parallelly)
})
source("code/WprevSeSp_SRS.R")

save_path <- path("//", "data", "bayerdm", "confidence_interval_imperfect_assay")
dir_create(save_path)

plan(multisession, workers = parallelly::availableCores())
cat("using", parallelly::availableCores(), "cores\n")
n_partitions <- 500
n_replications <- 10000
cat("session planned\n")


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
    prev = c(0.005, 0.05),
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

cat("experiments designed\n")

if (sjob == 1) {
  write_rds(experimental_design, path(save_path, "005_experimental_design.rds"))
}

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

  n_tested_for_prevalence <- n_groups * tests_per_group
  est_specificity <- rbinom(1, n_tested_for_specificity, specificity) / n_tested_for_specificity
  est_sensitivity <- rbinom(1, n_tested_for_sensitivity, sensitivity) / n_tested_for_sensitivity

  known_positive_counts <- integer(n_groups)
  known_positive_counts[groups_with_prev] <- rbinom(n = n_groups_with_prev, size = tests_per_group, prob = group_prev)
  known_negative_counts <- tests_per_group - known_positive_counts

  true_positive_counts <- rbinom(n = n_groups, size = known_positive_counts, prob = sensitivity)
  false_negative_counts <- known_positive_counts - true_positive_counts

  true_negative_counts <- rbinom(n = n_groups, size = known_negative_counts, prob = specificity)
  false_positive_counts <- known_negative_counts - true_negative_counts

  apparent_positive_counts <- true_positive_counts + false_positive_counts
  apparent_positive_counts <- 0
  apparent_prevalence <- sum(apparent_positive_counts / tests_per_group * weights)
  var_hat_p_hat <- sum((weights / tests_per_group)^2 * apparent_positive_counts)
  stdErrPrev <- sqrt(var_hat_p_hat)

  adjusted_prevalence <- prevAdj(AP = apparent_prevalence, sen = est_sensitivity, spec = est_specificity)
  # adjusted_positive_counts = apparent_positive_counts / apparent_prevalence * adjusted_prevalence
  # adjusted_var_hat_p_hat <- sum((weights / tests_per_group)^2 * adjusted_positive_counts)
  WprevSeSp_result <-
    WprevSeSp_original(AP = apparent_prevalence,
                       stdErrPrev = stdErrPrev,
                       nP = n_tested_for_prevalence,
                       Se = est_sensitivity,
                       nSe = n_tested_for_sensitivity,
                       Sp = specificity,
                       nSp = n_tested_for_specificity)
  w <- weights / tests_per_group
  WprevSeSp_gamma_result <- WprevSeSp_gamma(apparent_positive_counts = apparent_positive_counts,
                  w = w,
                  nP = n_tested_for_prevalence,
                  Se = est_sensitivity,
                  nSe = n_tested_for_sensitivity,
                  Sp = specificity,
                  nSp = n_tested_for_specificity)
 # Need to implement other method
  list(conf_int_WprevSeSp = as.vector(WprevSeSp_result$conf.int),
       conf_int_WprevSeSp_gamma_result = as.vector(WprevSeSp_gamma_result$conf.int))
}

plan(multisession, workers = parallelly::availableCores() / 2)
cat("starting results\n")

set.seed(200)
results <-
  experimental_design %>%
  filter(specificity == max(specificity),
         group_distribution == "high",
         prop_groups_with_prev == 0.75,
         n_groups == 50,
         prev == 0.005) %>%
  sample_n(5) %>%
  # filter(coef_var == min(coef_var)) %>%
  slice(rep(1:n(), 100)) %>%
  group_by(design) %>%
  mutate(replication = row_number()) %>%
  mutate(interval_results = pmap(.l = list(weights, n_groups, n_groups_with_prev, groups_with_prev, group_prev, tests_per_group, n_tested_for_sensitivity, n_tested_for_specificity, sensitivity, specificity),
                        .f = calculate_interval_results)) %>%
  select(design, replication, prev, interval_results) %>%
  unnest_wider(interval_results) %>%
  pivot_longer(cols = starts_with("conf_int_"),
               names_prefix = "conf_int_",
               names_to = "method",
               values_to = "conf_int") %>%
  unnest_wider(conf_int, names_sep = "_") %>%
  mutate(lower_error = conf_int_1 > prev,
         upper_error = prev > conf_int_2,
         covered = !(lower_error | upper_error)) %>%
  arrange(design, replication)

cat("results finished\n")

write_rds(x = results, file = path(save_path, str_c("results_raw", sprintf("%04d", sjob), "of", n_partitions, sep = "_"), ext = "rds"))

cat("summarizing results\n")
results_summary <-
  results %>%
  select(-starts_with("conf_int")) %>%
  group_by(design, method) %>%
  summarize(
    lower_error_freq = mean(lower_error),
    upper_error_freq = mean(upper_error),
    coverage = mean(covered),
    .groups = "drop")

write_rds(x = results_summary, file = path(save_path, str_c("results_summary", sprintf("%04d", sjob), "of", n_partitions, sep = "_"), ext = "rds"))
