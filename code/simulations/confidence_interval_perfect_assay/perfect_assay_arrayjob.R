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

save_path <- path("//", "data", "bayerdm", "confidence_interval_perfect_assay")
dir_create(save_path)

plan(multisession, workers = parallelly::availableCores())
cat("using", parallelly::availableCores(), "cores\n")
# n_partitions <- 100
n_partitions <- 250
n_replications <- 10000
cat("session planned\n")

# Agresti-Coull Method
AC_method <- function(p_hat,
                      var_hat_p_hat,
                      adjusted = F,
                      n_psu = NULL,
                      n_ssu = NULL,
                      n_strata = 1,
                      conf.level = 0.95) {
  if (adjusted & is.null(n_psu)) stop("Must specify n_psu if using adjusted method.")
  if (var_hat_p_hat == 0 & (is.null(n_psu) | is.null(n_ssu))) stop("Must specify n_psu and n_ssu if var_hat_p_hat is 0.")

  alpha <- 1 - conf.level
  z_val <- qnorm(1 - alpha / 2)
  c <- z_val^2 / 2

  if (adjusted) {
    design_df <- n_psu - n_strata
    t_val <- qt(1 - alpha / 2, design_df)
  }

  if (var_hat_p_hat > 0) {
    n_eff <- p_hat * (1 - p_hat) / var_hat_p_hat * ifelse(adjusted, (z_val / t_val)^2, 1)
  } else {
    n_eff <- n_psu * n_ssu
  }

  x_tilde <- p_hat * n_eff + c
  n_tilde <- n_eff + 2 * c
  p_tilde <- x_tilde / n_tilde

  ci <- p_tilde + c(-1, 1) * z_val * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  attr(ci, "conf.level") <- conf.level

  rval <- list(conf.int = ci, estimate = p_tilde, p_hat = p_hat, var_hat_p_hat = var_hat_p_hat)
  class(rval) <- "htest"
  rval
}

# Clopper-Pearson Method
CP_method <- function(p_hat,
                      var_hat_p_hat,
                      adjusted = F,
                      n_psu = NULL,
                      n_ssu = NULL,
                      n_strata = 1,
                      conf.level = 0.95) {

  if (var_hat_p_hat == 0 & (is.null(n_psu) | is.null(n_ssu))) stop("Must specify n_psu and n_ssu if var_hat_p_hat is 0.")
  alpha <- 1 - conf.level

  if (adjusted) {
    z_val <- qnorm(1 - alpha / 2)
    design_df <- n_psu - n_strata
    t_val <- qt(1 - alpha / 2, design_df)
  }

  if (var_hat_p_hat > 0) {
    n_eff <- p_hat * (1 - p_hat) / var_hat_p_hat * ifelse(adjusted, (z_val / t_val)^2, 1)
  } else {
    n_eff <- n_psu * n_ssu
  }

  n <- n_eff
  x <- p_hat * n_eff

  # nu_1 <- 2 * x
  # nu_2 <- 2 * (n - x + 1)
  # nu_3 <- 2 * (x + 1)
  # nu_4 <- 2 * (n - x)
  #
  # ci_l <- (nu_1 * qf(alpha / 2, nu_1, nu_2)) / (nu_2 + nu_1 * qf(alpha / 2, nu_1, nu_2))
  # ci_u <- (nu_3 * qf(1 - alpha / 2, nu_3, nu_4)) / (nu_4 + nu_3 * qf(1 - alpha / 2, nu_3, nu_4))

  ci_l <- qbeta(alpha / 2, x, n - x + 1)
  ci_u <- qbeta(1 - alpha / 2, x + 1, n - x)

  ci <- c(ci_l, ci_u)
  attr(ci, "conf.level") <- conf.level

  rval <- list(conf.int = ci, p_hat = p_hat, var_hat_p_hat = var_hat_p_hat)
  class(rval) <- "htest"
  rval
}

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
    prev = c(0.005, .05),
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
  write_rds(experimental_design, path(save_path, "experimental_design.rds"))
}

calculate_interval_results <- function(weights,
                                       n_groups,
                                       n_groups_with_prev,
                                       groups_with_prev,
                                       group_prev,
                                       tests_per_group) {
  positive_tests <- integer(n_groups)
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
plan(multisession, workers = parallelly::availableCores() / 2)
cat("starting results\n")
results <-
  experimental_design %>%
  filter(partition == sjob) %>% # this is where I would select the partition
  slice(rep(1:n(), each = n_replications)) %>%
  group_by(design) %>%
  mutate(replication = row_number()) %>%
  ungroup() %>%
  select(design, replication, everything()) %>%
  mutate(interval_results = future_pmap(.options = furrr_options(seed = 200),
                                        list(weights = weights,
                                             n_groups = n_groups,
                                             n_groups_with_prev = n_groups_with_prev,
                                             groups_with_prev = groups_with_prev,
                                             group_prev = group_prev,
                                             tests_per_group = tests_per_group),
                                        function(weights,
                                                 n_groups,
                                                 n_groups_with_prev,
                                                 groups_with_prev,
                                                 group_prev,
                                                 tests_per_group)
                                          calculate_interval_results(weights,
                                                                     n_groups,
                                                                     n_groups_with_prev,
                                                                     groups_with_prev,
                                                                     group_prev,
                                                                     tests_per_group))) %>%
  select(design, replication, prev, interval_results) %>%
  unnest_wider(interval_results) %>%
  pivot_longer(cols = starts_with("conf_int_"),
               names_prefix = "conf_int_",
               names_to = "method",
               values_to = "conf_int") %>%
  unnest_wider(conf_int, names_sep = "_") %>%
  mutate(lower_error = conf_int_1 > prev,
         upper_error = prev > conf_int_2,
         covered = !(lower_error | upper_error))
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
