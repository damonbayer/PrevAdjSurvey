suppressMessages({
  library(tidyverse)
  library(asht)
  library(parallel)
  library(foreach)
  library(doRNG)
  library(doParallel)
  library(fs)
  library(doMPI)
  library(future)
})

cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

n_replications <- 1

timings <- as.POSIXct(rep(NA, 9))

timings[1] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")
# Agresti-Coull Method
AC_method <- function(p_hat,
                      var_hat_p_hat,
                      adjusted = F,
                      n_psu = NULL,
                      n_ssu = NULL,
                      n_strata = 1,
                      conf.level = 0.95) {
  if (adjusted & is.null(n_psu)) stop("Must specify n_psu if using adjusted method.")

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

timings[2] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")

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


timings[3] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")

weight_candidates <-
  crossing(n_groups = c(50, 8000),
           # target_coef_var = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 100, 1000)) %>%
           target_coef_var = seq(0.01, 5, length.out = 500)) %>%
  mutate(design = 1:n()) %>%
  # slice(rep(1:n(), each = 100)) %>%
  group_by(design) %>%
  mutate(replication = row_number()) %>%
  ungroup()

set.seed(200)
weight_candidates$weights <-
  foreach(n_groups = weight_candidates$n_groups,
          target_coef_var = weight_candidates$target_coef_var) %dorng% {
            simulate_weights(n_groups, target_coef_var)
            }

timings[4] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")

experimental_design <-
  weight_candidates %>%
  mutate(obs_coef_var = map_dbl(weights, ~(sd(.) / mean(.)))) %>%
  mutate(coef_var_rel_error = abs(obs_coef_var - target_coef_var) / target_coef_var) %>%
  group_by(design) %>%
  filter(coef_var_rel_error == min(coef_var_rel_error)) %>%
  select(design, n_groups, coef_var = obs_coef_var, weights) %>%
  mutate(weights = map(weights, sort)) %>%
  mutate(tests_per_group = case_when(
    n_groups == 50 ~ 200,
    n_groups == 8000 ~ 1
  )) %>%
  crossing(
    prev = 0.005,
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
  mutate(design = 1:n())

timings[5] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")

results <-
  experimental_design %>%
  slice(rep(1:n(), each = n_replications)) %>%
  group_by(design) %>%
  mutate(replication = row_number()) %>%
  ungroup() %>%
  select(design, replication, everything())

write_rds(experimental_design, "//data/bayerdm/fixed_weights_experimental_design.rds")
rm(experimental_design)

timings[6] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")

results$interval_results <-
  foreach(weights = results$weights,
          n_groups = results$n_groups,
          n_groups_with_prev = results$n_groups_with_prev,
          groups_with_prev = results$groups_with_prev,
          group_prev = results$group_prev,
          n_groups = results$n_groups,
          tests_per_group = results$tests_per_group,
          .packages = "asht") %dorng% {
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

             list(result_wspoissonTest = result_wspoissonTest,
                  result_wspoissonTest_midp = result_wspoissonTest_midp,
                  result_AC = result_AC,
                  result_CP = result_CP)
           }

timings[7] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")

results <-
  results %>%
  select(design, replication, prev, interval_results)

timings[8] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")

write_rds(results, "//data/bayerdm/fixed_weights_results_raw_1.rds")

results_summary <-
  results %>%
  unnest_wider(interval_results) %>%
  pivot_longer(cols = starts_with("result"),
               names_prefix = "result_",
               names_to = "method",
               values_to = "htest") %>%
  mutate(conf.int = htest %>%
           map("conf.int")) %>%
  unnest_wider(conf.int, names_sep = "_") %>%
  mutate(lower_error = conf.int_1 > prev,
         upper_error = prev > conf.int_2,
         covered = !(lower_error | upper_error)) %>%
  select(-htest, -starts_with("conf.int")) %>%
  group_by(design, method) %>%
  summarize(
    lower_error_freq = mean(lower_error),
    upper_error_freq = mean(upper_error),
    coverage = mean(covered),
    .groups = "drop")

timings[9] <- Sys.time(); write_rds(timings, "//data/bayerdm/timings1.rds")
write_rds(results_summary, "//data/bayerdm/fixed_weights_results_summary_1.rds")

closeCluster(cl)
mpi.quit()
