library(tidyverse)
# raw_results <- read_csv("code/simulations/confidence_interval_imperfect_assay/combined_imperfect_assay_raw.csv")
# experimental_design <- readRDS("code/simulations/confidence_interval_imperfect_assay/experimental_design.rds")
#
# raw_results %>%
#   filter(conf_int_2 < conf_int_1)
#
# target_design <-
#   raw_results %>%
#   filter(method != "wspoissonTest") %>%
#   arrange(desc(lower_error_freq)) %>%
#   pull(design) %>%
#   unique() %>%
#   pluck(2)
#
# strange_design <- experimental_design %>% filter(design == target_design) %>% as.list()
#
# strange_design$weights <- strange_design$weights[[1]]
# strange_design$groups_with_prev <- strange_design$groups_with_prev[[1]]
# write_rds(strange_design, "code/strange_design.rds")

strange_design <- read_rds("code/strange_design.rds")
source("code/simulations/confidence_interval_imperfect_assay/WprevSeSp_SRS.R")
prevAdj <- function(AP, sen, spec) (AP + spec - 1) / (sen + spec - 1)


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
                       Sp = est_specificity,
                       nSp = n_tested_for_specificity)
  w <- weights / tests_per_group
  WprevSeSp_gamma_result <-
    WprevSeSp_gamma(apparent_positive_counts = apparent_positive_counts,
                    w = w,
                    nP = n_tested_for_prevalence,
                    Se = est_sensitivity,
                    nSe = n_tested_for_sensitivity,
                    Sp = est_specificity,
                    nSp = n_tested_for_specificity)

  # What if we treat it as perfect?
  result_wspoissonTest <- wspoissonTest(x = apparent_positive_counts, w = w, midp = F)

  list(conf_int_WprevSeSp = as.vector(WprevSeSp_result$conf.int),
       conf_int_WprevSeSp_gamma_result = as.vector(WprevSeSp_gamma_result$conf.int),
       conf_int_wspoissonTest = as.vector(result_wspoissonTest$conf.int))
}


set.seed(1000)
# Get intervals once
calculate_interval_results(
  weights = strange_design$weights, n_groups = strange_design$n_groups,
  n_groups_with_prev = strange_design$n_groups_with_prev,
  groups_with_prev = strange_design$groups_with_prev,
  group_prev = strange_design$group_prev,
  tests_per_group = strange_design$tests_per_group,
  n_tested_for_sensitivity = strange_design$n_tested_for_sensitivity,
  n_tested_for_specificity = strange_design$n_tested_for_specificity,
  sensitivity = strange_design$sensitivity,
  specificity = strange_design$specificity
)

n_reps <- 1000
set.seed(200)
results_original_g <-
  tibble(
    prev = strange_design$prev,
    interval_results = map(
      1:n_reps,
      ~ calculate_interval_results(
        weights = strange_design$weights, n_groups = strange_design$n_groups,
        n_groups_with_prev = strange_design$n_groups_with_prev,
        groups_with_prev = strange_design$groups_with_prev,
        group_prev = strange_design$group_prev,
        tests_per_group = strange_design$tests_per_group,
        n_tested_for_sensitivity = strange_design$n_tested_for_sensitivity,
        n_tested_for_specificity = strange_design$n_tested_for_specificity,
        sensitivity = strange_design$sensitivity,
        specificity = strange_design$specificity
      )
    )
  ) %>%
  unnest_wider(interval_results) %>%
  pivot_longer(
    cols = starts_with("conf_int_"),
    names_prefix = "conf_int_",
    names_to = "method",
    values_to = "conf_int"
  ) %>%
  unnest_wider(conf_int, names_sep = "_") %>%
  mutate(
    lower_error = conf_int_1 > prev,
    upper_error = prev > conf_int_2,
    covered = !(lower_error | upper_error)
  )

results_original_g %>%
  group_by(method) %>%
  summarize(mean(lower_error),
            mean(upper_error),
            mean(covered))

prevAdj <- function(AP, sen, spec) {
  # spec = 1 - phi_n
  phi_n <- 1 - spec
  phi_p <- sen
  theta_1 <- AP

  # if ((phi_n < phi_p) & (phi_p < theta_1)) {
  #   return(1)
  # } else if ((phi_p >= theta_1) & (theta_1 >= phi_n)) {
  #   return((theta_1 - phi_n) / (phi_p - phi_n))
  # } else {
  #   return(0)
  # }

  ((phi_n < phi_p) & (phi_p < theta_1)) + # (If this is TRUE, add 1)
    ((phi_p >= theta_1) & (theta_1 >= phi_n)) * ((theta_1 - phi_n) / (phi_p - phi_n)) # If this is true, run the normal adjustment
  # else, just returns 0
}

set.seed(200)
results_modified_g <-
  tibble(
    prev = strange_design$prev,
    interval_results = map(
      1:200,
      ~ calculate_interval_results(
        weights = strange_design$weights, n_groups = strange_design$n_groups,
        n_groups_with_prev = strange_design$n_groups_with_prev,
        groups_with_prev = strange_design$groups_with_prev,
        group_prev = strange_design$group_prev,
        tests_per_group = strange_design$tests_per_group,
        n_tested_for_sensitivity = strange_design$n_tested_for_sensitivity,
        n_tested_for_specificity = strange_design$n_tested_for_specificity,
        sensitivity = strange_design$sensitivity,
        specificity = strange_design$specificity
      )
    )
  ) %>%
  unnest_wider(interval_results) %>%
  pivot_longer(
    cols = starts_with("conf_int_"),
    names_prefix = "conf_int_",
    names_to = "method",
    values_to = "conf_int"
  ) %>%
  unnest_wider(conf_int, names_sep = "_") %>%
  mutate(
    lower_error = conf_int_1 > prev,
    upper_error = prev > conf_int_2,
    covered = !(lower_error | upper_error)
  )

results_modified_g %>%
  group_by(method) %>%
  summarize(mean(lower_error),
            mean(upper_error),
            mean(covered))
