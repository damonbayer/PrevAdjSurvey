library(tidyverse)
library(extraDistr)
library(asht)
library(cowplot)
library(furrr)
library(fs)
library(scales)
library(glue)

plan(multisession, workers = 7)

AC_method <- function(p_hat,
                      var_hat_p_hat,
                      adjusted = F,
                      n_psu = NULL,
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

  n_eff <- p_hat * (1 - p_hat) / var_hat_p_hat * ifelse(adjusted, (z_val / t_val)^2, 1)

  x_tilde <- p_hat * n_eff + c
  n_tilde <- n_eff + 2 * c
  p_tilde <- x_tilde / n_tilde

  ci <- p_tilde + c(-1, 1) * z_val * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  attr(ci, "conf.level") <- conf.level

  rval <- list(conf.int = ci, estimate = p_tilde, p_hat = p_hat, var_hat_p_hat = var_hat_p_hat)
  class(rval) <- "htest"
  return(rval)
}

CP_method <- function(p_hat,
                      var_hat_p_hat,
                      adjusted = F,
                      n_psu = NULL,
                      n_strata = 1,
                      conf.level = 0.95) {
  alpha <- 1 - conf.level

  if (adjusted) {
    z_val <- qnorm(1 - alpha / 2)
    design_df <- n_psu - n_strata
    t_val <- qt(1 - alpha / 2, design_df)
  }

  n_eff <- p_hat * (1 - p_hat) / var_hat_p_hat * ifelse(adjusted, (z_val / t_val)^2, 1)

  n <- n_eff <- p_hat * (1 - p_hat) / var_hat_p_hat
  x <- p_hat * n_eff

  nu_1 <- 2 * x
  nu_2 <- 2 * (n - x + 1)
  nu_3 <- 2 * (x + 1)
  nu_4 <- 2 * (n - x)

  ci_l <- (nu_1 * qf(alpha / 2, nu_1, nu_2)) / (nu_2 + nu_1 * qf(alpha / 2, nu_1, nu_2))
  ci_u <- (nu_3 * qf(1 - alpha / 2, nu_3, nu_4)) / (nu_4 + nu_3 * qf(1 - alpha / 2, nu_3, nu_4))

  ci <- c(ci_l, ci_u)
  attr(ci, "conf.level") <- conf.level

  rval <- list(conf.int = ci, p_hat = p_hat, var_hat_p_hat = var_hat_p_hat)
  class(rval) <- "htest"
  return(rval)
}

# n_replications <- 10000

simulation_design_8000_1 <-
  crossing(prevalence = 0.005,
           ICC = c(0.1),
           dirichlet_alpha = 10^(seq(-2, 2, 0.5)),
           n_weights = 8000,
           samples_per_weight = 1) %>%
  mutate(design = 1:n()) %>%
  select(design, everything())


simulation_design_200_50 <-
  crossing(prevalence = 0.005,
           ICC = c(0.0001, 0.025, 0.05, 0.1),
           dirichlet_alpha = 10^(seq(-2, 2, 0.5)),
           n_weights = 200,
           samples_per_weight = 50) %>%
  mutate(design = 1:n()) %>%
  select(design, everything())

set.seed(200)

generate_results <- function(simulation_design, n_replications) {
  experimental_results <-
    simulation_design %>%
    slice(rep(1:n(), each = n_replications)) %>%
    group_by(design) %>%
    mutate(replication = row_number()) %>%
    ungroup() %>%
    mutate(interval_results =
             future_pmap(.options = furrr_options(seed = 200),
                         list(prevalence, ICC, dirichlet_alpha, n_weights, samples_per_weight),
                         function(prevalence, ICC, dirichlet_alpha, n_weights, samples_per_weight) {
                           alpha <- (1 - ICC) / ICC * prevalence
                           beta <- (1 - ICC) / ICC * (1 - prevalence)

                           dat <- tibble(
                             weight = as.vector(rdirichlet(1, alpha = rep(dirichlet_alpha, n_weights))),
                             true_positivity = rbeta(n_weights, alpha, beta),
                             n_tested = samples_per_weight) %>%
                             mutate(sample_positives = rbinom(n = n(), size = n_tested, prob = true_positivity)) %>%
                             mutate(sample_positivity = sample_positives / n_tested)

                           mean_weight <- mean(dat$weight)
                           var_weight <- var(dat$weight)
                           coef_var <- sqrt(var_weight) / mean_weight

                           p_hat <- sum(dat$sample_positivity * dat$weight)
                           var_hat_p_hat <- sum((dat$weight / dat$n_tested)^2 * dat$sample_positives) # this is incorrect

                           result_wspoissonTest <-
                             wspoissonTest(x = dat$sample_positives,
                                           w = dat$weight / dat$n_tested,
                                           midp = F)

                           result_wspoissonTest_midp <-
                             wspoissonTest(x = dat$sample_positives,
                                           w = dat$weight / dat$n_tested,
                                           midp = T)

                           result_AC <- AC_method(p_hat = p_hat,
                                                  var_hat_p_hat = var_hat_p_hat,
                                                  adjusted = F)

                           # result_AC_adjusted <- AC_method(p_hat = p_hat,
                           #                        var_hat_p_hat = var_hat_p_hat,
                           #                        n_psu = n_weights,
                           #                        adjusted = T)

                           result_CP <- CP_method(p_hat = p_hat,
                                                  var_hat_p_hat = var_hat_p_hat,
                                                  adjusted = F)

                           # result_CP_adjusted <- CP_method(p_hat = p_hat,
                           #                        var_hat_p_hat = var_hat_p_hat,
                           #                        n_psu = n_weights,
                           #                        adjusted = T)

                           list(mean_weight = mean_weight,
                                var_weight = var_weight,
                                coef_var = coef_var,
                                result_wspoissonTest = result_wspoissonTest,
                                result_wspoissonTest_midp = result_wspoissonTest_midp,
                                result_AC = result_AC,
                                # result_AC_adjusted = result_AC_adjusted,
                                result_CP = result_CP
                                # result_CP_adjusted = result_CP_adjusted
                                )
                         })) %>%
    select(-n_weights, -samples_per_weight) %>%
    unnest_wider(interval_results) %>%
    pivot_longer(cols = starts_with("result"),
                 names_prefix = "result_",
                 names_to = "method",
                 values_to = "htest") %>%
    mutate(conf.int = htest %>%
             map("conf.int") %>%
             map(~set_names(., c("l", "u")))) %>%
    unnest_wider(conf.int, names_sep = "_") %>%
    mutate(lower_error = conf.int_l > prevalence,
           upper_error = prevalence > conf.int_u,
           covered = !(lower_error | upper_error)) %>%
    select(-htest, -starts_with("conf.int"))

  experimental_results
}


experimental_results_8000_1 <- generate_results(simulation_design = simulation_design_8000_1, n_replications = 10000)
experimental_results_200_50 <- generate_results(simulation_design = simulation_design_200_50, n_replications = 10000)

write_rds(experimental_results_8000_1, "experimental_results_8000_1.rds")
write_rds(experimental_results_200_50, "experimental_results_200_50.rds")

generate_plot <- function(experimental_results, n_weights = 0, samples_per_weight = 0) {
  alpha_to_coef_var_conversion <-
      experimental_results %>%
      group_by(design, replication) %>%
      slice(1) %>%
      distinct() %>%
      group_by(dirichlet_alpha) %>%
      summarize(avg_coef_var = scales::percent(mean(coef_var))) %>%
      mutate(avg_coef_var = avg_coef_var %>% fct_reorder(dirichlet_alpha) %>% fct_rev())

  result_plot <-
    experimental_results %>%
    left_join(alpha_to_coef_var_conversion) %>%
    group_by(method, avg_coef_var, ICC) %>%
    summarize(lower_error_freq = mean(lower_error),
              upper_error_freq = mean(upper_error),
              coverage = mean(covered),
              .groups = "drop") %>%
    pivot_longer(-c(method, avg_coef_var, ICC)) %>%
    mutate(name = fct_recode(name,
                             "Lower Error Frequency" = "lower_error_freq",
                             "Upper Error Frequency" = "upper_error_freq",
                             "Coverage" = "coverage"),
           method = fct_recode(method,
                               "Agresti-Coull (Unadjusted)" = "AC",
                               "Agresti-Coull (Adjusted)" = "AC_adjusted",
                               "Clopper-Pearson (Unadjusted)" = "CP",
                               "Clopper-Pearson (Adjusted)" = "CP_adjusted",
                               "wsPoisson" = "wspoissonTest",
                               "wsPoisson with mid-p" = "wspoissonTest_midp")) %>%
    ggplot(aes(avg_coef_var, value, group = method, color = method)) +
    facet_grid(name ~ ICC, scale = "free_y", labeller = labeller(ICC = label_both)) +
    geom_line() +
    geom_point() +
    geom_hline(data = tibble(
      name =  c("Lower Error Frequency", "Upper Error Frequency", "Coverage"),
      yintercept = c(0.025, 0.025, 0.95)),
      mapping = aes(yintercept = yintercept), linetype = "dashed") +
    scale_x_discrete(name = "Average Coefficient of Variation") +
    ggtitle("Properties of Confidence Intervals for Surveys with Perfect Assays",
            subtitle = glue("{comma(max(experimental_results$replication))} replications of {comma(n_weights)} groups of size {comma(samples_per_weight)}")) +
    scale_y_continuous(name = NULL, labels = scales::percent) +
    scale_color_discrete(name = "Method") +
    cowplot::theme_minimal_grid() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(ncol=2))

  result_plot
}


results_plot_8000_1 <- generate_plot(experimental_results = experimental_results_8000_1 %>% filter(str_ends(method, "adjusted", negate = T)), n_weights = 8000, samples_per_weight = 1)
results_plot_200_50 <- generate_plot(experimental_results = experimental_results_200_50 %>% filter(str_ends(method, "adjusted", negate = T)), n_weights = 200, samples_per_weight = 50)


save_plot(results_plot_8000_1, ncol = 1, nrow = 3, filename = path("figures", "results_plot_8000_1", ext = "pdf"), base_asp = 2.5)
save_plot(results_plot_200_50, ncol = 3, nrow = 3, filename = path("figures", "results_plot_200_50", ext = "pdf"), base_asp =  1.5)
