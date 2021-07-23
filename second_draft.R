library(tidyverse)
library(extraDistr)
library(asht)
library(cowplot)
AC_method <- function(p_hat, var_hat_p_hat, conf.level = 0.95) {
  alpha <- 1 - conf.level
  z_val <- qnorm(1 - alpha / 2)
  c <- z_val * z_val / 2

  n_eff <- p_hat * (1 - p_hat) / var_hat_p_hat

  x_tilde <- p_hat * n_eff + c
  n_tilde <- n_eff + 2 * c
  p_tilde <- x_tilde / n_tilde

  ci <- p_tilde + c(-1,1) * z_val * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  attr(ci, "conf.level") <- conf.level

  rval <- list(conf.int = ci, estimate = p_tilde, p_hat = p_hat, var_hat_p_hat = var_hat_p_hat)
  class(rval) <- "htest"
  return(rval)
}

CP_method <- function(p_hat, var_hat_p_hat, conf.level = 0.95) {
  alpha <- 1 - conf.level

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

# multipdplyr appears to be slower for this
# library(multidplyr)
# cluster <- new_cluster(parallel::detectCores() / 2 - 2)
# cluster_library(cluster, "asht")
# cluster_library(cluster, c("asht", "dplyr"))
# cluster_assign(cluster, AC_method = AC_method)



n_replications <- 2
set.seed(200)
tmp <-
crossing(prevalence = 0.005,
         ICC = c(0.01, 0.25, 0.9),
         dirichlet_alpha = 10^(-2:2),
         n_weights = 8000,
         samples_per_weight = 1) %>%
  mutate(design = 1:n()) %>%
  select(design, everything()) %>%
  mutate(data = pmap(list(prevalence, ICC, dirichlet_alpha, n_weights, samples_per_weight),
                     function(prevalence, ICC, dirichlet_alpha, n_weights, samples_per_weight) {
                       alpha <- (1 - ICC) / ICC * prevalence
                       beta <- (1 - ICC) / ICC * (1 - prevalence)

                       tibble(iteration = rep(1:n_replications, each = n_weights),
                              weight = as.vector(t(rdirichlet(n_replications, alpha = rep(dirichlet_alpha, n_weights)))),
                              true_positivity = rbeta(n_weights * n_replications, alpha, beta),
                              n_tested = samples_per_weight) %>%
                         mutate(sample_positives = rbinom(n = n(), size = n_tested, prob = true_positivity)) %>%
                         mutate(sample_positivity = sample_positives / n_tested)
                     })) %>%
  unnest(data) %>%
  group_by(design, prevalence, ICC, dirichlet_alpha, n_weights, samples_per_weight, iteration) %>%
  # partition(cluster) %>%
  summarize(
    var_weight = var(weight),
    var_true_positivity = var(true_positivity),
    result_wspoissonTest = list(
      wspoissonTest(
        x = sample_positives,
        w = weight / n_tested,
        midp = F
      )
    ),
    result_wspoissonTest_midp = list(
      wspoissonTest(
        x = sample_positives,
        w = weight / n_tested,
        midp = T)
    ),
    result_AC = list(
      AC_method(
        p_hat = sum(sample_positivity * weight),
        var_hat_p_hat = sum((weight / n_tested)^2 * sample_positives))
    ),
    result_CP = list(
      CP_method(p_hat = sum(sample_positivity * weight),
                var_hat_p_hat = sum((weight / n_tested)^2 * sample_positives))
    ),
    .groups = "drop") %>%
  # collect() %>%
  pivot_longer(cols = starts_with("result"),
               names_prefix = "result_",
               names_to = "method",
               values_to = "htest")

results <-
tmp %>%
  mutate(conf.int = htest %>%
           map("conf.int") %>%
           map(~set_names(., c("l", "u")))) %>%
  unnest_wider(conf.int, names_sep = "_") %>%
  mutate(lower_error = conf.int_l > prevalence,
         upper_error = prevalence > conf.int_u,
         covered = !(lower_error | upper_error)) %>%
  select(-htest, -starts_with("conf.int"))

# IT WOULD BE BETTER TO GENERATE DATA< DO TESTS AND TOSS DATA IN ONE GO

result_plot <-
  results %>%
  group_by(method, dirichlet_alpha, ICC) %>%
  summarize(lower_error_freq = mean(lower_error),
            upper_error_freq = mean(upper_error),
            coverage = mean(covered),
            .groups = "drop") %>%
  pivot_longer(-c(method, dirichlet_alpha, ICC)) %>%
  mutate(name = fct_recode(name,
                           "Lower Error Frequency" = "lower_error_freq",
                           "Upper Error Frequency" = "upper_error_freq",
                           "Coverage" = "coverage"),
         method = fct_recode(method,
                             "Agresti-Coull (Unadjusted)" = "AC",
                             "Clopper-Pearson (Unadjusted)" = "CP",
                             "wsPoisson" = "wspoissonTest",
                             "wsPoisson with mid-p" = "wspoissonTest_midp"),
         dirichlet_alpha = as_factor(dirichlet_alpha)) %>%
  ggplot(aes(dirichlet_alpha, value, group = method, color = method)) +
  facet_grid(name ~ ICC, scale = "free_y", labeller = labeller(ICC = label_both)) +
  geom_line() +
  geom_point() +
  geom_hline(data = tibble(
    name =  c("Lower Error Frequency", "Upper Error Frequency", "Coverage"),
    yintercept = c(0.025, 0.025, 0.95)),
    mapping = aes(yintercept = yintercept), linetype = "dashed") +
  # scale_x_discrete(name = latex2exp::TeX("\\overset{Dirichlet $\\alpha$}{smaller $\\Rightarrow$ larger variance of weights}")) +
  scale_x_discrete(name = TeX(r'(\overset{Dirichlet $\alpha$}{larger $\Rightarrow$ smaller variance of weights})')) +
  ggtitle("Properties of Confidence Intervals for Surveys with Perfect Assays", subtitle = "8000 subjects per study, replicated 1000 times each") +
  scale_y_continuous(name = NULL, labels = scales::percent) +
  scale_color_discrete(name = "Method") +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom")


save_plot(filename = "results_2021-07-23.pdf", plot = result_plot, ncol = 3, nrow = 3, base_height = 2.5)
