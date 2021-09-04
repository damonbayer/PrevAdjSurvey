library(tidyverse)
library(scales)
library(glue)
library(cowplot)
library(fs)

ed1 <- read_csv("code/simulations/confidence_interval_perfect_assay/experimental_design.csv")
ed2 <- read_csv("code/simulations/confidence_interval_perfect_assay/005_prev_experimental_design.csv")
sum1 <- read_csv("code/simulations/confidence_interval_perfect_assay/combined_fixed_weights_many_summary.csv")
sum2 <- read_csv("code/simulations/confidence_interval_perfect_assay/combined_fixed_weights_many_summary_005.csv")

results_summary <-
  left_join(sum1, ed1) %>%
  mutate(group_distribution = case_when( # Error because I labeled things wrong when setting up ex design
    group_distribution == "low" ~ "high",
    group_distribution == "high" ~ "low",
    group_distribution == "uniform" ~ "uniform"
  )) %>%
  bind_rows(left_join(sum2, ed2) %>%
              mutate(design = design + max(ed1$design))) %>%
  pivot_longer(cols = c(lower_error_freq, upper_error_freq, coverage)) %>%
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
                             "wsPoisson with mid-p" = "wspoissonTest_midp"),
         group_distribution = group_distribution %>%
           fct_relevel("high", "uniform", "low") %>%
           fct_relabel(str_to_title))

rm(ed1, ed2, sum1, sum2)


generate_plot <- function(name_to_plot, n_groups_to_plot, prev_to_plot){
  group_size <- if_else(n_groups_to_plot == 8000, 1, 200)

  generated_plot <-
    results_summary %>%
    filter(n_groups == n_groups_to_plot,
           name == name_to_plot,
           prev == prev_to_plot) %>%
    mutate(
      tests = 10000,
      successes = as.integer(value * tests),
      failures = tests - successes) %>%
    ggplot(aes(x = coef_var,
               y = successes / tests,
               color = method,
               group = method,
               # linetype = method,
               successes = successes,
               failures = failures)) +
    facet_grid(group_distribution ~ prop_groups_with_prev,
               scales = "free_y",
               labeller = labeller(
                 prop_groups_with_prev = function(x) x %>% as.numeric() %>% percent(accuracy = 1) %>% str_c("Weights with Non-0 Prev.\n", .),
                 group_distribution = function(x) str_c("Weights with Prev.\n", x, " Weights"))) +
    geom_hline(yintercept = case_when(name_to_plot == "Coverage" ~ 0.95,
                                      name_to_plot == "Upper Error Frequency" ~ 0.025,
                                      name_to_plot == "Lower Error Frequency" ~ 0.025),
               linetype = "dashed") +
    geom_point(size = 1, alpha = 0.1) +
    geom_smooth(
      method="glm",
      method.args=list(family="binomial"),
      formula = cbind(successes, failures) ~ splines::bs(x),
      se = F,
      alpha = 0.5, linetype = "dashed"
    ) +
    cowplot::theme_minimal_grid() +
    scale_x_continuous(name = "Weight Coefficient of Variation") +
    scale_y_continuous(name = name_to_plot, label = ~percent(., accuracy = 1)) +
    scale_color_discrete(name = "Method") +
    theme(legend.position = "bottom") +
    ggtitle(label = glue("{name_to_plot} Properties for Simulations with {percent(prev_to_plot, accuracy = 0.1)} Prevalence Among {comma(n_groups_to_plot)} Groups of {comma(group_size)}"),
            subtitle = "Each Point = 10,000 Replications")
  generated_plot
  }

generate_plot_2 <- function(name_to_plot, n_groups_to_plot, prev_to_plot){
  group_size <- if_else(n_groups_to_plot == 8000, 1, 200)

  generated_plot <-
    results_summary %>%
    filter(n_groups == n_groups_to_plot,
           name == name_to_plot,
           prev == prev_to_plot,
           method != "wsPoisson with mid-p",
           prop_groups_with_prev %in% c(0.05, 0.25, 0.75)) %>%
    mutate(
      tests = 10000,
      successes = as.integer(value * tests),
      failures = tests - successes) %>%
    ggplot(aes(x = coef_var,
               y = successes / tests,
               color = method,
               group = method,
               # linetype = method,
               successes = successes,
               failures = failures)) +
    facet_grid(group_distribution ~ prop_groups_with_prev,
               scales = "free_y",
               labeller = labeller(
                 prop_groups_with_prev = function(x) str_c("Groups with Non-0 Prev.\n", scales::comma(as.numeric(x) * n_groups_to_plot)),
                 group_distribution = function(x) str_c("Groups with Prev.\n", x, " Weights"))) +
    geom_hline(yintercept = case_when(name_to_plot == "Coverage" ~ 0.95,
                                      name_to_plot == "Upper Error Frequency" ~ 0.025,
                                      name_to_plot == "Lower Error Frequency" ~ 0.025),
               linetype = "dashed") +
    geom_point(size = 1, alpha = 0.25) +
    geom_smooth(
      method="glm",
      method.args=list(family="binomial"),
      formula = cbind(successes, failures) ~ splines::bs(x),
      se = F,
      alpha = 0.1, linetype = "dashed"
    ) +
    cowplot::theme_minimal_grid() +
    scale_x_continuous(name = "Weight Coefficient of Variation", labels = scales::percent) +
    scale_y_continuous(name = name_to_plot, label = ~percent(., accuracy = 1)) +
    scale_color_discrete(name = "Method") +
    theme(legend.position = "bottom") +
    ggtitle(label = glue("{name_to_plot} Properties for Simulations with {percent(prev_to_plot, accuracy = 0.1)} Prevalence Among {comma(n_groups_to_plot)} Groups of {comma(group_size)}"),
            subtitle = "Each Point = 10,000 Replications")
  generated_plot
}

generate_plot(name_to_plot = "Coverage", n_groups_to_plot = 50, prev_to_plot = 0.005)

results_summary %>%
  select(n_groups, name, prev) %>%
  distinct() %>%
  mutate(plot = pmap(
    list(name, n_groups, prev),
    ~generate_plot(name_to_plot = ..1,
                   n_groups_to_plot = ..2,
                   prev_to_plot = ..3))) %>%
  mutate(file_name = path("figures", str_c(str_replace_all(str_to_lower(name), "\\s", "_"), n_groups, str_replace_all(prev, "\\.", "_"), "full", sep = "_"), ext = "pdf")) %>%
  with(walk2(file_name, plot, ~save_plot(filename = .x, plot = .y, ncol = 5, nrow = 3, base_asp = (11 / 5) / (8.5 / 3), base_height = 4)))


results_summary %>%
  select(n_groups, name, prev) %>%
  distinct() %>%
  mutate(plot = pmap(
    list(name, n_groups, prev),
    ~generate_plot_2(name_to_plot = ..1,
                   n_groups_to_plot = ..2,
                   prev_to_plot = ..3))) %>%
  mutate(file_name = path("figures", str_c(str_replace_all(str_to_lower(name), "\\s", "_"), n_groups, str_replace_all(prev, "\\.", "_"), "reduced", sep = "_"), ext = "pdf")) %>%
  with(walk2(file_name, plot, ~save_plot(filename = .x, plot = .y, ncol = 3, nrow = 2, base_width = 3.5)))
