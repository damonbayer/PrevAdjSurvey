library(tidyverse)
library(scales)
library(glue)
library(cowplot)
results_summary <-
  read_csv("code/simulations/confidence_interval_perfect_assay/combined_fixed_weights_many_summary.csv") %>%
  left_join(read_csv("code/simulations/confidence_interval_perfect_assay/experimental_design.csv") %>%
              mutate(group_distribution = case_when( # Error because I labeled things wrong when setting up ex design
                group_distribution == "low" ~ "high",
                group_distribution == "high" ~ "low",
                group_distribution == "uniform" ~ "uniform"
              ))) %>%
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

name_to_plot <- "Coverage"
n_groups_to_plot <- 8000

generate_plot <- function(name_to_plot, n_groups_to_plot){
  group_size <- if_else(n_groups_to_plot == 8000, 1, 200)

  generated_plot <-
    results_summary %>%
    mutate(
      tests = 10000,
      successes = as.integer(value * tests),
      failures = tests - successes) %>%
    filter(n_groups == n_groups_to_plot,
           name == name_to_plot) %>%
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
      formula = cbind(successes, failures) ~ splines::ns(x, 3),
      se = F,
      alpha = 0.5, linetype = "dashed"
    ) +
    cowplot::theme_minimal_grid() +
    scale_x_continuous(name = "Weight Coefficient of Variation") +
    scale_y_continuous(name = name_to_plot, label = ~percent(., accuracy = 1)) +
    scale_color_discrete(name = "Method") +
    theme(legend.position = "bottom") +
    ggtitle(label = glue("{name_to_plot} Properties for Simulations with {comma(n_groups_to_plot)} Groups of {comma(group_size)}"),
            subtitle = "Each Point = 10,000 Replications")
  generated_plot
  }




results_summary %>%
  select(n_groups, name) %>%
  distinct() %>%
  mutate(plot = map2(name, n_groups, ~generate_plot(.x, .y))) %>%
  mutate(file_name = map2_chr(name, n_groups, ~path("figures", str_c(str_replace_all(str_to_lower(.x), "\\s", "_"), .y, sep = "_"), ext = "pdf"))) %>%
  # mutate(walk2(plot, file_name, ~save_plot(filename = .y, plot = .x, ncol = 5, nrow = 3, base_asp = (11 / 5) / (8.5 / 3), base_height = 8.5 / 3)))
mutate(walk2(plot, file_name, ~save_plot(filename = .y, plot = .x, ncol = 5, nrow = 3, base_asp = (11 / 5) / (8.5 / 3), base_height = 4)))



